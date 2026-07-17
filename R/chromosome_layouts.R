#' loading UCSC chromosome length files
#' @param version default "hg19", can also load "hg38"
#' @return tibgble of chromosome, total length, etc.
#' @export
load_chrom_info_file <- function(version = c("hg19", "hg38")) {
  version <- match.arg(version)

  chrom_files <- list(
    hg19 = "hg19_chromInfo.txt.gz",
    hg38 = "hg38_chromInfo.txt.gz"
  )

  chrom_info <- get_package_file_path(chrom_files[version]) |>
    vroom::vroom(
      col_names = c("chr", "total_length", "misc"),
      show_col_types = FALSE
    ) |>
    dplyr::filter(
      # remove the unnecessary chromosomes
      stringr::str_detect(chr, "_|M", negate = TRUE)
    ) |>
    dplyr::mutate(
      chr = stringr::str_replace(chr, "chr", "")
    ) |>
    dplyr::select(-misc)

  return(chrom_info)
}

#' add chromosome lengths to a cn dataframe
#'
#' Handles chromosomes specified with "chr#" or just #/X/Y.
#'
#' @param cn_df dataframe of copy number information
#' @param version string. hg19 (default) or hg38
#' @param chrom_col string. name of column with chromosome information.
add_chromosome_length <- function(cn_df, version = c("hg19", "hg38"), chrom_col = "chr") {
  chrom_info <- load_chrom_info_file(version = version)

  if (is_chr_used_in_chroms(cn_df[[chrom_col]])) {
    chrom_info <- dplyr::mutate(chrom_info, chr = stringr::str_c("chr", chr))
  }

  cn_df <- dplyr::left_join(
    cn_df,
    chrom_info,
    by = dplyr::join_by({{ chrom_col }} == chr)
  )

  return(cn_df)
}


#' add centromere information to reads by chromosome
#'
#' See dlptools::read_and_prep_ucsg_cenrtomeres() for details of file origins.
#'
#' @param cn_df dataframe of cn states for read bins or segments
#' @param centro_file NULL or string to the path if you download yourself
#' @param hg19 boolean to target hg19 for loading
#' @param hg38 boolean to target hg38 for loading
#' @return input table with centromere information added by chromosome
#' @export
add_centromere_locations <- function(
  cn_df, centro_file = NULL, version = c("hg19", "hg38")
) {
  version_choice <- match.arg(version)

  centros <- load_ucsc_centromeres(
    centro_file = centro_file, version = version_choice
  )

  cn_df <- dplyr::left_join(
    cn_df,
    centros,
    by = dplyr::join_by(chr == chrom)
  )

  return(cn_df)
}

#' add boolean if bin overlaps with a centromere.
#'
#' Can optionally specify a padding to mark locations of bins as being close
#' enough to a centromere. Often bins near centromeres are corrupt in their
#' state calls. In the past, we have filtered within 3 Mb of centromeres.
#'
#' IMPORTANT! This adds the padding to each side of the centromere. So if you
#' specify a padding of 3 Mb, it will be within 3 Mb of the start and 3 Mb of
#' the end of the centromere.
#'
#' @param reads_df tibble of read data
#' @param padding int of number of BP to add to each side of the centromere
#' @param bin_start_col Default: start. column name of the start of bins.
#' @param bin_end_col Default: end column name of the end of bins.
#' @param version default 'hg19', or choose 'hg38' for locations of centromeres.
#' @return input table, but with a boolean 'within_centro' column added (and
#' potentially other centromere information columns, if needed)
#' @export
mark_bins_overlapping_centromeres <- function(
  reads_df,
  padding = 0,
  bin_start_col = "start",
  bin_end_col = "end",
  version = c("hg19", "hg38")
) {
  version_choice <- match.arg(version)
  columns_check <- all(c("centro_start", "centro_end") %in% colnames(reads_df))

  if (!columns_check) {
    warning(paste0(
      "centromere locations not present, adding locations for ",
      version_choice,
      ". If that is not correct, set version to 'hg38')."
    ))

    reads_df <- add_centromere_locations(cn_df = reads_df)
  }

  reads_df <- reads_df |>
    dplyr::mutate(
      centromere_padding = padding,
      overlaps_centro = (
        .data[[bin_start_col]] <= centro_end + padding &
          .data[[bin_end_col]] >= centro_start - padding
      )
    )

  return(reads_df)
}


#' load and prep default centromerefiles from UCSC
#'
#' Cytoband files originally obtained from UCSC golden path genomes and parsed
#' to collect the total centromere spans.
#'
#' hg19:
#' https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
#' hg38:
#' https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
#'
#' TODO: can I stash the file if it's already loaded once?
#'
#' @param centro_file NULL or string to the path if you download yourself
#' @param hg19 boolean to target hg19 for loading
#' @param hg38 boolean to target hg38 for loading
#' @return tibble of parsed centromere spans
#' @export
load_ucsc_centromeres <- function(
  version = c("hg19", "hg38"),
  centro_file = NULL
) {
  version_choice <- match.arg(version)

  if (is.null(centro_file)) {
    default_centro_files <- c(
      # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
      hg19 = "hg19_cytoBand.txt.gz",
      # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
      hg38 = "hg38_cytoBand.txt.gz"
    )

    centro_file <- get_package_file_path(default_centro_files[version_choice])
  }

  # import file and filter to features of interest
  centros <- vroom::vroom(
    centro_file,
    col_names = c("chrom", "start", "end", "arm_full_name", "feat"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(feat == "acen") |>
    dplyr::mutate(
      # removes the 'chr' from the chromosome name, e.g., chr1 -> 1
      chrom = stringr::str_extract(chrom, "([A-Z]|[0-9]+)")
    )
  centros <- centros |>
    dplyr::group_by(chrom) |>
    dplyr::mutate(
      centro_start = min(start),
      centro_end = max(end),
      centro_span = centro_end - centro_start,
      arm = stringr::str_extract(arm_full_name, "[a-z]")
    )

  # pivot out to one line per chromosome of centromere measurements
  centros <- centros |>
    dplyr::select(-c(arm_full_name, feat)) |>
    tidyr::pivot_wider(names_from = arm, values_from = c(start, end))

  return(centros)
}


#' load UCSC gap files for telomeres
#'
#' Honestly, this function is sort of pointless. The telomeres in both genome
#' versions are marked as 10Kb for p and q on all chromosomes. Telomeres are
#' variable and better estimates of their size exist, e.g.,
#' https://www.nature.com/articles/s41467-024-48917-7
#'
#' This function just serves as a way to set up easy marking of where CNs occur
#' on their respective chromosome.
#'
#' The only exception to the 10Kb size is chr17. It is not listed in the gap
#' file for hg19. chr17 is known to have small chromosomes. Here, I set the p
#' to 3000 Kb and q to 5000 Kb, inferred from the article referenced above.
#'
#' @param version string. hg19 (default) or hg38
#' @return tibble of telomere information
#' @export
import_telos_file <- function(version = c("hg19", "hg38")) {
  version_choice <- match.arg(version)

  default_telo_files <- c(
    # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
    hg19 = "hg19_gap.txt.gz",
    # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
    hg38 = "hg38_gap.txt.gz"
  )

  telos_file <- get_package_file_path(default_telo_files[version_choice])


  # col names by inspection and checking against here by looking up gap table
  # data format description
  # https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2018561924_ZAtdZC9CFEJw8BKiYeUjd2ImvhS7&clade=mammal&org=Human&db=hg19&hgta_group=allTables&hgta_track=hg19&hgta_table=gap&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=primaryTable&hgta_outFileName=
  telos <- vroom::vroom(
    telos_file,
    col_names = c(
      "bin", "chrom", "telostart", "teloend", "ix", "n", "size", "feat_type", "bridge"
    )
  ) |>
    dplyr::filter(feat_type == "telomere") |>
    dplyr::mutate(
      chrom = dplyr::if_else(
        # removes the 'chr' from the chromosome name, e.g., chr1 -> 1
        # this is to facilitate matching with our reads/segs files that just
        # have 0-0/X/Y
        chrom %in% c("chrX", "chrY"),
        stringr::str_extract(chrom, "[A-Z]"),
        stringr::str_extract(chrom, "[0-9]+")
      )
    )

  # now mark p and q telos
  telos <- telos |>
    dplyr::mutate(
      arm = dplyr::if_else(telostart == 0, "p", "q")
    ) |>
    dplyr::select(arm, chr = chrom, telostart, teloend) |>
    tidyr::pivot_wider(names_from = arm, values_from = c(telostart, teloend))

  # special case for chr17
  # for hg19, chr17 doesn't have telomeres listed, so need to fake some
  # telomere data.
  # marking p as 3000 and q as 5000 following:
  # https://www.nature.com/articles/s41467-024-48917-7
  chr17_total_len <- load_chrom_info_file() |>
    dplyr::filter(chr == 17) |>
    dplyr::pull(total_length)

  chr17_telo <- data.frame(
    chr = "17",
    telostart_p = 0,
    teloend_p = 3000,
    telostart_q = chr17_total_len - 5000,
    teloend_q = chr17_total_len
  )

  telos <- dplyr::bind_rows(telos, chr17_telo) |>
    dplyr::relocate(chr, telostart_p, teloend_p)

  return(telos)
}

#' add telomere positions to a CN df
#'
#' See [dlptools::import_telos_file()] but really just a convenience function.
#' Telomere positions in the UCSC gap file are all just 10Kb, except chr17,
#' which I set to 3Kb and 5Kb, again see linked function for details.
#'
#' @param cn_df dataframe of CN information, segments or read bins. Just needs
#' a 'chr' column.
#' @param version string. hg19 (default) or hg38, not that it makes a
#' difference.
#' @return input dataframe with new columns of telomere information
#' @export
add_telomere_positions <- function(cn_df, version = c("hg19", "hg38")) {
  telos <- import_telos_file(version = version)

  cn_df <- dplyr::left_join(cn_df, telos, by = "chr")

  return(cn_df)
}


#' label CN segments based on relative chromosomal positions
#'
#' Uses centromere and telomere coordinates to label where a segment sits on a
#' chromosome, relative to telomeres and centromeres.
#'
#' Possible categories are:
#' 1. telomere bound (telo-bound) - segment touches a telomere
#' 2. centromere bound (centro-bound) - segment touches or crosses the
#' centromere
#' 3. arm (arm) - segment spans a whole are (*with conditions)
#' 4. whole chromosome (whole-chrom) - segment spans the entire chromosome
#' (*with conditions)
#' 5. intersitial (inter) - occuring within the chromosome, not touching the
#' centromere, telomeres, and not big enough to be an entire arm.
#'
#' You can set a min_bound_distance which reflects how close a feature needs to
#' be to be considered "touching". For example, we can considere a segment
#' telomere bound if within traditional 1 DLP bin by setting this distance to
#' 500k (default). This allows for some level of measurement error.
#'
#' Users can also set the proportion of the arm or whole chromosome a segment
#' needs to span to be considered either category. Default is spaning 90% of
#' either feature. Meaning, if the segment (end - start)/arm_length is at least
#' 90% of the arm_length, the segment is considered an "arm" spanning segment.
#'
#' This function runs several other functions including:
#' [dlptools::add_chromosome_length()],
#' [dlptools::add_centromere_locations()], and
#' [dlptools::add_telomere_positions()].
#'
#' @param segs_df dataframe of CN segments
#' @param min_bound_distance integer. Distance to adjacent feature to be
#' considered associated with that feature.
#' @param min_span_of_arm float. Proportion of arm to cover to be considered an
#' arm segment.
#' @param min_span_of_chrom float. Proportion of the chromosome to cover to be
#' considered a whole chromosome segment.
#' @param version string. hg19 (default) or hg38
#' @param acro_fix_whole_chrom. boolean. Whether to reset acrocentric
#' chromosome CN segments to "whole-chrom" if they span futher than the Q arm.
#' Honestly, probably not useful
#' @return input dataframe, annotated with segment scale information. Primary
#' column of interest being seg_span_event.
#' @export
mark_segs_chromosome_span <- function(
  segs_df,
  min_bound_distance = 5e5, # given scale of DLP, this should probably be 1 bin width, at least
  min_span_of_chrom = 0.9,
  min_span_of_arm = 0.9,
  version = c("hg19", "hg38"),
  acro_fix_whole_chrom = FALSE
) {
  event_labels <- c(
    arm = "arm", whole = "whole-chrom", telo = "telo-bound",
    centro = "centro-bound", inter = "inter"
  )

  segs_df <- segs_df |>
    add_chromosome_length(version = version) |>
    add_centromere_locations(version = version) |>
    add_telomere_positions(version = version) |>
    dplyr::mutate(
      # find closest telomere to the CN event
      telo_p_dist = start - teloend_p,
      telo_q_dist = telostart_q - end,
      telo_dist = dplyr::case_when(
        # in telo or touching at both ends
        telo_p_dist <= 0 & telo_q_dist <= 0 ~ 0,
        # in or touching at one end
        telo_p_dist <= 0 | telo_q_dist <= 0 ~ 0,
        # closer to p
        telo_p_dist <= telo_q_dist ~ telo_p_dist,
        # closer to q
        telo_p_dist > telo_q_dist ~ telo_q_dist,
      ),
      telo_bound = telo_dist <= min_bound_distance,

      # find where CN sits relative to centromere
      centro_p_dist = centro_start - end,
      centro_q_dist = start - centro_end,
      centro_dist = dplyr::case_when(
        # crossing centro somewhere
        centro_p_dist <= 0 & centro_q_dist <= 0 ~ 0,
        # on the p side
        centro_p_dist >= 0 & centro_q_dist < 0 ~ centro_p_dist,
        # on th q side
        centro_p_dist < 0 & centro_q_dist >= 0 ~ centro_q_dist,
      ),
      spans_centro = end > centro_end & start < centro_start,
      centro_bound = centro_dist <= min_bound_distance,

      # span of chromosome
      spans_chrom = (end - start) / total_length >= min_span_of_chrom,

      # marking of arm events.
      p_arm_length = centro_start - teloend_p,
      q_arm_length = telostart_q - centro_end,
      spans_arm = dplyr::case_when(
        end < centro_start ~ (end - start) / p_arm_length >= min_span_of_arm,
        start > centro_end ~ (end - start) / q_arm_length >= min_span_of_arm,
        .default = FALSE
      ),

      # trying to be simple and Shih et al. 2023 10.1038/s41586-023-06266-3
      # inspired
      seg_span_event = dplyr::case_when(
        spans_chrom ~ event_labels["whole"], # "whole-chrom",
        centro_bound & telo_bound & !spans_centro ~ event_labels["arm"],
        !telo_bound & !centro_bound & !spans_chrom & spans_arm ~ event_labels["arm"],
        telo_bound & !centro_bound & !spans_chrom ~ event_labels["telo"], # "telo-bound",
        centro_bound | spans_centro ~ event_labels["centro"], # "centro-bound"
        .default = event_labels["inter"]
      )
    )

  if (acro_fix_whole_chrom) {
    # these chromosomes have p arms >10Mb. Short, but not crazy. Might not need
    # this fix.
    acro_chroms <- c("13", "14", "15", "21", "22")

    segs_df <- segs_df |>
      dplyr::mutate(
        # basically span of the q arm
        acro_span = dplyr::if_else(
          chr %in% acro_chroms & start > centro_start,
          (end - start) / (total_length - centro_end),
          NA
        ),
        seg_span_event = dplyr::case_when(
          !(chr %in% acro_chroms) ~ seg_span_event,
          acro_span >= q_span ~ event_labels["whole"],
          .default = seg_span_event
        )
      )
  }

  segs_df <- dplyr::mutate(
    segs_df,
    seg_span_event = factor(
      seg_span_event,
      levels = unname(event_labels)
    )
  )

  return(segs_df)
}

#' break a chromosome up into intervals of a defined window size
#'
#' @param window_size integer. The size of window to split the chromosome into.
#' @param genome_version string. "hg19" (default) or "hg38"
#' @return list. Named by chromosome, vectors of window starts.
#' @export
create_chrom_window_intervals <- function(
  window_size = 1e7,
  genome_version = c("hg19", "hg38"),
  return_type = c("granges", "tibble")
) {
  genome_version <- match.arg(genome_version)

  suppressWarnings(
    chr_info <- load_chrom_info_file(version = genome_version)
  )

  chr_lens <- chr_info$total_length
  names(chr_lens) <- chr_info$chr

  genome_tiles <- GenomicRanges::tileGenome(
    seqlengths = chr_lens,
    tilewidth = window_size,
    cut.last.tile.in.chrom = TRUE
  )

  genome_tiles$window_name <- paste0("w", seq(1, length(genome_tiles)))

  if (return_type == "granges") {
    genome_tiles
  } else {
    tibble::as_tibble(genome_tiles)
  }
  # intervals <- purrr::map(chr_info$total_length, \(total_length) {
  #   max_end <- total_length + window_size
  #   intervals <- seq(1, max_end, window_size)
  #   intervals
  # })
  # names(intervals) <- chr_info$chr
  # return(intervals)
}


#' create a list of intervals spanning chromosome arms
#'
#' Splits a chromosome at the middle of the centromere. Sets up intervals for
#' splitting each chromosome arm.
#'
#' @param genome_version string. "hg19" (default) or "hg38"
#' @return list. Named by chromosome, vectors of how to break a chromsome into
#' intervals of arms.
#' @export
create_chrom_arm_intervals <- function(genome_version = c("hg19", "hg38")) {
  chrom_layouts <- suppressWarnings(
    load_ucsc_centromeres(version = genome_version)
  )
  chrom_lengths <- suppressWarnings(
    load_chrom_info_file(version = genome_version)
  )

  chrom_info <- dplyr::left_join(
    chrom_lengths,
    chrom_layouts,
    by = dplyr::join_by("chr" == "chrom")
  )

  intervals <- purrr::pmap(
    dplyr::select(
      chrom_info,
      total_length, centro_start, centro_end
    ),
    \(total_length, centro_start, centro_end) {
      middle_centro <- centro_start + round(((centro_end - centro_start) / 2))
      c(1, middle_centro, total_length + 1)
    }
  )
  names(intervals) <- chrom_info$chr

  return(intervals)
}

load_chromosome_arms_info <- function(
  genome_version = c("hg19", "hg38"),
  return_type = c("granges", "tibble")
) {
  centros <- suppressWarnings(
    load_ucsc_centromeres(version = genome_version)
  )
  total_lens <- suppressWarnings(
    load_chrom_info_file(version = genome_version)
  )

  p_arms <- dplyr::select(centros, chrom, end = centro_start) |>
    dplyr::mutate(
      start = 1, arm = "p"
    )

  q_arms <- dplyr::select(centros, chrom, start = centro_end) |>
    dplyr::left_join(total_lens, by = dplyr::join_by(chrom == chr)) |>
    dplyr::mutate(arm = "q") |>
    dplyr::rename(end = total_length)

  all_arms <- dplyr::bind_rows(p_arms, q_arms) |>
    dplyr::arrange(chrom) |>
    dplyr::relocate(chrom, start, end, arm)

  if (return_type == "granges") {
    GenomicRanges::GRanges(all_arms)
  } else {
    all_arms
  }
}
