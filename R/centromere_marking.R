#' add centromere information to reads by chromosome
#'
#' See read_and_prep_ucsg_cenrtomeres() for details of file origins.
#'
#' @param centro_file NULL or string to the path if you download yourself
#' @param hg19 boolean to target hg19 for loading
#' @param hg38 boolean to target hg38 for loading
#' @return input table with centromere information added by chromosome
#' @export
add_centromere_locations <- function(
    reads_df, centro_file = NULL, hg19 = FALSE, hg38 = FALSE) {
  centros <- read_and_prep_ucsg_cenrtomeres(
    centro_file = centro_file, hg19 = hg19, hg38 = hg38
  )

  reads_df <- dplyr::left_join(
    reads_df,
    centros,
    by = dplyr::join_by(chr == chrom)
  )

  return(reads_df)
}

#' mark edge of bin as within centromere or not
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
#' @param bin_locatation_column which column to use as the location of a bin.
#' popular choices are the start of a bin, middle of a bin, or and of a bin.
#' @param hg19 bool. Not needed unless this function is also adding the
#' centromere information to the data frame. See add_centromere_locations()
#' which will be called if the correct centromere columns are not present.
#' @param hg38 bool. Same as hg19...but for hg38.
#' @return input table, but with a boolean 'within_centro' column added (and
#' potentially other centromere information columns, if needed)
#' @export
bin_within_centromere <- function(
    reads_df,
    padding = 0,
    bin_location_column = "start",
    hg19 = FALSE,
    hg38 = FALSE) {
  version_choices <- c(hg19 = hg19, hg38 = hg38)
  columns_check <- all(c("centro_start", "centro_end") %in% colnames(reads_df))

  if (!columns_check && all(version_choices == FALSE)) {
    stop(paste0(
      "Centromere information columns not present, use",
      " add_centromere_locations() or call again with hg19=TRUE or hg38=TRUE"
    ))
  } else if (all(version_choices)) {
    stop("you gotta pick one, hg19=TRUE or hg38=TRUE")
  } else if (!columns_check && any(version_choices)) {
    warning("adding centromere information to the dataframe")

    reads_df <- add_centromere_locations(
      reads_df = reads_df,
      hg19 = version_choices[["hg19"]],
      hg38 = version_choices[["hg38"]]
    )
  }

  reads_df <- reads_df |>
    dplyr::mutate(
      within_centro = dplyr::between(
        .data[[bin_location_column]],
        centro_start - padding,
        centro_end + padding
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
read_and_prep_ucsg_cenrtomeres <- function(
    centro_file = NULL,
    hg19 = FALSE,
    hg38 = FALSE) {
  version_choices <- c(hg19 = hg19, hg38 = hg38)

  if (all(version_choices) && is.null(centro_file)) {
    stop("Need to set hg19 or hg38 to TRUE, but not both")
  } else if (all(version_choices == FALSE)) {
    stop("I don't know what you want me to do. hg19 or hg38? Set one to true.")
  } else {
    targ_version <- names(version_choices)[which(version_choices == TRUE)]
    print(targ_version)
  }

  if (is.null(centro_file)) {
    default_centro_files <- c(
      # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
      hg19 = "hg19_cytoBand.txt.gz",
      # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
      hg38 = "hg38_cytoBand.txt.gz"
    )

    centro_file <- get_package_file_path(default_centro_files[targ_version])
  }

  # import file and filter to features of interest
  centros <- vroom::vroom(
    centro_file,
    col_names = c("chrom", "start", "end", "arm_full_name", "feat")
  ) |>
    dplyr::filter(feat == "acen") |>
    dplyr::mutate(chrom = dplyr::if_else(
      # removes the 'chr' from the chromosome name, e.g., chr1 -> 1
      chrom %in% c("chrX", "chrY"),
      stringr::str_extract(chrom, "[A-Z]"),
      stringr::str_extract(chrom, "[0-9]+")
    ))
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
