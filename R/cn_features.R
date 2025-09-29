#' Extract CN features following Wu et al
#'
#'
#' This function extracts copy number features in the style of the paper:
#'
#'  Wu et al. 2025. Single-cell copy number alteration signature analysis
#'  reveals masked patterns and potential biomarkers for cancer. bioRxiv.
#'
#'  https://www.biorxiv.org/content/10.1101/2025.03.02.641098v1
#'
#' They employ 4 base features, that they cross for 90 categories:
#' 1. CN states: 5 bins, 0-1, 2, 3, 4, 5+
#' 2. segment size: 3 bins, <5 MB, 5-10Mb, 10Mb,
#' 3. segment shape: 3 bins, LL (low left, low right segment), HH (high left,
#' high right segment), OT (other)
#' 4. segment change: 2 bins: AA (difference between surrounding segments <=
#' some critical value), or BB
#'
#' Some issues:
#' * whole chromosome amplifications/losses not captured
#' * chromosomes need to have at least 3 changes to be truly reflected in these
#' categories. Those with fewer are backfilled based on hardcoded rules
#' * for segment change, the paper says they considered changes > 2 as BB, but
#' in their actual code is set to 1. This function follows their code, but can
#' be altered.
#'
#' @param segs_df dataframe. CN segments
#' @param sample_col string. Name of the column with cell_id/other sample name
#' @param state_bin_max int. Maximum CN to consider for bins. All CNs of this
#' value and higher are grouped together. Default of 5 follows paper.
#' @param bin_breaks floats, how to break up segment sizes. Bins will be one
#' more than breaks. Defaults follow paper. Default is < 5Mb, 5--10Mb, > 10Mb
#' specified as c(5e6, 10e6 + 1). Internally, base::cut() is used, so 2 splits
#' produces 3 bins.
#' @param annotate_input boolean. return input dataframe annotating each
#' segment with the feature categories it falls into.
#' @param return_matrix boolean. Return a cell-by-feature matrix of counts.
#' @param ... can pass change_split_val to alter critical value for AA/BB split
#' @return default return is a tibble of feature counts for each cell id.
#' @export
extract_wu_features <- function(segs_df, sample_col = "cell_id", state_bin_max = 5, bin_breaks = NA, annotate_input = FALSE, return_matrix = FALSE, ...) {
  # paper code: https://github.com/XSLiuLab/single-cell-CNA-signature/blob/main/code/divide_feature.R

  # default segment size breakpoints
  if (is.na(bin_breaks)) {
    bin_breaks <- c(5e6, 10e6 + 1)
  }

  segs_df <- add_wu_seg_state_bins(
    segs_df = segs_df,
    state_bin_max = state_bin_max,
    bin_breaks = bin_breaks
  )

  segs_df <- add_wu_change_shape(segs_df = segs_df, ...)

  if (annotate_input) {
    return(segs_df)
  }

  # they cross everything before counting for 90 features
  feat_count <- segs_df |>
    dplyr::filter(
      !dplyr::if_any(c(cn_bin, seg_bin, seg_change, seg_shape), is.na)
    ) |>
    dplyr::group_by(
      .data[[sample_col]], cn_bin, seg_bin, seg_change, seg_shape,
      .drop = FALSE
    ) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      feat_cat = stringr::str_c(
        seg_bin, seg_shape, cn_bin, seg_change,
        sep = ":"
      )
    )

  if (return_matrix) {
    feat_mtx <- make_cellid_matrix(feat_count, "feat_cat", "n")
    return(feat_mtx)
  }

  return(feat_count)
}

#' see dlptools::extract_wu_features
#' @return tibble
add_wu_change_shape <- function(segs_df, change_split_val = 1) {
  # their code seems to use 1, not 2, but paper mentions 2

  # segment change and shape requires some grouping
  segs_df <- segs_df |>
    dplyr::group_by(cell_id, chr) |>
    dplyr::mutate(
      next_state = dplyr::lead(state),
      prev_state = dplyr::lag(state),
      seg_change = dplyr::case_when(
        (
          abs(state - next_state) > change_split_val |
            abs(state - prev_state) > change_split_val ~ "BB"
        ),
        # backfill start
        is.na(prev_state) ~ dplyr::if_else(
          abs(state - next_state) <= change_split_val,
          "AA", "BB"
        ),
        is.na(next_state) ~ dplyr::if_else(
          abs(state - prev_state) <= change_split_val,
          "AA", "BB"
        ),
        .default = "AA"
      ),
      seg_change = factor(seg_change, levels = c("AA", "BB")),
      seg_shape = dplyr::case_when(
        state < next_state & prev_state > state ~ "LL",
        state > next_state & prev_state < state ~ "HH",
        # backfill start/end
        is.na(prev_state) ~ dplyr::case_when(
          state < next_state ~ "LL",
          state > next_state ~ "HH",
          .default = "OT"
        ),
        is.na(next_state) ~ dplyr::case_when(
          state < prev_state ~ "LL",
          state > prev_state ~ "HH", # they actually have a bug in their code
          # here and always compare to the second row, when they meant to
          # compare to the second to last row
          .default = "OT"
        ),
        .default = "OT"
      ),
      seg_shape = factor(seg_shape, levels = c("LL", "HH", "OT"))
    ) |>
    dplyr::ungroup()

  return(segs_df)
}


#' see dlptools::extract_wu_features
#' @return tibble
add_wu_seg_state_bins <- function(
    segs_df, state_bin_max = 5, bin_breaks = c(5e6, 10e6 + 1)) {
  # names of segment bins, "S": small, L: large, M#: medium and depends on
  # number of intermediate categories
  bin_names <- c(
    "S",
    paste0(rep("M", length(bin_breaks) - 1), 1:(length(bin_breaks) - 1)),
    "L"
  )

  # looks like they put CN 0 & 1 into one bin
  segs_df <- segs_df |>
    dplyr::mutate(
      # segment bins
      seg_span = end - start + 1,
      seg_bin = cut(
        seg_span,
        c(0, bin_breaks, Inf),
        right = FALSE,
        labels = bin_names
      ),
      # state bins
      cn_bin = cut(
        state,
        c(0, 2:state_bin_max, Inf),
        labels = c(
          "S0-1",
          paste0(
            "S",
            c(2:(state_bin_max - 1), paste0(state_bin_max, "+"))
          )
        ),
        right = FALSE
      )
    )

  return(segs_df)
}


#' thin wrapper around sigminer feature tally
#'
#' Reshapes/renames our typical dataframes into one sigminer is happy with
#' and runs sigminer::sig_tally with the "W" method.
#'
#' Features are based on the paper:
#'
#' Wang et al. Copy number signature analysis tool and its application in
#' prostate cancer reveals distinct mutational processes and clinical outcomes.
#' PLOS Genetics. 2021.
#'
#' https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009557
#'
#' Some issues, I, Ben Furman, have with it are:
#'
#' * any oscillations of CN count as an oscillation, including something like
#' 2 - 500 - 2. Oscillations are supposed to reflect potential chromothripsis,
#' which a pattern like that is not likely to be from.
#' * segment sizes are binned based on log10 values. This means there are a lot
#' of CN size bins < 1 Mb, and DLP starts at 500Kb (and segments that small
#' should probably be filtered anyway). Thus, many of the segment size bins are
#' not used with DLP. This isn't a sigminer issue, but an issue with the scale
#' of DLP data.
#'
#' @param segs_df dataframe of CN segments.
#' @return matrix of feature counts.
#' @export
extract_sigminer_wang_features <- function(segs_df) {
  sm_tally <- sigminer::read_copynumber(
    dplyr::select(segs_df, cell_id, chr, start, end, state),
    seg_cols = c("chr", "start", "end", "state"),
    samp_col = "cell_id"
  ) |>
    sigminer::sig_tally(method = "W")

  return(sm_tally$nmf_matrix)
}


#' Count changes of state relative to ploidy.
#'
#' Marks CN segments as a gain or loss, relative to the mode ploidy of the
#' sample. Internally using [dlptools::mark_cn_relative_to_ploidy]. See that
#' function for argument details.
#'
#' @param segs_df dataframe. CN segments
#' @param sample_col string. Name of the column with cell_id/other sample name
#' @param annotate_input boolean. return input dataframe annotating each
#' @param return_matrix boolean. Return a cell-by-feature matrix of counts.
#' @return tibble/dataframe of counts
#' @export
extract_ploidy_cn_feature <- function(
    segs_df = NA,
    sample_col = "cell_id",
    annotate_input = FALSE,
    return_matrix = FALSE) {
  segs_df <- mark_cn_relative_to_ploidy(
    in_df = segs_df,
    df_type = "segs",
    sample_col = sample_col
  )

  if (annotate_input) {
    return(segs_df)
  }

  feat_count <- segs_df |>
    dplyr::group_by(.data[[sample_col]], .drop = FALSE) |>
    dplyr::count(cn_v_ploidy) |>
    dplyr::ungroup()

  if (return_matrix) {
    feat_mtx <- make_cellid_matrix(feat_count, "cn_v_ploidy", "n")
    return(feat_mtx)
  }

  return(feat_count)
}


#' count the segment-span-on-chromosome event types.
#'
#' Critical to this function is [dlptools::mark_segs_chromosome_span()]. It is
#' important to read and understand that function and its arguments.
#'
#' This function basically just calls [dlptools::mark_segs_chromosome_span()]
#' and summarizes the results. Arguments can be passed to that underlying
#' function. Passing no arguments means you are happy with the defaults. See
#' [dlptools::mark_segs_chromosome_span()] to understand what the defaults are.
#'
#' @param segs_df dataframe of CN segments
#' @param sample_col string. Name of the column with cell_id/other sample name
#' @param annotate_input boolean. return input dataframe annotating each
#' segment.
#' @param return_matrix boolean. Return a cell-by-feature matrix of counts.
#' @return tibble/dataframe of counts
#' @export
extract_segment_position_feature <- function(
    segs_df,
    sample_col = "cell_id",
    annotate_input = FALSE,
    return_matrix = FALSE,
    ...) {
  segs_df <- mark_segs_chromosome_span(segs_df, ...)

  if (annotate_input) {
    return(segs_df)
  }

  seg_span_counts <- segs_df |>
    dplyr::group_by(.data[[sample_col]], .drop = FALSE) |>
    dplyr::count(seg_span_event) |>
    dplyr::ungroup()

  if (return_matrix) {
    feat_mtx <- make_cellid_matrix(seg_span_counts, "seg_span_event", "n")
    return(feat_mtx)
  }

  return(seg_span_counts)
}

#' sizes of the segments
#'
#' This function is as simple as it sounds, end - start.
#'
#' Used as a setup for extracting process based features *á la*:
#' * Macintyre et al. 2018
#' * Drews et al. 2018
#'
#' Really, something like this should be used to generate values that then you
#' define categories for to count occurrences.
#'
#' Related functions include: [dlptools::extract_extract_changepoint] and
#' [dlptools::extract_breakpoints]
#'
#' @param segs_df dataframe. copy number segments for samples.
#' @param sample_col string. Name of column with cell/sample names
#' @param return string. "values" (default) or "counts". Values are the
#' observed values for cells, counts are the counts of these values in
#' pre-determined categories.
#' @return dataframe. sample ids and all observed segment sizes.
#' @export
extract_segment_sizes <- function(
    segs_df,
    sample_col = "cell_id",
    return = c("values", "counts")) {
  return_type <- match.arg(return)

  confirm_cols_present(c("start", "end"), segs_df)

  seg_sizes <- segs_df |>
    dplyr::mutate(seg_size = end - start + 1) |>
    dplyr::select(.data[[sample_col]], seg_size)

  if (return_type == "values") {
    return(seg_sizes)
  } else if (return_type == "counts") {
    seg_counts <- cut_categories_and_count(
      df = seg_sizes,
      targ_col = "seg_size",
      sample_col = sample_col,
      breaks = c(0, 5e6, 10e6, 20e6, 50e6, 100e6, Inf),
      labels = paste0(
        "ss:",
        c("<5Mb", "5-10Mb", "10-20Mb", "20-50Mb", "50-100Mb", "100+Mb")
      )
    )
    return(seg_counts)
  }
}

#' extract copy state values.
#'
#' This function is really only here to offer the counting in pre-set
#' categories, and provide alignment with other feature extraction types. It is
#' mostly a pointless function and you should probably do this yourself.
#'
#' Also, it is a question if copy state should even be used in
#' signature fitting, as this can occassionally lead to the same process being
#' artifically split by ploidy of samples. See
#' \href{https://www.nature.com/articles/s41586-022-04789-9}{Drews et al.}
#' for a discussion of this idea.
#'
#' @param segs_df dataframe. copy number segments for samples.
#' @param sample_col string. Name of column with cell/sample names
#' @param sample_col string. Name of column with CN state values
#' @param return string. "values" (default) or "counts". Values are the
#' observed values for cells, counts are the counts of these values in
#' pre-determined categories.
#' @export
extract_cn <- function(
    segs_df,
    sample_col = "cell_id",
    state_col = "state",
    return = c("values", "counts")) {
  return_type <- match.arg(return)

  cn_values <- dplyr::select(segs_df, .data[[sample_col]], .data[[state_col]])

  if (return_type == "values") {
    return(cn_values)
  } else if (return_type == "counts") {
    cn_counts <- cut_categories_and_count(
      df = cn_values,
      targ_col = state_col,
      sample_col = sample_col,
      breaks = c(0:6, 7, Inf),
      labels = c(paste0("CN:", 0:6), "CN:7+")
    )
    return(cn_counts)
  }
}


#' extract the CN change between adjacent segments.
#'
#' Change points are the change in copy number state between adjacent segments.
#' If one segment is 4 and the adjacent segment is 1, the change point is 3.
#'
#' Change points are based on the difference to the "left" adjacent segment,
#' when moving from BP 1 to the end of a chromosome. So if there are 3
#' segments: 4 - 1 - 2, the change points would be:
#' |1 - 4| and |2 - 1| resulting in: 3, 1
#'
#' For the first segment on a chromosome,
#' \href{https://www.nature.com/articles/s41586-022-04789-9}{Drews et al.}
#' compared it to a hypothetical diploid. So if the first segment on a
#' chromosome is 5, the change point would be 5 - 2 = 3. That's fine if the
#' base genome is diploid, but doesn't work so well for other ploidies, or
#' cases where you don't want to assume a diploid base case.
#'
#' `first_seg_correction` provides options over what to do. "diploid" for Drews
#' solution, "cn_mode" to compare to sample ploidy estimate based on mode, or
#' "ignore" to not count anything for first segments.
#'
#' @param segs_df dataframe. Sample copy number segments.
#' @param first_seg_correction string. Default: 'ignore'. Options include
#' "diploid" or "cn_mode".
#' @param sample_col string. Name of column with cell/sample names
#' @param chrom_col string. Name of column with chromosome names
#' @param cn_col string. Name of column with segment copy number states.
#' @param return string. "values" (default) or "counts". Values are the
#' observed values for cells, counts are the counts of these values in
#' pre-determined categories.
#' @param ... can pass arguments to [dlptools::segs_to_reads]
#' @return dataframe. Sample IDs and the observed breakpoint counts per scope.
#' @export
extract_changepoint <- function(
    segs_df,
    first_seg_correction = c("ignore", "cn_mode", "diploid"),
    sample_col = "cell_id",
    chrom_col = "chr",
    cn_col = "state",
    return = c("values", "counts"),
    ...) {
  first_seg_correction <- match.arg(first_seg_correction)
  return_type <- match.arg(return)

  if (first_seg_correction == "cn_mode") {
    mode_ploidies <- segs_to_reads(segs_df, ...) |>
      mode_ploidy(
        sample_col = sample_col,
        cn_col = cn_col, chrom_col = chrom_col
      )

    segs_df <- dplyr::left_join(
      segs_df, mode_ploidies,
      by = sample_col
    ) |>
      dplyr::rename(
        first_seg_comp_val = mode_ploidy
      )
  } else if (first_seg_correction == "diploid") {
    segs_df <- dplyr::mutate(segs_df, first_seg_comp_val = 2)
  } else if (first_seg_correction == "ignore") {
    segs_df <- dplyr::mutate(segs_df, first_seg_comp_val = NA)
  }

  change_points <- segs_df |>
    dplyr::group_by(.data[[sample_col]], .data[[chrom_col]]) |>
    dplyr::mutate(
      n_segs = dplyr::n(),
      left_seg = dplyr::lag(.data[[cn_col]]),
      left_seg = dplyr::if_else(
        is.na(left_seg),
        first_seg_comp_val,
        left_seg
      ),
      cn_change = dplyr::if_else(
        # set whole chromosomes without CN events to 0 changepoints.
        n_segs == 1,
        0,
        abs(.data[[cn_col]] - left_seg)
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(cn_change)) |>
    dplyr::select(.data[[sample_col]], cn_change)

  if (return_type == "values") {
    return(change_points)
  } else if (return_type == "counts") {
    changepoint_cts <- cut_categories_and_count(
      df = change_points,
      targ_col = "cn_change",
      sample_col = sample_col,
      breaks = c(0:5, Inf),
      labels = c(paste0("CP:", 0:4), "CP:5+")
    )

    return(changepoint_cts)
  }
}

#' convenience function to extracting breakpoints per window size
#'
#' just calling the generic [dlptools::extract_breakpoints], with some
#' pre-loaded options. See that function for details.
#' @return dataframe. Sample IDs and the observed breakpoint counts observed in
#' the specified window size.
#' @export
extract_bp_per_window <- function(
    segs_df,
    window_size = 10e6, # 10Mb is standard for most of these papers
    sample_col = "cell_id",
    chrom_col = "chr",
    genome_version = c("hg19", "hg38"),
    return = c("values", "counts")) {
  extract_breakpoints(
    segs_df = segs_df,
    scope = "windows",
    window_size = window_size,
    genome_version = genome_version,
    sample_col = sample_col,
    chrom_col = chrom_col,
    return = return
  )
}

#' convenience function to extracting breakpoints per chromosome arm
#'
#' just calling the generic [dlptools::extract_breakpoints], with some
#' pre-loaded options. See that function for details.
#' @return dataframe. Sample IDs and the observed breakpoint counts on
#' chromosome arms.
#' @export
extract_bp_per_arm <- function(
    segs_df,
    sample_col = "cell_id",
    chrom_col = "chr",
    genome_version = c("hg19", "hg38"),
    return = c("values", "counts")) {
  extract_breakpoints(
    segs_df = segs_df,
    scope = "arms",
    genome_version = genome_version,
    sample_col = sample_col,
    chrom_col = chrom_col,
    return = return
  )
}

#' extracting counts of breakpoints per user-defined scope
#'
#' Counting the number of breakpoints (i.e., transitions between copy number
#' segments) per arm or per window (typically 10Mb).
#'
#' @param segs_df dataframe. Copy number segments for samples
#' @param scope string. "windows" (default) or "arms", i.e., what to target for
#' counting
#' @param genome_version string. "hg19" (default) or "hg38"
#' @param window_size integer. How big of a window to use, if extracting counts
#' per a "windows" scope. Most publications use 10Mb, which is the default
#' (10e6)
#' @param sample_col string. Name of column with cell/sample names
#' @param chrom_col string. Name of column with chromosome names
#' @param return string. "values" (default) or "counts". Values are the
#' observed values for cells, counts are the counts of these values in
#' pre-determined categories.
#' @return dataframe. Sample IDs and the observed breakpoint counts per scope.
#' @export
extract_breakpoints <- function(
    segs_df,
    scope = c("windows", "arms"),
    genome_version = c("hg19", "hg38"),
    window_size = 10e6, # 10Mb is standard for most of these papers
    sample_col = "cell_id",
    chrom_col = "chr",
    return = c("values", "counts")) {
  scope <- match.arg(scope)
  return_type <- match.arg(return)
  genome_version <- match.arg(genome_version)

  confirm_cols_present("end", segs_df)

  chr_info <- suppressWarnings(load_chrom_info_file(version = genome_version))

  # DLP can output bins longer than the chromosome to maintain it's 500Kb
  # bin size
  segs_df_p <- dplyr::left_join(
    segs_df, chr_info,
    by = setNames(chrom_col, "chr")
  ) |>
    dplyr::mutate(
      end = dplyr::if_else(end > total_length, total_length, end)
    )

  # set up chromosome intervals, could then generalize to arms?
  if (scope == "windows") {
    intervals <- create_chrom_window_intervals(
      window_size = window_size,
      genome_version = genome_version
    )
  } else if (scope == "arms") {
    intervals <- create_chrom_arm_intervals(
      genome_version = genome_version
    )
  }

  bp_counts <- segs_df %>%
    split(.[[sample_col]]) |>
    purrr::map(
      .x = _, \(x) split(x[["end"]], x[[chrom_col]])
    ) |>
    purrr::imap_dfr(\(samp, samp_name) {
      bps <- purrr::imap(samp, \(chrom_ends, chrom_name) {
        counts <- cut(chrom_ends, breaks = intervals[[chrom_name]]) |>
          table() |>
          as.vector()

        unname(counts)
      }) |>
        purrr::list_c()

      tibble::tibble(
        !!sample_col := samp_name,
        breakpoints = bps
      )
    })

  if (return_type == "values") {
    return(bp_counts)
  } else if (return_type == "counts") {
    count_cats <- list(
      windows = list(
        breaks = c(0:3, Inf),
        labels = c(paste0("wBP:", 0:2), "wBP:3+")
      ),
      arms = list(
        breaks = c(0:6, Inf),
        labels = c(paste0("armBP:", 0:5), "armBP:6+")
      )
    )

    bp_cts <- cut_categories_and_count(
      df = bp_counts,
      targ_col = "breakpoints",
      sample_col = sample_col,
      breaks = purrr::pluck(count_cats, scope, "breaks"),
      labels = purrr::pluck(count_cats, scope, "labels")
    )

    return(bp_cts)
  }
}

#' extract the legths of chains of osscilating copy number segments
#'
#' Osscilating segments of copy numbers take from of N-M-N, e.g., 3-4-3. The
#' length of these chains reflects the number of oscillations. Chains are
#' evaluated in 3-segment intervals (the minimum required to recognize an
#' oscillation). Chromosomse with < 3 segments cannot have oscillations, and
#' will receive a 0.
#'
#' Chromosomes are evaluated in 3-segments sliding windows, with a window
#' receiving a 1 if the 3-set is an oscillation, and a 0 if not. Adjacent
#' oscillation sets are summed, non-osscilating sets are left as 0.
#'
#' Examples of Chain : length
#' * 2-3-4: 0
#' * 3-4-3 : 1
#' * 3-4-3-4: 2
#' * 3-4-3-2: 1 0 (3-4-3 is a chain, 4-3-2 is not)
#' * 2-3-4-5: 0 0
#'
#' An early paper arguing for the use of oscillations to detect chromothripsis
#' outlined that oscillations should be within 1 or 2 CN values. Using N-M-N
#' notation from above, this would mean that Ns should be the same value and M
#' would be within 1 or 2 copy values of N. Some pubilications relax these
#' constraints, and for example allow M to be any value other than N. This
#' means that 3-300-3 is counted as an oscillation by several publications.
#' Granted, this is an oscillation, but not likely to be one indicative of
#' chromothripsis.
#'
#' Others allow some flexiblity in the matching of the first and second N,
#' allowing them to be within some margin of difference. Componding these two
#' things leads to fairly extreme senarios, with some methods allowing chains
#' like 5-300-6 to count as an oscillation.
#'
#' These behaviours can be replicated here, but isn't recommended. By default,
#' Ns must match exactly, and M must be within 2 CNs.
#'
#' @param segs_df dataframe. Copy number segments for samples.
#' @param middle_bound integer. How many copy number values away the middle
#' segment can be from the end segments to count as an oscillation. Default 2
#' (e.g., 3-4-3 and 3-5-3 both count as oscillations, but 3-7-3 would not).
#' @param ends_bound integer. How many CN values apart the two ends of a
#' 3-segment set are allow to be to count as an oscillation. Default is 0,
#' i.e., they must be the same copy number.
#' @param sample_col string. Name of column with cell/sample names
#' @param chrom_col string. Name of column with chromosome names
#' @param cn_col string. Name of column with segment copy number states.
#' @param return string. "values" (default) or "counts". Values are the
#' observed values for cells, counts are the counts of these values in
#' pre-determined categories.
#' @return dataframe. sample IDs and observed chain lengths.
#' @export
extract_oscillations <- function(
    segs_df,
    middle_bound = 2,
    ends_bound = 0,
    sample_col = "cell_id",
    chrom_col = "chr",
    cn_col = "state",
    return = c("values", "counts")) {
  return_type <- match.arg(return)

  # TODO: this is pretty gross. Is it really better that filtering,
  # splitting, or something else?
  res <- segs_df %>%
    split(.[[sample_col]]) |>
    purrr::map(
      .x = _, \(x) split(x[[cn_col]], x[[chrom_col]])
    ) |>
    purrr::map(\(samp, samp_name) {
      purrr::map(samp, \(chrom) {
        count_oscillations(chrom)
      })
    })

  sample_osc_counts <- res |>
    purrr::list_flatten(name_spec = "{outer}&&&{inner}") |>
    tibble::enframe() |>
    tidyr::unnest_longer(value) |>
    tidyr::separate_wider_delim(name, names = c(sample_col, "chrom"), "&&&") |>
    dplyr::select(dplyr::any_of(sample_col), OSc = value)

  if (return_type == "values") {
    return(sample_osc_counts)
  } else if (return_type == "counts") {
    osc_cts <- cut_categories_and_count(
      df = sample_osc_counts,
      targ_col = "OSc",
      sample_col = sample_col,
      breaks = c(0, 1, 4, 10, Inf),
      labels = c("osc:0", "osc:1-3", "osc:4-9", "osc:10+")
    )

    return(osc_cts)
  }
}
