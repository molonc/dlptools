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
#' 3. segment shape: 3 bins, LL (low left, low right segment), HH (high left, high right segment), OT (other)
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
#' @param segs_df a dataframe of CN segments
#' @param state_bin_max int. Maximum CN to consider for bins. All CNs of this
#' value and higher are grouped together. Default of 5 follows paper.
#' @param bin_breaks floats, how to break up segment sizes. Bins will be one
#' more than breaks. Defaults follow paper. Default is < 5Mb, 5--10Mb, > 10Mb
#' specified as c(5e6, 10e6 + 1). Internally, base::cut() is used, so 2 splits
#' produces 3 bins.
#' @param annotate_input boolean. return input dataframe annotating each segment with the feature categories it falls into.
#' @param return_matrix boolean. Return a cell-by-feature matrix of counts.
#' @param ... can pass change_split_val to alter critical value for AA/BB split
#' @return default return is a tibble of feature counts for each cell id.
#' @export
extract_wu_features <- function(segs_df, state_bin_max = 5, bin_breaks = NA, annotate_input = FALSE, return_matrix = FALSE, ...) {
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
      cell_id, cn_bin, seg_bin, seg_change, seg_shape,
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


#' Could changes of state relative to ploidy.
#'
#' Marks CN segments as a gain or loss, relative to the mode ploidy of the
#' sample.
#'
#' @param segs_df a dataframe of CN segments
#' @param annotate_input boolean. return input dataframe annotating each
#' @param return_matrix boolean. Return a cell-by-feature matrix of counts.
#' @return tibble/dataframe of counts
#' @export
ploidy_relative_cn_state <- function(
    segs_df,
    sample_col = "cell_id",
    annotate_input = FALSE,
    return_matrix = FALSE) {
  sample_mode_ploidy <- segs_df |>
    reads_to_segs() |>
    mode_ploidy(sample_col = sample_col)

  segs_df <- dplyr::left_join(
    segs_df,
    sample_mode_ploidy,
    by = sample_col
  )

  segs_df <- segs_df |>
    dplyr::mutate(
      cn_v_ploidy = dplyr::case_when(
        state < mode_ploidy ~ "ploidy-loss",
        state > mode_ploidy ~ "ploidy-gain",
        state == mode_ploidy ~ "ploidy-match"
      )
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
