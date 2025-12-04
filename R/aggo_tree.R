#' build tree with AGNES clustering
#'
#' Agglomerative hierarchical clustering is a method of clustering that can
#' produce a tree-like structure. Here, we construct a distance matrix based,
#' typically, on copy number states of bins, then feed that distance matrix
#' to AGNES clustering.
#'
#' This function will also do a preliminary tree cutting to provide groups
#' within the tree, which could be interpreted as copy-number clones of cells.
#' You can always re-cut the returned tree with [stats::cutree()], which is all
#' that is used internally here.
#'
#' There are 3 options of values to use for clustering:
#'
#' 1. states
#' 2. by_cn_change: changes from one state to another state are marked as 1,
#' and all other bins as 0. In essence, it's marking the breakpoints and
#' clustering based on shared breakpoints.
#' 3. by_ploidy_change: finds mode ploidy of each cell, and then bin states are
#' converted into their difference from the mode ploidy.
#'
#' Options 2 and 3 seem to work well to group WGD cells of lower ploidy clones
#' together into clades. The also circumvent the problem that a straight
#' distance of states doesn't always correspond to biological reality. For
#' example, from a state of 2, a state of 3 and and state of 4 can both be 1
#' mutational step away, the former being an amplification event and the latter
#' being a WGD event. Thus, bins of 3 and 4 are potentially equally distant
#' from a bin of 2, making distance based metrics using CN states dubious (see
#' example).
#'
#' @examples
#' mm <- matrix(c(2, 3, 4), byrow = TRUE)
#' rownames(mm) <- c(paste0("cell-", LETTERS[1:3]))
#' dist(mm)
#'
#' @param reads_df data.frame. bin based reads data
#' @param by_cn_change bool. cluster based on changes in CN state along
#' chromosomes.
#' @param by_ploidy_change bool. Use CN change from cell mode ploidy as the
#' feature to cluster sample by. This can help group clones that are WGD of
#' lower ploidy clones.
#' @param state_col str. Column to use for the clustering. Default: "state"
#' @param sample_col str. Column of the sample labels. Default: "cell_id"
#' @param cut_k int. level to cut the tree at.
#' @return list. Two elements. $phylo: the tree; $clones: tibble of clone ids
#' of tip labels based on tree cutting.
#' @export
build_aggo_tree <- function(
    reads_df,
    cut_k = 8,
    state_col = "state",
    sample_col = "cell_id",
    chrom_col = "chr",
    by_ploidy_change = FALSE,
    by_cn_change = FALSE) {
  if (by_ploidy_change) {
    reads_df <- reads_df |>
      mark_cn_relative_to_ploidy(
        sample_col = sample_col
      ) |>
      dplyr::mutate(
        ploidy_change = state - mode_ploidy
      )

    state_col <- "ploidy_change"
  } else if (by_cn_change) {
    reads_df <- reads_df |>
      dplyr::group_by(.data[[sample_col]], .data[[chrom_col]]) |>
      dplyr::mutate(
        prev_bin_state = dplyr::lag(.data[[state_col]]),
        bin_change = dplyr::case_when(
          state == prev_bin_state ~ 0,
          state != prev_bin_state ~ 1,
          .default = NA
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(bin_change))

    state_col <- "bin_change"
  }

  os_states_wm <- format_states_for_hm(
    states_df = reads_df,
    state_col = state_col
  )

  agnes_clust <- cluster::agnes(dist(os_states_wm), method = "ward")

  agnes_clones <- tibble::tibble(
    cell_id = rownames(os_states_wm),
    clone_id = cutree(agnes_clust, k = cut_k)
  )

  return(
    list(
      clones = agnes_clones,
      dendro = as.dendrogram(agnes_clust)
    )
  )
}
