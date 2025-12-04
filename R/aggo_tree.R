#' build tree with AGNES
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
#' @param reads_df data.frame. bin based reads data
#' @param by_ploidy_change. bool. Use CN change from cell mode ploidy as the
#' feature to cluster sample by. This can help group clones that are WGD of
#' lower ploidy clones.
#' @param state_col str. Column to use for the clustering. Default: "state"
#' @param sample_col str. Column of the sample labels. Default: "cell_id"
#' @param cut_k int. level to cut the tree at.
#' @return list. Two elements. $phylo: the tree; $clones: tibble of clone ids
#' of tip labels based on tree cutting.
#' @export
build_aggo_tree <- function(reads_df, cut_k = 8, state_col = "state", sample_col = "cell_id", by_ploidy_change = FALSE) {
  if (by_ploidy_change) {
    state_col <- "ploidy_change"
    reads_df <- reads_df |>
      mark_cn_relative_to_ploidy(
        sample_col = sample_col
      ) |>
      dplyr::mutate(
        ploidy_change = state - mode_ploidy
      )
  }

  # os_states_w <- convert_long_reads_to_wide(reads_df, state_col = state_col)

  # os_states_wm <- as.matrix(
  #   dplyr::select(os_states_w, -dplyr::all_of(sample_col))
  # )

  # rownames(os_states_wm) <- os_states_w[[sample_col]]

  # os_states_wms <- scale(os_states_wm) # needed?
  # os_states_dist <- dist(os_states_wms)


  os_states_wm <- format_states_for_hm(states_df = reads_df, state_col = state_col)

  agnes_clust <- cluster::agnes(dist(os_states_wm), method = "ward")

  agnes_clones <- tibble::tibble(
    # cell_id = os_states_w[[sample_col]],
    cell_id = rownames(os_states_wm), # do I know this is the right order?
    # oddly seems to be, though I don't know how.
    clone_id = cutree(agnes_clust, k = cut_k)
  )

  # agnes_phylo <- phylogram::as.phylo(as.dendrogram(agnes_clust))

  return(list(clones = agnes_clones, phylo = as.dendrogram(agnes_clust)))
}
