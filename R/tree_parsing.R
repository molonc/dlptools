#' resolve who the immediate parent nodes are for each tip.
#' @param phylogeny phylo object of tree of interest.
#' @return tibble of each tip and it's immediate parent node label.
#' @export
get_tip_parents <- function(phylogeny) {
  tip_ancestors <- phangorn::Ancestors(
    phylogeny,
    phylogeny$tip.label,
    type = "parent"
  )

  cell_node_ancestors <- tibble::tibble(
    cell_id = phylogeny$tip.label,
    parent_node = c(phylogeny$tip.label, phylogeny$node.label)[tip_ancestors]
  )

  return(cell_node_ancestors)
}

#' add parent node labels to a dataframe of cell ids with states
#' @param states_df bin based dataframe of cell ids and their states
#' @param phylogeny phylo object
#' @return input dataframe, but with parent node information added. Creates a
#' 'parent_node' column.
#' @export
add_tip_ancestors_to_df <- function(states_df, phylogeny) {
  # 1. get parent node of each tip
  cell_node_ancestors <- get_tip_parents(phylogeny)

  # 2. Add node information to the states DF
  states_df <- dplyr::left_join(
    states_df, cell_node_ancestors,
    by = "cell_id"
  )

  return(states_df)
}

#' read medicc tree to phylo object.
#'
#' Medicc forces the inclusion of a diploid ancestor for all cells and calls
#' that tip "diploid". This optionally drops that tip.
#'
#' @param tree_file file path to newick tree created by medicc
#' @param drop_dipliod boolean. Default TRUE. Remove diploid tip from tree.
#' @return phylo object
#' @export
read_medicc_tree <- function(tree_file, drop_diploid = TRUE) {
  phylo <- ape::read.tree(tree_file)

  if (drop_diploid) {
    phylo <- ape::drop.tip(phylo, "diploid")
  }

  return(phylo)
}
