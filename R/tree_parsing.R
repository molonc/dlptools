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
