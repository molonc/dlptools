#' get distance of states to siblings for a tree tip
#'
#' Computes a basic string distance between a tip and all siblings in
#' a tree. A tip can have one or more siblings, and the mean distance
#' to them will be returned.
#'
#' e.g., in tree (A, (B, C)), B has one sibling C, but for A both B & C are
#' siblings. So the mean distance of A and C and B is returned.
#'
#' This function is intended to be used through it's wrapper
#' dlptools::compute_all_tip_sibling_distances()
#'
#' @param tip label of a tree tip
#' @param states_string string of the states for a tree tip
#' @param tree the full tree of data
#' @param state_str_df dataframe of all tips and the states strings
#' @param tip_label_col name of the column in the state_str_df containing the
#' names of the tips on the tree
#' @importFrom rlang .data
#' @export
get_dist_to_sibs <- function(
    tip_label, states_string, tree, state_str_df, tip_label_col = "tip_label") {
  # pull out node numbers of the tips
  tip_nodes <- tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]]

  # get all siblings for the target tip
  sib_node <- phangorn::Siblings(tree, tip_label)

  # find the siblings
  if (sib_node %in% tip_nodes) {
    # if the siblings node is a tip node, get that label
    sib_tips <- tree$tip.label[sib_node]
  } else {
    # else it's an internal node, in which case grab the tip descendants
    sib_tips <- tree$tip.label[phangorn::Descendants(tree, sib_node, type = "tips")[[1]]]
  }

  # can now get distances
  sib_state_strs <- state_str_df |>
    dplyr::filter(.data[[tip_label_col]] %in% sib_tips) |>
    dplyr::pull(states_string)

  sib_dists <- stringdist::stringdist(states_string, sib_state_strs)
  return(mean(sib_dists, na.rm = TRUE))
}


#' convert a vector of states to letters
#'
#' The point is to then use measure string distance between cells. With raw
#' states double digit values would count as 2 characters and throw things off,
#' with standard R-functions to measure string differences.
#'
#' @param states vector of state values
#' @return vector of letters corresponding to those states
#' @export
map_states_to_letters <- function(states) {
  if (any(states > 52)) {
    warning("states higher than 52 observed, will be capped at 51.")
    states[states >= 52] <- 51
  }

  state_to_letter_map <- c(LETTERS, letters)
  names(state_to_letter_map) <- 0:(length(state_to_letter_map) - 1)

  return(state_to_letter_map[as.character(states)])
}


#' collapse cell states to strings
#'
#' Converts a long format states dataframe into 1 row per cell with all
#' states converted to a letter string.
#'
#' @param states_df a dataframe with states for each bin for cells (or whatever
#' the tree tips are)
#' @param states_col name of the column with state values
#' @param cell_id_col name of the column with the cell ids (tree tip names)
#' @importFrom rlang .data
#' @export
cell_states_to_strings <- function(
    states_df, states_col = "state", cell_id_col = "cell_id") {
  tip_state_strings <- states_df |>
    # need to remap states to letters or 2-digit values count as 2 characters
    # and throw the whole comparison off
    dplyr::mutate(
      state_letter = map_states_to_letters(.data[[states_col]])
    ) |>
    dplyr::group_by(.data[[cell_id_col]]) |>
    dplyr::summarise(
      states_string = stringr::str_c(state_letter, collapse = "")
    ) |>
    dplyr::rename(tip_label = {{ cell_id_col }})

  return(tip_state_strings)
}

#' measure string distances between sibling tips
#'
#' Basically, this function is useful for asking if one tree groups more
#' similar tips together better than another tree.
#'
#' For sibling tips, measure the distance between their states, treating
#' the states across the genome as a string and obtaining a string distance.
#'
#' States for a cell id are first converted to letters (to prevent double digit
#' states from counting as 2 characters) and then made into a single string
#' across the genome for each cell. I.e.,
#' 2 2 2 3 3 3 10 -> C C C D D D K
#' see map_states_to_letters() for details.
#'
#' Then for each tip, it's sister tip is found and the string distance is
#' measured. If the sister to a tip is a clade, the mean distance to all tips
#' in the clade are found. E.g., in tree (A, (B, C)) the sister to A is both
#' B & C. See get_dist_to_sibs() for details.
#'
#' Finally, a mean distance across all sibling clades is computed and returned.
#'
#' @param states_df long format read bin state data
#' @param tree phylo object to be checked
#' @param states_col name of the column containing state data
#' @param cell_id_col name of the column containing the tip labels
#' @importFrom rlang .data
#' @export
compute_tip_sibling_distances <- function(
    states_df, tree, states_col = "state", cell_id_col = "cell_id") {
  all_tip_state_strings <- cell_states_to_strings(
    states_df = states_df,
    states_col = states_col,
    cell_id_col = cell_id_col
  )

  non_redundant_tips <- get_tips_that_avoid_redundant_comps(tree)

  targ_tips <- dplyr::filter(
    all_tip_state_strings,
    tip_label %in% non_redundant_tips
  )

  all_dists <- purrr::pmap_dbl(
    targ_tips,
    get_dist_to_sibs,
    tree = tree,
    state_str_df = all_tip_state_strings
  )

  sib_dist <- mean(all_dists)

  return(sib_dist)
}

#' just a silly alias.
#'
#' see compute_tip_sibling_distances()
#'
#' @export
check_the_vibe <- function(
    states_df,
    tree,
    states_col = "state",
    cell_id_col = "cell_id") {
  compute_tip_sibling_distances(
    states_df = states_df,
    tree = tree,
    states_col = "state",
    cell_id_col = "cell_id"
  )
}

#' get tips labels that will avoid duplicate sibling comparisons
#'
#' I.e., in a tree: (A, (B, C)) we don't want to compare B to C
#' and then C to B, we only need to do one of those comparisons.
#'
#' @param tree a phylo object
#' @return vector of tips labels that will lead to non-redundant
#' @export
get_tips_that_avoid_redundant_comps <- function(tree) {
  sibs <- phangorn::Siblings(tree, tree$tip.label)

  distinct_comps <- purrr::imap(sibs, \(s, n) {
    if (is.null(s)) {
      s <- 0 # place holder
    }
    data.frame(node = n, sib = s)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      min_v = pmin(sib, node, na.rm = ),
      max_v = pmax(sib, node),
      sorted_sib_node = paste0(min_v, max_v)
    ) |>
    dplyr::distinct(sorted_sib_node, .keep_all = TRUE)

  # tip labels to use that avoid redundant comparisons
  non_redundant_labs <- tree$tip.label[distinct_comps$node]

  return(non_redundant_labs)
}
