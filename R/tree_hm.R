# for plotting the tree in a heatmap.


#' create a ggplot object of a phylogenetic tree
#'
#' Will make a plot of a tree that can be added to a heatmap object. The plot
#' itself is just a simple tree with no tip labels.
#'
#' @param phylo_tree a phylo object
#' @param clones (optional) adds clones as OTUs (TODO is to allow this to then
#' make colors)
#' @return ggplot object
#' @importFrom rlang .data
#' @export
make_tree_plot_obj <- function(
    phylo_tree, clones_df = NULL, clone_palette = NULL, color_clones = FALSE) {
  if (color_clones && is.null(clone_palette)) {
    clone_palette <- make_clone_palette(unique(clones_df$clone_id))
  }

  if (color_clones) {
    cell_clone_groups <- get_clone_members(clones_df)

    # puts tip labels into named groups, with clone id as the name of the list
    phylo_tree <- ggtree::groupOTU(phylo_tree, cell_clone_groups)

    tree_aes <- ggplot2::aes(x, y, color = group)
  } else {
    tree_aes <- ggplot2::aes(x, y)
  }

  # TODO: add colors like the original script
  tree_p <- ggplot2::ggplot(phylo_tree, tree_aes) +
    ggtree::geom_tree(size = 0.25) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ylim(0.5, length(phylo_tree$tip.label) + 0.5) +
    ggplot2::theme_void()

  if (!is.null(clone_palette)) {
    tree_p <- tree_p + ggplot2::scale_color_manual(values = clone_palette)
  }

  return(tree_p)
}

#' create list of list of cell_ids that belong to each clone id.
#'
#' @param clones_df a tibble/df with cell_id, clone_id columns
#' @return named list of lists, with clone_id as names
#' @export
#' @importFrom rlang .data
get_clone_members <- function(clones_df) {
  clone_members <- clones_df |>
    dplyr::select(.data$cell_id, .data$clone_id) |>
    base::split(f = base::as.factor(dplyr::pull(clones_df, .data$clone_id))) |>
    purrr::map(dplyr::pull, .data$cell_id)

  return(clone_members)
}

#' generate a color palette for clone labels
#'
#' @param clone_ids unique vector of clone labels
#' @return named vector of hex colors, names are clone labels
#' @export
make_clone_palette <- function(clone_ids) {
  clone_ids <- gtools::mixedsort(clone_ids)
  n_clones <- length(clone_ids)
  if (n_clones < 26) {
    clone_palette <- pals::alphabet() |> sample(size = n_clones)
  } else {
    clone_palette <- grDevices::colorRampPalette(pals::stepped3())(n_clones)
  }

  names(clone_palette) <- clone_ids

  return(clone_palette)
}

#' grab cell ids in the order that they are plotted
#'
#' This is to align state calls and other things. It will make a ggplot of the
#' tree the same as the heatmap code does, and then pull out the cell ID in the
#' order that they are plotted
#'
#' @param phylo_tree a phylo object of the tree being used in the heatmap
#' @return vector of cell_ids (or really, whatever are the tip labels)
#' @export
#' @importFrom rlang .data
cell_id_order_as_plotted <- function(phylo_tree) {
  tree_ggplot <- make_tree_plot_obj(phylo_tree)
  ggplot_tree_tip_order <- tree_ggplot$data |>
    dplyr::tibble() |>
    dplyr::filter(.data$isTip) |>
    dplyr::arrange(.data$y) |>
    dplyr::pull(.data$label) |>
    base::rev() # bottom up assignment

  return(ggplot_tree_tip_order)
}

#' create a complex heatmap object with a tree.
#'
#' Basically creates an empty heatmap with a tree added as a left side
#' annotation.
#'
#' @param phylo_tree a phylo object of a tree.
#' @param clones (optional) adds clones as OTUs, will eventually allow for color
#' @return a ComlexHeatmap heatmap object with a tree.
#' @export
make_corrupt_tree_heatmap <- function(
    phylo_tree, clones_df = NULL, color_clones = FALSE, clone_palette = NULL) {
  tree_ggplot <- make_tree_plot_obj(
    phylo_tree,
    clones_df = clones_df,
    color_clones = color_clones,
    clone_palette = clone_palette
  )

  tree_annot_func <- ComplexHeatmap::AnnotationFunction(
    fun = function(index) {
      grid::pushViewport(grid::viewport(height = 1))
      grid::grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[5]])
      grid::popViewport()
    },
    var_import = list(tree_ggplot = tree_ggplot),
    width = grid::unit(4, "cm"),
    which = "row"
  )
  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree = tree_annot_func, which = "row", show_annotation_name = FALSE
  )

  n_cells <- base::sum(tree_ggplot$data$isTip)
  tree_hm <- ComplexHeatmap::Heatmap(
    base::matrix(nc = 0, nr = n_cells),
    left_annotation = tree_annot
  )

  return(tree_hm)
}
