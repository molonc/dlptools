# This script is inspired from and directly lifts some code from other versions
# of this script. Primarily inspired by a script passed to me by Hoa Tran, and
# other versions given to me by Daniel Lai.
#
# this is left as a monolith module so it can easily be copied and modified if
# people need.


#' read a tree from a file in newick format
#'
#' @param tree_f a path to a tree file
#' @return phylo object
#' @export
import_tree <- function(tree_f) {
  s_tree <- ape::read.tree(tree_f)
  return(s_tree)
}

#' read a clones t(c)sv file
#'
#' @param clones_f a path to a t(c)sv file with cell_id and clone_id columns
#' @return tibble
#' @export
import_clones <- function(clones_f) {
  clones <- vroom::vroom(clones_f)
  # TODO: check for column names with a stop if wrong.
  if (!all(c("cell_id", "clone_id") %in% colnames(clones))) {
    stop("clones file requires columns: 'cell_id', 'clone_id'")
  }
  return(clones)
}

#' read an annotations file (t(c)sv)
#'
#' @param annotations_file path to file. Requires at least sample_id column
#' @return tibble
#' @export
import_annotations_df <- function(annotations_file) {
  anno_df <- vroom::vroom(annotations_file)

  if (!("sample_id" %in% colnames(anno_df))) {
    stop(
      "sample_id not in annotations df. Will not be able to add annotations"
    )
    return(NULL)
  }

  return(anno_df)
}

#' import wide format states data and flip to cell rows and bin columns.
#'
#' This is a file with the first column are bins (chrom_start_end) and the
#' other columns are cell ids. Internally, the function will flip this.
#'
#' @param states_f path to a state calls file
#' @return tibble with cell_id rows and bin columns
#' @export
import_wideformat_states_file <- function(states_f) {
  states <- vroom::vroom(states_f)

  # reset column names to account for index column
  base::colnames(states) <- c(
    "bin",
    base::colnames(states)[2:base::length(base::colnames(states))]
  )

  states_bin_wide <- convert_bins_to_columns_cells_to_rows(states)

  return(states_bin_wide)
}

#' flip bins from rows to columns and put cell ids as rows
#'
#' @param states_df tibble/dataframe with bins and rows and cell ids as columns
#' @return tibble of the inverse frame.
#' @export
#' @importFrom rlang .data
convert_bins_to_columns_cells_to_rows <- function(states_df) {
  cell_rows <- states_df |>
    tidyr::pivot_longer(cols = -c(.data$bin), values_to = "state") |>
    tidyr::pivot_wider(names_from = .data$bin, values_from = .data$state) |>
    dplyr::rename(cell_id = .data$name)

  return(cell_rows)
}

#' clean tree tip labels and drop any locus tips from sitka trees
#'
#' 'Locus tips' are from sitka and are locus values that end up on the
#' tip of trees. Also removes the "cell_" prefix from tip labels, which is
#' also a consequence of sitka.
#'
#' @param tree phylo object as read by ape::read.tree
#' @return phylo object cleaned of "cell_" notation
#' @export
format_sitka_tree <- function(tree) {
  locus_tips <- base::grep("locus", tree$tip.label, value = TRUE)
  tree <- ape::drop.tip(tree, locus_tips)

  tree$tip.label <- base::gsub("cell_", "", tree$tip.label)

  return(tree)
}

#' create a ggplot object of a phylogenetic tree
#'
#' Will make a plot of a tree that can be added to a heatmap object. The plot
#' itself is just a simple tree with no tip labels.
#'
#' @param phylo_tree a phylo object
#' @param clones (optional) adds clones as OTUs (TODO is to allow this to then
#' make colors)
#' @return ggplot object
#' @export
make_tree_plot_obj <- function(phylo_tree, clones = NULL) {
  if (!base::is.null(clones)) {
    cell_clone_groups <- get_clone_members(clones)
    # puts tip labels into named groups, with clone id as the name of the list
    phylo_tree <- ggtree::groupOTU(phylo_tree, cell_clone_groups)
  }

  # TODO: add colors like the original script
  tree_p <- ggplot2::ggplot(phylo_tree) +
    ggtree::geom_tree(size = 0.25) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ylim(0.5, length(phylo_tree$tip.label) + 0.5) +
    ggplot2::theme_void()

  return(tree_p)
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
make_corrupt_tree_heatmap <- function(phylo_tree, clones = NULL) {
  tree_ggplot <- make_tree_plot_obj(phylo_tree, clones = clones)

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
    purrr::map(dplyr::select, .data$cell_id)

  return(clone_members)
}


#' sort a table given a vector of cell_ids
#'
#' Typically used to sort a dataframed based on the plotted tip order to align
#' states heatmap/annotations/clone IDs to the plotted tree
#'
#' @param targ_df a table with cell_ids to sort
#' @param tip_order a vector of cell_ids in the desired order (e.g., pulled from
#' a ggplot of a tree)
#' @return table that has been sorted
#' @export
#' @importFrom rlang .data
sort_df_by_tip_order <- function(targ_df, tip_order) {
  sorted_df <- targ_df |>
    dplyr::arrange(base::match(.data$cell_id, tip_order))

  return(sorted_df)
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

#' prepare a states table for heatmap plotting
#'
#' Will take the states table (bin rows, cell ID columns) and sort it by the
#' plotted tree order, convert to a matrix (which complexheatmap needs), and
#' sort it by chromosome. It also converts state 11 to 11+, which is what it
#' really is.
#'
#' @param states table of bin rows, cell ID columns with state calls
#' @param tree_tips vector of the order that cell ids should be in (ideally
#' pulled from a ggplot tree object)
#' @return matrix of state calls.
#' @export
#' @importFrom rlang .data
format_states_for_hm <- function(states, tree_tips) {
  # order the states matrix same as the tree
  states <- sort_df_by_tip_order(states, tree_tips)

  # convert to matrix with cell_id rownames for ComplexHeatmap
  states_mat <- base::as.matrix(dplyr::select(states, -.data$cell_id))
  base::rownames(states_mat) <- states$cell_id

  # sort columns by chromosome and bin start_end
  states_mat <- states_mat[, gtools::mixedsort(base::colnames(states_mat))]
  base::class(states_mat) <- "character"
  states_mat[states_mat == "11"] <- "11+"
  return(states_mat)
}

#' extract chromosome name from a bin name (chr_start_end)
#'
#' an internal function for pulling information from a state matrix
#'
#' @param col_name string of chr_start_end (column name in the internally used
#' state dataframe)
#'
#' @return string of the chromosome
#' @export
pull_chr_from_col_name <- function(col_name) {
  chr_name <- stringr::str_split(col_name, pattern = "_", simplify = TRUE)[1]
  return(chr_name)
}

#' create a sorted factor vector of chromosomes
#'
#' This is mostly just for splitting the the states heatmap by chromosome.
#' Will be naturally sorted.
#'
#' @param states_mat matrix of cell_id named rows and bin (chr_start_end)
#' columns
#' @return factor of chromosomes with sorted levels
create_chromosome_column_fct <- function(states_mat) {
  chr_cols <- base::colnames(states_mat)
  chroms <- purrr::map_chr(chr_cols, pull_chr_from_col_name)
  chroms <- base::factor(
    chroms,
    levels = base::unique(gtools::mixedsort(chroms))
  )
  return(chroms)
}

#' get a row number of where to place each clone label.
#'
#' A single clone can appear on multiple clades of the tree or in one
#' monophyletic group. On the heatmap, we need to determine were to put a clone
#' ID label (e.g., A, B, C, etc), as we don't want a label for each cell id.
#' Typically, we label the largest clade with the clone ID, and then provide a
#' legend mapping clone colors to IDs. This function will find the "middle" cell
#' ID of the clade being labelled as the position to place the clone ID.
#'
#' @param clones table of cell_id, clone_id
#' @param only_largest_clone_group boolean. TRUE means only largest clone clade
#' will receive the clone ID label. FALSE means all clades wil get the clone
#' label.
#' @return tibble of clone_id, and a row_number to place it, which
#' ComplexHeatmap understands as the position to place the label on the left
#' side annotation.
#' @export
#' @importFrom rlang .data
get_clone_id_label_positions <- function(
    clones, only_largest_clone_group = TRUE) {
  # make a running number of identifying which clade groups for each clone ID
  # based on plotted row adjacency.
  clone_groups <- clones |>
    dplyr::mutate(row_num = dplyr::row_number()) |>
    dplyr::group_by(.data$clone_id) |>
    dplyr::mutate(
      con_group = base::cumsum(c(1, diff(.data$row_num) != 1))
    )

  if (only_largest_clone_group) {
    # get the largest group for each clone ID
    largest_clone_groups <- clone_groups |>
      dplyr::count(.data$con_group, name = "group_size") |>
      dplyr::slice_max(.data$group_size)

    clone_groups <- dplyr::inner_join(
      dplyr::ungroup(.data$clone_groups),
      largest_clone_groups,
      by = c("clone_id", "con_group")
    )
  }

  # find the middle of each clone group, which is where were will drop the
  # clone_id label in the heatmap
  clone_positions <- clone_groups |>
    dplyr::group_by(.data$clone_id, .data$con_group) |>
    dplyr::summarise(
      row_num = base::round(base::mean(.data$row_num))
    ) |>
    dplyr::ungroup()

  return(dplyr::select(clone_positions, clone_id, row_num))
}

#' extract sample ID from the typically formatted cell_ids
#'
#' expecting cell IDs as AT21350-A143952A-R10-C37 with the first position being
#' the sample ID.
#'
#' @param cell_id string of a cell_id
#' @return string of the sample ID contained within
#' @export
pull_sample_id_from_cell_id <- function(cell_id) {
  sample_id <- stringr::str_split(cell_id, pattern = "-", simplify = TRUE)[1]
  return(sample_id)
}

#' expand an annotations dataframe to be for each cell.
#'
#' Annotations are provided as per sample_id, this function expands those
#' annotations to each cell_id that they apply to. This is so ComplexHeatmap
#' can place the annotations for each cell in the heatmap.
#'
#' @param annos_df a table of sample_id and columns of the desired annotations.
#' @param cell_id_order the order that the cell ids are plotted in, typically
#' extracted from a ggplot of a tree (see [cell_id_order_as_plotted()])
#' @return tibble of annotations for each cell_id.
#' @export
make_ordered_cell_id_annos <- function(annos_df, cell_id_order) {
  if (base::is.null(annos_df)) {
    return(NULL)
  }
  # matching of annos to cells will all be based on sample_ids
  # df of cell_id and sample_id
  cell_sample_id <- purrr::map_chr(cell_id_order, pull_sample_id_from_cell_id)

  cell_sample_df <- dplyr::tibble(
    cell_id = cell_id_order,
    sample_id = cell_sample_id
  )

  ordered_cell_annos <- dplyr::left_join(
    cell_sample_df, annos_df,
    by = "sample_id"
  ) |>
    dplyr::select(-c("cell_id", "sample_id"))

  return(ordered_cell_annos)
}

#' create a ComlexHeatmap object of cell states
#'
#' This will create the typically DLP states heatmap of cells. It will add clone
#' labels and colors for each clade, and add any annotations given.
#'
#' @param states a table of cells states with a column of cell_ids, and other
#' columns of state bins (chrom_start_end).
#' @param cell_id_order a vector of cell_ids in the order that they should be
#' plotted in. Typically extracted from a ggplot of a tree (see
#' [cell_id_order_as_plotted()]).
#' @param clones a table with at least cell_id and clone_id labels
#' @param annos_df optional table with a `sample_id` column and any number of
#' annotations to add for those sample IDs.
#' @param labels_fontsize int of label size for the chromosomes on the X and
#' clone IDs on the y of the heatmap.
#' @return ComplexHeatmap::Heatmap object of states and any annotations.
#' @export
make_copynumber_heatmap <- function(
    states, cell_id_order, clones, annos_df = NULL, labels_fontsize = 8) {
  states_mat <- format_states_for_hm(states, cell_id_order)

  # set up a chromosome factor for column splits in the heatmap
  chroms <- create_chromosome_column_fct(states_mat)

  # build up annotations
  # first clones. Order same as tree
  clones_sort <- sort_df_by_tip_order(clones, cell_id_order)

  # TODO: add clone size to this for legend
  clone_label_pos <- get_clone_id_label_positions(
    clones_sort,
    only_largest_clone_group = FALSE
  )

  # now for any other annotations that are passed
  ordered_cell_annos <- make_ordered_cell_id_annos(
    annos_df, cell_id_order
  )

  left_annot <- ComplexHeatmap::HeatmapAnnotation(
    clones = clones_sort$clone_id,
    clone_label = ComplexHeatmap::anno_mark(
      labels = clone_label_pos$clone_id,
      at = clone_label_pos$row_num,
      link_width = grid::unit(2, "mm"),
      labels_gp = grid::gpar(fontsize = labels_fontsize)
    ),
    df = base::as.data.frame(ordered_cell_annos),
    na_col = "black",
    which = "row",
    annotation_legend_param = list(
      clones = list(nrow = 4)
    )
  )

  states_hm <- ComplexHeatmap::Heatmap(
    states_mat,
    col = CNV_COLOURS,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_quality = 5,
    column_split = chroms,
    column_title_side = "bottom",
    # might need to revisit when I size img
    column_title_gp = grid::gpar(fontsize = labels_fontsize),
    heatmap_legend_param = list(nrow = 4),
    left_annotation = left_annot
  )

  return(states_hm)
}

#' create the classic tree-state heatmap plot
#'
#' This will create a ComplexHeatmap::Heatmap with a tree on the left, and a
#' heatmap of states, with clone ID labels. It can also add annotations to the
#' left side, if provided. (not ideal set up) but you can also pass arguments
#' for [make_copynumber_heatmap()] and  [plot_combined_heatmap()] if needed. See
#' those functions for arguments not specified here.
#'
#' This function can either return a heatmap or dump it to a file, if a
#' file_name is provided.
#'
#' @param phylo_tree a phylo object of the tree to plot.
#' @param states a table of cells states with a column of cell_ids, and other
#' columns of state bins (chrom_start_end).
#' @param clones a table with at least cell_id and clone_id labels
#' @param annos_df optional table with a `sample_id` column and any number of
#' annotations to add for those sample IDs.
#' @param file_name optional string of where to save a png image of the heatmap.
#' @return either a ComplexHeatmap::Heatmap().draw() object or noting and writes
#' a file.
#' @export
create_tree_copynumber_heatmap <- function(
    phylo_tree, states_df, clones_df, annos_df = NULL, file_name = NULL, ...) {
  tree_hm <- make_corrupt_tree_heatmap(phylo_tree, clones_df)

  cell_id_plot_order <- cell_id_order_as_plotted(phylo_tree)

  copynumber_hm <- make_copynumber_heatmap(
    states_df, cell_id_plot_order, clones_df, annos_df, ...
  )

  total_hm <- tree_hm + copynumber_hm

  plot_combined_heatmap(total_hm = total_hm, file_name = file_name, ...)
}

#' generate a ComplexHeatmap::Heatmap image, either to console or file.
#'
#' Will return an image of the generated heatmap or dump the heatmap to a file.
#'
#' @param total_hm ComplexHeatmap::Heatmap of the combined tree and states.
#' @param file_name optional string of where to save a png image of the heatmap.
#' @return ComplexHeatmap::draw or nothing if a file is written.
#' @export
plot_combined_heatmap <- function(
    total_hm, file_name = NULL,
    png_height = 1600, png_width = 2800, png_res = 144) {
  if (!is.null(file_name)) {
    # open a file to dump it to
    grDevices::png(
      # strip any ext the user tries to give
      paste0(fs::path_ext_remove(file_name), ".png"),
      height = png_height, width = png_width,
      res = png_res
    )
  }
  ComplexHeatmap::draw(
    total_hm,
    padding = grid::unit(c(10, 2, 2, 2), "mm"),
    annotation_legend_side = "bottom",
    heatmap_legend_side = "bottom"
  )

  if (!is.null(file_name)) {
    dev.off()
  }
}
