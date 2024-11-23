# This script is inspired from and directly lifts some code from other versions
# of this script. Primarily inspired by a script passed to me by Hoa Tran, and
# other versions given to me by Daniel Lai, and the R package Signals.


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

#' format states for plotting in a heatmap
#'
#' @param states_df long format states information
#' @param state_col the column being plotted in the heatmap
#' @return tibble (or maybe just a matrix ready to go?)
#' @export
#' @importFrom rlang .data
format_states_for_hm <- function(states_df, state_col) {
  states_w <- states_df |>
    dplyr::select(cell_id, chr, start, end, state_col) |>
    convert_long_reads_to_wide(state_col = state_col)

  # convert to matrix
  states_mat <- base::as.matrix(dplyr::select(states_w, -.data$cell_id))
  rownames(states_mat) <- states_w$cell_id

  # sort columns by chromosome and bin start_end
  states_mat <- states_mat[, gtools::mixedsort(base::colnames(states_mat))]
  if (state_col != "BAF") {
    base::class(states_mat) <- "character"
  }

  if (state_col == "state") {
    states_mat[states_mat == "11"] <- "11+"
  }

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

#' creates a complex heatmap of a given matrix of states.
generate_state_hm <- function(
    states_mat, labels_fontsize = 8, plot_cols = STATE_COLORS, left_annot = NULL) {
  # set up a chromosome factor for column splits in the heatmap
  chroms <- create_chromosome_column_fct(states_mat)

  states_hm <- ComplexHeatmap::Heatmap(
    states_mat,
    col = plot_cols,
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
    na_col = "white",
    left_annotation = left_annot
  )

  return(states_hm)
}


#' generate a ComplexHeatmap::Heatmap image, either to console or file.
#'
#' Will return an image of the generated heatmap or dump the heatmap to a file.
#'
#' @param total_hm ComplexHeatmap::Heatmap of the combined tree and states.
#' @param file_name optional string of where to save a png image of the heatmap.
#' @return ComplexHeatmap::draw or nothing if a file is written.
#' @export
generate_hm_image <- function(
    hm, file_name = NULL,
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
    hm,
    padding = grid::unit(c(10, 2, 2, 2), "mm"),
    annotation_legend_side = "bottom",
    heatmap_legend_side = "bottom"
  )

  if (!is.null(file_name)) {
    dev.off()
  }
}

#' grab colors for various hm possibilities.
#'
#' Standard colors used by Signals and other people from the DLP world.
#' @export
fetch_heatmap_color_palette <- function(state_col, states_df) {
  color_choices <- list(
    "state" = STATE_COLORS,
    "A" = STATE_COLORS,
    "B" = STATE_COLORS,
    "BAF" = BAF_COLORS,
    "state_phase" = ASCN_PHASE_COLORS,
    "state_AS_phased" = ASCN_COLORS,
    "state_AS" = ASCN_COLORS
  )

  color_chosen <- purrr::pluck(
    color_choices, state_col,
    .default = STATE_COLORS
  )

  # check if the choice was ok
  plot_col_elements <- unique(states_df[[state_col]])

  if (state_col != "BAF" && length(plot_col_elements) > length(color_chosen)) {
    warning(
      paste0(
        "more elements that colors for ", state_col, " can plot them all.",
        " Defaulting to rainbow, but maybe don't plot this?"
      )
    )
    color_chosen <- grDevices::rainbow(length(plot_col_elements))
    names(color_chosen) <- gtools::mixedsort(plot_col_elements)
  }

  return(color_chosen)
}

#' builds the left-side annotations of the cells
build_left_annot <- function(anno_df = NULL, anno_cols_list = list()) {
  # and probably clones too
  if (!is.null(anno_df)) {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      df = as.data.frame(dplyr::select(anno_df, -c(cell_id))),
      na_col = "black",
      which = "row",
      col = anno_cols_list
    )
  } else {
    left_annot <- NULL
  }

  return(left_annot)
}

#' main hm building function
#'
#' anno_cols_list: list(Passage=c(`3`: #123456))
#' @export
plot_state_hm <- function(
    states_df, # long format data
    state_col, # column of data to plot
    phylogeny = NULL, # optional
    anno_df = NULL, # optional, can also specify columns in the dataframe
    anno_colors_list = list(), # for custom colors of annotations
    clones_df = NULL, # optional, can also specify columns in the dataframe
    anno_columns = NULL,
    file_name = NULL, # for direct saving to a file
    labels_fontsize = 8,
    ...) {
  # first, format the states for plotting
  states_mat <- format_states_for_hm(states_df, state_col)

  # deal with any annotations
  if (!is.null(anno_columns) && is.null(anno_df)) {
    anno_df <- dplyr::distinct(cell_id, anno_columns)
  }
  # ensure consistent ordering
  anno_df <- sort_df_by_cell_order(anno_df, rownames(states_mat))

  # deal with clones or consider it an annotation? They are a bit special
  # with the tree call out

  # deal with tree, and re-order states and annotations if so
  if (!is.null(phylogeny)) {
    cell_id_plot_order <- cell_id_order_as_plotted(phylogeny)
    states_mat <- states_mat[cell_id_plot_order, ]
    anno_df <- sort_df_by_cell_order(anno_df, cell_id_plot_order)

    tree_hm <- make_corrupt_tree_heatmap(phylogeny)
  } else {
    tree_hm <- NULL
  }


  # build left annotations, returns null if there is nothing
  left_annot <- build_left_annot(
    anno_df = anno_df,
    anno_cols_list = anno_colors_list
  )

  # determine plot colors for heatmap
  hm_colors <- fetch_heatmap_color_palette(state_col, states_df)

  # plot the heatmap
  state_hm <- generate_state_hm(
    states_mat,
    plot_cols = hm_colors,
    labels_fontsize = labels_fontsize,
    left_annot = left_annot
  )

  if (!is.null(tree_hm)) {
    total_hm <- tree_hm + state_hm
  } else {
    total_hm <- state_hm
  }

  generate_hm_image(total_hm, file_name = file_name, ...)
}
