#' naturally sort a column from a dataframe.
#'
#' Common plotting issue, used mixed sort so things like, e.g., chromosomes
#' get sorted properly chr1, chr2, [...], chr10
#'
#' @param df a dataframe that contains the column to sort
#' @param col the column to sort
#' @return factor of the column with sorted levels
#' @export
factor_column_mixedsort <- function(df, col) {
  # common plotting issue to set up a factor with a mixedsort
  factor(df[[col]], levels = gtools::mixedsort(unique(df[[col]])))
}

#' internal for checking the names of chromosome columns in a dataframe.
#'
#' used for error catching.
#'
#' @param df a dataframe being manipulated in some function
#' @param exp_chr_name a string of the chromosome column name the function is
#' expecting to see.
#' @return bool TRUE if all good, FALSE with a message otherwise.
chr_name_check <- function(df, exp_chr_name) {
  if (!(exp_chr_name %in% names(df))) {
    print(paste0(
      "expecting chromosome names to be ",
      exp_chr_name,
      " in the mask file. Change the corresponding function arg if it's",
      " different."
    ))
    return(FALSE)
  }
  return(TRUE)
}


#' convert long format reads to wide format
#'
#' A common manipulation with reads files for various analyses is to reshape
#' long format reads data (each row is a 500kb bin with state values for each
#' cell) to wide format, with chr_start_end rows and cell_id columns.
#'
#' minimal required columns for input are: chr,start,end,cell_id,state
#'
#' @param reads_df is the reads table to convert.
#' @return wide format table
#' @export
#' @importFrom rlang .data
convert_long_reads_to_wide <- function(reads_df, state_col = "state") {
  # takes a csv of: chr,start,end,cell_id,state
  # and coverts it to: chrom_start_end,state,state
  # with an index column of location and columns of states for each cell

  wide_states <- reads_df %>%
    dplyr::mutate(
      bin = base::paste(
        .data$chr,
        base::as.integer(.data$start),
        base::as.integer(.data$end),
        sep = "_"
      )
    ) %>%
    dplyr::select("bin", "cell_id", {{ state_col }}) %>%
    tidyr::pivot_wider(names_from = bin, values_from = {{ state_col }})

  return(wide_states)
}


#' sort a table given a vector of cell_ids
#'
#' Typically used to sort a dataframe based on the plotted tip order to align
#' states heatmap/annotations/clone IDs to the plotted tree
#'
#' @param targ_df a table with cell_ids to sort
#' @param cell_order a vector of cell_ids in the desired order (e.g., pulled
#' from a ggplot of a tree)
#' @return table that has been sorted
#' @export
#' @importFrom rlang .data
sort_df_by_cell_order <- function(targ_df, cell_order) {
  sorted_df <- targ_df |>
    dplyr::arrange(base::match(.data$cell_id, cell_order))

  return(sorted_df)
}


#' extract sample ID from the typically formatted cell_ids
#'
#' expecting cell IDs as AT21350-A143952A-R10-C37 with the first position being
#' the sample ID.
#'
#' @param cell_id string of a cell_id or vector of cell IDs
#' @return string of the sample ID contained within
#' @export
pull_sample_id_from_cell_id <- function(cell_id) {
  pull_info_from_cell_id(cell_id, sample_id = TRUE)
}

#' generic extractor of info contained in cell ids
#'
#' @param cell_id string or vector of cells id
#' @param library_id boolean to extract library IDs
#' @param sample_id boolean to extract sample IDs
#' @return vector of requested information
#' @export
pull_info_from_cell_id <- function(
    cell_id, library_id = FALSE, sample_id = FALSE) {
  if (library_id && sample_id) {
    stop("you gotta pick, sample or library")
  }
  if (library_id) {
    idx <- 2
  } else if (sample_id) {
    idx <- 1
  }
  cell_info <- stringr::str_split(
    cell_id,
    pattern = "-", simplify = TRUE
  )[, idx]

  return(cell_info)
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


#' get plotted values bounds
#'
#' min, max, median of a column to generate a color palette for
get_column_metrics <- function(vals, min_max = FALSE) {
  if (min_max) {
    min <- min(vals, na.rm = TRUE)
    median <- median(vals, na.rm = TRUE)
    max <- max(vals, na.rm = TRUE)
    metrics_span <- c(min = min, median = median, max = max)
  } else {
    metrics_span <- quantile(vals, c(0.25, 0.5, 0.75), na.rm = TRUE)
  }

  return(metrics_span)
}
