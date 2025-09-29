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
#' @return vector of the sample ID(s) contained within
#' @export
sample_from_cell <- function(cell_id) {
  pull_info_from_cell_id(cell_id, "sample")
}

#' extract library ID from the typically formatted cell_ids
#'
#' expecting cell IDs as AT21350-A143952A-R10-C37 with the second position being
#' the library ID.
#'
#' @param cell_id string of a cell_id or vector of cell IDs
#' @return vector of the library ID(s) contained within
#' @export
library_from_cell <- function(cell_id) {
  pull_info_from_cell_id(cell_id, "library")
}

#' generic extractor of info contained in cell ids
#'
#' @param cell_id string or vector of cells id
#' @param info "sample" or "library" to pull
#' @return vector of requested information
#' @export
pull_info_from_cell_id <- function(
    cell_id, info = c(NA, "sample", "library")) {
  info <- match.arg(info)
  info_indexes <- c(sample = 1, library = 2)

  if (is.na(info)) {
    stop("you gotta pick, info = 'sample' or 'library'")
  }

  cell_info <- stringr::str_split(
    cell_id,
    pattern = "-", simplify = TRUE
  )[, info_indexes[info]]

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

#' retrieve file path to package data files
#' @param file_name string of the file to target for loading
#' @return file path string
get_package_file_path <- function(file_name) {
  warning(paste0("loading package file: ", file_name), immediate. = TRUE)
  mask_f <- fs::path_package(
    "extdata", file_name,
    package = "dlptools"
  )

  if (is.null(mask_f)) {
    stop(past0("could not find package file: ", file_name))
  }

  return(mask_f)
}

#' pivot a cell_id dataframe to a wide matrix.
#'
#' Different functions may produce a dataframe of counts or something per
#' cell_id for various attributes. This converts it to a matrix with rownames
#' of cell_id and columns of some feature of interest from the dataframe that
#' can be pivoted out.
#'
#' @examples
#' cell_df <- data.frame(
#'   cell_id = c("cell_1", "cell_1", "cell_1"),
#'   some_col = c("A", "B", "C"),
#'   some_val = 1:3
#' )
#'
#' make_cellid_matrix(cell_df, name_col = "some_col", val_col = "some_val")
#'
#' #       A B C
#' # cell_1 1 2 3
#'
#' @param cell_df some dataframe with a column of cell_id and values for each
#' of those cells.
#' @param name_col the column to be pivoted out
#' @param val_col the column to use as values in the matrix
#' @return matrix
#' @export
make_cellid_matrix <- function(cell_df, name_col, val_col) {
  cell_mtx_init <- cell_df |>
    dplyr::select(cell_id, .data[[name_col]], .data[[val_col]]) |>
    tidyr::pivot_wider(
      id_cols = cell_id,
      names_from = .data[[name_col]],
      values_from = .data[[val_col]]
    )

  cell_mtx <- dplyr::select(cell_mtx_init, -cell_id) |> as.matrix()
  rownames(cell_mtx) <- dplyr::pull(cell_mtx_init, cell_id)
  return(cell_mtx)
}

#' find most common element in a vector
#'
#' For ties, just the first one in the list is returned.
#'
#' @param x vector of values
#' @return most common element of input vector
#' @export
cust_mode <- function(x, na_rm = FALSE) {
  if (na_rm) {
    x <- x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

#' detect if chromosomes labels include "chr"
#'
#' @param chrom_vec vector of chromosome labels, probably pulled from a
#' dataframe of cn information.
#' @return boolean. TRUE if "chr" found in the strings.
is_chr_used_in_chroms <- function(chrom_vec) {
  return(any(stringr::str_detect(unique(chrom_vec), "chr")))
}

#' confirm if columns exist in dataframe
#'
#' @param cols vector of strings. Names of columns you are searching for in the
#' dataframe.
#' @param df. The dataframe to look in for the columns
#' @return TRUE or raises a stop
confirm_cols_present <- function(cols, df) {
  if (!all(cols %in% colnames(df))) {
    col_names <- paste(cols, collapse = ", ")
    stop(paste0("need all columns: ", col_names))
  }
  return(TRUE)
}

#' count values that map to given categories
#'
#' @param df dataframe. Contains a column of values you want to count per sample
#' @param targ_col string. Column you want to put in categories and counted.
#' @param sample_col string. Name of the column with sample ids.
#' @param breaks vector of doubles. Bounds for the categories.
#' @param labels vector of strings. What to call the categories.
#' @return tibble
cut_categories_and_count <- function(
    df, targ_col, sample_col, breaks, labels) {
  if (length(breaks) != length(labels)) {
    stop("breaks and labels must be the same length.")
  }

  counts <- df |>
    dplyr::mutate(
      feat_cat = cut(
        .data[[targ_col]],
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE,
        right = FALSE
      )
    ) |>
    dplyr::group_by(.data[[sample_col]], feat_cat, .drop = FALSE) |>
    dplyr::count()

  return(counts)
}
