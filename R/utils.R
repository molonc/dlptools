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
convert_long_reads_to_wide <- function(reads_df) {
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
    dplyr::select("bin", "cell_id", "state") %>%
    tidyr::pivot_wider(names_from = bin, values_from = state)

  return(wide_states)
}

# no, I need this to go the other way
