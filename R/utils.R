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
