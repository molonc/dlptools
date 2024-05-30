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
