#' decide if a 3 segment chain counts as an oscillation
#'
#' See [dlptools::extract_oscillations] for a discussion of oscillations and
#' recognizing them in sample copy number segments.
#'
#' @param left_cn int. CN value of segment on left side of the chain
#' @param middle_cn int. CN value of segment middle of the chain
#' @param left_cn int. CN value of segment on right side of the chain
#' @param middle_bound int|Inf. Default Inf. How different the middle CN is
#' allowed to be from either side to count as an oscillation. See
#' [dlptools::extract_oscillations] for a discussion of this
#' @param ends_bound int. Default 0. How different the ends of a 3-segment set
#' are allowed to be to count as an oscillation
#' @return boolean. True if 3 CNs constitute an oscillation, False if not.
decide_if_oscillation <- function(
    left_cn, middle_cn, right_cn,
    middle_bound = Inf,
    ends_bound = 0) {
  # are the ends the same, or within a defined degree of difference
  ends_match <- abs(left_cn - right_cn) <= ends_bound
  # is the middle and within any boundary of difference
  middle_different <- (
    middle_cn != right_cn & abs(middle_cn - right_cn) <= middle_bound
  )
  is_osc <- ends_match & middle_different
  return(is_osc)
}

#' count the length of oscillating chains of CN.
#'
#' Given a vector of CN values, detect and summarize lengths of any oscillating
#' chains that are present in the sequence, along with 0s for non-oscillating
#' sets.
#'
#' A full discussion of this function can be found in the documentation of.
#' [dlptools::extract_oscillations].
#'
#' Examples:
#' * 3-2-3 = 1
#' * 1-2-3 = 0
#' * 3-2-3-1-2-3 = 1 0
#' * 3-2-3-2 = 2
#' * 3-2-3-2-1 = 2 0
#'
#' @param cn_vals vector of ints. Segment state values. Ideally from the
#' chromosome of 1 sample.
#' @param middle_bound int|Inf. Default 2. How different the middle CN is
#' allowed to be from either side to count as an oscillation. See
#' [dlptools::extract_oscillations] for a discussion of this. Default of 2
#' means the middle CN must be within 2 copy values of the flanking segments.
#' @param ends_bound int. Default 0. How different the ends of a 3-segment set
#' are allowed to be to count as an oscillation
#' @return vector of ints.
#' @export
count_oscillations <- function(
    cn_vals,
    middle_bound = 2,
    ends_bound = 0) {
  if (length(cn_vals) < 3) {
    return(0) # or NA?
  }

  osc_count <- 0
  total_osc <- c()
  for (i in 3:length(cn_vals)) {
    # since we start on position 3, we are effectively at the "right_cn"
    is_osc <- decide_if_oscillation(
      left_cn = cn_vals[i - 2],
      middle_cn = cn_vals[i - 1],
      right_cn = cn_vals[i],
      middle_bound = middle_bound,
      ends_bound = ends_bound
    )
    # this would be a vanilla comparison, which most papers take, the
    # defaults of the function reflect this case where basically the ends
    # need to match, and the middle can be anything.
    # is_osc = cn_vals[i-2] == cn_vals[i] && cn_vals[i-1] != cn_vals[i]

    is_end_of_seq <- i == length(cn_vals)

    if (is_osc) {
      # accumulate a count for it being within an oscillation
      osc_count <- osc_count + 1
    } else if (!is_osc && !is_end_of_seq) {
      # if we are not in an oscillation and not at the terminal end
      if (osc_count > 0) {
        # store any oscillation we have accumulated
        total_osc <- c(total_osc, osc_count)
        # reset the oscillation count, as we are not in an oscillation
        osc_count <- 0
      }
      # add 0 for the current 3-set segment not being an oscillation
      total_osc <- c(total_osc, 0)
    }

    # set up cases to catch the last 3-set segment
    if (is_end_of_seq) {
      if (is_osc) {
        # if we are ending on the oscillation, store the count
        total_osc <- c(total_osc, osc_count)
      } else if (!is_osc) {
        # ending not on an oscillation
        if (osc_count > 0) {
          # store any accumulated oscillation
          total_osc <- c(total_osc, osc_count)
        }
        # finally, store a 0 for the current 3-set segment not being an
        # oscillation
        total_osc <- c(total_osc, 0)
      }
    }
  }
  return(total_osc)
}
