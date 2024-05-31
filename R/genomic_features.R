#' add a mask column to a dataframe
#'
#' Determine if any start--end regions in the segs or reads dataframe overlaps
#' with regions in the mask file. Will add a 'mask' boolean column to the
#' dataframe to indicate whether the region should be masked.
#'
#' @param subject_df a table with chromosome, start, and end columns. Typically
#' the reads or segs dataframe created by [import_dlp_files()].
#' @param mask_f path to a mask file, see README.md ## Setup for details.
#' @param mask_chr_name the name of the chromosomes column in the mask file.
#' Default is "seqnames", which is the name in the provided file.
#' @param subject_chr_name name of the chromosome column in the reads/segs
#' table. Default is 'chr' which is what the import scripts call the column.
#' @export
mark_mask_regions <- function(
    subject_df, mask_f = NULL,
    mask_chr_name = "seqnames", subject_chr_name = "chr") {
  # in importing.R
  masks <- load_mask_file(mask_f = mask_f)

  # error checking of passed params
  # in utils.R
  c1 <- chr_name_check(masks, mask_chr_name)
  c2 <- chr_name_check(subject_df, subject_chr_name)
  if (!c1 || !c2) {
    return()
  }

  subject_grs <- GenomicRanges::GRanges(
    seqnames = subject_df[[subject_chr_name]],
    ranges = IRanges::IRanges(subject_df$start, subject_df$end)
  )

  mask_grs <- GenomicRanges::GRanges(
    seqnames = masks[[mask_chr_name]],
    ranges = IRanges::IRanges(masks$start, masks$end)
  )

  overlaps <- GenomicRanges::findOverlaps(subject_grs, mask_grs)

  # overlaps contain index of rows that are found in the mask
  subject_df["mask"] <- (
    base::seq_len(nrow(subject_df)) %in% base::data.frame(overlaps)$queryHits
  )

  return(subject_df)
}
