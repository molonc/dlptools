#' generic function to import dlp files
#'
#' Will create a dataframe of all files of a certain type. This is for the three
#' common dlp files that you might want to load: metrics, segments, or reads.
#'
#' @param dlp_sc_dir the parent directory where the dlp outputs are saved (see
#' the README.md for expected structure)
#' @param file_type options are: "metrics", "segs", "reads"
#' @return a tibble of the information from all dlp outputs combined
#' @export
#' @importFrom rlang .data
import_dlp_files <- function(
    dlp_sc_dir,
    file_type = c("metrics", "segs", "reads")) {
  glob_patterns <- c(
    metrics = "*/annotation/*metrics.csv.gz",
    segs = "*/hmmcopy/segments.csv.gz",
    reads = "*/hmmcopy/reads.csv.gz"
  )

  dlp_fs <- fs::dir_ls(
    path = dlp_sc_dir,
    glob = glob_patterns[file_type],
    recurse = TRUE
  )

  dlp_dat <- purrr::map_dfr(dlp_fs, vroom::vroom, id = "file_path")

  dlp_dat <- dplyr::mutate(
    dlp_dat,
    SCID = stringr::str_extract(.data$file_path, "SC-[0-9]+")
  )

  if ("chr" %in% names(dlp_dat)) {
    dlp_dat$chr <- factor_column_mixedsort(dlp_dat, "chr")
  }

  return(dlp_dat)
}


#' internal to control mask file loading.
#' @param mask_f path to the mask file to load.
load_mask_file <- function(mask_f) {
  masks <- vroom::vroom(mask_f)

  return(masks)
}
