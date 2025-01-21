read_medicc_profiles <- function(meddic_profiles_file) {
  med_profiles <- vroom::vroom(meddic_profiles_file) |>
    dplyr::rename(cell_id = sample_id) |>
    dplyr::mutate(
      chr = stringr::str_remove(chrom, "chr")
    )

  return(med_profiles)
}


profiles_to_foreground <- function(
    meddic_profiles_file, phylogeny, cn_type = c("total", "allele")) {
  cn_type <- match.arg(cn_type)

  med_profiles <- read_medicc_profiles(meddic_profiles_file)

  internal_profiles <- dplyr::filter(
    med_profiles,
    stringr::str_detect(cell_id, "internal")
  )

  tip_profiles <- dplyr::filter(
    med_profiles,
    !stringr::str_detect(cell_id, "internal|diploid")
  )

  if (cn_type == "total") {
    internal_profiles <- dplyr::select(
      internal_profiles,
      parent_node = cell_id,
      parent_state = total_cn,
      chrom, start, end
    )

    tip_profiles <- dplyr::rename(tip_profiles, state = total_cn)
  } else if (cn_type == "allele") {
    internal_profiles <- internal_profiles <- dplyr::select(
      internal_profiles,
      parent_node = cell_id,
      parent_a_state = cn_a,
      parent_b_state = cn_b,
      chrom, start, end
    )
    tip_profiles <- dplyr::rename(tip_profiles, state_a = cn_a, state_b = cn_b)
  }

  tip_profiles <- add_tip_ancestors_to_df(tip_profiles, phylogeny)

  tip_profiles <- dplyr::left_join(
    tip_profiles,
    internal_profiles,
    by = dplyr::join_by(parent_node, chrom, start, end)
  )

  summary_funcs <- list(
    "total" = purrr::partial(
      characterize_foreground_total_state_changes,
      state_col = "state",
      parent_state_col = "parent_state"
    ),
    "allele" = purrr::partial(
      characterize_foreground_allele_state_changes
      # TODO: actually implement this function
    )
  )

  tip_profiles <- summary_funcs[[cn_type]](tip_profiles)

  return(tip_profiles)
}
