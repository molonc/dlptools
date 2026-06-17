#' Calculate WGD states based on major allele CN
#'
#' Simple metric. If >50% of the genome has 2 or more copies of the "major"
#' allele ("A"; referred to as "FM2" in the paper), it's at least 1x WGD. If
#' there are 3 ("FM3") or more copies of A for >50% of the genome, then it's at
#' least 2x WGD. Otherwise, there is no support for WGD.
#'
#' Requires a [signals](https://shahcompbio.github.io/signals/) input, and is
#' based on the paper by
#' [McPherson et al. 2025](https://www.nature.com/articles/s41586-025-09240-3).
#'
#' Functions returns dataframe with "wgd_state" column, along with the fraction
#' estimates. Also calls on [weighted_ploidy()] and adds those results.
#'
#' This function will not use the sex chromosomes for the calculation.
#'
#' A potential issue with this approach is that larger chromosomes
#' disproportionately contribute to the result (e.g., sum of chr 1 and 2 is
#' more than many of the smaller chromosomes combined). Setting
#' `equalize_chromosomes` to TRUE will perform the calculations within each
#' chromosome, then assess if >50% of the chromosomes show FM2 >50%, or FM3 >50%
#' . This way, all chromosomes contribute equally to the inference of WGD.
#'
#' @param signals_df dataframe/tibble from signals, with A and B allele columns
#' @param cell_col string column name identifying cell IDs
#' @param start_col string column name identifying start of segments
#' @param end_col string column name identifying end of segments
#' @param chrom_col string column name for chromosomes
#' @param equalize_chromosomes boolean. See description.
#' @return tibble of per-sample ploidy estimates
#' @export
get_wgd_states <- function(
  signals_df,
  cell_col = "cell_id",
  start_col = "start",
  end_col = "end",
  chrom_col = "chr",
  equalize_chromosomes = FALSE
) {
  signals_df <- signals_df |>
    dplyr::filter(
      stringr::str_detect(
        stringr::str_to_lower(.data[[chrom_col]]),
        "x|y",
        negate = TRUE
      )
    )

  w_ploidy <- weighted_ploidy(
    signals_df,
    sample_col = cell_col,
    start_col = start_col,
    end_col = end_col,
    chrom_col = chrom_col,
    cn_col = "state" # expected for signals input.
  )

  if (!equalize_chromosomes) {
    wgd_state <- signals_df |>
      dplyr::group_by(.data[[cell_col]]) |>
      dplyr::mutate(
        span = .data[[end_col]] - .data[[start_col]],
      ) |>
      dplyr::summarise(
        mean_cn = mean(state),
        mode_cn = dlptools::cust_mode(state),
        total_genome = sum(span),
        frac_a_ge2 = sum(span[A >= 2]) / total_genome,
        frac_a_ge3 = sum(span[A >= 3]) / total_genome,
      )
  } else if (equalize_chromosomes) {
    wgd_state <- signals_df |>
      dplyr::group_by(.data[[cell_col]], .data[[chrom_col]]) |>
      dplyr::mutate(
        span = .data[[end_col]] - .data[[start_col]],
      ) |>
      dplyr::summarise(
        chr_mean_cn = mean(state),
        chr_mode_cn = dlptools::cust_mode(state),
        total_chrom = sum(span),
        chr_frac_a_ge2 = sum(span[A >= 2]) / total_chrom,
        chr_frac_a_ge3 = sum(span[A >= 3]) / total_chrom,
      ) |>
      dplyr::summarise(
        total_chroms = dplyr::n(),
        frac_a_ge2 = sum(chr_frac_a_ge2 > 0.5) / total_chroms,
        frac_a_ge3 = sum(chr_frac_a_ge3 > 0.5) / total_chroms,
        mean_ch = mean(chr_mean_cn),
        mode_cn = dlptools::cust_mode(chr_mode_cn)
      )
  }

  wgd_state <- wgd_state |>
    dplyr::mutate(
      wgd_state = dplyr::case_when(
        frac_a_ge3 > 0.5 ~ "2+ WGD",
        frac_a_ge2 > 0.5 ~ "1 WGD",
        .default = "no WGD"
      )
    ) |>
    dplyr::ungroup()


  wgd_res <- dplyr::left_join(
    dplyr::rename(w_ploidy, weighted_ploidy = ploidy),
    wgd_state,
    by = "cell_id"
  )

  return(wgd_res)
}
