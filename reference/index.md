# Package index

## All functions

- [`ASCN_COLORS`](https://molonc.github.io/dlptools/reference/ASCN_COLORS.md)
  : colors for signals results of allele balances
- [`ASCN_PHASE_COLORS`](https://molonc.github.io/dlptools/reference/ASCN_PHASE_COLORS.md)
  : ASCN phase colors
- [`CNV_COLOURS`](https://molonc.github.io/dlptools/reference/CNV_COLOURS.md)
  : standard colors used in dlp plots
- [`DEFAULT_CONTINUOUS_COLOR_RANGE`](https://molonc.github.io/dlptools/reference/DEFAULT_CONTINUOUS_COLOR_RANGE.md)
  : min, median, max for a continuous color range
- [`STATE_COLORS`](https://molonc.github.io/dlptools/reference/STATE_COLORS.md)
  : alias of standard colors used in dlp plots
- [`add_centromere_locations()`](https://molonc.github.io/dlptools/reference/add_centromere_locations.md)
  : add centromere information to reads by chromosome
- [`add_chromosome_length()`](https://molonc.github.io/dlptools/reference/add_chromosome_length.md)
  : add chromosome lengths to a cn dataframe
- [`add_missing_bins_for_cells()`](https://molonc.github.io/dlptools/reference/add_missing_bins_for_cells.md)
  : Add back chr, start, end of bins missing from cells
- [`add_telomere_positions()`](https://molonc.github.io/dlptools/reference/add_telomere_positions.md)
  : add telomere positions to a CN df
- [`add_tip_ancestors_to_df()`](https://molonc.github.io/dlptools/reference/add_tip_ancestors_to_df.md)
  : add parent node labels to a dataframe of cell ids with states
- [`add_wu_change_shape()`](https://molonc.github.io/dlptools/reference/add_wu_change_shape.md)
  : see dlptools::extract_wu_features
- [`add_wu_seg_state_bins()`](https://molonc.github.io/dlptools/reference/add_wu_seg_state_bins.md)
  : see dlptools::extract_wu_features
- [`bme_vec()`](https://molonc.github.io/dlptools/reference/bme_vec.md)
  : picks beginning, middle, and end of a vector to handle when vectors
  that are too long are passed
- [`build_aggo_tree()`](https://molonc.github.io/dlptools/reference/build_aggo_tree.md)
  : build tree with AGNES clustering
- [`build_left_annot()`](https://molonc.github.io/dlptools/reference/build_left_annot.md)
  : builds the left-side annotations of the cells
- [`cell_states_to_strings()`](https://molonc.github.io/dlptools/reference/cell_states_to_strings.md)
  : collapse cell states to strings
- [`characterize_foreground_allele_state_changes()`](https://molonc.github.io/dlptools/reference/characterize_foreground_allele_state_changes.md)
  : mark types of allele specific foreground state changes
- [`characterize_foreground_total_state_changes()`](https://molonc.github.io/dlptools/reference/characterize_foreground_total_state_changes.md)
  : mark the types of change between a parent node state and the tip.
- [`check_args()`](https://molonc.github.io/dlptools/reference/check_args.md)
  : confirms arguments are compatible for the plotting wrapper
- [`check_if_state_columns_are_valid()`](https://molonc.github.io/dlptools/reference/check_if_state_columns_are_valid.md)
  : internal check for using the characterization function
- [`check_or_fetch_clone_palette()`](https://molonc.github.io/dlptools/reference/check_or_fetch_clone_palette.md)
  : confirm pallete given has enough colors or generate one
- [`check_the_vibe()`](https://molonc.github.io/dlptools/reference/check_the_vibe.md)
  : just a silly alias.
- [`chr_name_check()`](https://molonc.github.io/dlptools/reference/chr_name_check.md)
  : internal for checking the names of chromosome columns in a
  dataframe.
- [`compare_two_cells()`](https://molonc.github.io/dlptools/reference/compare_two_cells.md)
  : internal workhorse of pairwise_bin_difference function
- [`compute_tip_sibling_distances()`](https://molonc.github.io/dlptools/reference/compute_tip_sibling_distances.md)
  : measure string distances between sibling tips
- [`confirm_cols_present()`](https://molonc.github.io/dlptools/reference/confirm_cols_present.md)
  : confirm if columns exist in dataframe
- [`convert_long_reads_to_wide()`](https://molonc.github.io/dlptools/reference/convert_long_reads_to_wide.md)
  : convert long format reads to wide format
- [`count_default_process_feats()`](https://molonc.github.io/dlptools/reference/count_default_process_feats.md)
  : wrapper around classic process based features
- [`count_oscillations()`](https://molonc.github.io/dlptools/reference/count_oscillations.md)
  : count the length of oscillating chains of CN.
- [`create_chrom_arm_intervals()`](https://molonc.github.io/dlptools/reference/create_chrom_arm_intervals.md)
  : create a list of intervals spanning chromosome arms
- [`create_chrom_window_intervals()`](https://molonc.github.io/dlptools/reference/create_chrom_window_intervals.md)
  : break a chromosome up into intervals of a defined window size
- [`create_chromosome_column_fct()`](https://molonc.github.io/dlptools/reference/create_chromosome_column_fct.md)
  : create a sorted factor vector of chromosomes
- [`create_expected_bins()`](https://molonc.github.io/dlptools/reference/create_expected_bins.md)
  : create bins of given size for a genome
- [`cust_mode()`](https://molonc.github.io/dlptools/reference/cust_mode.md)
  : find most common element in a vector
- [`cut_categories_and_count()`](https://molonc.github.io/dlptools/reference/cut_categories_and_count.md)
  : count values that map to given categories
- [`cut_equal_segments()`](https://molonc.github.io/dlptools/reference/cut_equal_segments.md)
  : cut segments to equal widths across cells/samples.
- [`decide_if_oscillation()`](https://molonc.github.io/dlptools/reference/decide_if_oscillation.md)
  : decide if a 3 segment chain counts as an oscillation
- [`does_palette_cover_plot_vals()`](https://molonc.github.io/dlptools/reference/does_palette_cover_plot_vals.md)
  : internal to confirm color choice was ok
- [`expand_length_to_bins()`](https://molonc.github.io/dlptools/reference/expand_length_to_bins.md)
  : For a given length, create bins of a given size
- [`extract_bp_per_arm()`](https://molonc.github.io/dlptools/reference/extract_bp_per_arm.md)
  : convenience function to extracting breakpoints per chromosome arm
- [`extract_bp_per_window()`](https://molonc.github.io/dlptools/reference/extract_bp_per_window.md)
  : convenience function to extracting breakpoints per window size
- [`extract_breakpoints()`](https://molonc.github.io/dlptools/reference/extract_breakpoints.md)
  : extracting counts of breakpoints per user-defined scope
- [`extract_changepoint()`](https://molonc.github.io/dlptools/reference/extract_changepoint.md)
  : extract the CN change between adjacent segments.
- [`extract_cn()`](https://molonc.github.io/dlptools/reference/extract_cn.md)
  : extract copy state values.
- [`extract_oscillations()`](https://molonc.github.io/dlptools/reference/extract_oscillations.md)
  : extract the legths of chains of osscilating copy number segments
- [`extract_ploidy_cn_feature()`](https://molonc.github.io/dlptools/reference/extract_ploidy_cn_feature.md)
  : Count changes of state relative to ploidy.
- [`extract_segment_position_feature()`](https://molonc.github.io/dlptools/reference/extract_segment_position_feature.md)
  : count the segment-span-on-chromosome event types.
- [`extract_segment_sizes()`](https://molonc.github.io/dlptools/reference/extract_segment_sizes.md)
  : sizes of the segments
- [`extract_sigminer_wang_features()`](https://molonc.github.io/dlptools/reference/extract_sigminer_wang_features.md)
  : thin wrapper around sigminer feature tally
- [`extract_wu_features()`](https://molonc.github.io/dlptools/reference/extract_wu_features.md)
  : Extract CN features following Wu et al
- [`factor_column_mixedsort()`](https://molonc.github.io/dlptools/reference/factor_column_mixedsort.md)
  : naturally sort a column from a dataframe.
- [`fetch_continuous_color_ramp()`](https://molonc.github.io/dlptools/reference/fetch_continuous_color_ramp.md)
  : internal function for setting up heatmap continuous range colors
  chooses defaults, unless overwritten by user.
- [`fetch_heatmap_color_palette()`](https://molonc.github.io/dlptools/reference/fetch_heatmap_color_palette.md)
  : grab colors for various hm possibilities.
- [`fill_state_from_neighbours()`](https://molonc.github.io/dlptools/reference/fill_state_from_neighbours.md)
  : infer missing data using the upstream neighbour bin
- [`find_nearest_neighbours()`](https://molonc.github.io/dlptools/reference/find_nearest_neighbours.md)
  : find the nearest neighbour of each cell
- [`find_outlier_cells()`](https://molonc.github.io/dlptools/reference/find_outlier_cells.md)
  : Find outlier cells using a beta distribution
- [`format_sitka_tree()`](https://molonc.github.io/dlptools/reference/format_sitka_tree.md)
  : clean tree tip labels and drop any locus tips from sitka trees
- [`format_states_for_hm()`](https://molonc.github.io/dlptools/reference/format_states_for_hm.md)
  : format states for plotting in a heatmap
- [`generate_state_hm()`](https://molonc.github.io/dlptools/reference/generate_state_hm.md)
  : Internal function. Creates a complex heatmap of a given matrix of
  states.
- [`get_clone_id_label_positions()`](https://molonc.github.io/dlptools/reference/get_clone_id_label_positions.md)
  : get a row number of where to place each clone label.
- [`get_column_metrics()`](https://molonc.github.io/dlptools/reference/get_column_metrics.md)
  : get plotted values bounds
- [`get_dist_to_sibs()`](https://molonc.github.io/dlptools/reference/get_dist_to_sibs.md)
  : get distance of states to siblings for a tree tip
- [`get_package_file_path()`](https://molonc.github.io/dlptools/reference/get_package_file_path.md)
  : retrieve file path to package data files
- [`get_tip_parents()`](https://molonc.github.io/dlptools/reference/get_tip_parents.md)
  : resolve who the immediate parent nodes are for each tip.
- [`get_tips_that_avoid_redundant_comps()`](https://molonc.github.io/dlptools/reference/get_tips_that_avoid_redundant_comps.md)
  : get tips labels that will avoid duplicate sibling comparisons
- [`get_wgd_states()`](https://molonc.github.io/dlptools/reference/get_wgd_states.md)
  : Calculate WGD states based on major allele CN
- [`import_clones()`](https://molonc.github.io/dlptools/reference/import_clones.md)
  : read a clones t(c)sv file
- [`import_telos_file()`](https://molonc.github.io/dlptools/reference/import_telos_file.md)
  : load UCSC gap files for telomeres
- [`import_tree()`](https://molonc.github.io/dlptools/reference/import_tree.md)
  : read a tree from a file in newick format
- [`is_chr_used_in_chroms()`](https://molonc.github.io/dlptools/reference/is_chr_used_in_chroms.md)
  : detect if chromosomes labels include "chr"
- [`is_not_numeric_lacks_order()`](https://molonc.github.io/dlptools/reference/is_not_numeric_lacks_order.md)
  : confirms a vector is not numeric and is unorderd
- [`library_from_cell()`](https://molonc.github.io/dlptools/reference/library_from_cell.md)
  : extract library ID from the typically formatted cell_ids
- [`load_chrom_info_file()`](https://molonc.github.io/dlptools/reference/load_chrom_info_file.md)
  : loading UCSC chromosome length files
- [`load_mask_file()`](https://molonc.github.io/dlptools/reference/load_mask_file.md)
  : internal to control mask file loading.
- [`load_ucsc_centromeres()`](https://molonc.github.io/dlptools/reference/load_ucsc_centromeres.md)
  : load and prep default centromerefiles from UCSC
- [`make_cell_focused_matrix()`](https://molonc.github.io/dlptools/reference/make_cell_focused_matrix.md)
  : mostly internal for rearranging the pairwise DF to focus on a cell.
- [`make_cellid_matrix()`](https://molonc.github.io/dlptools/reference/make_cellid_matrix.md)
  : pivot a cell_id dataframe to a wide matrix.
- [`make_clone_palette()`](https://molonc.github.io/dlptools/reference/make_clone_palette.md)
  : generate a color palette for clone labels
- [`map_states_to_letters()`](https://molonc.github.io/dlptools/reference/map_states_to_letters.md)
  : convert a vector of states to letters
- [`mark_bins_overlapping_centromeres()`](https://molonc.github.io/dlptools/reference/mark_bins_overlapping_centromeres.md)
  : add boolean if bin overlaps with a centromere.
- [`mark_cn_relative_to_ploidy()`](https://molonc.github.io/dlptools/reference/mark_cn_relative_to_ploidy.md)
  : Mark if CN states are gains or losses relative to cell ploidy
- [`mark_mask_regions()`](https://molonc.github.io/dlptools/reference/mark_mask_regions.md)
  : add a mask column to a dataframe
- [`mark_segs_chromosome_span()`](https://molonc.github.io/dlptools/reference/mark_segs_chromosome_span.md)
  : label CN segments based on relative chromosomal positions
- [`medicc_profiles_to_foreground()`](https://molonc.github.io/dlptools/reference/medicc_profiles_to_foreground.md)
  : medicc profiles files to foreground assessment
- [`mode_ploidy()`](https://molonc.github.io/dlptools/reference/mode_ploidy.md)
  : find the mode CN per chromosome, then mode across the chromosomes
- [`output_hm_image()`](https://molonc.github.io/dlptools/reference/output_hm_image.md)
  : generate a ComplexHeatmap::Heatmap image, either to console or file.
- [`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md)
  : metric of pairwise differences between two cells.
- [`plot_bg_state_highlight()`](https://molonc.github.io/dlptools/reference/plot_bg_state_highlight.md)
  : generate a plot of just the background states, making the foreground
  white.
- [`plot_cell_cn_profile()`](https://molonc.github.io/dlptools/reference/plot_cell_cn_profile.md)
  : plot copy number profile plot of a cell.
- [`plot_dlp_chip()`](https://molonc.github.io/dlptools/reference/plot_dlp_chip.md)
  : Plot metrics across the chip
- [`plot_fg_state_highlight()`](https://molonc.github.io/dlptools/reference/plot_fg_state_highlight.md)
  : generate a plot of just the foreground states, making the background
  white.
- [`plot_gc()`](https://molonc.github.io/dlptools/reference/plot_gc.md)
  : visualize GC correction for a cell
- [`plot_heatmap_of_tip_changes()`](https://molonc.github.io/dlptools/reference/plot_heatmap_of_tip_changes.md)
  : make a plot of just the foreground events
- [`plot_nnd_outlier_cells()`](https://molonc.github.io/dlptools/reference/plot_nnd_outlier_cells.md)
  : visualize cell nearest neighbour distance values
- [`plot_read_bins_basic()`](https://molonc.github.io/dlptools/reference/plot_read_bins_basic.md)
  : create a tile plot of read state calls.
- [`plot_state_hm()`](https://molonc.github.io/dlptools/reference/plot_state_hm.md)
  : main hm building function
- [`pull_chr_from_col_name()`](https://molonc.github.io/dlptools/reference/pull_chr_from_col_name.md)
  : extract chromosome name from a bin name (chr_start_end)
- [`pull_info_from_cell_id()`](https://molonc.github.io/dlptools/reference/pull_info_from_cell_id.md)
  : generic extractor of info contained in cell ids
- [`read_medicc_profiles()`](https://molonc.github.io/dlptools/reference/read_medicc_profiles.md)
  : ingest a medicc profiles file
- [`read_medicc_tree()`](https://molonc.github.io/dlptools/reference/read_medicc_tree.md)
  : read medicc tree to phylo object.
- [`reads_to_segs()`](https://molonc.github.io/dlptools/reference/reads_to_segs.md)
  : convert reads table to segments table
- [`rle_states()`](https://molonc.github.io/dlptools/reference/rle_states.md)
  : convert states to run length encoding
- [`sample_from_cell()`](https://molonc.github.io/dlptools/reference/sample_from_cell.md)
  : extract sample ID from the typically formatted cell_ids
- [`segs_to_reads()`](https://molonc.github.io/dlptools/reference/segs_to_reads.md)
  : split segments into bins of a requested size.
- [`sort_df_by_cell_order()`](https://molonc.github.io/dlptools/reference/sort_df_by_cell_order.md)
  : sort a table given a vector of cell_ids
- [`weighted_ploidy()`](https://molonc.github.io/dlptools/reference/weighted_ploidy.md)
  : calculate sample ploidy with a weighted CN mean
