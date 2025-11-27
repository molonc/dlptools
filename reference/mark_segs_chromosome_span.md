# label CN segments based on relative chromosomal positions

Uses centromere and telomere coordinates to label where a segment sits
on a chromosome, relative to telomeres and centromeres.

## Usage

``` r
mark_segs_chromosome_span(
  segs_df,
  min_bound_distance = 5e+05,
  min_span_of_chrom = 0.9,
  min_span_of_arm = 0.9,
  version = c("hg19", "hg38"),
  acro_fix_whole_chrom = FALSE
)
```

## Arguments

- segs_df:

  dataframe of CN segments

- min_bound_distance:

  integer. Distance to adjacent feature to be considered associated with
  that feature.

- min_span_of_chrom:

  float. Proportion of the chromosome to cover to be considered a whole
  chromosome segment.

- min_span_of_arm:

  float. Proportion of arm to cover to be considered an arm segment.

- version:

  string. hg19 (default) or hg38

- acro_fix_whole_chrom.:

  boolean. Whether to reset acrocentric chromosome CN segments to
  "whole-chrom" if they span futher than the Q arm. Honestly, probably
  not useful

## Value

input dataframe, annotated with segment scale information. Primary
column of interest being seg_span_event.

## Details

Possible categories are:

1.  telomere bound (telo-bound) - segment touches a telomere

2.  centromere bound (centro-bound) - segment touches or crosses the
    centromere

3.  arm (arm) - segment spans a whole are (\*with conditions)

4.  whole chromosome (whole-chrom) - segment spans the entire chromosome
    (\*with conditions)

5.  intersitial (inter) - occuring within the chromosome, not touching
    the centromere, telomeres, and not big enough to be an entire arm.

You can set a min_bound_distance which reflects how close a feature
needs to be to be considered "touching". For example, we can considere a
segment telomere bound if within traditional 1 DLP bin by setting this
distance to 500k (default). This allows for some level of measurement
error.

Users can also set the proportion of the arm or whole chromosome a
segment needs to span to be considered either category. Default is
spaning 90% of either feature. Meaning, if the segment (end -
start)/arm_length is at least 90% of the arm_length, the segment is
considered an "arm" spanning segment.

This function runs several other functions including:
[`add_chromosome_length()`](https://molonc.github.io/dlptools/reference/add_chromosome_length.md),
[`add_centromere_locations()`](https://molonc.github.io/dlptools/reference/add_centromere_locations.md),
and
[`add_telomere_positions()`](https://molonc.github.io/dlptools/reference/add_telomere_positions.md).
