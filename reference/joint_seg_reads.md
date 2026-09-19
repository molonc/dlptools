# Perform joint segmentation of read state calls

This function takes a cell bin state calls and uses recursive binary
segmentation to find as many breakpoints as it can. By default, it looks
for a maximum number of breakpoints equal to half the number of input
bins. E.g., 100 bins of data, will search for at most 50 breakpoints.

Joint segmentation happens across all cells for a single chromosome at a
time. This function is a wrapper for all chromosomes, calling
[joint_seg_chromosome](https://molonc.github.io/dlptools/reference/joint_seg_chromosome.md)
for each to do the joint segmentation.

This function was intended for use with the already segmented profiles
of hmmcopy, by leveraging the state calls. And the goal is to modify
them as little as possible to correct for slight variances among cells
in where the segment breakpoints are inferred.

## Usage

``` r
joint_seg_reads(
  reads_df,
  chrom_col = "chr",
  max_bps = NULL,
  bin_fraction = 0.5
)
```

## Arguments

- reads_df:

  dataframe. Binned read data.

- chrom_col:

  string. Name of the chromosome column

- max_bps:

  NULL/int. Can specify the maximum number of breakpoints to find

- bin_fraction:

  float. What fraction of the bins to target as the maximum number of
  segments.

## Value

dataframe. Cell segments that match the jointly found segments.
