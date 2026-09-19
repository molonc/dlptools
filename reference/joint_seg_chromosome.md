# Perform joint segmentation on a single chromosome

Perform joint segmentation on a single chromosome

## Usage

``` r
joint_seg_chromosome(
  chr_df,
  chrom_col = "chr",
  max_bps = NULL,
  bin_fraction = 0.5
)
```

## Arguments

- chr_df:

  dataframe. State calls of the DLP bins data.

- chrom_col:

  string. Name of the chromosome column

- max_bps:

  int. Can specify the maximum number of breakpoints to find

- bin_fraction:

  float. What fraction of the bins to target as the maximum number of
  segments.

## Value

dataframe. Cell segments that match the jointly found segments.

## Details

see
[joint_seg_reads](https://molonc.github.io/dlptools/reference/joint_seg_reads.md)
for explanation. But effectively, joint segmentation needs to happen per
chromosome. So this function does the segmentation, and the other is
just a wrapper for all chromosomes.
