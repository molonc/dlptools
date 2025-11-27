# create a list of intervals spanning chromosome arms

Splits a chromosome at the middle of the centromere. Sets up intervals
for splitting each chromosome arm.

## Usage

``` r
create_chrom_arm_intervals(genome_version = c("hg19", "hg38"))
```

## Arguments

- genome_version:

  string. "hg19" (default) or "hg38"

## Value

list. Named by chromosome, vectors of how to break a chromsome into
intervals of arms.
