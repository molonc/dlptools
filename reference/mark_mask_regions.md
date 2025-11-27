# add a mask column to a dataframe

Determine if any start–end regions in the segs or reads dataframe
overlaps with regions in the mask file. Will add a 'mask' boolean column
to the dataframe to indicate whether the region should be masked.

## Usage

``` r
mark_mask_regions(
  subject_df,
  mask_f = NULL,
  mask_chr_name = "seqnames",
  subject_chr_name = "chr"
)
```

## Arguments

- subject_df:

  a table with chromosome, start, and end columns. Typically the reads
  or segs dataframe created by
  [`import_dlp_files()`](https://molonc.github.io/dlptools/reference/import_dlp_files.md).

- mask_f:

  path to a mask file, see README.md \## Setup for details.

- mask_chr_name:

  the name of the chromosomes column in the mask file. Default is
  "seqnames", which is the name in the provided file.

- subject_chr_name:

  name of the chromosome column in the reads/segs table. Default is
  'chr' which is what the import scripts call the column.
