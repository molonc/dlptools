# generic function to import dlp files

Will create a dataframe of all files of a certain type. This is for the
three common dlp files that you might want to load: metrics, segments,
or reads.

## Usage

``` r
import_dlp_files(dlp_sc_dir, file_type = c("metrics", "segs", "reads"))
```

## Arguments

- dlp_sc_dir:

  the parent directory where the dlp outputs are saved (see the
  README.md for expected structure)

- file_type:

  options are: "metrics", "segs", "reads"

## Value

a tibble of the information from all dlp outputs combined
