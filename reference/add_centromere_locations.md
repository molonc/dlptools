# add centromere information to reads by chromosome

See dlptools::read_and_prep_ucsg_cenrtomeres() for details of file
origins.

## Usage

``` r
add_centromere_locations(cn_df, version = c("hg19", "hg38"))
```

## Arguments

- cn_df:

  dataframe of cn states for read bins or segments

- version:

  string. base centromeres on genome hg19 or hg38. Default: hg19

## Value

input table with centromere information added by chromosome
