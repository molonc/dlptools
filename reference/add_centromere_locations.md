# add centromere information to reads by chromosome

See dlptools::read_and_prep_ucsg_cenrtomeres() for details of file
origins.

## Usage

``` r
add_centromere_locations(
  cn_df,
  centro_file = NULL,
  version = c("hg19", "hg38")
)
```

## Arguments

- cn_df:

  dataframe of cn states for read bins or segments

- centro_file:

  NULL or string to the path if you download yourself

- hg19:

  boolean to target hg19 for loading

- hg38:

  boolean to target hg38 for loading

## Value

input table with centromere information added by chromosome
