# load UCSC gap files for telomeres

Honestly, this function is sort of pointless. The telomeres in both
genome versions are marked as 10Kb for p and q on all chromosomes.
Telomeres are variable and better estimates of their size exist, e.g.,
https://www.nature.com/articles/s41467-024-48917-7

## Usage

``` r
import_telos_file(version = c("hg19", "hg38"))
```

## Arguments

- version:

  string. hg19 (default) or hg38

## Value

tibble of telomere information

## Details

This function just serves as a way to set up easy marking of where CNs
occur on their respective chromosome.

The only exception to the 10Kb size is chr17. It is not listed in the
gap file for hg19. chr17 is known to have small chromosomes. Here, I set
the p to 3000 Kb and q to 5000 Kb, inferred from the article referenced
above.
