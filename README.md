# DLPTools

And R package for basic DLP+ data manipulation and plotting. Basically just a collection of the various functions we have been using to handle DLP data.

## Setup

A lot of the import functions are based on the idea that you have downloaded your dlp something like this:

```bash
cd where/to/save/my_dlp/

for ticket in "$@" ; 
do
  azcopy copy https://singlecellresults.blob.core.windows.net/results/${ticket}/results/annotation/ ${ticket} --recursive
  azcopy copy https://singlecellresults.blob.core.windows.net/results/${ticket}/results/hmmcopy/ ${ticket} --recursive
done
```

which produces a directory with this sort of structure:

```bash
├── SC-8382
│   ├── annotation
│   │   ├── metrics.csv.gz
│   │   ├── #[...other files ...]
│   └── hmmcopy
│   │   ├── #[...other files ...]
│       ├── reads.csv.gz
│       ├── reads.csv.gz.yaml
│       ├── reads_filtered.csv.gz
│       ├── segments.csv.gz
│       ├── segments.csv.gz.yaml
│       ├── segments_filtered.csv.gz
│       └── segs.tar.gz
├── SC-8408
│   ├── annotation
│   │   ├── metrics.csv.gz
│   │   ├── #[...other files ...]
│   └── hmmcopy
│   │   ├── #[...other files ...]
│       ├── reads.csv.gz
│       ├── reads.csv.gz.yaml
│       ├── reads_filtered.csv.gz
│       ├── segments.csv.gz
│       ├── segments.csv.gz.yaml
│       ├── segments_filtered.csv.gz
│       └── segs.tar.gz
├── SC-8650
```


As well, a lot of functions require external files that are too big to package effectively (TODO: maybe add downloads as a function?).

```bash
# centromere file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

# HG19 chromsome length information
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz

# information on where telomeres start and end
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
```

You'll also want a local copy of the masking file found in `meta_data/blacklist_2023.07.17.txt`, or your own version of regions to mask. This file was build by Daniel Lai, and has been used for many DLP analyses. It is masking rough approximations of centromeres and telomeres that are generally bad for DLP analyses. Other functions (*coming soon*) will add the explicit coordinates of centromeres and telomeres to dataframes of reads and segs.


## Main functions

File imports:

```R
library(dlptools)

# assuming you have the dlp files saved above, there are convenience functions 
# to load data:
my_dlp_dir <- "/projects/molonc/scratch/bfurman/projects/osteosarcoma/dlp/"
segs <- dlptools::import_dlp_files(my_dlp_dir, 'segs')
metrics <- dlptools::import_dlp_files(my_dlp_dir, 'metrics')
reads <- dlptools::import_dlp_files(my_dlp_dir, 'reads')
```


Adding features to reads/segs dataframes:

```R
# add a column called 'mask' to indicate that the reads or segs overlap with
# regions we want to exclude:
segs <- dlptools::mark_mask_regions(segs, "meta_data/blacklist_2023.07.17.txt")

reads <- dlptools::mark_mask_regions(reads, "meta_data/blacklist_2023.07.17.txt")
```


Plotting:

```R
# will create a simple plot with geom_tile of read states facetted by chromosome
dlptools::basic_tile_plot(reads)
```

## Development

```bash
# 1. clone the repo

git clone URL

# 2. install the pre-commit hooks
pre-commit install

# 3. code happy.
```


If you haven't worked on an R package before, [this](https://r-pkgs.org/whole-game.html) is a good read to get started.

To add functions, we can rely on `library(devtools)` to help with the metadata stuff and package construction. Steps include:

1. boot an R instance and set the working directory to wherever you cloned this repo

```R
setwd("/path/to/dlptools/")
library(devtools)
```

2. add your R function either within an existing `R/*R` module, or create a new module

```
vim R/my_new_module.R
```

```R
# add your function, with documentation in roxygen style! See other funcs for 
# examples.

#' the thing my function does
#'
# my function does this great thing
#' 
#' @param x a description of the arg for my function
#' @export
my_cool_func <- function(x) {
    return(x)
}

# also, try and use namespaces for functions that come from other packages, 
# e.g., purrr::map() instead of just map(), or dplyr::select() instead of just select()
# leads to fewer mixups as code bases grow and encounter redundantly named functions.
```

3. Then in your R instance, we can build the needed meta data for the new function.

```R
# searches the R/ directory for new functions and updates the NAMESPACE file.
# also creates .Rd files based on your function docstring.
devtools::document()
```

If your function includes code from another R package (even base R), you need to declare that and add it to the package metadata (unless it's already listed in the `Imports: ` section of the `DESCRIPTION` file).

From the R instance:

```R
usethis::use_package("purrr") # or whatever package your function needs
```


4. From the R instance, run a check that the package is still good

```R
devtools::check()
```

5. ideally we write unit tests 😬

https://r-pkgs.org/whole-game.html#use_testthat

6. commit your changes to the repo

```bash
git status # to see what is all changed

git add R/my_new_module.R
# and probably the other files that the package updated
git add NAMESPACE
git add man/my_cool_func.R
git commit
# and battle it out with the pre-commit hooks...the R ones are painfully slow

# write a GOOD commit message. For inspiration look at other commits and 
# checkout https://cbea.ms/git-commit/
# also try to include the type of thing being committed, i.e.,
# conventional commits style: https://www.conventionalcommits.org/en/v1.0.0/
```