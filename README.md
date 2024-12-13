# DLPTools

And R package for basic DLP+ data manipulation and plotting. Basically just a collection of the various functions we have been using to handle DLP data.

## Installing DLPtools

```r
library(devtools)
devtools::install_github("molonc/dlptools")
```

## Vignettes

See `vignettes/` for various overview (Rmd & htmls available) or browse the articles at the [package website](https://molonc.github.io/dlptools/).

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

3. To test your function, in the R instance you can load it and check it out

```R
# will load the whole package including your new function
devtools::load_all()

my_cool_func(x)
```

4. Then in your R instance, we can build the needed meta data for the new function.

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


5. From the R instance, run a check that the package is still good

```R
devtools::check()
```

6. ideally we write unit tests 😬

https://r-pkgs.org/whole-game.html#use_testthat

7. commit your changes to the repo

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