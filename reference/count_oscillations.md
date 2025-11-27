# count the length of oscillating chains of CN.

Given a vector of CN values, detect and summarize lengths of any
oscillating chains that are present in the sequence, along with 0s for
non-oscillating sets.

## Usage

``` r
count_oscillations(cn_vals, middle_bound = 2, ends_bound = 0)
```

## Arguments

- cn_vals:

  vector of ints. Segment state values. Ideally from the chromosome of 1
  sample.

- middle_bound:

  int\|Inf. Default 2. How different the middle CN is allowed to be from
  either side to count as an oscillation. See
  [extract_oscillations](https://molonc.github.io/dlptools/reference/extract_oscillations.md)
  for a discussion of this. Default of 2 means the middle CN must be
  within 2 copy values of the flanking segments.

- ends_bound:

  int. Default 0. How different the ends of a 3-segment set are allowed
  to be to count as an oscillation

## Value

vector of ints.

## Details

A full discussion of this function can be found in the documentation of.
[extract_oscillations](https://molonc.github.io/dlptools/reference/extract_oscillations.md).

Examples:

- 3-2-3 = 1

- 1-2-3 = 0

- 3-2-3-1-2-3 = 1 0

- 3-2-3-2 = 2

- 3-2-3-2-1 = 2 0
