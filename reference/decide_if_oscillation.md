# decide if a 3 segment chain counts as an oscillation

See
[extract_oscillations](https://molonc.github.io/dlptools/reference/extract_oscillations.md)
for a discussion of oscillations and recognizing them in sample copy
number segments.

## Usage

``` r
decide_if_oscillation(
  left_cn,
  middle_cn,
  right_cn,
  middle_bound = Inf,
  ends_bound = 0
)
```

## Arguments

- left_cn:

  int. CN value of segment on right side of the chain

- middle_cn:

  int. CN value of segment middle of the chain

- middle_bound:

  int\|Inf. Default Inf. How different the middle CN is allowed to be
  from either side to count as an oscillation. See
  [extract_oscillations](https://molonc.github.io/dlptools/reference/extract_oscillations.md)
  for a discussion of this

- ends_bound:

  int. Default 0. How different the ends of a 3-segment set are allowed
  to be to count as an oscillation

## Value

boolean. True if 3 CNs constitute an oscillation, False if not.
