# extract the legths of chains of osscilating copy number segments

Osscilating segments of copy numbers take from of N-M-N, e.g., 3-4-3.
The length of these chains reflects the number of oscillations. Chains
are evaluated in 3-segment intervals (the minimum required to recognize
an oscillation). Chromosomse with \< 3 segments cannot have
oscillations, and will receive a 0.

## Usage

``` r
extract_oscillations(
  segs_df,
  middle_bound = 2,
  ends_bound = 0,
  sample_col = "cell_id",
  chrom_col = "chr",
  cn_col = "state",
  return = c("values", "counts")
)
```

## Arguments

- segs_df:

  dataframe. Copy number segments for samples.

- middle_bound:

  integer. How many copy number values away the middle segment can be
  from the end segments to count as an oscillation. Default 2 (e.g.,
  3-4-3 and 3-5-3 both count as oscillations, but 3-7-3 would not).

- ends_bound:

  integer. How many CN values apart the two ends of a 3-segment set are
  allow to be to count as an oscillation. Default is 0, i.e., they must
  be the same copy number.

- sample_col:

  string. Name of column with cell/sample names

- chrom_col:

  string. Name of column with chromosome names

- cn_col:

  string. Name of column with segment copy number states.

- return:

  string. "values" (default) or "counts". Values are the observed values
  for cells, counts are the counts of these values in pre-determined
  categories.

## Value

dataframe. sample IDs and observed chain lengths.

## Details

Chromosomes are evaluated in 3-segments sliding windows, with a window
receiving a 1 if the 3-set is an oscillation, and a 0 if not. Adjacent
oscillation sets are summed, non-osscilating sets are left as 0.

Examples of Chain : length

- 2-3-4: 0

- 3-4-3 : 1

- 3-4-3-4: 2

- 3-4-3-2: 1 0 (3-4-3 is a chain, 4-3-2 is not)

- 2-3-4-5: 0 0

An early paper arguing for the use of oscillations to detect
chromothripsis outlined that oscillations should be within 1 or 2 CN
values. Using N-M-N notation from above, this would mean that Ns should
be the same value and M would be within 1 or 2 copy values of N. Some
pubilications relax these constraints, and for example allow M to be any
value other than N. This means that 3-300-3 is counted as an oscillation
by several publications. Granted, this is an oscillation, but not likely
to be one indicative of chromothripsis.

Others allow some flexiblity in the matching of the first and second N,
allowing them to be within some margin of difference. Componding these
two things leads to fairly extreme senarios, with some methods allowing
chains like 5-300-6 to count as an oscillation.

These behaviours can be replicated here, but isn't recommended. By
default, Ns must match exactly, and M must be within 2 CNs.

Optionally can summarize in counts of pre-defined categories of:

Oscillation chains of length: 0, 1-3, 4-9, 10+
