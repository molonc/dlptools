# generate a ComplexHeatmap::Heatmap image, either to console or file.

Will return an image of the generated heatmap or dump the heatmap to a
file.

## Usage

``` r
output_hm_image(
  hm,
  file_name = NULL,
  png_height = 1600,
  png_width = 2800,
  png_res = 144,
  custom_legend = list(),
  hm_title
)
```

## Arguments

- file_name:

  optional string of where to save a png image of the heatmap.

- png_height:

  optional height of png file

- total_hm:

  ComplexHeatmap::Heatmap of the combined tree and states.

## Value

ComplexHeatmap::draw or nothing if a file is written.
