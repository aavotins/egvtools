# Download and unpack vector templates (points/grids/gpkg)

Downloads vector templates archive (by default from Zenodo) and places
files into:

- Parquet grids into `grid_dir`

- Parquet points into `points_dir`

- GPKG into `gpkg_dir` The function auto-classifies files by filename
  patterns: "tikls\*.parquet" (grids), "pts\*.parquet" (points), and
  "vector_grids.gpkg" (GPKG).

## Usage

``` r
download_vector_templates(
  url = "https://zenodo.org/api/records/14277114/files-archive",
  grid_dir = "./Templates/TemplateGrids",
  points_dir = "./Templates/TemplateGridPoints",
  gpkg_dir = "./Templates",
  overwrite = FALSE,
  quiet = FALSE
)
```

## Source

Zenodo: https://zenodo.org/records/14277114

## Arguments

- url:

  Character. Default:
  "https://zenodo.org/api/records/14277114/files-archive"

- grid_dir:

  Character. Default "./Templates/TemplateGrids"

- points_dir:

  Character. Default "./Templates/TemplateGridPoints"

- gpkg_dir:

  Character. Default "./Templates"

- overwrite:

  Logical. Overwrite existing files? Default FALSE.

- quiet:

  Logical. Suppress progress messages? Default FALSE.

## Value

Invisibly returns a list of the three dirs.

## See also

[`download_raster_templates()`](https://aavotins.github.io/egvtools/reference/download_raster_templates.md),
[`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
download_vector_templates()
} # }
```
