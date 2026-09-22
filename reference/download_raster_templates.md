# Download and unpack raster templates

Downloads a ZIP from Zenodo (or a user URL), unpacks it into a
destination directory, and removes the archive after a successful
extraction. Writes to a tempfile first and then moves to final
destination for atomicity.

## Usage

``` r
download_raster_templates(
  url = "https://zenodo.org/api/records/14497070/files-archive",
  out_dir = "./Templates/TemplateRasters",
  overwrite = FALSE,
  quiet = FALSE
)
```

## Source

Zenodo: https://doi.org/10.5281/zenodo.14497070

## Arguments

- url:

  Character. Archive URL to download. Default: Zenodo (EGV raster
  templates) v2.0.0:
  "https://zenodo.org/api/records/14497070/files-archive".

- out_dir:

  Character. Destination directory for the unzipped rasters. Default:
  "./Templates/TemplateRasters".

- overwrite:

  Logical. Re-download and overwrite existing files/dirs? Default:
  FALSE.

- quiet:

  Logical. Suppress progress messages. Default: FALSE.

## Value

Invisibly returns a named list with `out_dir`, `files_written`.

## See also

[`download_vector_templates()`](https://aavotins.github.io/egvtools/reference/download_vector_templates.md),
[`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
download_raster_templates()
} # }
```
