# Create constant-background rasters from a directory of GeoTIFFs

For every raster in `in_dir`, creates a new raster where all **non-NA**
cells are set to `background_value` and NA cells are preserved. Outputs
go to `out_dir` (created if needed). Files are named with a prefix (by
default `nulls_` when `background_value == 0`, otherwise
`bg[background_value]_`) followed by the original name, e.g.
`bg10_my.tif`. Uses LZW compression and atomic writes. Large rasters can
stream to disk with `terra_todisk=TRUE`.

## Usage

``` r
create_backgrounds(
  in_dir,
  out_dir = NULL,
  background_value = 0,
  pattern = NULL,
  recursive = FALSE,
  out_prefix = NULL,
  NAflag = NULL,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  write_datatype = NULL,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = FALSE,
  force_gc = FALSE,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- in_dir:

  Directory containing input rasters.

- out_dir:

  Directory for outputs. If `NULL`, defaults to
  `file.path(in_dir,"backgrounds")`.

- background_value:

  Numeric value to assign to all non-NA pixels. Default `0`.

- pattern:

  Optional regex to filter filenames (applied after extension filter).

- recursive:

  Logical; search subdirectories in `in_dir`. Default `FALSE`.

- out_prefix:

  Optional filename prefix. If `NULL`, it defaults to:

  - `"nulls_"` when `background_value == 0`,

  - otherwise `"bg[background_value]_"` (e.g., `"bg10_"`, `"bg0.5_"`).

- NAflag:

  Optional NA flag for writing (passed to GDAL). Default `NULL` (auto).

- gdal_opts:

  Character vector of GDAL creation options (merged with tuned
  defaults). Default
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Optional terra datatype for writing (e.g., `"INT2S"`, `"FLT4S"`).
  Default `NULL` (auto as described above).

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra operations. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical or `NA`. If `TRUE`, prefer on-disk processing for the `ifel`
  step. Default `FALSE`.

- force_gc:

  Logical; call [`gc()`](https://rdrr.io/r/base/gc.html) at checkpoints.
  Default `FALSE`.

- overwrite:

  Overwrite existing outputs? Default `FALSE`.

- quiet:

  Suppress progress prints ([`cat()`](https://rdrr.io/r/base/cat.html))?
  Default `FALSE`.

## Value

Invisibly, a data.frame with columns: `in_file`, `out_file`, `n_cells`,
`elapsed_sec`.

## Details

- Reads each input with
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html).

- Applies a fast element-wise replacement with
  `terra::ifel(!is.na(x), background_value, NA)`.

- Writes with LZW compression to a temporary `._tmp.tif`, then
  atomically renames to the final file.

- If you do not set `write_datatype`/`NAflag`, the function picks
  sensible defaults:

  - If `background_value` is an integer and within `[-32768, 32767]`,
    then `INT2S` + `NAflag=-32768`.

  - Otherwise `FLT4S` + float `NAflag` (`-3.4028235e38`).

- Only `.tif` / `.tiff` files are processed (optionally refine with
  `pattern`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Classic "nulls_" backgrounds (fill 0)
create_backgrounds(
  in_dir = "./Templates/TemplateRasters",
  out_dir = "./Outputs/Nulls",
  background_value = 0,
  overwrite = TRUE,
  terra_todisk = TRUE
)

# Fill non-NA with 10, prefix becomes "bg10_"
create_backgrounds(
  in_dir = "./Templates/TemplateRasters",
  out_dir = "./Outputs/BG10",
  background_value = 10,
  overwrite = TRUE
)
} # }
```
