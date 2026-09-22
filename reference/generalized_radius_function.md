# Extract and Rasterize Buffered-Radius Summary Statistics

Extracts one summary statistic per buffered polygon from one or more
raster layers, transfers those values to a rasterization grid using a
user-defined key field, rasterizes the results to a common template, and
merges the per-tile outputs into final rasters.

The relationship between buffered point grids, radii, join keys, and
rasterization grids is controlled by `buffer_mode` and, for fully custom
workflows, `buffer_spec`. This avoids hard-coding particular grid names
or key fields in the extraction logic.

## Usage

``` r
generalized_radius_function(
  kvadrati_path,
  radii_path,
  reference_grid_path = NULL,
  template_path,
  input_layers,
  layer_prefixes,
  output_dir = "./Extracted_Layers",
  unlink_tiles = TRUE,
  n_workers = 1,
  radii = c("r500", "r1250", "r3000", "r10000"),
  fill_missing = TRUE,
  radius_mode = "sparse",
  IDW_weight = 2,
  extract_fun = "mean",
  future_max_size = 8 * 1024^3,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  write_datatype = NULL,
  NAflag = NULL,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = TRUE,
  quiet = FALSE,
  buffer_mode = NULL,
  buffer_spec = NULL,
  tikls100_path = NULL,
  crop_buffer = 1000
)
```

## Arguments

- kvadrati_path:

  Character. Directory containing tiled GeoParquet files for the base
  rasterization grid. The tile identifier is taken from the final
  underscore-separated component of each file name, before the file
  extension. For example, `grid_2434.parquet` has tile ID `2434`.

- radii_path:

  Character. Directory containing tiled buffered-polygon GeoParquet
  files. Fixed-radius files are expected to follow the naming convention
  `<buffer_grid>_r<radius>_<tileid>.parquet`, matching the output
  convention of
  [`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md).

- reference_grid_path:

  Character or `NULL`. Optional path to a GeoParquet rasterization grid
  used by the built-in `"sparse"` mode for coarser buffer grids. The
  grid must contain the relevant `key_field` values. It is not required
  for `"dense"` mode or when a custom `buffer_spec` supplies all
  required grid paths.

- template_path:

  Character. Path to the template raster defining final extent,
  resolution, alignment, mask, and CRS.

- input_layers:

  Named or unnamed character vector of raster paths from which summary
  statistics are extracted. Names are replaced by `layer_prefixes`.

- layer_prefixes:

  Character vector with one output prefix for each element of
  `input_layers`. Values must be unique.

- output_dir:

  Character. Output directory root. Default `"./Extracted_Layers"`.

- unlink_tiles:

  Logical. If `TRUE`, delete intermediate per-tile raster directories
  after successful mosaicking. Default `TRUE`.

- n_workers:

  Integer. Number of parallel workers. On SLURM systems this is capped
  at the detected CPU allocation. Default `1`.

- radii:

  Character or numeric vector selecting radii to process. Numeric values
  are interpreted as distances and converted to IDs such as `500` to
  `"r500"`. Character values may be supplied either as `"500"` or
  `"r500"`. In `"sparse"` and `"dense"` modes the default processes
  `r500`, `r1250`, `r3000`, and `r10000`. In `"specified"` mode, if
  `radii` is omitted, all rows of `buffer_spec` are processed.

- fill_missing:

  Logical. If `TRUE`, fill NA gaps within the template footprint using
  [`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html).
  Default `TRUE`.

- radius_mode:

  Deprecated compatibility argument corresponding to the former
  `radius_mode`. Prefer `buffer_mode`. If `buffer_mode` is `NULL`, this
  value is used. Default `"sparse"`.

- IDW_weight:

  Numeric. Inverse-distance weighting power passed to
  [`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html).
  Default `2`.

- extract_fun:

  A function or a single character string such as `"mean"`, passed to
  [`exactextractr::exact_extract()`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html).
  It must resolve to exactly one numeric value per polygon and raster
  layer.

- future_max_size:

  Numeric. Maximum size, in bytes, of globals exported to each future
  worker. Default `8 * 1024^3`.

- gdal_opts:

  Character vector of GDAL creation options used for raster writes. By
  default LZW compression, tiling, BigTIFF-if-needed, and 256-cell
  blocks are used. Threading is reduced to one GDAL thread per worker
  when parallel processing is active.

- write_datatype:

  Character or `NULL`. Optional terra datatype such as `"FLT4S"` or
  `"INT2S"`. Default `NULL` uses terra's default.

- NAflag:

  Numeric or `NULL`. Optional NA flag used when writing rasters. Default
  `NULL`.

- terra_memfrac:

  Numeric. Fraction passed to `terra::terraOptions(memfrac = ...)`.
  Default `0.7`.

- terra_tempdir:

  Character. Temporary directory used by terra and for intermediate
  mosaics. Default [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical. If `TRUE`, prefer on-disk terra operations. Default `TRUE`.

- quiet:

  Logical. If `TRUE`, suppress progress messages. Default `FALSE`.

- buffer_mode:

  Character or `NULL`. Preferred replacement for `radius_mode`. One of
  `"sparse"`, `"dense"`, or `"specified"`.

  - `"sparse"` reproduces the package's standard sparse setup: `pts100`
    at 500 and 1250 m, `pts300` at 3000 m, and `pts1000` at 10000 m. The
    500 and 1250 m results join directly to each tile by `id`; the 3000
    and 10000 m results use `reference_grid_path` with `rinda300` and
    `ID1km`, respectively.

  - `"dense"` uses `pts100` for all four standard radii and joins every
    result directly to the current tile by `id`.

  - `"specified"` uses `buffer_spec` and contains no built-in
    assumptions about buffer-grid names, radii, or key fields.

  If `NULL`, the value of the compatibility argument `radius_mode` is
  used.

- buffer_spec:

  A data.frame used when `buffer_mode = "specified"`. It must contain
  exactly one processing rule per radius and the columns:

  - `buffer_grid`: character prefix identifying buffered files. For
    example, `"pts250"` matches files such as
    `pts250_r2000_2434.parquet`.

  - `radius_m`: positive whole-number buffer radius used to construct
    the radius ID `r<radius_m>` and locate buffered files.

  - `key_field`: character name of the field that links each buffered
    polygon to the grid on which its extracted value is rasterized.

  - `grid_path`: character path to a GeoParquet rasterization grid, or
    `NA`/`"tile"` to use the current tile from `kvadrati_path`.

  Each radius must occur only once because each layer-radius combination
  produces one output raster.

- tikls100_path:

  Deprecated compatibility alias for `reference_grid_path`. If supplied,
  `reference_grid_path` must be `NULL`.

- crop_buffer:

  Numeric. Additional map-unit buffer around the bounding box of the
  largest available buffered polygons when cropping input rasters and
  the template. Default `1000`.

## Value

Invisibly returns a data.frame with one row per successfully written
layer-radius combination and columns:

- `layer`: output layer prefix;

- `radius`: radius ID such as `"r3000"`;

- `radius_m`: numeric radius;

- `buffer_grid`: buffered-file grid prefix;

- `output_path`: final raster path;

- `n_tiles_merged`: number of tile rasters included in the mosaic;

- `gaps_before`: number of NA cells inside the template before filling;

- `gaps_after`: number remaining after optional filling;

- `filter_size_used`: Whitebox filter width, or `NA` if not used; and

- `gap_filled`: logical indicating whether gap filling was successfully
  applied.

## Details

### Buffer specification

Internally all three modes are converted to the same four-column
specification: `buffer_grid`, `radius_m`, `key_field`, and `grid_path`.
Consequently, the extraction and rasterization code does not contain
radius-specific rules.

If `grid_path` is `NA`, empty, or `"tile"`, extracted values are joined
to the current tiled grid from `kvadrati_path`. Otherwise the referenced
GeoParquet file is loaded once per worker, filtered to the key values
present in the current buffer polygons, and used as the rasterization
grid.

The `key_field` must exist in both the buffered polygons and the
selected rasterization grid. It must uniquely identify buffered polygons
within a buffer file. It may be repeated in the rasterization grid; this
is what permits a statistic calculated for a coarser buffer center to be
assigned to multiple finer output cells.

### Workflow

1.  Resolve `buffer_mode` to a validated buffer specification.

2.  Discover base-grid tiles and buffered files using their file-name
    conventions.

3.  For each tile, read the available buffered polygons and identify the
    largest available radius.

4.  Crop the template and raster stack to the largest buffered extent
    plus `crop_buffer`.

5.  Run
    [`exactextractr::exact_extract()`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
    once per radius over the cropped raster stack.

6.  Join each extracted statistic to the rule's rasterization grid using
    `key_field`, rasterize with
    [`fasterize::fasterize()`](https://rdrr.io/pkg/fasterize/man/fasterize.html),
    and write per-tile GeoTIFFs.

7.  VRT-mosaic the per-tile rasters, project and mask them to
    `template_path`, optionally fill internal NA gaps, and write the
    final GeoTIFF.

### Parallel and memory behavior

Multisession futures are used rather than forking. Large reusable
objects, including the template, input raster stack, and external
rasterization grids, are cached once per worker through the internal
egvtools worker cache. Cache keys include the normalized source paths
(and raster layer names where relevant), preventing consecutive calls
with different inputs from reusing stale file-backed raster objects.
BLAS, OpenMP, and GDAL thread counts are constrained inside workers to
reduce oversubscription on HPC systems.

### Gap filling

After mosaicking and alignment, gaps are defined as NA cells inside the
non-NA template footprint. When `fill_missing = TRUE`, the maximum
distance into the gap mask is used to derive an odd Whitebox filter
width with a minimum of three cells. Gap counts before and after filling
are returned.

## See also

[`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard sparse workflow.
generalized_radius_function(
  kvadrati_path = "./Templates/TemplateGrids/lapas/",
  radii_path = "./Templates/TemplateGridPoints/lapas/",
  reference_grid_path =
    "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
  template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
  input_layers = c(
    Soils_txtSand = "./Rastri_100m/RAW/Soils_txtSand_cell.tif",
    Soils_txtSilt = "./Rastri_100m/RAW/Soils_txtSilt_cell.tif"
  ),
  layer_prefixes = c("Soils_txtSand", "Soils_txtSilt"),
  buffer_mode = "sparse",
  n_workers = 4
)

# Fully custom relationship between buffered files and rasterization grids.
spec <- data.frame(
  buffer_grid = c("pts100", "pts100", "pts300", "pts1000"),
  radius_m = c(500, 1250, 3000, 10000),
  key_field = c("id", "id", "rinda300", "ID1km"),
  grid_path = c(
    NA,
    NA,
    "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
    "./Templates/TemplateGrids/tikls100_sauzeme.parquet"
  ),
  stringsAsFactors = FALSE
)

generalized_radius_function(
  kvadrati_path = "./Templates/TemplateGrids/lapas/",
  radii_path = "./Templates/TemplateGridPoints/lapas/",
  template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
  input_layers = c(NDVI = "./Rastri_100m/NDVI_cell.tif"),
  layer_prefixes = "NDVI",
  buffer_mode = "specified",
  buffer_spec = spec,
  n_workers = 4
)

# A custom grid can join directly back to the current tile.
custom_spec <- data.frame(
  buffer_grid = "pts200",
  radius_m = 750,
  key_field = "id",
  grid_path = NA_character_,
  stringsAsFactors = FALSE
)
} # }
```
