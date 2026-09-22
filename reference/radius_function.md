# Extract and Rasterize Summary Statistics from Buffered Radii Using exactextractr

Extracts summary statistics from raster layers using buffered polygon
zones of multiple radii and rasterizes them onto a common template grid.
Designed for parallel execution.

## Usage

``` r
radius_function(
  kvadrati_path,
  radii_path,
  tikls100_path,
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
  quiet = FALSE
)
```

## Arguments

- kvadrati_path:

  Path to directory containing 1 ha grid-cell GeoParquet tiles (e.g.,
  `.../lapas/`).

- radii_path:

  Path to directory with buffered polygon radii tiles (`r500`, `r1250`,
  `r3000`, `r10000`).

- tikls100_path:

  Path to GeoParquet for the 100 m grid (must contain `rinda300` and
  `ID1km`).

- template_path:

  Path to the template raster (final alignment target; e.g.,
  `LV100m_10km.tif`).

- input_layers:

  Named character vector of covariate raster paths (values are file
  paths, names are the **layer prefixes** used in outputs).

- layer_prefixes:

  Character vector of names (same length/order as `input_layers`) to use
  as output layer prefixes.

- output_dir:

  Output directory root. Default `"./Extracted_Layers"`.

- unlink_tiles:

  Logical; delete per-tile rasters after merging. Default `TRUE`.

- n_workers:

  Number of parallel workers. Default `1`.

- radii:

  Character vector of radii to process (subset of: `"r500"`, `"r1250"`,
  `"r3000"`, `"r10000"`).

- fill_missing:

  If `TRUE`, run Whitebox IDW gap filling on the projected mosaic.

- radius_mode:

  `"sparse"` or `"dense"`. Controls which buffered files are used:
  *`"sparse"`* uses `pts100` for r500/r1250, `pts300` for r3000, and
  `pts1000/pts1km` for r10000; *`"dense"`* uses only `pts100` for all
  radii.

- IDW_weight:

  Numeric power for Whitebox IDW. Default `2`.

- extract_fun:

  Function or single string (e.g., `"mean"`) passed to
  [`exactextractr::exact_extract()`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
  that returns **one scalar per polygon**. If it returns multiple values
  per polygon, the function stops.

- future_max_size:

  Max globals per worker (bytes). Default `8 * 1024^3`.

- gdal_opts:

  GDAL creation options for writes. Default
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Optional terra datatype (e.g., `"FLT4S"`, `"INT2S"`). Default `NULL`
  (terra default).

- NAflag:

  Optional NA flag for writing. Default `NULL` (terra default).

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for heavy ops. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  If `TRUE`, prefer on-disk operations. Default `TRUE`.

- quiet:

  Suppress progress prints. Default `FALSE`.

## Value

Invisibly, a data.frame with columns:
`layer, radius, output_path, n_tiles_merged, gaps_before, gaps_after, filter_size_used, gap_filled`.

## Details

**Workflow**

1.  Discover inputs (1 ha tiles and buffered radii files) according to
    `radius_mode`.

2.  **Per tile (parallel):** read tile zones and available radii; crop
    to bbox + 1 km; **run
    [`exactextractr::exact_extract()`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
    once per radius** over the (cropped) raster stack, producing one
    value per polygon per layer; join to the appropriate 100 m grid (or
    its parent for r3000/r10000), and rasterize to the cropped template
    via `fasterize`, writing **per-tile GeoTIFFs**.

3.  For each `layer at radius`: VRT-merge tiles then **project to the
    full template grid** then mask to the template footprint.

4.  **Gap analysis & optional fill:** count NA gaps inside the template.
    If `fill_missing = TRUE` and gaps exist, estimate max gap width on
    the template
    ([`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
    on a fillable mask), set Whitebox
    `filter = 2 * ceil(max_gap / pixel_size)` with **min clamp = 3**,
    run
    [`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html),
    then mask again.

5.  **CRS guard & write:** set
    `terra::crs(out) <- terra::crs(template, proj = FALSE)` and write
    with LZW tiling via **atomic writes**.

6.  **Return:** a data.frame with
    `layer, radius, output_path, n_tiles_merged, gaps_before, gaps_after, filter_size_used, gap_filled`.

**Output parsing guard** Handles single vectors, data.frames
(`"mean.Layer"` or layer names), and lists (one element per feature:
numeric vectors or 1-row data.frames). Clear error otherwise.

**Per-worker cache** Large `sf` objects like a ~6.5M-feature `tikls100`
are **loaded once per worker** to avoid repeated I/O and pointer
invalidation issues in parallel sessions.

**CRS guard** Outputs copy the CRS string from the template using
`terra::crs(template, proj = FALSE)` to keep a Proj.4-style string for
maximal compatibility.

## Examples

``` r
if (FALSE) { # \dontrun{
radius_function(
  kvadrati_path  = "./Templates/TemplateGrids/lapas/",
  radii_path     = "./Templates/TemplateGridPoints/lapas/",
  tikls100_path  = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
  template_path  = "./Templates/TemplateRasters/LV100m_10km.tif",
  input_layers   = c(Soils_txtSand  = "./Rastri_100m/RAW/Soils_txtSand_cell.tif",
                     Soils_txtSilt  = "./Rastri_100m/RAW/Soils_txtSilt_cell.tif",
                     Soils_txtClay  = "./Rastri_100m/RAW/Soils_txtClay_cell.tif"),
  layer_prefixes = c("Soils_txtSand","Soils_txtSilt","Soils_txtClay"),
  output_dir     = "./Extracted_Layers",
  n_workers      = 4,
  radii          = c("r500","r1250","r3000","r10000"),
  radius_mode    = "sparse",
  extract_fun    = "mean",
  fill_missing   = TRUE,
  IDW_weight     = 2
)
} # }
```
