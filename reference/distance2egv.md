# Distance-to-class on the EGV grid

Computes Euclidean distance (in **map units**) from cells matching a set
of class values in an **input raster** to all cells of an **EGV
template** grid, then writes a Float32 GeoTIFF aligned to the template.
Designed to work with rasters produced by
[`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md).

## Usage

``` r
distance2egv(
  input,
  template_egv,
  values_as_one = NULL,
  project_to_template_input = FALSE,
  template_input = NULL,
  use_whitebox = TRUE,
  fill_gaps = FALSE,
  idw_weight = 2,
  filter_size = NULL,
  outlocation = "./Rastri_100m/RAW/",
  outfilename,
  layername,
  check_na = FALSE,
  plot_result = FALSE,
  plot_gaps = FALSE,
  NAflag = NULL,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  write_datatype = NULL,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = FALSE,
  force_gc = FALSE,
  quiet = FALSE
)
```

## Arguments

- input:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or file path. Raster prepared by
  [`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md).

- template_egv:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or path. Target EGV template (e.g., 100 m).

- values_as_one:

  `NULL`, numeric vector, and/or range strings in bracket notation (see
  Details). If `NULL`, non-`NA` cells are sources.

- project_to_template_input:

  Logical. If `TRUE`, reproject `input` to `template_input` **once at
  the start**. Default `FALSE`.

- template_input:

  When `project_to_template_input=TRUE`, a 10 m template used to define
  the target CRS/grid for that initial reprojection.

- use_whitebox:

  Logical; use Whitebox for distance (default `TRUE`), otherwise
  [`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html).

- fill_gaps:

  Logical. If `TRUE`, fill remaining gaps on the **template grid** using
  Whitebox `wbt_fill_missing_data()`. Default `FALSE`.

- idw_weight:

  IDW power for Whitebox fill. Default `2`.

- filter_size:

  Optional odd integer window for Whitebox. If `NULL`, auto:
  `2 * max_gap_width_cells`, with a **minimum of 3**.

- outlocation:

  Output directory. Default `"./Rastri_100m/RAW/"`.

- outfilename:

  Output filename (e.g., `"dist_class_100m.tif"`). **Required**.

- layername:

  Output band name. **Required**.

- check_na:

  Logical. If `TRUE`, report internal NA count on the template
  footprint. Default `FALSE`.

- plot_result:

  Logical; plot final result.

- plot_gaps:

  Logical; plot gap map (TRUE = gap). If both TRUE, side-by-side.

- NAflag:

  Optional NA flag for writing. Default `NULL` (omitted for Float32).

- gdal_opts:

  GDAL creation options (merged with tuned defaults). Default
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Terra datatype for writing. Default `NULL` will be coded as `"FLT4S"`.

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra ops. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical or `NA`. If `TRUE`, prefer on-disk for heavy steps. Default
  `FALSE`.

- force_gc:

  Logical; call [`gc()`](https://rdrr.io/r/base/gc.html) at checkpoints.
  Default `FALSE`.

- quiet:

  Suppress console prints. Default `FALSE`.

## Value

Invisibly, a list with:

- `path` (output file),

- `n_sources` (number of source cells on input grid),

- `n_na_final` (internal NA count on template footprint),

- `min_dist`, `max_dist`,

- `elapsed_sec`.

## Details

**Preparation & CRS:** Prepare inputs with the function
[`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md)
and the templates from https://zenodo.org/records/14497070 (download via
the function
[`download_raster_templates()`](https://aavotins.github.io/egvtools/reference/download_raster_templates.md)).
The function does **not** automatically reproject to the 100 m template.
It only warns if CRSs differ. For rare cases, set
`project_to_template_input = TRUE` and provide `template_input` (10 m)
to reproject the input once at the start. Distance is computed in the
(possibly reprojected) inputs map units.

**Selecting class values:** `values_as_one` accepts any-length vector
combining:

- numeric values (e.g., `500`, `c(610,620,630)`), and/or

- range strings in interval notation with inclusive/exclusive bounds:
  `"[400,600]"`, `"(400,600)"`, `"[400,600)"`, `"(400,600]"`. Example:
  `c("[600,700)", "500")`. If `values_as_one` is `NULL`, **non-NA**
  cells are considered sources.

**Distance engine:** Set `use_whitebox = TRUE` (default) to compute the
distance with
[`whitebox::wbt_euclidean_distance()`](https://rdrr.io/pkg/whitebox/man/wbt_euclidean_distance.html);
otherwise
[`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
is used.

**Auto alignment choice:** If input distance grid and the 100 m template
are same-CRS, bounding boxes match, and the template resolution is an
**integer multiple** of the input resolution (e.g., 10 m to 100 m), the
function uses `terra::aggregate(..., fun=mean)` (fast), followed by a
light `resample(..., "near")` to lock onto the template grid. Otherwise,
it falls back to `terra::resample(..., method="mean")`.

**Masking:** A single mask to `template_egv` is applied once **after**
alignment and **before** plotting/saving.

**Gap filling:** If `fill_gaps = TRUE`, gaps (cells where output is `NA`
but the template is not) are filled via
[`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html)
with IDW weight `idw_weight`. If `filter_size` is `NULL`, the function
uses `2 * max_gap_width_cells`, with a **minimum of 3** (no max clamp).
Whitebox temporaries are written with `COMPRESS=NONE` for speed; the
final GeoTIFF uses `gdal_opts`.

**Workflow**

1.  Load `input` and `template_egv`; optionally reproject `input` once
    if `project_to_template_input=TRUE` (using `template_input`).

2.  Build a **seeds** raster on the input grid: `1` where value matches
    `values_as_one` (or non-`NA` if `values_as_one=NULL`), `NA`
    elsewhere.

3.  Compute **distance** on the input grid via Whitebox
    (`wbt_euclidean_distance`) or
    [`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
    depending on `use_whitebox`.

4.  **Align** the distance raster to the template:

    - If perfectly nested (same CRS, matching extent, integer resolution
      ratio), do `aggregate(mean)` then `resample("near")`.

    - Else if same CRS but not nested, `resample(method="mean")`.

    - Else, `project(..., method="bilinear")`.

5.  Apply a **single final mask** to `template_egv`.

6.  Optionally **fill gaps** with Whitebox (if `fill_gaps=TRUE`) on the
    template grid.

7.  Optionally **plot** result and/or gap map (side-by-side if both
    requested).

8.  **Write atomically** to GeoTIFF with LZW compression and tuned GDAL
    options.

9.  Restore `terraOptions()` and sink state on exit (prevents stuck
    sinks).

Console safety: uses [`cat()`](https://rdrr.io/r/base/cat.html) for
progress and snapshots/restores sink state on exit so your console will
not remain "sunk" after interrupts.

## See also

[`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md)
for preparing the input raster,
[`input2egv()`](https://aavotins.github.io/egvtools/reference/input2egv.md)
for aggregation/resampling to the EGV grid,
[`downscale2egv()`](https://aavotins.github.io/egvtools/reference/downscale2egv.md),
[`landscape_function()`](https://aavotins.github.io/egvtools/reference/landscape_function.md),
[`radius_function()`](https://aavotins.github.io/egvtools/reference/radius_function.md),
and the template downloaders:
[`download_raster_templates()`](https://aavotins.github.io/egvtools/reference/download_raster_templates.md),
[`download_vector_templates()`](https://aavotins.github.io/egvtools/reference/download_vector_templates.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Distance from classes 1 and [10,20) using Whitebox, with side-by-side plots
distance2egv(
  input         = "./TestejuPakotni_early/Forests_StandAge_bg0_b.tif",
  template_egv  = "./Templates/TemplateRasters/LV100m_10km.tif",
  values_as_one = c(1, "[10,20)"),
  outlocation   = "./Rastri_100m/RAW/",
  outfilename   = "dist_age_100m.tif",
  layername     = "dist_age",
  use_whitebox  = TRUE,
  plot_result   = TRUE,
  plot_gaps     = TRUE,
  terra_todisk  = TRUE
)

# Same with terra::distance() and gap filling
distance2egv(
  input         = "./Masks/roads_mask_10m.tif",
  template_egv  = "./Templates/TemplateRasters/LV100m_10km.tif",
  values_as_one = NULL,  # non-NA are sources
  use_whitebox  = FALSE,
  fill_gaps     = TRUE,
  idw_weight    = 2,
  outlocation   = "./Rastri_100m/RAW/",
  outfilename   = "dist_roads_100m.tif",
  layername     = "dist_roads",
  terra_todisk  = TRUE
)
} # }
```
