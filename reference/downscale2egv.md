# Downscale & Align a Raster to a Template, with optional gap fill & IDW smoothing

Aligns a source raster to a template grid (CRS, resolution, extent),
masks to the template footprint, and optionally: (1) fills NoData gaps
using WhiteboxTools' IDW-based `fill_missing_data`, and (2) applies
**IDW smoothing** to reduce blockiness from low-resolution inputs.

**Mass preservation:** IDW smoothing (and most simple smoothers) are
**not mass-preserving**. If conservation of totals/means matters
(globally or by zones), post-smoothing rescaling is recommended (see
*Details*).

If `interpolation_method = "auto"`, integer/factor inputs use `"near"`;
continuous inputs use `"bilinear"`.

## Usage

``` r
downscale2egv(
  template_path,
  grid_path,
  rawfile_path,
  out_path,
  file_name,
  layer_name,
  interpolation_method = "auto",
  buffer_m = 10000,
  check_na = FALSE,
  fill_gaps = TRUE,
  idw_weight = 2,
  filter_size_cells = "auto",
  plot_gaps = FALSE,
  plot_result = FALSE,
  smooth = FALSE,
  smooth_radius_km = 10,
  smooth_agg_factor = 10,
  smooth_power = 0.5,
  smooth_epsilon = 0,
  smooth_nmax = 50,
  smooth_force = FALSE,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
  write_datatype = NULL,
  NAflag = NULL,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = TRUE,
  force_gc = FALSE,
  return_visible = FALSE,
  quiet = FALSE
)
```

## Arguments

- template_path:

  Character path to a template GeoTIFF **or** a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- grid_path:

  Character path to a GeoParquet grid **or** an `sf` object (used only
  for bbox).

- rawfile_path:

  Character path to the input raster **or** a `SpatRaster`.

- out_path:

  Character. Directory for the output GeoTIFF (created if missing).

- file_name:

  Character. Output file name (e.g., `"result.tif"`).

- layer_name:

  Character. Band name to set on the output.

- interpolation_method:

  `"auto"` (default), or `"bilinear"`, `"near"`, `"cubicspline"`,
  `"cubic"`.

- buffer_m:

  Numeric. Buffer distance (meters) applied to the grid bbox before
  cropping. Default `10000`.

- check_na:

  Logical. If `TRUE`, report NA gap count (inside template). Default
  `FALSE`.

- fill_gaps:

  Logical. If `TRUE`, fill NA gaps via Whitebox. Default `TRUE`.

- idw_weight:

  Numeric. IDW power used by Whitebox gap filling. Default `2`.

- filter_size_cells:

  Integer or `"auto"`. Window size (in **cells**) for gap filling.
  Default `"auto"`.

- plot_gaps, plot_result:

  Logical flags for diagnostics/plots (side-by-side if both).

- smooth:

  Logical. Enable **IDW** smoothing. Default `FALSE`.

- smooth_radius_km:

  Numeric. IDW smoothing radius (km). Default `10`.

- smooth_agg_factor:

  Integer. Aggregate factor before point creation. Default `10`.

- smooth_power:

  Numeric. IDW power. Default `0.5`.

- smooth_epsilon:

  Numeric. Small smoothing term (if supported by your `terra`). Default
  `0`.

- smooth_nmax:

  Integer. Max neighbors (if supported by your `terra`). Default `50`.

- smooth_force:

  Logical. Allow smoothing for integer/factor rasters (not recommended).
  Default `FALSE`.

- gdal_opts:

  Character vector. GDAL creation options; merged with tuned defaults:
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Character or `NULL`. Optional terra datatype (e.g., `"FLT4S"`).

- NAflag:

  Numeric/integer NoData value to write; `NULL` to omit (default).

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra operations. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical or `NA`. If `TRUE`, prefer on-disk processing for heavy ops.
  Default `TRUE`.

- force_gc:

  Logical; call [`gc()`](https://rdrr.io/r/base/gc.html) at checkpoints.
  Default `FALSE`.

- return_visible:

  Logical. If `TRUE`, return the summary visibly; otherwise invisibly.
  Default `FALSE`.

- quiet:

  Suppress progress prints ([`cat()`](https://rdrr.io/r/base/cat.html))?
  Default `FALSE`.

## Value

A **data.frame** with columns: `output`, `gap_count`,
`max_gap_distance`, `filter_size_cells`, `smoothed`,
`smoothing_radius_used`, `elapsed_sec`.

## Details

**Workflow**

1.  **Inputs**: accepts paths or in-memory objects (`SpatRaster` /
    `sf`).

2.  **Fast crop**: transforms the grid bbox to the template CRS, expands
    by `buffer_m`, crops the raw raster early to minimise work.

3.  **Method auto-choice**: `"near"` for integer/factor rasters,
    otherwise `"bilinear"`.

4.  **Align**: `terra::project(raw_crop, template, method)` then
    `terra::mask(..., template)`.

5.  **Gaps**: optional internal NA counts on the template footprint
    only.

6.  **Gap fill (optional)**: Whitebox `wbt_fill_missing_data`, with
    window auto-sized from the widest gap (or user `filter_size_cells`).

7.  **Smoothing (optional)**: IDW via
    [`terra::interpIDW()`](https://rspatial.github.io/terra/reference/interpIDW.html)
    on points derived from an aggregated version of the filled raster;
    masked back to template.

8.  **Plotting**: `plot_gaps`, `plot_result` (side-by-side if both).

9.  **Write**: atomic write (tmp + rename), LZW, tiling,
    `BIGTIFF=IF_SAFER`, threaded; output CRS **forced to template
    WKT/EPSG** so `crs(out) == crs(template)` is TRUE.

10. **Sink safety & memory**: snapshots/restore sinks (prevents stuck
    console if interrupted); `terraOptions(memfrac/tempdir/todisk)` set
    then restored; optional [`gc()`](https://rdrr.io/r/base/gc.html).

**IDW smoothing** Aggregate to points to IDW back to template. Controls:

- `smooth_radius_km` search radius (km; converted to template units),

- `smooth_agg_factor` aggregation factor before point creation,

- `smooth_power`, `smooth_epsilon`, `smooth_nmax`. Set
  `smooth_force=TRUE` to allow smoothing of integer/factor rasters (not
  recommended).

**Mass preservation tips** The smoothed result will generally *not*
preserve totals/means.

- **Global**: scale smoothed raster to match global sum.

- **Zonal**: compute per-zone factors on the original coarse grid and
  multiply.

## See also

terra::project(), terra::mask(), terra::interpIDW(),
whitebox::wbt_fill_missing_data(), sfarrow::st_read_parquet()

## Examples

``` r
if (FALSE) { # \dontrun{
out_dir <- tempdir()
df <- downscale2egv(
  template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
  grid_path     = "./Templates/TemplateGrids/grid_10km.parquet",
  rawfile_path  = "./MyCoarse/coarse_indicator.tif",
  out_path      = out_dir,
  file_name     = "indicator_egv.tif",
  layer_name    = "indicator",
  fill_gaps     = TRUE,
  smooth        = TRUE,
  smooth_radius_km = 10,
  plot_result   = TRUE
)
print(df)
} # }
```
