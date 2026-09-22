# Compute landscape-level metrics per zone (tiled), merge, analyze gaps, optionally fill with IDW, and write rasters

Computes a landscapemetrics metric (default `"lsm_l_shdi"`), optionally
with extra `lm_args`, that yields one value per zone and **per input
layer**. Runs tile-by-tile (by `tile_field`), writes per-tile rasters,
merges to final per-layer GeoTIFF(s), then performs **gap analysis** (NA
count within the template footprint and optional maximum gap width) and
**optional IDW gap filling** via WhiteboxTools. Returns a compact
**data.frame** with per-layer stats and timing.

## Usage

``` r
landscape_function(
  landscape,
  zones = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
  id_field = "rinda500",
  tile_field = "tks50km",
  template = "./Templates/TemplateRasters/LV500m_10km.tif",
  out_dir,
  out_filename,
  out_layername,
  what = "lsm_l_shdi",
  lm_args = NULL,
  buffer_m = 1000,
  rasterize_engine = c("fasterize", "terra"),
  n_workers = 1,
  future_max_size = 4 * 1024^3,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
  write_datatype = NULL,
  NAflag = NA_real_,
  keep_tiles = FALSE,
  skip_existing = TRUE,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = TRUE,
  force_gc = FALSE,
  quiet = FALSE,
  report_gaps = TRUE,
  report_gap_size = TRUE,
  fill_gaps = FALSE,
  idw_weight = 2,
  filter_size_cells = "auto",
  plot_result = FALSE,
  plot_gaps = FALSE
)
```

## Arguments

- landscape:

  SpatRaster or path(s) to .tif (class labels). Multi-layer supported.

- zones:

  sf polygons or path to sfarrow GeoParquet. Default
  `"./Templates/TemplateGrids/tikls500_sauzeme.parquet"`.

- id_field:

  Zone id field (default `"rinda500"`).

- tile_field:

  Tiling field (default `"tks50km"`).

- template:

  SpatRaster or path defining target grid/CRS/cellsize. Default
  `"./Templates/TemplateRasters/LV500m_10km.tif"`.

- out_dir:

  Output directory (created if missing).

- out_filename:

  Character vector: final filename(s) (no dir); one per input layer,
  matching order.

- out_layername:

  Character vector: final layer name(s); one per input layer, matching
  order.

- what:

  Landscapemetrics metric (default `"lsm_l_shdi"`).

- lm_args:

  Named list of extra args for
  [`landscapemetrics::sample_lsm()`](https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html)
  compatible with `what`.

- buffer_m:

  Numeric buffer (m) for tile bbox before crop (default `1000`).

- rasterize_engine:

  `"fasterize"` (default) or `"terra"`.

- n_workers:

  Integer; `1` = sequential.

- future_max_size:

  Max globals per worker (bytes); default `4 * 1024^3` (~4 GiB).

- gdal_opts:

  GDAL creation options (merged with tuned defaults:
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`).

- write_datatype:

  terra datatype string; if `NULL`, defaults to `"FLT4S"`.

- NAflag:

  Numeric NA value to write; if `NA`/`NULL`, a sensible default is
  chosen (`FLT* to -9999`, `INT2S to -32768`, `INT4S to -2147483648`).

- keep_tiles:

  Logical; keep temp tiles dir (debug). Default `FALSE`.

- skip_existing:

  Logical; skip tiles/finals that already exist. Default `TRUE`.

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra ops. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical or `NA`. If `TRUE`, prefer on-disk. Default `TRUE`.

- force_gc:

  Logical; call [`gc()`](https://rdrr.io/r/base/gc.html) at key
  checkpoints. Default `FALSE`.

- quiet:

  Suppress progress prints ([`cat()`](https://rdrr.io/r/base/cat.html))?
  Default `FALSE`.

- report_gaps:

  Logical; print/compute `gap_count`. Default `TRUE`.

- report_gap_size:

  Logical; compute/print maximum gap width (expensive). Default `TRUE`.

- fill_gaps:

  Logical; run Whitebox IDW fill if gaps exist. Default `FALSE`.

- idw_weight:

  Numeric; IDW power for Whitebox fill. Default `2`.

- filter_size_cells:

  Integer or `"auto"`; Whitebox window size (cells). Default `"auto"`.

- plot_result:

  Logical; plot final per-layer raster(s). Default `FALSE`.

- plot_gaps:

  Logical; plot gap map(s) (1 = NA inside template). Default `FALSE`.

## Value

A **data.frame** (one row per layer) with:

- `layer_name`, `output_path`,

- `tiles_written`, `tiles_skipped_existing`, `tiles_skipped_empty_crop`,
  `tiles_skipped_empty_join`,

- `merge_skipped`, `gap_count`, `max_gap_distance`,
  `filter_size_cells_used`, `gap_filled`,

- `n_tiles`, `n_zones`, `n_layers`, `elapsed_sec`. Attributes:
  `"tile_dir"` (path or `NULL`), `"run_params"` (list).

## Details

**Workflow**

1.  **Inputs & CRS**

    - Accept paths or in-memory (`SpatRaster`, `sf`). Zones and
      landscape are aligned to the **template** CRS/grid.

2.  **Tiling**

    - Split `zones` by `tile_field`; process each tile independently.

3.  **Per-tile (all layers)**

    - Quick extent crop (buffered by `buffer_m`),
      [`landscapemetrics::sample_lsm()`](https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html)
      per layer, join back by `id_field`, rasterize (engine:
      `"fasterize"` or `"terra"`), mask to template crop, and write tile
      raster. Skips empty crops/joins cleanly.

4.  **Merge per layer**

    - VRT over all tile rasters to final GeoTIFF with tuned GDAL options
      (LZW, tiling, BIGTIFF=IF_SAFER, threads).

    - CRS is **forced to the template WKT** so
      `crs(output) == crs(template)` string-matches.

5.  **Gap analysis & optional fill (per layer)**

    - `gap_count`: number of NA cells **inside** the template footprint.

    - If `report_gap_size=TRUE` or `filter_size_cells="auto"`, compute
      **max gap width** (distance to nearest non-NA).

    - If `fill_gaps=TRUE` and gaps exist, fill with
      [`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html)
      using `idw_weight` and `filter_size_cells` (or auto from the max
      gap width; odd, \>=3).

    - Plot result and/or gaps if requested (side-by-side if both).

6.  **Return**

    - One row per output layer with paths, tile counters, gap stats,
      fill parameters, and elapsed seconds.

## See also

landscapemetrics::sample_lsm, terra::writeRaster, fasterize::fasterize,
whitebox::wbt_fill_missing_data

## Examples

``` r
if (FALSE) { # \dontrun{
res_tbl <- landscape_function(
  landscape      = "./Rastri_10m/Ainava_KopejaiDaudzveidibai.tif",
  zones          = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
  id_field       = "rinda500",
  tile_field     = "tks50km",
  template       = "./Templates/TemplateRasters/LV500m_10km.tif",
  out_dir        = "./out/lsm/",
  out_filename   = "Landscape_diversity.tif",
  out_layername  = "Landscape_diversity",
  what           = "lsm_l_shdi",
  rasterize_engine = "fasterize",
  n_workers      = 8,
  fill_gaps      = TRUE,
  plot_gaps      = TRUE,
  plot_result    = TRUE
)
print(res_tbl)

#' ## --- Total edge length per zone (ignore outside map boundary) ---
## Using a binary "water vs other" landscape, compute lsm_l_te per 500 m zone.
rez_edges <- landscape_function(
  landscape        = "./Rastri_10m/Ainava_vienk_mask.tif",
  zones            = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
  id_field         = "rinda500",
  tile_field       = "tks50km",
  template         = "./Templates/TemplateRasters/LV500m_10km.tif",
  out_dir          = "./",
  out_filename     = "edges_water.tif",
  out_layername    = "edges_water",
  what             = "lsm_l_te",
  lm_args          = list(count_boundary = FALSE),
  rasterize_engine = "terra",
  n_workers        = 8,
  future_max_size  = 2 * 1024^3,
  report_gaps      = TRUE,
  plot_result      = TRUE,
  plot_gaps        = TRUE
)
rez_edges
} # }
```
