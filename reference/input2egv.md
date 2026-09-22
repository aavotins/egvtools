# Convert an input raster to an EGV-aligned raster

Align a fine-resolution input raster to a (coarser) EGV template,
optionally cover missing values and/or fill gaps (IDW via Whitebox), and
write the result to disk. Designed for large runs: fast gap counting
(inside template footprint only), optional filling, tuned GDAL write
options, and controlled `terra` memory/temp behavior.

## Usage

``` r
input2egv(
  input,
  egv_template,
  summary_function = "average",
  missing_job = c("CoverInput", "CoverOutput", "FillOutput", "CoverInput-CoverOutput",
    "CoverInput-FillOutput", "none"),
  input_bg = NULL,
  input_template = NULL,
  output_bg = NULL,
  idw_weight = 2,
  is_categorical = FALSE,
  check_alignment = TRUE,
  outlocation = "./Rastri_100m/RAW/",
  outfilename,
  layername,
  plot_gaps = FALSE,
  plot_final = FALSE,
  NAflag = NULL,
  gdal_opts = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  write_datatype = NULL,
  terra_memfrac = 0.7,
  terra_tempdir = tempdir(),
  terra_todisk = TRUE,
  force_gc = FALSE,
  return_visible = FALSE,
  quiet = FALSE
)
```

## Arguments

- input:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or character path (fine resolution input).

- egv_template:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or path (coarser EGV template).

- summary_function:

  Character resample method. Default `"average"`; `"mean"` is mapped to
  `"average"`. If `is_categorical = TRUE`, `"near"` is enforced.

- missing_job:

  One of `"CoverInput"`, `"CoverOutput"`, `"FillOutput"`,
  `"CoverInput-CoverOutput"`, `"CoverInput-FillOutput"`, `"none"`.
  Default `"none"`.

- input_bg:

  Background for CoverInput: `SpatRaster`, path, or **numeric** (+
  `input_template`).

- input_template:

  Template used when `input_bg` is numeric (builds a background raster).

- output_bg:

  Background for CoverOutput: `SpatRaster`, path, or **numeric** (+
  `egv_template`).

- idw_weight:

  Numeric power for Whitebox IDW. Default `2`.

- is_categorical:

  Logical; if `TRUE`, use `"near"` and skip WBT fill. Default `FALSE`.

- check_alignment:

  Logical; if `TRUE`, cheaply checks CRS/res; warns if mismatched.
  Default `TRUE`.

- outlocation:

  Output directory. Default `"./Rastri_100m/RAW/"`.

- outfilename:

  Output filename (e.g., `"layer.tif"`). **Required**.

- layername:

  Layer name to assign before writing. **Required**.

- plot_gaps:

  Logical; plot the **initial** gap map (guarded by
  [`interactive()`](https://rdrr.io/r/base/interactive.html)). Default
  `FALSE`.

- plot_final:

  Logical; plot the final raster (guarded by
  [`interactive()`](https://rdrr.io/r/base/interactive.html)). Default
  `FALSE`.

- NAflag:

  Optional numeric NA flag for writing. Default `NULL` (terra default).

- gdal_opts:

  Character vector of GDAL creation options (merged with tuned
  defaults). Default
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Optional terra datatype (e.g., `"FLT4S"`, `"INT2S"`). Default `NULL`.

- terra_memfrac:

  Memory fraction for `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra/Whitebox operations. Default:
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical (set `terraOptions(todisk=...)` for this run). Default `TRUE`.

- force_gc:

  Logical; if `TRUE`, call [`gc()`](https://rdrr.io/r/base/gc.html) at
  checkpoints. Default `FALSE`.

- return_visible:

  Logical; if `TRUE`, return data.frame is visible; otherwise invisible.
  Default `FALSE`.

- quiet:

  Logical; suppress console messages. Default `FALSE`.

## Value

A **single-row `data.frame`** with:

- `path` (character): written file path

- `method` (character): resample method used

- `missing_job` (character)

- `n_gaps_initial` (integer; after resample/CoverInput, before
  output-stage ops)

- `n_gaps_final` (integer; after CoverOutput/FillOutput)

- `max_gap_dist` (numeric; `NA` if not computed)

- `filter_w` (integer; `NA` if not used)

- `elapsed_sec` (numeric)

## Details

**Workflow (high level)**

1.  Load `input` and `egv_template`. Optionally warn on CRS/resolution
    mismatch (`check_alignment`).

2.  If `missing_job` contains **CoverInput**, build `input_bg`
    (numeric + `input_template`, or raster/path) and
    [`terra::cover()`](https://rspatial.github.io/terra/reference/cover.html)
    the input.

3.  **Resample/aggregate** the input to the template with
    `summary_function` (`"mean"` is linked to `"average"`). If
    `is_categorical=TRUE` the method is forced to `"near"`.

4.  Count **initial gaps** = `is.na(aligned) & !is.na(template)`.

5.  If `missing_job` contains **CoverOutput**, apply
    [`terra::cover()`](https://rspatial.github.io/terra/reference/cover.html)
    using `output_bg`.

6.  If `missing_job` contains **FillOutput** **and** `!is_categorical`,
    estimate a **maximum gap width** (via
    [`terra::distance()`](https://rspatial.github.io/terra/reference/distance.html)
    on a fillable mask) and set the Whitebox **filter width** as
    `ceil(max_gap/pixel_size)*2`.

7.  Mask to the template footprint, set the final layer name, optionally
    plot initial gaps and/or result, and write with LZW/tiling (adds a
    `PREDICTOR` appropriate to datatype if missing).

8.  Return a single-row **data.frame** with path, method, gap counts,
    max gap, filter size, and elapsed time.

**Missing value workflows** (`missing_job`):

- `"CoverInput"`: cover the **input** before resampling (with
  `input_bg`).

- `"CoverOutput"`: cover the **resampled** raster (with `output_bg`).

- `"FillOutput"`: if gaps remain after resampling, compute a max-gap
  distance (via
  [`terra::distance`](https://rspatial.github.io/terra/reference/distance.html))
  and call
  [`whitebox::wbt_fill_missing_data()`](https://rdrr.io/pkg/whitebox/man/wbt_fill_missing_data.html)
  with a filter width approx. twice the maximum gap (in pixels). Only a
  minimum clamp (`>= 3`) is applied.

- Combinations: `"CoverInput-CoverOutput"`, `"CoverInput-FillOutput"`.

- `"none"`: do nothing about gaps.

**Resampling**:

- Default `summary_function = "average"` (area-preserving). `"mean"` is
  mapped to `"average"`.

- If `is_categorical = TRUE`, method is forced to `"near"`.

**I/O and stability**:

- Writes are tiled/ threaded with `LZW` compression. Adds `PREDICTOR=2`
  (floats) or `PREDICTOR=3` (ints) if not provided.

- `terraOptions(memfrac, tempdir, todisk)` are set for the call and
  **restored** on exit.

- Console messages use [`cat()`](https://rdrr.io/r/base/cat.html) via a
  `say()` helper and the call is **sink-safe**.

## See also

[`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md),
[`downscale2egv()`](https://aavotins.github.io/egvtools/reference/downscale2egv.md),
[`distance2egv()`](https://aavotins.github.io/egvtools/reference/distance2egv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
input2egv(
  input            = "./TestejuPakotni_early/Forests_StandAge_bg0.tif",
  egv_template     = "./Templates/TemplateRasters/LV100m_10km.tif",
  summary_function = "average",
  missing_job      = "CoverInput-FillOutput",
  input_bg         = 0,                     # numeric background + input_template
  input_template   = "./Templates/TemplateRasters/LV10m_10km.tif",
  outlocation      = "./Rastri_100m/RAW/",
  outfilename      = "StandAge_100m.tif",
  layername        = "stand_age",
  force_gc         = TRUE,
  quiet            = FALSE
)
} # }
```
