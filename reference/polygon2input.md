# Rasterize polygons to a template grid, optionally restrict & cover gaps, then write GeoTIFF

Rasterizes polygon/multipolygon **sf** data to a raster aligned to a
**template** GeoTIFF. Rasterization targets a `raster::RasterLayer`
built from the template (so grids normally match). Projection is
optional (`project_mode`). Missing values are counted **only** over
valid template cells. You may optionally restrict the result with a
raster mask (`restrict_to`) using numeric values or **bracketed range
strings** (e.g., `"(0,5]"`, `"[10,)"`). Remaining `NA` cells can be
filled by covering with a background raster (`background_raster`) or a
constant (`background_value`). For large rasters, heavy steps
(projection/mask/cover) can stream to disk via `terra_todisk=TRUE`.

## Usage

``` r
polygon2input(
  vector_data,
  template_path,
  out_path = "./",
  file_name,
  value_field = NULL,
  constant_value = 1,
  fun = "max",
  value_type = c("categorical", "continuous"),
  project_mode = c("auto", "never", "always"),
  prepare = TRUE,
  restrict_to = NULL,
  restrict_values = NULL,
  background_raster = NULL,
  background_value = NULL,
  check_na = TRUE,
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
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- vector_data:

  sf polygons/multipolygons to rasterize.

- template_path:

  Path to the template raster (.tif).

- out_path:

  Output directory. Default `"./"`.

- file_name:

  Output filename, e.g., `"my_input.tif"`.

- value_field:

  Character or `NULL`. Attribute to burn; if `NULL`, uses
  `constant_value`.

- constant_value:

  Numeric scalar burned when `value_field=NULL`. Default `1`.

- fun:

  Aggregation in `fasterize` for overlaps: one of `"max"`, `"sum"`,
  `"first"`, `"last"`, `"min"`, `"count"`. Default `"max"`.

- value_type:

  Guides resampling **if projection happens**: `"categorical"` with
  method `"near"`, `"continuous"` with `"bilinear"`. Default
  `"categorical"`.

- project_mode:

  `"auto"` (default; project only if misaligned), `"never"`, or
  `"always"`.

- prepare:

  Logical. If `TRUE`, run `st_make_valid()` and transform to template
  CRS. Default `TRUE`.

- restrict_to:

  Optional raster mask (path or `terra` `SpatRaster`).

- restrict_values:

  Optional values to keep from `restrict_to`: numbers (scalar/vector)
  and/or **range strings** using bracket inclusivity, e.g. `"(0,5]"`,
  `"[10,)"`, `"(-inf,0)"`. Multiple entries are OR-ed. If `NULL`, keep
  non-NA cells of `restrict_to`.

- background_raster:

  Optional path/`SpatRaster` used to cover remaining NAs.

- background_value:

  Numeric constant to fill where result is NA **if** no
  `background_raster`. Default `NULL` (no covering).

- check_na:

  Logical. If `TRUE`, report NA counts before/after. Default `TRUE`.

- plot_result:

  Logical. Plot final raster (after all processing). Default `FALSE`.

- plot_gaps:

  Logical. Plot NA gaps (only within template footprint). If both plot
  flags are `TRUE`, plots are side-by-side. Default `FALSE`.

- NAflag:

  Optional NA flag for writing (passed to GDAL). Default `NULL` (auto if
  needed).

- gdal_opts:

  Character vector of GDAL creation options (merged with tuned
  defaults). Default
  `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.

- write_datatype:

  Optional terra datatype for writing (e.g., `"FLT4S"`, `"INT2S"`).
  Default `NULL` (auto).

- terra_memfrac:

  `terraOptions(memfrac=...)`. Default `0.7`.

- terra_tempdir:

  Temp dir for terra operations. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- terra_todisk:

  Logical or `NA`. If `TRUE`, prefer on-disk processing for heavy ops
  (project/mask/cover). If `FALSE`, prefer memory. If `NA`, leave
  session default unchanged. Default `FALSE`.

- force_gc:

  Logical; call [`gc()`](https://rdrr.io/r/base/gc.html) at checkpoints.
  Default `FALSE`.

- overwrite:

  Overwrite output? Default `FALSE`.

- quiet:

  Suppress progress prints ([`cat()`](https://rdrr.io/r/base/cat.html))?
  Default `FALSE`.

## Value

Invisibly, a list with `out_file`, `n_cells`, `n_na_initial`,
`n_na_final`, `elapsed_sec`, and `crs`.

## Details

**Workflow**

1.  (Optional) `prepare=TRUE` ensures valid polygons and transforms to
    the template CRS
    ([`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html) +
    [`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)).
    Set `prepare=FALSE` if you guarantee this.

2.  Rasterize with
    [`fasterize::fasterize()`](https://rdrr.io/pkg/fasterize/man/fasterize.html)
    onto a `raster::RasterLayer` built from the template.

3.  Depending on `project_mode`, **project** to the template grid:

    - `"auto"` (default): checks CRS/resolution/extent/rows/cols;
      projects **only if needed**.

    - `"never"`: skips projection.

    - `"always"`: forces projection (`"near"` for
      `value_type="categorical"`, `"bilinear"` for `"continuous"`).

4.  Apply a **template mask** (keeps template extent and `NA`
    footprint).

5.  (Optional) **Restrict** by `restrict_to` (path or in-memory
    `SpatRaster`):

    - If `restrict_values` is `NULL`, keep **non-NA** cells of
      `restrict_to`.

    - If supplied, keep cells matching **numbers** (scalar/vector)
      and/or **ranges** using bracket syntax: `"(a,b)"` open; `"[a,b]"`
      closed; mix ends like `"(a,b]"`. Use `-inf`/`+inf` for unbounded,
      e.g., `"[10,)"`, `"(-inf,0)"`. Multiple entries are OR-ed.

6.  (Optional) **Cover** remaining `NA` cells:

    - with `background_raster` (aligned to template; projected if
      needed), or

    - with a **constant** using `background_value` (fast and
      memory-light).

7.  A **final mask** to the template is applied **before plotting and
    saving**.

8.  NA counts are computed with
    [`terra::global()`](https://rspatial.github.io/terra/reference/global.html)
    on raster masks (memory-safe).

**Datatype & NAflag auto-chooser** If you do **not** provide
`write_datatype`/`NAflag`, the function picks:

- `value_type="categorical"` to `write_datatype="INT2S"`,
  `NAflag=-32768`

- `value_type="continuous"` to `write_datatype="FLT4S"`, `NAflag`
  omitted You can always override via arguments.

**Performance & stability**

- For very large rasters (e.g., ~1.3B cells), set `terra_todisk=TRUE`
  and consider a fast SSD for `terra_tempdir`. This streams big
  operations to disk and avoids "vector memory limit" errors.

- Projection is often unnecessary here because rasterization targets the
  template grid; `"auto"` will detect alignment and skip it.

- Writes are **atomic**: a temporary file is written and then moved into
  place.

## Range string syntax for restrict_values

Use `(a,b)` for open interval, `[a,b]` for closed; mix ends like
`(a,b]`. Use `-inf`/`+inf` (or `inf`) for unbounded, e.g. `"[10,)"`,
`"(-inf,0)"`. Supply multiple strings to OR them, e.g.
`c("(0,5]","[10,15)")`.

## See also

[`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md),
[`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic: burn constant "1", continuous output, no covering
polygon2input(
vector_data = my_polys_sf,
template_path = "./Templates/TemplateRasters/template10m.tif",
out_path = "./Outputs",
file_name = "mask_const1.tif",
value_field = NULL,
constant_value = 1,
value_type = "continuous",
project_mode = "auto",
prepare = FALSE,
check_na = TRUE,
plot_result = TRUE,
overwrite = TRUE
)

# Restrict to classes 1 and 2, plus (10,20], then fill remaining NAs with 0
polygon2input(
vector_data = my_polys_sf,
template_path = "./Templates/TemplateRasters/template10m.tif",
out_path = "./Outputs",
file_name = "mask_restricted_bg0.tif",
value_field = "attr",
value_type = "categorical",
restrict_to = "./mask_classes.tif",
restrict_values = c(1, 2, "(10,20]"),
background_value = 0,
terra_todisk = TRUE, # stream big ops to disk
terra_tempdir = tempdir(), # or a fast SSD scratch
check_na = TRUE,
plot_result = TRUE,
plot_gaps = TRUE,
overwrite = TRUE
)
} # }
```
