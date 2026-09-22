# Create buffered tiles from point layers

Buffers point parquet files and writes **per-tile** outputs named:

- Fixed radius: `<base>_r<radius>_<tileid>.parquet`

- Per-feature radii (`radius_field`):
  `<base>_varradius_<tileid>.parquet`

Where:

- **base** = prefix of input filename before the first underscore (e.g.,
  "pts100" from `pts100_sauzeme.parquet`)

- **tileid** = unique value of `split_field` (default "tks50km") inside
  the input file

The `buffer_mode` determines how radii are assigned:

- **"dense"**: Buffers the best-matching `pts100*.parquet` (prefers
  `pts100_sauzeme.parquet`) for each tile by `radii_dense` (default:
  500, 1250, 3000, 10000 m).

- **"sparse"**: Uses a file to radius mapping. Default mapping:

  - "pts100_sauzeme.parquet" to c(500, 1250)

  - "pts300_sauzeme.parquet" to 3000

  - "pts1000_sauzeme.parquet" to 10000 You can override this via
    `mapping_sparse` (named list or data.frame with columns `file`,
    `radius`).

- **"specified"**: You provide `points_path` and either `buffer_radius`
  (uniform, one or more fixed radii) **or** `radius_field` (numeric
  column with per-feature radii in meters). Outputs are still split by
  `split_field`.

## Usage

``` r
tiled_buffers(
  in_dir = "./Templates/TemplateGridPoints",
  out_dir = "./Templates/TemplateGridPoints/lapas",
  buffer_mode = "dense",
  radii_dense = c(500, 1250, 3000, 10000),
  mapping_sparse = list(pts100_sauzeme.parquet = c(500, 1250), pts300_sauzeme.parquet =
    3000, pts1000_sauzeme.parquet = 10000),
  points_path = NULL,
  buffer_radius = NULL,
  radius_field = NULL,
  split_field = "tks50km",
  n_workers = max(1L, parallel::detectCores()),
  future_max_mem_gb = 4,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Source

Zenodo grids/points example: https://zenodo.org/records/14277114

## Arguments

- in_dir:

  Character. Directory containing input points parquet files (default
  "./Templates/TemplateGridPoints"). Used in "dense" and "sparse" modes.

- out_dir:

  Character. Output directory for buffered tiles (default
  "./Templates/TemplateGridPoints/lapas").

- buffer_mode:

  Character. One of "dense", "sparse", "specified". Default "dense".

- radii_dense:

  Numeric vector of radii (m) used when `buffer_mode = "dense"`. Default
  c(500, 1250, 3000, 10000).

- mapping_sparse:

  Named list **or** data.frame describing file to radii for
  `buffer_mode = "sparse"`. Default:
  `list("pts100_sauzeme.parquet" = c(500, 1250), "pts300_sauzeme.parquet" = 3000, "pts1000_sauzeme.parquet" = 10000)`.
  If a data.frame is supplied, it must have columns `file` and `radius`.

- points_path:

  Character. Path to a single parquet file for
  `buffer_mode = "specified"`.

- buffer_radius:

  Numeric vector. Used in "specified" when you want fixed radii. Ignored
  if `radius_field` is provided.

- radius_field:

  Character or NULL. Column name in `points_path` that provides
  per-feature radii (meters) for "specified". If given, `buffer_radius`
  is ignored.

- split_field:

  Character. Field in the point data that defines tiles. Default
  "tks50km".

- n_workers:

  Integer. Parallel workers. Default `max(1L, parallel::detectCores())`.

- future_max_mem_gb:

  Numeric. Max size of exported globals per worker (GiB). Sets
  `options(future.globals.maxSize = future_max_mem_gb * 1024^3)`.
  Default 4.

- overwrite:

  Logical. Overwrite existing outputs? Default FALSE.

- quiet:

  Logical. Suppress messages? Default FALSE.

## Value

Invisibly returns a data.frame with columns: `input`, `tileid`, `mode`,
`radius_m` (NA if `radius_field` is used), `radius_field`, `out_file`,
`wrote`.

## Details

- Uses sfarrow for parquet I/O and sf for buffering.

- Jobs are created **per tile** (unique `split_field` value) and (when
  applicable) per radius.

- Workers open data from file paths (keeps RAM low).

- Files are written atomically (temp file then move).

- **Safety**: jobs are de-duplicated by content **and** by predicted
  output path, so concurrent workers never write to the same file.

## See also

[`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tiled_buffers(buffer_mode = "sparse", split_field = "tks50km")
} # }
```
