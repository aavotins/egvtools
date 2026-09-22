# Tile a vector grid into parquet tiles

Reads a grid parquet (e.g., `tikls100_sauzeme.parquet`) and writes
per-tile parquet files into `out_dir`. Output filenames are formed as:
**`[base]_[tilename].parquet`**, where:

- **base** = the substring of the input filename before the first
  underscore (e.g., `"tikls100"` from `tikls100_sauzeme.parquet`);

- **tilename** = the value from `tile_field` (e.g., `2434`), or an
  auto-generated chunk id (e.g., `tile_00001`) if `tile_field` is `NULL`
  or not present.

Tiling can be controlled either by an existing `tile_field` (preferred)
or by automatic chunking using `chunk_size`.

## Usage

``` r
tile_vector_grid(
  grid_path = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
  out_dir = "./Templates/TemplateGrids/lapas",
  tile_field = "lapa",
  chunk_size = 50000L,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Source

Zenodo (grid example): https://zenodo.org/records/14277114

## Arguments

- grid_path:

  Character. Path to the source grid parquet (default
  `"./Templates/TemplateGrids/tikls100_sauzeme.parquet"`).

- out_dir:

  Character. Where to write tiles (default
  `"./Templates/TemplateGrids/lapas"`).

- tile_field:

  Character or `NULL`. If provided and exists in data (e.g., `"lapa"`),
  tiles are split by unique values of this field. If `NULL` or not
  found, uses `chunk_size`. Default `"lapa"`.

- chunk_size:

  Integer. Number of rows per tile when `tile_field` is `NULL` or
  missing. Default `50000L`.

- overwrite:

  Logical. Overwrite existing tile files? Default `FALSE`.

- quiet:

  Logical. Suppress messages? Default `FALSE`.

## Value

Invisibly returns a `data.frame` with columns: `tile_id`, `n_rows`,
`path`, `wrote`.

## Details

Uses sfarrow for Parquet I/O. CRS and geometry are preserved. Designed
for stable, idempotent outputs, deterministic filenames, and clean
skipping of existing tiles.

## See also

[`download_vector_templates()`](https://aavotins.github.io/egvtools/reference/download_vector_templates.md),
[`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Split by an existing field (e.g., "lapa")
tile_vector_grid(
  grid_path = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
  out_dir   = "./Templates/TemplateGrids/lapas",
  tile_field = "lapa"
)

# Fallback to chunking if tile_field is NULL or missing
tile_vector_grid(
  grid_path  = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
  out_dir    = "./Templates/TemplateGrids/lapas",
  tile_field = NULL,
  chunk_size = 25000
)
} # }
```
