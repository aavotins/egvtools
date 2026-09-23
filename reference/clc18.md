# Example Corine Land Cover (2018) polygons

A reduced Corine Land Cover (2018) dataset used in spatial examples.

## Usage

``` r
clc18
```

## Format

An `sf` object containing polygon geometries and grid identifiers:

- `OBJECTID`: unique identificator of geometry.

- `code_18`: Corine Land Cover third level class' identificator.

- `Remark`: redundant field from original database.

- `Area_Ha`: redundant field from original database.

- `Shape_Length`: redundant field from original database.

- `Shape_Area`: redundant field from original database.

- `geom`: `sf geometry` field.

## Source

Derived from a reduced section of the Latvian 500-m harmonization grid.
