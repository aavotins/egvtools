#' Example 100-m grid points
#'
#' Point centres of a reduced 100-m grid used in spatial examples.
#'
#' @format An `sf` object containing point geometries and grid identifiers:
#'
#' - `id`: identificator of 100 m grid cell.
#'
#' - `yes`: indicator for location in terrestrial Latvia (all the cells contain value `1`).
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `rinda300`: identificator of 300 m grid cell.
#'
#' - `ID1km`: identificator of 1 km grid cell.
#'
#' - `rinda500`: identificator of 500 m grid cell.
#'
#' - `geom`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 100-m harmonization grid.
"pts100_sauzeme"

#' Example 1000-m grid points
#'
#' Point centres of a reduced 1000-m grid used in spatial examples.
#'
#' @format An `sf` object containing point geometries and grid identifiers:
#'
#' - `ID1km`: identificator of 1 km grid cell.
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `geometry`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 1000-m harmonization grid.
"pts1000_sauzeme"

#' Example 300-m grid points
#'
#' Point centres of a reduced 300-m grid used in spatial examples.
#'
#' @format An `sf` object containing point geometries and grid identifiers:
#'
#' - `rinda300`: identificator of 300 m grid cell.
#'
#' - `X`: EPSG:3059 x-coordinate of centroid.
#'
#' - `Y`: EPSG:3059 y-coordinate of centroid.
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `x`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 300-m harmonization grid.
"pts300_sauzeme"

#' Example 500-m grid points
#'
#' Point centres of a reduced 500-m grid used in spatial examples.
#'
#' @format An `sf` object containing point geometries and grid identifiers:
#'
#' - `rinda500`: identificator of 500 m grid cell.
#'
#' - `X`: EPSG:3059 x-coordinate of centroid.
#'
#' - `Y`: EPSG:3059 y-coordinate of centroid.
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `x`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 500-m harmonization grid.
"pts500_sauzeme"

#' Example 100-m grid polygons
#'
#' A reduced 100-m grid used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `id`: identificator of 100 m grid cell.
#'
#' - `yes`: indicator for location in terrestrial Latvia (all the cells contain value `1`).
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `rinda300`: identificator of 300 m grid cell.
#'
#' - `ID1km`: identificator of 1 km grid cell.
#'
#' - `rinda500`: identificator of 500 m grid cell.
#'
#' - `geom`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 100-m harmonization grid.
"tikls100_sauzeme"

#' Example 1000-m grid polygons
#'
#' A reduced 1000-m grid used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `yes`: indicator for location in terrestrial Latvia (all the cells contain value `1`).
#'
#' - `ID1km`: identificator of 1 km grid cell.
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `geometry`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 1000-m harmonization grid.
"tikls1km_sauzeme"

#' Example 300-m grid polygons
#'
#' A reduced 300-m grid used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `rinda300`: identificator of 300 m grid cell.
#'
#' - `x`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 300-m harmonization grid.
"tikls300_sauzeme"

#' Example 500-m grid points
#'
#' Point centres of a reduced 500-m grid used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `rinda500`: identificator of 500 m grid cell.
#'
#' - `tks50`: identificator of topographic map's 50 km page.
#'
#' - `x`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 500-m harmonization grid.
"tikls500_sauzeme"

#' Example 50-km map tiles
#'
#' A reduced 50-km grid used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `NOSAUKUMS`: name of topographic map's 50 km page.
#'
#' - `NUMURS`: identificator of topographic map's 50 km page.
#'
#' - `Shape_Length`: redundant field from original database.
#'
#' - `Shape_Area`: redundant field from original database.
#'
#' - `Shape`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 500-m harmonization grid.
"tks93_50km"

#' Example Corine Land Cover (2018) polygons
#'
#' A reduced Corine Land Cover (2018) dataset used in spatial examples.
#'
#' @format An `sf` object containing polygon geometries and grid identifiers:
#'
#' - `OBJECTID`: unique identificator of geometry.
#'
#' - `code_18`: Corine Land Cover third level class' identificator.
#'
#' - `Remark`: redundant field from original database.
#'
#' - `Area_Ha`: redundant field from original database.
#'
#' - `Shape_Length`: redundant field from original database.
#'
#' - `Shape_Area`: redundant field from original database.
#'
#' - `geom`: `sf geometry` field.
#'
#' @source Derived from a reduced section of the Latvian 500-m harmonization grid.
"clc18"

