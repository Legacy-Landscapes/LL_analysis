project_for_sf <- function(func, sf_data, ...) {
  # sf can not correctly project longlat format, so we project to
  # mercator, simplify and then project back to longlat
  crs <- sf::st_crs(sf_data)
  sf_data <- sf::st_transform(sf_data, 3857)
  sf_data <- func(sf_data, ...)
  sf_data <- sf::st_transform(sf_data, crs)
  return(sf_data)
}

simplify_polygons <- function(sf_data, tolerance = 0.1) {
  return(
    project_for_sf(
      sf::st_simplify,
      sf_data,
      preserveTopology = TRUE,
      dTolerance = tolerance
    )
  )
}

load_worldmap <- function(filename) {
  worldmap <- sf::st_read(filename, layer = "ne_50m_admin_0_countries")
  worldmap <- simplify_polygons(worldmap)
  return(worldmap)
}

load_realmmap <- function(filename) {
  worldmap <- sf::st_read(filename, layer = "realms")
  worldmap <- simplify_polygons(worldmap)
  return(worldmap)
}

