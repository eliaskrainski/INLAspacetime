#' Spatial and Spatio-Temporal Models using INLA
#'
#' This package main purpose is to provide user friendly functions
#' to fit temporal, spatial and spatio-temporal models using the
#' INLA software available at www.r-inla.org
#'
#' @return opens the Vignettes directory on a browser
#' @export
#'
#' @import methods
#' @import Matrix
#' @importFrom sp CRS

INLAspacetime <- function() {
  print("Welcome to the INLAspacetime package!")
##  utils::browseVignettes("INLAspacetime")
  utils::browseURL("https://eliaskrainski.github.io/INLAspacetime")
}
