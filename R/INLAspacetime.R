#' Functions for easier code to fit of some models with INLA.
#'
#' This package main purpose is to provide user friendly functions
#' to fit temporal, spatial and spatio-temporal models using the
#' INLA software available at www.r-inla.org
#'
#' @return opens the Vignettes directory on a browser
#' @export
#'
#' @importFrom methods as is

INLAspacetime <- function() {
  print("Welcome to the INLAspacetime package!")
  utils::browseVignettes("INLAspacetime")
}
