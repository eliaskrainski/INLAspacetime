#' Spatial and Spatio-Temporal Models using INLA
#'
#' This package main purpose is to provide user friendly functions
#' to fit temporal, spatial and space-time models using the
#' INLA software available at www.r-inla.org as well the
#' inlabru package available
#' @return opens the Vignettes directory on a browser
#' @useDynLib INLAspacetime, .registration = TRUE
#' @import methods
#' @import Matrix
#' @import fmesher
#' @export
INLAspacetime <- function() {
  print("Welcome to the INLAspacetime package!")
  utils::browseURL("https://eliaskrainski.github.io/INLAspacetime")
}
.onAttach <- function(...) {
  packageStartupMessage(
   "see more on https://eliaskrainski.github.io/INLAspacetime"
  )
}
