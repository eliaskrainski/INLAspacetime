#' Internal util functions
#'
#' Function used internally to compute the area of a triangle.
#' @rdname polyUtils
#' @aliases Heron
#' @description
#' This computes the area of a triangle given its three coordinates.
#' @param x,y coordinate vectors.
#' @section Warning: Internal functions, not exported.
#' @return the area of a triangle
Heron <- function(x, y) {
  ### function to compute the area of a triangle
  aa <- sqrt((x[2]-x[1])^2 + (y[2]-y[1])^2)
  bb <- sqrt((x[3]-x[2])^2 + (y[3]-y[2])^2)
  cc <- sqrt((x[1]-x[3])^2 + (y[1]-y[3])^2)
  s <- 0.5*(aa+bb+cc)
  sqrt(s*(s-aa)*(s-bb)*(s-cc))
}
#' Function to compute the area of a polygon.
#' @rdname polyUtils
#' @aliases Area
#' @return the area of a general polygon
Area <- function(x, y) {
  n <- length(x)
  stopifnot(length(y)==n)
  abs(0.5*sum(x[1:n]*y[c(2:n,1)]-
              y[1:n]*x[c(2:n,1)]))
}
#' Function to compute the stiffness contribution from a triangle.
#' @rdname polyUtils
#' @aliases Stiffness
#' @return the stiffness matrix for a triangle
Stiffness <- function(x, y) {
  d <- rbind(c(x[3]-x[2], x[1]-x[3], x[2]-x[1]),
             c(y[3]-y[2], y[1]-y[3], y[2]-y[1]))
  crossprod(d)/4
}
