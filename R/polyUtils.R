#' Internal util functions for polygon properties.
#'
#' Function used internally to compute the area of a triangle.
#' @rdname polyUtils
#' @aliases Heron
#' @description
#' This computes the area of a triangle given its three coordinates.
#' @param x,y coordinate vectors.
#' @section Warning: Internal functions, not exported.
#' @return the area of a 2d triangle
Heron <- function(x, y) {
  ### function to compute the area of a triangle
  aa <- sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2)
  bb <- sqrt((x[3] - x[2])^2 + (y[3] - y[2])^2)
  cc <- sqrt((x[1] - x[3])^2 + (y[1] - y[3])^2)
  s <- 0.5 * (aa + bb + cc)
  sqrt(s * (s - aa) * (s - bb) * (s - cc))
}
#' Function to compute the area of a polygon.
#' @rdname polyUtils
#' @aliases Area
#' @return the area of a 2d polygon
Area <- function(x, y) {
  n <- length(x)
  stopifnot(length(y) == n)
  abs(0.5 * sum(x[1:n] * y[c(2:n, 1)] -
    y[1:n] * x[c(2:n, 1)]))
}
#' Function to compute area of a S2 triangle
#' @rdname polyUtils
#' @aliases s2trArea
#' @param tr the triangle coordinates
#' @param R the radius of the spherical domain
#' @return the area of a triangle in S2
s2trArea <- function(tr, R = 1) {
  costh <- R * R +
    sum(tr[1, ] * tr[2, ]) +
    sum(tr[2, ] * tr[3, ]) +
    sum(tr[3, ] * tr[1, ])
  sinth <- (
    (tr[2,2]*tr[3,3] - tr[2,3]*tr[3,2]) * tr[1,1] +
    (tr[2,3]*tr[3,1] - tr[2,1]*tr[3,3]) * tr[1,2] +
    (tr[2,1]*tr[3,2] - tr[2,2]*tr[3,1]) * tr[1,3] ) / R
  return((2 * R * R) * atan2(sinth, costh))
}
#' Function to compute flat area approximation
#' @rdname polyUtils
#' @aliases flatArea
#' @return the area of a triangle
flatArea <- function(tr) {
  e0 <- tr[3, ] - tr[2, ]
  e1 <- tr[1, ] - tr[3, ]
  s0 <- e0[2] * e1[3] - e0[3] * e1[2]
  s1 <- e0[3] * e1[1] - e0[1] * e1[3]
  s2 <- e0[1] * e1[2] - e0[2] * e1[1]
  return(sqrt(s0 * s0 + s1 * s1 + s2 * s2) * 0.5)
}
#' Function to compute the stiffness contribution from a triangle.
#' @rdname polyUtils
#' @aliases Stiffness
#' @return the stiffness matrix for a triangle
Stiffness <- function(tr) {
  e <- cbind(tr[3, ] - tr[2, ],
             tr[1, ] - tr[3, ],
             tr[2, ] - tr[1, ])
  return(0.25 * crossprod(e))
}
