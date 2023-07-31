#' Spacetime precision matrix.
#'
#' To build the the precision matrix for a spacetime model
#' given the temporal and the spatial meshes.
#'
#' @param tmesh a mesh object over the time domain.
#' @param smesh a mesh object over the spatial domain.
#' @param model a string identifying the model.
#' So far we have the following models:
#' '102', '121', '202' and '220' models.
#' @return a (sparse) precision matrix,
#' as in Lindgren et. al. (2023)
#' @export
stModel.precision <-
  function(smesh, tmesh, model, theta, verbose = FALSE) {

    stopifnot(model %in% c("102", "121", "202", "220"))

    stopifnot(any(substr(smesh$manifold, 1, 1) %in% c("S", "R")))
    Rmanifold <- (substr(smesh$manifold, 1, 1) == "R") + 0L

    dimension <- as.integer(substr(smesh$manifold, 2, 2))
    stopifnot(dimension > 0)

    alphas <- as.integer(strsplit(model, "")[[1]])
    alpha <- alphas[3] + alphas[2] * (alphas[1] - 0.5)
    nu.s <- alpha - dimension / 2
    nu.t <- min(alphas[1] - 0.5, nu.s / alphas[2])

    if (verbose) {
      print(c(alphas = alphas))
      print(c(alpha = alpha, nu.s = nu.s, nu.t = nu.t))
    }

    log.C.t <- lgamma(alphas[1] - 0.5) - lgamma(alphas[1]) - 0.5 * log(4 * pi)
    cc <- c(
      c1 = 0.5 * log(8 * nu.s),
      c2 = -0.5 * log(8 * nu.t),
      c3 = NA
    )
    if (Rmanifold) {
      log.C.Rd <- lgamma(alpha - (dimension * 0.5)) - lgamma(alpha) -
        (dimension * 0.5) * log(4 * pi)
      cc[3] <- 0.5 * (log.C.t + log.C.Rd)
      if (verbose) {
        cat("R manifold, cc[3] = ", cc[3], "\n")
      }
    } else {
      log.C.S2.part <- -log(4 * pi) ## S1???
      cc[3] <- 0.5 * (log.C.t + log.C.S2.part
      ) ## c3 part for S2 (S1???), to be completed in C
      if (verbose) {
        cat("S manifold, cc[3] = ", cc[3], "\n")
      }
    }
    if (verbose) {
      print(c(cc = cc))
    }

    mm <- stModel.matrices(smesh, tmesh, model, constr = FALSE)
    n <- smesh$n * tmesh$n
    nm <- ncol(mm$TT)
    stopifnot(nm == length(mm$bb))
    jmm <- pmatch(paste0("M", 1:nm), names(mm))
    stopifnot(length(jmm[complete.cases(jmm)]) == nm)

    lgammas <- params2gammas(
      theta, alphas[1], alphas[2], alphas[3], smanifold = smesh$manifold)

    params <- double(nm)
    for(i in 1:nm) {
      a1 <- lgammas[1] * mm$TT[1, i]
      a2 <- lgammas[1] * mm$TT[2, i]
      params[i] <- exp(2 * (lgammas[3] + a1 + a2) * mm$bb[i])
    }

    lmats <- upperPadding(mm[jmm], relative = FALSE)

    val <- sparseMatrix(
      i = lmats$graph@i + 1L,
      j = lmats$graph@j + 1L,
      x = drop(lmats$xx %*% params)
    )

    return(forceSymmetric(val))

  }

#' Define the spacetime model matrices.
#'
#' This function computes all the matrices needed to build
#' the precision matrix for spatio-temporal model,
#' as in Lindgren et. al. (2023)
#'
#' @param tmesh a mesh object over the time domain.
#' @param smesh a mesh object over the spatial domain.
#' @param model a string identifying the model.
#' So far we have the following models:
#' '102', '121', '202' and '220' models.
#' @param constr logical to indicate if the integral of the field
#' over the domain is to be constrained to zero. Default value is FALSE.
#' @details
#' See the paper and the vignette('stspdemodels') for details.
#' @return a list containing needed objects for model definition.
#' 1. 'manifold' to spedify the a string with the model identification
#' 2. a length three vector with the constants `c1`, `c2` and `c3`
#' 3. the vector `d`
#' 4. the matrix `T`
#' 5. the model matrices `M_1`, ..., `M_m`
#' @export
stModel.matrices <-
  function(smesh, tmesh, model, constr = FALSE) {
    stopifnot(inherits(smesh, "inla.mesh"))
    stopifnot(inherits(tmesh, "inla.mesh.1d"))
    stopifnot(nchar(model) == 3)
    stopifnot(model %in% c("102", "121", "202", "220"))

    uM <- function(m) { ### extract the upper
      m <- INLA::inla.as.dgTMatrix(m)
      i.u <- which(m@i <= m@j)
      m@i <- m@i[i.u]
      m@j <- m@j[i.u]
      m@x <- m@x[i.u]
      return(m)
    }

    val <- list()

    if (tmesh$cyclic) {
      if (model == "102") {
        val$bb <- rep(c(1, 2, 1), 2)
        val$TT <- rbind(rep(2:0, 2), rep(0:1, each = 3))
      }
      if (model == "121") {
        val$bb <- c(1, 3, 3, 1, 1, 1)
        val$TT <- rbind(
          c(3:0, 1:0),
          rep(0:1, c(4, 2))
        )
      }
      if (model == "202") {
        val$bb <- rep(c(1, 2, 1), 3)
        val$TT <- rbind(rep(2:0, 3), rep(0:2, each = 3))
      }
      if (model == "220") {
        val$bb <- c(1, 4, 6, 4, 1, 1, 2, 1, 1)
        val$TT <- rbind(
          c(4:0, 2:0, 0),
          rep(0:2, c(5, 3, 1))
        )
      }
    } else {
      if (model == "102") {
        val$bb <- rep(c(1, 2, 1), 3)
        val$TT <- rbind(rep(2:0, 3), rep(0:2 / 2, each = 3))
      }
      if (model == "121") {
        val$bb <- c(1, 3, 3, 1, 1, 2, 1, 1, 1)
        val$TT <- rbind(
          c(3:0, 2:0, 1:0),
          rep(0:2 / 2, c(4, 3, 2))
        )
      }
      if (model == "202") {
        val$bb <- rep(c(1, 2, 1), 5)
        val$TT <- rbind(rep(2:0, 5), rep(0:4 / 2, each = 3))
      }
      if (model == "220") {
        val$bb <- c(1, 4, 6, 4, 1, 1, 3, 3, 1, 1, 2, 1, 1, 1, 1)
        val$TT <- rbind(
          c(4:0, 3:0, 2:0, 1:0, 0),
          rep(0:4 / 2, c(5, 4, 3, 2, 1))
        )
      }
    }

    tfe <- INLA::inla.mesh.fem(tmesh, order = 2)
    sfe <- INLA::inla.mesh.fem(smesh, order = 4)

    ns <- nrow(sfe$c0)
    nt <- nrow(tfe$c0)
    h <- mean(diff(tmesh$loc))

    if (model %in% c("102", "121")) {
      J0 <- tfe$c0
      if (!tmesh$cyclic) {
        J1 <- Matrix::sparseMatrix(
          i = c(1, nt), j = c(1, nt), x = c(0.5, 0.5)
        )
      }
      J2 <- tfe$g1

      val$M1 <- uM(kronecker(J0, sfe$c0))
      val$M2 <- uM(kronecker(J0, sfe$g1))
      val$M3 <- uM(kronecker(J0, sfe$g2))

      if (model == "102") {
        if (!tmesh$cyclic) {
          val$M4 <- uM(kronecker(J1, sfe$c0))
          val$M5 <- uM(kronecker(J1, sfe$g1))
          val$M6 <- uM(kronecker(J1, sfe$g2))
        }

        val$M7 <- uM(kronecker(J2, sfe$c0))
        val$M8 <- uM(kronecker(J2, sfe$g1))
        val$M9 <- uM(kronecker(J2, sfe$g2))
      } else {
        val$M4 <- uM(kronecker(J0, sfe$g3))

        if (!tmesh$cyclic) {
          val$M5 <- uM(kronecker(J1, sfe$c0))
          val$M6 <- uM(kronecker(J1, sfe$g1))
          val$M7 <- uM(kronecker(J1, sfe$g2))
        }

        val$M8 <- uM(kronecker(J2, sfe$c0))
        val$M9 <- uM(kronecker(J2, sfe$g1))
      }
    }

    if (model %in% c("202", "220")) {
      Jm <- Jmatrices(tmesh)

      val$M1 <- uM(kronecker(Jm$J0, sfe$c0))
      val$M2 <- uM(kronecker(Jm$J0, sfe$g1))
      val$M3 <- uM(kronecker(Jm$J0, sfe$g2))

      if (model == "202") {
        if (!tmesh$cyclic) {
          val$M4 <- uM(kronecker(Jm$J1, sfe$c0))
          val$M5 <- uM(kronecker(Jm$J1, sfe$g1))
          val$M6 <- uM(kronecker(Jm$J1, sfe$g2))
        }

        val$M7 <- uM(kronecker(Jm$J2, sfe$c0))
        val$M8 <- uM(kronecker(Jm$J2, sfe$g1))
        val$M9 <- uM(kronecker(Jm$J2, sfe$g2))

        if (!tmesh$cyclic) {
          val$M10 <- uM(kronecker(Jm$J3, sfe$c0))
          val$M11 <- uM(kronecker(Jm$J3, sfe$g1))
          val$M12 <- uM(kronecker(Jm$J3, sfe$g2))
        }

        val$M13 <- uM(kronecker(Jm$J4, sfe$c0))
        val$M14 <- uM(kronecker(Jm$J4, sfe$g1))
        val$M15 <- uM(kronecker(Jm$J4, sfe$g2))
      } else {
        val$M4 <- uM(kronecker(Jm$J0, sfe$g3))
        val$M5 <- uM(kronecker(Jm$J0, sfe$g4))

        if (!tmesh$cyclic) {
          val$M6 <- uM(kronecker(Jm$J1, sfe$c0))
          val$M7 <- uM(kronecker(Jm$J1, sfe$g1))
          val$M8 <- uM(kronecker(Jm$J1, sfe$g2))
          val$M9 <- uM(kronecker(Jm$J1, sfe$g3))
        }

        val$M10 <- uM(kronecker(Jm$J2, sfe$c0))
        val$M11 <- uM(kronecker(Jm$J2, sfe$g1))
        val$M12 <- uM(kronecker(Jm$J2, sfe$g2))

        if (!tmesh$cyclic) {
          val$M13 <- uM(kronecker(Jm$J3, sfe$c0))
          val$M14 <- uM(kronecker(Jm$J3, sfe$g1))
        }

        val$M15 <- uM(kronecker(Jm$J4, sfe$c0))
      }
    }

    stopifnot(length(val$bb) == ncol(val$TT))
    stopifnot(ncol(val$TT) == (length(val) - 2L))
    names(val)[3:length(val)] <- paste0("M", 1:ncol(val$TT))

    if (length(unique(diff(tmesh$loc))) > 1) {
      warning("Edge correction for irregular are not OK yet!")
    }

    if (constr) {
      val$extraconstr <- list(
        A = kronecker(tfe$c0@x, sfe$va[, 1]), e = 0
      )
      val$extraconstr$A <-
        matrix(val$extraconstr$A / sum(val$extraconstr$A), 1)
    }

    return(val)
  }

#' The 2nd order temporal matrices with boundary correction
#' @param tmesh Temporal mesh
#' @return return a list of temporal finite element method matrices
#' for the supplied mesh.
#' @export
Jmatrices <- function(tmesh) {
  tfe <- INLA::inla.mesh.fem(tmesh, 2L)
  nt <- nrow(tfe$g1)
  h <- mean(diff(tmesh$loc))

  J0 <- tfe$c0
  J0[2, 2] <- tfe$c0[2, 2] / 2
  J0[nt - 1, nt - 1] <- tfe$c0[nt - 1, nt - 1] / 2

  if (tmesh$cyclic) {
    J1 <- NULL
  } else {
    J1 <- Matrix::sparseMatrix(
      i = c(1, 1, 2, nt - 1, nt - 1, nt),
      j = c(1, 2, 2, nt - 1, nt, nt),
      x = c(5, -1, 5, 5, -1, 5) / 4
    )
  }

  J2 <- tfe$g1 * 2
  J2[2, 2] <- tfe$g1[2, 2]
  J2[1, 2] <- J2[2, 1] <- tfe$g1[1, 2]
  J2[nt, nt - 1] <- J2[nt - 1, nt] <- tfe$g1[nt, nt - 1]
  J2[nt - 1, nt - 1] <- tfe$g1[nt - 1, nt - 1]

  if (tmesh$cyclic) {
    J3 <- NULL
  } else {
    J3 <- Matrix::sparseMatrix(
      i = c(1, 1, 2, nt - 1, nt - 1, nt),
      j = c(1, 2, 2, nt - 1, nt, nt),
      x = c(2, -2, 2, 2, -2, 2) / (h^2)
    )
  }

  J4 <- tfe$g2
  J4[1, 1] <- tfe$g2[1, 1] / 3
  J4[nt, nt] <- tfe$g2[nt, nt] / 3
  J4[1, 2] <- J4[2, 1] <- tfe$g2[1, 2] / 2
  J4[nt - 1, nt] <- J4[nt, nt - 1] <- tfe$g2[nt - 1, nt] / 2
  J4[2, 2] <- tfe$g2[2, 2] * 5 / 7
  J4[nt - 1, nt - 1] <- tfe$g2[nt - 1, nt - 1] * 5 / 7

  return(list(J0 = J0, J1 = J1, J2 = J2, J3 = J3, J4 = J4))
}
