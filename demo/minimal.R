library(fmesher)
library(INLAspacetime)

## sample data
set.seed(1)
n <- 5
dataf <- data.frame(
  s1   = runif(n, -1, 1),
  s2   = runif(n, -1, 1),
  time = runif(n, 1, 4),
  y    = rnorm(n, 0, 1))

## mesh over space and time
smesh <- fm_mesh_2d(
  loc = cbind(0,0),
  max.edge = 5,
  offset = 2)
tmesh <- fm_mesh_1d(
  loc = 0:5)

if(require(INLA)) {

  ## define the cgeneric model
  stmodel <- stModel.define(
    smesh = smesh, ## spatial mesh
    tmesh = tmesh, ## temporal mesh
    model = '121', ## model, see the paper
    control.priors = list(
      prs = c(1, 0.1), ## P(spatial range < 1) = 0.1
      prt = c(5, 0), ## temporal range fixed to 5
      psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
    )
  )

  ## projector/mapper matrix
  Aproj <- inla.spde.make.A(
    mesh = smesh,
    loc = cbind(dataf$s1, dataf$s2),
    group = dataf$time,
    group.mesh = tmesh
  )

  ## control for the likelihood parameter (fixed)
  ctrl.lik <- list(
    hyper = list(
      prec = list(
        initial = 10,
        fixed = TRUE)
    )
  )

  ## model definition
  fmodel <- y ~ f(st, model = stmodel, A.local = Aproj)
  dataf$st <- NA

  ## model fitting
  fit <- inla(
    formula = fmodel,
    data = dataf,
    control.family = ctrl.lik)

  ## summary for the fixed effects
  fit$summary.fixed

  ## summary of the internal parameters
  fit$summary.hyperpar

}

