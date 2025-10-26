
library(ggplot2)
library(sf)
library(INLA)
library(INLAspacetime)
stopifnot(packageVersion("INLAspacetime")>'0.1.9.2')

mesh <- fm_mesh_2d(
    loc = cbind(0, 0),
    offset = c(7,3),
    max.edge = c(5, 10)/10,
    n = 4
)

mesh$n

plot(mesh)

## stationary model
smodel <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(1, 0.1),
    prior.sigma = c(1, 0.5)
)

## Non-stationary model

## triangle centers
ce.tri <- cbind(   
  mesh$loc[mesh$graph$tv[,1], 1:2] +
  mesh$loc[mesh$graph$tv[,2], 1:2] +
  mesh$loc[mesh$graph$tv[,3], 1:2])/3

## barriers
barrs <- list(
    cbind(c(0, 10, 10, 0, 0),
          c(-1, -1, 1, 1, -1)),
    cbind(c(-2, 0, 0, -8, -11, -2),
          c(0, 0, 1, 10, 10, 0)),
    cbind(c(0, 0, -2, -11, -8, 0),
          c(-1, 0, 0, -10, -10, -1))
    )

barrier1 <- st_sfc(
    st_multipolygon(
        list(st_polygon(barrs[1]))))
barrier2 <- st_sfc(
    st_multipolygon(
        list(st_polygon(barrs[2]))))
barrier3 <- st_sfc(
    st_multipolygon(
        list(st_polygon(barrs[3]))))

## triangles in the barrier
tri.ids <- list(
    unlist(
        fmesher::fm_contains(
                     x = barrier1,
                     y = mesh,
                     type = 'centroid'
                 )
    ),
    unlist(
        fmesher::fm_contains(
                     x = barrier2,
                     y = mesh,
                     type = 'centroid'
                 )
    ),
    unlist(
        fmesher::fm_contains(
                     x = barrier3,
                     y = mesh,
                     type = 'centroid'
                 )
    )
)

ggplot() + theme_minimal() +
    geom_sf(data = barrier1, fill = rgb(1,.5,.5,.5)) +
    geom_sf(data = barrier2, fill = rgb(.5,1,.5,.5)) +
    geom_sf(data = barrier3, fill = rgb(.5,.5,1,.5)) 

plot(mesh)
for(i in 1:length(tri.ids)) {
    points(ce.tri[tri.ids[[i]], ], pch = 19, col = i)
}

### define the cgeneric barrier model
bmodel <- barrierModel.define(
    mesh = mesh,
    barrier.triangles = tri.ids,
    prior.range = c(1, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = c(1, 1, 1) ## as stationary
    ##, useINLAprecomp = FALSE
)

## model fit with no data (posterior = prior)
sfit0 <- inla(
    y ~ 0 + f(i, model = smodel),
    data = data.frame(y = NA, i = 1:mesh$n),
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)
bfit0 <- inla(
    y ~ 0 + f(i, model = bmodel),
    data = data.frame(y = NA, i = 1:mesh$n),
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)

rbind(sfit0$mode$theta,
      bfit0$mode$theta)

spmargs0 <- lapply(sfit0$internal.marginals.hyperpar, function(m)
    inla.tmarginal(exp, m))
bpmargs0 <- lapply(bfit0$internal.marginals.hyperpar, function(m)
    inla.tmarginal(exp, m))

## visualize
par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1,0), las = 1)
for(i in 1:2) {
    plot(spmargs0[[i]], type = "l",
         xlab = c("range", "sigma")[i], ylab = "Density")
    lines(bpmargs0[[i]], lty = 2)
}

## simulate data and fit
## set parameter values
range <- 5
sigma <- 2

Qs <- inla.spde2.precision(smodel, theta = log(c(range, sigma)))
xx <- inla.qsample(n = 1, Q = Qs)

## model fit
sfit <- inla(
    y ~ 0 + f(i, model = smodel),
    data = data.frame(y = xx[,1], i = 1:mesh$n),
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)
bfit <- inla(
    y ~ 0 + f(i, model = bmodel),
    data = data.frame(y = xx[,1], i = 1:mesh$n),
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)

rbind(sfit$mode$theta,
      bfit$mode$theta)

spmargs <- lapply(sfit$internal.marginals.hyperpar, function(m)
                 inla.tmarginal(exp, m))
bpmargs <- lapply(bfit$internal.marginals.hyperpar, function(m)
                 inla.tmarginal(exp, m))

## visualize
par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1,0), las = 1)
for(i in 1:2) {
    plot(spmargs[[i]], type = "l",
         xlab = c("range", "sigma")[i], ylab = "Density")
    lines(bpmargs[[i]], lty = 2)
}

