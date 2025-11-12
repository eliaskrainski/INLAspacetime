
library(ggplot2)
library(sf)
library(INLA)
library(INLAspacetime)
stopifnot(packageVersion("INLAspacetime")>'0.1.9.2')

s <- 10
bnd <- st_sfc(list(st_polygon(
    list(cbind(c(-1, 1, 1, -1, -1),
               c(-1, -1, 1, 1, -1)) * s
         ))))

mesh <- fm_mesh_2d(
    loc = fm_hexagon_lattice(
        st_buffer(bnd, s/10),
        edge_len = s/20),
    offset = s/3, 
    max.edge = s/5
)

mesh$n

plot(mesh)

## stationary model
smodel <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(s/5, 0.1),
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
    cbind(c(0, 2, 2, 0, 0) * s,
          c(-1, -1, 1, 1, -1) * s/10),
    cbind(c(0, 0, -s/5, -1.8*s, -1.5*s, 0),
          c(-s/10, 0, 0, -1.5*s, -1.5*s, -s/10)),
    cbind(c(-s/5, 0, 0, -1.5*s, -1.8*s, -s/5),
          c(0, 0, s/10, 1.5*s, 1.5*s, 0))
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
rfracs <- c(1,1,1)##c(.1, .5, .7)
bmodel <- barrierModel.define(
    mesh = mesh,
    barrier.triangles = tri.ids,
    prior.range = c(s/5, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = rfracs
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
      bfit0$mode$theta) ## equal if all(rfracs == 1)

spmargs0 <- lapply(sfit0$internal.marginals.hyperpar, function(m)
    inla.tmarginal(exp, m))
bpmargs0 <- lapply(bfit0$internal.marginals.hyperpar, function(m)
    inla.tmarginal(exp, m))

## visualize (equal if all(rfracs == 1))
par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1,0), las = 1, bty = "n")
for(i in 1:2) {
    plot(bpmargs0[[i]], type = "l",
         xlim = range(bpmargs0[[i]][, 1],
                      spmargs0[[i]][, 1]),
         ylim = range(bpmargs0[[i]][, 2],
                      spmargs0[[i]][, 2]),
         xlab = c("range", "sigma")[i], ylab = "Density", log = 'x')
    lines(spmargs0[[i]], lty = 2)
}
legend("top", bty = "n", lty = 1:2,
       c("barrier", "stationary"))

## simulate data and fit
## set parameter values
range <- s
sigma <- 2

## barrier FEM
system.time(bfem <- mesh2fem.barrier(mesh, tri.ids))

## barrier model precision
Qb <- inla.barrier.q(bfem, ranges = range * c(1, rfracs), sigma = sigma)
xx <- inla.qsample(n = 1, Q = Qb)

all.equal(Qb, inla.spde2.precision(smodel, theta = log(c(range, sigma))))

ylm <- xlm <- c(-1, 1) * s
gproj <- fm_evaluator(mesh = mesh, xlim = xlm, ylim = ylm)
z <- fm_evaluate(gproj, field = xx[, 1])

par(mfrow = c(1, 1), mar = c(0,0,0,0))
image(gproj$x, gproj$y, z)

ggplot() + theme_minimal() +
    geom_tile(aes(x = gproj$lattice$loc[, 1],
                  y = gproj$lattice$loc[, 2],
                  fill = as.vector(z))) +
    geom_sf(data = barrier1, fill = "transparent") + 
    geom_sf(data = barrier2, fill = "transparent") + 
    geom_sf(data = barrier3, fill = "transparent") + 
    xlab("") + ylab("") + 
    xlim(xlm) + ylim(ylm) + 
    scale_fill_distiller(palette = "RdBu")
rfracs

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
par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1,0), las = 1, bty = "n")
for(i in 1:2) {
    plot(bpmargs[[i]], type = "l",
         xlim = range(bpmargs[[i]][, 1],
                      spmargs[[i]][, 1]),
         ylim = range(bpmargs[[i]][, 2],
                      spmargs[[i]][, 2]),
         xlab = c("range", "sigma")[i], ylab = "Density", log = 'x')
    lines(spmargs[[i]], lty = 2)
    abline(v = c(range, sigma)[i], lwd = 2, lty = 3)
}
legend("top", bty = "n", lty = 1:2,
       c("barrier", "stationary"))

