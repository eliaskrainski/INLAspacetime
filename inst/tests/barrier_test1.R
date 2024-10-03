
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

## Stationary FEM
system.time(sfem <- fm_fem(mesh))

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

### barrier FEM
system.time(bfem <- mesh2fem.barrier(mesh, tri.ids))

## check the structure matrices
all.equal(sfem$g1,
          Reduce("+", bfem$D))

all.equal(sfem$c1, bfem$I)

all.equal(sfem$c0@x,
          Reduce("+", bfem$C))

### define the cgeneric barrier model
bmodel <- barrierModel.define(
    mesh = mesh,
    barrier.triangles = tri.ids,
    prior.range = c(1, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = c(1, 1, 1), ## as stationary
    useINLAprecomp = FALSE
)

## model parameters
range <- runif(1, 3, 10)
sigma <- runif(1, .3, 3)

## get the precision
ifit <- inla(
    y ~ 0 + f(i, model = bmodel),
    data = data.frame(y = NA, i = 1:mesh$n),
    control.mode = list(
        theta = c(10, log(c(range, sigma))),
        fixed = TRUE),
    control.compute = list(config = TRUE)
)

## the upper part of the prior
Qub <- inla.as.sparse(
    ifit$misc$configs$config[[1]]$Qprior
)

## the non-stationary precision matrix
Qb <- inla.as.sparse(
    sparseMatrix(
        i = Qub@i + 1L, 
        j = Qub@j + 1L,
        x = Qub@x,     useINLAprecomp = FALSE
        symmetric = TRUE,
        repr = "T"
    )
)
Vb <- inla.qinv(Qb)

## The stationary precision matrix
range2 <- range^2
sigma2 <- sigma^2
k2 <- 8/(range2)
t2 <- 1/(4*pi*k2*sigma2)
Qs <- inla.as.sparse(
    t2 * ((k2^2)*sfem$c1 + ## C1: closer comparison with Hakkon's
          2*k2*sfem$g1 + sfem$g2))
Vs <- inla.qinv(Qs)

## compare
summary(diag(Vs))
summary(diag(Vb))

summary(Qs@x)
summary(Qb@x)

sum(diag(Vb)) / sum(diag(Vs))

