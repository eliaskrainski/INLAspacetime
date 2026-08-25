## multiple barriers and multiple domain

library(ggplot2)
library(sf)
library(fmesher)
library(INLA)
library(INLAspacetime)
stopifnot(packageVersion("INLAspacetime")>'0.1.14.902')

## Define two domains
s <- 10
bnds <- list(
    bnd1 = st_sfc(list(st_polygon(
        list(cbind(c(1, 3, 3, 1, 1),
                   c(0, 0, 3, 3, 0)) * s
             )))),
    bnd2 = st_sfc(list(st_polygon(
        list(cbind(c(-3, -1, -1, -3, -3),
                   c(0, 0, 3, 3, 0)) * s
             ))))
)

plot(st_union(bnds[[1]], bnds[[2]]))

## build a mesh over each one and collect 
mesh <- fm_collect(lapply(bnds, function(bnd)
    fm_mesh_2d(
        loc = fm_hexagon_lattice(
            st_buffer(bnd, s/10),
            edge_len = s/20),
        offset = s/3, 
        max.edge = s/5)
    ))

str(mesh,2)

sapply(mesh[[1]], function(x) x$n)

bbm <- fm_bbox(mesh)

par(mfrow=c(1,1), mar = c(0,0,0,0))
plot(mesh[[1]]$bnd1,
     xlim = range(sapply(bbm, function(x) range(x[[1]]))))
plot(mesh[[1]]$bnd2, add = TRUE)

## Stationary FEM
system.time(sfem <- fm_fem(mesh))

## Barrier definition
## Set of barriers poly 
barrs <- lapply(
    list(
        st_polygon(list(
            cbind(c(-5,-5,  5, 5, -5),
                  c( 5,5.5, -8,-8.5, 5)) * s)),
        st_polygon(list(
            cbind(c(-5,-5,  5, 5, -5),
                  c(-8.5,-8,5.5,5,-8.5)) * s))
        ), function(p)
            st_sfc(st_multipolygon(list(p))))

ggplot() + theme_minimal() +
    geom_sf(data = barrs[[1]], fill = rgb(1,.5,.1,.5)) +
    geom_sf(data = barrs[[2]], fill = rgb(.1,.5,1,.5)) 

## triangles in the barrier
tri.ids <- barrier_mesh_centroids(mesh, barrs)

## triangle centorids (to be visualized)
ce.tri <- lapply(mesh$fun_spaces, fm_centroids)
str(ce.tri)

par(mfrow=c(1,1), mar = c(0,0,0,0))
plot(mesh[[1]]$bnd1,
     xlim = range(sapply(bbm, function(x) range(x[[1]]))))
plot(mesh[[1]]$bnd2, add = TRUE)
for(d in 1:length(tri.ids)) {
    for(b in 1:length(tri.ids[[d]])) {
        points(ce.tri[[d]][tri.ids[[d]][[b]], ], pch = 19, col = b)
    }
}

### barrier FEM for a fm_collect() object
system.time(bfem <- collect2fem.barrier(mesh, tri.ids))

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
    prior.range = c(s/5, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = c(1, 1)
)

## model parameters
range <- s
sigma <- 2

## get the precision
ifit <- inla(
    y ~ 0 + f(i, model = bmodel),
    verbose = !TRUE,
    data = data.frame(y = NA, i = 1:bmodel$f$n), 
    control.mode = list(
        theta = c(10, log(c(range, sigma))),
        fixed = TRUE),
    control.compute = list(config = TRUE)
)

summary(ifit$summary.random$i$sd)

## the upper part of the prior
Qub <- inla.as.sparse(
    ifit$misc$configs$config[[1]]$Qprior
)

## the non-stationary precision matrix
Qb <- inla.as.sparse(
    sparseMatrix(
        i = Qub@i + 1L, 
        j = Qub@j + 1L,
        x = Qub@x,
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
    t2 * ((k2^2)*sfem$c1 + ## C1: closer comparison with Hakkon's choice
          2*k2*sfem$g1 + sfem$g2))
Vs <- inla.qinv(Qs)

## compare
summary(diag(Vs))
summary(diag(Vb))

summary(Qs@x)
summary(Qb@x)

sum(diag(Vb)) / sum(diag(Vs))

