## multiple barriers and multiple domain

library(ggplot2)
library(sf)
library(fmesher)
library(INLA)
library(INLAspacetime)
stopifnot(packageVersion("INLAspacetime")>='0.1.14.905')

## Define two domains
plfn <- function(a=0, b=0, s=1, r=0.7) {
    return(cbind(c(-1, 1, 1, -1, -1)*r + a,
                 c(-1, -1, 1, 1, -1) + b)*s)
}

s <- 50
bnds <- list(
    bnd1 = st_sfc(list(st_polygon(list(plfn(-1,1,s))))),
    bnd2 = st_sfc(list(st_polygon(list(plfn(1,1,s)))))
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
barrs <- list(
    bnd1 = st_sfc(list(st_polygon(list(plfn(-4,5,s/4,1))))),
    bnd2 = st_sfc(list(st_polygon(list(plfn(4,5,s/4,1))))),
    bnd3 = st_sfc(list(st_polygon(list(
        cbind(c(seq(-1, 1, 0.1), seq(1,-1,-0.1), -1)*1.5,
              c(1.4-cos(seq(-1,1,0.1)),
                1.2-cos(seq(1,-1,-0.1)),
                1.4-cos(-1)))*s))))
)

ggplot() + theme_minimal() +
    geom_sf(data = bnds[[1]], fill = rgb(1,.7,.5)) +
    geom_sf(data = bnds[[2]], fill = rgb(.5,.7,1)) +
    geom_sf(data = barrs[[1]], fill = rgb(1,.3,.1,.5)) +
    geom_sf(data = barrs[[2]], fill = rgb(.1,.3,1,.5)) +
    geom_sf(data = barrs[[3]], fill = rgb(0.5,1,0.5,.5)) 

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
bmodel_c <- barrierModel.define(
    mesh = mesh,
    barrier.triangles = tri.ids,
    prior.range = c(s/5, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = rep(1, length(barrs)),
    constr = TRUE ## one per domain
)

bmodel_c
sapply(bfem$C, sum)
rowSums(bmodel_c$f$extraconstr$A)

## model parameters
range <- s
sigma <- 2

## get the precision
fit_c <- inla(
    y ~ 0 + f(i, model = bmodel_c),
    verbose = !TRUE,
    data = data.frame(y = NA, i = 1:bmodel_c$f$n), 
    control.mode = list(
        theta = c(10, log(c(range, sigma))),
        fixed = TRUE),
    control.compute = list(config = TRUE)
)

## now without constraints, for comparison
bmodel <- barrierModel.define(
    mesh = mesh,
    barrier.triangles = tri.ids,
    prior.range = c(s/5, 0.1),
    prior.sigma = c(1, 0.5),
    range.fraction = rep(1, length(barrs))
)

fit <- inla(
    y ~ 0 + f(i, model = bmodel),
    verbose = !TRUE,
    data = data.frame(y = NA, i = 1:bmodel$f$n), 
    control.mode = list(
        theta = c(10, log(c(range, sigma))),
        fixed = TRUE),
    control.compute = list(config = TRUE)
)

summary(fit_c$summary.random$i$sd)
summary(fit$summary.random$i$sd)

## the upper part of the prior
Qub <- inla.as.sparse(
    fit$misc$configs$config[[1]]$Qprior
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
    t2 * ((k2^2)*sfem$c1 + ## C1: closer comparison 
          2*k2*sfem$g1 + sfem$g2))
Vs <- inla.qinv(Qs)

## compare
summary(diag(Vs))
summary(diag(Vb))

summary(Qs@x)
summary(Qb@x)

sum(diag(Vb)) / sum(diag(Vs))

