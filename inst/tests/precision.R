library(INLA)
library(INLAspacetime)

rxy <- c(10, 5) ## size of spatial domain
nt <- 10 ## number of time points

(r0 <- mean(rxy))

## setup rectangle for spatial domain
bb <- rbind(
    x = c(0, rxy[1]),
    y = c(0, rxy[2]))

domain <- cbind(
    x = bb[c(1, 3, 3, 1, 1)],
    y = bb[c(2, 2, 4, 4, 2)])

## spatial mesh
smesh <- inla.mesh.2d(###fm_mesh_2d(
    loc.domain = domain,
    offset = r0 / c(5, 2), ##r0 / c(40, 3),
    max.edge = r0 / c(5, 2), ##c(r0 / c(20, 5),
    cutoff = r0 / 10)##r0 / 40)
smesh$n

if(FALSE)
    plot(smesh)

## temporal mesh
tmesh <- inla.mesh.1d(
    loc = 1:nt)

## model parameters
params <- c(
    rs = r0 / 3, ## spatial range
    rt = nt / 2, ## temporal range
    sigma.u = 1) ## standard deviation
params

system.time(
    stPrecision <- stModel.precision(
        smesh = smesh,
        tmesh = tmesh,
        model = "121",
        theta = log(params)
    )
)

image(stPrecision)

system.time(
    cmodel <- stModel.define(
        smesh = smesh,
        tmesh = tmesh,
        model = "121",
        control.priors = list(
            prt = c(1, 0.1),
            prs = c(1, 0.1),
            psigma = c(1, 0.1)
        ),
        verbose = TRUE,
        debug = !TRUE,
        useINLAprecomp = FALSE
    )
)

library(corGraphs)

str(cmodel$f)

system.time(
    Qij <- cgeneric_get(cmodel, cmd = "graph", theta = log(params))
)

system.time(
    qq <- cgeneric_get(cmodel, theta = log(params), cmd = "Q")
)

str(Qij)
str(qq)

sapply(Qij, range)
str(Qij)

Q <- sparseMatrix(
    i = Qij[[1]] + 1L,
    j = Qij[[2]] + 1L,
    x = qq
)

str(stPrecision)
str(Q)

all.equal(stPrecision@x,
          Q@x)

plot(stPrecision@x,
          Q@x)

n <- 100
xyt <- cbind(runif(n, bb[1, 1], bb[1, 2]),
             runif(n, bb[2, 1], bb[2, 2]),
             sort(sample(1:nt, n, replace = TRUE)))

A <- inla.spde.make.A(
    mesh = smesh,
    loc = xyt[, 1:2], ### spatial locations
    group = xyt[, 3], ### temporal locations
    group.mesh = tmesh)

image(A)
