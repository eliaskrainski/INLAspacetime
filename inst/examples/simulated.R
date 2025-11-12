
library(INLA)
library(INLAspacetime)
library(inlabru)
library(sf)

rxy <- c(7, 5) ## size of spatial domain
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
smesh <- fm_mesh_2d(
    loc.domain = domain,
    offset = r0 / c(50, 3),
    max.edge = r0 / c(20, 5),
    cutoff = r0 / 20)
smesh$n

if(FALSE)
    plot(smesh)

## temporal mesh
tmesh <- fm_mesh_1d(
    loc = 1:nt)

## model parameters
params <- c(
    rs = r0 / 3, ## spatial range
    rt = nt / 2, ## temporal range
    sigma.u = 1) ## standard deviation
params

## build 'cgeneric' model
stmodel <- stModel.define(
    smesh = smesh,
    tmesh = tmesh,
    model = "121",
    control.priors = list(
        prs = c(r0 / 10, 0.05),
        prt = c(nt / 10, 0.05),
        psigma = c(2, 0.05)))

## build the precision matrix
qq <- stModel.precision(smesh, tmesh, '121', log(params))

if(FALSE)
    image(qq)

## sample
xx <- inla.qsample(n = 1, Q = qq)

ns <- 1000 ## number locations in space
(nd <- nt * ns) ## number of observations

(sigma.e <- 1/sqrt(9))
error <- rnorm(nd, 0, sigma.e)

## data spacetime locations
dataf <- data.frame(
    xloc = runif(nd, bb[1, 1], bb[1, 2]),
    yloc = runif(nd, bb[2, 1], bb[2, 1]),
    tloc = sort(sample(1:nt, nd, replace = TRUE)))

## project the sample to spacetime data locations
A.d <- inla.spde.make.A(
    mesh = smesh,
    loc = cbind(dataf$xloc, dataf$yloc),
    group = dataf$tloc,
    group.mesh = tmesh)
dataf$y <- drop(A.d %*% xx + error)

############################################################
## model fit with inlabru

## likelihood setup
mlike <- like(
    y ~ .,
    data = dataf,
    control.family = list(
        hyper = list(
            prec = list(
                prior = "pc.prec",
                param = c(0.5, 0.05)))))

## linear predictor definition
mcomps <- ~ Intercept(1) +
    spacetime(list(space = cbind(xloc, yloc),
                   time = tloc),
              model = stmodel)

## modelfit
result <- bru(
    mcomps,
    mlike,
    options = list(
        num.threads = '5',
##        control.mode = list(
  ##          theta = log(c(2 / sigma.e^2, params / 2)),
    ##        restart = TRUE),
        verbose = TRUE)
)
result$cpu.used

## compare with truth
cbind(true = c(1 / sigma.e^2, log(params)),
      result$summary.hyper[, c(1, 2)])
