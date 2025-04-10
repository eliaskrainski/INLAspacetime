
library(INLA)
library(INLAspacetime)
library(inlabru)
library(sf)

## some INLA setup
if(Sys.info()["user"]=="eliask") {
    inla.setOption(
        inla.call = "remote",
        num.threads = "8:1",
##        pardiso.license = "~/.pardiso.lic",
        safe = FALSE
        )
}

rxy <- c(10, 5) ## size of spatial domain
nt <- 30 ## number of time points
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
    offset = r0 / c(40, 3),
    max.edge = r0 / c(20, 5),
    cutoff = r0 / 40)
smesh$n

if(FALSE)
    plot(smesh)

## temporal mesh
tmesh <- fm_mesh_1d(
    loc = 1:nt)

## model parameters
params <- c(
    rs = r0 / 3, ## spatial range
    rt = nt / 3, ## temporal range
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

## collect params from the defined model
.gammas <- params2gammas(log(params), 1, 2, 1)
.b <- stmodel$f$cgeneric$data$doubles$bb
.tt <- matrix(stmodel$f$cgeneric$data$matrices$tt[-(1:2)], 2)
.dxx <- stmodel$f$cgeneric$data$matrices$xx
.d <- .b * exp(2 * (.gammas[1] * .tt[1, ] + .gammas[2] * .tt[2, ]))

## build precision
qq <- sparseMatrix(
    i = stmodel$f$cgeneric$data$ints$ii + 1L,
    j = stmodel$f$cgeneric$data$ints$jj + 1L,
    x = drop(t(matrix(.dxx[-c(1,2)], .dxx[1], byrow = TRUE)) %*%
             .d) * exp(2 * .gammas[3]),
    symmetric = TRUE)
rm(.gammas, .b, .tt, .dxx, .d)

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
##        control.mode = list(
  ##          theta = log(c(2 / sigma.e^2, params / 2)),
    ##        restart = TRUE),
        verbose = TRUE)
)

## compare with truth
cbind(true = c(1 / sigma.e^2, log(params)),
      result$summary.hyper[, c(1, 2)])
