
library(sf)
library(fmesher)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(sf)

rxy <- c(7, 5) ## size of spatial domain
(r0 <- mean(rxy))
tlim <- c(0, 20)

## setup rectangle for spatial domain
bb <- rbind(
    x = c(0, rxy[1]),
    y = c(0, rxy[2]))
domain <- cbind(
    x = bb[c(1, 3, 3, 1, 1)],
    y = bb[c(2, 2, 4, 4, 2)])

## data spacetime locations
set.seed(1)
nd <- 3000
dataf <- data.frame(
    xloc = runif(nd, bb[1, 1], bb[1, 2]),
    yloc = runif(nd, bb[2, 1], bb[2, 1]),
    tloc = sort(runif(nd, tlim[1], tlim[2])))
print(t(sapply(dataf, summary)))

## model parameters
params <- c(
    rs = r0 / 2, ## spatial range
    rt = diff(tlim) / 1, ## temporal range
    sigma.u = 1) ## standard deviation
params

## noise parameter
(sigma.e <- 1/sqrt(exp(10))) ## will be fixed

nt <- 11    ## temporal resolution
xr <- r0/12 ## spatial resolution

## spatial mesh
smesh <- fm_mesh_2d(
    loc = fm_hexagon_lattice(
        bnd = st_buffer(
            sf::st_sfc(sf::st_polygon(list(domain))),
            dist = xr*2),
        edge_len = xr), 
    max.edge = xr * 3,
    offset = xr * 5,
    cutoff = xr / 2)

cat("Number of spatial mesh nodes:", smesh$n, "\n")

if(FALSE) {
    
    plot(smesh)
    lines(domain, col = 4)
    
}

## temporal mesh
tmesh <- fm_mesh_1d(
    loc = seq(tlim[1], tlim[2], length = nt))

cat("Number of time points:", nt, "\n")

## build the precision matrix
t0 <- Sys.time()
qq <- stModel.precision(smesh, tmesh, '220', log(params))
cat("Precision built\n")
print(Sys.time()-t0)

if(FALSE)
    image(qq)

## sample
set.seed(2)
zz <- rnorm(nrow(qq))
t0 <- Sys.time()
xx <- inla.qsolve(qq, matrix(zz, ncol = 1))[,1]
cat("Field simulation finnished\n")
print(Sys.time()-t0)
print(summary(xx))

set.seed(3)
error <- rnorm(nd, 0, sigma.e)
print(summary(error))

## project the sample to spacetime data locations
A.d <- inla.spde.make.A(
    mesh = smesh,
    loc = cbind(dataf$xloc, dataf$yloc),
    group = dataf$tloc,
    group.mesh = tmesh)
dataf$y <- drop(A.d %*% xx + error)

cat("Outcome summary:\n")
print(summary(dataf$y))

############################################################
## model fit with inlabru

## build 'cgeneric' model
t0 <- Sys.time()
stmodel <- stModel.define(
    smesh = smesh,
    tmesh = tmesh,
    model = "121",
    control.priors = list(
        prs = c(r0 / 3, 0.05),
        prt = c(nt / 5, 0.05),
        psigma = c(2, 0.05)))
cat("Model built finnish\n")
print(Sys.time()-t0)

## likelihood setup
mlike <- bru_obs(
    y ~ .,
    data = dataf,
    control.family = list(
        hyper = list(
            prec = list(initial = 10, fixed = TRUE)
        )
    )
)

## linear predictor definition
mcomps <- ~ Intercept(1) +
    spacetime(list(space = cbind(xloc, yloc),
                   time = tloc),
              model = stmodel)

## modelfit
t0 <- Sys.time()
result <- bru(
    mcomps,
    mlike,
    options = list(
        verbose = !TRUE)
)
cat("Model fit finished\n")
print(Sys.time()-t0)
print(result$cpu.used)

## compare with truth
print(cbind(true = log(params), 
            result$summary.hyper[, c(1, 2)]))

cor(result$summary.random$spacetime$mean, xx)
sd(result$summary.random$spacetime$mean)
