### fit one of the spacetime models (Lindgren et. al. 2022)
### for the piemonte dataset (Cameleti et. al. 2012)

### packages
library(INLA)
library(INLAspacetime)
library(inlabru)

### overall INLA setup
inla.setOption(
    inla.mode='compact',
    smtp='pardiso',
    pardiso.license='~/.pardiso.lic')

ctri <- list(
    int.strategy='ccd',
    parallel.linesearch=TRUE)

ctrc <- list(
    waic=TRUE,
    dic=TRUE,
    cpo=TRUE)

### the data filenames
u0 <- paste0(
    'http://inla.r-inla-download.org/',
    'r-inla.org/case-studies/Cameletti2012/')
coofl <- 'coordinates.csv'
datafl <- 'Piemonte_data_byday.csv'
bordersfl <- 'Piemonte_borders.csv'

### get the domain borders
if(!file.exists(bordersfl))
    download.file(paste0(u0, bordersfl), bordersfl)
dim(pborders <- read.csv(bordersfl))

### get the coordinates
if(!file.exists(coofl))
    download.file(paste0(u0, coofl), coofl)
dim(locs <- read.csv(coofl))

### get the dataset
if(!file.exists(datafl))
    download.file(paste0(u0, datafl), datafl)
dim(pdata <- read.csv(datafl))

head(pdata)

### prepare and select time
range(pdata$Date <- as.Date(pdata$Date, '%d/%m/%y'))
pdata$time <- as.integer(difftime(
    pdata$Date, min(pdata$Date), units='days'))+1

### define a temporal mesh
nt <- max(pdata$time)
tmesh <- fm_mesh_1d(1:nt, degree=1)
fm_dof(tmesh)

### mesh as in Cameleti et. al. 2012
smesh <- fm_mesh_2d(
    cbind(locs[,2], locs[,3]),
    loc.domain=pborders,
    max.edge=c(50, 300),
    offset=c(10, 140),
    cutoff=5,
    min.angle=c(26, 21))
fm_dof(smesh)

### visualize
par(mfrow=c(1,1), mar=c(0,0,1,0))
plot(smesh, asp=1)
points(locs[,2:3], pch=8, col='red')
lines(pborders, lwd=2, col='green4')

### prepare the covariates
xnames <- c('A', ###'UTMX', 'UTMY',
            'WS', 'TEMP', 'HMIX', 'PREC', 'EMI')
xmean <- colMeans(pdata[, xnames])
xsd <- sapply(pdata[xnames], sd)

### prepare the data (st loc, scale covariates and log PM10)
dataf <- data.frame(pdata[c('UTMX', 'UTMY', 'time')],
                    scale(pdata[xnames], xmean, xsd),
                    y=log(pdata$PM10))
str(dataf)

### likelihood precision prior
lkprec <- list(
    prec=list(prior='pcprec', param=c(1, 0.1)))

### likelihood data object
lhood <- like(
    formula = y ~ .,
    family="gaussian",
    control.family = list(
        hyper = lkprec),
    data=dataf)


### define the data Model
M <- ~ -1 + Intercept(1) + A + WS + TEMP + HMIX + PREC + EMI +
    field(list(space = cbind(UTMX, UTMY), time=time),
          model=stmodel)

### initial theta (higher prec and ranges and lower sigma)
theta.ini <- list('102'=c(4, 5.5, 5, 0),
                  '121'=c(4, 7.5, 12, 2))
theta.ini

### fit two first order in time models (separable and non-separable)
models <- c('102', '121')
results <- vector('list', 2)
names(results) <- models

for(m in 1:2) {

### define the spacetime model
    stmodel <- stModel.define(
        smesh, tmesh, models[m],
        control.priors=list(
            prs=c(70, 0.5),
            prt=c(50, 0.5),
            psigma=c(20, 0.5)),
        constr = TRUE) ## no need but helps

### fit
    results[[m]] <-
        bru(M,
            lhood,
            options = list(
                verbose=TRUE,
                control.mode=list(
                    theta=theta.ini[[m]],
                    restart=TRUE),
                control.inla=ctri,
                control.compute=ctrc))

}

### time, nfn and theta
sapply(results, function(r) r$cpu)
sapply(results, function(r) r$misc$nfunc)
sapply(results, function(r) unname(r$mode$theta))

### compare with Table 2 of Cameletti et. al. 2013
## (UTMX and UTMY were not included here)
round(results$'102'$summary.fixed[, c(1,2,3,4,5)], 4)

### with the non-separable
round(results$'121'$summary.fixed[, c(1,2,3,4,5)], 4)

### user parametrization marginals
marginals <- lapply(results, function(r)
    list(sigma.e=inla.tmarginal(
             function(x) exp(-x/2),
             r$internal.marginals.hyperpar[[1]]),
         srange=inla.tmarginal(
             function(x) exp(x),
             r$internal.marginals.hyperpar[[2]]),
         trange=inla.tmarginal(
             function(x) exp(x),
             r$internal.marginals.hyperpar[[3]]),
         sigma.u=inla.tmarginal(
             function(x) exp(x),
             r$internal.marginals.hyperpar[[4]])))

### user interpretable parameters summary
margz <- lapply(marginals, lapply, function(m)
    unlist(inla.zmarginal(m, silent=TRUE)))
lapply(lapply(margz, data.frame), function(x)
    round(t(x), 2))

### Compare ing "102" with Table 3 in Cameletti et. al. 2022
### Consider the temporal range relation and `a` as `a = exp(-2/r_s)`
c(s2e=margz$'102'$sigma.e[1]^2,
  s2w=margz$'102'$sigma.u[1]^2,
  rho=margz$'102'$srange[1],
  a=exp(-sqrt(8*0.5)/margz$'102'$trange[1]))

### fit statistics
round(t(sapply(results, stats.inla, y=log(pdata$PM10),
               fsummarize=function(x) mean(x, na.rm=TRUE))), 4)

