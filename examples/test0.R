
library(INLAspacetime)
library(inlabru)

inla.setOption(
    smtp='pardiso',
    inla.mode='compact',
    pardiso.license='~/.pardiso.lic')

smesh <- inla.mesh.2d(cbind(0,0), max.edge=5, offset=2)

tmesh <- inla.mesh.1d(0:5)

n <- 5
dataf <- data.frame(
    s1=runif(n, -1, 1),
    s2=runif(n, -1, 1),
    time=runif(n, 1, 4),
    y=rnorm(n, 0, 1))
str(dataf)

#### mappter from the model domain to the data
stmapper <- bru_mapper_multi(
    list(space = bru_mapper(smesh),
         time = bru_mapper(tmesh, indexed=TRUE)))

### define the data Model
M <- ~ -1 + Intercept(1) + 
    field(list(space = cbind(s1, s2), time=time),
          mapper=stmapper, model=stmodel)

### likelihood precision prior
lkprec <- list(prec=list(initial=10, fixed=TRUE))

### define the spacetime model
stmodel <- stModel.define(
    smesh, tmesh, '121', 
    control.priors=list(
        prs=c(1, 0.0),
        prt=c(5, 0.0),
        psigma=c(1, 0.5)))

result <- 
    bru(M, 
        like(formula = y ~ ., 
             family="gaussian",
             control.family = list(
                 hyper = lkprec), 
             data=dataf),
        options = list(
            verbose=TRUE)
        )

result$summary.hyperpar
