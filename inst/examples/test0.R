library(INLA)
library(INLAspacetime)
library(inlabru)

### We need a plain test-example...
inla.setOption(smtp='taucs', inla.mode='compact')

smesh <- fm_mesh_2d(cbind(0,0), max.edge=5, offset=2)
tmesh <- fm_mesh_1d(0:5)

n <- 5
dataf <- data.frame(
    s1=runif(n, -1, 1),
    s2=runif(n, -1, 1),
    time=runif(n, 1, 4),
    y=rnorm(n, 0, 1))
str(dataf)

### define the data Model
M <- ~ -1 + Intercept(1) +
    field(list(space = cbind(s1, s2), time=time),
          model=stmodel)

### define the spacetime model
stmodel <- stModel.define(
    smesh, tmesh, '121',
    control.priors=list(
        prs=c(1, NA), ## fix spatial range to 1
        prt=c(5, 0.0),## fix temporal range to 5
        psigma=c(1, 0.5)))

### print number of non-zeros in Q_u
cat("Number of non-zeros in Q_u:",
    stmodel$f$cgeneric$data$matrices$xx[2], "\n")

### likelihood precision prior
lkprec <- list(prec=list(initial=10, fixed=TRUE))

### fit
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
