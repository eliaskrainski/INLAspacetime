
library(INLAspacetime)
library(INLA)

alpha <- 2
param0 <- c(0.7, 0.3)
theta0 <- log(param0)

(kappa0 <- sqrt(8 * (alpha-2/2))/param0[1])
log(kappa0^2)

globe <- !FALSE

if(globe) {
    if(FALSE){ ## check

        k = 0:50
        plot(function(x) sapply(x, function(xx) log(sum((2*k+1)/(4*pi*(xx^2+k*(k+1))^alpha)))), 0.01, 100, log = 'x')
        plot(function(x) -2*log(x^2)-log(4*pi), 0.01, 0.5, add = TRUE, col=2, lty = 2)
        plot(function(x) -log(alpha-1)-2*(alpha-1)*log(x)-log(4*pi), 2, 100, add = TRUE, col=2, lty = 2)
        plot(function(x) -log(4*pi) + sapply(x, function(xx) .C("CSphere_gamma_alpha", log(xx)*2, as.integer(alpha), c=double(1))$c), 0.01, 100, col = 5, lty = 2, add = TRUE)
        
    }
    mesh <- fm_rcdt_2d_inla(globe = 10)
} else {
    mesh <- fm_mesh_2d(
        loc = c(0, 0),
        max.edge = param0[1]/5,
        offset = param0[1]*2,
        n = 7,
        cutoff = param0[1]/10
    )
}
mesh$n

plot(mesh)

spde <- inla.spde2.pcmatern(
    mesh = mesh,
    alpha = alpha,
    prior.range = c(param0[1], 0.5),
    prior.sigma = c(param0[2], 0.5)
)

Q0 <- inla.spde.precision(spde = spde, theta = theta0)

cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
ccomp <- list(config = TRUE)
cmode <- list(theta = theta0, fixed = TRUE)

dataf <- data.frame(
    spatial = 1:spde$n.spde,
    y = NA)

fit0 <- inla(
    formula = y ~ 0 + f(spatial, model = spde), 
    data = dataf,
    control.family = cfam,
    control.compute = ccomp,
    control.mode = cmode
)

csspde <- cgeneric_sspde(
    mesh = mesh,
    alpha = 2,
    control.priors = list(
        prange = c(param0[1], 0.5),
        psigma = c(param0[2], 0.5)
    ), useINLAprecomp=FALSE
)

Q1 <- inla.as.sparse(graphpcor::precision(csspde, theta = theta0))

fit1 <- inla(
    formula = y ~ 0 + f(spatial, model = csspde),
    data = dataf,
    control.family = cfam,
    control.compute = ccomp,
    control.mode = cmode
)

Q0i <- inla.as.sparse(forceSymmetric(fit0$misc$configs$config[[1]]$Q))
Q1i <- inla.as.sparse(forceSymmetric(fit1$misc$configs$config[[1]]$Q))

c(all.equal(Q0, Q0i),
  all.equal(Q1, Q1i),
  all.equal(Q0, Q1))

param0[2]
summary(apply(inla.qsample(100, Q0), 2, sd))
summary(apply(inla.qsample(100, Q1), 2, sd))
