
library(fmesher)
library(INLA)
library(INLAspacetime)

pclrange <- function(lrange, lam, d = 2, logdens = FALSE) {
    dh <- 0.5 * d
    out <- log(lam * dh) -dh * lrange - lam * exp(-dh * lrange) 
    if(logdens)
        return(out)
    return(exp(out))
}

mesh <- fm_mesh_2d(
    loc = cbind(0, 0),
    offset = c(50, 50),
    max.edge = c(10, 30),
    n = 4,
    cutoff = 5
)
mesh$n

plot(mesh)

spde1 <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(50, 0.95),
    prior.sigma = c(1, NA)
)
spde2 <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(5, 0.05),
    prior.sigma = c(1, NA)
)
spde1b <- cgeneric(
    model = "sspde", mesh = mesh, alpha = 2,
    control.priors = list(
        prange = c(50, 0.95),
        psigma = c(1, NA)
    ), useINLAprecomp = FALSE)
spde2b <- cgeneric(
    model = "sspde", mesh = mesh, alpha = 2,
    control.priors = list(
        prange = c(5, 0.05),
        psigma = c(1, NA)
    ), useINLAprecomp = FALSE)

fit1 <- inla(
    y ~ 0 + f(s, model = spde1),
    data = list(y=NA, s=NA),
    control.family = list(
        hyper = list(prec = list(initial = 20, fixed = TRUE))
    )
)
fit2 <- inla(
    y ~ 0 + f(s, model = spde2),
    data = list(y=NA, s=NA),
    control.family = list(
        hyper = list(prec = list(initial = 20, fixed = TRUE))
    )
)
fit1b <- inla(
    y ~ 0 + f(s, model = spde1b),
    data = list(y=NA, s=NA),
    control.family = list(
        hyper = list(prec = list(initial = 20, fixed = TRUE))
    )
)
fit2b <- inla(
    y ~ 0 + f(s, model = spde2b),
    data = list(y=NA, s=NA),
    control.family = list(
        hyper = list(prec = list(initial = 20, fixed = TRUE))
    )
)

par(mfrow = c(2, 1), mar = c(4,4,1,1), mgp = c(2,1,0), las = 1, bty = 'n')
plot(fit1$marginals.hyperpar[[1]], type = 'l', log = 'x')
lines(fit2$marginals.hyperpar[[1]], col = 2, lty = 2, lwd = 2)
plot(function(x) pclrange(log(x), -log(0.95)*50)/x,
     0.1, 200, n = 10001, add = TRUE, ylab = '',
     col = 5, lty = 3, lwd = 3)
plot(function(x) pclrange(log(x), -log(0.05)*5)/x,
     0.1, 200, n = 10001, add = TRUE, ylab = '',
     col = 6, lty = 3, lwd = 3)

lines(fit1b$marginals.hyperpar[[1]], col = 3, lty = 2, lwd = 2)
lines(fit2b$marginals.hyperpar[[1]], col = 4, lty = 2, lwd = 2)

plot(fit1$internal.marginals.hyperpar[[1]], type = 'l', xlim = c(-10, 15))
lines(fit2$marginals.hyperpar[[1]], col = 2, lty = 2, lwd = 2)
lines(fit1b$marginals.hyperpar[[1]], col = 3, lty = 2, lwd = 2)
lines(fit2b$marginals.hyperpar[[1]], col = 4, lty = 2, lwd = 2)
plot(function(x) pclrange(x, -log(0.95)*50),
     -10, 15, n = 10001, add = TRUE, ylab = '',
     col = 5, lty = 3, lwd = 3)
plot(function(x) pclrange(x, -log(0.05)*5),
     -10, 15, n = 10001, add = TRUE, ylab = '',
     col = 6, lty = 3, lwd = 3)
