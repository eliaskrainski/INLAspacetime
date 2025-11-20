library(fmesher)
library(INLAspacetime)
library(INLA)

## common
alpha <- 2
dataf <- data.frame(s = NA, y = NA) ## fake data
cfam.pfix <- list( ## fix nugget to "zero"
    hyper = list(prec = list(
                     initial = 20,
                     fixed = TRUE)))
## prior for the log range
pclrange <- function(lrange, lam, d = 2, logdens = FALSE) {
    dh <- 0.5 * d
    out <- log(lam * dh) -dh * lrange - lam * exp(-dh * lrange)
    if(logdens)
        return(out)
    return(exp(out))
}

### test 1: zeros.rm = TRUE
mesh <- fm_mesh_2d(
    loc = cbind(
        rep(0:7, 8),
        rep(0:7, each = 8))*10,
    offset = 30, 
    max.edge = 25, 
    n = 4,
    cutoff = 5
)
mesh$n

plot(mesh)

spde1 <- inla.spde2.pcmatern(
    mesh = mesh, alpha = alpha, 
    prior.range = c(5, 0.05),
    prior.sigma = c(1, NA)
)
Q1 <- forceSymmetric(
    inla.spde.precision(spde = spde1, theta = log(c(50,1))))

spde2 <- cgeneric_sspde(
    mesh = mesh, alpha = alpha,
    control.priors = list(
        prange = c(5, 0.05),
        psigma = c(1, NA)
    ), useINLAprecomp = FALSE)
Q2 <- prec(spde2, theta = log(c(50)))

all.equal(Q1, Q2)

cmode <- list(theta = c(20, log(50)), fixed = TRUE)
fit1 <- inla(
    y ~ 0 + f(s, model = spde1), data = dataf, 
    control.mode = cmode)
fit2 <- inla(
    y ~ 0 + f(s, model = spde2), data = dataf, 
    control.mode = cmode)

q1 <- prec(fit1)
q2 <- prec(fit2)

all.equal(q1, q2)
all.equal(Q1, q2)

### however...
idd <- which(!(paste(q1@i, q1@j) %in%
               paste(q2@i, q2@j)))
length(idd)

qd <- Sparse(sparseMatrix(
    i = q1@i[idd],
    j = q1@j[idd],
    x = q1@x[idd]
))
summary(qd@x)

image(qd)

### test 2: prior
fit1b <- inla(
    y ~ 0 + f(s, model = spde1), data = dataf, 
    control.family = cfam.pfix
)
fit2b <- inla(
    y ~ 0 + f(s, model = spde2), data = dataf, 
    control.family = cfam.pfix
)

grep("fn-calls=", fit1b$logfile, value = TRUE)[1]
grep("fn-calls=", fit2b$logfile, value = TRUE)[1]

par(mfrow = c(1, 2), mar = c(4,4,1,1),
    mgp = c(3,1,0), las = 1, bty = 'n')
plot(function(x) pclrange(x, -log(0.05)*5),
     -2, 8, n = 10001, lwd = 4,
     xlab = 'log(range)', ylab = 'Density')
lines(fit1b$internal.marginals.hyperpar[[1]],
      lwd = 2, col = 2, lty = 2)
lines(fit2b$internal.marginals.hyperpar[[1]],
      lwd = 3, col = 3, lty = 3)
plot(function(x) pclrange(log(x), -log(0.05)*5)/x,
     0.1, 200, n = 10001, lwd = 4,
     xlab = 'range', ylab = 'Density')
lines(fit1b$marginals.hyperpar[[1]],
      lwd = 2, col = 2, lty = 2)
lines(inla.tmarginal(exp, fit2b$marginals.hyperpar[[1]]),
      lwd = 3, col = 3, lty = 3)

## test 3: sphere (constant in cgeneric_sspde accounts for sphere)

if(FALSE){ ## check the constant for shere

    Fab <- function(gs2,a,b,alpha=1) {
        res <- 0.0
        for(k in a:b) {
            res <- res + (2*k+1)/((gs2 + k*(k+1))^alpha)
        }
    }
    IaInf <- 
    
    k = 0:50
    par(mfrow = c(1, 1), mar = c(4,4,1,1))
    plot(function(x) sapply(x, function(xx)
        log(sum((2*k+1)/(4*pi*(xx^2+k*(k+1))^alpha)))),
        0.01, 100, log = 'x', lwd = 8)
    plot(function(x) -2*log(x^2)-log(4*pi),
         0.01, 0.5, add = TRUE, lwd = 3, col = 2, lty = 2)
    plot(function(x) -log(alpha-1)-2*(alpha-1)*log(x)-log(4*pi),
         2, 100, add = TRUE, lwd = 3, col = 2, lty = 2)
    plot(function(x) -log(4*pi) + sapply(x, function(xx)
        .C("CSphere_gamma_alpha",
           as.double(log(xx)*2),
           as.double(alpha),
           const = double(1))$const),
        0.01, 100, lwd = 2, col = 5, lty = 3, add = TRUE)
    
}

mesh <- fm_rcdt_2d_inla(globe = 10)
mesh$n

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(mesh)

param0 <- c(2.5, 0.3)
theta0 <- log(param0)

(kappa0 <- sqrt(8 * (alpha-2/2))/param0[1])
log(kappa0^2)

spde1 <- inla.spde2.pcmatern(
    mesh = mesh,
    alpha = alpha,
    prior.range = c(param0[1], 0.5),
    prior.sigma = c(param0[2], 0.5)
)

Q1 <- forceSymmetric(inla.as.sparse(
    inla.spde.precision(spde = spde1, theta = theta0)))

fit1 <- inla(
    formula = y ~ 0 + f(s, model = spde1), 
    data = dataf,
    control.family = cfam.pfix
)

spde2 <- cgeneric_sspde(
    mesh = mesh,
    alpha = alpha,
    control.priors = list(
        prange = c(param0[1], 0.5),
        psigma = c(param0[2], 0.5)
    ), useINLAprecomp = FALSE
)

Q2 <- forceSymmetric(prec(spde2, theta = theta0))

all.equal(Q1, Q2)

fit2 <- inla(
    formula = y ~ 0 + f(s, model = spde2),
    data = dataf,
    control.family = cfam.pfix)

grep("fn-calls=", fit1$logfile, value = TRUE)[1]
grep("fn-calls=", fit2$logfile, value = TRUE)[1]

q1 <- prec(fit1)
q2 <- prec(fit2)

all.equal(q1, q2)

## posterior (=prior)
par(mfrow = c(1, 2), mar = c(4,4,1,1),
    mgp = c(3,1,0), las = 1, bty = 'n')
plot(function(x) pclrange(x, -log(0.05)*5),
     -2, 8, n = 10001, lwd = 4,
     xlab = 'log(range)', ylab = 'Density')
lines(fit1b$internal.marginals.hyperpar[[1]],
      lwd = 2, col = 2, lty = 2)
lines(fit2b$internal.marginals.hyperpar[[1]],
      lwd = 3, col = 3, lty = 3)
plot(function(x) pclrange(log(x), -log(0.05)*5)/x,
     0.1, 200, n = 10001, lwd = 4,
     xlab = 'range', ylab = 'Density')
lines(fit1b$marginals.hyperpar[[1]],
      lwd = 2, col = 2, lty = 2)
lines(inla.tmarginal(exp, fit2b$marginals.hyperpar[[1]]),
      lwd = 3, col = 3, lty = 3)
