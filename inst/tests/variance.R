library(parallel)
library(INLA)
library(INLAspacetime)

lpars0 <- list(
    srange = seq(0, 2, 1),
    trange = seq(1.0, 3.0, 1),
    sigma = seq(-3, 3, 1))

sapply(lpars0, function(b) range(exp(b)))

str(upars <- expand.grid(lpars0))

nt <- 20
tmesh <- inla.mesh.1d(1:nt)

ldomain <- cbind(
    x = c(0, 1, 1, 0, 0) * 5,
    y = c(0, 0, 1, 1, 0) * 5
)
mean(apply(ldomain, 2, function(x) diff(range(x))))

smesh <- inla.mesh.2d(
    loc.domain = ldomain, 
    offset = c(1, 2),
    max.edge = c(0.4, 2))
smesh$n

c(nt, smesh$n, nt*smesh$n)

## plot(smesh)

models <- c("102", "121", "202", "220")

mm1 <- lapply(models, function(m) {
    stModel.matrices(
        smesh = smesh,
        tmesh = tmesh,
        model = m)
})

lapply(mm1, function(x) x$TT)

inla.setOption(num.threads = "1:1")

t0 <- Sys.time()
vv <- t(simplify2array(mclapply(1:nrow(upars), function(i) {
    cat(i, "")
    v <- sapply(models, function(m)
        mean(log(diag(inla.qinv(stModel.precision(
            smesh = smesh,
            tmesh = tmesh,
            model = m,
            theta = c(upars$srange[i],
                      upars$trange[i],
                      upars$sigma[i])
        ))))))
    if((i%%10)==0) cat(':', i, "\n")
    return(v)
}, mc.cores = 6L)))
print(Sys.time()-t0)

summary(vv)

cor(upars$sigma, vv)

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(2, 1, 0), las = 1)
for(k in 1:4) {
    plot(exp(vv[, k]/2), exp(upars$sigma), log = "xy")
    abline(0:1)
}

summary(exp(vv))

mm <- lapply(1:4, function(j) {
    lm(vv[, j]/2 ~ as.matrix(upars))
})

lapply(mm, summary)

sapply(mm, coef)

exp(sapply(mm, coef)[1, ])

##params2gammas
