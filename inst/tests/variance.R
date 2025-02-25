library(parallel)
library(INLA)
library(INLAspacetime)

sphere <- !FALSE

outer(1:2, c(3,1,0.3), function(a,b) sqrt(8*a)/b)

lpars0 <- list(
    srange = log(c(1, 3, 10)), 
    trange = log(c(5,20)),
    sigma = seq(-3, 3, 1))

outer(1:2, exp(lpars0$srange), function(a,b) sqrt(8*a)/b)

sapply(lpars0, exp)

sapply(lpars0, function(b) range(exp(b)))

str(upars <- expand.grid(lpars0))

nt <- 20
tmesh <- inla.mesh.1d(1:nt)

if(sphere){
    smesh <- inla.mesh.create(globe = 10)
} else {
    smesh <- inla.mesh.2d(
        loc.domain = c(0,0), 
        offset = 2,
        max.edge = 0.2,
        n = 4
    )
}
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
        mean(log(diag(inla.qinv(
            stModel.precision(
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
