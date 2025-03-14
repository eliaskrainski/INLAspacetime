## quote form the spacetime SPDE paper, page 31 at
## https://www.idescat.cat/sort/sort481/48.1.1.Lindgren-etal.pdf
## "the computational costs of the sep-arable and non-separable
##  models are similar, as the sparsity structure of posterior
##  precisions, given irregularly spaced observations in
##  generalised latent Gaussian models, is only marginally affected
##  by the non-separability, and can even be more sparse in the
##  non-separable cases; the separable precision neighbourhood
##  structures are space-time prisms, whereas the non-separable
##  neighbourhood structures are double-cones."

library(INLA)
library(INLAspacetime)

sresol <- 50
nt <- 10

nlocs <- 10000 ## number of data locations in space
ndata <- nt * nlocs ## number of observations

sphere <- !FALSE

if(sphere) {
  srange <- 0.3
  srange.u <- 0.5
} else {
  srange <- 3
  srange.u <- 5
}
ssigma <- 2
trange.u <- nt/3
sigma.u <- 1
sigma.e <- 0.1

theta.s <- log(c(srange, ssigma))
thetau.fix <- log(c(srange.u, trange.u, sigma.u))
theta.fix <- c(-2 * log(sigma.e), theta.s, thetau.fix)    

##for(sresol in c(50, 100, 200, 500, c(1,2,4,8,16)*1000)) {
if(sphere) {
  bb <- rbind(c(0, 2), c(0, 2)) * srange.u
  locs <- inla.mesh.map(
    loc = cbind(runif(nlocs, -180, 180), 
                runif(nlocs, -90, 90)),
    projection = "longlat", 
    inverse = TRUE
  )
  smesh <- inla.mesh.create(globe = round(0.32 * sqrt(sresol))) 
} else {
  bb <- rbind(c(0, 2), c(0, 2)) * srange.u
  pol <- matrix(bb[c(1,3,3,1,1, 2,2,3,3,2)], ncol = 2)
  locs <- cbind(
    x = runif(nlocs, bb[1, 1], bb[1, 2]),
    y = runif(nlocs, bb[2, 1], bb[2, 2]))
  sr <- 7 * srange.u / (sresol^0.625)
  smesh <- inla.mesh.2d(
    loc.domain = locs, 
    cutoff = 1 * sr,
    offset = c(3, 7) * sr,
    max.edge = c(2, 5) * sr
  )
}
(ns <- smesh$n)
cat('sresol = ', sresol, ', # spatial nodes', ns, '\n')
##}

if(FALSE)
  plot(smesh, asp = 1)

B <- cbind(
  b0 = 1,
  rep(sin(locs[, 1] - diff(bb[1,])), nt),
  rep(sin(locs[, 2] - diff(bb[2,])), nt),
  rep(sin(2*pi * 1:nt/nt), each = nlocs),
  rep(cos(2*pi * 1:nt/nt), each = nlocs)
); colnames(B) <- paste0("b", 0:(ncol(B)-1))

beta <- c(30, -1, 1, -0.5, -0.3)

## spatial
smodel <- inla.spde2.pcmatern(
  mesh = smesh,
  prior.range = c(srange.u/2, 0.01),
  prior.sigma = c(sigma.u*3, 0.01)
)

Q.s <- inla.spde2.precision(
  spde = smodel,
  theta = theta.s)

s0 <- inla.qsample(
  n = 1,
  Q = Q.s)[, 1]

A.s <- inla.spde.make.A(
  mesh = smesh,
  loc = kronecker(matrix(1, nt, 1), locs)
)

ssim <- drop(as.matrix(A.s %*% s0))

tmesh <- inla.mesh.1d(
  loc = 1:nt)

A.st <- inla.spde.make.A(
  mesh = smesh,
  loc = kronecker(matrix(1, nt), locs),
  group = rep(1:nt, each = nlocs),
  group.mesh = tmesh
)

p0 <- ns * (nt-5)
pn <- ns * 5

imodels <- c('102', '121', '202', '220')

QQ <- vector('list', length(imodels))

for(i in 1:length(imodels)) {
    
    model.id <- imodels[i]
    cat('Computing for model', model.id, '\n')

    stmodel <- stModel.define(
        smesh = smesh,
        tmesh = tmesh,
        model = model.id,
        control.priors = list(
            prs = c(srange.u/3, 0.01),
            prt = c(trange.u/3, 0.01),
            psigma = c(3 * sigma.u, 0.01)
        )
    )
    
    Qu <- stModel.precision(
        smesh = smesh,
        tmesh = tmesh,
        model = model.id,
        theta = thetau.fix
    )
    
    u0 <- inla.qsample(
        n = 1,
        Q = Qu)[, 1]
    
    usim <- drop(as.matrix(A.st %*% u0))

    ysim <- drop(B %*% beta) + ssim + usim + 
        rnorm(ndata, mean = 0, sd = sigma.e)
    
    sdata <- inla.stack(
        data = list(y = ysim),
        effects = list(
            as.data.frame(B),
            spatial = 1:smodel$n.spde,
            field = 1:stmodel$f$n),
        A = list(1, A.s, A.st)
    )
    
    imodal <- inla(
        y ~ 0 + b0 + b1 + b2 + b3 + b4 +
            f(field, model = stmodel) + f(spatial, model = smodel),
        data = inla.stack.data(sdata),
        control.fixed = list(prec.intercept = 1, prec = rep(1, 5)),
        control.predictor = list(
            A = inla.stack.A(sdata)),
        control.mode = list(theta = theta.fix, fixed = TRUE),
        control.compute = list(config = TRUE),
        ##    inla.call = "remote",
        verbose = FALSE
    )


    QQ[[i]] <- imodal$misc$configs$config[[1]]$Q
}


library(gridExtra)

ii <- list((1):(pn),
           ceiling(nt/2)*ns + (1):(pn),
           (p0+1):(p0+pn))

nnz <- sapply(QQ, function(x) length(x@x)*2-nrow(x))
nnz
nnzr <- nnz/sapply(QQ,nrow)

grid.arrange(
    image(QQ[[1]][ii[[1]], ii[[1]]],
          main = paste0(imodels[1], ': nnz/r = ', format(nnzr[1], digits = 2))),
    image(QQ[[2]][ii[[1]], ii[[1]]],
          main = paste0(imodels[2], ': nnz/r = ', format(nnzr[2], digits = 2))),
    image(QQ[[3]][ii[[1]], ii[[1]]],
          main = paste0(imodels[3], ': nnz/r = ', format(nnzr[3], digits = 2))),
    image(QQ[[4]][ii[[1]], ii[[1]]],
          main = paste0(imodels[4], ': nnz/r = ', format(nnzr[4], digits = 2))),
    ncol = 2)


