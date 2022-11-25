### fit the four models for the piemonte dataset

library(INLAspacetime)

inla.setOption(
    inla.mode='experimental',
    num.threads='8:-1',
    smtp='pardiso', 
    inla.call='remote',
    pardiso.license='~/.pardiso.lic')

ctri <- list(
    int.strategy='ccd',
    parallel.linesearch=TRUE)

ctrc <- list(
    config=TRUE,
    waic=TRUE,
    dic=TRUE,
    cpo=TRUE)

if(!any(ls()=='iverbose'))
    iverbose <- !FALSE

### get the dataset
source('piemonte_data.R')

if(FALSE) {

    system.time(wdat <- reshape(pdata[c('Station.ID', 'Date', 'PM10')],
                                direction='wide', timevar='Station.ID', idvar='Date'))
    dim(wdat)
    wdat[1:3, 1:5]

    png('~/github/slides/figures/piemonte_stlines.png', 1200, 700, res=100)
    par(mfrow=c(1,2), mar=c(0,0,0,0))
    plot(pborders, type='l', lwd=2, col=3, asp=1, las=2, axes=FALSE)
    points(locs[,2:3], pch=8, cex=2)
    plot(pborders, type='l', lwd=2, col=3, asp=1, las=2, axes=FALSE)
    stlines(wdat[,-1], SpatialPoints(locs[,2:3]))
    dev.off()

}

table(pdata$time)
if(!any(ls()=='nt'))
    nt <- 182 ### max number of time points to be used
pdata <- pdata[pdata$time<=nt,]

### define a temporal mesh
tmesh <- inla.mesh.1d(1:nt, degree=1)
tmesh$n

if(FALSE) {
    
    bnd <- inla.nonconvex.hull(
        rbind(cbind(locs[,2], locs[,3]),
              cbind(pborders[,1], pborders[,2])), 
        convex=10, concave=50, resolution=100)
    
    locs0 <- inla.mesh.2d(    
        boundary=bnd, 
        max.edge=10,
        cutoff=5)$loc[,1:2]
    
    smesh <- inla.mesh.2d(
        rbind(cbind(locs[,2], locs[,3]), locs0),
        max.edge=150, offset=150, cutoff=15)

} else {

### mesh as in Cameleti et. al. 2012
    smesh <- inla.mesh.2d(
        cbind(locs[,2], locs[,3]),
        loc.domain=pborders, 
        max.edge=c(50, 300), 
        offset=c(10, 140), 
        cutoff=5, 
        min.angle=c(26, 21))

}

smesh$n

par(mfrow=c(1,1), mar=c(0,0,1,0))
plot(smesh, asp=1)
points(locs[,2:3], pch=8, col='red')
lines(pborders, lwd=2, col='green4')

### prepare the covariates
xnames <- c('A', ###'UTMX', 'UTMY',
            'WS', 'TEMP', 'HMIX', 'PREC', 'EMI')
xmean <- colMeans(pdata[, xnames])
xsd <- sapply(pdata[xnames], sd)

### define the projector matrix
A <- inla.spde.make.A(
    smesh, as.matrix(pdata[c('UTMX', 'UTMY')]),
    group=pdata$time, group.mesh=tmesh)
dim(A)
smesh$n * tmesh$n
stopifnot(sum(A)==nrow(pdata))

dsstack <- inla.stack(
    tag='e',
    data=list(y=log(pdata$PM10)),
    effects=list(
        data.frame(b0=1,
                   scale(pdata[xnames], xmean, xsd)),
        spacetime=1:(tmesh$n*smesh$n)),
    A=list(1, A))

diff(range(pdata$UTMX))
diff(range(pdata$UTMY))

sd(pdata$PM10, na.rm=TRUE)

### define the prior parameters
prs <- c(70, 0.5)
prt <- c(50, 0.5)
psigma <- c(20, 0.5)

### likelihood precision prior
lkprec <- list(
    prec=list(prior='pcprec', param=c(1, 0.1)))

stConstr <- list(
    A=matrix(diag(kronecker(
        Diagonal(tmesh$n, colSums(inla.mesh.fem(tmesh)$c1)),
        Diagonal(smesh$n, colSums(inla.mesh.fem(smesh)$c1)))), nrow=1), e=0)

### define the model formula     
ff <- update(
    y~0+b0+f(spacetime, model=cmodel, extraconstr=stConstr), 
    paste('.~.+', paste(xnames, collapse='+')))

theta.ini <- c(3.5, log(70), log(50), log(1))
theta.ini

models <- c('102', '121', '202', '220')

results <- vector('list', 4L); names(results) <- models

for(model in 1:4) {

    cat('Running model', models[model], '\n')
    
    cmodel <- stModel.define(
        smesh, tmesh, models[model], 
        control.priors=list(
            prs=prs, prt=prt, psigma=psigma))    

### fit the separable model
    results[[model]] <- inla(
        ff,
        data=inla.stack.data(dsstack),
        control.predictor=list(
            A=inla.stack.A(dsstack)), 
        control.family=list(hyper=lkprec),
        control.fixed=list(prec=c(b0=1)),
        control.mode=list(theta=theta.ini, restart=TRUE),
        verbose=iverbose,
        control.inla=ctri,
        control.compute=ctrc)

}

sapply(results, function(r) r$cpu)
sapply(results, function(r) r$misc$nfunc)
sapply(results, function(r) unname(r$mode$theta))
sapply(results, function(r)
    stats.inla(r, y=log(pdata$PM10), fsummarize=function(x) mean(x, na.rm=TRUE)))

detach("package:INLAspacetime", unload=TRUE)
library(INLAspacetime)
