
library(INLA)

if(!file.exists('cgeneric-demo.so')) {
## generate the so    
    icall <- inla.getOption('inla.call')
    isub <- regexpr('/INLA/', icall)
    cgpath <- paste0(substr(icall, 1, isub+5), 'cgeneric/')
    print(cgpath)
    system(paste0('gcc -c -fpic -I', cgpath, ' ',
                  cgpath, 'cgeneric-demo.c -o cgeneric-demo.o'))
    system('gcc -shared cgeneric-demo.o -o cgeneric-demo.so ')
}

(nn <- 4^(3:8))
(nsim <- 5)
(nths <- paste0('4:', 1:1))
(cmod <- c('generic0', 'cgeneric'))
times <- as.data.frame(cbind(
    n0=rep(nn, each=nsim*2*length(nths)),
    n=NA, model=NA, nth=NA,
    nfn=NA, cpu=NA))

inla.setOption(
    inla.mode='experimental')

f1 <- y ~ -1 + f(i, model='generic0', Cmatrix=Rmat0, constr=FALSE,
                 hyper=list(theta=list(param=c(1,1))))

kk <- 0
for(n0 in nn) {
    (dd <- round(sqrt(n0)*c(0.8, 1.25)))
### define space domain as a grid
    grid <- SpatialGrid(GridTopology(c(0,0), c(1, 1), dd))
    (n <- nrow(xy <- coordinates(grid)))
### build a spatial neighborhood list
    system.time(jj <- lapply(1:n, function(i)
        which(sqrt((xy[i,1]-xy[,1])^2 + (xy[i,2]-xy[,2])^2)==1)))
### build the spatial adjacency matrix
    nneigh <- sapply(jj, length)
    agraph <- sparseMatrix(
        i=rep(1:n, nneigh),
        j=unlist(jj))
    Rmat0 <- Diagonal(n, colSums(agraph)) - agraph
    ijx <- list(i=rep(1:n, nneigh+1),
                j=unlist(lapply(1:n, function(i)
                    c(i, jj[[i]]))),
                x=unlist(lapply(1:n, function(i)
                    c(nneigh[i], rep(-1, length(jj[[i]]))))))
    upp <- ijx$j>=ijx$i
    Rmat <- sparseMatrix(
        i=ijx$i[upp], j=ijx$j[upp], x=ijx$x[upp],
        dims=c(n,n), repr='T')
    xx <- inla.qsample(n=nsim, Q=Rmat0+Diagonal(n=n,x=1e-3),
                       constr=list(A=matrix(1,1,n), e=0))
    for(s in 1:nsim) {
        cat(s, '')
        dtest <- list(y=xx[,s]+rnorm(n), i=1:n)
        for(nth in nths) {
            r0 <- inla(f1, data=dtest, num.threads=nth)
            kk <- kk+1
            times$s[kk] <- s
            times$n[kk] <- n
            times$model[kk] <- cmod[1]
            times$nth[kk] <- nth
            times$nfn[kk] <- r0$misc$nfunc
            times$cpu[kk] <- r0$cpu[['Total']]
            model1 <- inla.cgeneric.define(
                model = "inla_cgeneric_generic0_model",
                shlib = "cgeneric-demo.so",
                n=n, Cmatrix=Rmat, debug=FALSE)
            r1 <- inla(y ~ -1+f(i, model=model1),
                       data=dtest, num.threads=nth)
            kk <- kk+1
            times$s[kk] <- s
            times$n[kk] <- n
            times$model[kk] <- cmod[2]
            times$nth[kk] <- nth
            times$nfn[kk] <- r1$misc$nfunc
            times$cpu[kk] <- r1$cpu[['Total']]
        }
    }
    cat('\n')
}

times$model <- factor(times$mode, cmod)
times$case <- factor(
    paste(times$model, times$nth), 
    paste(rep(cmod, each=length(nths)), nths))

if (require(ggplot2)) {

ggplot(times) +
    geom_boxplot(aes(x=factor(n), y=nfn, fill=case))

ggplot(times) +
    geom_boxplot(aes(x=factor(n), y=cpu, fill=case)) +
    scale_y_log10()

ggplot(times) +
    geom_boxplot(aes(x=factor(n), y=nfn/cpu, fill=case)) +
    scale_y_log10()

head(times,2)

### differences
nn0 <- sort(unique(times$n))
agg.n <- as.data.frame(t(sapply(nn0, function(n)
    colMeans(times[times$n==n &
                 times$case==times$case[1], ][c('nfn', 'cpu')]))))
agg.n

ii.a <- pmatch(times$n, nn0, duplicates.ok=TRUE)

times$dcpu <- times$cpu-agg.n$cpu[ii.a]
times$rcpu <- times$cpu/agg.n$cpu[ii.a]

ggplot(times) +
    geom_boxplot(aes(x=factor(n), y=dcpu, fill=case)) +
    geom_abline(slope=0, intercept=0, color='red', linetype='dashed')+
    xlab('n (time in seconds)') +
    ylab('CPU: cgeneric - generic0')

ggplot(times) +
    geom_boxplot(aes(x=factor(n), y=rcpu, fill=model)) +
    geom_abline(slope=0, intercept=1, color='red', linetype='dashed')+
    xlab('n (time in seconds)') +
    ylab('CPU: cgeneric / generic0') +
    facet_wrap(~nth)

}
