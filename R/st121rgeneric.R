#' Define the spacetime model with alpha_t=1, alpha_s=2 and alpha_E=1.
#'
#' This function define the object set for using the 'rgeneric'
#' model tool in INLA.
#'
#' @param cmd see the 'rgeneric.define' in INLA.
#' @param theta hyperparameters of the model.
#'  Note that theta is log(range_s, range_t, sigma) 
#' @param args arguments (for example, n: size of the model).
#' @param ... objects for internal use.
#' @details
#' See the paper.
#' @return an rgeneric INLA model.
#' @export
spacetime121 <-
    function(cmd=c("graph", "Q", "mu", "initial",
                   "log.norm.const", "log.prior", "quit"),
             theta = NULL, args = NULL, ...)
{
    
    stopifnot(length(obj$lambdas) == 3)
    stopifnot(inherits(obj$mesh.time, "inla.mesh.1d"))
    stopifnot(inherits(obj$mesh.space, "inla.mesh"))
    
    ## Set fixed.theta if it exists
    fixed.theta = obj$fixed.theta
    
    inla.st121.make.objects <- function(smesh, tmesh) {
        
        uM <- function(m) { ### extract the upper
            m <- inla.as.dgTMatrix(m)
            i.u <- which(m@i<=m@j)
            m@i <- m@i[i.u]
            m@j <- m@j[i.u]
            m@x <- m@x[i.u]
            return(m)
        }
        
        result <- inla.mesh.fem(
            smesh, order=3)[
            c('c0', 'g1', 'g2', 'g3')]
        tfe <- inla.mesh.fem(tmesh, order = 2)
        
        result$M0 <- tfe$c0
        N = nrow(tfe$c0)
        
        result$M1 <- sparseMatrix(
            i=c(1,N), j=c(1,N),
            x=1*(tmesh$boundary!='cyclic'),
            dims=c(N, N))
        result$M2 <- tfe$g1
        
        result$M0c0 <- uM(result$M0 %x% result$c0)
        result$M0g1 <- uM(result$M0 %x% result$g1)
        result$M0g2 <- uM(result$M0 %x% result$g2)
        result$M0g3 <- uM(result$M0 %x% result$g3)
        
        result$M1c0 <- uM(result$M1 %x% result$c0)
        result$M1g1 <- uM(result$M1 %x% result$g1)
        result$M1g2 <- uM(result$M1 %x% result$g2)
        
        result$M2c0 <- uM(result$M2 %x% result$c0)
        result$M2g1 <- uM(result$M2 %x% result$g1)

        mnams <- c('M0c0', 'M0g1', 'M0g2', 'M0g3',
                   'M1c0', 'M1g1', 'M1g2', 'M2c0', 'M2g1')
        result$graph <- Reduce(
            '+', lapply(result[mnams], function(m) {
                m@x <- m@x*0.0 + 1.0
                return(m)
            }))
        result$graph@x <- rep(1.0, length(result$graph@x))
        
        c1t <- gamma(1-1/2)/(gamma(1)*(4*pi)^(1/2))
        ## constant for the marginal variance depends on spatial manifold
        if(smesh$manifold=='S2') {
            c1s <- gamma(2-2/2)/(gamma(2)*(4*pi)^(3/2))
        } else {
            if(smesh$manifold=='R2') {
                c1s <- gamma(2-2/2)/(gamma(2)*(4*pi)^(2/2))
            } else {
                stop(paste(smesh$manifold, 'spatial manifold not supported'))
            }
        } 
        result$lc12 <- log(c1t * c1s)

        return(result)
        
    }

    envir = parent.env(environment())
    if (!exists("init.cache", envir = envir)) {
        ## initialize the cache
        assign('lobj',
               inla.st121.make.objects(
                   obj$mesh.space,
                   obj$mesh.time),
               envir=envir)
        assign("init.cache", TRUE, envir = envir)
        ## print(ls(envir = envir))
    }

    interpret.theta <- function(n, theta) {
        ## Assume the imput theta is log(range_s, range_t, sigma)
        if (is.null(fixed.theta)) {
            itheta <- theta
        } else {
            stopifnot(length(fixed.theta)<4)
            itheta = fixed.theta
            itheta[is.na(fixed.theta)] = theta
        }

        ## Define the log-gamma's: log(g_t, g_s, g_e)
        theta.gam = rep(NA, 3)
        theta.gam[2] = 1.03972077083992 - itheta[1]
        theta.gam[1] = itheta[2] -0.693147180559945 + 2*theta.gam[2]
        theta.gam[3] = 0.5*(lobj$lc12 - theta.gam[1]) - itheta[3] - theta.gam[2]
        return(theta.gam)
    }

    graph <- function(n, theta) {
      return(lobj$graph)
    }

    Q <- function(n, theta) {

        theta.gam <- interpret.theta(n, theta)
        gt <- exp(theta.gam[1]) ## squared \gamma_t
        gs <- exp(theta.gam[2]) ## \gamma_s
        ge2 <- exp(2*theta.gam[3]) ## squared \gamma_e
        
        gt2 <- gt^2
        gt.2 <- 2*gt
        gs2 <- gs^2
        gs4 <- gs2*gs2
        gs6 <- gs2*gs4
        gs4.3 <- 3*gs4
        gs2.2 <- 2*gs2
        gs2.3 <- 3*gs2

        qq <- (gs6*lobj$M0c0 + gs4.3*lobj$M0g1 +
               gs2.3*lobj$M0g2 + lobj$M0g3 +
               gt*(gs4*lobj$M1c0 + gs2.2*lobj$M1g1 + lobj$M1g2) +
               gt2*(gs2*lobj$M2c0 + lobj$M2g1))*ge2
      
      return(inla.as.sparse(qq)@x)

    }

    mu <- function(n, theta) return(numeric(0))

    log.norm.const <- function(n, theta) return(numeric(0))

    log.prior <- function(n, theta) {
        ## pc priors lambdas for the corresponding thetas 
        lambdas = obj$lambdas
        if (is.null(fixed.theta)) {
            ## log prior value
            ##i = 2 # log(rs) # pc prior for range in d=2
            val = log(lambdas[1]) - theta[1] - lambdas[1]*exp(-theta[1]) 
            ##i = 1 # log(rt) # pc prior for range in d=1
            val = val + log(lambdas[2]) - 0.5*theta[2] -
                lambdas[2]*exp(-0.5*theta[2]) +log(0.5) 
            ##i = 3 # sigma # pc prior for sigma
            val = val + log(lambdas[3]) + theta[3] - lambdas[3]*exp(theta[3]) 
        } else {
            ## In this case the prior will not be properly re-scaled
            ## But that makes no difference to INLA (except mlik etc.)
            stopifnot(length(fixed.theta)<4)
            theta.interpret = fixed.theta
            theta.interpret[is.na(fixed.theta)] = theta
            theta = theta.interpret
        }        
        
        return(val)
    }

    initial <- function(n, theta) {
      if (is.null(fixed.theta)) {
        return(c(0,0,0))
      } else {
        return(rep(0, sum(is.na(fixed.theta))))
      }
    }

    quit <- function(n, theta) return(invisible())

    cmd <- match.arg(cmd)
    val <- do.call(
      cmd, args=list(n=as.integer(args$n),
                     theta=theta))

    return(val)

  }

