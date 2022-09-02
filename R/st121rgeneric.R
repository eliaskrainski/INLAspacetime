#' Define the spacetime model with alpha_t=1, alpha_s=2 and alpha_E=1.
#'
#' This function define the object set for using the 'rgeneric'
#' model tool in INLA.
#'
#' @param cmd see the 'rgeneric.define' in INLA.
#' @param theta hyperparameters of the model.
#' @param args arguments (for example, n: size of the model).
#' @param ... objects for internal use.
#' @details
#' See the paper.
#' @return an rgeneric INLA model.
#' @export
st121rgeneric <-
  function(cmd=c("graph", "Q", "mu", "initial",
                 "log.norm.const", "log.prior", "quit"),
           theta = NULL, args = NULL, ...)
  {

    ### check mesh objects
    stopifnot(inherits(obj$mesh.time, "inla.mesh.1d"))
    stopifnot(inherits(obj$mesh.space, "inla.mesh"))

    ### check priors
    stopifnot(length(obj$prior.rs)==2)
    stopifnot(length(obj$prior.rt)==2)
    stopifnot(length(obj$prior.sigma)==2)
    stopifnot(is.numeric(obj$prior.rs))
    stopifnot(is.numeric(obj$prior.rt))
    stopifnot(is.numeric(obj$prior.sigma))

    ### compute lambdas
    prior.ipar <- c(obj$prior.rs[1], obj$prior.rt[1], obj$prior.sigma[1])
    prior.p <- c(obj$prior.rs[2], obj$prior.rt[2], obj$prior.sigma[2])
    lambdas <- c(-log(prior.p)/prior.ipar)

    ## Set fixed.theta if it exists
    is.fixed.theta = is.na(prior.p)
    fixed.theta <- as.numeric(rep(NA, 3))
    fixed.theta[is.fixed.theta] <- log(prior.ipar[is.fixed.theta])

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
      ##     stopifnot(abs(tfe$c0[1,1] -
      ##                 0.5*tfe$c0[2,2])<1e-3)
      N = nrow(tfe$c0)
      result$M1 <- sparseMatrix(
        i=c(1,N), j=c(1,N),
        x=0.5*(tmesh$boundary!='cyclic'),
        dims=c(N, N))
      result$M2 <- uM(tfe$g1) ### do not need the lower

      result$M0c0 <- uM(result$M0 %x% result$c0)
      result$M0g1 <- uM(result$M0 %x% result$g1)
      result$M0g2 <- uM(result$M0 %x% result$g2)
      result$M0g3 <- uM(result$M0 %x% result$g3)

      result$M1c0 <- uM(result$M1 %x% result$c0)
      result$M1g1 <- uM(result$M1 %x% result$g1)
      result$M1g2 <- uM(result$M1 %x% result$g2)

      result$M2c0 <- uM(result$M2 %x% result$c0)
      result$M2g1 <- uM(result$M2 %x% result$g1)

      result$graph <- sparseMatrix(
        i=c(result$M0c0@i, result$M0g1@i,
            result$M0g2@i, result$M0g3@i,
            result$M1c0@i, result$M1g1@i, result$M0g2@i,
            result$M2c0@i, result$M2g1@i)+1,
        j=c(result$M0c0@j, result$M0g1@j,
            result$M0g2@j, result$M0g3@j,
            result$M1c0@j, result$M1g1@j, result$M0g2@j,
            result$M2c0@j, result$M2g1@j)+1, x=1)
      result$graph@x <- rep(1.0, length(result$graph@x))

      ## c1 constant for the marginal variance depends on spatial manifold
      ##        if(smesh$manifold=='S2') {
      ##          result$c1 <- gamma(0.5) * gamma(1) / (gamma(1)*gamma(2) * 4 * pi^0.5)
      ##    } else {
      result$c1 <- gamma(0.5) * gamma(1) / (gamma(1)*gamma(2) * 8 * pi^1.5)
      ##  }

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
      if (is.null(fixed.theta)) {
        ## Assume the input theta is log(rs, rt, sigma)
        theta.interpret = theta
      } else {
        stopifnot(length(fixed.theta)<4)
        theta.interpret = fixed.theta
        theta.interpret[is.na(fixed.theta)] = theta
      }
      ## Note that theta is log(range_s, range_t, sigma)
      alpha_t = 1; alpha_s = 2; alpha_e = 1;
      alpha = alpha_e + alpha_s*(alpha_t-1/2);
      nu.s = alpha -1; nu.t = alpha_t - 0.5;
      ## old&davids:
      ##c1 = gamma(alpha_t - 1/2)*gamma(alpha-1)/(gamma(alpha_t) *gamma(alpha) *4*sqrt(pi))
      ## New
      ##c1 = gamma(alpha_t - 1/2)*gamma(alpha-1)/(gamma(alpha_t) *gamma(alpha) *8*pi^1.5)
      c1 <- lobj$c1

      ## Define the log-gamma's
      theta.gam = rep(NA, 3)
      theta.gam[2] = 0.5*log(8*nu.s) - theta.interpret[1]
      theta.gam[1] = theta.interpret[2] - 0.5*log(8*(alpha_t-1/2)) + alpha_s * theta.gam[2]
      theta.gam[3] = 0.5*(log(c1) - 2*theta.gam[1] - (alpha-1)*theta.gam[2] - theta.interpret[3]*2)
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

      ##        qq <- (kronecker(gt^2*lobj$M2,
      ##                       gs^2*lobj$c0 + lobj$g1) +
      ##           kronecker(lobj$M0,
      ##                   gs^6*lobj$c0 +
      ##                 gs^4*lobj$g1 +
      ##               gs^2*lobj$g2) + lobj$c0g3 +
      ##   kronecker(2*gt*lobj$M1,
      ##           gs^4*lobj$c0 +
      ##         2*gs^2*lobj$g1 +
      ##       lobj$g2)) * ge2

      gt2 <- gt^2
      gs2 <- gs^2
      gs4 <- gs2*gs2
      gs4.3 <- 3*gs4
      gs2.2 <- 2*gs2
      gs2.3 <- 3*gs2
      gt.2 <- 2*gt
      gs6 <- gs2*gs4
      qq <- (gt2*(gs2*lobj$M2c0 + lobj$M2g1) +
               (gs6*lobj$M0c0 + gs4.3*lobj$M0g1 +
                  gs2.3*lobj$M0g2 + lobj$M0g3) +
               gt.2*(gs4*lobj$M1c0 +
                       gs2.2*lobj$M1g1 + lobj$M1g2))*ge2

      ##cat('nnz(Q) =', length(qq@x), '\n')
      return(inla.as.sparse(qq)@x)

    }

    mu <- function(n, theta) return(numeric(0))

    log.norm.const <- function(n, theta) return(numeric(0))

    log.prior <- function(n, theta) {
      ## pc priors with 3 lambdas
      if (is.null(fixed.theta)) {
        ## Assume the input theta is log(rt, rs, sigma)
        #theta = theta
      } else {
        ## In this case the prior will not be properly re-scaled
        ## But that makes no difference to INLA (except mlik etc.)
        stopifnot(length(fixed.theta)<4)
        theta.interpret = fixed.theta
        theta.interpret[is.fixed.theta] = theta
        theta = theta.interpret
      }

      ## log prior value
      ##i = 1 # log(rt) # pc prior for range in d=1
      val = log(lambdas[1]) - lambdas[1] * exp(-0.5*theta[1]) +
        log(0.5) - 0.5*theta[1]
      ##i = 2 # log(rs) # pc prior for range in d=2
      val = val + log(lambdas[2]) - lambdas[2] * exp(-theta[2]) +
        -theta[2]
      ##i = 3 # sigma # pc prior for sigma
      val = val + log(lambdas[3]) - lambdas[3] * exp(theta[3]) +
        theta[3]

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

