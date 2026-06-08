
library(INLA)

args(inla.matern.cov)

c(inla.matern.cov(nu=1, kappa = 1, x = 1, d = 1, corr = TRUE),
  inla.matern.cov(nu=1, kappa = 1, x = 1, d = 2, corr = TRUE),
  inla.matern.cov(nu=1, kappa = 1, x = 1, d = 1, corr = FALSE),
  inla.matern.cov(nu=1, kappa = 1, x = 1, d = 2, corr = FALSE))

library(INLAspacetime)

(nnu <- length(nu <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 10)))
range <- c(1, 5)

nsub <- 1001
a = 1/nsub; b = 10
x <- c(.Machine$double.eps,
       seq(0, b, length=nsub)[-1])

par(mfrow = c(2, 4), mar = c(4,4,0.5,0.5),
    mgp = c(1.5,0.5,0), bty = 'n', las = 1)
for(i in 1:nnu) {
    plot(x, 
         inla.matern.cov(
             nu[i], kappa = sqrt(8 * nu[i])/range[1],
             x, d = 2, corr = TRUE), type = "l", lwd = 2,
         xlab = "Distance", ylab = "Correlation")
    lines(x, .C("cWMatern", as.integer(nsub), 1.0, 
                sqrt(8*nu[i])/range[1], nu[i],
                as.double(x), r=double(nsub), 
                PACKAGE = "INLAspacetime")$r,
          lty = 2, col = 3, lwd = 2)
    lines(x, 
          inla.matern.cov(
              nu[i], kappa = sqrt(8 * nu[i])/range[2],
              x, d = 2, corr = TRUE), col = 2, lwd = 2)
    lines(x, .C("cWMatern", as.integer(nsub), 1.0,
                sqrt(8*nu[i])/range[2], nu[i],
                as.double(x), r=double(nsub), 
                PACKAGE = "INLAspacetime")$r,
          lty = 2, col = 3, lwd = 2)
    legend("topright", paste("range = ", range), bty = "n",
           title = as.expression(bquote(nu == .(nu[i]))))
    abline(h = c(0.9, 0.1, 0.139), v = range,
           lty = 2, col = gray(0.5))
}
