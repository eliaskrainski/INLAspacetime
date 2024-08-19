
library(INLAspacetime)

a1 <- c(-1.7, -1.9)
a2 <- 0.95

ar2cov(a1[1], a2, k=5)

length(a1s <- seq(-1.999, -1.5, 0.001))

library(microbenchmark)

microbenchmark(
    R = ar2cov(a1s[1:10], 0.999999, k = 50),
    C = ar2cov(a1s[1:10], 0.999999, k = 50, TRUE)
)

microbenchmark(
    R = ar2cov(a1s[1:10], 0.999999, k = 500),
    C = ar2cov(a1s[1:10], 0.999999, k = 500, TRUE)
)

microbenchmark(
    R = ar2cov(a1s, 0.999999, k = 50),
    C = ar2cov(a1s, 0.999999, k = 50, TRUE)
)

microbenchmark(
    R = ar2cov(a1s, 0.999999, k = 500),
    C = ar2cov(a1s, 0.999999, k = 500, TRUE)
)
