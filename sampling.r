# PRML 9.2 Mixture of Gaussians
#   There are two Gaussians.
library(MASS)

pi1 <- 0.3;
pi2 <- 0.7;

mu1 <- c(5, 5);
mu2 <- c(3, 7);

sigma1 <- matrix(c(1,0,0,1), 2, 2);
sigma2 <- matrix(c(1,0,0,2), 2, 2);

ancestralSampling1 <- function(n) {
  x <- runif(1);
  if (x[1] <= pi1) {
    mvrnorm(1, mu1, sigma1);
  } else {
    mvrnorm(1, mu2, sigma2);
  }
}

ancestralSampling1();

lapply(1:10, ancestralSampling1);

