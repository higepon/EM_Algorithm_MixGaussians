# PRML 9.2 Mixture of Gaussians
#   There are two Gaussians.
library(mvtnorm)

pi <- list(0.7, 0.3);
mu <- list(c(6, 7), c(1, 1));

sigma <- list(matrix(c(7,0,0,7), 2, 2), matrix(c(10,3,3,10), 2, 2));

ancestralSampling1 <- function(n) {
  x <- runif(1);
  if (x[1] <= pi[[1]]) {
    rmvnorm(n=1, mu[[1]], sigma[[1]]);
  } else {
    rmvnorm(n=1, mu[[2]], sigma[[2]]);
  }
}

ancestralSampling <- function(n) {
  t(sapply(1:n, ancestralSampling1))
}

responsibility <- function(xn, k, K, pi, mu, sigma) {
  a <- pi[[k]] * dmvnorm(xn, mu[[k]], sigma[[k]]);
  b <- sum(sapply(1:K, function(j) { pi[[j]] * dmvnorm(xn, mu[[j]], sigma[[j]]) }));
  a / b;
}

# Estep works
Estep <- function(xx, pi, mu, sigma) {
  K <- length(mu);
  apply(xx, 1, function(...) {
    sapply(1:K, function(k) {
      responsibility(list(...)[[1]], k, K, pi, mu, sigma);
    });
  });
}

nK <- function(xx, k, K, pi, mu, sigma) {
  sum(apply(xx, 1, function(...) { x = list(...)[[1]]; responsibility(x, k, K, pi, mu, sigma); }));
}

muNew <- function(xx, k, K, pi, mu, sigma) {
#  colSums(apply(xx, 1, function(...) { x <- list(...)[[1]]; responsibility(x, k, K, pi, mu, sigma) * x; })) / nK(xx, k, K, pi, mu, sigma);
  apply(xx, 1, function(...) { x <- list(...)[[1]]; responsibility(x, k, K, pi, mu, sigma) * x; });
}

xx <- matrix(c(1, 2, 3, 4, 5, 6),ncol=2, byrow=TRUE);
Estep(xx, pi, mu, sigma);
##            [,1]      [,2]      [,3]
## [1,] 0.08630052 0.5957016 0.9313342
## [2,] 0.91369948 0.4042984 0.0686658

nK(xx, 1, 2, pi, mu, sigma)
## [1] 1.613336

muNew(xx, 1, 2, pi, mu, sigma)

## muNew(list(c(1, 2), c(3, 4), c(5, 6)), 1, 2, pi, mu, sigma)
## nK(list(c(1, 2), c(3, 4), c(5, 6)), 1, 2, pi, mu, sigma)


## sigmaNew <- function(xx, k, K, mu, muNew, sigma) {

#plot(ancestralSampling(1000));

## Estep(list(c(1, 2), c(3, 4), c(5, 6)), pi, mu, sigma);
## > Estep(list(c(1, 2), c(3, 4), c(5, 6)), pi, mu, sigma);
##            [,1]      [,2]      [,3]
## [1,] 0.08630052 0.5957016 0.9313342
## [2,] 0.91369948 0.4042984 0.0686658


## (dmvnorm(c(3,4), mu[[2]], sigma[[2]]) * pi[[2]]) / (dmvnorm(c(3,4), mu[[1]], sigma[[1]]) * pi[[1]] + dmvnorm(c(3,4), mu[[2]], sigma[[2]]) * pi[[2]])
