# PRML 9.2 Mixture of Gaussians
#   There are two Gaussians.
library(mvtnorm)
library(RUnit)

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
  rowSums(apply(xx, 1, function(x) { responsibility(x, k, K, pi, mu, sigma) * x; })) / nK(xx, k, K, pi, mu, sigma);
}


sigmaNew <- function(xx, k, K, pi, mu, sigma, muKNew) {
  matrix(rowSums(apply(xx, 1, function(x) { responsibility(x, k, K, pi, mu, sigma) * ((x - muKNew) %*% t(x - muKNew)); })) / nK(xx, k, K, pi, mu, sigma), ncol=2);
}

  xx <- matrix(c(1, 2, 3, 4, 5, 6),ncol=2, byrow=TRUE);
gammaNk <- Estep(xx, pi, mu, sigma);
gammaNk
##            [,1]      [,2]      [,3]
## [1,] 0.08630052 0.5957016 0.9313342
## [2,] 0.91369948 0.4042984 0.0686658


input.data <- function() {
  pi <- list(0.7, 0.3);
  mu <- list(c(6, 7), c(1, 1));
  sigma <- list(matrix(c(7,0,0,7), 2, 2), matrix(c(10,3,3,10), 2, 2));
  xx <- matrix(c(1, 2, 3, 4, 5, 6),ncol=2, byrow=TRUE);
  list(pi, mu, sigma, xx);
}


muKNew <- muNew(xx, 1, 2, pi, mu, sigma)
## [1] 4.047560 5.047560

sigmaNew(xx, 1, 2, pi, mu, sigma, muKNew)
##          [,1]     [,2]
## [1,] 1.425674 1.425674
## [2,] 1.425674 1.425674

## sigmaNew
## > (gammaNk[1, 1] * (c(1, 2) - muKNew)  %*% t((c(1, 2) - muKNew)) +
## + gammaNk[1, 2] * (c(3, 4) - muKNew)  %*% t((c(3, 4) - muKNew)) + 
## + gammaNk[1, 3] * (c(5, 6) - muKNew)  %*% t((c(5, 6) - muKNew))) / 1.613336
##          [,1]     [,2]
## [1,] 1.425674 1.425674
## [2,] 1.425674 1.425674
## > 


#plot(ancestralSampling(1000));

## Estep(list(c(1, 2), c(3, 4), c(5, 6)), pi, mu, sigma);
## > Estep(list(c(1, 2), c(3, 4), c(5, 6)), pi, mu, sigma);
##            [,1]      [,2]      [,3]
## [1,] 0.08630052 0.5957016 0.9313342
## [2,] 0.91369948 0.4042984 0.0686658


## (dmvnorm(c(3,4), mu[[2]], sigma[[2]]) * pi[[2]]) / (dmvnorm(c(3,4), mu[[1]], sigma[[1]]) * pi[[1]] + dmvnorm(c(3,4), mu[[2]], sigma[[2]]) * pi[[2]])

## Unit Tests
test.nK <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  checkEqualsNumeric(nK(xx, 1, 2, pi, mu, sigma), 1.613336, tolerance = 0.0001);
}

test.Estep <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaNk <- Estep(xx, pi, mu, sigma);
  checkEqualsNumeric(gammaNk, matrix(c(0.08630052, 0.5957016, 0.9313342, 0.91369948, 0.4042984, 0.0686658), nrow=2, byrow=T), tolerance=0.0001);
}