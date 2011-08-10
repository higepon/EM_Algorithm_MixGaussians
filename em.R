# PRML 9.2 Mixture of Gaussians
#   There are two Gaussians.
library(mvtnorm)
library(RUnit)

responsibility <- function(xn, k, K, pi, mu, sigma) {
  a <- pi[[k]] * dmvnorm(xn, mu[[k]], sigma[[k]]);
  b <- sum(sapply(1:K, function(j) { pi[[j]] * dmvnorm(xn, mu[[j]], sigma[[j]]) }));
  a / b;
}

Estep <- function(xx, pi, mu, sigma) {
  K <- length(mu);
  apply(xx, 1, function(x) {
    sapply(1:K, function(k) {
      responsibility(x, k, K, pi, mu, sigma);
    });
  });
}

nK <- function(xx, k, gammaKn) {
  ret <- 0;
  N <- nrow(xx);
  for(n in 1:N) {
    ret <- ret + gammaKn[k, n];
  }
  ret;
}

muNew <- function(xx, k, gammaKn, nK) {
  N <- nrow(xx);
  sum = c(0, 0);
  for(n in 1:N) {
    sum <- sum + gammaKn[k, n] * xx[n,];
  }
  sum / nK;
}

sigmaNew <- function(xx, k, gammaKn, muKNew, nK) {
  N <- nrow(xx);
  sum = c(0, 0, 0, 0);
  for(n in 1:N) {
    sum <- sum + gammaKn[k, n] * ((xx[n,] - muKNew) %*% t(xx[n,] - muKNew));
  }
  matrix(sum / nK, ncol=2);

  ## matrix(rowSums(apply(xx, 1, function(x) { gammaKn[k, n] * ((x - muKNew) %*% t(x - muKNew)); })) / nK(xx, k, gammaKn), ncol=2);
}

piNew <- function(xx, k, gammaKn, nK) {
  nK / length(xx) * length(xx[1,]);
}

Mstep <- function(xx, gammaKn) {
  K <- length(xx[1,]);
  nKList <- lapply(1:K, function(k) { nK(xx, k, gammaKn); });
  piNext <- lapply(1:K, function(k) { piNew(xx, k, gammaKn, nKList[[k]]); });
  muNext <- lapply(1:K, function(k) { muNew(xx, k, gammaKn, nKList[[k]]); });
  sigmaNext <- lapply(1:K, function(k) { sigmaNew(xx, k, gammaKn, muNext[[k]], nKList[[k]]); });
  list(piNext, muNext, sigmaNext);
}

input.data <- function() {
  pi <- list(0.7, 0.3);
  mu <- list(c(6, 7), c(1, 1));
  sigma <- list(matrix(c(7,0,0,7), 2, 2), matrix(c(10,3,3,10), 2, 2));
  xx <- matrix(c(1, 2, 3, 4, 5, 6),ncol=2, byrow=TRUE);
  list(pi, mu, sigma, xx);
}

## Unit Tests
test.nK <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaKn <- Estep(xx, pi, mu, sigma);
  checkEqualsNumeric(nK(xx, 1, gammaKn), 1.613336, tolerance = 0.0001);
}

test.Estep <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaNk <- Estep(xx, pi, mu, sigma);
  checkEqualsNumeric(gammaNk, matrix(c(0.08630052, 0.5957016, 0.9313342,
                                       0.91369948, 0.4042984, 0.0686658), nrow=2, byrow=T), tolerance=0.0001);
}

test.Mstep <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaNk <- Estep(xx, pi, mu, sigma);
  checkEquals(Mstep(xx, gammaNk),
              list(list(0.5377787, 0.4622212), # pi
                   list(c(4.047560, 5.047560), c(1.781199, 2.781199)), # mu
                   list(matrix(c(1.425674, 1.425674,
                                 1.425674, 1.425674), byrow=T, ncol=2),
                        matrix(c(1.348276, 1.348276,
                                 1.348276, 1.348276), byrow=T, ncol=2))),
              tolerance=0.0001);
}

test.muKNew <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaKn <- Estep(xx, pi, mu, sigma);
  muKNew <- muNew(xx, 1, gammaKn, nK(xx, 1, gammaKn));
  checkEqualsNumeric(muKNew, c(4.047560, 5.047560), tolerance=0.0001);
}

test.sigmaNew <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaKn <- Estep(xx, pi, mu, sigma);
  muKNew <- muNew(xx, 1, gammaKn, nK(xx, 1, gammaKn));
  checkEqualsNumeric(sigmaNew(xx, 1, gammaKn, muKNew, nK(xx, 1, gammaKn)),
                     matrix(c(1.425674, 1.425674,
                              1.425674, 1.425674), byrow=T, ncol=2), tolerance=0.0001);
## sigmaNew
## > (gammaNk[1, 1] * (c(1, 2) - muKNew)  %*% t((c(1, 2) - muKNew)) +
## + gammaNk[1, 2] * (c(3, 4) - muKNew)  %*% t((c(3, 4) - muKNew)) +
## + gammaNk[1, 3] * (c(5, 6) - muKNew)  %*% t((c(5, 6) - muKNew))) / 1.613336
##          [,1]     [,2]
## [1,] 1.425674 1.425674
## [2,] 1.425674 1.425674
## >
}

test.piKNew <- function() {
  input <- input.data();
  pi <- input[[1]];
  mu <- input[[2]];
  sigma <- input[[3]];
  xx <- input[[4]];
  gammaKn <- Estep(xx, pi, mu, sigma);
  k <- 1;
  checkEqualsNumeric(piNew(xx, k, gammaKn, nK(xx, k, gammaKn)), 0.5377787, tolerance=0.0001);
}

# runTestFile("./em.R")

## E, M step


true_pi <- list(0.7, 0.3);
true_mu <- list(c(8, 7), c(1, 1));
true_sigma <- list(matrix(c(4,0,0,3), 2, 2), matrix(c(2,3,3,7), 2, 2));

ancestralSampling1 <- function(n) {
  x <- runif(1);
  if (x[1] <= true_pi[[1]]) {
    rmvnorm(n=1, true_mu[[1]], true_sigma[[1]]);
  } else {
    rmvnorm(n=1, true_mu[[2]], true_sigma[[2]]);
  }
}

ancestralSampling <- function(n) {
  t(sapply(1:n, ancestralSampling1))
}

xx <- ancestralSampling(1000)
#plot(xx);
start_pi <- list(rnorm(1, 0.5), rnorm(1, 0.5));
start_mu <- list(c(rnorm(1, 10), rnorm(1, 10)), c(rnorm(1, 10), rnorm(1, 10)));
start_sigma <- list(matrix(c(1,0,0,1), 2, 2), matrix(c(1,0,0,1), 2, 2));

gammaKn = Estep(xx, start_pi, start_mu, start_sigma)
v <- Mstep(xx, gammaKn);
v[[1]]
for(n in 1:30) {
  gammaKn = Estep(xx, v[[1]], v[[2]], v[[3]]);
  v <- Mstep(xx, gammaKn);
  cat(sprintf("n=%d mu=v%s\n", n, v[[1]]));
}


