% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{simulate_alphas_FOHM}
\alias{simulate_alphas_FOHM}
\title{Generate attribute trajectories under the first order hidden Markov model}
\usage{
simulate_alphas_FOHM(Omega, alpha0s, L)
}
\arguments{
\item{Omega}{A 2^K-by-2^K \code{matrix} of transition probabilities from row pattern to column pattern}

\item{alpha0s}{An N-by-K \code{matrix} of subjects' initial attribute patterns.}

\item{L}{An \code{int} of number of time points}
}
\value{
An N-by-K-by-L \code{array} of attribute patterns of subjects at each time point.
}
\description{
Based on the initial attribute patterns and probability of transitioning between different patterns,
create cube of attribute patterns of all subjects across time.
}
\examples{
N = length(Test_versions)
J = nrow(Q_matrix)
K = ncol(Q_matrix)
L = nrow(Test_order)
Jt = J/L
TP <- TPmat(K)
Omega_true <- rOmega(TP)
class_0 <- sample(1:2^K, N, replace = L)
Alphas_0 <- matrix(0,N,K)
for(i in 1:N){
  Alphas_0[i,] <- inv_bijectionvector(K,(class_0[i]-1))
}
Alphas <- simulate_alphas_FOHM(Omega_true, Alphas_0,L)
}
