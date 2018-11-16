#' Generate a random p x p correlation matrix
#'
#' @title Generate a random p x p correlation matrix
#' @param p A scalar positive integer denoting the size of the correlation matrix
#' @section Details:
#' This function implements the algorithm by Pourahmadi and Wang [1] for generating a random p x p correlation matrix. 
#' Briefly, the idea is to represent the correlation matrix using Cholesky factorization and p(p-1)/2 hyperspherical coordinates
#' (i.e., angles), sample the angles form a particular distribution and then convert to the standard correlation matrix form. The  
#' angles are sampled from a distribution with probability density function sin^k(theta) (0 < theta < pi, k >= 1) using the efficient
#' sampling algorithm described in [2].
#' 
#' @return A random p x p correlation matrix
#' 
#' @note     To cite this package please reference: 
#'
#' Makalic, E. & Schmidt, D. F.
#' An efficient algorithm for sampling from sin^k(x) for generating random correlation matrices
#' arXiv:1809.05212, 2018 \url{https://arxiv.org/abs/1809.05212}
#' 
#' A MATLAB-compatible implementation of the sampler in this package can be obtained from:
#' 
#' \url{https://au.mathworks.com/matlabcentral/fileexchange/68810-randcorr}
#' 
#' @references 
#' 
#'  [1] Mohsen Pourahmadi and Xiao Wang,
#'  Distribution of random correlation matrices: Hyperspherical parameterization of the Cholesky factor,
#'  Statistics & Probability Letters, Volume 106, Pages 5-12, 2015.
#'  
#'  [2] Enes Makalic and Daniel F. Schmidt
#'  An efficient algorithm for sampling from sin^k(x) for generating random correlation matrices,
#'  arXiv:1809.05212, 2018.
#'  
#'   
#' @examples 
#' # -----------------------------------------------------------------
#' # Example 1: Generate a 5x5 correlation matrix
#' C = randcorr(5)
#' 
#' # Example 2: Generate a 1000x1000 correlation matrix
#' C = randcorr(1000)
#'
#'
#'
#' @seealso \code{\link{randcorr.sample.sink}}
#' @export
randcorr <- function(p)
{
  # Check inputs
  if (length(p)>1 || p<2 || p%%1)
  {
    stop("p must be a scalar integer greater than 1.")
  }
  
  # Step 1 - generate angles theta from PDF (sin(theta))^k, k>=1, 0<theta<pi
  e = matrix(1,p,1)
  theta = matrix(0,p,p)
  for(j in 1 : (p-1)) {
    theta[(j+1):p,j] = randcorr.sample.sink( (p-j)*e[(j+1):p] )
  }
  
  # Step 2 - construct lower triangular Cholesky factor
  L = matrix(1,p,p)
  for (i in 2 : p) {
    L[i,2:i] = cumprod( sin(theta[i,1:i-1]) )
  }
  R = cos(theta)
  R[upper.tri(R)] = 0
  L = L * R
  
  # Form correlation matrix
  C = L %*% t(L)

  return(C)
}



#' Sample from the (unnormalized) distribution sin(x)^k, 0 < x < pi, k >= 1
#'
#' @title Sample from the (unnormalized) distribution sin(x)^k, 0 < x < pi, k >= 1
#' @param k The \code{k} parameter of the distribution. If this is a vector, the function draws a random variate for every entry in \code{k}. 
#' @section Details:
#' This code generates samples from the sin(x)^k distribution using the specified vector \code{k}.
#' @return A vector of samples with length equal to the length of \code{k}
#' 
#' @references 
#' 
#'  Enes Makalic and Daniel F. Schmidt
#'  An efficient algorithm for sampling from sin^k(x) for generating random correlation matrices,
#'  arXiv:1809.05212, 2018.
#' 
#' @examples 
#' # -----------------------------------------------------------------
#' # Example 1: Draw a random variate from sin(x), 0<x<pi
#' x = randcorr.sample.sink(1)
#' 
#' # Example 2: Draw a million random variate from sin^3(x), 0<x<pi
#' x = randcorr.sample.sink( matrix(3, 1e6,1) )
#' mean(x)
#' var(x)
#' 
#'
#' @seealso \code{\link{randcorr}}
#' @export
randcorr.sample.sink <- function(k)
{
  N = length(k)
  logconst = 2*log(pi/2)
   
  # Sampling loop - vectorized
  x = matrix(0,N,1)
  accept = rep(FALSE,N)
  while(!all(accept))
  {
    # index of samples that need to be accepted 
    ix = !accept
    T = sum(ix)
    
    # Beta(k+1,k+1) rng
    x[ix] = pi * stats::rbeta(T,k[ix]+1,k[ix]+1)
    
    # Check acceptance
    accept[ix] = log(stats::runif(T))/k[ix] < logconst + log(sin(x[ix])) - log(x[ix]) - log(pi-x[ix])
  }
  return(x)
}
