#' This package contains a function to generate a random p x p correlation matrix.
#' This function implements the algorithm by Pourahmadi and Wang [1] for generating a random p x p correlation matrix. 
#' Briefly, the idea is to represent the correlation matrix using Cholesky factorization and p(p-1)/2 hyperspherical 
#' coordinates (i.e., angles), sample the angles from a particular distribution and then convert to the standard correlation 
#' matrix form. The angles are sampled from a distribution with a probability density function proportional to sin^k(theta) (0 < theta < pi, k >= 1) 
#' using the efficient sampling algorithm described in [2].
#' 
#' 
#' For usage, see the examples in \code{\link{randcorr}} and \code{\link{randcorr.sample.sink}}.
#' 
#' @title The randcorr package
#' @docType package
#' @author Daniel Schmidt \email{daniel.schmidt@@monash.edu} 
#' 
#' Faculty of Information Technology, Monash University, Australia
#'
#' Enes Makalic \email{emakalic@@unimelb.edu.au}
#' 
#' Centre for Epidemiology and Biostatistics, The University of Melbourne, Australia
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
#' @keywords correlation matrix, distribution
#' @seealso \code{\link{randcorr}}, \code{\link{randcorr.sample.sink}}
#' @name randcorr-package
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
# @exportPattern "^[[:alpha:]]+"
NULL