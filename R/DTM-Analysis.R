#' Calculation of the Empirical DTM-Function for all Points in a Sample
#'
#' @param sample A n x d matrix that contains the sample points
#' @param m The mass parameter of the DTM-function
#' @param r Parameter in the definition of the general DTM-function
#' @return
#' The function returns a vector that contains the DTM-function evaluated at every point in the input sample.
#' @references{
#'   \insertAllCited{}
#' }
#' @importFrom TDA dtm
#' @importFrom Rdpack reprompt
#' @noRd
DTM.function.helper <- function(sample, m, r = 2){
  dtmP_2 <- (dtm(X = sample, Grid = sample, m0 = m, r = r))^r
  return(dtmP_2 )
}

#' Calculation of the Empirical DTM-Function
#'
#' This function calculates the (r'th power of the) empirical DTM-function (see \insertCite{chazal2016rates;textual}{DTMdemo} for a definition) for all points in a given sample.  It is a wrapper for \code{\link[TDA]{dtm}}.
#' @param sample either a n x d numeric coordinate matrix  or a list of such matrices. Here, n corresponds to the size of the sample and d to the dimension of the data.
#' @param m  the mass parameter of the (empirical) DTM-function. It is a double in (0,1].
#' @param r a double in the definition of the general (empirical) DTM-function. By definition r>=1 and the default value is 2.
#' @return
#' Let \code{sample} be a n x d matrix. Then, the function returns a vector that contains the value of the (r'th power of the) empirical DTM-function evaluated at each point encoded in \code{sample}.
#'
#' If \code{sample} is a list of input matrices, then a list containing the return of each input matrix is returned.
#' @examples
#' n <- 1000
#'
#' square <- matrix(0,n,2)
#' square[,1] <- runif(n,0,1)
#' square[,2] <- runif(n,0,1)
#'
#' disc <- matrix(0,n,2)
#' radius <- 0.5
#' r <- radius*sqrt(runif(n,0,1))
#' theta <- runif(n,0,1)*2*pi
#' disc[,1] <- r*cos(theta)
#' disc[,2] <- r*sin(theta)
#'
#' sample <- list(square, disc)
#'
#' DTM.function(sample, m = 1)
#' @references{
#'   \insertAllCited{}
#' }
#' @importFrom Rdpack reprompt
#' @export
DTM.function <- function(sample, m, r = 2){
  if(is.list(sample)){
    res <- lapply(sample,DTM.function.helper, m = m, r = r)
  } else {
    res <- DTM.function.helper(sample, m, r)
  }
  return(res)
}

#' Estimation of a DTM-Density
#'
#' @param sample a n x d numeric coordinate matrix.
#' @param m  the mass parameter of the corresponding (empirical) DTM-function. It is a double in (0,1].
#' @param r a double in the definition of the general (empirical) DTM-function. By definition r>=1.
#' @param kernel a character string giving the smoothing kernel to be used. This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "biweight".
#' @param bw 	the smoothing bandwidth to be used. The kernels are scaled such that this is the standard deviation of the smoothing kernel. bw can also be a character string giving a rule to choose the bandwidth. See \code{\link[stats]{bw.nrd}}.
#' @param full_return a boolean specifying whether the output of the function should also include the DTM-function evaluated at each point of the input sample.
#' @param ... additional arguments passed to the function "\code{\link[stats]{density}}"
#' @return
#' Returns an estimator for the DTM-density proposed in \insertCite{proksch2022small;textual}{DTMdemo} for the given sample (i.e., an object of class "\code{\link[stats]{density}}").
#'
#' If \code{full_return} is \code{TRUE} a list with two components is returned. It contains the the DTM-function evaluated at each sample point ("\code{dtm}") and the corresponding kernel density estimator ("\code{kde}").
#' @references{
#'   \insertAllCited{}
#' }
#' @importFrom TDA dtm
#' @importFrom stats density
#' @noRd
DTM.density.helper <- function(sample, m, r = 2, kernel = "biweight", bw = "SJ", full_return = F, ...){
  dtmP_2 <- (dtm(X = sample, Grid = sample, m0 = m, r = r))^r
  dtm.density <- density(dtmP_2, kernel = kernel, bw = bw, ...)
  if(full_return){
    res <- list("dtm" = dtmP_2, "kde" = dtm.density)
  } else {
    res <- dtm.density
  }
  return(res)
}

#' Estimation of DTM-Densities
#'
#' Calculate the DTM-densities estimators proposed in \insertCite{proksch2022small;textual}{DTMdemo} for all input samples.
#' @param sample either a n x d numeric coordinate matrix  or a list of such matrices. Here, n corresponds to the size of the sample and d to the dimension of the data.
#' @param m  the mass parameter of the corresponding (empirical) DTM-functions. It is a double in (0,1].
#' @param r a double in the definition of the general (empirical) DTM-function. By definition \code{r >= 1} and the default is \code{r = 2}.
#' @param kernel a character string giving the smoothing kernel to be used. This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "biweight".
#' @param bw 	the smoothing bandwidth to be used. The kernels are scaled such that this is the standard deviation of the smoothing kernel. bw can also be a character string giving a rule to choose the bandwidth. See \code{\link[stats]{bw.nrd}}. By default \code{bw = "SJ"}.
#' @param full_return a boolean specifying whether the output of the function should also include the DTM-function evaluated at each point of the input sample (disabled by default).
#' @param ... additional arguments passed to the function \code{\link[stats]{density}}.
#' @details
#'  Given a n x d coordinate matrix, the (r'th power of the) empirical DTM-function is calculated at each sample point encoded in said matrix. To this end, the function \code{\link[TDA]{dtm}} is used. The resulting vector is used to estimate the DTM-density of the input sample via a kernel density estimator (see \insertCite{proksch2022small;textual}{DTMdemo} for details).
#'
#'  If the input is a list of coordinate matrices, this procedure is repeated for every entry of the list.
#' @return
#' Let \code{sample} be a n x d matrix. Then, the function returns an estimator for the DTM-density proposed in \insertCite{proksch2022small;textual}{DTMdemo} based the given sample (i.e., an object of class "\code{\link[stats]{density}}"). If \code{full_return} is \code{TRUE}, a list with two components will be returned. It contains the the DTM-function evaluated at each sample point ("\code{dtm}") and the corresponding kernel density estimator ("\code{kde}").
#'
#' If \code{sample} is a list of input matrices, then a list containing the return of each input matrix is returned.
#' @examples
#' n <- 1000
#'
#' square <- matrix(0,n,2)
#' square[,1] <- runif(n,0,1)
#' square[,2] <- runif(n,0,1)
#'
#' disc <- matrix(0,n,2)
#' radius <- 0.5
#' r <- radius*sqrt(runif(n,0,1))
#' theta <- runif(n,0,1)*2*pi
#' disc[,1] <- r*cos(theta)
#' disc[,2] <- r*sin(theta)
#'
#' sample <- list(square, disc)
#'
#' den.list <- DTM.density(sample, m = 1)
#' plot(den.list[[1]], ylim = c(0, 5), main = "")
#' lines(den.list[[2]], col = "red")
#' @references{
#'   \insertAllCited{}
#' }
#' @importFrom Rdpack reprompt
#' @export

DTM.density <- function(sample, m, r = 2, kernel = "biweight", bw = "SJ", full_return = F, ...){
  if(is.list(sample)){
    res <- lapply(sample, DTM.density.helper, m = m, r = r, kernel = kernel, bw = bw, full_return = full_return, ...)
  } else {
    res <- DTM.density.helper(sample, m = m, r = r, kernel = kernel, bw = bw, full_return = full_return, ...)
  }
  return(res)
}


#' Calculate the L1-distance between Two Kernel Density Estimates.
#'
#' @param den.one,den.two two objects with class "\code{\link[stats]{density}}".
#' @param lower,upper the limits of integration. Can be infinite.
#' @param ... additional arguments passed to the function \code{\link[stats]{integrate}}.
#' @return
#' The function returns a matrix that contains an estimate for the absolute integrated difference between den.one and den.two from lower to upper.
#' @importFrom stats integrate
#' @importFrom stats approxfun
#' @noRd
L1.dist <- function(den.one, den.two, lower, upper, ...){
  f  <- approxfun(den.one$x, den.one$y)
  g  <- approxfun(den.two$x, den.two$y)
  f2 <- function(v) ifelse(is.na(f(v)), 0, f(v))
  f3 <- function(v) ifelse(is.na(g(v)), 0, g(v))
  f4 <- function(v) abs(f2(v) - f3(v))
  res <- integrate(f4, lower, upper, ...)$value
  return(res)
}

#' L1-distance between Kernel Density Estimators
#'
#' Calculate the L1-distance distance matrix for a list of objects with class "\code{\link[stats]{density}}".
#' @param density.list a list of objects with class "\code{\link[stats]{density}}".
#' @param lower,upper the limits of integration. Can be infinite.
#' @param ... additional arguments passed to the function "\code{\link[stats]{integrate}}."
#' @return
#' The function returns an estimate for the absolute integrated difference matrix for the elements in \code{density.list} (integration from \code{lower} to \code{upper}).
#' @examples
#' n = 1000
#'
#' square <- matrix(0,n,2)
#' square[,1] <- runif(n,0,1)
#' square[,2] <- runif(n,0,1)
#'
#' disc <- matrix(0,n,2)
#' radius <- 0.5
#' r <- radius*sqrt(runif(n,0,1))
#' theta <- runif(n,0,1)*2*pi
#' disc[,1] <- r*cos(theta)
#' disc[,2] <- r*sin(theta)
#'
#' sample <- list(square, disc)
#'
#' den.list <- DTM.density(sample, m = 1)
#' L1.density.dist(den.list, 0, 1)
#' @importFrom Rdpack reprompt
#' @export
L1.density.dist <- function(density.list, lower, upper, ...){
  k <- length(density.list)
  dist_mat <- outer(seq(k), seq(k),
                   Vectorize(function(i, j) if(i>j){0}else{L1.dist(density.list[[i]], density.list[[j]], lower =lower, upper = upper, ...)})
  )
  return(dist_mat + t(dist_mat))
}