#' @title Closure
#'
#' @description
#' \code{clo} divides each row of \code{X} by its row sum
#'
#' @details If \code{check} is \code{TRUE} then this function will stop if
#' there are any negative or \code{NA} values in \code{X}
#' @param X A matrix or dataframe of positive numeric values
#' @param check A logical scalar
#' @return A version of \code{X} where each row has been scaled so they sum to 1.
#' @examples
#' X <- matrix(1:12, nrow=3)
#' x <- clo(X)
#' rowSums(x)
#' @export
propr.clo <- function(X, check=FALSE){
  if(check){
    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }
  return(sweep(X, 1, rowSums(X), "/"))
}

#######################################################################################
#' @title Centred logratio transformation
#'
#' @description
#' \code{clr} takes the log of each row of X and centres it (i.e., subtracts the mean).
#'
#' @details If \code{check} is \code{TRUE} then this function will stop if
#' there are any negative or \code{NA} values in \code{X}
#' @param X A matrix or dataframe of positive numeric values
#' @param check A logical scalar
#' @return The logarithm of \code{X} where each row has been shifted to have mean 0.
#' \deqn{\mathrm{clr}(x) = \log x_i - \sum_{i=1}^D \log x_i}
#' @examples
#' X <- matrix(1:12, nrow=3)
#' x <- clr(X)
#' rowSums(x)             # Pretty close to zero
#' apply(exp(x), 1, prod) # The row products of exp(x) will be 1
#' @export
propr.clr <- function(X, check=FALSE){
  if(check){
    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }
  logX <- log(X)
  return(sweep(logX, 1, rowMeans(logX), "-"))
}

#######################################################################################
#' @title Variance of logratios
#'
#' @description
#' \code{vlr} returns a matrix where element (i,j) is
#' the variance (over rows) of the log of the ratios of column i and j.
#'
#' @details If \code{check} is \code{TRUE} then this function will stop if
#' there are any negative or \code{NA} values in \code{X}.
#' @param X A matrix or dataframe of positive numeric values
#' @param check A logical scalar
#' @return The symmetric matrix
#' \eqn{\mathrm{Var}{\log(X_i/X_j)}}{Var(log(X_i/X_j))} where \eqn{X_i} and \eqn{X_j}
#' denote \emph{columns} \eqn{i} and \eqn{j} of \eqn{X}.
#' @examples
#' N <- 10 # Number of observations
#' # Make a data frame with columns a and b roughly proportional
#' # and columns c and d roughly proportional
#' X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
#'                 c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
#' round(vlr(X),2)
#' @export
propr.vlr <- function(X, check=FALSE){
  if(check){
    if(any(X < 0))    stop("negative values found")
    if(any(is.na(X))) stop("NA values found")
  }
  logX <- log(X)
  Cov    <- stats::var(logX)  ## Note the avoidance of compositions::var
  D      <- ncol(logX)
  VarCol <- matrix(rep(diag(Cov), D), ncol = D)
  return(-2 * Cov + VarCol + t(VarCol))
}


#######################################################################################
#' @title Symmetric phi statistic
#'
#' @description
#' \code{phisym} returns a matrix where element (i,j) is
#' the symmetric phi statistic between columns i and j of \code{X}.
#' @details \code{X} should be the result of a centred logratio transformation
#' @param X A matrix or dataframe
#' @return TBA.
#' @examples
#' N <- 10 # Number of observations
#' # Make a data frame with columns a and b roughly proportional
#' # and columns c and d roughly proportional
#' X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
#'                 c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
#' round(phisym(clr(X)),2)
#' @export
propr.phisym <- function (X)
{
  Cov    <- stats::var(X)
  tmp    <- 2 * Cov / outer(diag(Cov), diag(Cov), "+")
  return((1-tmp)/(1+tmp))
}

#######################################################################################
#' @title Expected value of phi from Dirichlet log-ratio distributions
#'
#' @description
#' default returns dataframe of the lower-triangle of symmetrical phi metric,
#' alternatively returns matrix of the summetrical phi metric
#' in either case, the value of phi is the expected value of a number of Dirichlet
#' Monte-Carlo replicates of the data. This reduces the problem of
#' 0-count and low-count features being highly variable because their
#' values range wildly and so the expected value is always large
#' @details requires aldex.clr function from ALDEx2
#' param aldex.clr is an S3 object from the aldex.clr function
#' we ignore all the other measures that are used for trouble-shooting phi
#' the sma.df function in particular is very time and memory intensive
#' @examples
#' # use a count table where the samples are by column, features by row
#' x <- aldex.clr(count.table, return="df")
#' # if return = df, returns a dataframe of the expected value of the lower
#' triangle of the propr.phisym function.
#' # if return = mat, returns the symmetric matrix
#' The number of Dirichlet Monte-Carlo replicates is
#' obtained from the aldex.clr object

propr.aldex.phi <- function(aldex.clr, return="df"){

	# calculate expected value of phi
	# a single high phi value will push the component out of consideration
	# a median is right out for memory considerations

	# get first value
	sym.phi <- propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
	    function(y){y[,1]})))

	# sum the rest of the values as we proceed through the DIR MC instances
	for(i in 2:numMCInstances(aldex.clr)){
		#print(i)
		sym.phi <- sym.phi + propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
		    function(y){y[,i]})))
	}
	##### Done ALDEx2 stuff

	# make indices of the correct size
	lt <- which(col(sym.phi)<row(sym.phi), arr.ind=FALSE)
	lt.ind <- which(col(sym.phi)<row(sym.phi), arr.ind=TRUE)

	# dataframe to hold the info,
	# data is a set of vectors where only the lower triangle is kept, because the matrix
	#    is symmetrical
	# this is needed so subsequent subset function works properly
	sma.df <- data.frame(row=factor(rownames(sym.phi)[lt.ind[,"row"]]),
		col=factor(colnames(sym.phi)[lt.ind[,"col"]]))

	#save the lower triangle as an expected value
	sma.df$phi <- sym.phi[lt] /  numMCInstances(aldex.clr)

	if(return=="df") return(sma.df)
	if(return=="mat") return(sym.phi /  numMCInstances(aldex.clr))
}

#######################################################################################
#######################################################################################
#' @title Standardised Major Axis fits of pairs of columns
#'
#' @description
#' \code{sma} returns a list whose elements are matrices whose elements (i,j)
#' relate to the Standardised Major Axis fits of columns i and j of \code{X}
#' @details \strong{Note:} \code{X} should be the result of a centred logratio transformation
#' @param X A matrix or dataframe
#' @return A list of three elements \code{b}, \code{p} and \code{r2}.
#' @examples
#' N <- 100 # Number of observations
#' # Make a data frame with columns a and b roughly proportional
#' # and columns c and d unrelated
#' a <- seq(from=5, to=15, len=N)
#' b <- a * rnorm(N, 1, 0.1)
#' c <- rnorm(N,10)
#' d <- rnorm(N,10)
#' X <- data.frame(a, b, c, d)
#' pairs(X)
#' pairs(clr(X)) # Note the spurious correlation between variables c and d
#' sma(clr(X))
#' @export
propr.sma <- function(X){
  X.cor <- stats::cor(X, use="pairwise.complete.obs")
  X.var <- stats::cov(X, use="pairwise.complete.obs")
  X.sd  <- sqrt(diag(X.var))

  # Following the approach of Warton et al. Biol. Rev. (2006), 81, pp. 259-291
  # r.rf2 = cor(X+Y, X-Y)^2
  #       = (var(X) - var(Y))^2  /  ((var(X) + var(Y))^2 - 4cov(X,Y)^2)
  r.rf2   <-
    (outer(diag(X.var), diag(X.var), "-")^2 ) /
    (outer(diag(X.var), diag(X.var), "+")^2 - 4 * X.var^2 )

  # At this point the diagonal of r.rf2 will be 0/0 = NaN. The correlation should be 0
  diag(r.rf2) <- 0
  res.dof     <- nrow(X) - 2
  F           <- r.rf2/(1 - r.rf2) * res.dof

  list(b=sign(X.cor) * outer(X.sd, X.sd, "/"), # slope = sign(s_xy) s_y/s_x
       p=1 - pf(F, 1, res.dof),                # p-value of the test that b = 1
       r2=X.cor^2)                             # the squared correlation coefficient
}



#######################################################################################
#' @title Pairwise proportionality, slope and other statistics of columns
#'
#' @description
#' \code{phiDF} returns a dataframe of various statistics for each pair of columns in
#' \code{X}.
#' @details \strong{Note:} \code{X} should be the result of a centred logratio transformation
#' @param X A matrix or dataframe
#' @return A list of three elements \code{b}, \code{p} and \code{r2}.
#' @examples
#' N <- 100 # Number of observations
#' # Make a data frame with columns a and b roughly proportional
#' # and columns c and d unrelated
#' a <- seq(from=5, to=15, len=N)
#' b <- a * rnorm(N, 1, 0.1)
#' c <- rnorm(N,10)
#' d <- rnorm(N,10)
#' X <- data.frame(a, b, c, d)
#' pairs(X)
#' pairs(clr(X)) # Note the spurious correlation between variables c and d
#' phiDF(X)
#' # Note that phi and phisym are related to the slope and r^2 values:
#' with(phiDF(X),
#'   all.equal(
#'     phisym,
#'     (1 + b^2 - 2*b*sqrt(r2))/(1 + b^2 + 2*b*sqrt(r2))
#'   )
#' )
#' with(phiDF(X),
#'   all.equal(
#'     phi,
#'     (1 + b^2 - 2*b*sqrt(r2))
#'   )
#' )
#' @export
propr.phiDF <- function(X){
  X.clr          <- clr(X)
  X.sma          <- sma(X.clr)
  X.vlr          <- vlr(X)
  X.clr.var      <- apply(X.clr, 2, var)  # The variance of each column
  X.phi          <- sweep(X.vlr, 2, X.clr.var, FUN="/")
  X.phisym       <- phisym(X.clr)
  lt             <- which(col(X.sma$b)<row(X.sma$b),arr.ind=FALSE)
  lt.ind         <- which(col(X.sma$b)<row(X.sma$b),arr.ind=TRUE)
  result         <- data.frame(
    row=factor(rownames(X.sma$b)[lt.ind[,"row"]]),
    col=factor(colnames(X.sma$b)[lt.ind[,"col"]])
  )
  result$b       <- X.sma$b[lt]
  result$p       <- X.sma$p[lt]
  result$r2      <- X.sma$r2[lt]
  result$vlr     <- X.vlr[lt]
  result$phi     <- X.phi[lt]
  result$phisym  <- X.phisym[lt]
  return(result)
}

