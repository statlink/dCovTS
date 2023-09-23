TstarBoot_boot <- function(x, type, testType, p, b, parallel = FALSE) {

  if ( is.matrix(x) ) {
    if ( dim(x)[2] > 1 ) {
      stop( 'Univariate time series only.' )
    } else  x <- x[, 1]
  }
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if ( (b == 0)  |  missing(b) )  stop( 'b must be grater than 0.' )

  n <- length(x)
  MaxLag <- n - 1
  Atilde0 <- ATilde( Rfast::vecdist(x) )
  dvarx <- mean(Atilde0 * Atilde0)
  Wtstar <- Rfast::matrnorm(b, n) ## matrix of Z variables with b/2 rows and n - k columns

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) ) 
    require(doParallel, quiet = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- parallel::makePSOCKcluster( parallel::detectCores() )
    doParallel::registerDoParallel(cl)

    k <- 1:MaxLag
    d <- foreach(k = k, .combine = cbind, 
         .export = c("kernelFun", "crossDist", "ATilde", "bootCov", "bootCor"), .packages = "Rfast") %dopar% {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) < 1e-16 ) {
        return( numeric(b) )
      } else {
        cross <- crossDist(x, lags = k)
        Atilde <- ATilde( cross$A[[ 1 ]] )
        Btilde <- ATilde( cross$B[[ 1 ]] )
        if ( testType == "covariance" ) {
          return( bootCov(Atilde, Btilde, k, b) )
        } else  return( bootCor(Atilde, Btilde, k, b, dvarx) )
      } ##  end  if ( kern == 0 )
    }  ## end foreach
    parallel::stopCluster(cl)
    
  } else {

    d <- matrix( nrow = b, ncol = MaxLag)
    for ( k in 1:MaxLag ) {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) < 1e-16 ) {
        return( numeric(b) )
      } else {
        cross <- crossDist(x, lags = k)
        Atilde <- ATilde( cross$A[[ 1 ]] )
        Btilde <- ATilde( cross$B[[ 1 ]] )
        com <- Atilde * Btilde  

        if ( testType == "covariance" ) {
           d[, k] <- kern^2 * Rfast::rowsums( Wtstar[, 1:(n - k)] %*% com * Wtstar[, 1:(n - k)] ) / (n - k)  ## this is in fact a Mahalanobis distance
        } else  d[, k] <- kern^2 * Rfast::rowsums( Wtstar[, 1:(n - k)] %*% com * Wtstar[, 1:(n - k)] ) / dvarx / (n - k) 
      }  ##  end  if ( kern == 0 )
    }  ## end  for ( k in 1:MaxLag )
 
  } ##  end if ( parallel )

  Rfast::rowsums(d)
}