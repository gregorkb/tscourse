#'Computes the sample autocovariance and autocorrelation functions.
#'
#' @param x a vector containing time series data.
#' @param max.lag the largest lag at which to compute the autocovariance and autocorrelation functions.
#' @return a list containing the autocovariance at lags zero through \code{max.lag}.
#' This function computes the sample autocovariance and autocorrelation functions based on the input data.
sample.acf <- function(x,max.lag=12)
{

	n <- length(x)
	x.bar <- mean(x)
	
	gamma.hat <- numeric(max.lag+1)
	
	for(h in 0:min(max.lag,n-1))
	{
		
		gamma.hat[h+1] <- 0
		
		for(t in 1:(n-h))
		{
		
			gamma.hat[h+1] <- gamma.hat[h+1] + (x[t] - x.bar)*(x[t+h] - x.bar)
		
		}
		
	}
	
	gamma.hat <- gamma.hat / n
	rho.hat <- gamma.hat / gamma.hat[1]
	
	
	output <- list( gamma.hat = gamma.hat,
					rho.hat = rho.hat,
					lags = 0:max.lag)
					
	return(output)
	
}

#'Performs the Durbin-Levinson algorithm for one-step-ahead prediction.
#'
#' @param X a vector containing time series data.
#' @param gamma.0 the value of the autocovariance function at lag zero
#' @param gamma.n a vector containing the values of the autocovariance function at lags \code{1} through \code{length(X)}
#' @return a list containing the one-step-ahead predictions,the values of the partial autocorrelation function at lags \code{1} through \code{length(X)}, the MSPEs of the predictions, and the matrix \code{Phi}.
#' This function performs the Durbin-Levinson algorithm for one-step-ahead prediction.  Data are centered before the prediction are computed and then mean is added back to the predictions. 
DL.1step <- function(X, gamma.0, gamma.n) 
{
	
	X.cent <- X - mean(X)
	
    n <- length(X)
    alpha <- numeric(n)
    v <- numeric(n + 1)
    Phi <- matrix(0,n+1,n)
    
    v[1] <- gamma.0
    Phi[2,1] <- gamma.n[1]/gamma.0
    alpha[1] <- Phi[2,1]
    
    for (k in 1:(n - 1))
    {

        v[k + 1] <- v[k] * (1 - Phi[1+k,1]^2)
        Phi[1+k+1,1] <- (gamma.n[k + 1] - sum(gamma.n[1:k] * Phi[1+k,1:k]))/v[k + 1]
        Phi[1+k+1,(k+1):2]	<- Phi[1+k,k:1] - Phi[1+k+1,1] * Phi[1+k,1:k]
        alpha[k + 1] <- Phi[1+k+1,1]

    }
    
	v[n+1] <- v[n] * (1 - Phi[n+1,1]^2)
	
	X.pred <- as.numeric(Phi %*% X.cent) + mean(X)

    output <- list(	X.pred = X.pred, 
    					alpha = alpha, 
    					v = v, 
    					Phi = Phi)

    return(output)

}

#'Performs the innovations algorithm for h-step-ahead prediction.
#'
#' @param X a vector containing time series data.
#' @param h the number of steps ahead for which to make predictions.
#' @param K the covariance matrix of the random variables X_1,\dots,X_{n+h}.
#' @return a list containing the predicted values as well as the MSPEs of the predictions and the matrix \code{Theta}.
#' This function performs the innovations algorithm for one-step-ahead and h-step-ahead prediction. Data are centered before the prediction are computed and then mean is added back to the predictions. 
innov.hstep<- function(X,h,K){
		
		
	X.cent <- X - mean(X)	
	n <- length(X)
	v <- numeric(n+h)
	X.pred <- numeric(n+h)
	Theta <- matrix(0,n+h,n+h)

	v[1] <- K[1,1]
	X.pred[1] <- 0
	
	Theta[1+1,1] <-  K[2,1] / v[1]
	v[2] <- K[2,2] - Theta[1+1,1]^2*v[1]
	X.pred[2] <- Theta[1+1,1]*X[1]
	
	for(k in 2:n)
	{
		
		Theta[1+k,k] <- K[k+1,1] / v[1]
						
		for(j in 1:(k-1))
		{
			
			Theta[1+k,k-j] <- (K[k+1,j+1]-sum(Theta[1+j,j:1]*Theta[1+k,k:(k-j+1)]*v[1:j]))/v[j+1]
				
		}
		
		v[k+1] <- K[k+1,k+1] - sum( Theta[1+k,k:1]^2 * v[1:k] )
		
		X.pred[k+1] <- sum( Theta[1+k,1:k] *(X.cent[k:1] - X.pred[k:1]) )
		
	}

	if( h > 1)
	{
		
		for(k in (n+1):(n+h-1))
		{
			
			Theta[1+k,k] <- K[k+1,1] / v[1]
			
			for(j in 1:(k-1))
			{
				
				Theta[1+k,k-j]<-(K[k+1,j+1]-sum(Theta[1+j,j:1]*Theta[1+k,k:(k-j+1)]*v[1:j]))/v[j+1]
					
			}
		
			v[k+1] <- K[k+1,k+1] - sum( Theta[1+k,(k-n+1):k]^2 * v[n:1] )
			
			X.pred[k+1] <- sum( Theta[1+k,(k-n+1):k] *(X.cent[n:1] - X.pred[n:1]) )	
						
		}

	}

	for( k in 1:(n+h-1))
	{
		
		Theta[1+k,1:k] <- Theta[1+k,k:1] # switch order of indices for export
		
	}
		
		
	X.pred <- X.pred + mean(X)	 # add mean back to predictions
	
	output <- list( X.pred = X.pred,
					v = v,
					Theta = Theta[1:n,1:n])
					
	return(output)

}

#' Gives coefficients of causal representation of a causal ARMA(p,q) model.
#'
#' @param phi a vector with autoregressive coefficients.
#' @param theta a vector the moving average coefficients.
#' @param trun the number of terms to keep in the causal representation of the model
#' @return a vector containing the first \code{trun+1} moving average coefficients in the causal representation of the ARMA(p,q) model
#' This function recursively computes the coefficients in the MA(Inf) representation of the causal ARMA(p,q) model.
ARMAtoMAinf <- function(phi=NULL,theta=NULL,trun=500)
{
   
  	if(length(phi)==0)
  	{
  		q <- length(theta)	
  		psi <- numeric(trun+1)
  		psi[1:(q+1)] <- c(1,theta)
  		
  	} else if(length(phi)>0)
  	{
  		# check to see if the time series is causal:
		minroot <- min(Mod(polyroot(c(1,-phi))))
		if( minroot < 1)
			stop("The ARMA process specified is not causal.")
		
		p <- length(phi)
		q <- length(theta)
		
		# set theta_j = 0 for j > q
		theta.0 <- c(theta,rep(0,trun-q))
		
		# set psi_j = 0 for j < 0
		psi.0 <- numeric(trun+p)
		psi.0[p] <- 1 # this is psi_0
		
		for(j in 1:trun)
		{
			
			psi.0[p+j] <- theta.0[j] + sum( phi[1:p] * psi.0[(p+j-1):j] )
			
		}	
		
		# take away zeroes at beginning
	  	psi <- psi.0[p:(p+trun)]
	
	}
	
	return(psi)
	
}

#' Generates data from an ARMA(p,q) model with iid Normal innovations.
#'
#' @param phi a vector with autoregressive coefficients.
#' @param theta a vector the moving average coefficients.
#' @param sigma the white noise variance.
#' @param n the length of the realization to be generated.
#' @param trun the number of terms to keep in the causal representation of the model, which is used to generate the data.
#' @return a length-\code{n} realization of the time series.
#'
#' This function generates a length-\code{n} realization from a causal invertible ARMA(p,q) model with iid Normal innovations.
get.ARMA.data <- function(phi=NULL,theta=NULL,sigma=1,n,trun=500)
{
	
	# check to see if the time series is causal:
	if(length(phi) > 0)
	{
	minroot <- min(Mod(polyroot(c(1,-phi))))
	if( minroot < 1)
		stop("The ARMA process specified is not causal.")
	}
	# check to see if the time series is invertible
	if(length(theta)>0)
	{
	minroot <- min(Mod(polyroot(c(1,theta))))
	if( minroot < 1)
		stop("The ARMA process specified is not invertible.")
	}
	
	psi <- ARMAtoMAinf(phi,theta,trun)	
	
	Z <- rnorm(n+trun,0,sigma)
	X <- numeric(n)
	
	for( t in 1:n)
	{
		ind <- trun + t:(t-trun)
		X[t] <- sum( psi * Z[ind] )
	}

	return(as.ts(X))

}

#' Computes the autocovariance function of a causal ARMA(p,q) process.
#'
#' @param phi a vector with autoregressive coefficients.
#' @param theta a vector the moving average coefficients.
#' @param sigma the white noise variance.
#' @param max.lag the number of lags at which to return the value of the autocovariance function.
#' @param trun the order at which to truncate the MA(Inf) representation of the causal ARMA(p,q) process.
#' @return A vector containing the values of the autocovariance function at lags \code{0} through \code{max.lag}.
ARMAacvf <- function(phi=NULL,theta=NULL,sigma=1,max.lag=12,trun=500)
{
	# check to see if the time series is causal:
	if(length(phi)>0)
	{
	minroot <- min(Mod(polyroot(c(1,-phi))))
	if( minroot < 1)
		stop("The ARMA process specified is not causal.")
	}
	
	psi <- ARMAtoMAinf(phi,theta,trun)
		
	gamma.0 <- sigma^2 * sum(psi^2)
	ARMAacvf <- numeric(max.lag+1)
	ARMAacvf[1] <- gamma.0
	
	for(h in 1:max.lag)
	{
		
		ARMAacvf[1+h] <- sigma^2 * sum( psi[1:(trun-h)] * psi[(1+h):trun])
		
	}
	
	return(ARMAacvf)
	
}

#' Computes h-step-ahead predictions from an ARMA(p,q) model
#'
#' @param X a vector containing time series data.
#' @param h the number of steps ahead for which to make predictions.
#' @param phi a vector with autoregressive coefficients.
#' @param theta a vector the moving average coefficients.
#' @param sigma the white noise variance.
#' @return a list containing the predicted values as well as the MSPEs of the predictions and the AIC and BIC.
#' This function builds a matrix of autocovariances for the ARMA(p,q) model using the MA(inf) representation of the process. It then runs the innovations algorithm on this matrix of autocovariances.
ARMA.hstep <- function(X,h,phi,theta,sigma)
{
	
	X.cent <- X - mean(X)
	n <- length(X)
	gamma.hat <- ARMAacvf(phi,theta,sigma,max.lag=n+h)
	gamma.0 <- gamma.hat[1]
	gamma.n <- gamma.hat[-1]
	
	K <- matrix(0,n+h,n+h)
	for(j in 1:(n+h-1))
		for(i in (j+1):(n+h))
		{
			
			K[i,j] <- c(gamma.0,gamma.n)[1+abs(i-j)]
			
		}
	
	K <- K + t(K) + diag(rep(gamma.0),n+h)
	
	innov.hstep.out <- innov.hstep(X.cent,h,K) # this part is inefficient. Speed up someday with 5.3.9 of B&D Theory
	X.pred <- innov.hstep.out$X.pred + mean(X)
	v <- innov.hstep.out$v
	
	p <- length(phi)
	q <- length(theta)
	ll <- -(n/2)*log(2*pi) - (1/2) * sum( log( v[1:n] )) - (1/2) * sum( (X - X.pred[1:n])^2/v[1:n])
	aic <- -2*ll + 2 * ( p + q + 1 ) # plus 1 for the variance 
	bic <- -2*ll + log(n) * ( p + q + 1 )
		
	output <- list(	X.pred = X.pred,
					v = v,
					aic = aic,
					bic = bic)
					
	return(output)
	
}



#' Find the spectral density function of an ARMA(p,q) process.
#'
#' @param phi a vector with autoregressive coefficients.
#' @param theta a vector the moving average coefficients.
#' @param sigma the white noise variance.
#' @param plot if \code{TRUE} the spectral density is plotted
#' @return a list containing values of the spectral density function and the corresponding frequencies
ARMAtoSD <- function(phi=NULL,theta=NULL,sigma,plot=TRUE)
{

	lambda <- seq(-pi,pi,by=2*pi/1e5)
	
	p <- length(phi)
	q <- length(theta)
	
	theta.pol <- rep(1,length(lambda))
	if(q > 0)
	{
		for(k in 1:q)
		{
			theta.pol <- theta.pol + theta[k] * exp(-1i*k*lambda)
		}
	}
	
	phi.pol <- rep(1,length(lambda))
	if( p > 0)
	{
		for(k in 1:p)
		{
			phi.pol <- phi.pol - phi[k] * exp(-1i*k*lambda)
		}
	}
	
	f <- sigma^2/(2*pi) * Mod(theta.pol)^2/Mod(phi.pol)^2
	
	if(plot==TRUE)
	{
		plot(f[lambda>=0]~lambda[lambda>=0],type="l",ylab="spectral density",xlab="lambda")
	}
	
	output <- list(	f = f,
			lambda = lambda)
					
	return(output)
	
}

#' Find the moving average representation of a time series based on the spectral density.
#'
#' @param f a vector with evaluations of the spectral density
#' @param lambda the frequencies to which the values of f correspond. Should be evenly spaced and dense between \code{-pi} and \code{pi}.
#' @param trun the number of terms to keep in the moving average representation of the time series
#' @param tol controls the accuracy. Smaller values require more computation time.
#' @return a vector containing the coefficients of the moving average representation of the time series.
#' This is based on the work of Pourahmadi (1984).
#' @references Pourahmadi, M. (1984). Taylor expansion of and some applications. \emph{The American Mathematical Monthly}, 91(5), 303-307.
SDtoMAinf <- function(f,lambda,trun=500,tol=1e-4)
{
		
	delta <- 2*pi/(length(lambda)-1) # assume equally spaced lambdas
	a <- numeric(trun)
	for(k in 1:trun)
	{
		# integrate over the log of the spectral density
		a[k] <- 1/(2*pi) * sum( log(f) * exp(-1i*(k-1)*lambda)) * delta
		
		if(abs(a[k]) < tol) break
		
	}
	
	c <- numeric(trun)
	c[1] <- 1
	for(k in 0:(trun-2))
	{
		
		c[k+2] <- sum( (1 - c(0:k) / (k+1) ) * a[k:0 + 2] * c[0:k + 1] )	
	
	}
	
	c <- ifelse(abs(Re(c)) > tol,Re(c),0)
	
	return(c)
	
}

#' Compute the periodogram.
#'
#' @param X a vector of data
#' @param plot if TRUE a plot of the periodogram is generated
#' @return a list containing the periodogram ordinates, the Fourier frequencies, and the cumulative periodogram
#'
#' This function computes the periodogram at the Fourier frequencies.
pgram <- function(X,plot=FALSE)
{
	
	n <- length(X)
	lambda <- (-floor((n-1)/2):floor(n/2))/n*2*pi
	
	E <- matrix(NA,n,n)
	for(i in 1:n)
	{
		E[i,] <- 1/sqrt(n) * exp(1i*i*lambda)
	}		
	
	# compute discrete Fourier transform of X		
	D <- as.vector(t(Conj(E)) %*% X)
	I <- Mod(D)^2
	
	if(plot == TRUE)
	{
		plot(I[lambda>0]~lambda[lambda>0],type="o")
	}
	
	# compute cumulative periodogram	
	Y0 <- cumsum(I[lambda < 0]) / sum(I[lambda < 0])
	Y <- Y0[-length(Y0)]	

	output <- list( I = I,
					lambda = lambda,
					Y = Y)
	
	return(output)
	
}


#' Perform Bartlett's test for whether a time series is iid Normal.
#'
#' @param X a vector containing time series data.
#' @param plot if TRUE a plot of the cumulative periodogram is generated, on which the Kolmogorov-Smirnov bounds are displayed.
#' @return the p value of Bartlett's test
#'
#' This function performs Bartlett's test for whether the time series is iid Normal.
WNtest.Bartlett <- function(X,plot=FALSE)
{
	
	Y <- pgram(X)$Y
	q <- length(Y) + 1
	
	# compute test statistic
	dev <-sqrt(q-1) * abs(Y - c(1:(q-1))/(q-1)) 
	B <- max(dev)
	
	# get p-value
	j <- c(-100:100)
	pval <- 1 - sum( (-1)^j * exp(-2*B^2*j^2) )
		
	if(plot == TRUE)
	{
		plot(Y,
				main=paste("Bartlett test: p val = ",round(pval,3),sep=""),
				xlab="frequencies",
				ylab="cumulative periodogram",
				col = ifelse( dev > 1.36,"red","black" ),
				pch = ifelse( dev > 1.36,19,1 ))
		abline(-1.36/sqrt(q-1),1/(q-1))
		abline(+1.36/sqrt(q-1),1/(q-1))
		abline(-1.63/sqrt(q-1),1/(q-1),lty=3)
		abline(+1.63/sqrt(q-1),1/(q-1),lty=3)

	}
	
	return(pval)
	
}

#' Perform the Ljung-Box test for whether a time series is iid.
#'
#' @param X a vector containing time series data.
#' @param k number of lags on which to base the test statistic
#' @param nparms number of parameters in the fitted model; the degrees of freedom of the chi-squared distribution from which the p value is computed will be set equal to \code{k - nparms}.
#' @return the p value of the Ljung-Box test
#'
#' This function performs the Ljung-Box test for whether the time series is iid.
WNtest.LB <- function(X,k=1,nparms=0)
{
	
	n <- length(X)
	
	if(k > n-1 ) 
	  stop("k must be less than n")
	if(k <= nparms)
	  stop("k must be greater than nparms")
	
	rho.k <- sample.acf(X)$rho.hat[2:(k+1)]
	Q.LB <- n*(n+2) * sum( rho.k^2 / (n - 1:k))
	pval <- 1 - pchisq(Q.LB,k-nparms)
	
	return(pval)
	
}

#' Perform the Dickey-Fuller unit-root test
#'
#' @param X a vector containing time series data.
#' @param p the autoregressive order on which to base the test.
#' @return a list with logicals indicating rejection or failure to reject at the 0.05 and 0.01 significance levels as well as the value of the Dickey-Fuller test statistic.
#'
#' This function performs the Dickey-Fuller unit-root test.
unitroottest.DF <- function(X,p=1)
{
	
	if(p==1)
	{
		
		Xdiff <- diff(X)
		Xlag <- X[1:(n-1)]
		lm.out <- lm(Xdiff ~ Xlag)
		DF <- summary(lm.out)$coef[2,3]
		
	} else if(p > 1){
		
    Xdiff <- diff(X)
	  
    XX <- matrix(NA,n-p,p)
    XX[,1] <- X[p:(n-1)]
	  
    for(j in 1:(p-1))
	  {
	    
	    XX[,1+j] <- Xdiff[(p-j):(n-1-j)]
	    
	  }
	 
	  lm.out <- lm(Xdiff[p:(n-1)] ~ XX)
		DF <- summary(lm.out)$coef[2,3]
	  	
	}
	
	reject.01 <- FALSE
	reject.05 <- FALSE
	if(DF < -3.43)
	{
		reject.01 <- TRUE
	}
	if(DF < -2.86)
	{
		reject.05 <- TRUE
		
	}
	
	output <- list( reject.05 = reject.05,
					        reject.01 = reject.01,
					        DF = DF)
					
	return(output)
	
}

#' Evaluate the Parzen window
#'
#' @param x a real number.
#' @return the value of the Parzen window function.
parzen <- function(x)
{
	
	if( abs(x) < 1/2){
		
		return( 1 - 6 * x^2 + 6 * abs(x)^3 )
		
	} else if( (abs(x) >= 1/2) & (abs(x) <= 1)){
		
		return( 2*(1 - abs(x))^3 )
		
	} else {
		
		return( 0 )
		
	}
	
}

#' Compute a lag-window estimator of the spectral density
#'
#' @param X a numeric vector containing the observed data
#' @param L the number of lags at which to truncate the autocovariance function.
#' @param window the lag window function to be used. Default is the Parzen window.
#' @param nlambda the number of values between \code{-pi} and \code{pi} for which the estimate of the spectral density should be computed.
#' @param plot if \code{TRUE} then a plot of the estimated spectral density is produced.
#' @return a list containing values of the estimated spectral density function and the corresponding frequencies.
SDlagWest <- function(X,L,window=parzen,nlambda=1e4,plot=TRUE)
{
	
	lambda <- seq(-pi,pi,by=2*pi/nlambda)

	gamma.hat.0toL <- sample.acf(X,max.lag = L)$gamma.hat
	gamma.hat <- c(gamma.hat.0toL[(L+1):2],gamma.hat.0toL)
	
	E.lambda <- matrix(NA,2*L+1,length(lambda))
	gamma.hat.weighted <- numeric(2*L+1)
	for(j in (-L):L)
	{
		
		E.lambda[j+L+1,] <- exp( -1i * j * lambda )
		gamma.hat.weighted[j+L+1] <- window(abs(j/L)) * gamma.hat[j+L+1]
		
	}
	
	f.hat.L <- as.numeric( Re(t(Conj(E.lambda)) %*% gamma.hat.weighted )) / ( 2*pi )
	
	if(plot==TRUE)
	{
		plot(f.hat.L[lambda>=0]~lambda[lambda>=0],type="l",ylab="estimate of spectral density",xlab="lambda")
	}
	
	output <- list( f.hat = f.hat,
					lambda = lambda)
					
	return(output)
	
}

