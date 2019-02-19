#'Computes the sample autocovariance and autocorrelation functions.
#'
#' @param x a vector containing time series data.
#' @param max.lag the largest lag at which to compute the autocovariance and autocorrelation functions.
#' @return a list containing the autocovariance at lags zero through \code{max.lag}
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
#' @return a list containing the one-step-ahead predictions,the values of the partial autocorrelation function at lags \code{1} through \code{length(X)}, the MSPEs of the predictions, and the matrix \code{Phi} for which \code{Phi%*%X} returns the predictions.
#' This function performs the Durbin-Levinson algorithm for one-step-ahead prediction
DL.1step <- function(X, gamma.0, gamma.n) 
{
	
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
	
	X.pred <- as.numeric(Phi %*% X)

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
#' @return a list containing the predicted values as well as the MSPEs of the predictions and the matrix \code{Theta} such that \code{Theta %*% (X - X.pred) = X.pred}.
#' This function performs the innovations algorithm for one-step-ahead and h-step-ahead prediction
innov.hstep<- function(X,h,K){
	  
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
		
		X.pred[k+1] <- sum( Theta[1+k,1:k] *(X[k:1] - X.pred[k:1]) )
		
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
			
			X.pred[k+1] <- sum( Theta[1+k,(k-n+1):k] *(X[n:1] - X.pred[n:1]) )	
						
		}

	}

	for( k in 1:(n+h-1))
	{
		
		Theta[1+k,1:k] <- Theta[1+k,k:1] # switch order of indices for export
		
	}
	
	
	output <- list( X.pred = X.pred,
					v = v,
					Theta = Theta[1:n,1:n])
					
	return(output)

}

