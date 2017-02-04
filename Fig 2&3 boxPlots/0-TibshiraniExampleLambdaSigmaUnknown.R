# AMP algorithm which takes MLE for sigma.sq based on current parameter estimate in each iteration, i.e. sample variance of estimated residuals
BGAMP <- function(y, A, psi, theta = 1) {

# Set hyperparameters for the Bernoulli-Gaussian / spike-and-slab prior
dzeta = 0 # Mean of the normal


# Initialize lambda (prior inclusion probability)
lambda <- NULL
lambda[1] <- a/(a+b) # Prior mean


# Initialization for the AMP algorithm
n.iter = 10000 # Maximum number of iterations for AMP algorithm
n = length(A[1,]) # Number of predictors
m = length(A[,1]) # Sample size
x.hat = matrix(lambda[1]*dzeta, nrow = n.iter, ncol = n) # Initialize vector of parameter estimates at prior mean
tau.x = matrix(lambda[1]^2*psi, nrow = n.iter, ncol = n) # Initialize vector of variances at prior variance
s = matrix(0, nrow = n.iter, ncol = m) # Note that row indexing starts at 0 in Andersen's thesis
z = matrix(0, nrow = n.iter, ncol = m)
tau.p = matrix(0, nrow = n.iter, ncol = m)
p.hat = matrix(0, nrow = n.iter, ncol = m)
tau.z = matrix(0, nrow = n.iter, ncol = m)
tau.s = matrix(0, nrow = n.iter, ncol = m)
tau.r = matrix(0, nrow = n.iter, ncol = n)
r.hat = matrix(0, nrow = n.iter, ncol = n)
pi = matrix(0, nrow = n.iter, ncol = n)
counter <- 0

# The AMP algorithm
for(k in 1:(n.iter-1)) {
	counter <- counter+1
	#print(paste("Counter:", counter))
	if(counter == n.iter-1) print("Warning: Maximum number of iterations reached")
	
	# Step 1
	for(aa in 1:m) {
		z[k,aa] <- A[aa,]%*%x.hat[k,]
		tau.p[k,aa] <- A[aa,]^2%*%tau.x[k,]
		p.hat[k,aa] <- z[k,aa]-tau.p[k,aa]*s[k,aa]
	}
	
	#print(paste("z:", z[k,1]))
	#print(z[k,m])
	
	# E-M step for sigma.sq
	#sigma.sq <- var(y-z[k,])
	
	# Set sigma.sq equal to rough estimate of posterior mode
	sigma.sq <- (2*b.0+sum((y-z[k,])^2))/(2*a.0+m)
	#print(paste("sigma.sq:", sigma.sq))
	
	# Step 2
	for(aa in 1:m) {
		tau.z[k,aa] <- tau.p[k,aa]*sigma.sq/(tau.p[k,aa]+sigma.sq) # Var(z_a|y,p.hat[k,a],tau.p[k,a])
		z.hat <- tau.z[k,aa]*(y[aa]/sigma.sq+p.hat[k,aa]/tau.p[k,aa]) # E(z.a|y,p.hat[k,a],tau.p[k,a])
		
		# Dampened:
		s[k+1,aa] <- (1-theta)*s[k,aa] + theta*(z.hat-p.hat[k,aa])/tau.p[k,aa]
		tau.s[k,aa] <- (1-tau.z[k,aa]/tau.p[k,aa])/tau.p[k,aa]
	}
	
	# Step 3
	for(i in 1:n) {
		tau.r[k,i] <- 1/A[,i]^2%*%tau.s[k,]
		r.hat[k,i] <- x.hat[k,i]+tau.r[k,i]*A[,i]%*%s[k+1,]
	}
	#print(paste("tau.r:", tau.r[k,1]))
	#print(tau.r[k,n])
	
	#print(paste("r.hat:", r.hat[k,1]))
	#print(r.hat[k,n])
	
	# Step 4
	for(i in 1:n) {
		pi[k,i] <- 1/(1+(1-lambda[k])*dnorm(0, r.hat[k,i], sqrt(tau.r[k,i]))/(lambda[k]*dnorm(0, dzeta-r.hat[k,i], sqrt(psi+tau.r[k,i]))))
		nu <- 1/(1/psi+1/tau.r[k,i])
		gamma <- (dzeta/psi+r.hat[k,i]/tau.r[k,i])*nu
		
		# Dampened:
		x.hat[k+1,i] <- (1-theta)*x.hat[k+1,i] + theta*pi[k,i]*gamma # E(x_i|y,r.hat[k,i],tau.r[k,i])
		tau.x[k+1,i] <- pi[k,i]*(gamma^2-pi[k,i]*gamma^2+nu) # Var(x_i|y,r.hat[k,i],tau.r[k,i])
	}
	
	# Update lambda with a rough estimate of the posterior mean
	lambda[k+1] <- (a+sum(pi[k,]))/(a+b+n)
	#print(paste("lambda:", lambda[k+1]))
	
	if(is.na(sd(x.hat[k,]-x.hat[k+1,]) < 1E-9)) return(BGAMP(y, A, psi, theta/2)) # If AMP fails to converge, add dampening
	else if(sd(x.hat[k,]-x.hat[k+1,]) > 1) {
		print("Hey!")
		return(BGAMP(y, A, psi, theta/2))
	}
	else if(sd(x.hat[k,]-x.hat[k+1,]) < 1E-9) break # Stop iterating if AMP has sufficiently converged
}

return(list(mean=x.hat[counter,], var=tau.x[counter,], incProb=pi[counter,], lambda=lambda[counter], sigma.sq=sigma.sq, theta=theta))

}


# Log of marginal likelihood with sigma.sq ~ IG(a.0, b.0), beta ~ N(0, sigma.sq*Sigma.beta) and the proportionality constant (2pi)^(-n/2) omitted
l.marginal.likelihood <- function(y, X, Sigma.beta, a.0, b.0) {
	Lambda.n <- t(X)%*%X+solve(Sigma.beta)
	mu.n <- solve(Lambda.n)%*%t(X)%*%y
	
	a.n <- a.0+length(y)/2
	b.n <- as.numeric(b.0+( t(y)%*%y - t(mu.n)%*%Lambda.n%*%mu.n )/2)
	
	return( a.0*log(b.0)+lgamma(a.n) - (log(det(Sigma.beta))+log(det(Lambda.n)))/2 - a.n*log(b.n)-lgamma(a.0) )
}


library(pracma) # For gramSchmidt
library(Matrix) # For bdiag
BCR <- function(y, X, lambda, psi, Sigma.beta) {
	p <- length(X[1,]) # Number of predictors
	m <- length(Sigma.beta[1,])
	
	rep <- 10 # Number of BCR models averaged
	incProb.rep <- matrix(nrow = rep, ncol = p)
	l.rep <- matrix(nrow = rep, ncol = p)
	l.sum <- NULL
	
	for(r in 1:rep) {
		#print(r)
		theta <- runif(1, min = .1, max = .9)
		Theta <- matrix(0, nrow = p, ncol = m)
		while(qr(Theta)$rank < m) {
			for(i in 1:p) for(j in 1:m) {
				u <- runif(1)
				if(u < theta^2)						Theta[i,j] <- -sqrt(1/theta)
				else if(u < 1-2*theta+2*theta^2)	Theta[i,j] <- sqrt(1/theta)
			}
		}
		Theta <- gramSchmidt(Theta)$Q
		
		for(j in 1:p) {
			l.0 <- l.marginal.likelihood(y, X[,-j]%*%Theta[-j,], Sigma.beta, a.0, b.0)
			l.j <- l.marginal.likelihood(y, cbind(X[,-j]%*%Theta[-j,], X[,j]), bdiag(Sigma.beta, psi), a.0, b.0)
			
			incProb.rep[r,j] <- 1 / ( 1+( (1-lambda) / lambda ) * exp( l.0 - l.j ) )
			
			l.rep[r,j] <- log(1-lambda)+l.0+log(1+exp(log(lambda)+l.j-log(1-lambda)-l.0))
			if(r == 1) l.sum[j] <- l.rep[1,j]
			else l.sum[j] <- l.sum[j]+log(1+exp(l.rep[r,j]-l.sum[j]))
		}
	}
	
	model.weights <- exp(t(l.rep)-l.sum)
	# print(model.weights) # Print this to check whether specific model dominates
	
	incProb <- NULL
	for(j in 1:p) incProb[j] <- model.weights[j,]%*%incProb.rep[,j]
	
	return(list(incProb=incProb))
}


BCR.2 <- function(y, X, psi, m = 4) {
	p <- length(X[1,]) # Number of predictors
	k.max <- 5 # Number of iterations
	
	# Initialize lambda (prior inclusion probability)
	lambda <- NULL
	lambda[1] <- a/(a+b) # Prior mean
	
	for(k in 1:k.max) {
		# print(paste("progress: ", k/k.max))
		# Update lambda with a rough estimate of the posterior mean
		lambda[k+1] <- (a+sum(BCR(y, X, lambda[k], psi, diag(x = psi, m))$incProb))/(a+b+p)
		#print(paste("lambda:", lambda[k+1]))
	}
	#print(lambda)
	return(list(incProb=BCR(y, X, lambda[k.max+1], psi, diag(x = psi, m))$incProb, lambda=lambda))
}



# Generate data similar to first example in Tibshirani's lasso paper
gendata <- function(n, p, rho, SNR) {
	beta <- rep(0,p)
	beta[1] <- 3
	beta[2] <- 1.5
	beta[3] <- 2
	
	sigma.sq <- sum(beta^2)/SNR
	
	X <- matrix(rnorm(n*p), nrow = n, ncol = p)
	for(j in 2:p) X[,j] <- rho*X[,j-1]+sqrt(1-rho^2)*X[,j]
	X <- scale(X)
	
	y <- X%*%beta+rnorm(n, mean = 0, sd = sqrt(sigma.sq))
	
	lambda <- 3/p
	
	return(list(beta=beta, lambda=lambda, X=X, y=scale(y), a=7*lambda, b=7*(1-lambda)))
}


# Hyperparameters for prior on sigma.sq
a.0 <- 3
b.0 <- 1


psi <- 10
n <- 100
p <- 7

repeats <- 200
AMPincProb <- matrix(nrow = repeats, ncol = p)
BCRincProb <- matrix(nrow = repeats, ncol = p)

pb <- txtProgressBar(min = 0, max = repeats, style = 3)

for(rho in c(0, .2, .8, .5, .7)) {
	if(is.element(rho, c(0, .2))) SNR <- 1 else SNR <- 10
	
	setTxtProgressBar(pb, 0)
	for(r in 1:repeats) {
	
	d <- gendata(n,p,rho,SNR)
	
	# Hyperparameters for beta prior on lambda
	a <- d$a
	b <- d$b
	
	AMPresult <- BGAMP(d$y, d$X, psi)
	BCRresult <- BCR.2(d$y, d$X, psi)
	AMPincProb[r,] <- AMPresult$incProb
	BCRincProb[r,] <- BCRresult$incProb
	
	setTxtProgressBar(pb, r)
	}
	
	save(AMPincProb, BCRincProb, file = paste0('rho', 10*rho, 'SNR', SNR))
}