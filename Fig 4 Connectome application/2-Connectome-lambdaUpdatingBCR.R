# AMP algorithm which takes MLE for sigma.sq based on current parameter estimate in each iteration, i.e. sample variance of estimated residuals
BGAMP <- function(y, A, psi, theta) {

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
	print(paste("Counter:", counter))
	if(counter == n.iter-1) print("Warning: Maximum number of iterations reached")
	
	# Step 1
	for(a in 1:m) {
		z[k,a] <- A[a,]%*%x.hat[k,]
		tau.p[k,a] <- A[a,]^2%*%tau.x[k,]
		p.hat[k,a] <- z[k,a]-tau.p[k,a]*s[k,a]
	}
	
	#print(paste("z:", z[k,1]))
	#print(z[k,m])
	
	# E-M step for sigma.sq
	#sigma.sq <- var(y-z[k,])
	
	# Set sigma.sq equal to rough estimate of posterior mode
	sigma.sq <- (2*b.0+sum((y-z[k,])^2))/(2*a.0+m)
	print(paste("sigma.sq:", sigma.sq))
	
	# Step 2
	for(a in 1:m) {
		tau.z[k,a] <- tau.p[k,a]*sigma.sq/(tau.p[k,a]+sigma.sq) # Var(z_a|y,p.hat[k,a],tau.p[k,a])
		z.hat <- tau.z[k,a]*(y[a]/sigma.sq+p.hat[k,a]/tau.p[k,a]) # E(z.a|y,p.hat[k,a],tau.p[k,a])
		
		# Dampened:
		s[k+1,a] <- (1-theta)*s[k,a] + theta*(z.hat-p.hat[k,a])/tau.p[k,a]
		tau.s[k,a] <- (1-tau.z[k,a]/tau.p[k,a])/tau.p[k,a]
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
	print(paste("lambda:", lambda[k+1]))
	
	if(sd(x.hat[k,]-x.hat[k+1,]) < 1E-13) break # Stop iterating if AMP has sufficiently converged
}

return(list(mean=x.hat[counter,], var=tau.x[counter,], incProb=pi[counter,]))

}




# Log of marginal likelihood with sigma.sq ~ IG(a.0, b.0), beta ~ N(0, sigma.sq*Sigma.beta) and the proportionality constant (2pi)^(-n/2) omitted
l.marginal.likelihood <- function(y, X, Sigma.beta, a.0, b.0) {
	Lambda.n <- t(X)%*%X+solve(Sigma.beta)
	mu.n <- solve(Lambda.n)%*%t(X)%*%y
	
	a.n <- a.0+length(y)/2
	b.n <- as.numeric(b.0+( t(y)%*%y - t(mu.n)%*%Lambda.n%*%mu.n )/2)
	
	return( a.0*log(b.0)+lgamma(a.n) - (log(det(Sigma.beta))+log(det(Lambda.n)))/2 - a.n*log(b.n)-lgamma(a.0) )
}


library(R.utils) # for intToBin()
# Log of marginal likelihood with spike-and-slab prior
l.marginal.likelihood.BG <- function(y, X, lambda, psi) {
	p = length(X[1,]) # Number of predictors
	n = length(y) # Sample size
	
	a.n <- a.0+n/2
	l.sum <- p*log(1-lambda) + a.0*log(b.0) + lgamma(a.n) - a.n*log(b.0+as.numeric(t(y)%*%y)/2 ) - lgamma(a.0) # The marginal likelihood of the model without predictors is proportional to this
	
	for(k in 1:(2^p-1)) {
		gam <- rep(0, p)
		gam[min(p-ceiling(log2(k+1)-1),p):p] <- as.numeric(strsplit(intToBin(k), split = "")[[1]])
		X.gam <- X[,which(gam==1)]
		s.gam <- sum(gam)
		l.add <- s.gam*log(lambda)+ (p-s.gam)*log(1-lambda) + l.marginal.likelihood(y, X.gam, psi*diag(s.gam), a.0, b.0)
		l.sum <- l.sum + log(1+exp(l.add-l.sum))
	}
	
	return(l.sum)
}


exact <- function(y, A, lambda, psi) {
	incProb <- NULL
	l.marginalLikelihood <- l.marginal.likelihood.BG(y, A, lambda, psi)
	
	for(j in 1:length(A[1,])) {
		l.marginalLikelihoodj <- l.marginal.likelihood.BG(y, A[,-j], lambda, psi)
		incProb[j] <- 1-(1-lambda)*exp(l.marginalLikelihoodj-l.marginalLikelihood)
	}
	
	return(list(incProb=incProb))
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
		print(r)
		theta <- runif(1, min = .1, max = .9)
		Theta <- matrix(0, nrow = p, ncol = m)
		for(i in 1:p) for(j in 1:m) {
			u <- runif(1)
			if(u < theta^2)						Theta[i,j] <- -sqrt(1/theta)
			else if(u < 1-2*theta+2*theta^2)	Theta[i,j] <- sqrt(1/theta)
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



screeningMethod.BCR <- function(y, A, lambda, psi) { # No screening this time as n > p
	n = length(A[,1]) # Sample size
	p = length(A[1,]) # Number of predictors
	
	mean <- rep(0,p)
	incProb <- rep(0,p)
	
	p.star <- round(3*lambda*p)
	if(p >= length(y) & p.star < p) {
		# (Ridge) HOLP screening (Sam Wang) to get superset of three times the average size of nonzero predictors
		x.hat <- t(A)%*%solve(diag(n)+A%*%t(A))%*%y
		supset <- order(-abs(x.hat))[1:p.star]
		lambda.star <- lambda*(p/p.star)
	}
	else {
		supset <- 1:p
		lambda.star <- lambda
	}
	
	incProb[supset] <- BCR(y, A[,supset], lambda.star, psi, diag(x = psi, round(p.star*lambda.star)))$incProb
	
	return(list(incProb=incProb))
}



BCR.2 <- function(y, X, psi, m) {
	p <- length(X[1,]) # Number of predictors
	k.max <- 20 # Number of iterations
	
	# Initialize lambda (prior inclusion probability)
	lambda <- NULL
	lambda[1] <- a/(a+b) # Prior mean
	
	for(k in 1:k.max) {
		print(paste("progress: ", k/k.max))
		# Update lambda with a rough estimate of the posterior mean
		lambda[k+1] <- (a+sum(BCR(y, X, lambda[k], psi, diag(x = psi, m))$incProb))/(a+b+p)
		print(paste("lambda:", lambda[k+1]))
	}
	print(lambda)
	return(BCR(y, X, lambda[k.max+1], psi, diag(x = psi, m)))
}




load("connectome_data")
y <- scale(y) # Standardize y
n <- length(y)

subset <- NULL
i=1
for(j in 1:70^2) {
	if(!all(connect[,j] == rep(0,n))) {
		subset[i] <- j
		i <- i+1
	}
}

con <- scale(connect[,subset]) # Demean columns and make L2 column norm equal to one

# Hyperparameters for prior on sigma.sq
a.0 <- 3
b.0 <- 1


# Hyperparameters for beta prior on lambda
a <- 2
b <- 5



BCRincProb <- BCR.2(y, con, 1, 20)$incProb
save(BCRincProb, file="BCRResult")