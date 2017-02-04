BGAMP <- function(y, A, lambda, psi, sigma.sq, theta = 1) {

# Set hyperparameters for the Bernoulli-Gaussian / spike-and-slab prior
dzeta = 0 # Mean of the normal

# Initialization for the AMP algorithm
n.iter = 10000 # Maximum number of iterations for AMP algorithm
n = length(A[1,]) # Number of predictors
m = length(A[,1]) # Sample size
x.hat = matrix(lambda*dzeta, nrow = n.iter, ncol = n) # Initialize vector of parameter estimates at prior mean
tau.x = matrix(lambda^2*psi, nrow = n.iter, ncol = n) # Initialize vector of variances at prior variance
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
	if(counter == n.iter-1) print("Warning: Maximum number of iterations reached")
	
	# Step 1
	for(a in 1:m) {
		z[k,a] <- A[a,]%*%x.hat[k,]
		tau.p[k,a] <- A[a,]^2%*%tau.x[k,]
		p.hat[k,a] <- z[k,a]-tau.p[k,a]*s[k,a]
	}
	
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
	
	# Step 4
	for(i in 1:n) {
		pi[k,i] <- 1/(1+(1-lambda)*dnorm(0, r.hat[k,i], sqrt(tau.r[k,i]))/(lambda*dnorm(0, dzeta-r.hat[k,i], sqrt(psi+tau.r[k,i]))))
		nu <- 1/(1/psi+1/tau.r[k,i])
		gamma <- (dzeta/psi+r.hat[k,i]/tau.r[k,i])*nu
		
		# Dampened:
		x.hat[k+1,i] <- (1-theta)*x.hat[k+1,i] + theta*pi[k,i]*gamma # E(x_i|y,r.hat[k,i],tau.r[k,i])
		
		tau.x[k+1,i] <- pi[k,i]*(gamma^2-pi[k,i]*gamma^2+nu) # Var(x_i|y,r.hat[k,i],tau.r[k,i])
	}
	
	if(is.na(sd(x.hat[k,]-x.hat[k+1,]) < 1E-9)) return(BGAMP(y, A, lambda, psi, sigma.sq, theta/2)) # If AMP fails to converge, add dampening
	else if(sd(x.hat[k,]-x.hat[k+1,]) < 1E-9) break # Stop iterating if AMP has sufficiently converged
}

return(list(mean=x.hat[counter,], var=tau.x[counter,], incProb=pi[counter,]))

}



l.marginal.likelihood.lm <- function(y, X, psi, sigma.sq) {
	if(length(dim(X)) == 2) Phi <- t(X)%*%X + (sigma.sq/psi)*diag(length(X[1,]))
	else Phi <- t(X)%*%X + sigma.sq/psi
	return( ( t(y)%*%X%*%solve(Phi)%*%t(X)%*%y / (2*sigma.sq) )  - ( length(X[1,])* (log(psi) - log(sigma.sq) ) + log(det(Phi)) )/2 )
}


library(pracma) # For gramSchmidt
library(Matrix) # For bdiag
BCR <- function(y, X, lambda, psi, sigma.sq, m = 3) {
	p <- length(X[1,]) # Number of predictors
	
	rep <- 10 # Number of BCR models averaged
	incProb.rep <- matrix(nrow = rep, ncol = p)
	l.rep <- matrix(nrow = rep, ncol = p)
	l.sum <- NULL
	
	for(r in 1:rep) {
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
			l.0 <- l.marginal.likelihood.lm(y, X[,-j]%*%Theta[-j,], psi, sigma.sq)
			l.j <- l.marginal.likelihood.lm(y, cbind(X[,-j]%*%Theta[-j,], X[,j]), psi, sigma.sq)
			
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





library(R.utils) # for intToBin()
marginal.likelihood <- function(y, X, lambda, psi, sigma.sq) {
	p = length(X[1,]) # Number of predictors
	n = length(y) # Sample size
	sum <- (1-lambda)^p # The marginal likelihood of the model without predictors is proportional to this
	
	for(k in 1:(2^p-1)) {
		gam <- rep(0, p)
		gam[min(p-ceiling(log2(k+1)-1),p):p] <- as.numeric(strsplit(intToBin(k), split = "")[[1]])
		X.gam <- X[,which(gam==1)]
		s.gam <- sum(gam)
		sum <- sum + lambda^s.gam*(1-lambda)^(p-s.gam)* marginal.likelihood.lm(y, X.gam, psi, sigma.sq)
	}
	
	return(sum)
}


post.mean <- function(y, X, lambda.star, psi, sigma.sq) {
	p = length(X[1,]) # Number of predictors
	n = length(y) # Sample size
	mean <- rep(0, p)
	
	for(k in 1:(2^p-1)) {
		gam <- rep(0, p)
		gam[min(p-ceiling(log2(k+1)-1),p):p] <- as.numeric(strsplit(intToBin(k), split = "")[[1]])
		X.gam <- X[,which(gam==1)]
		Phi <- t(X.gam)%*%X.gam + (sigma.sq/psi)*diag(sum(gam))
		mean[gam==1] <- mean[gam==1] + prod(lambda.star[gam==1])*prod(1-lambda.star[gam==0])*solve(Phi)%*%t(X.gam)%*%y
	}
	
	return(mean)
}


marginal.likelihood.lm <- function(y, X, psi, sigma.sq) {
	if(length(dim(X)) == 2) Phi <- t(X)%*%X + (sigma.sq/psi)*diag(length(X[1,]))
	else Phi <- t(X)%*%X + sigma.sq/psi
	return( exp( t(y)%*%X%*%solve(Phi)%*%t(X)%*%y / (2*sigma.sq) )  * det(( psi/sigma.sq )*Phi)^(-.5) )
}


exact <- function(y, A, lambda, psi, sigma.sq) {
	incProb <- NULL
	marginalLikelihood <- marginal.likelihood(y, A, lambda, psi, sigma.sq)
	
	for(j in 1:length(A[1,])) {
		marginalLikelihoodj <- marginal.likelihood(y, A[,-j], lambda, psi, sigma.sq)
		incProb[j] <- 1-(1-lambda)*marginalLikelihoodj/marginalLikelihood
	}
	
	return(list(mean=post.mean(y, A, incProb, psi, sigma.sq), incProb=incProb))
}


normalApprox.AMP <- function(y, X, lambda, psi, sigma.sq) {
	n <- length(y) # Sample size
	est.post.mean <- NULL
	est.post.inc.prob <- NULL

for(j in 1:length(X[1,])) {
	eigen.decom <- eigen(diag(n)-X[,j]%*%t(X[,j])/(t(X[,j])%*%X[,j])[1,1], symmetric = T)
	Q.j <- diag(sqrt(eigen.decom$values[-n]))%*%t(eigen.decom$vectors[,-n])
	AMP <- BGAMP(Q.j%*%y, Q.j%*%X[,-j], lambda, psi, sigma.sq)
	mu.star <- t(X[,j])%*%X[,-j]%*%AMP$mean
	sigma.sq.star <- t(X[,j])%*%( X[,-j]%*%diag(AMP$var)%*%t(X[,-j]) + sigma.sq*diag(length(y)))%*%X[,j]
	sigma.sq0 <- 1/( (t(X[,j])%*%X[,j])^2/sigma.sq.star + 1/psi )
	mean0 <- sigma.sq0*(t(X[,j])%*%y-mu.star)*(t(X[,j])%*%X[,j])/sigma.sq.star
	C0 <- dnorm(t(X[,j])%*%y, mean = mu.star, sd = sqrt((t(X[,j])%*%X[,j])^2*psi+sigma.sq.star))
	est.post.inc.prob[j] <- C0/( C0+dnorm(t(X[,j])%*%y, mean=mu.star, sd=sqrt(sigma.sq.star)) )
	est.post.mean[j] <- est.post.inc.prob[j]*mean0
}
	
	return(list(mean=est.post.mean, incProb=est.post.inc.prob))
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
	
	return(list(beta=beta, lambda=lambda, X=X, y=y, sigma.sq=sigma.sq))
}



MSEf <- function(rho) {
	if(exists("lambda")) rm(beta, lambda, X, y, sigma.sq, pos = 2L)
	attach(gendata(n = 100, p = 12, rho = rho, SNR = 2))
	
	psi = 10*sigma.sq
	
	ex <- exact(y, X, lambda, psi, sigma.sq)$incProb
	AMP <- normalApprox.AMP(y, X, lambda, psi, sigma.sq)$incProb
	BCRincProb <- BCR(y, X, lambda, psi, sigma.sq, m = 5)$incProb
	
	return( c( mean((ex-AMP)^2), mean((ex-BCRincProb)^2) ) )
}

n.rep <- 100 # Number of repeats per rho
g.rho <- 0:9/10 # Grid of rho values
n.rho <- length(g.rho)

MSEs <- array(dim = c(n.rho, n.rep, 2))

for(rho in g.rho) for(r in 1:n.rep) {
	print(paste(rho, r))
	MSEs[which(g.rho == rho), r, ] <- MSEf(rho)
}

save(MSEs, file = "MSEn100p12SNR2rep100")