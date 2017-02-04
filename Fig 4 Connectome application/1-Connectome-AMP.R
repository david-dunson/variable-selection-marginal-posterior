# AMP algorithm which takes MLE for sigma.sq based on current parameter estimate in each iteration, i.e. sample variance of estimated residuals
BGAMP <- function(y, A, psi, theta = 1) {
print(paste("theta:", theta))

# Set hyperparameters for the Bernoulli-Gaussian / spike-and-slab prior
dzeta = 0 # Mean of the normal


# Initialize lambda (prior inclusion probability)
lambda <- NULL
lambda[1] <- a/(a+b) # Prior mean


# Initialization for the AMP algorithm
n.iter = 742 # Maximum number of iterations for AMP algorithm
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
sigma.sq <- NULL

convergence <- NULL

# The AMP algorithm
for(k in 1:(n.iter-1)) {
	counter <- counter+1
	if((counter %% 10) == 0) print(paste("Counter:", counter))
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
	if(k <= 713) sigma.sq[k] <- (2*b.0+sum((y-z[k,])^2))/(2*a.0+m)
	else sigma.sq[k] <- sigma.sq[k-1]
	#print(paste("sigma.sq:", sigma.sq))
	
	# Step 2
	for(aa in 1:m) {
		tau.z[k,aa] <- tau.p[k,aa]*sigma.sq[k]/(tau.p[k,aa]+sigma.sq[k]) # Var(z_a|y,p.hat[k,a],tau.p[k,a])
		z.hat <- tau.z[k,aa]*(y[aa]/sigma.sq[k]+p.hat[k,aa]/tau.p[k,aa]) # E(z.a|y,p.hat[k,a],tau.p[k,a])
		
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
	if(k <= 713) lambda[k+1] <- (a+sum(pi[k,]))/(a+b+n)
	else lambda[k+1] <- lambda[k]
	#print(paste("lambda:", lambda[k+1]))
	
	# for debugging
	#print(x.hat[k+1,1])
	
	th <- 1E-5
	convergence[counter] <- sd(x.hat[k,]-x.hat[k+1,])
	
	print(convergence[counter])
	
	if(is.na(convergence[counter] < th)) {
		print("Diverged!")
		break
	}
	#else if(sum(abs(x.hat[k,]-x.hat[k+1,])) < th) break # Stop iterating if AMP has sufficiently converged
}

return(list(mean=x.hat[counter,], var=tau.x[counter,], incProb=pi[counter,], lambda=lambda, sigma.sq=sigma.sq, convergence=convergence))

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



AMPresult <- BGAMP(y, con, 1, .0001)
plot(log(AMPresult$convergence))

save(AMPresult, file = 'AMPresult')