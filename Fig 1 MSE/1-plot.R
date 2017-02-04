load("MSEn100p12SNR2rep100")

n.rep <- 100
n.rho <- 10

AMP1 <- NULL
AMP2 <- NULL
AMP3 <- NULL

BCR1 <- NULL
BCR2 <- NULL
BCR3 <- NULL


for(r in 1:n.rho) {
	AMP <- sort(MSEs[r,,1])
	AMP1[r] <- quantile(AMP, .2)
	AMP2[r] <- mean(AMP)
	AMP3[r] <- quantile(AMP, .8)
	
	BCR <- sort(MSEs[r,,2])
	BCR1[r] <- quantile(BCR, .2)
	BCR2[r] <- mean(BCR)
	BCR3[r] <- quantile(BCR, .8)
}


# Create plot

w <- 4
h <- 2
margins <- c(4,4,0,0)


rho <- 0:9/10
rho1 <- rho-.01
rho2 <- rho+.01
pdf(file = "MSE.pdf", width = w, height = h)
par(mar=margins)
plot(rho1, AMP3, col = 'red', ylim = c(0,.1), type = 'n', ylab = 'MSE', xlab = 'rho')
segments(x0 = rho1, y0 = AMP1, y1 = AMP3, col = 'red')
points(rho1, AMP2, col = 'red', pch = 20)

segments(x0 = rho2, y0 = BCR1, y1 = BCR3, col = 'blue')
points(rho2, BCR2, col = 'blue', pch = 20)
dev.off()