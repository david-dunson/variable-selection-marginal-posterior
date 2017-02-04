w <- 2
h <- 3
margins <- c(4,4,0.2,0)


load("rho0SNR1")

pdf(file = "AMPrho0.pdf", width = w, height = h)
par(mar=margins)
boxplot(AMPincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()

pdf(file = "BCRrho0.pdf", width = w, height = h)
par(mar=margins)
boxplot(BCRincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()


load("rho2SNR1")

pdf(file = "AMPrho2.pdf", width = w, height = h)
par(mar=margins)
boxplot(AMPincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()

pdf(file = "BCRrho2.pdf", width = w, height = h)
par(mar=margins)
boxplot(BCRincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()


load("rho8SNR10")

pdf(file = "AMPrho8.pdf", width = w, height = h)
par(mar=margins)
boxplot(AMPincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()

pdf(file = "BCRrho8.pdf", width = w, height = h)
par(mar=margins)
boxplot(BCRincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()


load("rho5SNR10")

pdf(file = "AMPrho5.pdf", width = w, height = h)
par(mar=margins)
boxplot(AMPincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()

pdf(file = "BCRrho5.pdf", width = w, height = h)
par(mar=margins)
boxplot(BCRincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()


load("rho7SNR10")

pdf(file = "AMPrho7.pdf", width = w, height = h)
par(mar=margins)
boxplot(AMPincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()

pdf(file = "BCRrho7.pdf", width = w, height = h)
par(mar=margins)
boxplot(BCRincProb, outline = F, xlab = "parameter index", ylab = "Posterior inclusion probability")
dev.off()