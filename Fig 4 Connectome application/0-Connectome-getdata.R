data <- read.csv("Covariates_Responces.csv", header = T, sep="\t")
n <- length(data$URSI)

library(R.matlab) # For readMat()

connect <- matrix(nrow = n, ncol = 70^2)
for(i in 1:n) connect[i,] <- as.vector(readMat(paste(
	"http://openconnecto.me/data/public/MR/MIGRAINE_v1_0/MRN_111/small_graphs/MRN-114_",
	data$URSI[i], "_small_graph.mat", sep = ""
))$fibergraph)

y <- as.numeric(gsub(",", ".", levels(data$CCI)))[data$CCI]

# Delete 105th observation as it has CCI missing
y <- y[-105]
connect <- connect[-105,]

save(y, connect, file="connectome_data")