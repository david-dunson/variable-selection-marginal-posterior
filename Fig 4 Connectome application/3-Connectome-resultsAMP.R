load("connectome_data")
n <- length(y)
subset <- NULL
i=1
for(j in 1:70^2) {
	if(!all(connect[,j] == rep(0,n))) {
		subset[i] <- j
		i <- i+1
	}
}

load("AMPresult")

BCRincProb <- AMPresult$incProb # For convenience, we use BCRincProb as variable with results both when we use BCR as well as AMP
BCR.a <- rep(0, 70^2)
BCR.a[subset] <- BCRincProb
BCR.A <- matrix(BCR.a,70,70)
BCR.A <- BCR.A+t(BCR.A)



library(qgraph)


V<-70
Coord_Brain<-matrix(,V,2)

Coord_Brain[1,]<-c(-3,9.5)
Coord_Brain[2,]<-c(-7.1,8.9)
Coord_Brain[3,]<-c(-0.8,16.4)
Coord_Brain[4,]<-c(-4.7,14.2)
Coord_Brain[5,]<-c(-1.3,12.1)
Coord_Brain[6,]<-c(-1.1,1.9)
Coord_Brain[7,]<-c(-3.9,15.3)
Coord_Brain[8,]<-c(-4.3,8)
Coord_Brain[9,]<-c(-5,4.2)
Coord_Brain[10,]<-c(-6.4,10.8)
Coord_Brain[11,]<-c(-1.1,7.6)
Coord_Brain[12,]<-c(-3.3,1.5)
Coord_Brain[13,]<-c(-3,19.7)
Coord_Brain[14,]<-c(-1.7,5.5)
Coord_Brain[15,]<-c(-0.9,20.2)
Coord_Brain[16,]<-c(-7.6,11.3)
Coord_Brain[17,]<-c(-3.8,11.4)
Coord_Brain[18,]<-c(-1.5,13.5)
Coord_Brain[19,]<-c(-6.3,16)
Coord_Brain[20,]<-c(-5.2,21)
Coord_Brain[21,]<-c(-6.2,18.3)
Coord_Brain[22,]<-c(-1.2,3.3)
Coord_Brain[23,]<-c(-5.7,9.4)
Coord_Brain[24,]<-c(-0.7,10.8)
Coord_Brain[25,]<-c(-5.3,12.1)
Coord_Brain[26,]<-c(-1.3,4.8)
Coord_Brain[27,]<-c(-0.9,19.1)
Coord_Brain[28,]<-c(-4.5,19.5)
Coord_Brain[29,]<-c(-1.5,16.8)
Coord_Brain[30,]<-c(-2.6,4.2)
Coord_Brain[31,]<-c(-6.8,12.7)
Coord_Brain[32,]<-c(-6.8,8.2)
Coord_Brain[33,]<-c(-0.7,23)
Coord_Brain[34,]<-c(-3.5,17.3)
Coord_Brain[35,]<-c(-5.4,10.6)
Coord_Brain[36:70,]<-Coord_Brain[1:35,]
Coord_Brain[36:70,1]<--Coord_Brain[36:70,1]

LAB<-c("0L","1L","2L","3L","4L","5L","6L","7L","8L","9L","10L","11L","12L","13L","14L","15L","16L","17L","18L","19L","20L","21L","22L","23L","24L","25L","26L","27L","28L","29L","30L","31L","32L","33L","34L","0R","1R","2R","3R","4R","5R","6R","7R","8R","9R","10R","11R","12R","13R","14R","15R","16R","17R","18R","19R","20R","21R","22R","23R","24R","25R","26R","27R","28R","29R","30R","31R","32R","33R","34R")
SHA<-c("square","triangle","square","circle","square","rectangle","triangle","triangle","diamond","triangle","square","rectangle","circle","rectangle","circle","triangle","square","circle","circle","circle","circle","rectangle","diamond","square","circle","diamond","square","circle","circle","diamond","triangle","diamond","circle","triangle","triangle")
SHA<-c(SHA,SHA)


SHA_NEW<-SHA
SHA_NEW[SHA_NEW=="triangle"]<-"square"
SHA_NEW[SHA_NEW=="diamond"]<-"square"
SHA_NEW[SHA_NEW=="rectangle"]<-"square"
SHA_NEW[SHA_NEW=="circle"]<-"circle"



diag(BCR.A)<-NA
BCR.A0 <- BCR.A
BCR.A0[is.na(BCR.A0)] <- 0
l.b <- sort(as.vector(BCR.A0), decreasing = T)[60]
BCR.A <- BCR.A0 >= l.b


val2col<-function(z, zlim, col = heat.colors(100), breaks){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 colorlevels <- col[((as.vector(z)-breaks[1])/(range(breaks)[2]-range(breaks)[1]))*(length(breaks)-1)+1] # assign colors to heights for each point
 colorlevels
}

col.vec <- val2col(as.vector(BCR.A0)[as.vector(BCR.A)])
col.vec.2 <- rep(0,4900)
col.vec.2[as.vector(BCR.A)] <- col.vec
colors <- matrix(col.vec.2, 70, 70)



# Create plot

w <- 10
h <- 10
margins <- c(0,0,0,0)

pdf(file = "Plots/AMPnew.pdf", width = w, height = h)
par(mar=margins)
qgraph(BCR.A, bg= "gray",layout=Coord_Brain,labels=LAB,shape=SHA_NEW, edge.color = colors, color = "white", pin = c(3,.5), width = w, height = h)
dev.off()
