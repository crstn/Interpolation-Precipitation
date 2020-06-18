# the sp library has all the interpolation stuff in R
library(sp)
library(sf)
library(dismo)
library(deldir)
library(gstat)
library(fields)

# set the working directory in R Studio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the stations
stations <- read.table("niederschlagsdaten/nieder_stationen.txt", 
                       header = TRUE, sep = ";", dec = ".", fileEncoding="latin1")
stations #check

# drop the columns we don't need:
stations <- stations[c('Name', 'east', 'north')]
stations #check

# load the daily precipitation data
precip <- read.table("Niederschlagsdaten/nieder_tageswerte.txt", 
                     header = TRUE, sep = ";", dec = ".", fileEncoding="latin1")
precip # you get the idea...

# we will only look at June 5, 2020 - that day has a nice distribution of different 
# values across NRW
datum = "2020-06-05"
precip <- precip[(precip$Datum == datum),]
precip 

# only keep the value and the station name
precip <- precip[c('Name', 'Tagessumme')] 
precip

# join on the station name to get the coordinates
precip <- merge(stations, precip, by="Name")
precip

# write that out, so we can look it in a GIS etc.
write.csv(precip, paste(datum, ".csv", sep="")) 

# Ready to go. Let's turn the data into a SpatialPointsDataFrame
preciPoints <- SpatialPoints(precip[,2:3], proj4string=CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) #UTM 32N
preciPoints <- SpatialPointsDataFrame(preciPoints, precip[,c(1,4)])

# check if we have duplicate points
zerodist(preciPoints)
# two stations at the same location - that will cause trouble later, let#s remove them both:
preciPoints <- preciPoints[-c(32,77),]

# check the maximum value
max(preciPoints$Tagessumme)

# read NRW shape file to plot the state boundary
NRW.sf <- st_read("NRW/NRW.shp")
NRW.sp <- as(NRW.sf, "Spatial")
NRW.pols <- list("sp.polygons", NRW.sp, fill = "lightgray")

spplot(preciPoints, 'Tagessumme',sp.layout=NRW.pols, 
       pch=20, cex=1, cuts=15, colorkey=TRUE)


# nearest neighbor interpolation

# just to demonstrate, make a Delaunay triangulation:
vtess=deldir(preciPoints@coords[,1],preciPoints@coords[,2])
plot(vtess, wlines='triang')
plot(vtess, wlines='both')


# create voronoi polygons
voro <- voronoi(preciPoints)
plot(voro)

# clip to NRW
voro.nrw <- intersect(voro, aggregate(NRW.sp))
spplot(voro.nrw, 'Tagessumme', main="Nearest neighbor interpolation")

# rasterize:

# blank raster
r <- raster(NRW.sf, res=2500)

vr <- rasterize(voro.nrw, r, 'Tagessumme')
spplot(vr, main="Nearest Neighbor")


# IDW

gs <- gstat(formula=Tagessumme~1, locations=preciPoints, nmax=7, set = list(idp = 1))
idw <- interpolate(r, gs)
## [inverse distance weighted interpolation]
idwr <- mask(idw, NRW.sf)
spplot(idwr)

# Spline
m <- Tps(coordinates(preciPoints), preciPoints$Tagessumme, lambda=0.000001)
tps <- interpolate(r, m)
tps <- mask(tps, NRW.sf)
spplot(tps, main="Spline mit lambda=0.000001")


# Semivariogram

gs <- gstat(formula=Tagessumme~1, locations=preciPoints)
v <- variogram(gs, width=10000, cutoff=200000)
plot(v, main="Semivariogramm", pch=20, col='red')

# fit exponential model
fve <- fit.variogram(v, vgm("Exp"))
fve
plot(variogramLine(fve, 200000), type='l', ylim=c(0,22), main="Semivariogramm (exponential)")
points(v[,2:3], pch=20, col='red')

# fit spherical model
fvsph <- fit.variogram(v, vgm("Sph"))
fvsph
plot(variogramLine(fvsph, 200000), type='l', ylim=c(0,22), main="Semivariogramm (spherical)")
points(v[,2:3], pch=20, col='red')

# using the exponential model for a kriging interpolation
k <- gstat(formula=Tagessumme~1, locations=preciPoints, model=fve)
# predicted values - first we need a grid:
g <- as(r, 'SpatialGrid')
kp <- predict(k, g)
spplot(kp)

# clip to NRW
ok <- brick(kp)
ok <- mask(ok, NRW.sf)
names(ok) <- c('prediction', 'variance')
spplot(ok)

# Cross-validation

# define root mean square error
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=Tagessumme~1, locations=train, nmax=nmx, set=list(idp=idp))
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$Tagessumme, p)
}

# randomly split into test (20%) and training data (80%)
set.seed(20200618)
i <- sample(nrow(preciPoints), 0.2 * nrow(preciPoints))
tst <- preciPoints[i,]
trn <- preciPoints[-i,]
opt <- optim(c(8, .5), f1, test=tst, train=trn)

nfolds <- 20
k <- kfold(preciPoints, nfolds)
spline.rmse <- kriging.rmse <- idw.rmse <- rep(NA, nfolds)

for (i in 1:nfolds) {
  test <- preciPoints[k!=i,]
  train <- preciPoints[k==i,]
  m <- gstat(formula=Tagessumme~1, 
             locations=train, 
             nmax=opt$par[1], 
             set=list(idp=opt$par[2]))
  p1 <- predict(m, newdata=test, debug.level=0)$var1.pred
  idw.rmse[i] <-  RMSE(test$Tagessumme, p1)
  m <- gstat(formula=Tagessumme~1, locations=train, model=fve)
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  kriging.rmse[i] <-  RMSE(test$Tagessumme, p2)
  m <- Tps(coordinates(train), train$Tagessumme)
  p3 <- predict(m, coordinates(test))
  spline.rmse[i] <-  RMSE(test$Tagessumme, p3)
  w <- c(idw.rmse[i], kriging.rmse[i], spline.rmse[i])
}

rmi <- mean(idw.rmse)
rmk <- mean(kriging.rmse)
rmt <- mean(spline.rmse)
rms <- c(rmi, rmt, rmk)
rms


boxplot(data.frame(idw.rmse, spline.rmse, kriging.rmse), 
        ylab="RMSE",
        labels=c("IDW", "Spline", "Kriging"))


# plot them all together
s <- stack(vr, idwr, tps, ok[[1]])
names(s) <- c('Nearest Neighbor', 'IDW', 'Spline', 'Kriging')
spplot(s)


