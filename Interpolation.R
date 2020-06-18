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
r <- raster(NRW, res=2500)

vr <- rasterize(voro.nrw, r, 'Tagessumme')
spplot(vr, main="Nearest Neighbor")


# IDW

gs <- gstat(formula=Tagessumme~1, locations=preciPoints, nmax=7, set = list(idp = 1))
idw <- interpolate(r, gs)
## [inverse distance weighted interpolation]
idwr <- mask(idw, NRW)
spplot(idwr)

# Spline
m <- Tps(coordinates(preciPoints), preciPoints$Tagessumme, lambda=0.000001)
tps <- interpolate(r, m)
tps <- mask(tps, NRW)
spplot(tps, main="Spline mit lambda=0.000001")




# define root mean square error
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}





