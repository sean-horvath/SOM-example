library(kohonen)
library(raster)
library(rasterVis)
library(maptools)
library(RColorBrewer)


# Arctic Melt Onset SOM ---------------------------------------------------

# EASE2.0 Coordinate Reference System
ease <- '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
# read in data
gph <- read.csv('data/gph_data.csv',header=T)

# convert data frame into spatial object, then raster stack
coordinates(gph) <- ~lon+lat
gridded(gph)=TRUE
gphstack <- stack(gph)
data(wrld_simpl)
crs(gphstack) <- crs(wrld_simpl)

# re-project raster to EASE crs (twice to fix some missing points)
gphstack <- projectRaster(gphstack,
                          crs=CRS(ease))
gphstack <- projectRaster(gphstack,
                          crs=CRS(ease))

# find locations with NA values and remove
all_vals <- getValues(gphstack)
index <- complete.cases(all_vals)
vals <- all_vals[index,]

# Perform SOM (timing it just for fun)
begin <- Sys.time()
SOM <- som(scale(t(vals)),
           grid = somgrid(5,5,"rectangular"))
end <- Sys.time()
end-begin

png(filename='images/som_plot.png')
plot(SOM)
dev.off()

# Plot number of observations in each node
colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

png(filename='images/counts_plot.png')
plot(SOM, type = "counts",
     palette.name = colors,
     heatkey = TRUE)
dev.off()

# Display composited GPH data into nodes
all_vals <- t(all_vals)
rlist <- list()
for(i in 1:25){
  tmp_df <- all_vals[which(SOM$unit.classif==i),]
  vals <- colMeans(tmp_df,na.rm=T)
  rlist[[i]] <- setValues(raster(gphstack,layer=1),vals)
}

rstack <- stack(rlist)

# Create spatial polygons object to map land masses
out <- crop(wrld_simpl, extent(-180, 180, 60, 83.57027))
wrld_arctic <- spTransform(out,
                        CRS(ease))

# Crop for nicer maps
e <- extent(-3700000,3700000,-3700000,3700000)
rstack <- crop(rstack,e)

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'RdBu')),
                        axis.line=list(col='transparent'))

# Plot GPH data onto SOM map
png(filename='images/GPH Master SOM.png')
levelplot(rstack,margin=F,
          par.settings=mapTheme,
          main='GPH @ 500hpa Master SOM Map',
          at=seq(4850,5750,by=100),
          layout=c(5,5),scales=list(draw=FALSE),
          names.attr=c(1,2,3,4,5,rep('',20)),
          between=list(x=0,y=0)) +
  layer(sp.lines(wrld_arctic, col="grey20",lwd=0.2),under=F)
dev.off()


# Map surface air temperature to Master SOM -------------------------------

temps <- read.csv('data/temperature_data.csv',header=T)

# convert data frame into spatial object, then raster stack
coordinates(temps) <- ~lon+lat
gridded(temps)=TRUE
tempsstack <- stack(temps)
crs(tempsstack) <- crs(wrld_simpl)
# re-project raster to EASE crs
tempsstack <- projectRaster(tempsstack,
                          crs=CRS(ease))
tempsstack <- projectRaster(tempsstack,
                          crs=CRS(ease))

# get values from raster stack and transpose so columns are locations
tempvals <- getValues(tempsstack)
tempvals <- t(tempvals)

# composite temperature days for each node
rlist <- list()
for(i in 1:25){
  tmp_df <- tempvals[which(SOM$unit.classif==i),]
  vals <- colMeans(tmp_df,na.rm=T)
  rlist[[i]] <- setValues(raster(tempsstack,layer=1),vals)
}

rstack <- stack(rlist)
rstack <- crop(rstack,e)

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'RdBu')),
                        axis.line=list(col='transparent'))

png(filename='images/Temperature Mapped to SOM.png')
levelplot(rstack,margin=F,
          par.settings=mapTheme,
          main='Temperature @ 925hpa SOM Map',
          at=seq(240,295,by=5),
          layout=c(5,5),scales=list(draw=FALSE),
          names.attr=c(1,2,3,4,5,rep('',20)),
          between=list(x=0,y=0)) +
  layer(sp.lines(wrld_arctic, col="grey20",lwd=0.2),under=F)
dev.off()



