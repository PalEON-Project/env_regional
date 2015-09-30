# ----------------------------------------------
# Script to extract 
# Christine Rollinson, crollinson@gmail.com
# Original: 29 September, 2015
#
# --------------
# Mask Proceedure (loop over each file)
# --------------
# 1) Open raw nitrogen files
# 2) Crop to paleon domain extent
# 3) Mask to make sure we have consistent land/water coverage
# 4) Convert to .nc file for consistent convention
#
# NOTE: We're going to maintain the original file structure for
#       what we're distributing, but I'll also write a set of 
#       utils that will organize into separate time series by
#       grid cell ID (i.e. ED) or by year (MsTMIP)
# --------------
#
# ----------------------------------------------


# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4); library(raster); library(rgdal)

paleon.mask <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc"
nitrogen.path <- "~/Dropbox/PalEON_CR/env_regional/env_drivers_raw/nitrogen/NACP_MSTMIP_MODEL_DRIVER/data/"
nitrogen.out <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/nitrogen"

# Create the driver folder if it doesn't already exist
if(!dir.exists(nitrogen.out)) dir.create(nitrogen.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export mask
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)

# The "native" nitrogen driver goes 1860-2050
nitrogen.years <- 1860:2050

# 1) load the raw data
nhx.nc <- nc_open(file.path(nitrogen.path, "mstmip_driver_global_hd_nitrogen_nhx_v1.nc4"))
noy.nc <- nc_open(file.path(nitrogen.path, "mstmip_driver_global_hd_nitrogen_noy_v1.nc4"))

nhx.name     <- nhx.nc$var[[5]]$name
nhx.longname <- nhx.nc$var[[5]]$longname
nhx.unit     <- nhx.nc$var[[5]]$units
# nhx.lat      <- ncvar_get(nhx.nc, "lat")
# nhx.lon      <- ncvar_get(nhx.nc, "lon")
# nhx.dat      <- ncvar_get(nhx.nc, nhx.name)[,,which(nitrogen.years<=2010)]
dim(nhx.dat) # dim=c(lon, lat, time)

noy.name     <- noy.nc$var[[5]]$name
noy.longname <- noy.nc$var[[5]]$longname
noy.unit     <- noy.nc$var[[5]]$units
# noy.lat      <- ncvar_get(noy.nc, "lat")
# noy.lon      <- ncvar_get(noy.nc, "lon")
# noy.dat      <- ncvar_get(noy.nc, nhx.name)[,,which(nitrogen.years<=2010)]
dim(noy.dat) # dim=c(lon, lat, time)

nc_close(nhx.nc); nc_close(noy.nc)

# Taking the spatial data, but ignoring the years beyond the paleon temporal domain
# Note: there will be some weird message about an invalid CRS, but it's fine
nhx.raw <- stack(file.path(nitrogen.path, "mstmip_driver_global_hd_nitrogen_nhx_v1.nc4"))[[which(nitrogen.years<=2010)]]
noy.raw <- stack(file.path(nitrogen.path, "mstmip_driver_global_hd_nitrogen_noy_v1.nc4"))[[which(nitrogen.years<=2010)]]
# nhx.raw
# noy.raw
# plot(nhx.raw[[1:6]]) # Just to double check that things got read in okay
# plot(noy.raw[[1:6]]) # Just to double check that things got read in okay

# 2) crop raw data (extents need to be same before masking)
nhx.crop <- crop(nhx.raw, extent(paleon))
noy.crop <- crop(noy.raw, extent(paleon))
# plot(nhx.crop[[1]]) 
# plot(paleon, add=T, alpha=0.8)

# plot(noy.crop[[1]]) 
# plot(paleon, add=T, alpha=0.8)

# 3) mask raw data
nhx.mask <- mask(nhx.crop, paleon)
noy.mask <- mask(noy.crop, paleon)
# plot(nhx.mask[[1]])
# plot(paleon, add=T, alpha=0.4)
# plot(nhx.mask[[1:6]])

# plot(noy.mask[[1]])
# plot(paleon, add=T, alpha=0.3)
# plot(noy.mask[[1:6]])

# 6) save as netcdf
writeRaster(nhx.mask, file.path(nitrogen.out, "paleon_nhx.nc"), format="CDF", overwrite=T, varname=nhx.name, varunit=nhx.unit, longname=nhx.longname, zname="time", zunit="years since 1860")
writeRaster(noy.mask, file.path(nitrogen.out, "paleon_noy.nc"), format="CDF", overwrite=T, varname=noy.name, varunit=noy.unit, longname=noy.longname, zname="time", zunit="years since 1860")
# ----------------------------------------------

# ----------------------------------------------
# Create some animations of nitrogen deposition
# 1) Total N Deposition
# 2) NHx
# 3) NOy
# ----------------------------------------------
# Load Spatial & Graphing Libraries
library(raster); library(ggplot2); library(grid); 
library(animation)
library(car)

# A couple settings for running on my local mac only; comment these out elsewhere
path.to.convert <- "/opt/local/bin/convert"
ani.options(convert=path.to.convert)

# Get state outlines
paleon.states <- map_data("state")

# get a spatial file for each PFT
nhx <- stack(file.path(nitrogen.out, "paleon_nhx.nc"))
noy <- stack(file.path(nitrogen.out, "paleon_noy.nc"))

plot(nhx[[1:6]])
plot(noy[[1:6]])

# Convert each layer to a dataframe
# Extract just the lat/lon
lat.lon <- data.frame(rasterToPoints(gcrop))[,1:2]
names(lat.lon)     <- c("lon", "lat")

# Extract the N deposition data & code in years
nhx.df <- stack(data.frame(rasterToPoints(nhx))[,3:(nlayers(nhx)+2)])
noy.df <- stack(data.frame(rasterToPoints(noy))[,3:(nlayers(noy)+2)])
years <- as.numeric(substr(nhx.df[,2],2,nchar(as.vector(nhx.df[,2]))))-1+1860

# Making a single Nitrogen deposition data frame
nitrogen <- data.frame(lon=lat.lon$lon, lat=lat.lon$lat, year=years, nhx=nhx.df[,1], noy=noy.df[,1])
nitrogen$Total <- nitrogen$nhx + nitrogen$noy
summary(nitrogen)

# -----------------------
# 1) Total N Deposition
# -----------------------
saveGIF(
for(i in min(nitrogen$year):max(nitrogen$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= nitrogen[nitrogen$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill=Total)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(nitrogen$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(nitrogen$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "blue2"), name="MgN/m2/yr", limits=range(nitrogen$Total)) +
    ggtitle(paste0("Total N Deposition : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(nitrogen.out,"NitrogenDeposition_Total_1860-2010.gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="NitrogenDeposition_Total_1860-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------

# -----------------------
# 2) NHx Deposition
# -----------------------
saveGIF(
for(i in min(nitrogen$year):max(nitrogen$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= nitrogen[nitrogen$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill=nhx)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(nitrogen$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(nitrogen$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "blue2"), name="MgN/m2/yr", limits=range(nitrogen$nhx)) +
    ggtitle(paste0("NHx Deposition : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(nitrogen.out,"NitrogenDeposition_NHx_1860-2010.gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="NitrogenDeposition_NHx_1860-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------

# -----------------------
# 3) NOy Deposition
# -----------------------
saveGIF(
for(i in min(nitrogen$year):max(nitrogen$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= nitrogen[nitrogen$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill=noy)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(nitrogen$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(nitrogen$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "blue2"), name="MgN/m2/yr", limits=range(nitrogen$noy)) +
    ggtitle(paste0("NHx Deposition : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(nitrogen.out,"NitrogenDeposition_NOy_1860-2010..gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="NitrogenDeposition_NOy_1860-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------
# ----------------------------------------------

