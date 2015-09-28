# ----------------------------------------------
# Script to extract 
# Christine Rollinson, crollinson@gmail.com
# Original: 28 September, 2015
#
# --------------
# Mask Proceedure (loop over each file)
# --------------
# 1) Open raw potential veg file
# 2) Crop to paleon domain extent
# 3) resample to make sure grid cells line up
# 4) Mask to make sure we have consistent land/water coverage
# 5) Convert to .nc file for consistent convention
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

paleon.mask <- "~/Dropbox/PalEON CR/env_regional/phase2_env_drivers_v1/domain_mask/paleon_domain.nc"
lulcc.path <- "~/Dropbox/PalEON CR/env_regional/env_drivers_raw/lulcc/LAND_USE_HARMONIZATION_V1/data/"
lulcc.out <- "~/Dropbox/PalEON CR/env_regional/phase2_env_drivers_v1/lulcc"

# Create the driver folder if it doesn't already exist
if(!dir.exists(lulcc.out)) dir.create(lulcc.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export mask
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)

# Set up file loop
lulcc.files <- dir(lulcc.path)
for(i in 1:length(lulcc.files)){
print(paste0("-- Loading File: ", lulcc.files[i]))	

# 1) load the raw data
# Assume that the long/lat is WGS84 and lines up with PalEON
nc.name <- strsplit(lulcc.files[i], "[.]")[[1]]
lu.code <- substr(nc.name[2],4, nchar(nc.name[2]))

# lulcc.raw <- stack(file.path(lulcc.path, lulcc.files), crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
lulcc.nc <- nc_open(file.path(lulcc.path, lulcc.files[i]))
summary(lulcc.nc$var)
var.name     <- lulcc.nc$var[[1]]$name
var.longname <- lulcc.nc$var[[1]]$longname
var.unit     <- lulcc.nc$var[[1]]$units
nc_close(lulcc.nc)

lulcc.raw <- stack(file.path(lulcc.path, lulcc.files[i]))
# lulcc.raw
# plot(lulcc.raw[[1:10]]) # Just to double check that things got read in okay

# 2) crop raw data (extents need to be same before masking)
lulcc.crop <- crop(lulcc.raw, extent(paleon))
# plot(lulcc.crop[[1]]) 
# plot(paleon, add=T, alpha=0.8)

# 3) mask raw data
lulcc.mask <- mask(lulcc.crop, paleon)
# plot(lulcc.mask[[1]])
# plot(paleon, add=T, alpha=0.3)
# plot(lulcc.mask[[1:6]])

# 6) save as netcdf
writeRaster(lulcc.mask, file.path(lulcc.out, paste0("paleon_lulcc_", lu.code, ".nc")), format="CDF", overwrite=T, varname=var.name, varunit=var.unit, longname=var.longname, zname="time", zunit="years since 1500 + 1")
} # End file loop
# ----------------------------------------------

# ----------------------------------------------
# Create some animations of landcover characters
# 1) Landcover Type
# 2) Wood Harvest
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

# List of codes of LU cover (not transitions)
lu.class <- c("gothr", "gsecd", "gcrop", "gpast")
harvest.f<- c("gfsh1", "gfsh2", "gfsh3", "gfvh1", "gfvh2")
harvest.b<- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2")


# -----------------------
# 1) Doing the dominant PFT 
# -----------------------
# get a spatial file for each PFT
gcrop <- stack(file.path(lulcc.out, dir(lulcc.out, "gcrop")))
gothr <- stack(file.path(lulcc.out, dir(lulcc.out, "gothr")))
gpast <- stack(file.path(lulcc.out, dir(lulcc.out, "gpast")))
gsecd <- stack(file.path(lulcc.out, dir(lulcc.out, "gsecd")))
# gurbn <- stack(file.path(lulcc.out, dir(lulcc.out, "gurbn"))) # Note: We're not using urban coverage

# Extract just the lat/lon
lat.lon <- data.frame(rasterToPoints(gcrop))[,1:2]
names(lat.lon)     <- c("lon", "lat")

# extract just the fraction of each landcover
# note: +2 is because it adds the lat/lon columns at the front
gcrop.df <- stack(data.frame(rasterToPoints(gcrop))[,3:(nlayers(gcrop)+2)])
gothr.df <- stack(data.frame(rasterToPoints(gothr))[,3:(nlayers(gothr)+2)])
gpast.df <- stack(data.frame(rasterToPoints(gpast))[,3:(nlayers(gpast)+2)])
gsecd.df <- stack(data.frame(rasterToPoints(gsecd))[,3:(nlayers(gsecd)+2)])

years <- as.numeric(substr(gcrop.df[,2],2,nchar(as.vector(gcrop.df[,2]))))-1+1500
land.use <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gothr=gothr.df[,1], gsecd=gsecd.df[,1], gcrop=gcrop.df[,1], gpast=gpast.df[,1])

# Find the dominant land use
land.use$LandUse <- as.factor(apply(land.use[, lu.class], 1, FUN=function(x){which(x==max(x))}))
levels(land.use$LandUse) <- c("Primary", "Secondary", "Cropland", "Pasture")
summary(land.use)

colors.lu <- c("darkgreen", "green3", "darkgoldenrod2", "darkorange2")

saveGIF(
for(i in min(land.use$year):max(land.use$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= land.use[land.use$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= LandUse)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(land.use$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(land.use$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_manual(values= colors.lu) +
    ggtitle(paste0("Dominant Land Use : ", i)) +
    guides(fill=guide_legend(direction="horizontal", nrow=2, title="lulcc") )+
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(lulcc.out,"Dominant_LandUse_1500-2005.gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="Dominant_LandUse_1500-2005.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------

# -----------------------
# 2) Doing total area harvest
# -----------------------
harvest.b<- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2")
# get a spatial file for each PFT
gsbh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh1")))
gsbh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh2")))
gsbh3 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh3")))
gvbh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gvbh1")))
gvbh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gvbh2"))) 

# Extract just the lat/lon
lat.lon <- data.frame(rasterToPoints(gsbh1))[,1:2]
names(lat.lon)     <- c("lon", "lat")

# extract just the fraction of each landcover
# note: +2 is because it adds the lat/lon columns at the front
df1 <- stack(data.frame(rasterToPoints(gsbh1))[,3:(nlayers(gsbh1)+2)])
df2 <- stack(data.frame(rasterToPoints(gsbh2))[,3:(nlayers(gsbh2)+2)])
df3 <- stack(data.frame(rasterToPoints(gsbh3))[,3:(nlayers(gsbh3)+2)])
df4 <- stack(data.frame(rasterToPoints(gvbh1))[,3:(nlayers(gvbh1)+2)])
df5 <- stack(data.frame(rasterToPoints(gvbh2))[,3:(nlayers(gvbh2)+2)])

years <- as.numeric(substr(df1[,2],2,nchar(as.vector(df1[,2]))))-1+1500
harvest.bm <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gsbh1=df1[,1], gsbh2=df2[,1], gsbh3=df3[,1], gvbh1=df4[,1], gvbh2=df5[,1])
harvest.bm$TotalHarvest <- rowSums(harvest.bm[,harvest.b])
summary(harvest.bm)

# Put Carbon on a Better scale
# assume 0.5 x 0.5 - degree roughly equal to 50 x 50 km
# 50 x 50 km = 5e4*5e4 m2
harvest.bm$TotalHarvest_MgHa <- harvest.bm$TotalHarvest*1e-3/((5e4*5e4)/10000)
summary(harvest.bm)

saveGIF(
for(i in min(harvest.bm$year):max(harvest.bm$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data=harvest.bm[harvest.bm$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= TotalHarvest_MgHa)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(harvest.bm$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(harvest.bm$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(name="Forest Harvest, Mg Ha-1", colours=c("gray50","darkgreen"), limits=c(0,quantile(harvest.bm$TotalHarvest_MgHa,0.95))) +
    ggtitle(paste0("Total Forest Harvest : ", i)) +
    # guides(fill=guide_legend(direction="horizontal", title="lulcc") )+
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(lulcc.out,"ForestHarvest_Biomass_1500-2005.gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="ForestHarvest_Biomass_1500-2005.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# dev.off()
# -----------------------

# -----------------------
# 3) Doing total harvested area
# -----------------------
harvest.f<- c("gfsh1", "gfsh2", "gfsh3", "gfvh1", "gfvh2")
# get a spatial file for each PFT
gfsh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gfsh1")))
gfsh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gfsh2")))
gfsh3 <- stack(file.path(lulcc.out, dir(lulcc.out, "gfsh3")))
gfvh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gfvh1")))
gfvh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gfvh2"))) 

# Extract just the lat/lon
lat.lon <- data.frame(rasterToPoints(gsbh1))[,1:2]
names(lat.lon)     <- c("lon", "lat")

# extract just the fraction of each landcover
# note: +2 is because it adds the lat/lon columns at the front
df1 <- stack(data.frame(rasterToPoints(gfsh1))[,3:(nlayers(gfsh1)+2)])
df2 <- stack(data.frame(rasterToPoints(gfsh2))[,3:(nlayers(gfsh2)+2)])
df3 <- stack(data.frame(rasterToPoints(gfsh3))[,3:(nlayers(gfsh3)+2)])
df4 <- stack(data.frame(rasterToPoints(gfvh1))[,3:(nlayers(gfvh1)+2)])
df5 <- stack(data.frame(rasterToPoints(gfvh2))[,3:(nlayers(gfvh2)+2)])

years <- as.numeric(substr(df1[,2],2,nchar(as.vector(df1[,2]))))-1+1500
harvest.area <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gfsh1=df1[,1], gfsh2=df2[,1], gfsh3=df3[,1], gfvh1=df4[,1], gfvh2=df5[,1])
harvest.area$TotalHarvest <- rowSums(harvest.area[,harvest.f])
summary(harvest.area)

# png(file.path(lulcc.out, "lulcc_PotentialVegetation_Dominant.png"), height=600, width=800)
saveGIF(
for(i in min(harvest.area$year):max(harvest.area$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= harvest.area[harvest.area$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= TotalHarvest)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(harvest.area$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(harvest.area$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(name="Forest Harvest, Fractional Area", colours=c("gray50","darkgreen"), limits=c(0,quantile(harvest.area$TotalHarvest,0.95))) +
    ggtitle(paste0("Total Forest Harvest : ", i)) +
    # guides(fill=guide_legend(direction="horizontal", title="lulcc") )+
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
)
# }, movie.name=file.path(lulcc.out,"ForestHarvest_Area_1500-2005.gif"), interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
# Note: Locally, movie.name doesn't like file.path, so use this setting
}, movie.name="ForestHarvest_Area_1500-2005.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------



# ----------------------------------------------
