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
# Corrections & gap filling
# 6) Make sure land cover types sum to 1 -- if it doesn't redistribute remaining fraction
# 7) Fill land cover type 850-1500 as all primary forest
# 8) Set all 850-1500 disturbance rates to 0
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

paleon.mask <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc"
lulcc.path <- "~/Desktop/Research/PalEON_CR/env_regional/env_drivers_raw/lulcc/LAND_USE_HARMONIZATION_V1/data/"
lulcc.out <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/lulcc"

# Create the driver folder if it doesn't already exist
if(!dir.exists(lulcc.out)) dir.create(lulcc.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export raw masked data
# 1) Open raw potential veg file
# 2) Crop to paleon domain extent
# 3) resample to make sure grid cells line up
# 4) Mask to make sure we have consistent land/water coverage
# 5) Convert to .nc file for consistent convention
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

# 4) Conservatively remove some outliers
# Outlier criteria: >=2 times the 99.9 percentile
# If outlier is >= 2 x 0.999%˚ for a year, make it equal the 99.9%˚ for that year

for(i in 1:nlayers(lulcc.mask)){
	layer.max  <- quantile(values(lulcc.mask[[i]]),0.999, na.rm=T)
	fill.cells <- which(values(lulcc.mask[[i]])>2*layer.max)
	print(paste0("     Layer ", i, ";  n = ", length(fill.cells)))
	if(length(fill.cells)>0){
    	lulcc.mask[[i]][fill.cells] <- 2*layer.max
    }
}

# 5) save as netcdf
writeRaster(lulcc.mask, file.path(lulcc.out, paste0("paleon_lulcc_", lu.code, ".nc")), format="CDF", overwrite=T, varname=var.name, varunit=var.unit, longname=var.longname, zname="time", zunit="years since 1500 + 1")
} # End file loop
# ----------------------------------------------

# ----------------------------------------------
# Doing some gapfilling & error removing
# 6) Make sure land cover types sum to 1 -- if it doesn't redistribute remaining fraction
# 7) Fill land cover type 850-1500 as all primary forest
# 8) Set all 850-1500 disturbance rates to 0
# ----------------------------------------------

# Load the mask
paleon <- raster(paleon.mask)

lu.class    <- c("gothr", "gsecd", "gcrop", "gpast") # Fraction land in each class
harvest.f   <- c("gfsh1", "gfsh2", "gfsh3", "gfvh1", "gfvh2") # Fraction area harvested
harvest.b   <- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2") # Biomass harvested
second.char <- c("gssma", "gssmb") # mean age & biomass of secondary land
transitions <- c("gflcp", "gflcs", "gflpc", "gflps", "gflsc", "gflsp", "gflvc", "gflvp") # Transitions among LU types

# ---------------------------------
# Re-distributing land cover types so that fractions sum to 1
# ---------------------------------
# Open land cover types
gcrop <- stack(file.path(lulcc.out, dir(lulcc.out, "gcrop")))
gothr <- stack(file.path(lulcc.out, dir(lulcc.out, "gothr")))
gpast <- stack(file.path(lulcc.out, dir(lulcc.out, "gpast")))
gsecd <- stack(file.path(lulcc.out, dir(lulcc.out, "gsecd")))


# Creating a new data frame that we fill with the assumption of no land cover change 
#   assume gothr (primary land cover type) = 1 (no human influence before 1500)
#   gcrop, gpast, gsecd = 0
gcrop2 <- gothr2 <- gpast2 <- gsecd2 <- brick(ncols=ncol(gcrop), nrows=nrow(gcrop), nl=length(850:2010),
                                              xmn=xmin(gcrop), xmx=xmax(gcrop), 
                                              ymn=ymin(gcrop), ymx=ymax(gcrop),
                                              crs=projection(gcrop))
names(gcrop2) <- names(gothr2) <- names(gpast2) <- names(gsecd2) <- paste0("X", 850:2010)
gcrop2 <- setValues(gcrop2, 0)
gpast2 <- setValues(gpast2, 0)
gsecd2 <- setValues(gsecd2, 0)
gothr2 <- setValues(gothr2, 1)

# Sum land cover types
gsum <- gcrop + gothr + gpast + gsecd
gsum
# plot(gsum[[1]])
# plot(gothr[[1]])

# Adding in the new layer normalized by the sum of the land cover types
# if land cover fraction sums to 1, nothing should change, otherwise it'll be normalized
for(i in 1:nlayers(gsum)){
	print(paste0("    Filling Year : ", i-1+1500))
	gcrop2[[650+i]] <- gcrop[[i]]/gsum[[i]]
	gothr2[[650+i]] <- gothr[[i]]/gsum[[i]]
	gpast2[[650+i]] <- gpast[[i]]/gsum[[i]]
	gsecd2[[650+i]] <- gsecd[[i]]/gsum[[i]]
}

# Assume no change in land use after 2005 
# (not great, but better than trying to make sure harvest rates line up etc)
for(i in (nlayers(gcrop2)-4):nlayers(gcrop2)){
	print(paste0("    Filling Year : ", i-1+850))
	gcrop2[[i]] <- gcrop[[nlayers(gcrop)]]/gsum[[nlayers(gsum)]]
	gothr2[[i]] <- gothr[[nlayers(gothr)]]/gsum[[nlayers(gsum)]]
	gpast2[[i]] <- gpast[[nlayers(gpast)]]/gsum[[nlayers(gsum)]]
	gsecd2[[i]] <- gsecd[[nlayers(gsecd)]]/gsum[[nlayers(gsum)]]
}

# redoing the mask to make sure everything is okay
gcrop2 <- mask(gcrop2, paleon)
gothr2 <- mask(gothr2, paleon)
gpast2 <- mask(gpast2, paleon)
gsecd2 <- mask(gsecd2, paleon)

# plot(gothr2[[1157]])
# plot(gsecd2[[1157]])
# plot(gcrop2[[1157]])
# plot(gpast2[[1157]])

writeRaster(gcrop2, file.path(lulcc.out, paste0("paleon_lulcc_gcrop.nc")), format="CDF", overwrite=T, varname="prop_crop", varunit="fraction area", longname="Proportion of landcover in crops", zname="time", zunit="years since 850")
writeRaster(gpast2, file.path(lulcc.out, paste0("paleon_lulcc_gpast.nc")), format="CDF", overwrite=T, varname="prop_pasture", varunit="fraction area", longname="Proportion of landcover in pasture", zname="time", zunit="years since 850")
writeRaster(gothr2, file.path(lulcc.out, paste0("paleon_lulcc_gothr.nc")), format="CDF", overwrite=T, varname="prop_primary", varunit="fraction area", longname="Proportion of landcover in primary landcover", zname="time", zunit="years since 850")
writeRaster(gsecd2, file.path(lulcc.out, paste0("paleon_lulcc_gsecd.nc")), format="CDF", overwrite=T, varname="prop_secd", varunit="fraction area", longname="Proportion of landcover in secondary landcover", zname="time", zunit="years since 850")
# ---------------------------------


# ---------------------------------
# Set all other transitions, harvest, etc to 0 for the time outside of the 
# actual landcover driver
# ---------------------------------
harvest.f   <- c("gfsh1", "gfsh2", "gfsh3", "gfvh1", "gfvh2") # Fraction area harvested
harvest.b   <- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2") # Biomass harvested
second.char <- c("gssma", "gssmb") # mean age & biomass of secondary land
transitions <- c("gflcp", "gflcs", "gflpc", "gflps", "gflsc", "gflsp", "gflvc", "gflvp") # Transitions among LU types

lulcc.files <- c(harvest.f, harvest.b, second.char, transitions)
for(i in 1:length(lulcc.files)){
	print(paste0("-- Loading File: ", lulcc.files[i]))	
	# 1) load the original masked
	# Assume that the long/lat is WGS84 and lines up with PalEON
	lu.code <- lulcc.files[i]

	# Extract the variable information
	lulcc.nc <- nc_open(file.path(lulcc.out, paste0("paleon_lulcc_", lulcc.files[i], ".nc")))
	summary(lulcc.nc$var)
	var.name     <- lulcc.nc$var[[1]]$name
	var.longname <- lulcc.nc$var[[1]]$longname
	var.unit     <- lulcc.nc$var[[1]]$units
	nc_close(lulcc.nc)

	# Load the original data
	lulcc.raw <- stack(file.path(lulcc.out, paste0("paleon_lulcc_", lulcc.files[i], ".nc")))

	# Set up a blank brink with the right dimensions
	lulcc2 <- brick(ncols=ncol(lulcc.raw), nrows=nrow(lulcc.raw), nl=length(850:2010),
                                              xmn=xmin(lulcc.raw), xmx=xmax(lulcc.raw), 
                                              ymn=ymin(lulcc.raw), ymx=ymax(lulcc.raw),
                                              crs=projection(lulcc.raw))
	names(lulcc2) <- paste0("X", 850:2010)
	lulcc2 <- setValues(lulcc2, 0)
	for(i in 1:nlayers(lulcc.raw)){
		lulcc2[[650+i]] <- lulcc.raw[[i]]
	}

	# if we're dealing with the secondary age file, we just need to keep adding 1 
	# to the age
	# (not great, but better than trying to make sure harvest rates line up etc)
	age <- 1
	if(lu.code == "gssma"){
		for(i in (nlayers(lulcc2)-4):nlayers(gcrop2)){
			lulcc2[[i]] <- lulcc2[[nlayers(gcrop)]]+age
			age <- age+1
		}
	}

	# redoing the mask to make sure everything is okay
	lulcc2 <- mask(lulcc2, paleon)

	writeRaster(lulcc2, file.path(lulcc.out, paste0("paleon_lulcc_", lu.code, ".nc")), format="CDF", overwrite=T, varname=var.name, varunit=var.unit, longname=var.longname, zname="time", zunit="years since 850 + 1")
                                              
}
# ---------------------------------
# ----------------------------------------------








# ----------------------------------------------
# ----------------------------------------------
# Create some animations of landcover characters
# 1) Landcover Type
# 2) Wood Harvest
# ----------------------------------------------
# ----------------------------------------------
# Load Spatial & Graphing Libraries
library(raster); library(ggplot2); library(grid); 
library(animation)
library(car)

# A couple settings for running on my local mac only; comment these out elsewhere
path.to.convert <- "/opt/ImageMagick/bin/convert"
ani.options(convert=path.to.convert)

# Get state outlines
paleon.states <- map_data("state")

# List of codes of LU cover (not transitions)
lu.class    <- c("gothr", "gsecd", "gcrop", "gpast") # Fraction land in each class
harvest.f   <- c("gfsh1", "gfsh2", "gfsh3", "gfvh1", "gfvh2") # Fraction area harvested
harvest.b   <- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2") # Biomass harvested
second.char <- c("gssma", "gssmb") # mean age & biomass of secondary land
transitions <- c("gflcp", "gflcs", "gflpc", "gflps", "gflsc", "gflsp", "gflvc", "gflvp") # Transitions among LU types
# -----------------------
# 1) Doing the dominant  & Fractional Land Uses 
# -----------------------
lu.class    <- c("gothr", "gsecd", "gcrop", "gpast") # Fraction land in each class

# get a spatial file for each land use type
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

years <- as.numeric(substr(gcrop.df[,2],2,nchar(as.vector(gcrop.df[,2]))))-1+850
land.use <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gothr=gothr.df[,1], gsecd=gsecd.df[,1], gcrop=gcrop.df[,1], gpast=gpast.df[,1])

# Find the dominant land use
land.use$LandUse <- as.factor(apply(land.use[, lu.class], 1, FUN=function(x){which(x==max(x))}))
levels(land.use$LandUse) <- c("Primary", "Secondary", "Cropland", "Pasture")
land.use$LU.sum <- rowSums(land.use[,c("gothr", "gsecd", "gcrop", "gpast")])
summary(land.use)

#               primary     secondary     cropland        pasture
colors.lu <- c("darkgreen", "green3", "darkgoldenrod2", "darkorange2")

# ---------
# Dominant land use
# ---------
saveGIF(
# for(i in min(land.use$year):max(land.use$year)){
for(i in 1475:max(land.use$year)){
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
}, movie.name="Dominant_LandUse_1475-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# ---------

# ---------
# Fractional Land Usage
# ---------
saveGIF(
# for(i in min(land.use$year):max(land.use$year)){
for(i in 1475:max(land.use$year)){
print(paste0("Graphing Year : ", i))

plot.prim <- ggplot(data= land.use[land.use$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= gothr)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(land.use$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(land.use$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "darkgreen"), limits=c(0,1), name="Frac. Primary") +
    ggtitle(paste0("Fraction Primary : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())

plot.secd <- ggplot(data= land.use[land.use$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= gsecd)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(land.use$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(land.use$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "green3"), limits=c(0,1), name="Frac. Secondary") +
    ggtitle(paste0("Fraction Primary : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())

plot.crop <- ggplot(data= land.use[land.use$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= gcrop)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(land.use$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(land.use$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "darkgoldenrod2"), limits=c(0,1), name="Frac. Crop") +
    ggtitle(paste0("Fraction Primary : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())

plot.past <- ggplot(data= land.use[land.use$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= gpast)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(land.use$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(land.use$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray50", "darkorange2"), limits=c(0,1), name="Frac. Pasture") +
    ggtitle(paste0("Fraction Primary : ", i)) +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
    
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
print(plot.prim, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.secd, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plot.crop, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.past, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
}, movie.name="Fraction_LandUse_1475-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=800, ani.width=1000)# dev.off()
# ---------

# -----------------------

# -----------------------
# 2) Doing total biomass harvest
# -----------------------
harvest.b<- c("gsbh1", "gsbh2", "gsbh3", "gvbh1", "gvbh2")
# get a spatial file for each PFT
gsbh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh1")))
gsbh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh2")))
gsbh3 <- stack(file.path(lulcc.out, dir(lulcc.out, "gsbh3")))
gvbh1 <- stack(file.path(lulcc.out, dir(lulcc.out, "gvbh1")))
gvbh2 <- stack(file.path(lulcc.out, dir(lulcc.out, "gvbh2"))) 
gvbh2[[1]]

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

years <- as.numeric(substr(df1[,2],2,nchar(as.vector(df1[,2]))))-1+850
harvest.bm <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gsbh1=df1[,1], gsbh2=df2[,1], gsbh3=df3[,1], gvbh1=df4[,1], gvbh2=df5[,1])
harvest.bm$TotalHarvest <- rowSums(harvest.bm[,harvest.b])
summary(harvest.bm)

# Put Carbon on a Better scale
# assume 0.5 x 0.5 - degree roughly equal to 50 x 50 km
# 50 x 50 km = 5e4*5e4 m2
harvest.bm$TotalHarvest_MgHa <- harvest.bm$TotalHarvest*1e-3/((5e4*5e4)/10000)
summary(harvest.bm)

summary(harvest.bm[harvest.bm$TotalHarvest_MgHa>2,])
dim(harvest.bm[harvest.bm$TotalHarvest_MgHa>2.3,])
quantile(harvest.bm$TotalHarvest_MgHa, 0.999)


saveGIF(
# for(i in min(harvest.bm$year):max(harvest.bm$year)){
for(i in 1475:max(harvest.bm$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data=harvest.bm[harvest.bm$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= TotalHarvest_MgHa)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(harvest.bm$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(harvest.bm$lat), expand=c(0,0), name="Latitude (degrees)") +
    # scale_fill_gradientn(name="Forest Harvest, Mg/Ha/yr", colours=c("gray50","darkgreen"), limits=c(0,quantile(harvest.bm$TotalHarvest_MgHa,0.95))) +
    scale_fill_gradientn(name="Forest Harvest, Mg/Ha/yr", colours=c("gray50","darkgreen"), limits=range(harvest.bm$TotalHarvest_MgHa)) +
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
}, movie.name="ForestHarvest_Biomass_1475-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)
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

years <- as.numeric(substr(df1[,2],2,nchar(as.vector(df1[,2]))))-1+850
harvest.area <- data.frame(lat=lat.lon$lat, lon=lat.lon$lon, year=years, gfsh1=df1[,1], gfsh2=df2[,1], gfsh3=df3[,1], gfvh1=df4[,1], gfvh2=df5[,1])
harvest.area$TotalHarvest <- rowSums(harvest.area[,harvest.f])
summary(harvest.area)

# png(file.path(lulcc.out, "lulcc_PotentialVegetation_Dominant.png"), height=600, width=800)
saveGIF(
# for(i in min(harvest.area$year):max(harvest.area$year)){
for(i in 1475:max(land.use$year)){
# for(i in 2000:2005){
print(paste0("Graphing Year : ", i))
print(
ggplot(data= harvest.area[harvest.area$year==i,]) +
    geom_raster(aes(x=lon, y=lat, fill= TotalHarvest)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(harvest.area$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(harvest.area$lat), expand=c(0,0), name="Latitude (degrees)") +
    # scale_fill_gradientn(name="Forest Harvest, Fractional Area/yr", colours=c("gray50","darkgreen"), limits=c(0,quantile(harvest.area$TotalHarvest,0.95))) +
    scale_fill_gradientn(name="Forest Harvest, Fractional Area/yr", colours=c("gray50","darkgreen"), limits=range(harvest.area$TotalHarvest)) +
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
}, movie.name="ForestHarvest_Area_1475-2010.gif", interval=0.3, nmax=10000, autobrowse=F, autoplay=F, ani.height=600, ani.width=800)# dev.off()
# -----------------------





# ----------------------------------------------
