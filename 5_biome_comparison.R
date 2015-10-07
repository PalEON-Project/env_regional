# ----------------------------------------------
# Script to extract 
# Christine Rollinson, crollinson@gmail.com
# Original: 16 September, 2015
#
# --------------
# Mask Proceedure
# --------------
# 1) Open raw potential veg file
# 2) Calculate percent PFT over a 5˚window
# 3) Crop to paleon domain extent
# 4) resample to make sure grid cells line up
# 5) Get rid of sparse water values by using whatever had the greatest weight above
# 6) Mask to make sure we have consistent land/water coverage
# 7) Convert to .nc file for consistent convention
# --------------
#
# ----------------------------------------------


# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4); library(raster); library(rgdal)
library(ggplot2); library(grid)
library(car)

paleon.mask <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc"
biom.path <- "~/Desktop/Research/PalEON_CR/env_regional/env_drivers_raw/biome/ISLSCP_II_POT_VEG_961/data/potential_veg_hd.asc"
mstmip.path <- "~/Desktop/Research/PalEON_CR/env_regional/env_drivers_raw/biome/NACP_MSTMIP_MODEL_DRIVER/data"
biom.out <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/biome"

# Create the driver folder if it doesn't already exist
if(!dir.exists(biom.out)) dir.create(biom.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Dominant Biome Type
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)
# Get state outlines
paleon.states <- map_data("state")

biome.rf <- raster(biom.path, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
biome.rf


biome.mstmip <- nc_open(file.path(mstmip.path, "mstmip_driver_global_hd_biome_potveg_v1.nc4"))
summary(biome.mstmip$var)

biome.type <- ncvar_get(biome.mstmip, "biome_type")
biome.lat <- ncvar_get(biome.mstmip, "lat")
biome.lon <- ncvar_get(biome.mstmip, "lon")
summary(biome.type)
dim(biome.type)
plot(biome.rf)
plot(biome.mstmip)

biome.raster <- raster(t(biome.type), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat), crs=projection(biome.rf))
plot(biome.raster)


rf.crop <- crop(biome.rf, paleon)
plot(rf.crop)

ms.crop <- crop(biome.raster, paleon)
plot(ms.crop)

rf.df <- data.frame(rasterToPoints(rf.crop))
names(rf.df) <- c("lon", "lat", "pot_veg")
rf.df$pot_veg <- as.factor(rf.df$pot_veg)
rf.df$pft <- as.factor(rf.df$pot_veg)
levels(rf.df$pft) <- c("NA", "conifer, temp", "decid., temp", "conifer, boreal", "mixed", "savanna", "grass")
rf.df <- rf.df[!(rf.df$pot_veg==0),]
summary(rf.df)


ms.df <- data.frame(rasterToPoints(ms.crop))
names(ms.df) <- c("lon", "lat", "pot_veg")
ms.df$pot_veg <- as.factor(ms.df$pot_veg)
ms.df$pft <- ms.df$pot_veg
# 0, 1, 5, 9, 14, 19, 23, 27, 37, 38, 41
levels(ms.df$pft) <- c("NA", "tree: needle, evg", "tree: broad, decid,", "tree: mixed, mixed", "tree/shrub: broad, decid", "tree/grass: needle, evg", "tree/grass: broad, decid", "tree/grass: mixed, mixed", "shrubs", "shrubs/grass", "grass")
ms.df <- ms.df[!(ms.df$pot_veg==0),]
summary(ms.df)

#                   conif.temp         decid.temp      conif.bor    mixed            grass
biome.colors1 <- c("darkolivegreen4", "springgreen3", "darkgreen", "darkseagreen3", "darkorange2", "goldenrod2")

#
biome.colors2 <- c("darkolivegreen4", "springgreen3", "darkseagreen3", "darkorange2", "darkorange2", "darkorange2", "darkorange2", "darkorange2", "goldenrod2", "goldenrod2")

png("Ramankutty_Foley.png", width=800, height=400)
ggplot(data=rf.df) +
	geom_raster(aes(x=lon, y=lat, fill=pft))	 +
	geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
	ggtitle("Ramankutty & Foley 2010") +
    scale_x_continuous(limits=range(rf.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(rf.df$lat), expand=c(0,0), name="Latitude (degrees)") +
	scale_fill_manual(values=biome.colors1) +
    coord_equal(ratio=1) +
	# theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()

png("Synmap.png", width=800, height=400)
ggplot(data=ms.df) +
	geom_raster(aes(x=lon, y=lat, fill=pft))	 +
	geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
	ggtitle("Synmap (MsTMIP)") +
    scale_x_continuous(limits=range(ms.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(ms.df$lat), expand=c(0,0), name="Latitude (degrees)") +
	scale_fill_manual(values=biome.colors2) +
    coord_equal(ratio=1) +
	# guide_legend(nrow=2) +
	# theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()
# ----------------------------------------------



# ----------------------------------------------
# Fractional Biome coverage (Synmap )
# ----------------------------------------------
# ----------------------------------------------
# Dominant Biome Type
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)
# Get state outlines
paleon.states <- map_data("state")

biome.mstmip <- nc_open(file.path(mstmip.path, "mstmip_driver_global_hd_biome_potveg_v1.nc4"))
summary(biome.mstmip$var)

biome.frac <- ncvar_get(biome.mstmip, "biome_frac")
biome.lat <- ncvar_get(biome.mstmip, "lat")
biome.lon <- ncvar_get(biome.mstmip, "lon")
# summary(biome.type)
dim(biome.frac)
# plot(biome.rf)
# plot(biome.mstmip)


biome.raster <- raster(t(biome.frac[,,1]), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat))
plot(biome.raster)
biome.raster <- stack(biome.raster)
biome.raster

for(i in 2:dim(biome.frac)[3]){
	biome.raster[[i]] <- raster(t(biome.frac[,,i]), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat))
}
names(biome.raster) <- paste0("X", 1:48)
biome.raster
plot(biome.raster)

ms.crop <- crop(biome.raster, paleon)
ms.crop
plot(ms.crop[[1:12]])
plot(ms.crop[[13:24]])
plot(ms.crop[[25:36]])
plot(ms.crop[[37:48]])

# ms.df <- data.frame(rasterToPoints(ms.crop))
# names(ms.df) <- c("lon", "lat", "pot_veg")
# ms.df$pot_veg <- as.factor(ms.df$pot_veg)
# ms.df$pft <- ms.df$pot_veg
# # 0, 1, 5, 9, 14, 19, 23, 27, 37, 38, 41
# levels(ms.df$pft) <- c("NA", "tree: needle, evg", "tree: broad, decid,", "tree: mixed, mixed", "tree/shrub: broad, decid", "tree/grass: needle, evg", "tree/grass: broad, decid", "tree/grass: mixed, mixed", "shrubs", "shrubs/grass", "grass")
# ms.df <- ms.df[!(ms.df$pot_veg==0),]
# summary(ms.df)

# #
# biome.colors2 <- c("darkolivegreen4", "springgreen3", "darkseagreen3", "darkorange2", "darkorange2", "darkorange2", "darkorange2", "darkorange2", "goldenrod2", "goldenrod2")


# png("Synmap.png", width=800, height=400)
ggplot(data=ms.df) +
	geom_raster(aes(x=lon, y=lat, fill=pft))	 +
	geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
	ggtitle("Synmap (MsTMIP)") +
    scale_x_continuous(limits=range(ms.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(ms.df$lat), expand=c(0,0), name="Latitude (degrees)") +
	scale_fill_manual(values=biome.colors2) +
    coord_equal(ratio=1) +
	# guide_legend(nrow=2) +
	# theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
# dev.off()
# ----------------------------------------------
