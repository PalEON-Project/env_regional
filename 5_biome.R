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

paleon.mask <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc"
biom.path <- "~/Dropbox/PalEON_CR/env_regional/env_drivers_raw/biome/ISLSCP_II_POT_VEG_961/data/potential_veg_hd.asc"
biom.out <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/biome"

# Create the driver folder if it doesn't already exist
if(!dir.exists(biom.out)) dir.create(biom.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export mask
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)

# 1) load the raw data
# Assume that the long/lat is WGS84 and lines up with PalEON
biome.raw <- raster(biom.path, crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
biome.raw
# plot(biome.raw) # Just to double check that things got read in okay

# 2) Calculate the percent coverage of each PFT in a given cell based on a 5˚ x 5˚ window
# 2.1) First, put each PFT in its own layer
pfts <- stack(biome.raw)
projection(pfts) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
pfts$pft4  <- biome.raw; pfts$pft4 [!pfts$pft4 == 4] <- 0
pfts$pft5  <- biome.raw; pfts$pft5 [!pfts$pft5 == 5] <- 0
pfts$pft6  <- biome.raw; pfts$pft6 [!pfts$pft6 == 6] <- 0
pfts$pft8  <- biome.raw; pfts$pft8 [!pfts$pft8 == 8] <- 0
pfts$pft9  <- biome.raw; pfts$pft9 [!pfts$pft9 == 9] <- 0
pfts$pft10 <- biome.raw; pfts$pft10[!pfts$pft10==10] <- 0
# plot(pfts)

# 2.2) assign an importance on a scale of 0-1; note: becuase the cells are currently
#    being evaluated as numeric, the weight needs to be 1/pftID
f.pfts <- stack(biome.raw)
projection(f.pfts) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
f.pfts$pft4  <- focal(pfts$pft4, w=matrix(1/4,nrow=11, ncol=11), fun=mean, na.rm=T)
f.pfts$pft5  <- focal(pfts$pft5, w=matrix(1/5,nrow=11, ncol=11), fun=mean, na.rm=T)
f.pfts$pft6  <- focal(pfts$pft6, w=matrix(1/6,nrow=11, ncol=11), fun=mean, na.rm=T)
f.pfts$pft8  <- focal(pfts$pft8, w=matrix(1/8,nrow=11, ncol=11), fun=mean, na.rm=T)
f.pfts$pft9  <- focal(pfts$pft9, w=matrix(1/9,nrow=11, ncol=11), fun=mean, na.rm=T)
f.pfts$pft10 <- focal(pfts$pft10, w=matrix(1/10,nrow=11, ncol=11), fun=mean, na.rm=T)
# plot(f.pfts)

# 2.3) relative the 0-1 scale
rf.pfts <- stack(biome.raw)
projection(rf.pfts) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
rf.pfts$pft4 <- f.pfts$pft4/sum(f.pfts[[2:7]])
rf.pfts$pft5 <- f.pfts$pft5/sum(f.pfts[[2:7]])
rf.pfts$pft6 <- f.pfts$pft6/sum(f.pfts[[2:7]])
rf.pfts$pft8 <- f.pfts$pft8/sum(f.pfts[[2:7]])
rf.pfts$pft9 <- f.pfts$pft9/sum(f.pfts[[2:7]])
rf.pfts$pft10 <- f.pfts$pft10/sum(f.pfts[[2:7]])
# plot(rf.pfts)

# 3) crop raw data (extents need to be same before masking)
biom.crop <- crop(rf.pfts, extent(paleon))
# plot(biom.crop[[1]])
# plot(biom.crop[[2:7]], zlim=c(0,1))

# plot(biom.crop[[1]])
# plot(paleon, add=T, alpha=0.8)

# lets make sure everythign sums to 1
check <- sum(biom.crop[[2:7]])
check

# 4) Resample -- It would probably be okay to leave this out, but just in case...
biom.resamp <- resample(biom.crop, paleon)
# plot(biom.resamp[[1]])
# plot(paleon, add=T, alpha=0.8)

# 5) mask raw data
biom.mask <- mask(biom.resamp, paleon)
# plot(biom.mask[[1]])
# plot(paleon, add=T, alpha=0.8)

# plot(biom.mask[[2:7]])



# 6) Get rid of water values in the first layer of potential veg
#    Since we've done the fraction of each PFT, lets just fill it with whichever
#    is greatest
# A vector of PFTs in order
pft.list <- c(4,5,6,8,9,10)

# 6.1) Find the pft values at each of the water cells
pft.max <- biom.mask[[2:7]][biom.mask$potential_veg_hd==0]
# 6.2) find out which pft has the highest value is the max 
pft.replace <- apply(pft.max, 1, function(x,...){pft.list[which(x == max(x))]})

biom.final <- biom.mask
biom.final$potential_veg_hd[biom.final$potential_veg_hd==0] <- pft.replace
# plot(biom.final[[1]])
# plot(biom.final[[2:7]])

# 7) convert to format for .nc
# # Saving as an array so that all the proper variables names get written to netcdf
writeRaster(biom.final[[1]], file.path(biom.out, "biome_potential_vegtype_dominant.nc"), format="CDF", overwrite=T, varname="potential_veg_dom", varunit="cateogrical", longname="Vegetation Type", zname="vegtype_names", zunit="cateogrical")
writeRaster(biom.final[[2:7]], file.path(biom.out, "biome_potential_vegtype_relative.nc"), format="CDF", overwrite=T, varname="VegType_Fraction", varunit="fraction", longname="Vegetation Type", zname="Veg_Type")
# ----------------------------------------------

# ----------------------------------------------
# Open & view the potential veg map
# ----------------------------------------------
# Load Spatial & Graphing Libraries
library(raster); library(ggplot2); library(grid)
library(car)

# Get state outlines
paleon.states <- map_data("state")

# 1) Doing the dominant PFT 
biome.main <- raster(file.path(biom.out, "biome_potential_vegtype_dominant.nc"))
plot(biome.main)

biome.df <- data.frame(rasterToPoints(biome.main))
names(biome.df)     <- c("lon", "lat", "biome")
# biome.df$biome      <- biome.df$biome) # making the mask categorical
# Adding in the text code for biome (note: only the relevant ones)
biome.df$biome.name   <- as.ordered(biome.df$biome)
levels(biome.df$biome.name) <- c(" 4 = Temp. Needle. Evg.", " 5 = Temp. Decid.", " 6 = Boreal Needle. Evg", " 8 = Mixed", " 9 = Savanna", "10 = Grassland")
summary(biome.df)

biome.colors <- c("darkolivegreen4", "springgreen3", "darkgreen", "darkseagreen3", "darkorange2", "goldenrod2")

png(file.path(biom.out, "Biome_PotentialVegetation_Dominant.png"), height=600, width=800)
ggplot(data=biome.df) +
    geom_raster(aes(x=lon, y=lat, fill=biome.name)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(biome.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(biome.df$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_manual(values=biome.colors) +
    ggtitle("PalEON Potential Vegetation ('Biome' Driver)") +
    guides(fill=guide_legend(direction="horizontal", nrow=2, title="Biome") )+
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()

# --------------------
# 2) Doing the relative composition of each PFT
biome.pft <- stack(file.path(biom.out, "biome_potential_vegtype_relative.nc"))
biome.df <- data.frame(rasterToPoints(biome.pft))
names(biome.df)     <- c("lon", "lat", "pft04", "pft05", "pft06", "pft08", "pft09", "pft10")
summary(biome.df)

biome.df2 <- stack(biome.df[,3:ncol(biome.df)])
names(biome.df2) <- c("Fraction", "PFT")
biome.df2[,c("lon", "lat")] <- biome.df[,c("lon", "lat")]
biome.df2$biome.name   <- as.ordered(biome.df2$PFT)
levels(biome.df2$biome.name) <- c(" 4 = Temp. Needle. Evg.", " 5 = Temp. Decid.", " 6 = Boreal Needle. Evg", " 8 = Mixed", " 9 = Savanna", "10 = Grassland")
summary(biome.df2)

png(file.path(biom.out, "Biome_PotentialVegetation_Fraction.png"), height=600, width=1000)
ggplot(data=biome.df2) + facet_wrap(~biome.name) +
    geom_raster(aes(x=lon, y=lat, fill=Fraction)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(biome.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(biome.df$lat), expand=c(0,0), name="Latitude (degrees)") +
    # scale_fill_manual(values=biome.colors) +
    ggtitle("PalEON Potential Vegetation ('Biome' Driver)") +
    # guides(fill=guide_legend(direction="horizontal", title="Fraction") )+
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()
