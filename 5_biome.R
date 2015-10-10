# ----------------------------------------------
# Script to extract 
# Christine Rollinson, crollinson@gmail.com
# Original: 16 September, 2015
#
# --------------
# Mask Proceedure
# --------------
# 1) Open & extract .nc data
# 2) Convert to spatial domain
# 3) Fill in gaps where MsTMIP doesn't cover
# 4) Crop to paleon domain extent
# 5) Resample to make sure grid cells line up & gapfill
# 6) Mask to make sure we have consistent land/water coverage
# 7) Aggregate up to PalEON Biomes & PFTs
# 8) Convert to .nc file for consistent convention
# --------------
#
# ----------------------------------------------


# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4); library(raster); library(rgdal)
library(car)

setwd("~/Desktop/Research/PalEON_CR/env_regional/")
paleon.mask <- "env_paleon/domain_mask/paleon_domain.nc"
biom.path <- "env_drivers_raw/biome/NACP_MSTMIP_MODEL_DRIVER/data"
biom.out <- "env_paleon/biome"

# Create the driver folder if it doesn't already exist
if(!dir.exists(biom.out)) dir.create(biom.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export masked data
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)

# 1) load the raw data
# Note: because of how the .nc file is set up, we have to do some extra steps to make it a spatial
#       file (which just makes my life easier, but could be avoided)
biome.file <- nc_open(file.path(biom.path, "mstmip_driver_global_hd_biome_potveg_v1.nc4"))
summary(biome.file$var)

biome.lat  <- ncvar_get(biome.file, "lat")
biome.lon  <- ncvar_get(biome.file, "lon")
biome.type <- ncvar_get(biome.file, "biome_type")
biome.frac <- ncvar_get(biome.file, "biome_frac")

nc_close(biome.file)

dim(biome.type)
dim(biome.frac)


# 2) Translating the nc arrays into rasters to make visualization & masking easier
# biome type first
type.raster <- raster(t(biome.type), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat), crs=projection(paleon))
type.raster[type.raster==0] <- NA
# plot(type.raster)

# PFT fraction next;
frac.raster <- stack(raster(t(biome.frac[,,1]), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat), crs=projection(paleon)))
for(i in 2:dim(biome.frac)[3]){
	frac.raster[[i]] <- raster(t(biome.frac[,,i]), xmn=min(biome.lon), xmx=max(biome.lon), ymn=min(biome.lat), ymx=max(biome.lat))
}
names(frac.raster) <- paste0("X", 0:(dim(biome.frac)[3]-1))
# plot(frac.raster[[2]])

# 3) create the smoothed datasets to fill in MsTMIP Gaps
#    Note: right now, we're going to skip over the fraction PFT and see how it pans out
type.fill <- focal(type.raster, w=matrix(1,nrow=3, ncol=3), fun=modal, na.rm=T)

# frac.fill <- frac.raster
# for(i in 1:nlayers(frac.raster)){
	# frac.fill[[i]] <- focal(frac.raster[[i]], w=matrix(1,nrow=3, ncol=3), fun=mean, na.rm=T)
# }
# plot(frac.fill[[2]])

# 4) crop the region
type.crop <- crop(type.raster, paleon)
type.fill <- crop(type.fill, paleon)

frac.crop <- crop(frac.raster, paleon)
# frac.fill <- crop(frac.fill, paleon)

# plot(type.crop)
# plot(paleon, add=T, alpha=0.5)

# plot(type.fill)
# plot(paleon, add=T, alpha=0.5)

# 5) Resample type & fraction; gapfill biome type
# Note: the biome type needs the nearest neighbor method to keep it discreet
type.resamp <- resample(type.crop, paleon, method="ngb")
type.fill   <- resample(type.fill, paleon, method="ngb")

type.resamp[is.na(type.resamp)] <- type.fill[is.na(type.resamp)]

frac.resamp <- resample(frac.crop, paleon, method="ngb")
# frac.fill   <- resample(frac.fill, paleon)

# 6) Mask the files
type.biome  <- mask(type.resamp, paleon)
# plot(type.biome)
# plot(paleon, add=T, alpha=0.4)

frac.biome  <- mask(frac.resamp, paleon)
# drop the water
frac.biome  <- dropLayer(frac.biome,1)
# plot(frac.biome[[2]])
# plot(paleon, add=T, alpha=0.4)

# 8) Reclassify SYNMAP codes into ones that make sense for PalEON
# 8a) Biome Reclassificaiton following Kevin Schaefer's Crosswalk for SibCASA
# 8b) PFT Crosswalk following Poulter et al 2015.

# --------------------
# 8a) recoding the SYNMAP Biomes to PalEON (See Synmap_Crosswalk.xlsx for more details)
# PalEON Code    Description                                       SYNMAP CODES   
#      1         Broadleaf-Deciduous Forest                        5
#      2         Mixed Broadleaf-Deciduous & Needleleaf-Evergreen  6, 9
#      3         Needleleaf-Evergreen Forest                       1, 19
#      4         Savanna & Shrublands                              10, 14, 23, 27, 37, 38, 40
#      5         Grasslands                                        31, 32, 33, 36, 39, 41, 42, 44, 47
# --------------------

type.values <- as.factor(values(type.biome))
summary(type.values)

type.paleon <- recode(type.values, " '5'='1'; 
                                     '6'='2';  '9'='2';
                                     '1'='3'; '19'='3';
                                    '10'='4'; '14'='4'; '23'='4'; '27'='4'; '37'='4'; '38'='4'; '40'='4';
                                    '31'='5'; '32'='5'; '33'='5'; '36'='5'; '39'='5'; '41'='5'; '42'='5'; '44'='5'; '47'='5'")
summary(type.paleon)

paleon.biome <- setValues(type.biome, type.paleon)
plot(paleon.biome)
# --------------------

# --------------------
# 8b) Reorganizing SYNMAP biomes into common model PFTs following
#     Poulter et al. 2015
#     PFT   Code    Form     Leaf       Phenology
#      1    TBrEv   Tree     Broad      Evergreen
#      2    TBrDe   Tree     Broad      Deciduous
#      3    TNeEv   Tree     Needle     Evergreen
#      4    TNeDe   Tree     Needle     Deciduous
#      5    SBrEv   Shrub    Broad      Evergreen
#      6    SBrDe   Shrub    Broad      Deciduous
#      7    SNeEv   Shrub    Needle     Evergreen
#      8    SNeDe   Shrub    Needle     Deciduous
#      9    Grass   Grass      -            - 
#     10    Bare    Bare       -            - 
#
# Steps
# 8.1) Reclassify SYNMAP biomes into LCCS biomes
# 8.2) Convert LCCS biomes to PFTs
# --------------------
# Load Crosswalk dataframes
biome2pft <- read.csv("Poulter2015_Table2_LCCS_Crosswalk.csv")
biome2pft$LCCS_Code <- as.factor(paste0("X", biome2pft$LCCS_Code))
biome2pft[is.na(biome2pft)] <- 0
summary(biome2pft)

synmap2lccs <- read.csv("SYNMAP_LCCS_Crosswalk.csv", na.strings="-")
# synmap2lccs$LCCS_Code1 <- as.factor(paste0("X", synmap2lccs$LCCS_Code1))
# synmap2lccs$LCCS_Code2 <- as.factor(ifelse(!is.na(synmap2lccs$LCCS_Code2), paste0("X", synmap2lccs$LCCS_Code2), NA))
summary(synmap2lccs)


# Create some blank rasters
# LCCS Biomes
lccs.biomes <- brick(ncols=ncol(frac.biome), nrows=nrow(frac.biome), nl=nrow(biome2pft),
                     xmn=xmin(frac.biome), xmx=xmax(frac.biome), 
                     ymn=ymin(frac.biome), ymx=ymax(frac.biome),
                     crs=projection(frac.biome))
lccs.biomes <- setValues(lccs.biomes, 0)
names(lccs.biomes) <- biome2pft$LCCS_Code
lccs.biomes <- mask(lccs.biomes, paleon)
# plot(lccs.biomes)

paleon.pfts <- brick(ncols=ncol(frac.biome), nrows=nrow(frac.biome), nl=ncol(biome2pft)-2,
                     xmn=xmin(frac.biome), xmx=xmax(frac.biome), 
                     ymn=ymin(frac.biome), ymx=ymax(frac.biome),
                     crs=projection(frac.biome))
paleon.pfts <- setValues(paleon.pfts, 0)
paleon.pfts <- mask(paleon.pfts, paleon)
names(paleon.pfts) <- names(biome2pft)[3:ncol(biome2pft)]


# 8.1) Convert SYNMAP Biomes to LCCS (used in Poulter 2015)
for(i in 1:nrow(synmap2lccs)){
	if(!is.na(synmap2lccs[i,"LCCS_Code2"])){
		biome1 <- paste0("X", synmap2lccs[i,"LCCS_Code1"])
		biome2 <- paste0("X", synmap2lccs[i,"LCCS_Code2"])
		lccs.biomes[[biome1]] <- lccs.biomes[[biome1]] + frac.biome[[i]]*0.5
		lccs.biomes[[biome2]] <- lccs.biomes[[biome2]] + frac.biome[[i]]*0.5
	} else {
		biome1 <- paste0("X", synmap2lccs[i,"LCCS_Code1"])
		lccs.biomes[[biome1]] <- lccs.biomes[[biome1]] + frac.biome[[i]]
	}
}
lccs.biomes

plot(lccs.biomes, zlim=c(0,1))


# 8.2) Convert LCCS Biomes to PFTs following Poulter 2015
for(i in 1:nrow(biome2pft)){
	for(j in names(paleon.pfts)){
		paleon.pfts[[j]] <- paleon.pfts[[j]] + lccs.biomes[[i]]*biome2pft[i,j]*.01
	}
}
paleon.pfts
plot(paleon.pfts, zlim=c(0,1))

# Normalizing so that all PFT coverage sums to 1
paleon.sum <- sum(paleon.pfts)
paleon.sum

for(i in 1:nlayers(paleon.pfts)){
	paleon.pfts[[i]] <- paleon.pfts[[i]]/paleon.sum
}
paleon.pfts
names(paleon.pfts) <- names(biome2pft)[3:(ncol(biome2pft))]
plot(paleon.pfts, zlim=c(0,1))
# --------------------

# Writing as netcdf rasters
writeRaster(paleon.biome, file.path(biom.out, "biome_potential_vegtype_biome.nc"), format="CDF", overwrite=T, varname="biome", varunit="categorical", longname="Biome Classification")

writeRaster(paleon.pfts, file.path(biom.out, "biome_potential_vegtype_pft_fraction.nc"), format="CDF", overwrite=T, varname="pft_frac", varunit="fraction", longname="PFT Fractional Coverage (area)", zname="pft", zunit="categorical")
# ----------------------------------------------









# ----------------------------------------------
# Graphing Biome & PFT coverage
# ----------------------------------------------
library(ggplot2); library(grid); library(raster)
paleon.states <- map_data("state")

# ------------------
# 1) Biome classification
# ------------------
paleon.biome <- raster(file.path(biom.out, "biome_potential_vegtype_biome.nc"))
plot(paleon.biome)

biome.df <- data.frame(rasterToPoints(paleon.biome))
names(biome.df) <- c("lon", "lat", "biome")
biome.df$biome <- as.factor(biome.df$biome)
levels(biome.df$biome) <- c("1. Forest, Hardwood", "2. Forest, Mixed", "3. Forest, Conifer", "4. Savanna & Shrublands", "5. Grasslands")
summary(biome.df)

biome.colors <- c("green3", "darkseagreen3", "darkgreen", "darkorange2", "goldenrod2")

png(file.path(biom.out, "Biome_PotentialVegetation_Biome.png"), height=600, width=1000)
ggplot(biome.df) +
	geom_raster(aes(x=lon, y=lat, fill=biome))	 +
	geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
	ggtitle("Dominant Biome Type (modified from SYNMAP)") +
    scale_x_continuous(limits=range(biome.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(biome.df$lat), expand=c(0,0), name="Latitude (degrees)") +
	scale_fill_manual(values=biome.colors, name="Biome") +
    coord_equal(ratio=1) +
	# theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()
# ------------------


# ------------------
# 1) pft classification
#     PFT   Code    Form     Leaf       Phenology
#      1    TBrEv   Tree     Broad      Evergreen
#      2    TBrDe   Tree     Broad      Deciduous
#      3    TNeEv   Tree     Needle     Evergreen
#      4    TNeDe   Tree     Needle     Deciduous
#      5    SBrEv   Shrub    Broad      Evergreen
#      6    SBrDe   Shrub    Broad      Deciduous
#      7    SNeEv   Shrub    Needle     Evergreen
#      8    SNeDe   Shrub    Needle     Deciduous
#      9    Grass   Grass      -            - 
#     10    Bare    Bare       -            - 
# ------------------
paleon.pft <- stack(file.path(biom.out, "biome_potential_vegtype_pft_fraction.nc"))
names(paleon.pft)[1:9] <- paste0("X0",1:9)
paleon.pft
plot(paleon.pft, zlim=c(0,1))

pft.df1 <- data.frame(rasterToPoints(paleon.pft))
names(pft.df1)[1:2] <- c("lon", "lat")

pft.df <- stack(pft.df1[,3:ncol(pft.df1)])
names(pft.df) <- c("fraction", "PFT")

pft.df[,c("lon", "lat")] <- pft.df1[,c("lon", "lat")]
levels(pft.df$PFT) <- c("1. Tree, Broad Evg", "2. Tree, Broad Decid", "3. Tree, Needle Evg", "4. Tree, Needle Decid", "5. Shrub, Broad Evg", "6. Shrub, Broad Decid", "7. Shrub, Needle Evg", "8. Shrub, Needle Decid", "Grass (C3)", "10. Bare/Non-Vegetated")
summary(pft.df)

png(file.path(biom.out, "Biome_PotentialVegetation_PFT_Fraction.png"), height=600, width=1000)
ggplot(pft.df) + facet_wrap(~PFT) +
	geom_raster(aes(x=lon, y=lat, fill=fraction))	 +
	geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
	ggtitle("PFT Fraction (modified from SYNMAP)") +
    scale_x_continuous(limits=range(pft.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(pft.df$lat), expand=c(0,0), name="Latitude (degrees)") +
	scale_fill_gradientn(colours=c("gray35", "green3")) +
    coord_equal(ratio=1) +
	# theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()
# ------------------
# ----------------------------------------------

