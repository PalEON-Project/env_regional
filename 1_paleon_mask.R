# ----------------------------------------------
# Script to generate the PalEON domain mask (land only)
# Christine Rollinson, crollinson@gmail.com
# Original: 16 September, 2015
#
# --------------
# Mask Proceedure
# --------------
# 1) open met file & convert to binary presence/absence (NAs)
# 2) export to ncdf
# 3) Plot and save an image of the mask
#
# Note: The met data was originally masked to precipf_0850_01_01
#       inside the met processing scripts due to inconsistencies 
#       in cell presence/absence throughout the data.  See
#       https://github.com/PalEON-Project/met_regional/blob/master/
#       5_rewrite_timestamps.R
# Note: A decent amount of coastline is missing because this mask is
#       ultimately derived from the CCSM4 met data, which is at a 
#       coarser spatial resolution than we are running the models at.
# --------------
#
# ----------------------------------------------


# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4)
# dir.met <- "/projectnb/dietzelab/paleon/met_regional/phase2_met_drivers.v1/precipf/precipf_0850_01_01.nc"
path.met <- "/projectnb/dietzelab/paleon/met_regional/bias_corr/corr_timestamp_v2/precipf/precipf_0850_01.nc" # temporary until all qa/qc is done
path.env <- "/projectnb/dietzelab/paleon/env_regional/phase2_env_drivers_v1/domain_mask"
env.name <- "paleon_domain.nc"

# Create the driver folder if it doesn't already exist
if(!dir.exists(path.env)) dir.create(path.env)

# dummy fill values
fillv <- 1e30

# ----------------------------------------------


# ----------------------------------------------
# Create & export mask
# ----------------------------------------------
# Get met data
nc.met <- nc_open(path.met)
met.mask <- ncvar_get(nc.met, "precipf")[,,1] # Note, we only need 1 layer for this to work
lat      <- ncvar_get(nc.met, "lat")
lon      <- ncvar_get(nc.met, "lon")
nc_close(nc.met)

# Convert non-NA values to 1
met.mask[!is.na(met.mask)] <- 1

# Export as a netcdf file
# # Set up dims & var defs
nc_variable_long_name=paste('PalEON Domain Spatial Mask', sep='')
nc_variable_units='binary, unitless'
dimY <- ncdim_def( "lat", "longitude: degrees", lat )
dimX <- ncdim_def( "lon", "latitude: degrees", lon )
nc_var  <- ncvar_def("domain",nc_variable_units, list(dimX,dimY), fillv, longname=nc_variable_long_name, prec="float")

# Create & fill the file
domain.mask <- nc_create(file.path(path.env, env.name), nc_var)
ncatt_put(domain.mask, '0', 'description', 'PalEON spatial domain mask for MIP Phase 2')
ncvar_put(domain.mask, nc_var, met.mask)
nc_close(domain.mask)
# ----------------------------------------------

# ----------------------------------------------
# Open & view the mask
# ----------------------------------------------
# Load Spatial & Graphing Libraries
library(raster); library(ggplot2); library(grid)

# Get state outlines
paleon.states <- map_data("state")

domain <- raster(file.path(path.env, env.name))

domain.df <- data.frame(rasterToPoints(domain))
names(domain.df) <- c("lon", "lat", "mask")
domain.df$mask   <- as.factor(domain.df$mask) # making the mask categorical
summary(domain.df)

png(file.path(path.env, "paleon_domain.png"), height=400, width=800)
ggplot(data=domain.df) +
    geom_raster(aes(x=lon, y=lat), fill="gray50") +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(domain.df$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(domain.df$lat), expand=c(0,0), name="Latitude (degrees)") +
    ggtitle("PalEON Spatial Domain Mask") +
    coord_equal(ratio=1) +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank())
dev.off()