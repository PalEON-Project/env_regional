# ----------------------------------------------
# Script to extract 
# Christine Rollinson, crollinson@gmail.com
# Original: 29 September, 2015
#
# --------------
# Mask Proceedure
# --------------
# 1) Open soil files
# 2) create smoothed product to fill in missing cells
# 3) Crop to paleon domain extent
# 4) Resample to half-degree resolution
# 5) Fill in values that are missing from the paleon mask
# 6) Mask to make sure we have consistent land/water coverage
# 7) Convert to .nc file for consistent convention
#
# Note: PalEON is using the North America data product because the documentation says
#       it's more accurate.  There are a few key data products missing, noticeably texture
#       and USDA soil code.  If someone needs them, we'll parse it from the global product
# --------------
#
# ----------------------------------------------


# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4); library(raster); library(rgdal)

paleon.mask <- "~/Desktop/Research/PalEON CR/env_regional/phase2_env_drivers_v1/domain_mask/paleon_domain.nc"
soil.path <- "~/Desktop/Research/PalEON CR/env_regional/env_drivers_raw/soil/NACP_MSTMIP_MODEL_DRIVER/data/mstmip_driver_na_qd_soil_v1.nc4"
soil.out <- "~/Desktop/Research/PalEON CR/env_regional/phase2_env_drivers_v1/soil"

# Create the driver folder if it doesn't already exist
if(!dir.exists(soil.out)) dir.create(soil.out)

# dummy fill values
fillv <- 1e30
# ----------------------------------------------


# ----------------------------------------------
# Create & export mask
# ----------------------------------------------
# Load the mask
paleon <- raster(paleon.mask)

# ------------------------
# 1) load the raw data: 
# ------------------------
#    Note: unlike other datasets, .nc4 isn't working as a raster, so we'll have to do a 
#          some extra steps
# Assume that the long/lat is WGS84 and lines up with PalEON
soil.na  <- nc_open(soil.path)
soil.vars <- data.frame(name=names(soil.na$var)[!names(soil.na$var) %in% c("crs", "lat_bnds", "lon_bnds")], longname=NA, units=NA)

soil.dat <- list()
soil.dat$lat <- ncvar_get(soil.na, "lat")
soil.dat$lon <- ncvar_get(soil.na, "lon")


for(v in soil.vars$name){
	soil.dat[[v]] <- t(ncvar_get(soil.na, v))
	soil.vars[soil.vars$name==v, "longname"] <- soil.na$var[[v]]$longname
	soil.vars[soil.vars$name==v, "units"   ] <- soil.na$var[[v]]$units
}
summary(soil.dat)
soil.vars

nc_close(soil.na)

# Create a single stack that we'll add to
soil <- stack(raster(soil.dat[[3]], xmn=min(soil.dat$lon), xmx=max(soil.dat$lon), ymn=min(soil.dat$lat), ymx=max(soil.dat$lat), crs=projection(paleon)))
names(soil) <- names(soil.dat)[3]
# soil
# plot(soil[[1]])

for(v in names(soil.dat)[4:length(soil.dat)]){
	soil[[v]] <- raster(soil.dat[[v]], xmn=min(soil.dat$lon), xmx=max(soil.dat$lon), ymn=min(soil.dat$lat), ymx=max(soil.dat$lat), crs=projection(paleon))
}
soil
plot(soil)
# ------------------------

# ------------------------
# 2) Create smoothed product to fill in the gaps
# ------------------------
soil.smooth <- soil
for(v in names(soil.smooth)){
  soil.smooth[[v]] <- focal(soil[[v]], w=matrix(1, nrow=3, ncol=3), fun=mean, na.rm=T)
}
# soil.smooth
# plot(soil.smooth[[1]])
# plot(soil[[1]])
# ------------------------

# ------------------------
# 2) crop raw data (extents need to be same before masking)
# ------------------------
soil.crop   <- crop(soil       , extent(paleon))
smooth.crop <- crop(soil.smooth, extent(paleon))
# soil.crop
# plot(soil.crop)
# plot(soil.crop[[1]])
# plot(smooth.crop[[1]])
# ------------------------

# ------------------------
# 3) resample soil data to match paleon mask
# ------------------------
soil.resamp <- resample(soil.crop, paleon, method="bilinear", na.rm=T)
soil.resamp

smooth.resamp <- resample(smooth.crop, paleon, method="bilinear", na.rm=T)
# smooth.resamp
# plot(soil.resamp)
# plot(soil.resamp[[1]])
# plot(smooth.resamp[[1]])

# plot(paleon)
# plot(soil.resamp[[1]], add=T)
# plot(soil.resamp[[1]], add=F)

# plot(smooth.resamp[[1]])
# plot(paleon)
# plot(smooth.resamp[[1]], add=T)
# ------------------------


# ------------------------
# 4) mask raw data
# ------------------------
soil.mask <- mask(soil.resamp, paleon)
# soil.mask
# plot(paleon)
# plot(soil.mask)

smooth.mask <- mask(smooth.resamp, paleon)
# smooth.mask

# plot(soil.mask[[1]])
# plot(smooth.mask[[1]])
# ------------------------

# ------------------------
# 5) replace missing values in mask with those from smooth.resamp
# ------------------------
soil.fill <- soil.mask
for(v in names(soil.fill)){
  soil.fill[[v]][is.na(soil.mask[[v]])] <- smooth.mask[[v]][is.na(soil.mask[[v]])]
}
plot(soil.mask[[1]])
plot(soil.fill[[1]])
# ------------------------


# ------------------------
# 6) save as netcdf: 
# ------------------------
# NOTE: Because the metadata is so varied & important in this one, it's important we reconvert it 
#       to a data frame and write the other way


# ----------
# 6a) Write each layer separately (easier mapping)
# ----------
for(v in names(soil.fill)){
  writeRaster(soil.fill[[v]], file.path(soil.out, paste0("paleon_soil_", v, ".nc")), format="CDF", overwrite=T, varname=v, varunit=soil.vars[soil.vars$name==v,"units"], longname=soil.vars[soil.vars$name==v,"longname"])
}
# ----------


# ----------
# 6b) Write single .nc file that mirrors MsTMIP format
# ----------
soil.list <- list()
soil.list$lon <- seq(xmin(soil.fill) + res(soil.fill)[1]/2, xmax(soil.fill), by=res(soil.fill)[1])
soil.list$lat <- seq(ymax(soil.fill) - res(soil.fill)[2]/2, ymin(soil.fill), by=-res(soil.fill)[2])

for(v in soil.vars$name){
	soil.list[[v]] <- t(as.array(soil.fill[[v]])[,,1])
}
summary(soil.list)
dim(soil.list[[3]])

dimY <- ncdim_def( "lat", "longitude: degrees", soil.list$lat)
dimX <- ncdim_def( "lon", "longitude: degrees", soil.list$lon)

nc_var <- list()
for(v in soil.vars$name){
nc_var[[v]]  <- ncvar_def(v, soil.vars[soil.vars$name==v,"units"], list(dimX,dimY), fillv, longname=soil.vars[soil.vars$name==v,"longname"], prec="double")
}
summary(nc_var)
summary(soil.list)

soil.final <- nc_create(file.path(soil.out, "paleon_soil.nc"), nc_var)
# ncatt_put(soil.final, '0', 'description', 'PalEON soil drivers extracted from MsTMIP North America')
for(v in names(nc_var)){
ncvar_put(soil.final, nc_var[[v]], soil.list[[v]])
}
nc_close(soil.final)
# ----------
# ------------------------
# ----------------------------------------------

# ----------------------------------------------
# Create some maps of soil characteristics
# ----------------------------------------------
# Load Spatial & Graphing Libraries
library(raster); library(ggplot2); library(grid); 
library(car)

# Get state outlines
paleon.states <- map_data("state")

# get a spatial file for each PFT
# soil.nc <- nc_open(file.path(soil.out, "paleon_soil.nc"))
# summary(soil.nc$var)
dir.soil <- dir(soil.out, ".nc")
dir.soil <- dir.soil[!dir.soil=="paleon_soil.nc"]

soil.raster <- stack(file.path(soil.out, dir.soil))
names(soil.raster) <- c("s_cec", "s_clay", "s_gravel", "s_oc", "s_ph", "s_ref_bulk", "s_sand", "s_silt", "soil_depth", "t_cec", "t_clay", "t_gravel", "t_oc", "t_ph", "t_ref_bulk", "t_sand", "t_silt")
soil.raster

# # Cheat way of just getting these done
# pdf(file.path(soil.out, "Soil_Properties_Ugly.pdf"))
# plot(soil.raster)
# dev.off()

soil <- data.frame(rasterToPoints(soil.raster))
names(soil) <- c("lon", "lat", "s_cec", "s_clay", "s_gravel", "s_oc", "s_ph", "s_ref_bulk", "s_sand", "s_silt", "soil_depth", "t_cec", "t_clay", "t_gravel", "t_oc", "t_ph", "t_ref_bulk", "t_sand", "t_silt")
summary(soil)


soil2 <- data.frame(lat=soil$lat, lon=soil$lon, layer=substr(stack(soil[,c("t_cec", "s_cec")])[,2],1,1),
                    cec      = stack(soil[,c("t_cec"     , "s_cec"     )])[,1], 
                    clay     = stack(soil[,c("t_clay"    , "s_clay"    )])[,1], 
                    gravel   = stack(soil[,c("t_gravel"  , "s_gravel"  )])[,1], 
                    oc       = stack(soil[,c("t_oc"      , "s_oc"      )])[,1], 
                    ph       = stack(soil[,c("t_ph"      , "s_ph"      )])[,1], 
                    ref_bulk = stack(soil[,c("t_ref_bulk", "s_ref_bulk")])[,1], 
                    sand     = stack(soil[,c("t_sand"    , "s_sand"    )])[,1], 
                    silt     = stack(soil[,c("t_silt"    , "s_silt"    )])[,1],
                    depth    = soil$soil_depth)
soil2$layer <- as.factor(ifelse(soil2$layer=="t", 1, 2))
levels(soil2$layer) <- c("topsoil", "subsoil")
summary(soil2)

# -----------------------
# 1) Save all the figures to 1 pdf
# -----------------------
# pdf(file.path(soil.out, "Soil_Prperties.pdf"))
# Soil Depth
png(file.path(soil.out, "soil_depth.png"), height=800, width=800)
print( ggplot(data= soil2[soil2$layer=="topsoil",]) + 
    geom_raster(aes(x=lon, y=lat, fill=depth)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="cm") +
    ggtitle("Soil Depth") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(panel.background=element_blank()) )
dev.off()

# OC
png(file.path(soil.out, "soil_oc.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=oc)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="% Weight") +
    ggtitle("Organic Carbon") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# CEC
png(file.path(soil.out, "soil_cec.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=cec)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="meq/100g") +
    ggtitle("Cation Exchange Capacity") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# pH
png(file.path(soil.out, "soil_ph.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=ph)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("red4", "red2", "gray50", "blue2"), name="-log(H+)", values=c(0,3.5/max(soil2$ph),7/max(soil2$ph),1)) +
    ggtitle("pH") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# ref_bulk
png(file.path(soil.out, "soil_ref_bulk.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=ref_bulk)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="g/cm3") +
    ggtitle("Reference Bulk Density") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# Sand
png(file.path(soil.out, "soil_sand.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=sand)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="% weight") +
    ggtitle("Sand Fraction") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# Silt
png(file.path(soil.out, "soil_silt.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=silt)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="% weight") +
    ggtitle("Silt Fraction") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# Clay
png(file.path(soil.out, "soil_clay.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=clay)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="% weight") +
    ggtitle("Clay Fraction") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()

# Gravel
png(file.path(soil.out, "soil_gravel.png"), height=800, width=800)
print( ggplot(data= soil2[,]) + facet_grid(layer~.) +
    geom_raster(aes(x=lon, y=lat, fill=gravel)) +
    geom_path(data=paleon.states, aes(x=long, y=lat, group=group)) +
    scale_x_continuous(limits=range(soil$lon), expand=c(0,0), name="Longitude (degrees)") +
    scale_y_continuous(limits=range(soil$lat), expand=c(0,0), name="Latitude (degrees)") +
    scale_fill_gradientn(colours=c("gray80", "blue3"), name="% volume") +
    ggtitle("Gravel Content") +
    coord_equal(ratio=1) +
	theme(legend.position="bottom") +
    theme(axis.text.x =element_text(color="black", size=rel(2)),
          axis.text.y =element_text(color="black", size=rel(2)), 
          axis.title.x=element_text(size=rel(2)),  
          axis.title.y=element_text(size=rel(2)),
          plot.title  =element_text(size=rel(2))) +
    theme(strip.text  = element_text(size=rel(1.5))) +
    theme(panel.background=element_blank()) )
dev.off()
# dev.off()
# -----------------------
# ----------------------------------------------
