# ----------------------------------------------
# convert PalEON Phase 1 MIP annual CO2 concentrations to monthly seasonal files based on MsTMIP
# Original Script: Jaclyn Hatala Matthes
#                  17 March, 2014
# Updated  Script: Christine Rollinson, crollinson@gmail.com
#                  30 September, 2015
#
# --------------
# CO2 Correction Proceedure
# --------------
# 1) Download and crop MsTMIP CO2 driver
# 2) Calculate annual MsTMIP annual mean, 
#    create monthly adjustment value
# 3) use 1700 seasonal variability for 0850-1700
#    use calculated variability for 1700-2010
# --------------
#
# ----------------------------------------------

# ----------------------------------------------
# Load libaries, Set up Directories, etc
# ----------------------------------------------
library(ncdf4); library(raster); library(rgdal)
paleon.mask <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/domain_mask/paleon_domain.nc"

co2.bjorn  <- "~/Desktop/Research/PalEON_CR/env_regional/env_drivers_raw/co2/bjorn/paleon_co2_mix.nc"
co2.mstmip <- "~/Desktop/Research/PalEON_CR/env_regional/env_drivers_raw/co2/NACP_MSTMIP_MODEL_DRIVER/data/mstmip_driver_global_hd_co2_v1.nc4"
outpath    <- "~/Desktop/Research/PalEON_CR/env_regional/env_paleon/co2/"

# Create the driver folder if it doesn't already exist
if(!dir.exists(soil.out)) dir.create(soil.out)

# dummy fill values
fillv <- 1e30

# Load the paleon mask
paleon <- raster(paleon.mask)
paleon
# ----------------------------------------------

# ----------------------------------------------
pl.nc   <- nc_open(co2.bjorn)
pl.time <- ncvar_get(pl.nc,"time")
pl.co2  <- ncvar_get(pl.nc,"CO2air")
nc_close(pl.nc)

ms.nc   <- nc_open(co2.mstmip)
ms.time <- ncvar_get(ms.nc,"time") #monthly, starting in 01-1700
ms.lon  <- ncvar_get(ms.nc,"lon")
ms.lat  <- ncvar_get(ms.nc,"lat")
ms.co2  <- ncvar_get(ms.nc,"CO2")
ms.yr   <- 1700:2010
nc_close(ms.nc)
length(ms.yr)*12; length(ms.time)

#loop over mstmip time, crop to PalEON domain and get monthly variability
dom.lon <- which(ms.lon > xmin(paleon) & ms.lon < xmax(paleon))
dom.lat <- which(ms.lat > ymin(paleon) & ms.lat < ymax(paleon))

co2.mon.var <- vector()
for(t in 1:(length(ms.time)/12)){
  t.ind <- ((t-1)*12+1):(t*12)
  co2.avg <- apply(ms.co2[dom.lon,dom.lat,t.ind],c(1,2),mean) #annual mean
  for(m in 1:12){
    co2.var <- co2.avg - ms.co2[dom.lon,dom.lat,(min(t.ind)+m-1)] 
    co2.mon.var[min(t.ind)+m-1] <- mean(co2.var)
  }
}
plot(seq(1700,2011-1/12,by=1/12),co2.mon.var,main="MsTMIP CO2 monthly variability",xlab="time",ylab="CO2 variability [ppm]", type="l")

#add mstmip seasonal cycle to PalEON CO2 record
# Note: here is where we get rid of the extra first 850 years of bjorn's record
co2.mon.new <- vector()
for(y in 850:2010){
  for(m in 1:12){
    if(y<1700){
      co2.mon.new[(y-850)*12+m] <- pl.co2[y] + co2.mon.var[m]
    } else {
      co2.mon.new[(y-850)*12+m] <- pl.co2[y] + co2.mon.var[(y-1700)*12+m]
    }
  }
}
# plot(seq(850,2011-1/12,by=1/12),co2.mon.new,main="PalEON CO2 monthly variability",xlab="time",ylab="CO2 [ppm]", type="l")
# plot(co2.mon.new[1:(100*12)], type="l")
# plot(co2.mon.new[(length(co2.mon.new)-100*12):length(co2.mon.new)], type="l")
# length(co2.mon.new)

#format monthly netCDF file for output
# Specify time units for this year and month
nc_time_units <- paste('months since 0850-01-01 00:00:00', sep='')
nc.time       <- (seq(850,2011-1/12,by=1/12)-850)*12
time          <- ncdim_def("time",nc_time_units,nc.time,unlim=TRUE)

data <- co2.mon.new 
nc_variable_long_name <- 'Average monthly atmospheric CO2 concentration [ppm]'
nc_variable_units='ppm'

fillv   <- 1e+30
nc_var  <- ncvar_def('co2',nc_variable_units,
                    list(time), fillv, longname=nc_variable_long_name,prec="double")

ofname  <- paste(outpath,"paleon_monthly_co2.nc",sep="")
newfile <- nc_create( ofname, nc_var ) # Initialize file 

ncatt_put( newfile, nc_var, time, 'monthly')
ncatt_put( newfile, 0, 'description',"PalEON annual CO2 with MsTMIP seasonal CO2 variability imposed")
ncvar_put(newfile, nc_var, data) # Write netCDF file

nc_close(newfile)  

#format netCDF file for output
# Specify time units for this year and month
nc_time_units <- paste('years since 0850-01-01 00:00:00', sep='')
nc.time       <- 850:2010-850
time          <- ncdim_def("time",nc_time_units,nc.time,unlim=TRUE)

data <- pl.co2[850:2010] 
nc_variable_long_name <- 'Average annual atmospheric CO2 concentration [ppm]'
nc_variable_units='ppm'

fillv   <- 1e+30
nc_var  <- ncvar_def('co2',nc_variable_units,
                        list(time), fillv, longname=nc_variable_long_name,prec="double")

ofname  <- paste(outpath,"paleon_annual_co2.nc",sep="")
newfile <- nc_create( ofname, nc_var ) # Initialize file 

ncatt_put( newfile, nc_var, time, 'yearly')
ncatt_put( newfile, 0, 'description',"PalEON annual CO2 concentrations")
ncvar_put(newfile, nc_var, data) # Write netCDF file

nc_close(newfile) 
