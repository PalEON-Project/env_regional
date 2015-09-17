Information about PalEON MIP Phase 2 (Regional) Environmental Drivers
Christy Rollinson, crollinson@gmail.com
17 September 2015

This directory contains the environmental drivers for phase 2 MIP runs (regional runs, 
formerly phase 1b) & any scripts used to generate them

There are 6 Environmental Drivers.  
1) Land Mask
   — File Path/Name: domain_mask/paleon_domain.nc
   — File Format: netcdf
   — Source: met drivers, precipf_0850_01_01.nc (first time step)
   — Generation Scripts: 1_paleon_mask.R
   — File Description: This is a spatial file that all data (drivers & outputs) should 
     match spatial.  Resolution: 0.5-degree, bounding box (xmin, xmax, ymin, ymax): -100, -60, 35, 50
   — Additional Notes: The base mask for the met drivers was originally created from 
     precipf_0850_01.nc and all met files were masked using this layer during met
     processing (5_rewrite_timestamps.R, https://github.com/PalEON-Project/met_regional).
     Because this met originated with CCSM4 data, which was coarser than the 0.5-degree 
     resolution we are currently working with, a decent amount of coastline is missing
     in the east and along the upper peninsula of Michigan.

2) CO2
   — File Format:
   — Source:
   — Generation Scripts:
   — File Description:
   — Additional Notes:

3) Land-Use
   — Source:
   — Generation Scripts:
   — File Format:
   — File Description:
   — Additional Notes:

4) Soil
   — Source:
   — Generation Scripts:
   — File Format:
   — File Description:
   — Additional Notes:

5) Biome
   — Source:
   — Generation Scripts:
   — File Format:
   — File Description:
   — Additional Notes:

6) Nitrogen Deposition
   — Source:
   — Generation Scripts:
   — File Format:
   — File Description:
   — Additional Notes:

More information on the drivers and their use can be found in the Phase 2 protocol.