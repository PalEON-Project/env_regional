Information about PalEON MIP Phase 2 (Regional) Environmental Drivers
Christy Rollinson, crollinson@gmail.com
17 September 2015

This directory contains the environmental drivers for phase 2 MIP runs (regional runs, 
formerly phase 1b) & any scripts used to generate them.  Note for BU paleon users: all 
raw data sets (pre-processing) can be found in the env_drivers_raw folder in this 
repository. All other users can obtain these data via links provided or drivers and raw
data will be made available upon request.

There are 6 Environmental Drivers.  
1) Land Mask
   — File Path/Name: domain_mask/paleon_domain.nc
   — File Format: netcdf
   — Processing Scripts: 1_paleon_mask.R
   — File Description: This is a spatial file that all data (drivers & outputs) should 
     match spatial.  Resolution: 0.5-degree, bounding box (xmin, xmax, ymin, ymax): -100, -60, 35, 50
   — Source: PalEON met drivers phase2 v1, precipf_0850_01_01.nc (first time step)
   — Additional Notes: The base mask for the met drivers was originally created from 
     precipf_0850_01.nc and all met files were masked using this layer during met
     processing (5_rewrite_timestamps.R, https://github.com/PalEON-Project/met_regional).
     Because this met originated with CCSM4 data, which was coarser than the 0.5-degree 
     resolution we are currently working with, a decent amount of coastline is missing
     in the east and along the upper peninsula of Michigan.

2) CO2
   — File Path/Name:
   — File Format:
   — Processing Script:
   — File Description:
   — Citation:
   — Web Link:
   — Date Accessed: 
   — Additional Notes: 

3) Land-Use
   — File Path/Name:
   — File Format:
   — Processing Script:
   — File Description:
   — Citation:
   — Web Link:
   — Date Accessed: 
   — Additional Notes: 

4) Soil
   — File Path/Name:
   — File Format:
   — Processing Script:
   — File Description:
   — Citation:
   — Web Link:
   — Date Accessed: 
   — Additional Notes: 

5) Biome
   — File Path/Name:
   — File Format:
   — Processing Script:
   — File Description:
   — Citation: Ramankutty, N. and J.A. Foley. 2010. ISLSCP II Potential Natural Vegetation Cover. 
               In Hall, Forest G., G. Collatz, B. Meeson, S. Los, E. Brown de Colstoun, and D. 
               Landis (eds.). ISLSCP Initiative II Collection. Data set. Available on-line 
               [http://daac.ornl.gov/] from Oak Ridge National Laboratory Distributed Active Archive 
               Center, Oak Ridge, Tennessee, U.S.A. doi:10.3334/ORNLDAAC/961
   — Web Link: http://dx.doi.org/10.3334/ORNLDAAC/961 (http://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=961)
   — Date Accessed: 17 September, 2015
   — Additional Notes: We are using the half-degree “potential_veg” product (not potential_veg_diffs)
     Codes & Biome Types are as follows: (taken from 0_potential_veg_readme.txt)
       Class #     Biome Type
       -------     ----------
          0        Water Bodies (Water bodies OR missing data for
                   "potential_veg_diffs_XX.asc" files)
          1        Tropical Evergreen Forest/Woodland
          2        Tropical Deciduous Forest/Woodland
          3        Temperate Broadleaf Evergreen Forest/Woodland
          4        Temperate Needleleaf Evergreen Forest/Woodland
          5        Temperate Deciduous Forest/Woodland
          6        Boreal Evergreen Forest/Woodland
          7        Boreal Deciduous Forest/Woodland
          8        Mixed Forest
          9        Savanna
         10        Grassland/Steppe
         11        Dense Shrubland
         12        Open Shrubland
         13        Tundra
         14        Desert
         15        Polar desert/Rock/Ice
         16        No Data over Land (not included in "potential_veg_diffs_XX.asc" 
                   files)

6) Nitrogen Deposition
   — File Path/Name:
   — File Format:
   — Processing Script:
   — File Description:
   — Citation:
   — Web Link:
   — Date Accessed: 
   — Additional Notes: 

More information on the drivers and their use can be found in the Phase 2 protocol.