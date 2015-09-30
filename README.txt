Information about PalEON MIP Phase 2 (Regional) Environmental Drivers
Christy Rollinson, crollinson@gmail.com
17 September 2015

This directory contains the environmental drivers for phase 2 MIP runs (regional runs, 
formerly phase 1b) & any scripts used to generate them.  Note for BU paleon users: all 
raw data sets (pre-processing) can be found in the env_drivers_raw folder in this 
repository. All other users can obtain these data via links provided or drivers and raw
data will be made available upon request.

There are 6 Environmental Drivers.  File paths associated each driver are nested within
the current version of env drivers.
Current Env Driver: phase2_env_drivers_v1

All drivers have a consistent spatial extent, resolution, and coordinate reference system
Spatial Extent (xmin, xmax, ymin, ymax) = -100, -60, 35, 50
Spatial Resolution = 0.5 x 0.5 degree
CRS (R format) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

See specific driver information for temporal extent and resolution

1) Land Mask
   — File Path/Name: domain_mask/paleon_domain.nc
   — File Format: netcdf, dim=[lat,lon]
   — Units: binary (presence/absence)
   — Temporal Extent: static (one time)
   - Temporal Resolution: static (one time)
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
   — Units: variable (see “Additional Notes” below)
   — Temporal Extent: 
   - Temporal Resolution: 
   — Processing Script:
   — File Description:
   — Citation:
   — Web Link:
   — Date Accessed: 
   — Additional Notes: 

3) Land-Use
   — File Path/Name:lulcc/paleon_lulcc_*.nc
   — File Format: netcdf, dim=[time,lat,lon]
   — Units: variable (see “Additional Notes” below)
   — Temporal Extent: 1500-2005
   - Temporal Resolution: annual
   — Processing Script: 3_lulcc.R
   — File Description:Spatial time series describing fractional land use types, fraction land use 
     transition, and harvest statistics (area, biomass).  These are modified Hurtt and differ greatly
     from MsTMIP drivers due to land use classifications.
   — Citation: Chini, L.P., G.C. Hurtt, and S. Frolking. 2014. Harmonized Global Land Use for Years 
               1500-2100, V1. Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge 
               National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA. 
               http://dx.doi.org/10.3334/ORNLDAAC/1248
   — Web Link: http://dx.doi.org/10.3334/ORNLDAAC/1248
   — Date Accessed: 28 September, 2015
   — Additional Notes: PalEON drivers are from LUHa.v1 files only, which are 1500-2005, with no urban 
     land. Raw files are provided in a format similar to the original.  Utility scripts will be generated
     to help reformat drivers for particular models. 

     Variable codes are as follows:
     Variable	Description
     Land Cover- Units=fraction (%) of each grid cell
       gcrop   cropland
       gothr   primary land
       gpast   pasture
       gsecd   secondary land
       gurbn   urban land
     Secondary land characteristics                     Units
       gssma   mean age of secondary land               years
       gssmb   mean biomass density of secondary land   kgC/m2
     Wood harvest data	                                                                                Units
       gfsh1   fraction of each grid cell that had wood harvested from mature secondary forested land   fraction (%) of each grid cell
       gfsh2   fraction of each grid cell that had wood harvested from mature secondary forested land   fraction (%) of each grid cell
       gfsh3   fraction of each grid cell that had wood harvested from secondary non-forested land      fraction (%) of each grid cell
       gfvh1   fraction of each grid cell that had wood harvested from primary forested land            fraction (%) of each grid cell
       gfvh2   fraction of each grid cell that had wood harvested from primary non-forested land        fraction (%) of each grid cell
       gsbh1   mature secondary forest biomass harvested                                                kgC
       gsbh2   young secondary forest biomass harvested                                                 kgC
       gsbh3   secondary non-forest biomass harvested                                                   kgC
       gvbh1   primaryforest biomass harvested                                                          kgC
       gvbh2   primary non-forest biomass harvested                                                     kgC
     Secondary land characteristics- Units=fraction (%) of each grid cell
       gflcp   transition from cropland to pasture
       gflcs   transition from cropland to secondary land
       gflcu   transition from cropland to urban land
       gflpc   transition from pasture to cropland
       gflps   transition from pasture to secondary land
       gflpu   transition from pasture to urban land
       gflsc   transition from secondary land to cropland
       gflsp   transition from secondary land to pasture
       gflsu   transition from secondary land to urban land
       gfluc   transition from urban land to crop
       gflup   transition from urban land to pasture
       gflus   transition from urban land to secondary land
       gflvc   transition from primary land to cropland
       gflvp   transition from primary land to pasture
       gflvu   transition from primary land to urban land

4) Soil
   — File Path/Name: soil/paleon_soil.nc (all variables), soil/paleon_soil_*.nc (individual files)
   — File Format: netcdf dim=[lat,lon]
   — Units: variable (see “Additional Notes” below)
   — Temporal Extent: static (one time)
   - Temporal Resolution: static (one time)
   — Processing Script: 4_soil.R
   — File Description: There are two file options for soil drivers:  1) single .nc file that mirror
     the format of the MsTMIP driver; 2) individual .nc files for each driver.  The data in these
     formats are identical.
   — Citation: Wei, Y., S. Liu, D.N. Huntzinger, A.M. Michalak, N. Viovy, W.M. Post, C.R. Schwalm, K. 
               Schaefer, A.R. Jacobson, C. Lu, H. Tian, D.M. Ricciuto, R.B. Cook, J. Mao, and X. Shi. 
               2014. NACP MsTMIP: Global and North American Driver Data for Multi-Model Intercomparison. 
               Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory 
               Distributed Active Archive Center, Oak Ridge, Tennessee, USA. 
               http://dx.doi.org/10.3334/ORNLDAAC/1220

   — Web Link: http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Date Accessed: 29 September, 2015
   — Additional Notes: PalEON soil drivers were aggregated up from MsTMIP drivers for North America. In
     the process of aggregating from quarter to half-degree resolution, we had a few cells that were 
     missing in MsTMIP and which were filled with the mean of the surrounding adjacent cells.
     MsTMIP data source: Unified North American Soil Database (UNASD) [U.S. General Soil Map (STATSGO2) 
     + Soil Landscapes of Canada v3.2 and v2.2 + HWSD v1.1]; topsoil=0-30 cm depth, subsoil=30-100 cm 
     depth; See highlighted section on page 18 of NACP_MsTMIP_Model_Driver.pdf in soil folder for more 
     information.

     Units & descriptions are as follows:
     Soil Property    Description                            Units in UNASD
      soil_code        soil mapping unit code                 code
      ref_depth        reference soil depth                   code
      roots            obstacles to roots (Europe only)       NA
      il               impermeable layer (Europe only)        NA
      t_cec_clay       topsoil CEC (clay)                     meq/100g
      t_clay           topsoil clay fraction                  % weight
      t_gravel         topsoil gravel content                 % volume
      t_oc             topsoil organic carbon                 % weight
      t_ph_h20         topsoil pH (H2O)                       -log(H+)
      t_ref_bulk       topsoil bulk density                   g/cm3
      t_sand           topsoil sand fraction                  % weight
      t_silt           topsoil silt fraction                  % weight
      t_usda_tex       topsoil USDA texture classification    name
      s_cec_clay       subsoil CEC (clay)                     meq/100g
      s_clay           subsoil clay fraction                  % weight
      s_gravel         subsoil gravel content                 % volume
      s_oc             subsoil organic carbon                 % weight
      s_ph_h20         subsoil pH (H2O)                       -log(H+)
      s_ref_bulk       subsoil bulk density                   g/cm3
      s_sand           subsoil sand fraction                  % weight
      s_silt           subsoil silt fraction                  % weight
      s_usda_tex       subsoil USDA texture classification    name



5) Biome
   — File Path/Name: biome/biome_potential_vegtype_dominant.nc; biome/biome_potential_vegtype_relative.nc
   — File Format: netcdf; dim=[PFT, lat, lon]
   — Units: categorical
   — Temporal Extent: static (one time)
   - Temporal Resolution: static (one time)
   — Processing Script: 5_biome.R
   — File Description: This file contains the potential dominant vegetation type (*_dominant.nc) and 
     estimate fraction of each vegetation type (*_relative.nc) in a grid cell.  Fraction of each vegetation
     type in a grid cell was estimated using a 5 x 5 degree smoothing.
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
   — File Path/Name: nitrogen/paleon_nhx.nc; nitrogen/paleon_noy.nc
   — File Format: netcdf, dim=[time,lat,lon]
   — Units: mgN/m2/yr
   — Temporal Extent: 1860-2010
   - Temporal Resolution: annual
   — Processing Script: 6_nitrogen.R
   — File Description: Annual nitrogen inputs from spatially extracted from the MsTMIP drivers.
   — Citation: Wei, Y., S. Liu, D.N. Huntzinger, A.M. Michalak, N. Viovy, W.M. Post, C.R. Schwalm, K. 
               Schaefer, A.R. Jacobson, C. Lu, H. Tian, D.M. Ricciuto, R.B. Cook, J. Mao, and X. Shi. 
               2014. NACP MsTMIP: Global and North American Driver Data for Multi-Model Intercomparison. 
               Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory 
               Distributed Active Archive Center, Oak Ridge, Tennessee, USA. 
               http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Web Link: http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Date Accessed: 29 September, 2015
   — Additional Notes: See highlighted section on page 16 of NACP_MsTMIP_Model_Driver.pdf in nitrogen 
     folder for description of MsTMIP methods.  Because there is no methodological difference between
     quarter-degree North America data and half-degree global data, PalEON drivers were extracted from
     the 0.5-degree native global product.
