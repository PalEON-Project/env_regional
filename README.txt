Information about PalEON MIP Phase 2 (Regional) Environmental Drivers
Christy Rollinson, crollinson@gmail.com
10 October 2015

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
   — File Path/Name: co2/paleon_*_co2.nc
   — File Format: netcdf, dim=[time]
   — Units: ppm
   — Temporal Extent: 850-2010
   - Temporal Resolution: annual, monthly
   — Processing Script: 2_seasonal_co2.R
   — File Description: This file provides the CO2 concentration for the entire PalEON domain.  The 
     entire domain uses the same CO2 concentration.  Two options exist: 1) 1 value per year (original 
     splicing by Bjorn Brooks in 2013; 2) monthly CO2 with variation added in from MsTMIP.
   — Citation: Original CO2: PMIP-3, NOAA Mona Loa Observatory (MLO)

               MsTMIP Seasonal Variability
               Wei, Y., S. Liu, D.N. Huntzinger, A.M. Michalak, N. Viovy, W.M. Post, C.R. Schwalm, K. 
               Schaefer, A.R. Jacobson, C. Lu, H. Tian, D.M. Ricciuto, R.B. Cook, J. Mao, and X. Shi. 
               2014. NACP MsTMIP: Global and North American Driver Data for Multi-Model Intercomparison. 
               Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory 
               Distributed Active Archive Center, Oak Ridge, Tennessee, USA. 
               http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Web Link: MsTMIP: http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Date Accessed: MsTMIP: 30 September, 2015
   — Additional Notes: We did not redo original script by Bjorn Brooks that spliced PMIP-3 and NOAA
     MLO data nor replace PMIP-3 with CMIP CO2.  

3) Land-Use
   — File Path/Name:lulcc/paleon_lulcc_*.nc
   — File Format: netcdf, dim=[time,lat,lon]
   — Units: variable (see “Additional Notes” below)
   — Temporal Extent: 0850-2010 (PalEON filled, native extent = 1500-2005)
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
     land.  LULCC was filled to the PalEON temporal domain according to the following scheme:
     1) 0850 - 1499: All land in primary vegetation cover,  no disturbances or land use transitions
        This is essentially assuming no significant human influence prior to 1500.  For models relying
        on the land cover portion, this may lead to abrupt transition to crops in the Ohio River Valley
        and could cause issues since we did not calculate the area moving from primary forest to cropland.
        If this is an issue, please let us know and we’ll fix it.
     2) 2005 - 2010: No more land cover transitions or harvest; land use designation stays as it was in 
        2005 (secondary land will keep aging).
                      
     Raw files are provided in a format similar to the original.  Utility scripts will be generated to help
     reformat drivers for particular models. 

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
   — File Path/Name: biome/biome_potential_vegtype_biome.nc; biome/biome_potential_vegtype_pft_fraction.nc
   — File Format: netcdf; dim=[PFT, lat, lon]
   — Units: categorical
   — Temporal Extent: static (one time)
   - Temporal Resolution: static (one time)
   — Processing Script: 5_biome.R (see 5_biome_comparison.R for some ugly script deciding whether to use
     SYNMAP (MsTMIP) or Ramankutty & Foley (PalEON original))
   — File Description: This file contains the dominant biome type (biome_potential_vegtype_biome.nc) and 
     estimated corresponding fraction of common model plant functional types (*_pft_fraction.nc) in each 
     grid cell. Original biome classification is from SYNMAP (MsTMIP biome driver). see additional notes
     for details on biome & PFT translation for PalEON.
   — Citation: Wei, Y., S. Liu, D.N. Huntzinger, A.M. Michalak, N. Viovy, W.M. Post, C.R. Schwalm, K. 
               Schaefer, A.R. Jacobson, C. Lu, H. Tian, D.M. Ricciuto, R.B. Cook, J. Mao, and X. Shi. 
               2014. NACP MsTMIP: Global and North American Driver Data for Multi-Model Intercomparison. 
               Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory 
               Distributed Active Archive Center, Oak Ridge, Tennessee, USA. 
               http://dx.doi.org/10.3334/ORNLDAAC/1220

               SYNMAP
               Jung, M., Henkel, K., Herold, M., and Churkina, G. 2006. Exploiting synergies of global land 
               cover products for carbon cycle modeling. Remote Sens. Environ. 101: 534-553. 
               DOI:10.1016/j.rse.2006.01.020

               Biome Crosswalk
               Poulter, B. N. MacBean. A. Hartley, I Khlystova, O. Arino, R. Betts, S. Bontemps, M. Boettcher,
               C. Brockmann, P. Defourny, S. Hagemann, M. Herold, G. Kirches, C. Lamarche, D. Lederer, C. 
               Ottle, M. Peters, and P. Peylin. 2015. Plant functional type classification for earth system 
               models: results from the European Space Agency’s Land Cover Climate Change Initiative.  
               Geoscientific Model Development 8: 2315-2328. doi:10.5194/gmd-8-2315-2015
   — Web Link: http://dx.doi.org/10.3334/ORNLDAAC/1220
   — Date Accessed: 28 September, 2015
   — Additional Notes: PalEON biome driver is derived from SYNMAP (Jung et al. 2006) which is used by 
     MsTMIP.  Crosswalk between SYNMAP and PalEON/model biomes generally followed previous work by K. Schaefer
     for SibCASA.  SYNMAP-PFT crosswalk was modified from Poulter et al. 2015.  In both cases, crops were 
     equated with grassland and urban area ignored. Mosaic tree-shrub or tree-grass/crop were treated as the 
     corresponding LCCS “open” forest types (<15% forest cover) treating shrub and grass as equivalent.  
     Current product keeps the bare ground cover, and models that do not have bare ground should reallocate 
     that proportion accordingly.

     Biome Classifications
     Code    Description                                          SYNMAP CODES   
      1         Broadleaf-Deciduous Forest                        5
      2         Mixed Broadleaf-Deciduous & Needleleaf-Evergreen  6, 9
      3         Needleleaf-Evergreen Forest                       1, 19
      4         Savanna & Shrublands                              10, 14, 23, 27, 37, 38, 40
      5         Grasslands                                        31, 32, 33, 36, 39, 41, 42, 44, 47

     PFT Classifications
     PFT   Code    Form     Leaf       Phenology
      1    TBrEv   Tree     Broad      Evergreen
      2    TBrDe   Tree     Broad      Deciduous
      3    TNeEv   Tree     Needle     Evergreen
      4    TNeDe   Tree     Needle     Deciduous
      5    SBrEv   Shrub    Broad      Evergreen
      6    SBrDe   Shrub    Broad      Deciduous
      7    SNeEv   Shrub    Needle     Evergreen
      8    SNeDe   Shrub    Needle     Deciduous
      9    Grass   Grass      -            - 
     10    Bare    Bare       -            - 



6) Nitrogen Deposition
   — File Path/Name: nitrogen/paleon_nhx.nc; nitrogen/paleon_noy.nc
   — File Format: netcdf, dim=[time,lat,lon]
   — Units: mgN/m2/yr
   — Temporal Extent: 0850-2010 (PalEON filled, native extent = 1860-2010)
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
     
     Temporal Gapfilling Protocol: Assume 1860 is indicative of background N deposition rates.  850 - 
     1859 receive the 1860 N deposition rates.
