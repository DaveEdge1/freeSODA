# cside2lipd

This script provides the starting point for transforming data from a netcdf file to a tileserver. The tiles are not currently overlaid at the correct position on the map. This may be related to the 'origins'

1. Grab data - (https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/files?subdataset=cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_202311)
2. Create a geotiff (netcdf2geotiff.R)
3. Use gdal2tiles in OSGeo4W shell to create tiles (gdal2tiles -z2-5 -sEPSG:3031 -praster --xyz -x C:\Users\dce25\Documents\RProjects\freeSODA\jan_ice_concentration3031.tif C:\Users\dce25\Documents\RProjects\freeSODA\tiles3031)
4. Push tiles to presto server (scp -r tiles3031 root@143.198.98.66:/root/presto/query/public/sea_ice_conc_jan)
