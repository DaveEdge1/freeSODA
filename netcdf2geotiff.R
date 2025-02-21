library(ncdf4)
library(reshape2)
library(raster)
library(ggplot2)

#Grab data
#https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/files?subdataset=cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_202311
oceanDat <- ncdf4::nc_open("C:/Users/dce25/Downloads/mercatorglorys12v1_gl12_mean_1993_2016_01.nc")
oceanDat$natts
oceanDat$ndims

#grab lat and lon, subset lat to below 30-degrees S
latitude_deg_n <- oceanDat$dim$latitude$vals
longitude_deg_e <- oceanDat$dim$longitude$vals

#subset to 1-degree res
lat_indices <- which(latitude_deg_n %% 1 == 0)[which(latitude_deg_n %% 1 == 0) <= max(which(latitude_deg_n<=-30))]
lon_indices <- which(longitude_deg_e %% 1 == 0)
latitude_deg_n <- latitude_deg_n[lat_indices]
longitude_deg_e<- longitude_deg_e[lon_indices]

#grab variables, subsetting to southern ocean at 1-degree res
sea_ice_thickness <- ncdf4::ncvar_get(oceanDat, "sithick")[lon_indices,lat_indices]
dim(sea_ice_thickness)
sea_floor_potential_temperature <- ncdf4::ncvar_get(oceanDat, "bottomT")[lon_indices,lat_indices]
dim(sea_floor_potential_temperature)
mixed_layer_thickness <- ncdf4::ncvar_get(oceanDat, "mlotst")[lon_indices,lat_indices]
dim(mixed_layer_thickness)
ice_concentration <- ncdf4::ncvar_get(oceanDat, "siconc")[lon_indices,lat_indices]
dim(ice_concentration)

#rownames(ice_concentration) <- make.names(longitude_deg_e)
#colnames(ice_concentration) <- make.names(abs(latitude_deg_n))
ice_concentration <- cbind.data.frame(ice_concentration, "lon"=longitude_deg_e)
#ice_concentration[1:10,1:10,1]
ice_concentration2 <- melt(ice_concentration, id.vars = "lon")
head(ice_concentration2)
levels(ice_concentration2$variable) <- latitude_deg_n
head(ice_concentration2)
colnames(ice_concentration2)[2] <- "lat"

#ice_concentration3 <- ice_concentration2[seq(1,length(ice_concentration2$lon),by=12),]

#p1 <- ggplot(data = ice_concentration2, mapping = aes(x=lat, y=lon, color=value)) + geom_tile()
#ggsave('jan_ice_concentration.png', plot = p1)

r = rasterFromXYZ(ice_concentration2, crs = ('+proj=longlat'))
#crs(r) <- "epsg:3857"
#r_ext <- extent(-180,180,-90,90)
#r2 <- crop(r,r_ext)
plot(r)
r_polar = raster::projectRaster(r, crs = CRS('epsg:3031'))
plot(r_polar)
raster::writeRaster(r_polar, "jan_ice_concentration.tiff")

terra::writeRaster(r_polar, "Data/jan_ice_concentration.tif", filetype = "GTiff", overwrite = TRUE)
