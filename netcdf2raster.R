library(ncdf4)
oceanDat <- ncdf4::nc_open("C:/Users/dce25/Downloads/cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_1740090063083.nc")
oceanDat$natts
oceanDat$ndims

sea_ice_thickness <- ncdf4::ncvar_get(oceanDat, "sithick")
dim(sea_ice_thickness)
sea_floor_potential_temperature <- ncdf4::ncvar_get(oceanDat, "bottomT")
dim(sea_floor_potential_temperature)
mixed_layer_thickness <- ncdf4::ncvar_get(oceanDat, "mlotst")
dim(mixed_layer_thickness)
ice_concentration <- ncdf4::ncvar_get(oceanDat, "siconc")
dim(ice_concentration)
#salinity <- ncdf4::ncvar_get(oceanDat, "so")
#dim(salinity)
#oceanDat$var$so$dim[[3]]$units
#oceanDat$var$so$dim[[3]]$vals
#temperature <- ncdf4::ncvar_get(oceanDat, "thetao")
#dim(temperature)
#oceanDat$var$so$dim[[3]]$units
#oceanDat$var$so$dim[[3]]$vals

latitude_deg_n <- oceanDat$dim$latitude$vals
longitude_deg_e <- oceanDat$dim$longitude$vals

# df1=data.frame(lon=longitude_deg_e,
#                lat=latitude_deg_n,
#                value=ice_concentration[,1]
#                  )
library(reshape2)
library(raster)
#rownames(ice_concentration) <- make.names(longitude_deg_e)
#colnames(ice_concentration) <- make.names(abs(latitude_deg_n))
ice_concentration <- cbind.data.frame(ice_concentration[,,1], "lon"=longitude_deg_e)
#ice_concentration[1:10,1:10,1]
ice_concentration2 <- melt(ice_concentration, id.vars = "lon")
head(ice_concentration2)
levels(ice_concentration2$variable) <- latitude_deg_n
head(ice_concentration2)
colnames(ice_concentration2)[2] <- "lat"

library(ggplot2)

ggplot(data = ice_concentration2, mapping = aes(x=lat, y=lon, color=value)) + geom_tile()

r = rasterFromXYZ(ice_concentration2)
plot(r)

projection(r)="+init=epsg:4326"



unique(ice_concentration2$variable)
