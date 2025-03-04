library(ncdf4)
library(reshape2)
library(terra)
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
lat_indices <- which(latitude_deg_n %% 1 == 0)[which(latitude_deg_n %% 1 == 0) <= max(which(latitude_deg_n<=-50))]
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
summary(ice_concentration2$lat)
ice_concentration2$lat <- as.numeric(as.character(ice_concentration2$lat))

#ice_concentration3 <- ice_concentration2[seq(1,length(ice_concentration2$lon),by=12),]

#p1 <- ggplot(data = ice_concentration2, mapping = aes(x=lat, y=lon, color=value)) + geom_tile()
#ggsave('jan_ice_concentration.png', plot = p1)

saveRDS(ice_concentration2, "jan_ice_concentration.RDS")

#ice_concentration2[is.na(ice_concentration2)] <- 0
r = raster::rasterFromXYZ(ice_concentration2,)
# crs(r)
# crs(r) <- crs("4326")
# crs(r)
# #r_ext <- extent(-180,180,-90,90)
# #r2 <- crop(r,r_ext)
# summary(r)
# plot(r)
# extent(r)
# r*254
# #r_polar = raster::projectRaster(r, crs = CRS('epsg:3031'))
# #plot(r_polar)
# #raster::NAvalue(r_polar) <- 9999
# raster::writeRaster(r*254,'jan_ice_concentration_latlon.tif',format='GTiff',overwrite=TRUE, datatype='INT1U')


# crsAntartica <-  leafletCRS(
#   crsClass = 'L.Proj.CRS',
#   code = 'EPSG:3031',
#   proj4def = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#   resolutions = c(48310.14147851562, 24155.07073925781, 12077.535369628906, 6038.767684814453, 3019.3838424072264, 1509.6919212036132, 754.8459606018066, 377.4229803009033, 188.71149015045165, 94.35574507522583, 47.17787253761291, 23.588936268806457, 11.794468134403228, 5.897234067201614, 2.948617033600807, 1.4743085168004035, 0.7371542584002018),
#   origin = c(-12367396.2185, 12367396.2185)#,
#   #bounds =  list( c(-4194304, -4194304), c(4194304, 4194304) )
# )

r2 <- terra::rast(r)
#CRS(r2, "epsg:4326")
plot(r2)
summary(r2)
#r2 <- r2*255
#r2 <- terra::subst(r2, NA, 0)
#plot(r2)

# terra::writeRaster(r2,
#                    "Data/jan_ice_concentration.tif",
#                    filetype = "GTiff",
#                    overwrite = TRUE)


summary(r2)
crs(r2) <- "epsg:4326"
r_polar2 <- terra::project(r2, y = "EPSG:3031")
#terra::origin(r_polar2) <- c(-8929486, 6564912)
#r_polar2 <- terra::subst(r_polar2, NA, 0)
plot(r_polar2)
r_polar2 <- r_polar2*254
r2 <- terra::subst(r2, NA, 254)
#class(r_polar2)
summary(r_polar2)
print(r_polar2)
terra::writeRaster(r_polar2,
                   "jan_ice_concentration3031.tif",
                   filetype = "GTIFF",
                   overwrite = TRUE,
                   datatype='INT1U')


# terra::describe("jan_ice_concentration3031.tif")
# grep("NoData Value", terra::describe("jan_ice_concentration3031.tif"), value=TRUE)

library(tiler)
library(raster)
tiler_options(osgeo4w="C:\\OSGeo4W")
#tile_dir <- file.path(tempdir(), "tiles")
map <- system.file("maps/map_wgs84.tif", package = "tiler")
(r <- raster(map))
plot(r)
tile(map, "C:/Users/dce25/Documents/R Projects/freeSODA/tiles", "0-3")
list.files("C:/Users/dce25/Documents/R Projects/freeSODA/tiles")
tiler_options()

