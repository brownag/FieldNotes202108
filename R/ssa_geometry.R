# get geometry for all SSAs
library(sf)
library(soilDB)

areasym <- SDA_query("SELECT DISTINCT areasymbol FROM legend")

ssafeats <- fetchSDA_spatial(areasym$areasymbol,
                             by.col = "areasymbol", 
                             geom.src = "sapolygon", 
                             chunk.size = 1)

sf::st_write(sf::st_as_sf(ssafeats), dsn="E:/Geodata/all_ssas.shp")

