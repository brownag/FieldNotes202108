# sda level of horizon detail
# 
library(soilDB)
library(dplyr)

# number of horizons and components in ssurgo
q <- SDA_query("SELECT COUNT(DISTINCT chkey), COUNT(DISTINCT component.cokey) FROM chorizon 
                INNER JOIN component ON component.cokey = chorizon.cokey
                INNER JOIN mapunit ON mapunit.mukey = component.mukey
                INNER JOIN legend ON legend.lkey = mapunit.lkey
                WHERE areasymbol != 'US'
                GROUP BY areasymbol")

# averages about 4 layer per
mean(q$V1/q$V2)

q2 <- SDA_query("SELECT legend.areasymbol, legend.areaname, 
                        COUNT(DISTINCT chkey), COUNT(DISTINCT component.cokey),
                        distlegendmd.cordate FROM chorizon 
                INNER JOIN component ON component.cokey = chorizon.cokey
                INNER JOIN mapunit ON mapunit.mukey = component.mukey
                INNER JOIN legend ON legend.lkey = mapunit.lkey
                INNER JOIN distlegendmd ON legend.areasymbol = distlegendmd.areasymbol
                WHERE legend.areasymbol != 'US'
                GROUP BY legend.areasymbol, legend.areaname, distlegendmd.cordate")

res <- q2 %>% 
  dplyr::group_by(areasymbol) %>% 
  dplyr::filter(!is.na(cordate)) %>% 
  arrange(desc(cordate)) %>% 
  dplyr::slice(1) 

na.cordate <- q2 %>% 
  dplyr::filter(is.na(cordate)) %>% 
  select(areasymbol) %>%  
  unique() %>% `[[`(1) %>% sort()

res[which.max(res$V3 / res$V4),]

res$year <- as.numeric(strftime(as.Date(res$cordate, "%m/%d/%Y %H:%M:%S"), "%Y"))
# q3 <- soilDB:::.SDA_query_FOR_JSON_AUTO(
#   "SELECT areasymbol, areaname, mapunit.mukey, 
#                     COUNT(DISTINCT chkey) AS count_chkey, 
#                     COUNT(DISTINCT component.cokey) AS count_cokey 
#                   FROM chorizon 
#                 INNER JOIN component ON component.cokey = chorizon.cokey
#                 INNER JOIN mapunit ON mapunit.mukey = component.mukey
#                 INNER JOIN legend ON legend.lkey = mapunit.lkey
#                 WHERE areasymbol != 'US'
#                 GROUP BY mapunit.mukey")

yrclasses <- cut(res$year, seq(1939, 2019, 10))
levels(yrclasses) <- paste0(seq(1940, 2019, 10), "s")
plot((res$V3 / res$V4) ~ res$year, 
     xlab = "Correlation Date (distlegendmd.cordate)", 
     ylab = "Average # of Layers per Component (whole SSA)")

plot((res$V3 / res$V4) ~ yrclasses, 
     xlab = "Correlation Date (distlegendmd.cordate)", 
     ylab = "Average # of Layers per Component (whole SSA)")

summary(lm((res$V3 / res$V4) ~ factor(yrclasses)))

library(sf)

ssas <- st_read("E:/Geodata/all_ssas.shp")

st_crs(ssas) <- st_crs(4326)
conus <- state.abb[!state.abb %in% c("HI","AK")]

usssas <- vapply(conus, function(x) stringi::stri_match_all(ssas$areasymbol, regex = x),
                 FUN.VALUE = vector('list', nrow(ssas))) %>% data.frame()
idx <- do.call('c', apply(usssas, 2, function(x) which(!is.na(x))))
ssas2 <- merge(ssas, res, by = "areasymbol", all.x = TRUE)
ssas2 <- st_transform(ssas2, 6350)
ssas2$ratio <- ssas2$V3 / ssas2$V4
plot(ssas2[idx, "ratio"], main = 'Ratio of Component Layers to Components')

plot(subset(ssas2, areasymbol %in% na.cordate))
# #raster::plot(fasterize::fasterize(ssas2[idx,], raster::raster(res=1000), field = "ratio"))
# coord_sf(crs = st_crs(2163), xlim = c(-2500000, 2500000), ylim = c(-2300000, 730000))

# ldm snapshot level of detail
ldm <- DBI::dbConnect(RSQLite::SQLite(), "E:/Geodata/soils/LDM-compact.sqlite")
DBI::dbListTables(ldm)
ldm_hz <- DBI::dbReadTable(ldm, "layer")
length(unique(ldm_hz$layer_key)) / length(unique(ldm_hz$pedon_key))

# pedon snapshot level of detail
mrp <- DBI::dbConnect(RSQLite::SQLite(), "C:/Geodata/soils/NASIS-data_patched.sqlite")
DBI::dbListTables(mrp)
mrp_hz <- DBI::dbReadTable(mrp, "phorizon")
length(unique(mrp_hz$phiid))
length(unique(mrp_hz$peiidref))
prop.table(table(!is.na(mrp_hz$bounddistinct)))
prop.table(table(!is.na(mrp_hz$boundtopo)))
