library(aqp)
library(soilDB)

res <- fetchKSSL(mlra = c("18","22A"), returnMorphologicData = TRUE)
res <- res$SPC
sub <- aqp::subset(res, taxorder %in% c("Ultisols", "Alfisols"))
coordinates(sub) <- ~ x + y
proj4string(sub) <- "EPSG:4326"

bdy <- sf::st_as_sf(fetchSDA_spatial("CA649", by.col = "areasymbol", geom.src = "sapolygon"))
bdy <- sf::st_transform(bdy, sf::st_crs(32610))
bdy_buf <- sf::st_buffer(bdy, 3000)
pts <- sf::st_transform(sf::st_as_sf(as(sub, 'SpatialPointsDataFrame')), sf::st_crs(32610))

plot(bdy_buf$geometry)
plot(pts$geometry, add = TRUE)

y <- sf::st_intersection(pts, bdy_buf)
mysub <- subset(sub, profile_id(sub) %in% y$pedon_key)
aqp::plot(mysub, label = "pedon_id")

mysub[1,]$bs82
mysub[1,]$estimated_oc
mysub[1,]$pedon_id

mysub[4,]$bs82
mysub[4,]$estimated_oc
mysub[4,]$taxonname
plot(mysub[c(1,4),], color="estimated_oc")

# bound.lut <- c('V'=0.5,'A'=2,'C'=5,'G'=15,'D'=45) / 3
horizons(mysub)$boundsd <- 2.5 #bound.lut[toupper(substr(mysub$bounddistinct,0,1))]

c_stocks <-  function(x) {
  x %>%
    trunc(0, 100) %>%
    transform(thk = hzn_bot - hzn_top) %>%
    mutate_profile(
      c_stock_partial = (estimated_oc / 100) * ((db_13b * 1e6) / 1000) * (thk / 100),
      c_stock_pfrags = c_stock_partial * (1 - frags / 100),
      cml_stock = cumsum(ifelse(
        is.na(c_stock_partial), 0, c_stock_partial
      )),
      cml_stock_frags = cumsum(ifelse(is.na(c_stock_pfrags), 0, c_stock_pfrags)),
      c_stock = sum(c_stock_partial, na.rm = TRUE),
      c_stock_frags = sum(c_stock_pfrags, na.rm = TRUE)
    )
}

perturb_c_stocks <- function(x, n = 10000, boundary.attr) {
  perturb(x, n = n, boundary.attr = boundary.attr) %>% 
    c_stocks()
}

f <- fetchNASIS(SS=FALSE)
f.sub <- subset(f, pedon_id == "S09CA109002")
mysub_1 <- mysub[1,]
mysub_1$bounddistinct <- f.sub[2,]$bounddistinct[1:5]

# constant
mysub_1_sim <- perturb_c_stocks(mysub_1[1,], boundary.attr = "boundsd")
# varying
bound.lut <- c('V'=0.5,'A'=2,'C'=5,'G'=15,'D'=45) / 3
mysub_1$boundsd2 <- bound.lut[toupper(substr(mysub_1$bounddistinct,0,1))]
mysub_1$boundsd2[is.na(mysub_1$boundsd2)] <- 0
mysub_1_sim2 <- perturb_c_stocks(mysub_1[1,], boundary.attr = "boundsd2")

# constant SD
plot(density(mysub_1_sim$c_stock_frags, bw = 0.05), xlim = c(5, 10), ylim = c(0, 1.3),
     main = paste0("Simulated 0-100cm Organic Carbon Stocks, kg/m^2\n(", mysub[1,]$pedon_id, ")"),
     sub = "Assuming 2.5cm standard deviation in horizon boundaries")
abline(v = c_stocks(mysub[1,])$c_stock_frags, lty=2)
legend("topleft", c("Simulated","Source"), lty=c(1,2))

# boundary distinctness based SD
plot(density(mysub_1_sim2$c_stock_frags, bw = 0.05), xlim = c(5, 10), ylim = c(0, 1.3), 
     main = paste0("Simulated 0-100cm Organic Carbon Stocks, kg/m^2\n(", mysub[1,]$pedon_id, ")"),
     sub = "Assuming standard deviation proportional to boundary distinctness")
abline(v = c_stocks(mysub[1,])$c_stock_frags, lty=2)
legend("topleft", c("Simulated","Source"), lty=c(1,2))

prop.table(table(mysub_1_sim$c_stock_frags > c_stocks(mysub[1,])$c_stock_frags))
# ```{r}
# library(aqp)
# library(soilDB)
# 
# res <- fetchKSSL(mlra = "22A", returnMorphologicData = TRUE)
# res <- res$SPC
# sub <- aqp::subset(res, taxorder %in% c("Ultisols", "Alfisols"))
# coordinates(sub) <- ~ x + y
# proj4string(sub) <- "EPSG:4326"
# 
# bdy <- sf::st_as_sf(fetchSDA_spatial("CA649", by.col = "areasymbol", geom.src = "sapolygon"))
# bdy <- sf::st_transform(bdy, sf::st_crs(32610))
# bdy_buf <- sf::st_buffer(bdy, 3000)
# pts <- sf::st_transform(sf::st_as_sf(as(sub, 'SpatialPointsDataFrame')), sf::st_crs(32610))
# 
# plot(bdy_buf$geometry)
# plot(pts$geometry, add = TRUE)
# 
# y <- sf::st_intersection(pts, bdy_buf)
# aqp::plot(subset(sub, profile_id(sub) %in% y$pedon_key))
# ```
