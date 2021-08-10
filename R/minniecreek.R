library(aqp)
library(soilDB)
library(magrittr)

# get data from NASIS as SoilProfileCollection
# 2021CA0436001,2021CA0436004,2021CA0436005,2021CA0436016

f <- fetchNASIS()

# create a data.frame (site data denormalized to (repeated at) horizon level)
h <- as(f, 'data.frame')

# filter horizons to keep only those with sand, silt and clay populated
# h.sub <- h[complete.cases(h[,c("sand","silt","clay")]),]
h.sub <- h

## install pedotransfR packge to estimate Ksat from sand, silt, clay using Saxton (1986) PTF
# remotes::install_github('ncss-tech/pedotransfR')
ldx <-  complete.cases(h.sub[,c("sand","clay")])
h$ksat <- NA
h.sub$ksat[ldx] <- pedotransfR::saxton_1986_Ksat(h.sub$sand[ldx], h.sub$clay[ldx])

# promote to SoilProfileCollection
depths(h.sub) <- peiid ~ hzdept + hzdepb

# for this demo, we will use only pedons correlated as Minniecreek
p <- subset(h.sub, taxonname == "Minniecreek")
p$genhz <- NULL

# apply boundary standard deviations based on [min. transition thickness] / 3
#   - "approximates" the Empirical rule, in that 99.7% of data falls within 3 SD of source (mean=0)
bound.lut <- c(
  'V' = 0.5,
  'A' = 2,
  'C' = 5,
  'G' = 15,
  'D' = 45
) / 3

# first letter of class is the label
p$boundsd <- bound.lut[toupper(substr(p$bounddistinct, 0, 1))]

# if the boundary is not populated, variance is 0
p$boundsd[is.na(p$boundsd)] <- 0

# iterate over each source profile
#  - perturb it 100 times
#  - truncate to [0,100]
#  - calculate geometric mean ksat (for instance)
#  
res <- combine(profileApply(p, function(x) {
    res2 <- perturb(x, boundary.attr = "boundsd", n = 100) %>% 
      trunc(0, 100) %>% 
      mutate_profile(ksat_mean_geom = exp(log(ksat)))
    profile_id(res2) <- paste(profile_id(x), profile_id(res2))
    res2
  }))

# density plot of geometric mean ksat [0,100]
plot(density(res$ksat_mean_geom, na.rm = TRUE))

# add vertical lines showing values estimated for source layers only
lines(density(p$ksat, na.rm=T))
abline(v = p$ksat, lty = 2)

# combine perturbed profiles with source profiles
res <- combine(p, res)

# calculate generalized horizon labels
p$genhz <- generalize.hz(p$hzname, 
                           new = c("Oi","A", "BA", "Bt"), 
                           pat = c("O","^A[1-9]?$", 
                                   "^BA|^AB", 
                                   "^[2-9]?B[^A]*t"))

res$genhz <- generalize.hz(res$hzname, 
                           new = c("A", "BA", "Bt"), 
                           pat = c("^A[1-9]?$", 
                                   "^BA|^AB", 
                                   "^[2-9]?B[^A]*t"))

# set generalized horizon as SPC designation name column
hzdesgnname(res) <- "genhz"
site(p)<- ~pedon_id
par(mar=c(1,1,4,1))
# minniecreek_genhz.png
plotSPC(p, color="genhz", label="pedon_id", id.style = "side", cex.names = 0.8)

# helper function to calculate thickness of horizons matching a pattern within a profile
hzDesgnThickness <- function(p, pattern) {
  hzd <- horizonDepths(p)
  maxDepthOf(p, pattern, top = FALSE)[[hzd[2]]] - minDepthOf(p, pattern)[[hzd[1]]]
}

# calculate thickness of A horizons
ahzthk <- hzDesgnThickness(res, "^A$")
bahzthk <- hzDesgnThickness(res, "^BA$")
bthzthk <- hzDesgnThickness(res, "^Bt$")

# total thickness is the sum of the above
totalthk <- ahzthk + ifelse(is.na(bahzthk), 0, bahzthk) + bthzthk

# 

# make base density plot using A horizon data
plot(density(ahzthk, bw=1, kernel = "rectangular"), xlim = c(0, 100), 
     col = "green",
     ylim = c(0, 0.3),
     main = "Observed v.s. Simulated Horizon Thickness Distribution")

# add lines for BA, Bt and total
lines(density(bahzthk, bw=1, kernel = "rectangular"), col = "red")
lines(density(bthzthk, bw=1, kernel = "rectangular"), col = "blue")
lines(density(totalthk, bw=1, kernel = "rectangular"), col = "black")

# add source pedon horizon thicknesses as vertical lines
abline(v = ahzthk[profile_id(res) %in% profile_id(p)], lty = 2, col = "green")
abline(v = bahzthk[profile_id(res) %in% profile_id(p)], lty = 2, col = "red")
abline(v = bthzthk[profile_id(res) %in% profile_id(p)], lty = 2, col = "blue")
abline(v = totalthk[profile_id(res) %in% profile_id(p)], lty = 2, col = "black")

# add explanation
legend("topright", legend = c("A [observed]",
                              "A [simulated]",
                              "BA [observed]",
                              "BA [simulated]",
                              "Bt [observed]",
                              "Bt [simulated]",
                              "Total [observed]",
                              "Total [simulated]"),
       col = c("green", "green", "red", "red", "blue", "blue", "black", "black"),  
       lty = c(2, 1, 2, 1, 2, 1, 2, 1))
