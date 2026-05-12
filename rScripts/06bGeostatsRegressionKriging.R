## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(733)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(gstat)


## -----------------------------------------------------------------------------
n <- 50
b0 <- 4
b1 <- 0.6
b2 <- 3
originalDat <- data.frame(x1=rnorm(n), x2=rnorm(n), epsilon = rnorm(n,sd = 2))
originalDat$y <- b0 + b1 * originalDat$x1 + b2 * originalDat$x2 + originalDat$epsilon


## -----------------------------------------------------------------------------
lm1 <- lm(y~x1 + x2,data = originalDat)
summary(lm1)


## -----------------------------------------------------------------------------
newDat <- data.frame(x1=rnorm(n),x2=rnorm(n))
head(newDat)


## -----------------------------------------------------------------------------
newDat$yhat <- predict(object = lm1,newdata = newDat)
head(newDat)


## -----------------------------------------------------------------------------
data(meuse.all)
str(meuse.all)


## -----------------------------------------------------------------------------
meuse.grid <- readRDS("../data/meuse.grid.Rds")
str(meuse.grid)


## -----------------------------------------------------------------------------
meusePoints_sf <- st_as_sf(meuse.all,coords = c("x","y"), crs=28992)
meuseGrid_sf <- st_as_sf(meuse.grid,coords = c("x","y"), crs=28992)
meusePoints_sf
meuseGrid_sf


## -----------------------------------------------------------------------------
# get the relative distance from meuseGrid_sf and add it to meusePoints_sf
meuseGrid_dist_sf <- meuseGrid_sf %>% select(dist)
# join the two objects by the nearest point
meusePoints_sf <- st_join(meusePoints_sf,meuseGrid_dist_sf, join = st_nearest_feature)
# note that we now have relative dist for each point now in the column `dist`
meusePoints_sf


## -----------------------------------------------------------------------------
meusePoints_sf <- meusePoints_sf %>% mutate(sqrt_dist = sqrt(dist))
meuseGrid_sf <- meuseGrid_sf %>% mutate(sqrt_dist = sqrt(dist))



## -----------------------------------------------------------------------------
meusePoints_sf <- meusePoints_sf %>% drop_na(om)


## -----------------------------------------------------------------------------
meusePoints_sf <- meusePoints_sf %>% mutate(soil = factor(soil))
meuseGrid_sf <- meuseGrid_sf %>% mutate(soil = factor(soil))


## -----------------------------------------------------------------------------
om_lm <- lm(om~sqrt_dist+soil, data=meusePoints_sf)
summary(om_lm)
anova(om_lm)


## -----------------------------------------------------------------------------
meuseGrid_sf$omHat <- predict(om_lm, newdata=meuseGrid_sf)

ggplot() + geom_sf(data = meuseGrid_sf,mapping = aes(color=omHat))


## ----echo=FALSE---------------------------------------------------------------
w <- spdep::knn2nb(spdep::knearneigh(meusePoints_sf,k=4))
foo <- spdep::moran.test(residuals(om_lm),spdep::nb2listw(w))


## -----------------------------------------------------------------------------
sf_2_rast <-function(sfObject,variable2get = 1){
  # coerce sf to a data.frame
  tmp <- sfObject[,variable2get] %>% st_drop_geometry()
  dfObject <- data.frame(st_coordinates(sfObject),
                         z=tmp)
  # coerce data.frame to SpatRaster
  rastObject <- rast(dfObject,crs=crs(sfObject))
  
  return(rastObject)
}

meuseGrid_rast <- sf_2_rast(meuseGrid_sf,variable2get = "omHat")
ggplot() + geom_spatraster(data = meuseGrid_rast,
                           mapping = aes(fill=omHat)) +
  scale_fill_terrain_c() +
  labs(fill = "omHat")
  


## -----------------------------------------------------------------------------
omVar <- variogram(om~1, meusePoints_sf)
plot(omVar,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "% Soil Organic Matter")


## -----------------------------------------------------------------------------
omGstat <- gstat(id = "omModel", formula = om ~ sqrt_dist + soil, 
                 data = meusePoints_sf)
omGstat_obsVariogram <- variogram(omGstat)
plot(omGstat_obsVariogram,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Model Residuals")


## -----------------------------------------------------------------------------
omGau_fittedVariogram <- fit.variogram(omGstat_obsVariogram, vgm(psill = 6, model = "Gau", range = 500, nugget = 4))

plot(omGstat_obsVariogram, omGau_fittedVariogram,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Model Residuals")


## -----------------------------------------------------------------------------
# Update the gstat object with the variogram:
omGstat_w_variogram <- gstat(omGstat, id="omModel", model=omGau_fittedVariogram)
# And predict
omHat_sf <- predict(omGstat_w_variogram, newdata = meuseGrid_sf)
omHat_sf


## -----------------------------------------------------------------------------
omHat_rast <- sf_2_rast(omHat_sf,variable2get = "omModel.pred")
omHat_rast


## -----------------------------------------------------------------------------
ggplot() + 
  geom_spatraster(data = omHat_rast,
                           mapping = aes(fill=omModel.pred)) +
  scale_fill_terrain_c() +
  labs(fill = "Organic Matter (%)")



## -----------------------------------------------------------------------------
omHatVar_rast <- sf_2_rast(omHat_sf,variable2get = "omModel.var")
# take the square root of variance to create standard deviation (units match then)
omHatVar_rast <- omHatVar_rast %>% mutate(omModel.var.sqrt = sqrt(omModel.var))
ggplot() + geom_spatraster(data = omHatVar_rast,
                           mapping = aes(fill=omModel.var.sqrt)) +
  scale_fill_terrain_c() +
  labs(fill = "Organic Matter Variance (%)")



## ----eval=FALSE,echo=FALSE,message=FALSE, warning=FALSE-----------------------
# # the rk
# out <- gstat.cv(object=omGstat_w_variogram,nfold=10)
# #r2
# cor(out$observed, out$observed - out$residual)^2
# cor(out$observed, out$omModel.pred)^2 # same
# 
# # the lm
# out2 <- gstat.cv(omGstat,nfold=10)
# #r2
# cor(out2$observed, out2$observed - out2$residual)^2
# cor(out2$observed, out2$omModel.pred)^2 # same
# 
# # oh! by fold. this is better prob
# 
# out@data %>% group_by(fold) %>% summarise(r2=cor(observed, omModel.pred)^2) %>% summarise(mean(r2))
# 
# out2@data %>% group_by(fold) %>% summarise(r2=cor(observed, omModel.pred)^2) %>% summarise(mean(r2))
# 
# # same as lm?
# library(caret)
# trainControl <- trainControl(method = "cv", number = 10)
# fit <- train(om~sqrt_dist+soil,data=meusePoints_sf, trControl = trainControl, method = "lm")
# fit

