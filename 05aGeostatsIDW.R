## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)
knitr::purl("05aGeostatsIDW.qmd", output = "05aGeostatsIDW.R", documentation = 1)



## ----message=FALSE------------------------------------------------------------
library(tidyverse)
library(sf)
library(gstat)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
foo <- data.frame(d = 1:100, w = (1:100)^-1)
foo %>% ggplot(mapping = aes(x=d,y=w)) + geom_line() +
  labs(x="Distance",y="Weight")


## -----------------------------------------------------------------------------
foo %>% ggplot(mapping = aes(x=d,y=w)) + geom_line() +
  labs(x="Distance",y="Weight") +
  scale_y_log10()


## -----------------------------------------------------------------------------
foo <- rbind(data.frame(d = 1:100, w = (1:100)^0, p = "0"),
             data.frame(d = 1:100, w = (1:100)^-0.5, p = "0.5"),
             data.frame(d = 1:100, w = (1:100)^-1, p = "1"),
             data.frame(d = 1:100, w = (1:100)^-1.5, p = "1.5"),
             data.frame(d = 1:100, w = (1:100)^-2, p = "2"),
             data.frame(d = 1:100, w = (1:100)^-2.5, p = "2.5"))
foo %>% ggplot(mapping = aes(x=d,y=w,color=p)) + geom_line() +
  labs(x="Distance",y="Weight")


## -----------------------------------------------------------------------------
foo %>% ggplot(mapping = aes(x=d,y=w,color=p)) + geom_line() +
  labs(x="Distance",y="Weight") +
  scale_y_log10()


## -----------------------------------------------------------------------------
foo <- data.frame(x = c(1,3,1,4,5),
                  y = c(5,4,3,5,1),
                  z = c(100,105,105,100,115))
foo
p1 <- foo %>% ggplot() + 
  geom_point(aes(x=x,y=y,size=z)) +
  lims(x=c(0,6),y=c(0,6))
p1


## -----------------------------------------------------------------------------
p1 + geom_point(aes(x=2,y=4),color="red",size=10,shape=0) +
  geom_point(aes(x=2,y=4),color="red",size=6,shape=63)


## ----echo=FALSE---------------------------------------------------------------
d2s0 <- as.matrix(dist(cbind(c(foo$x,2),c(foo$y,4))))[1:5,6]


## ----echo=FALSE---------------------------------------------------------------
p <- 2
w <- d2s0^-p


## ----echo=FALSE---------------------------------------------------------------
zhat_s0 <- sum(w*foo$z)/sum(w)


## ----warning=FALSE------------------------------------------------------------
# distance to missing point s0
d2s0 <- as.matrix(dist(cbind(c(foo$x,2),c(foo$y,4))))[1:5,6]
# power
p <- 2
# weights
w <- d2s0^-p
# and the estimation itself
zhat_s0 <- sum(w*foo$z)/sum(w)
zhat_s0
# add it to the plot
p1 + geom_point(aes(x=2,y=4,size=zhat_s0))


## -----------------------------------------------------------------------------
data(meuse.all)
glimpse(meuse.all)
class(meuse.all)
meuse.all$logLead <- log(meuse.all$lead)
# or for the tidyverse fans this is the same output
meuse.all <- meuse.all %>% mutate(logLead = log(lead))
# make into sf
meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

class(meuse_sf) # note change in class from data.frame to sf and data.frame

p2 <- ggplot(data = meuse_sf) +
  geom_sf(aes(fill=logLead), size=4, 
          shape = 21, color="white",alpha=0.8)+
  scale_fill_continuous(type = "viridis",name="log(ppm)") + 
  labs(title="Lead concentrations")
p2


## -----------------------------------------------------------------------------
meuse.grid <- readRDS("data/meuse.grid.Rds")
head(meuse.grid)


## -----------------------------------------------------------------------------
meuse_grid_sf <- st_as_sf(meuse.grid, 
                          coords = c("x","y"), 
                          crs = st_crs(meuse_sf))
meuse_grid_sf


## -----------------------------------------------------------------------------
idw_p2_model <- gstat(formula=logLead~1, 
                      locations = meuse_sf,
                      set=list(idp = 2))
logLeadIDW_p2_sf <- predict(idw_p2_model,meuse_grid_sf)
logLeadIDW_p2_sf


## -----------------------------------------------------------------------------
sf_2_rast <-function(sfObject,variableIndex = 1){
  # coerce sf to a data.frame
  dfObject <- data.frame(st_coordinates(sfObject),
                         z=as.data.frame(sfObject)[,variableIndex])
  # coerce data.frame to SpatRaster
  rastObject <- rast(dfObject,crs=crs(sfObject))
  
  names(rastObject) <- names(sfObject)[variableIndex]
  
  return(rastObject)
}

logLeadIDW_p2_rast <- sf_2_rast(logLeadIDW_p2_sf)
logLeadIDW_p2_rast
# and plot
ggplot() +
  geom_spatraster(data=logLeadIDW_p2_rast, mapping = aes(fill=var1.pred),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") + 
  labs(title="Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
ggplot() +
  geom_spatraster_contour_filled(data=logLeadIDW_p2_rast,
                                 breaks = seq(from=3.5, to=6.5,by=0.25),
                                 alpha = 0.9) +
  scale_fill_discrete(name="log(ppm)",na.value = "transparent") + 
  labs(title="Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
obs <- meuse_sf$logLead
preds <- extract(logLeadIDW_p2_rast, meuse_sf) %>% pull(var1.pred)
rsq <- cor(obs,preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse
ggplot() +
  geom_abline(slope=1,intercept = 0) +
  geom_point(aes(x=obs,y=preds)) + 
  coord_cartesian() + 
  labs(x="Observed Values",
       y="Predicted Values",
       title="Lead log(ppm)")


## -----------------------------------------------------------------------------
n <- nrow(meuse_sf)
rows4test <- sample(x = 1:n,size = n*0.25)
meuseTest <- meuse_sf[rows4test,]
meuseTrain <- meuse_sf[-rows4test,]

# note that we build the model with meuseTrain
idw_p2_model <- gstat(formula=logLead~1, 
                      locations = meuseTrain,
                      set=list(idp = 2))
logLeadIDW_p2_sf <- predict(idw_p2_model,meuse_grid_sf)

logLeadIDW_p2_rast <- sf_2_rast(logLeadIDW_p2_sf)

# and plot
ggplot() +
  geom_spatraster(data=logLeadIDW_p2_rast, mapping = aes(fill=var1.pred),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") + 
  labs(title="Lead concentrations", subtitle = "IDW with p=2") +
  theme_minimal()


## -----------------------------------------------------------------------------
# note use of meuseTest here
obs <- meuseTest$logLead
preds <- extract(logLeadIDW_p2_rast, meuseTest) %>% pull(var1.pred)
rsq <- cor(obs,preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse

ggplot() +
  geom_abline(slope=1,intercept = 0) +
  geom_point(aes(x=obs,y=preds)) + 
  coord_fixed(ratio=1, xlim = range(preds,obs),ylim = range(preds,obs)) +
  labs(x="Observed Values",
       y="Predicted Values",
       title="Lead log(ppm)")


## -----------------------------------------------------------------------------
rmseNULL <- sqrt(mean((mean(meuseTrain$logLead) - obs)^2))
rmseNULL
1 - (rmse / rmseNULL)


## -----------------------------------------------------------------------------
idw_p3_model <- gstat(formula=logLead~1, 
                        locations = meuseTrain,
                        set=list(idp = 3))

logLeadIDW_p3_sf <- predict(idw_p3_model,meuse_grid_sf)

logLeadIDW_p3_rast <- sf_2_rast(logLeadIDW_p3_sf)

# and plot
ggplot() +
  geom_spatraster(data=logLeadIDW_p3_rast, mapping = aes(fill=var1.pred),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") + 
  labs(title="Lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
# note use of meuseTest here
obs <- meuseTest$logLead
preds <- extract(logLeadIDW_p3_rast, meuseTest) %>% pull(var1.pred)
rsq <- cor(obs,preds)^2
rmse <- sqrt(mean((preds - obs)^2))
rsq
rmse
ggplot() +
  geom_abline(slope=1,intercept = 0) +
  geom_point(aes(x=obs,y=preds)) + 
  coord_fixed(ratio=1, xlim = range(preds,obs),ylim = range(preds,obs)) +
  labs(x="Observed Values",
       y="Predicted Values",
       title="Lead log(ppm)")


## -----------------------------------------------------------------------------
1 - (rmse / rmseNULL)


## ----eval=FALSE, echo=FALSE---------------------------------------------------
# #Optimizing p 3.058529
# 
# f = function(idp, formula, data,...){
#   res <- sum(krige.cv(formula,data,set=list(debug=0,idp=idp),nfold=10,...)$residual**2)
#   res
# }
# optimize(f, interval=c(0.01,4), formula=log(lead)~1, data=meuse_sf)


## -----------------------------------------------------------------------------
# precip point data
prcpCA <- readRDS("data/prcpCA.rds")
# empty grid to interpolate into
gridCA <- readRDS("data/gridCA.rds")

# make as sf
prcpCA <- prcpCA %>% st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(value = 3310)

# simple map -- see postscript below for a fancy map
prcpCA %>% ggplot() + 
  geom_sf(aes(fill=ANNUAL,size=ANNUAL),color="white",
          shape=21,alpha=0.8) + 
  scale_fill_continuous(type = "viridis",name="mm") + 
  labs(title="Total Annual Precipitation") +
  scale_size(guide="none")


## ----echo=FALSE, eval=FALSE---------------------------------------------------
# gridCA <- gridCA %>% st_as_sf(coords = c("X", "Y")) %>%
#   st_set_crs(value = 3310)
# 
# idw_p2.5_model <- gstat(formula=ANNUAL~1,locations=prcpCA,
#                         set=list(idp = 2.5))
# 
# prcpIDW_p2.5_sf <- predict(idw_p2.5_model,gridCA)
# 
# prcpIDW_p2.5_rast <- sf_2_rast(prcpIDW_p2.5_sf)
# 
# # and plot
# ggplot() +
#   geom_spatraster(data=prcpIDW_p2.5_rast, mapping = aes(fill=var1.pred),alpha=0.8) +
#   scale_fill_continuous(type = "viridis",name="mm",na.value = "transparent") +
#   labs(title="Total Annual Precip") +
#   theme_minimal()
# 
# 
# library(tmap)
# 
# tmap_mode("view")
# tm_shape(prcpIDW_p2.5_rast) +
#   tm_raster(col="var1.pred",alpha = 0.5) + #title = "TAP (mm)"
#   tm_shape(prcpCA) +
#   tm_symbols(col="deeppink",alpha=0.5,size = "ANNUAL")
# 

