## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)


## ----message=FALSE------------------------------------------------------------
library(sf)
library(gstat)
library(tidyverse)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
# load
data(meuse.all)
meuse.grid <- readRDS("../data/meuse.grid.Rds")

# make a variable to work with
meuse.all$logLead <- log(meuse.all$lead)
# make into sf
meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

meuse_grid_sf <- st_as_sf(meuse.grid, 
                          coords = c("x","y"), 
                          crs = st_crs(meuse_sf))
meuse_grid_sf

p1 <- ggplot(data = meuse_sf) +
  geom_sf(aes(fill=logLead), size=4, 
          shape = 21, color="white",alpha=0.8)+
  scale_fill_continuous(type = "viridis",name="log(ppm)") + 
  labs(title="Lead concentrations")
p1


## -----------------------------------------------------------------------------
leadVar <- variogram(logLead~1, meuse_sf)
plot(leadVar,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Lead concentrations (log(ppm))")


## -----------------------------------------------------------------------------
# note our initial estimates for the partial sill, range, and nugget. 
sph.model <- vgm(psill=0.5, model="Sph", range=750, nugget=0.05)
sph.fit <- fit.variogram(object = leadVar, model = sph.model)
sph.fit # look at the fitted values
plot(leadVar,model=sph.fit,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Lead concentrations (log(ppm))",
     sub="Points: Empirical, Line: Spherical Model")


## -----------------------------------------------------------------------------
# note our initial estimates for the sill, range, and nugget
exp.model <- vgm(psill=0.5, model="Exp", range=750, nugget=0.05)
exp.fit <- fit.variogram(object = leadVar, model = exp.model)
plot(leadVar,model=exp.fit,pch=20,cex=1.5,col="black",
     ylab=expression("Semivariance ("*gamma*")"),
     xlab="Distance (m)", main = "Lead concentrations (log(ppm))",
     sub="Points: Empirical, Line: Exponential Model")


## -----------------------------------------------------------------------------
foo <- data.frame(x = c(1,3,1,4,5),
                  y = c(5,4,3,5,1),
                  z = c(100,105,105,100,115))
foo
p2 <- ggplot() + 
  geom_point(data=foo,aes(x=x,y=y,size=z)) + 
  lims(x=c(0,6),y=c(0,6))
p2


## -----------------------------------------------------------------------------
p2 <- p2 + 
  geom_point(aes(x=2,y=4),color="red",size=10,shape=0) +
  geom_point(aes(x=2,y=4),color="red",size=6,shape=63)
p2


## ----echo=FALSE---------------------------------------------------------------
d2s0 <- as.matrix(dist(cbind(c(foo$x,2),c(foo$y,4))))[1:5,6]


## ----echo=FALSE---------------------------------------------------------------
g <- 0 + 13.5*d2s0


## ----echo=FALSE---------------------------------------------------------------
d <- dist(cbind(foo$x,foo$y))


## ----echo=FALSE---------------------------------------------------------------
G <- 0 + 13.5*d


## ----echo=FALSE---------------------------------------------------------------
#round(solve(G),3)
lambda <- solve(G) %*% g
zhat <- sum(lambda[,1] * foo$z)


## -----------------------------------------------------------------------------
leadVar <- variogram(logLead~1, meuse_sf)
leadModel <- vgm(psill=0.6, model="Sph", range=750, nugget=0.05)
leadFit <- fit.variogram(object = leadVar, model = leadModel)
leadGstat <- gstat(formula = logLead~1, locations = meuse_sf, 
                   model = leadFit)
leadKrige_sf <- predict(leadGstat,newdata = meuse_grid_sf)
leadKrige_sf


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

leadKrige_rast <- sf_2_rast(leadKrige_sf)

leadKrige_rast

# and plot
ggplot() +
  geom_spatraster(data=leadKrige_rast, mapping = aes(fill=var1.pred),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") + 
  labs(title="Lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
leadKrige_sf$var1.var.sqrt <- sqrt(leadKrige_sf$var1.var)

leadKrige_rast <- sf_2_rast(leadKrige_sf,variableIndex = 4)
leadKrige_rast

# and plot
ggplot() +
  geom_spatraster(data=leadKrige_rast, mapping = aes(fill=var1.var.sqrt ),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") + 
  labs(title="Variance of lead concentrations") +
  theme_minimal()


## -----------------------------------------------------------------------------
vgm()


## -----------------------------------------------------------------------------
show.vgms(models = c("Exp", "Mat", "Gau", "Sph"))


## -----------------------------------------------------------------------------
leadKrigeLOOCV_sf <- krige.cv(formula = logLead~1, 
                           locations = meuse_sf, 
                           model = leadFit, verbose = FALSE)
leadKrigeLOOCV_sf
# CV R2
cor(leadKrigeLOOCV_sf$observed,leadKrigeLOOCV_sf$var1.pred)^2


## -----------------------------------------------------------------------------
# precip point data
prcpCA <- readRDS("../data/prcpCA.rds")
# empty grid to interpolate into
gridCA <- readRDS("../data/gridCA.rds")

prcpCAsf <- prcpCA %>% st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(value = 3310)

prcpCAsf %>% ggplot() + 
  geom_sf(aes(fill=ANNUAL,size=ANNUAL),color="white",
          shape=21,alpha=0.8) + 
  scale_fill_continuous(type = "viridis",name="mm") + 
  labs(title="Total Annual Precipitation") +
  scale_size(guide="none")


## -----------------------------------------------------------------------------
library(automap)
leadVar <- autofitVariogram(formula = logLead~1,input_data = meuse_sf)
summary(leadVar)
plot(leadVar)


## -----------------------------------------------------------------------------
leadAutoKrigeLOOCV <- autoKrige.cv(formula = logLead~1, input_data = meuse_sf,
                                   verbose = c(FALSE,FALSE))
summary(leadAutoKrigeLOOCV)
# R2
cor(leadAutoKrigeLOOCV$krige.cv_output$observed,leadAutoKrigeLOOCV$krige.cv_output$var1.pred)^2

