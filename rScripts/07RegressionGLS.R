## ----message=FALSE------------------------------------------------------------
library(automap)
library(gstat)
library(nlme)
library(spdep)
library(ncf)
library(tidyverse)


## ----echo=FALSE---------------------------------------------------------------
set.seed(3515) #234 #104


## ----message=FALSE, warning=FALSE---------------------------------------------
# We will do n points
n <- 75

# CSR. Our mapped points.
easting <- runif(n,0,100)
northing <- runif(n,0,100)
points <- cbind(easting,northing)

# Generate the spatially autocorrelated error term, epsilon.
# Build neighbor objects with an inclusive distance threshold
dnb <- dnearneigh(x = points, d1 = 0, d2 = 150)
dsts <- nbdists(nb = dnb, coords = points)
# And follow an IDW process
p <- 2.25 # IDW power
idw <- lapply(dsts, function(x) x^-p)
# Weights matrix as a list
wList <- nb2listw(neighbours = dnb, glist=idw, style="W") 

# Compute SAR generating operator via a weights matrix
inv <- spatialreg::invIrW(x = wList, rho = 0.75)

# Calculate epsilon and rescale as z scores
epsilon <- inv %*% rnorm(n) 
epsilon <- scale(epsilon[,1])[,1]


## -----------------------------------------------------------------------------
moran.test(epsilon, wList)


## -----------------------------------------------------------------------------
# Generate x as noise.
x <- rnorm(n)
# Generate y as a function of x and epsilon (which is spatially autocorrelated)
B0 <- 3 # intercept
B1 <- 0.5 # slope
y <- B0 + B1*x + epsilon

# Make a data.frame
dat <- data.frame(easting,northing,y,x)


## -----------------------------------------------------------------------------
ggplot(data = dat, mapping = aes(x=easting,y=northing,
                                 size=x, fill=x)) + 
  geom_point(shape=21) +scale_fill_viridis_c() +
  guides(fill=guide_legend(), size = guide_legend()) +
  coord_fixed() +
  ggtitle("x")


## -----------------------------------------------------------------------------
ggplot(data = dat, mapping = aes(x=easting,y=northing,
                                 size=y, fill=y)) +
  geom_point(shape=21) + scale_fill_viridis_c() + 
  guides(fill=guide_legend(), size = guide_legend()) + 
  coord_fixed() + ggtitle("y")



## -----------------------------------------------------------------------------
# Fit a model using gls
glsNaive <- gls(y~x, dat)
summary(glsNaive)


## -----------------------------------------------------------------------------
# Add residuals to dat -- think epsilon
dat$glsNaiveResids <- residuals(glsNaive,type="normalized")
ggplot(data = dat, mapping = aes(x=easting,y=northing,
                                 size=glsNaiveResids, fill=glsNaiveResids)) +
  geom_point(shape=21) +
  scale_fill_gradient2() +
  guides(fill=guide_legend(), size = guide_legend()) +
  coord_fixed()


## -----------------------------------------------------------------------------
# Test them for spatial autocorrelation using a correlogram.
residsI <- spline.correlog(x=dat$easting, y=dat$northing,
                           z=dat$glsNaiveResids, resamp=50, quiet=TRUE)
plot(residsI,xlim=c(0,max(dist(points))/3))


## -----------------------------------------------------------------------------
datSF <- dat %>% st_as_sf(coords = c("easting","northing"))
# A variogram of the residuals
plot(autofitVariogram(glsNaiveResids~1, input_data = datSF, model = c("Gau","Sph","Exp")))


## -----------------------------------------------------------------------------
csSpatial <- corSpatial(form=~easting+northing,nugget=TRUE, type = "spherical")
glsUpdated <- update(glsNaive,correlation=csSpatial)
summary(glsUpdated)


## -----------------------------------------------------------------------------
# add the residuals to the spatial dat object
dat$glsResids <- residuals(glsUpdated,type="normalized")
# Map the residuals from gls
ggplot(data = dat, mapping = aes(x=easting,y=northing,
                                 size=glsResids, fill=glsResids)) +
  geom_point(shape=21) +
  scale_fill_gradient2() +
  guides(fill=guide_legend(), size = guide_legend()) +
  coord_fixed()


## -----------------------------------------------------------------------------
residsI <- spline.correlog(x=dat$easting, y=dat$northing,
                           z=dat$glsResids, resamp=50, quiet=TRUE)
plot(residsI,xlim=c(0,max(dist(points))/3))


## ----eval=FALSE,echo=FALSE----------------------------------------------------
# 
# library(nlme)
# library(ncf)
# library(sf)
# library(tidyverse)
# library(automap)
# 
# dat <- read.csv("birdDiv.csv")
# 
# glsNaive <- gls(birdDiv~plantDiv, dat)
# 
# summary(glsNaive)
# 
# dat$glsNaiveResids <- residuals(glsNaive,type="normalized")
# 
# residsIglsNaive <- spline.correlog(x=dat$UTME, y=dat$UTMN,xmax=8e04,
#                             z=dat$glsNaiveResids, resamp=50, quiet=TRUE)
# plot(residsIglsNaive)
# 
# # note that you can choose what models you fit in autofit vgm
# datSF <- dat %>% st_as_sf(coords=c("UTME","UTMN"))
# aVariogram <- autofitVariogram(glsNaiveResids~1,
#                                input_data = datSF,
#                                model = c("Gau","Exp","Sph"),
#                                verbose = TRUE)
# plot(aVariogram)
# # what is the max distance to look at for autocorrelation?
# # let's use 1/3 of the max point to point distance
# max(dist(cbind(dat$UTME, y=dat$UTMN)))/3 # about 84km
# 
# plot(spline.correlog(x=dat$UTME, y=dat$UTMN,xmax=84e03,
#                      z=dat$glsNaiveResids, resamp=100, quiet=TRUE))
# # def short distance ac
# 
# # try a bunch
# csGaus <- corSpatial(form=~UTME+UTMN,nugget=TRUE,type="gaus")
# csGausNoNugget <- corSpatial(form=~UTME+UTMN,nugget=FALSE,type="gaus")
# 
# glsGau <- update(glsNaive,correlation=csGaus)
# glsGauNoNug <- update(glsNaive,correlation=csGausNoNugget)
# 
# 
# csSph <- corSpatial(form=~UTME+UTMN,nugget=TRUE,type="sph")
# csSphNoNug <- corSpatial(form=~UTME+UTMN,nugget=FALSE,type="sph")
# 
# glsSph <- update(glsNaive,correlation=csSph)
# glsSphNoNug <- update(glsNaive,correlation=csSphNoNug)
# 
# csExp <- corSpatial(form=~UTME+UTMN,nugget=TRUE,type="exp")
# csExpNoNug <- corSpatial(form=~UTME+UTMN,nugget=FALSE,type="exp")
# 
# glsExp <- update(glsNaive,correlation=csExp) # nope
# glsExpNoNug <- update(glsNaive,correlation=csExpNoNug)
# 
# AIC(glsNaive,glsGau,glsGauNoNug,glsSphNoNug,glsExpNoNug)
# 
# # why no work? There is no sill!
# anova(glsNaive,glsGau)
# plot(Variogram(glsGau))
# plot(Variogram(glsGauNoNug))
# plot(Variogram(glsSph))
# plot(Variogram(glsSphNoNug))
# plot(Variogram(glsExp))
# plot(Variogram(glsExpNoNug))
# 
# 
# dat$glsSphResids <- residuals(glsSph,type="normalized")
# dat$glsGauResids <- residuals(glsGau,type="normalized")
# 
# 
# residsIglsSphResids <- spline.correlog(x=dat$UTME, y=dat$UTMN,xmax=8e04,
#                             z=dat$glsSphResids, resamp=50, quiet=TRUE)
# residsIglsGauResids <- spline.correlog(x=dat$UTME, y=dat$UTMN,xmax=8e04,
#                             z=dat$glsGauResids, resamp=50, quiet=TRUE)
# 
# par(mfcol=c(1,2))
# plot(residsIglsGauResids,main="GLS w/ corGaus Good!")
# plot(residsIglsSphResids,main="GLS w/ corSper Bad!")
# 

