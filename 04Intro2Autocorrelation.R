## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(1984)
knitr::purl("04Intro2Autocorrelation.qmd", output = "04Intro2Autocorrelation.R", documentation = 1)

## ----message=FALSE------------------------------------------------------------
library(sf)
library(tidyverse)
library(tmap)
library(gstat)
library(ncf)
library(spdep)
library(tmap)
birds_sf <- readRDS("data/birdRichnessMexico.rds")
tmap_mode("view")
tm_shape(birds_sf) + 
  tm_symbols(col="nSpecies", alpha = 0.7)


## -----------------------------------------------------------------------------
birdsDF <- data.frame(st_coordinates(birds_sf),nSpecies=birds_sf$nSpecies)


## -----------------------------------------------------------------------------
max(dist(st_coordinates(birds_sf)))/3

hist(birds_sf$nSpecies) # fine

birdVar <- variogram(nSpecies~1, data = birds_sf)
plot(birdVar)
birdVar <- variogram(nSpecies~1, data = birds_sf, alpha=c(0,90))
plot(birdVar)
birdsDF <- data.frame(st_coordinates(birds_sf),nSpecies=birds_sf$nSpecies)
zI <- correlog(x=birdsDF[,1],y=birdsDF[,2],
               z =birdsDF$nSpecies,increment = 50,resamp = 1000)
plot(foo) # yucky!
plot(foo,xlim=c(0,1e3))
abline(h=0)
pBirdsI <-data.frame(I = zI$correlation, 
                     d = zI$mean.of.class,
                     p = zI$p) %>%
  filter(d <= 1000) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),linewidth=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=3,pch=21) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (km)",y="Moran's I",
       caption = "Crit value of p<0.001")
pBirdsI


## -----------------------------------------------------------------------------
transect <- data.frame(
  x = seq(0, 700, by = 100),
  z = c(10, 12, 11, 13, 5, 4, 6, 5)
)


## -----------------------------------------------------------------------------
ggplot(transect, aes(x = x, y = z)) +
  geom_line() +
  geom_point(size = 4) +
  labs(x = "Location (m)", y = "z")


## -----------------------------------------------------------------------------
# Step 1: compute all pairwise distances between stations.
# dist() returns a distance object; we convert it to a plain matrix
# so we can index into it like a spreadsheet.
Dmat <- as.matrix(dist(transect$x))

# Step 2: we only want each pair once (i,j is the same pair as j,i),
# so we pull out just the upper triangle of the distance matrix.
# upper.tri() returns a logical matrix -- TRUE above the diagonal, FALSE elsewhere.
DmatUpperTri <- upper.tri(Dmat)

# Step 3: find the row and column positions of those TRUE cells.
# Each row of idx is one pair: idx[,1] is the i station, idx[,2] is the j station.
idx <- which(DmatUpperTri, arr.ind = TRUE)

# Step 4: build the pairs data frame.
# For each pair, grab the z value at station i, the z value at station j,
# and the distance between them from the matrix.
pairs_df <- data.frame(
  zi   = transect$z[idx[,1]],  # z at station i
  zj   = transect$z[idx[,2]],  # z at station j
  dist = Dmat[DmatUpperTri]     # distance between i and j
)
pairs_df


## -----------------------------------------------------------------------------
# compute r for each group to use as labels in the plot
cor_labels <- pairs_df %>%
  filter(dist <= 100 | dist > 200) %>%
  mutate(group = ifelse(dist <= 100, "Close (≤ 100m)", "Far (≥ 200m)")) %>%
  group_by(group) %>%
  summarise(r = round(cor(zi, zj), 2)) %>%
  mutate(label = paste0("r = ", r))

pairs_df %>%
  filter(dist <= 100 | dist > 200) %>%
  mutate(group = ifelse(dist <= 100, "Close (≤ 100m)", "Far (≥ 200m)")) %>%
  ggplot(aes(x = zi, y = zj)) +
  geom_point(size = 3) +
  geom_text(data = cor_labels, aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
  facet_wrap(~group) +
  labs(x = expression(z[i]), y = expression(z[j]))


## -----------------------------------------------------------------------------
data(meuse.all)
class(meuse.all)
head(meuse.all)


## -----------------------------------------------------------------------------
meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)


## -----------------------------------------------------------------------------
meuse_sf$log_lead <- log(meuse_sf$lead)


## -----------------------------------------------------------------------------
ggplot(data = meuse_sf) + 
  geom_sf(mapping = aes(fill=lead,size=lead),shape=21,alpha=0.6) +
  scale_fill_continuous(type = "viridis",name="ppm")


## -----------------------------------------------------------------------------
tmap_mode('view')
tm_shape(meuse_sf) +
  tm_bubbles(
    size        = "lead",
    fill        = "lead",
    fill.scale  = tm_scale_intervals(values = "viridis"),
    fill_alpha  = 0.7,
    fill.legend = tm_legend(title = "Lead (ppm)"),
    size.legend = tm_legend(title = "Lead (ppm)"),
    id          = "lead"
  ) +
  tm_scalebar(position = c("left", "bottom"))


## -----------------------------------------------------------------------------
leadVarCloud <- variogram(log_lead~1, locations = meuse_sf, cloud = TRUE)
plot(leadVarCloud,pch=20,cex=1.5,col="black",alpha=0.1,
     ylab=expression(Semivariance~(gamma)),
     xlab="Distance (m)",main = "Log lead concentrations")


## -----------------------------------------------------------------------------
leadVar <- variogram(log_lead~1, locations = meuse_sf, cloud = FALSE)
plot(leadVar,pch=20,cex=1.5,col="black",
     ylab=expression(Semivariance~(gamma)),
     xlab="Distance (m)", main = "Log lead concentrations")


## -----------------------------------------------------------------------------
summary(lm(log_lead~1,data=meuse_sf))


## -----------------------------------------------------------------------------
mean(meuse_sf$log_lead)
sd(meuse_sf$log_lead)/sqrt(nrow(meuse_sf))


## ----echo=FALSE,fig.width=9,fig.height=3--------------------------------------
n <- 20

chess <- expand.grid(x = 1:n, y = 1:n) %>%
  mutate(z = ifelse((x + y) %% 2 == 0, 1, 0),
         pattern = "Moran's I ≈ -1\nChessboard")

set.seed(42)
random <- expand.grid(x = 1:n, y = 1:n) %>%
  mutate(z = runif(n^2),
         pattern = "Moran's I ≈ 0\nRandom noise")

gradient <- expand.grid(x = 1:n, y = 1:n) %>%
  mutate(z = x / n,
         pattern = "Moran's I ≈ 1\nGradient")

bind_rows(chess, random, gradient) %>%
  mutate(pattern = factor(pattern, levels = c(
    "Moran's I ≈ -1\nChessboard",
    "Moran's I ≈ 0\nRandom noise",
    "Moran's I ≈ 1\nGradient"))) %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~pattern) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))


## ----warning=FALSE------------------------------------------------------------
# distance matrix 
D <- dist(st_coordinates(meuse_sf))
# inverse distance matrix
W <- as.matrix(1/D)
# convert a the weights matrix to a weights list object
# so spdep is happy
wList <- mat2listw(W)
# calculate I
moran.test(meuse_sf$log_lead,wList)


## ----warning=FALSE------------------------------------------------------------
x <- as.vector(as.matrix(dist(st_coordinates(meuse_sf))))
y <- as.vector(W)
ggplot() + geom_line(aes(x=x[x>0],y=y[x>0])) +
  labs(x="Distance (m)",y="Inverse distance (aka W)") +
  lims(x=c(0,1500))


## ----warning=FALSE------------------------------------------------------------
D <- as.matrix(dist(st_coordinates(meuse_sf)))
# make an empty weights matrix
W <- matrix(0,ncol=ncol(D),nrow=nrow(D))
# set some values to 1 using a logical mask of D<100
W[D<100] <- 1
# convert a the weights matrix to a weights list object
# so spdep is happy
wList <- mat2listw(W)
# calculate I
moran.test(meuse_sf$log_lead,wList)


## ----warning=FALSE------------------------------------------------------------
# make an empty weights matrix of the right dimensions
W <- matrix(0,ncol=ncol(D),nrow=nrow(D))
# set some values to 1
W[D>500 & D<=1000] <- 1
# convert a the weights matrix to a weights list object
# so spdep is happy
wList <- mat2listw(W)
# calculate I
moran.test(meuse_sf$log_lead,wList)


## ----warning=FALSE------------------------------------------------------------
distanceInterval <- 100
distanceVector <- seq(0,1500,by=distanceInterval)
n <- length(distanceVector)
D <- as.matrix(dist(st_coordinates(meuse_sf)))
# make an object to hold results
res <- data.frame(midBin=rep(NA,n-1),I=rep(NA,n-1))
for(i in 2:n){
  W <- matrix(0,ncol=ncol(D),nrow=nrow(D))
  # set some values to 1
  W[D >= distanceVector[i-1] & D < distanceVector[i]] <- 1
  # convert a the weights matrix to a weights list object
  # so spdep is happy
  wList <- mat2listw(W)
  # calculate I
  res$I[i-1] <- moran.test(meuse_sf$log_lead,wList,zero.policy=TRUE)$estimate[1]
  # centered distance bin
  res$midBin[i-1] <- distanceVector[i] - distanceInterval/2
}
ggplot(data=res, mapping = aes(x=midBin,y=I)) + 
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_line() + geom_point(size=3) +
  labs(x="Distance (m)",y="Moran's I")



## -----------------------------------------------------------------------------
W <- knn2nb(knearneigh(meuse_sf,k=8))
moran.test(meuse_sf$log_lead,nb2listw(W))


## ----warning=FALSE------------------------------------------------------------
n <- 7
res <- data.frame(k=2^(1:n),I=rep(NA,n))
for(i in 1:n){
  W <- knn2nb(knearneigh(meuse_sf,k=2^i))
  res$I[i] <- moran.test(meuse_sf$log_lead,nb2listw(W))$estimate[1]
}
ggplot(data=res, mapping = aes(x=k,y=I)) + 
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_line() + geom_point(size=3) +
  labs(x="K neighbors",y="Moran's I")


## -----------------------------------------------------------------------------
meuseX <- st_coordinates(meuse_sf)[,1]
meuseY <- st_coordinates(meuse_sf)[,2]
meuseLead <- meuse_sf$log_lead


## -----------------------------------------------------------------------------
leadI <- correlog(x=meuseX, y=meuseY, z=meuseLead,
                  increment=100, resamp=100, quiet=TRUE)
plot(leadI,xlim=c(0,1500))
abline(h=0,lty="dashed")


## -----------------------------------------------------------------------------
plot(leadI)
abline(h=0,lty="dashed")


## -----------------------------------------------------------------------------
max(meuseX) - min(meuseX)
max(meuseY) - min(meuseY)


## -----------------------------------------------------------------------------
meuseD <- dist(cbind(meuseX,meuseY))
max(meuseD)


## -----------------------------------------------------------------------------
data.frame(n=leadI$n,
           I = leadI$correlation, 
           d = leadI$mean.of.class,
           p = leadI$p) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),size=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=5,pch=21) +
  lims(x=c(0,1500),y=c(-0.6,0.6)) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (m)",y="Moran's I",
       title = "Autocorrelation of log(Lead)",subtitle = "Crit value of p<0.01")


## ----echo=FALSE, message=FALSE------------------------------------------------
library(fields)
library(plotly)

n <- 100
nVec <- 0:(n-1)

###################################################################
#
## Sine Wave 10
#
###################################################################

amp <- 0.5
nBumps <- 10
nBumps2 <- n/nBumps
z <- amp * sin(2 * pi / nBumps2 * nVec)
Z <- matrix(rep(z,n),nrow=n, ncol=n)
Z <- Z + t(Z)
eps <- matrix(abs(rnorm(n^2,sd=0.1)),nrow=n, ncol=n)
Z <- Z + eps
# make flat
Z[Z<0] <- 0
foo <- expand.grid(x=nVec,y=nVec)
foo$z <- as.vector(Z)

p10BumpsMap <- ggplot(foo,aes(x=x,y=y,fill=z)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  coord_fixed() + 
  labs(x="Easting (m)", y="Northing (m)") + 
  theme_minimal() +
  theme(legend.position = "none")


fooMat <- matrix(foo$z,n,n)
p10BumpsPersp <- plot_ly(z = fooMat) %>% 
  add_surface() %>%
  hide_colorbar() %>%
  layout(title = "",
         scene = list(
           xaxis = list(title = "Easting (m)"),
           yaxis = list(title = "Northing (m)"),
           camera = list(eye = list(x = 1.95, y = -1.25, z = 1.25))
         ))

samps2get <- sample(n^2,size = 500)
bar <- foo[samps2get,]

zI <- correlog(x=bar$x, y=bar$y, z=bar$z, increment = 1, resamp = 200, quiet = TRUE)

p10BumpsI <- data.frame(I = zI$correlation, 
                        d = zI$mean.of.class,
                        p = zI$p) %>%
  filter(I < 1 & d <= 50) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),size=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=3,pch=21) +
  lims(x=c(0,50))+ #,y=c(-0.6,0.6)) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (m)",y="Moran's I",
       caption = "Crit value of p<0.01")

###################################################################
#
## Sine Wave 5
#
###################################################################

nBumps <- 5
nBumps2 <- n/nBumps
z <- amp * sin(2 * pi / nBumps2 * nVec)
Z <- matrix(rep(z,n),nrow=n, ncol=n)
Z <- Z + t(Z)
eps <- matrix(abs(rnorm(n^2,sd=0.1)),nrow=n, ncol=n)
Z <- Z + eps
# make flat
Z[Z<0] <- 0
foo <- expand.grid(x=nVec,y=nVec)
foo$z <- as.vector(Z)

p5BumpsMap <- ggplot(foo,aes(x=x,y=y,fill=z)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  coord_fixed() + 
  labs(x="Easting (m)", y="Northing (m)") + 
  theme_minimal() +
  theme(legend.position = "none")

fooMat <- matrix(foo$z,n,n)
p5BumpsPersp <- plot_ly(z = fooMat) %>% 
  add_surface() %>%
  hide_colorbar() %>%
  layout(title = "",
         scene = list(
           xaxis = list(title = "Easting (m)"),
           yaxis = list(title = "Northing (m)"),
           camera = list(eye = list(x = 1.95, y = -1.25, z = 1.25))
         ))

samps2get <- sample(n^2,size = 500)
bar <- foo[samps2get,]
zI <- correlog(x=bar$x, y=bar$y, z=bar$z, increment = 1, resamp = 200, quiet = TRUE)

p5BumpsI <- data.frame(I = zI$correlation, 
                       d = zI$mean.of.class,
                       p = zI$p) %>%
  filter(I < 1 & d <= 50) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),size=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=3,pch=21) +
  lims(x=c(0,50))+ #,y=c(-0.6,0.6)) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (m)",y="Moran's I",
       caption = "Crit value of p<0.01")


###################################################################
#
## Random bumps 
#
###################################################################

foo <- expand.grid(x=nVec,y=nVec)
foo$z <- NA

# pin corners
foo$z[foo$x == 0 & foo$y == 0] <- 0
foo$z[foo$x == 0 & foo$y == 99] <- 0
foo$z[foo$x == 99 & foo$y == 0] <- 0
foo$z[foo$x == 99 & foo$y == 99] <- 0

# fill in random zeros
ran0 <- sample(1:nrow(foo),size=1000)
foo$z[ran0] <- 0

# make a random hillS
# small
foo$z[foo$x %in% 73:77 & foo$y %in% 73:77] <- runif(25,0.5,1)
foo$z[foo$x %in% 53:57 & foo$y %in% 43:47] <- runif(25,0.5,1)

#bigger
foo$z[foo$x %in% 21:30 & foo$y %in% 33:42] <- runif(100,0.25,1)
foo$z[foo$x %in% 11:15 & foo$y %in% 13:22] <- runif(50,0.25,1)

#INTERPOLATE
fooTPS <- Tps(x = foo[,1:2], Y = foo$z,verbose=FALSE)
fooTPS <- c(predict(object=fooTPS, x = cbind(foo[,1:2])))
foo$z <- fooTPS
eps <- abs(rnorm(n^2,sd=0.1))
foo$z <- foo$z + eps

pRanBumpsMap <- ggplot(foo,aes(x=x,y=y,fill=z)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  coord_fixed() + 
  labs(x="Easting (m)", y="Northing (m)") + 
  theme_minimal() +
  theme(legend.position = "none")

fooMat <- matrix(foo$z,n,n)
pRanBumpsPersp <- plot_ly(z = fooMat) %>% 
  add_surface() %>%
  hide_colorbar() %>%
  layout(title = "",
         scene = list(
           xaxis = list(title = "Easting (m)"),
           yaxis = list(title = "Northing (m)"),
           camera = list(eye = list(x = 1.95, y = -1.25, z = 1.25))
         ))

samps2get <- sample(n^2,size = 500)
bar <- foo[samps2get,]

zI <- correlog(x=bar$x, y=bar$y, z=bar$z, increment = 1, resamp = 200, quiet = TRUE)

pRanBumpsI <- data.frame(I = zI$correlation, 
                         d = zI$mean.of.class,
                         p = zI$p) %>%
  filter(I < 1 & d <= 50) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),size=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=3,pch=21) +
  lims(x=c(0,50))+ #,y=c(-0.6,0.6)) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (m)",y="Moran's I",
       caption = "Crit value of p<0.01")


###################################################################
#
## Random bumps on a gradient
#
###################################################################

# add a gradient to foo
for(i in 1:n-1){
  foo$z[foo$x==i] <- foo$z[foo$x==i] + i*0.005
}

pGradBumpsMap <- ggplot(foo,aes(x=x,y=y,fill=z)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  coord_fixed() + 
  labs(x="Easting (m)", y="Northing (m)") + 
  theme_minimal() +
  theme(legend.position = "none")

fooMat <- matrix(foo$z,n,n)
pGradBumpsPersp <- plot_ly(z = fooMat) %>% 
  add_surface() %>%
  hide_colorbar() %>%
  layout(title = "",
         scene = list(
           xaxis = list(title = "Easting (m)"),
           yaxis = list(title = "Northing (m)"),
           camera = list(eye = list(x = 1.95, y = -1.25, z = 1.25))
         ))


samps2get <- sample(n^2,size = 500)
bar <- foo[samps2get,]
zI <- correlog(x=bar$x, y=bar$y, z=bar$z, increment = 1, resamp = 200, quiet = TRUE)

pGradBumpsI <-data.frame(I = zI$correlation, 
                         d = zI$mean.of.class,
                         p = zI$p) %>%
  filter(I < 1 & d <= 50) %>%
  mutate(Significant = p < .01) %>%
  ggplot() +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_path(aes(x = d, y = I,group = 1, color=Significant),size=1) + 
  geom_point(aes(x = d, y = I, fill=Significant),size=3,pch=21) +
  lims(x=c(0,50))+ #,y=c(-0.6,0.6)) +
  scale_fill_manual(values = c("grey","darkgreen")) +
  scale_color_manual(values = c("grey","darkgreen")) +
  labs(x="Distance (m)",y="Moran's I",
       caption = "Crit value of p<0.01")



## ----echo=FALSE---------------------------------------------------------------
p5BumpsMap


## ----echo=FALSE---------------------------------------------------------------
p5BumpsPersp


## ----echo=FALSE---------------------------------------------------------------
p5BumpsI


## ----echo=FALSE---------------------------------------------------------------
p10BumpsMap


## ----echo=FALSE---------------------------------------------------------------
p10BumpsPersp


## ----echo=FALSE---------------------------------------------------------------
p10BumpsI


## ----echo=FALSE---------------------------------------------------------------
pRanBumpsMap


## ----echo=FALSE---------------------------------------------------------------
pRanBumpsPersp


## ----echo=FALSE---------------------------------------------------------------
pRanBumpsI


## ----echo=FALSE---------------------------------------------------------------
pGradBumpsMap


## ----echo=FALSE---------------------------------------------------------------
pGradBumpsPersp


## ----echo=FALSE---------------------------------------------------------------
pGradBumpsI


## ----message=FALSE------------------------------------------------------------

## -----------------------------------------------------------------------------
# continuous
leadI <- spline.correlog(x=meuseX, y=meuseY, z=meuse_sf$log_lead, 
                         resamp=100, xmax=1500, quiet=TRUE)
plot(leadI)

