## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)


## ----message=FALSE------------------------------------------------------------
library(terra)
library(sf)
library(tidyverse)
library(ggnewscale)
library(tidyterra)


## -----------------------------------------------------------------------------
mtbDEM <- rast("../data/mtbDEM.tif")
mtbDEM
# and a quick plot
mtbDEM_df <- as.data.frame(mtbDEM, xy = TRUE)
ggplot() + 
  geom_raster(data = mtbDEM_df , mapping = aes(x = x, y = y,
                                               fill = elev)) + 
  scale_fill_gradientn(colours = terrain.colors(100)) +
  labs(x="Easting (m)", y = "Northing (m)") +
  coord_equal()


## -----------------------------------------------------------------------------
ggplot() +
  geom_spatraster(data = mtbDEM) +
  scale_fill_hypso_c(palette = "usgs-gswa2")


## -----------------------------------------------------------------------------
ggplot() +
  geom_spatraster(data = mtbDEM) +
  scale_fill_hypso_c(palette = "usgs-gswa2") +
  coord_sf(datum = 32610)


## -----------------------------------------------------------------------------
mtbSlope <- terrain(mtbDEM, "slope", unit="radians")
mtbAspect <- terrain(mtbDEM, "aspect", unit="radians")
mtbHillshade <- shade(mtbSlope, mtbAspect, angle = 40, direction = 270)
mtbHillshade[mtbHillshade < 0] <- 0
names(mtbHillshade) <- "value" # these are unitless
mtbHillshade
# and a quick plot
ggplot() +
  geom_spatraster(data = mtbHillshade)


## -----------------------------------------------------------------------------
p1 <- ggplot() +
  geom_spatraster(data = mtbHillshade) +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0.1,
                                            end = 0.9), guide = "none") +
  new_scale_fill() +
  geom_spatraster(data = mtbDEM) +
  scale_fill_hypso_c(name = "Elevation (m)", 
                     palette = "usgs-gswa2",alpha = 0.6) +
  theme_minimal()
p1


## -----------------------------------------------------------------------------
chairs <- st_read("../data/mtbChairLines.shp")
p1 + geom_sf(data=chairs)


## -----------------------------------------------------------------------------
buildings <- read.csv("../data/mtbLodges.csv")
buildings


## ----echo=FALSE---------------------------------------------------------------
# buildings
buildings <- buildings[,1:3] # guard against extra col
buildings <- st_as_sf(buildings, coords = c("X", "Y"), crs = st_crs(chairs))

# make midpoints for lines
st_line_midpoints <- function(sf_lines = NULL) {
  
  g <- st_geometry(sf_lines)
  
  g_mids <- lapply(g, function(x) {
    
    coords <- as.matrix(x)
    
    # this is just a copypaste of View(maptools:::getMidpoints):
    get_mids <- function (coords) {
      dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid <- sum(dist)/2
      dist_cum <- c(0, cumsum(dist))
      end_index <- which(dist_cum > dist_mid)[1]
      start_index <- end_index - 1
      start <- coords[start_index, ]
      end <- coords[end_index, ]
      dist_remaining <- dist_mid - dist_cum[start_index]
      mid <- start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    
    mids <- st_point(get_mids(coords))
  })
  
  out <- st_sfc(g_mids, crs = st_crs(sf_lines))
  out <- st_sf(out)
}
chairMids <- st_line_midpoints(chairs)
names(chairMids) <- "geometry"
st_geometry(chairMids) <- "geometry"
chairMids$id <- c("Chair 1",
                  "Chair 2",
                  "Chair 3",
                  "Chair 6",
                  "Chair 4",
                  "Chair 5",
                  "Chair 7",
                  "Chair 8")

# sneaky point object
points4labs <- rbind(chairMids,buildings)
points4labs$type <- c(rep("Chair",8),rep("Building",3))


## ----echo=FALSE,fig.height = 8, fig.width = 8, warning=FALSE------------------
#devtools::install_github("yutannihilation/ggsflabel")
ggplot() +
  geom_spatraster(data = mtbHillshade) +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0.1,
                                            end = 0.9), guide = "none") +
  new_scale_fill() +
  geom_spatraster(data = mtbDEM) +
  scale_fill_hypso_c(name = "Elevation (m)", 
                     palette = "wiki-2.0_hypso",alpha = 0.6) +
  geom_spatraster_contour(data=mtbDEM,breaks=c(seq(1e3,1.6e3,by=100))) +
  geom_sf(data=buildings,shape=20,size=3) + 
  geom_sf(data=chairs,color="blue") + 
  ggsflabel::geom_sf_label_repel(data=points4labs,
                                 aes(label = id,color=type),
                                 alpha=0.8) +
  scale_color_manual(values=c("black","blue")) +
  labs(x=NULL, y=NULL,
       title = "Mt Baker Ski Area",
       caption = "Elevation contours in 100 m intervals") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid=element_blank(),legend.position = "none",
        text = element_text(size = 14))

