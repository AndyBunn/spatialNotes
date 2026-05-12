## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
set.seed(184)


## -----------------------------------------------------------------------------
#| warning: false
#| message: false
library(fields)
library(tidyverse)
library(gstat)
library(sf)
library(terra)
library(tidyterra)


## -----------------------------------------------------------------------------
#| warning: false
#| message: false
data(meuse.all)
meuse.all$logLead <- log(meuse.all$lead)
meuse_sf <- st_as_sf(meuse.all, coords = c("x", "y")) %>%
  st_set_crs(value = 28992)

meuse.grid <- readRDS("../data/meuse.grid.Rds")


## -----------------------------------------------------------------------------
ggplot(data = meuse_sf) +
  geom_sf(aes(fill=logLead), size=4,
          shape = 21, color="white",alpha=0.8)+
  scale_fill_continuous(type = "viridis",name="log(ppm)") +
  labs(title="Lead concentrations")


## -----------------------------------------------------------------------------
logLeadTPSmodel <- Tps(x = meuse.all[,2:3], Y = meuse.all$logLead)
logLeadTPSmodel


## -----------------------------------------------------------------------------
# Predict the model over all the coordinates in meuse.grid
logLeadPreds <- c(predict(object=logLeadTPSmodel, x = meuse.grid[,1:2]))
# Store in a data.frame with the x,y coordinates
logLeadTPS <- data.frame(x = meuse.grid[,1],
                         y = meuse.grid[,2],
                         logLead=logLeadPreds)
# And into SpatRaster
logLeadTPS_rast <- rast(logLeadTPS, crs=crs(meuse_sf))

# Plot
ggplot() +
  geom_spatraster(data=logLeadTPS_rast, mapping = aes(fill=logLead),alpha=0.8) +
  scale_fill_continuous(type = "viridis",name="log(ppm)",na.value = "transparent") +
  labs(title="Lead concentrations", subtitle = "TPS") +
  theme_minimal()


## -----------------------------------------------------------------------------
obs <- meuse.all$logLead
preds <- extract(logLeadTPS_rast, meuse_sf) %>% pull(logLead)
rsq <- cor(obs,preds)^2
rmse <- sqrt(mean((preds - obs)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
ggplot() +
  geom_abline(slope=1,intercept = 0) +
  geom_point(aes(x=obs,y=preds)) +
  coord_fixed(ratio=1, xlim = range(preds,obs),ylim = range(preds,obs)) +
  labs(x="Observed Values",
       y="Predicted Values",
       title="Lead log(ppm)")


## -----------------------------------------------------------------------------
n <- nrow(meuse.all)
rows4test <- sample(x = 1:n,size = n*0.2)
meuseTest <- meuse.all[rows4test,]
meuseTrain <- meuse.all[-rows4test,]

# Note meuseTrain here
logLeadTPSmodel <- Tps(x = meuseTrain[,2:3], Y = meuseTrain$logLead)
logLeadTPSmodel

# Predict the model over all the coordinates in meuse.grid
logLeadPreds <- c(predict(object=logLeadTPSmodel, x = meuse.grid[,1:2]))
# Store in a data.frame with the x,y coordinates
logLeadTPS <- data.frame(x = meuse.grid[,1],
                         y = meuse.grid[,2],
                         logLead=logLeadPreds)
# And into SpatRaster
logLeadTPS_rast <- rast(logLeadTPS, crs=crs(meuse_sf))

# Look at skill on withheld data. Note meuseTest:
obs <- meuseTest$logLead
preds <- extract(logLeadTPS_rast, meuseTest[,2:3]) %>% pull(logLead)
rsq <- cor(obs,preds)^2
rmse <- sqrt(mean((preds - obs)^2))
c(rsq = rsq, rmse = rmse)


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
1 - (rmse / rmseNULL)

