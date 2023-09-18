#' ---
#' title: "AirPollution"
#' author: "Lenka Sefcakova"
#' date: "2023-02-10"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

#' 
#' ## London Air Pollution and population density
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# Load libraries and import data
library(rgeoboundaries)
library(ggplot2)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)
library(terra)
library(gstat)
library(INLA)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#load the data
filtered_data <- read_csv(file = "data/filtered_data.csv")
#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# Summarize the data - compute mean for each location in case of multiple observations 
all_sites <- filtered_data %>%
  group_by(site) %>%
  summarize(mean = mean(pm10, na.rm = TRUE), latitude = first(latitude), longitude = first(longitude), site_type = first(site_type)) %>% drop_na()
head(all_sites)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# Download the boundaries for UK 
uk_boundary <- geoboundaries(country = "GBR")

# Plot the map with observations and their respective values
png(file="obs.png",width=600, height=350)
ggplot(data = uk_boundary) +
  geom_sf() +
  geom_point(data = all_sites, aes(x = longitude, y = latitude, color = mean), size = 3) +
  scale_color_gradient(name = "Level of Air pollution PM10", low = "blue", high = "red") +
  ggtitle("Mean level of air pollution in the UK in 2020") +
  labs(x = "Longitude", y = "Latitude")
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
c(max(all_sites$mean),min(all_sites$mean))

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# make a grid of points for interpolatioin
grid <- st_make_grid(uk_boundary, what = "centers", n = c(100, 100)) %>%  st_sf()
coop <- st_filter(grid,uk_boundary)
#plot(coop)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
png(file="histogram.png",width=600, height=350)
hist(all_sites$mean)
dev.off()


#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#create a sf object with all observations as points 
sites <-st_as_sf(all_sites,coords = c("longitude","latitude"))
#sites

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#we compute predictions for voronoi,knn and inverse distance weighting saving the results to later apply them in ensemble and show ensemble has best performance when measuring RMSE

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# Compute Voronoi and plot the results
v <- terra::voronoi(x = terra::vect(sites), bnd = uk_boundary)
v <- st_as_sf(v)
st_crs(v) = "EPSG::4258"
st_crs(coop) = "EPSG::4258"
p <- st_intersection(v, coop)
p1 <- p$mean

p$x <- st_coordinates(p)[,1]
p$y <- st_coordinates(p)[,2]
png(file="voronoi.png",width=500, height=500)
ggplot(p, aes(x, y)) + geom_tile(aes(fill = mean)) + scale_fill_viridis(option = "C") +ggtitle('Voronoi model') +labs(fill = "PM10")+ coord_fixed()
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#K nearest neighbors max neighbors 10, uniform weighting of neighbors 
res <- gstat(formula = mean ~ 1, locations = sites, nmax = 25, set = list(idp = 0))

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#obtain predictions
resp <- predict(res, coop)
p2 <- resp$var1.pred
resp$x <- st_coordinates(resp)[,1]
resp$y <- st_coordinates(resp)[,2]

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#plot nearest neighbors predictions 
png(file="knn.png",width=500, height=500)
ggplot(resp, aes(x, y)) + geom_tile(aes(fill = var1.pred)) + scale_fill_viridis(option = "C") + ggtitle('KNN model') +labs(fill = "PM10")+ coord_fixed()
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#Inverse Distance Weighting, using all points for prediction with less information coming from distant observations
library(gstat)
res <- gstat(formula = mean ~ 1, locations = sites,
             nmax = nrow(sites), # use all the neighbors locations
             set = list(idp = 7)) # beta = 1 

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
resp <- predict(res, coop)
p3 <- resp$var1.pred
resp$x <- st_coordinates(resp)[,1]
resp$y <- st_coordinates(resp)[,2]

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#Inverse Distance Weighting plot
png(file="idw.png",width=500, height=500)
ggplot(resp, aes(x, y)) + geom_tile(aes(fill = var1.pred)) + scale_fill_viridis(option = "C")+ggtitle('Inverse Distance Weighting model') +labs(fill = "PM10") + coord_fixed()
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#apply uniform wieghts on ensembles 
weights <- c(1/4, 3/8, 3/8)
p4 <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
# Plot ensemble method 
resp <- data.frame(
x = st_coordinates(coop)[, 1],
y = st_coordinates(coop)[, 2],
value = p4)

png(file="ensemble.png",width=500, height=500)
ggplot(resp, aes(x, y)) + geom_tile(aes(fill = value)) + scale_fill_viridis(option = "C") +ggtitle('Ensemble model') +labs(fill = "PM10")+ coord_fixed()
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)

# Function to calculate the RMSE
RMSE <- function(observed, predicted) {
sqrt(mean((observed - predicted)^2))
}

# Split data in 10 sets
kf <- dismo::kfold(nrow(sites), k = 20) # k-fold partitioning

# Vectors where RMSE values obtained with each of the methods will be stored
rmse1 <- rep(NA, 20) # Voronoi
rmse2 <- rep(NA, 20) # Nearest neighbors
rmse3 <- rep(NA, 20) # IDW
rmse4 <- rep(NA, 20) # Ensemble

for(k in 1:20) {
# Split data in test and train
test <- sites[kf == k, ]
train <- sites[kf != k, ]
# Voronoi
v <- terra::voronoi(x = terra::vect(train), bnd = uk_boundary)
v <- st_as_sf(v)
p1 <- st_intersection(v, test)$mean
rmse1[k] <- RMSE(test$mean, p1)
# Nearest neighbors
nn <- gstat(formula = mean ~ 1, locations = train, nmax = 10, set = list(idp = 0))
p2 <- predict(nn, test)$var1.pred
rmse2[k] <- RMSE(test$mean, p2)
# IDW
gs <- gstat(formula = mean ~ 1, locations = train,set = list(idp = 1))
p3 <- predict(gs, test)$var1.pred
rmse3[k] <- RMSE(test$mean, p3)
# Ensemble (weights are inverse RMSE so lower RMSE has higher weight)
w <- 1/c(rmse1[k], rmse2[k], rmse3[k])
weights <- w/sum(w)
p4 <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
rmse4[k] <- RMSE(test$mean, p4)
}

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#data.frame(voronoi = rmse1, near.neigh = rmse2, IDW = rmse3, ensemble = rmse4)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
library(xtable)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#data.frame(voronoi = mean(rmse1), near.neigh = mean(rmse2), IDW = mean(rmse3), ensemble = mean(rmse4))


#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
coo <- cbind(all_sites$longitude, all_sites$latitude)
bnd <- inla.nonconvex.hull(coo)
mesh <- inla.mesh.2d(boundary = bnd,loc = coo, max.edge = c(0.9, 10), cutoff = 0.3) #create a mesh to apply inla on with specified values 

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#mesh$n

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#plot(mesh)
#points(coo, col = "red")

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
spde <- inla.spde2.matern(mesh = mesh, alpha = 2) 

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
indexs <- inla.spde.make.index("s", spde$n.spde)
#lengths(indexs)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
A <- inla.spde.make.A(mesh = mesh, loc = coo) #create projection matrix A of site observations

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
library(raster)
r <- getData(name = 'alt', country = 'GBR', mask = TRUE) #retrieve raster points of UK 

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
ra <- aggregate(r, fact = 5, fun = mean) 

dp <- rasterToPoints(ra)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
coop <- dp[, c("x", "y")] #create grid with points from retrieved raster data of uk

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)#create projection matrix Ap this thime over our grid coop

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#plot(coop, asp = 1)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = all_sites$mean),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)
# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
formula <- y ~ 0 + b0 + f(s, model = spde) #define model formula with no intercept and b0 as fixed effects

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
res <- inla(formula, #call inla 
  data = inla.stack.data(stk.full),
  control.predictor = list(
    compute = TRUE,
    A = inla.stack.A(stk.full)
  )
)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
index <- inla.stack.index(stk.full, tag = "pred")$data 

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
dpm <- rbind( #create maps of the predicted values
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_mean, variable = "pred_mean"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ll, variable = "pred_ll"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ul, variable = "pred_ul"
  )
)
dpm$variable <- as.factor(dpm$variable)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
png(file="inla.png",width=1200, height=300) #save resulted plot
ggplot(dpm) + geom_tile(aes(east, north, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) + scale_fill_viridis(option = "C",name="PM10")+
  theme_bw()+ggtitle('INLA predictions')
dev.off()

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------
#summary(res)

#' 
## ---------------------------------------------------------------------------------------------------------------------------------------


