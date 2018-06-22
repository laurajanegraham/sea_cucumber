---
title: 'Sea Cucumber: spatial hurdle model'
author: "Laura Graham"
output: 
  html_document: 
    keep_md: true
---



## Aim

In this analysis we want to create spatial estimates of sea cucumber abundance within the < 30m depth region of the Gulf of California. Additionally, we want to create maps of habitat suitability and diversity for this same region. Habitat suitability is defined as any substrate type that is not sand. Diversity is the Shannon evenness of the 4 substrate types (sand, cobble, boulder, bedrock). It should be noted that because these will only be based on spatial structure in the measured data, we are making the rather restrictive assumption that the sampled points are representative of the full region. 

## Data

We have data on sea cucumber density in a number of locations across a number of years. For each of the sites we also have information about substrate type. At present we are modelling sea cucumber abundance (sc_density*100). 

### Load data


```r
# substrate data - calculate from this a % suitable and Shannon evenness for each transect
# (NB 0 = sand = not suitable; all other substrate types are suitable)
substr_dat <- read_csv("data/substrate_data.csv") %>% 
  mutate(date = dmy(date)) %>% 
  group_by(site_name, date, transect) %>%
  mutate(n_transect = max(row_number())) %>% 
  group_by(site_name, date, transect, substrate_code) %>% 
  summarise(count = n(), 
            n_transect = mean(n_transect),
            p = count / n_transect,
            plogp = p*log(p, base = exp(1))) %>% 
  mutate(p_suitable = case_when(substrate_code == 0 ~ 0, TRUE ~ p)) %>% 
  group_by(site_name, date, transect) %>% 
  summarise(shannon = (-(sum(plogp)))/log(4, base = exp(1)),
            p_suitable = sum(p_suitable)) %>% 
  group_by(site_name, transect) %>% 
  summarise(shannon = mean(shannon), 
            p_suitable = mean(p_suitable)) %>% 
  mutate(p = (p_suitable - min(p_suitable)) / (max(p_suitable) - min(p_suitable)),
         p_logit = log((p + 0.1) / (1 - p + 0.1), base = exp(1)))
```

```
## Parsed with column specification:
## cols(
##   site_name = col_character(),
##   latitude = col_double(),
##   longitude = col_double(),
##   date = col_character(),
##   transect = col_integer(),
##   substrate_code = col_integer()
## )
```

```
## Warning: 3007 failed to parse.
```

```r
# alternative substr_dat - this one a column per substrate with percentages (or perhaps counts)
substr_dat2 <- read_csv("data/substrate_data.csv") %>% 
  mutate(yday = yday(dmy(date))) %>% 
  select(-site_name, -date) %>% 
  mutate(idx = group_indices(., latitude, longitude, yday, transect))
```

```
## Parsed with column specification:
## cols(
##   site_name = col_character(),
##   latitude = col_double(),
##   longitude = col_double(),
##   date = col_character(),
##   transect = col_integer(),
##   substrate_code = col_integer()
## )
```

```
## Warning: 3007 failed to parse.
```

```r
# sea cucumber abundance and some environmental variables
dat <- read_csv("data/sea_cucumber_data.csv", na = c("", "NA", "ND")) %>% 
  # get the substrate type
  inner_join(substr_dat) %>% 
  # select the appropriate variables and omit and remaining NAs
  select(latitude, longitude, site_name, date, air_temp, diss_o2, 
         sc_density, shannon, p_suitable, p_logit) %>% 
  na.omit %>% 
  # create derived variables for day of year and y/z for the hurdle model and fix longitude (given as E)
  mutate(date = dmy(date), 
         yday = yday(date),
         longitude = -longitude,
         sc_abundance = round(sc_density*100),
         z = case_when(sc_abundance == 0 ~ 0, 
                       TRUE ~ 1),
         y = case_when(sc_abundance == 0 ~ NaN,
                       TRUE ~ sc_abundance) %>% as.integer) %>% 
  select(-date)
```

```
## Parsed with column specification:
## cols(
##   site_name = col_character(),
##   latitude = col_double(),
##   longitude = col_double(),
##   date = col_character(),
##   air_temp = col_double(),
##   wind = col_double(),
##   salinity = col_double(),
##   diss_o2 = col_double(),
##   transect = col_integer(),
##   depth = col_double(),
##   toc = col_double(),
##   sc_density = col_double(),
##   longitude2 = col_double()
## )
```

```
## Joining, by = c("site_name", "transect")
```

```r
# narrow version of dataset for plotting
dat_narrow <- dat %>% gather(covariate, value, -latitude, -longitude, -yday, -z, -y, -sc_abundance)
```

We have 724 observations of sea cucumber densities. NB to avoid complications, we will be modelling as abundances (density * 100) and then converting the predictions back to densities. 

### Data exploration

What do the sea cucumber densities look like. 


```r
ggplot(data = dat, aes(x = sc_abundance)) + 
  geom_histogram(binwidth = 1)
```

![](sea_cucumber_spatial_hurdle_files/figure-html/explore_distributions-1.png)<!-- -->

```r
dat %>% group_by(z) %>% summarise(count = n()) %>% rename(occurrence = z) %>% kable()
```



 occurrence   count
-----------  ------
          0     285
          1     439

Non-zero density ranges from 0.01 to 0.22; we will model the abundances (*100) as a negative binomial hurdle model. Occurrence will be fitted with a binomial model. 

### Spatial data

We only want to model and predict within the < 30m depth zone, so let's load this data, plot and create some prediction points. 


```r
# set the boundary for analysis
# currently fixing these boundaries so that we don't get the wide bit - revisit
bounds <- c(xmin = -114, ymin = 28, xmax = -112.5, ymax = 29.51)

# load in the depth data and crop to analysis region
depth_sf <- st_read("data/30m_isobath_baja_california.kml") %>% 
  st_cast("POLYGON") %>% st_buffer(0) %>% st_crop(bounds)
```

```
## Reading layer `Isobata Baja California' from data source `C:\Users\lg1u16\GIT_PROJECTS\SIDEPROJ\sea_cucumber\data\30m_isobath_baja_california.kml' using driver `KML'
## Simple feature collection with 33 features and 2 fields
## geometry type:  MULTILINESTRING
## dimension:      XYZ
## bbox:           xmin: -118.3721 ymin: 27.96276 xmax: -112.7483 ymax: 32.52787
## epsg (SRID):    4326
## proj4string:    +proj=longlat +datum=WGS84 +no_defs
```

```
## Warning in st_buffer.sfc(st_geometry(x), dist, nQuadSegs): st_buffer does
## not correctly buffer longitude/latitude data
```

```
## dist is assumed to be in decimal degrees (arc_degrees).
```

```
## although coordinates are longitude/latitude, st_intersection assumes that they are planar
```

```
## Warning: attribute variables are assumed to be spatially constant
## throughout all geometries
```

```r
# load in the land shapefile and crop to analysis region
map_sf <- st_read("data/bc_shp/bc_municipio.shp") %>% 
  st_crop(bounds)
```

```
## Reading layer `bc_municipio' from data source `C:\Users\lg1u16\GIT_PROJECTS\SIDEPROJ\sea_cucumber\data\bc_shp\bc_municipio.shp' using driver `ESRI Shapefile'
## Simple feature collection with 5 features and 4 fields
## geometry type:  MULTIPOLYGON
## dimension:      XY
## bbox:           xmin: -118.4076 ymin: 28 xmax: -112.6542 ymax: 32.71865
## epsg (SRID):    4326
## proj4string:    +proj=longlat +ellps=WGS84 +no_defs
```

```
## although coordinates are longitude/latitude, st_intersection assumes that they are planar
```

```
## Warning: attribute variables are assumed to be spatially constant
## throughout all geometries
```

```r
depth_sp <- as(depth_sf, "Spatial")
map_sp <- as(map_sf, "Spatial")

# create a prediction surface
ras <- raster(extent(depth_sp), crs = crs(depth_sp))

res(ras) <- 0.001
depth_ras <- rasterize(depth_sp, ras)
depth_ras <- mask(depth_ras, map_sp, inverse = T)

pred_dat <- as.data.frame(depth_ras, xy = TRUE) %>% 
  na.omit %>% select(x, y)

# just for the poster
dat_sml <- dat %>% filter(longitude >= bounds[1], 
                          longitude <= bounds[3], 
                          latitude >= bounds[2], 
                          latitude <= bounds[4])
ggplot() + 
  geom_sf(data=map_sf, fill = "grey", colour = "grey") +
  geom_point(data = dat_sml %>% arrange(sc_abundance), 
             aes(x = longitude, y = latitude, colour = sc_abundance), size = 3) +
  scale_colour_viridis_c(name = "Sea Cucumber Abundance") + 
  theme(axis.title = element_blank(),
        legend.position = "bottom")
```

![](sea_cucumber_spatial_hurdle_files/figure-html/spatial_data-1.png)<!-- -->

## INLA Modelling

We will compare 4 different modelling approaches. All approaches are hurdle models where occurrence is modelled using the binomial distribution and abundance is modelled using the negative binomial distribution: 

1. null: Intercept-only model
2. habitat: Proportion of suitable habitat (non-sand substrate) as a fixed non-spatial covariate
3. spatial: Spatial structure included as a random field (SPDE Matern correlation)
4. joint: Proportion of suitable habitat is estimated using a random spatial field and Gaussian distribution; this is then used in both the occurrence and abundance part of the hurdle model as a fixed spatial factor. This allows us to estimate both habitat suitability and abundance at the same time. 

### INLA Mesh

We need a spatial mesh to create the spatial random field. This will be restricted to the convex hull around the prediction points. 


```r
# Create mesh
pred_dat <- pred_dat %>% as.matrix
ch <- inla.nonconvex.hull(rbind(pred_dat, select(dat, x = longitude, y = latitude)) %>% as.matrix, convex = -0.05)
```

```
## Warning: 'cBind' is deprecated.
##  Since R version 3.2.0, base's cbind() should work fine with S4 objects
```

```r
# NB this mesh could be smaller - when we've got a good analysis going, need to change it. 
mesh <- inla.mesh.2d(boundary = ch,
                     #offset = 0.1,
                     max.edge = c(0.1))

# get locs of points and plot the mesh with the pred and obs points
xy <- cbind(dat$longitude, dat$latitude)
plot(mesh)
#points(pred_dat[,1:2], col = 1, pch = 16, cex = 1)
points(xy, col = 2, pch = 16, cex = 1)
```

![](sea_cucumber_spatial_hurdle_files/figure-html/create_mesh-1.png)<!-- -->

```r
mesh$n
```

```
## [1] 326
```

```r
# Associate observation locations with mesh vertices
A <- inla.spde.make.A(mesh, loc = xy)
spde <- inla.spde2.matern(mesh, alpha = 2)

A_pred <- inla.spde.make.A(mesh, loc = pred_dat)
```

### INLA Stacks

To model in INLA we need to create data stacks of fitting and prediction data. 


```r
# create the stacks
nobs = nrow(dat)

stack_y <- inla.stack(tag = "est.y", 
                      data = list(amount = dat$y, # single model
                                  alldata = cbind(dat$y, NA)), # joint model
                      A = list(A, A, 1, 1),
                      effects = list(
                        y_field = 1:spde$n.spde,
                        yc_field = 1:spde$n.spde,
                        p_suitable = dat$p_suitable,
                        y_intercept = rep(1, nobs)))
                      
stack_z <- inla.stack(tag = "est.z", 
                      data = list(occurrence = dat$z, # single model
                                  alldata = cbind(NA, dat$z)), # joint model
                      A = list(A, 1, 1),
                      effects = list(
                        z_field = 1:spde$n.spde,
                        p_suitable = dat$p_suitable,
                        z_intercept = rep(1, nobs)))

stack_yh <- inla.stack(tag = "est.y", 
                       data = list(amount = dat$y, # single model
                                   alldata = cbind(dat$y, NA, NA)), # joint model
                       A = list(A, A, A, 1),
                       effects = 
                         list(y_field = 1:spde$n.spde,
                              yc_field = 1:spde$n.spde,
                              yh_field = 1:spde$n.spde,
                              y_intercept = rep(1, nobs)))

stack_zh <- inla.stack(tag = "est.z", 
                       data = list(occurrence = dat$z, # single model
                                   alldata = cbind(NA, dat$z, NA)), # joint model
                       A = list(A, A, 1),
                       effects = list(
                         z_field = 1:spde$n.spde,
                         zh_field = 1:spde$n.spde,
                         z_intercept = rep(1, nobs)))

stack_hab <- inla.stack(tag = "est.hab",
                        data = list(p_suitable = dat$p_suitable, # single model
                                    alldata = cbind(NA, NA, dat$p_suitable)), # joint model
                        A = list(A, 1), 
                        effects = list(
                          hab_field = 1:spde$n.spde,
                          hab_intercept = rep(1, nobs)))

stack_yz <- inla.stack(stack_y, stack_z)
stack_all <- inla.stack(stack_yh, stack_zh, stack_hab)
```

### INLA Fitting

Finally, we're going to run the model. We will use a hurdle model approach (binomial for occupancy, negbin for abundance). To do this, we are following Chapter Six of the [INLA SPDE Tutorial](https://folk.ntnu.no/fuglstad/Lund2016/Session6/spde-tutorial.pdf). 

Starting off as intercept only models


```r
# 1. Null (intercept only)
f_null <- alldata ~ -1 + z_intercept + y_intercept
m_null <- inla(f_null, family = c("nbinomial", "binomial"),
                      data = inla.stack.data(stack_yz),
                      control.predictor = list(A = inla.stack.A(stack_yz)),
                      control.compute = list(dic = TRUE),
                      control.inla = list(strategy = "laplace"))

# 2. Habitat only (fixed effect)
f_habitat <- alldata ~ -1 + z_intercept + y_intercept + p_suitable

m_habitat <- inla(f_habitat, family = c("nbinomial", "binomial"),
                  data = inla.stack.data(stack_yz),
                  control.predictor = list(A = inla.stack.A(stack_yz)),
                  control.compute = list(dic = TRUE),
                  control.inla = list(strategy = "laplace"))

# 3. Space only (random spatial field)
f_spatial <- alldata ~ -1 + z_intercept + y_intercept + 
  f(z_field, model = spde) + 
  f(y_field, model = spde) + 
  f(yc_field, copy = "z_field", fixed = FALSE)

m_spatial <- inla(f_spatial, 
               family = c("nbinomial", "binomial"), 
               data = inla.stack.data(stack_yz),
               control.predictor = list(A = inla.stack.A(stack_yz)),
               control.compute = list(dic = TRUE, config = TRUE),
               control.inla = list(strategy = "laplace"))

# 4. Habitat and space (random spatial field, joint habitat model)
f_joint <- alldata ~ -1 + z_intercept + y_intercept + hab_intercept + 
  f(z_field, model = spde) + 
  f(zh_field, copy = "hab_field", fixed = TRUE) + 
  f(y_field, model = spde) + 
  f(yc_field, copy = "z_field", fixed = FALSE) + 
  f(yh_field, copy = "hab_field", fixed = TRUE) + 
  f(hab_field, model = spde)

m_joint <- inla(f_joint, 
               family = c("nbinomial", "binomial", "gaussian"), 
               data = inla.stack.data(stack_all),
               control.predictor = list(A = inla.stack.A(stack_all)),
               control.compute = list(dic = TRUE, config = TRUE),
               control.inla = list(strategy = "laplace"))

mods <- list(null = m_null, habitat = m_habitat, spatial = m_spatial, joint = m_joint)
save(mods, file = "results/model_out.Rda")
```

## Model Comparison

### DIC Comparison


```r
load("results/model_out.Rda")

# DIC table
map_dfr(mods, function(x) {
  tibble(dic = x$dic$local.dic, family = x$dic$family) %>% 
    group_by(family) %>% 
    summarise(dic = sum(dic)) %>% 
    filter(family %in% c(1, 2)) %>% 
    mutate(measure = c("Abundance", "Occurrence"))
}, .id = "model") %>% 
  select(-family) %>% 
  spread(measure, dic) %>% 
  arrange(Abundance) %>% 
  kable
```



model      Abundance   Occurrence
--------  ----------  -----------
joint       2106.861     927.3106
spatial     2109.251     925.4781
habitat     2151.477     975.4491
null        2162.974     973.2475

Joint model outperforms all others for both abundance and occurrence, based on DIC. NB not much improvement between spatial and joint (which is the spatial model that includes habitat), but the habitat part of the model requires improvement. 

### Predicted vs Observed


```r
idy <- which(dat[,"y"] > 0)
idz <- which(!is.na(dat$z))

pred_vals <- map_dfr(mods, function(x) {
  df <- tibble(predicted = x$summary.fitted.values$mean, 
             pred_sd = x$summary.fitted.values$sd)
  df[idy, "measure"] <- "Density"
  df[idy, "observed"] <- (dat$y/100) %>% na.omit
  df[idz + nrow(dat), "measure"] <- "Occurrence"
  df[idz + nrow(dat), "observed"] <- dat$z %>% na.omit
  return(df)
}, .id = "model") %>% na.omit %>% 
  mutate(model = factor(model, levels = c("null", "habitat", "spatial", "joint")),
         predicted = case_when(measure == "Density" ~ predicted / 100,
                               TRUE ~ predicted),
         pred_sd = case_when(measure == "Density" ~ pred_sd / 100,
                               TRUE ~ pred_sd))

corrs <- pred_vals %>% 
  filter(measure == "Density") %>% 
  select(model, observed, predicted) %>% 
  group_by(model) %>% 
  nest() %>% 
  mutate(corr = map_dbl(data, function(x) cor.test(x$observed, x$predicted)$estimate)) %>% 
  select(model, corr) %>% 
  unnest()

density <- ggplot(pred_vals %>% filter(measure == "Density"), aes(x = observed, y = predicted)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = predicted - pred_sd, ymax = predicted + pred_sd), width = 0) + 
  annotate("text", x = 0.05, y = 0.15, label = paste0("r = ", round(corrs$corr, 3))) +
  facet_wrap(~ model, ncol = 2) + 
  geom_smooth(method = "lm")

occurrence <- ggplot(pred_vals %>% filter(measure == "Occurrence"), aes(x = as.factor(observed), y = predicted)) + 
  geom_boxplot() + 
  facet_wrap(~ model, ncol = 2) + 
  labs(x = "observed")

plot_grid(density, occurrence)
```

![](sea_cucumber_spatial_hurdle_files/figure-html/pred_obs-1.png)<!-- -->

```r
save_plot(density, filename = "figures/density_performance.jpg")
save_plot(density, filename = "figures/occurrence_performance.jpg")
```

Both spatial and joint outperform the null and habitat only models, with the joint (habitat and spatial) model performing best (again, only marginally). No models are predicting the full range of abundance, and the high predicted values carry a lot of uncertainty. 

## Prediction of the response

All predictions use the joint model.


```r
# get the id's for each field
samples <- inla.posterior.sample(n = 1000, result = mods$joint)

ids <- lapply(c("y_intercept", "y_field", "yc_field", "yh_field",
                "z_intercept", "z_field", "zh_field", 
                "hab_field", "hab_intercept"),
              function(x) grep(x, rownames(samples[[1]]$latent), fixed = TRUE))

predict_y <- function(s) exp(s$latent[ids[[1]], 1] + 
                               s$latent[ids[[2]], 1] + 
                               s$latent[ids[[3]], 1] + 
                               s$latent[ids[[4]], 1])

predict_z <- function(s) 1/ (1 + exp(-(s$latent[ids[[5]], 1] + 
                                        s$latent[ids[[6]], 1] + 
                                         s$latent[ids[[7]], 1])))

predict_hab <- function(s) s$latent[ids[[8]], 1] + s$latent[ids[[9]], 1]

pred_y <- sapply(samples, predict_y)
pred_z <- sapply(samples, predict_z)
pred_hab <- sapply(samples, predict_hab)

projgrid <- inla.mesh.projector(mesh, loc = pred_dat)

pred_vals <- as.tibble(pred_dat) %>% 
  mutate(ymean = inla.mesh.project(projgrid, field = rowMeans(pred_y)),
         ysd = inla.mesh.project(projgrid, field = apply(pred_y, 1, sd)),
         zmean = inla.mesh.project(projgrid, field = rowMeans(pred_z)),
         zsd = inla.mesh.project(projgrid, field = apply(pred_z, 1, sd)),
         habmean = inla.mesh.project(projgrid, field = rowMeans(pred_hab)),
         habsd = inla.mesh.project(projgrid, field = apply(pred_hab, 1, sd)))
```

### Abundance predictions


```r
ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") + 
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = ymean)) + 
  scale_fill_viridis_c(name = "Abundance (mean)") + 
  labs(x = "", y = "") + 
  theme(legend.position = "bottom")
```

![](sea_cucumber_spatial_hurdle_files/figure-html/abund_plot-1.png)<!-- -->

```r
ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") +
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = ysd/100)) + 
  scale_fill_viridis_c(name = "Abundance (SD)") + 
  labs(x = "", y = "") + 
  theme(legend.position = "bottom")
```

![](sea_cucumber_spatial_hurdle_files/figure-html/abund_plot-2.png)<!-- -->

### Occurrence predictions


```r
zmean_plot <- ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") +
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = zmean)) + 
  scale_fill_viridis_c(name = "Probability of\noccurrence\n(mean)") + 
  labs(x = "", y = "")

zsd_plot <- ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") +
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = zsd)) + 
  scale_fill_viridis_c(name = "Probability of\noccurrence\n(SD)") + 
  labs(x = "", y = "")

plot_grid(zmean_plot, zsd_plot)
```

![](sea_cucumber_spatial_hurdle_files/figure-html/occurrence_plot-1.png)<!-- -->

```r
save_plot(plot_grid(zmean_plot, zsd_plot), 
          filename = "figures/occurrence.jpg", 
          base_width = 15, 
          base_height = 10)
```

### Habitat proportion prediction


```r
ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") +
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = habmean)) + 
  scale_fill_viridis_c(name = "Habitat proportion (mean)") + 
  labs(x = "", y = "") + 
  theme(legend.position = "bottom")
```

![](sea_cucumber_spatial_hurdle_files/figure-html/habitat_plot-1.png)<!-- -->

```r
ggplot() + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") +
  geom_raster(data = pred_vals, aes(x = x, y = y, fill = habsd)) + 
  scale_fill_viridis_c(name = "Habitat proportion (SD)") + 
  labs(x = "", y = "") + 
  theme(legend.position = "bottom")
```

![](sea_cucumber_spatial_hurdle_files/figure-html/habitat_plot-2.png)<!-- -->

## To do

1. Model habitat as multinomial: see [INLA FAQ](http://www.r-inla.org/faq#TOC-I-have-multinomial-data-but-I-cannot-find-the-multinomial-likelihood-isn-t-it-available-) on MN distribution. 

2. Look into spatio-temporal mesh. Reasoning behind this is that date is a key factor affecting density, so to make more accurate predictions, need specific mesh for specific times of year. Speak to Luis about what would be most useful for this. 

3. Expand to the full study site (currently cutting off some points for ease of plotting for poster)