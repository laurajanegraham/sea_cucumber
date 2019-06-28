---
title: 'Sea Cucumber: spatial hurdle model'
author: "Laura Graham"
output: 
  html_document: 
    df_print: kable
    highlight: tango
    keep_md: yes
    theme: flatly
params: 
  range0: 100
---



## Aim

In this analysis we want to create spatial estimates of sea cucumber abundance within the < 30m depth region of the Gulf of California. 

## Data


```r
# sea cucumber abundance and some environmental variables
dat <- read_csv("data/sea_cucumber_data.csv", na = c("", "NA", "ND")) %>% 
  select(site_name, transect, latitude, longitude, date, sc_density) %>% 
  mutate(longitude = -longitude,
         transect = paste0(site_name, transect),
         date = dmy(date),
         month = str_sub(date, 1, 7),
         z = case_when(sc_density == 0 ~ 0, 
                       TRUE ~ 1),
         y = case_when(sc_density == 0 ~ NaN,
                       TRUE ~ sc_density)) %>% 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) %>% 
  st_transform(crs_ea_mex)
```

We have 1107 observations of sea cucumber densities. These are from 118 sites, collected between 2014-10-16 and 2015-09-12.

### Data exploration

What do the sea cucumber densities look like. 


```r
ggplot(data = dat, aes(x = sc_density)) + 
  geom_histogram(binwidth = 0.01)
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/explore_distributions-1.png)<!-- -->

Non-zero density ranges from 10^{-4} to 0.0058; we will model the abundances (*100) as a negative binomial hurdle model. Occurrence will be fitted with a binomial model. 

### Spatial data

We only want to model and predict within the < 30m depth zone, so let's load this data, plot and create some prediction points. 


```r
# set the boundary for analysis
# currently fixing these boundaries so that we don't get the wide bit - revisit
bounds <- c(xmin = -114.5, ymin = 28, xmax = -112.5, ymax = 30)

# load in the depth data and crop to analysis region
depth_sf <- st_read("data/30m_isobath_baja_california.kml", quiet = TRUE) %>% 
  st_crop(bounds) %>% 
  st_cast("POLYGON") %>% 
  st_zm() %>% 
  st_transform(crs_ea_mex) %>% 
  st_buffer(0) %>% 
  filter(Name != 0) %>% 
  st_cast("MULTIPOLYGON")

# load in the land shapefile and crop to analysis region
map_sf <- st_read("data/bc_shp/bc_municipio.shp", quiet = TRUE) %>% 
  st_crop(bounds) %>% 
  st_transform(crs_ea_mex)

# remove the intersection of land and depth
depth_sf <- st_difference(depth_sf, st_combine(map_sf)) %>% st_cast("MULTIPOLYGON")

# get the estates and quotas
f <- list.files("data/quota_shp/", pattern = ".shp", full.names = TRUE)

quotas <- st_read("data/quota_estates.shp", quiet = TRUE)

ggplot() + geom_sf(data = quotas, aes(fill = Estate))
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/spatial_data-1.png)<!-- -->

```r
quotas_sp <- as(quotas, "Spatial")

study_site <- ggplot() + 
  geom_sf(data=map_sf, fill = "grey", colour = "grey") +
  geom_sf(data = quotas) +
  geom_sf(data = dat %>% arrange(sc_density), 
             aes(colour = sc_density)) +
  scale_colour_viridis_c(name = "Sea Cucumber\nDensity") + 
  theme(axis.title = element_blank(),
        panel.grid.major=element_line(colour="transparent")) 

ggsave(plot = study_site, filename = "figures/study_site.jpg", dpi = 300, width = 140, units = "mm")
```


```r
ggplot() +
  geom_sf(data = dat %>% select(-month), colour = "lightgrey") + 
  geom_sf(data = dat %>% arrange(sc_density), aes(colour = sc_density)) + 
  scale_colour_viridis_c(name = "Sea Cucumber Abundance") + 
  facet_wrap(~month)
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/spatiotemporal_data-1.png)<!-- -->

Date structure in the sampling, will need to account for. 

## INLA Modelling

We will compare 2 different modelling approaches: one with a spatial random field (spatial) and one without (null). Both approaches are hurdle models where occurrence is modelled using the binomial distribution and density is modelled using the gamma distribution. In both cases, we include the following additional random effects: site (iid), site:transect (iid), date (rw2). 

### INLA Mesh

We need a spatial mesh to create the spatial random field. This will be restricted to the convex hull around the prediction points, and such that land forms a barrier. The range value for the PC prior is $\rho =$ 500 and the sd value is `\sigma = 30`. Both have low probability (0.05). 


```r
# get boundary
bound <- dat %>% 
  st_union %>% 
  st_convex_hull %>% 
  st_buffer(dist = 20000) %>% 
  st_difference(st_combine(map_sf))

bound_sp <- as(bound, "Spatial")
bound_sp <- gSimplify(bound_sp, tol = 100)

# Create mesh
coords_dat <- st_coordinates(dat)
mesh <- inla.mesh.2d(boundary = as.inla.mesh.segment(bound_sp),
                     #offset = 0.1,
                     max.edge = 5000)

max_edge = 10000
bound_outer = 50000
mesh = inla.mesh.2d(boundary = bound_sp,
                     loc=coords_dat,
                    max.edge = c(1,5)*max_edge,
                    cutoff = max_edge/10,
                    offset = c(max_edge, bound_outer))

# get locs of points and plot the mesh with the pred and obs points
plot(mesh, lwd = 0.5)
points(coords_dat, col = "black", pch = 16, cex = 1)
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/create_mesh-1.png)<!-- -->

```r
mesh$n
```

```
## [1] 1399
```

```r
# get the triangles which are in/outside of the boundary
tl = length(mesh$graph$tv[,1])

pos_tri <- sapply(1:tl, function(t) {
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  return(colMeans(temp)[c(1,2)])
})

pos_tri = SpatialPoints(t(pos_tri), proj4string = crs(bound_sp))

# - compute the triangle positions
normal = over(bound_sp, pos_tri, returnList=T)
# - checking which mesh triangles are inside the normal area
normal = unlist(normal)
barrier_triangles = setdiff(1:tl, normal)

poly_barrier = inla.barrier.polygon(mesh, barrier_triangles)

# Associate observation locations with mesh vertices
A <- inla.spde.make.A(mesh, loc = coords_dat)

spde <- inla.barrier.pcmatern(mesh, 
                              barrier.triangles = barrier_triangles, 
                              prior.range = c(range0, 0.05), 
                              prior.sigma = c(30, 0.05))
```

### INLA Stacks

To model in INLA we need to create data stacks of fitting and prediction data. 


```r
# create the stacks
nobs = nrow(dat)
dat <- dat %>% mutate(site_name = as.numeric(as.factor(site_name)),
                      transect = as.numeric(as.factor(transect)),
                      date = as.numeric(date))

stack_y <- inla.stack(tag = "est.y", 
                      data = list(alldata = cbind(dat$y, NA), link = 1), 
                      A = list(A, A, 1, 1, 1, 1),
                      effects = list(
                        y_field = 1:mesh$n,
                        yc_field = 1:mesh$n,
                        y_intercept = rep(1, nobs),
                        site = dat$site_name,
                        transect = dat$transect,
                        date = dat$date))

stack_z <- inla.stack(tag = "est.z", 
                      data = list(alldata = cbind(NA, dat$z), link = 2), 
                      A = list(A, 1, 1, 1, 1),
                      effects = list(
                        z_field = 1:mesh$n,
                        z_intercept = rep(1, nobs),
                        site = dat$site_name,
                        transect = dat$transect,
                        date = dat$date))

stack_yz <- inla.stack(stack_y, stack_z)
```

### INLA Fitting

Finally, we're going to run the model. We will use a hurdle model approach (binomial for occupancy, negbin for abundance). To do this, we are following Chapter Six of the [INLA SPDE Tutorial](https://folk.ntnu.no/fuglstad/Lund2016/Session6/spde-tutorial.pdf). 


```r
f_null <- alldata ~ -1 + z_intercept + y_intercept + 
  f(date, model = "rw2") + 
  f(transect, model = "iid") + 
  f(site, model = "iid")

m_null <- inla(f_null, family = c("gamma", "binomial"),
                data = inla.stack.data(stack_yz),
                control.predictor = list(A = inla.stack.A(stack_yz)),
                control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))

# 3. Space only (random spatial field)
f_spatial <- alldata ~ -1 + z_intercept + y_intercept + 
  f(date, model = "rw2") +
  f(site, model = "iid") +
  f(transect, model = "iid") + 
  f(z_field, model = spde) +
  f(y_field, model = spde) + 
  f(yc_field, copy = "z_field", fixed = FALSE)

m_spatial <- inla(f_spatial, family = c("gamma", "binomial"),
                  data = inla.stack.data(stack_yz),
                  control.predictor = list(A = inla.stack.A(stack_yz)),
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))

mods <- list(null = m_null,
             spatial = m_spatial)

save(mods, file = fname)
```

## Model Comparison

### DIC/WAIC Comparison


```r
load(fname)

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
spatial    -3362.256     1289.612
null       -3312.765     1323.099

```r
# WAIC table
map_dfr(mods, function(x) {
  tibble(waic = x$waic$local.waic, family = x$dic$family) %>% 
    group_by(family) %>% 
    summarise(waic = sum(waic)) %>% 
    filter(family %in% c(1, 2)) %>% 
    mutate(measure = c("Abundance", "Occurrence"))
}, .id = "model") %>% 
  select(-family) %>% 
  spread(measure, waic) %>% 
  arrange(Abundance) %>% 
  kable
```



model      Abundance   Occurrence
--------  ----------  -----------
spatial    -3381.616     1288.858
null       -3326.761     1324.217

### Predicted vs Observed


```r
dat <- dat %>% st_set_geometry(NULL)
idy <- which(dat[,"y"] > 0)
idz <- which(!is.na(dat$z))

pred_vals <- map_dfr(mods, function(x) {
  df <- tibble(predicted = x$summary.fitted.values$mean, 
             pred_sd = x$summary.fitted.values$sd)
  df[idy, "measure"] <- "Density"
  df[idy, "observed"] <- (dat$y) %>% na.omit
  df[idz + nrow(dat), "measure"] <- "Occurrence"
  df[idz + nrow(dat), "observed"] <- dat$z %>% na.omit
  return(df)
}, .id = "model") %>% na.omit

corrs <- pred_vals %>% 
  filter(measure == "Density") %>% 
  select(model, observed, predicted) %>% 
  group_by(model) %>% 
  nest() %>% 
  mutate(corr = map_dbl(data, function(x) cor.test(x$observed, x$predicted)$estimate)) %>% 
  select(model, corr) %>% 
  unnest()

density <- ggplot(pred_vals %>% filter(measure == "Density"), 
                  aes(x = observed, y = predicted)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = predicted - pred_sd, ymax = predicted + pred_sd), width = 0) + 
  annotate("text", x = 0.3, y = 0.4, label = paste0("r = ", round(corrs$corr, 3))) +
  facet_wrap(~ model, nrow = 1) + 
  geom_smooth(method = "lm")

occurrence <- ggplot(pred_vals %>% filter(measure == "Occurrence"), 
                     aes(x = as.factor(observed), y = predicted)) + 
  geom_boxplot() + 
  facet_wrap(~ model, nrow = 1) + 
  labs(x = "observed")

density + occurrence + plot_layout(nrow = 2)
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/pred_obs-1.png)<!-- -->

### Model summary


```r
summary(mods$spatial)
```

```
## 
## Call:
##    c("inla(formula = f_spatial, family = c(\"gamma\", \"binomial\"), 
##    data = inla.stack.data(stack_yz), ", " control.compute = list(dic 
##    = TRUE, waic = TRUE, config = TRUE), ", " control.predictor = 
##    list(A = inla.stack.A(stack_yz)))") 
## Time used:
##     Pre = 1.42, Running = 920, Post = 4.53, Total = 926 
## Fixed effects:
##               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
## z_intercept  0.362 0.165      0.027    0.364      0.682  0.369   0
## y_intercept -3.548 0.096     -3.747   -3.544     -3.368 -3.538   0
## 
## Random effects:
##   Name	  Model
##     date RW2 model
##    site IID model
##    transect IID model
##    z_field RGeneric2
##    y_field RGeneric2
##    yc_field Copy
## 
## Model hyperparameters:
##                                                     mean       sd
## Precision parameter for the Gamma observations     5.927 5.56e-01
## Precision for date                             26553.492 1.45e+04
## Precision for site                                 8.660 2.81e+00
## Precision for transect                             4.771 7.18e-01
## Theta1 for z_field                                -0.142 1.78e-01
## Theta2 for z_field                                 8.342 4.30e-01
## Theta1 for y_field                                -1.618 5.43e-01
## Theta2 for y_field                                 6.746 6.66e-01
## Beta for yc_field                                  0.503 8.70e-02
##                                                0.025quant  0.5quant
## Precision parameter for the Gamma observations      4.895     5.907
## Precision for date                               7486.545 23663.710
## Precision for site                                  4.563     8.179
## Precision for transect                              3.537     4.709
## Theta1 for z_field                                 -0.499    -0.138
## Theta2 for z_field                                  7.494     8.342
## Theta1 for y_field                                 -2.777    -1.575
## Theta2 for y_field                                  5.385     6.769
## Beta for yc_field                                   0.335     0.501
##                                                0.975quant      mode
## Precision parameter for the Gamma observations      7.081     5.874
## Precision for date                              62987.772 17853.808
## Precision for site                                 15.469     7.313
## Precision for transect                              6.354     4.579
## Theta1 for z_field                                  0.199    -0.126
## Theta2 for z_field                                  9.187     8.342
## Theta1 for y_field                                 -0.664    -1.419
## Theta2 for y_field                                  7.999     6.851
## Beta for yc_field                                   0.677     0.496
## 
## Expected number of effective parameters(stdev): 362.98(20.95)
## Number of equivalent replicates : 4.90 
## 
## Deviance Information Criterion (DIC) ...............: -2072.64
## Deviance Information Criterion (DIC, saturated) ....: 2292.26
## Effective number of parameters .....................: 362.76
## 
## Watanabe-Akaike information criterion (WAIC) ...: -2092.76
## Effective number of parameters .................: 274.73
## 
## Marginal log-Likelihood:  876.99 
## Posterior marginals for the linear predictor and
##  the fitted values are computed
```

## Prediction of the response


```r
# get the id's for each field
samples <- inla.posterior.sample(n = 10000, result = mods$spatial)

ids <- lapply(c("y_intercept", "y_field", "yc_field",
                "z_intercept", "z_field"),
              function(x) grep(x, rownames(samples[[1]]$latent), fixed = TRUE))

predict_y <- function(s) exp(s$latent[ids[[1]], 1] + 
                               s$latent[ids[[2]], 1] + 
                               s$latent[ids[[3]], 1])

predict_z <- function(s) 1/ (1 + exp(-(s$latent[ids[[4]], 1] + 
                                        s$latent[ids[[5]], 1])))

pred_y <- sapply(samples, predict_y)
pred_z <- sapply(samples, predict_z)

r <- raster(st_sf(bound), res = 100)
out <- fasterize::fasterize(st_sf(bound), r)
locs <- as.data.frame(out, xy = TRUE) %>% na.omit %>% select(-layer) %>% 
  filter(x > -910000)
proj = inla.mesh.projector(mesh, loc = as.matrix(locs))

locs <- locs %>% 
  mutate(ymean = inla.mesh.project(proj, field = rowMeans(pred_y)),
         ylci = inla.mesh.project(proj, field = apply(pred_y, 1, quantile, 0.025)),
         yuci = inla.mesh.project(proj, field = apply(pred_y, 1, quantile, 0.975)),
         zmean = inla.mesh.project(proj, field = rowMeans(pred_z)),
         zlci = inla.mesh.project(proj, field = apply(pred_z, 1, quantile, 0.025)),
         zuci = inla.mesh.project(proj, field = apply(pred_z, 1, quantile, 0.975)))
```

### Density predictions


```r
# rasterize the boundary area to get locations for plot
ymean_plot <- ggplot() + 
  geom_raster(data = locs, aes(x = x, y = y, fill = ymean)) + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") + 
  scale_fill_viridis_c(name = "Density\n(mean)", option = "plasma") + 
  labs(x = "", y = "") 

yci_plot <- ggplot() + 
  geom_raster(data = locs, aes(x = x, y = y, fill = yuci - ylci)) + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") + 
  scale_fill_viridis_c(name = "Density\n(95% CI)", option = "plasma") + 
  labs(x = "", y = "") 

ymean_plot + yci_plot + plot_annotation(tag_levels = "a", tag_suffix = ")")
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/density_plot-1.png)<!-- -->

```r
ggsave(plot = ymean_plot, filename = "~/Google Drive/SIDEPROJ/Papers/sea cucumber/density_preds.jpg", dpi = 300, width = 80, height = 60, units = "mm")
```

### Occurrence predictions


```r
zmean_plot <- ggplot() + 
  geom_raster(data = locs, aes(x = x, y = y, fill = zmean)) + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") + 
  scale_fill_viridis_c(name = "Occurence\n(mean)", option = "plasma") + 
  labs(x = "", y = "") 

zci_plot <- ggplot() + 
  geom_raster(data = locs, aes(x = x, y = y, fill = zuci - zlci)) + 
  geom_sf(data = map_sf, colour = "lightgrey", fill = "lightgrey") + 
  scale_fill_viridis_c(name = "Occurrence\n(95% CI)", option = "plasma") + 
  labs(x = "", y = "") 

zmean_plot + zci_plot + plot_annotation(tag_levels = "a", tag_suffix = ")")
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/occ_plot-1.png)<!-- -->
## Comparison against quotas

We use the shapes from Luis to make predictions of units within the estates for comparison against granted quotas. 


```r
density_mean <- locs %>% 
  select(x, y, ymean) %>% 
  mutate(ymean = ymean * 100 * 100) %>% 
  rasterFromXYZ(crs = crs_ea_mex)

density_lci <- locs %>% 
  select(x, y, ylci) %>% 
  mutate(ylci = ylci * 100 * 100) %>% 
  rasterFromXYZ(crs = crs_ea_mex)

density_uci <- locs %>% 
  select(x, y, yuci) %>% 
  mutate(yuci = yuci * 100 * 100) %>% 
  rasterFromXYZ(crs = crs_ea_mex)

density <- stack(list(mean = density_mean, 
                 lci = density_lci, 
                 uci = density_uci))

# for the null models
ymean <- exp(mods$null$summary.fixed[2, 1]) * 100 * 100
ylci <- exp(mods$null$summary.fixed[2, 3]) * 100 * 100
yuci <- exp(mods$null$summary.fixed[2, 5]) * 100 * 100
nulldensity_mean <- density_mean
nulldensity_lci <- density_lci
nulldensity_uci <- density_uci
values(nulldensity_mean) <- ifelse(is.na(values(nulldensity_mean)), NA, ymean)
values(nulldensity_lci) <- ifelse(is.na(values(nulldensity_lci)), NA, ylci)
values(nulldensity_uci) <- ifelse(is.na(values(nulldensity_uci)), NA, yuci)

nulldensity <- stack(list(nullmean = nulldensity_mean,
                          nulllci = nulldensity_lci,
                          nulluci = nulldensity_uci))
```


```r
# multiply by depth raster so we only predict in < 30m depth
depth_ras <- fasterize::fasterize(depth_sf, density[[1]])
density <- density * depth_ras
names(density) <- c("mean", "lci", "uci")
nulldensity <- nulldensity * depth_ras
names(nulldensity) <- c("nullmean", "nulllci", "nulluci")

# get predictions
quota_preds <- raster::extract(density, quotas_sp, fun = sum, na.rm = TRUE) %>% as_tibble
quota_nullpreds <- raster::extract(nulldensity, quotas_sp, fun = sum, na.rm = TRUE) %>% as_tibble
quota_area <- raster::extract(density, quotas_sp) %>% 
  map_dfr(function(x) {
    tibble(suitable_area = x %>% na.omit %>% nrow)
  })

bind_cols(quotas, quota_preds, quota_nullpreds, quota_area) %>% 
  st_set_geometry(NULL) %>% 
  mutate(uniform = 0.3*10^4*total_area*0.1,
         smean = mean*0.1,
         slci = lci*0.1,
         nmean = nullmean*0.1,
         nlci = nulllci*0.1) %>% 
  arrange(Estate) %>% 
  select(Estate, total_area, suitable_area, # info on estates
         mean, lci, uci, # spatial model estimates
         nullmean, nulllci, nulluci, # non-spatial model estimates
         uniform, smean, slci, nmean, nlci) %>% # quotas
  write_csv("results/quotas.csv") %>% 
  kable
```



Estate    total_area   suitable_area       mean         lci         uci    nullmean    nulllci     nulluci      uniform      smean        slci       nmean       nlci
-------  -----------  --------------  ---------  ----------  ----------  ----------  ---------  ----------  -----------  ---------  ----------  ----------  ---------
1           2730.271             791   257437.9    82113.88    613908.3    258715.1   228531.3    291927.9     819081.4   25743.79    8211.388    25871.51   22853.13
2           8276.371             900   294784.3    98401.80    686625.3    294366.1   260023.0    332155.6    2482911.4   29478.43    9840.180    29436.61   26002.30
3          17363.769            2344   920247.7   385411.43   1881234.3    766660.3   677215.3    865080.8    5209130.6   92024.77   38541.143    76666.03   67721.53
4          14141.100            1119   418349.2   145578.82    962401.3    365995.2   323295.2    412980.1    4242329.9   41834.92   14557.882    36599.52   32329.52
5           6577.154            1201   424387.5   162456.22    911240.9    392815.3   346986.2    443243.2    1973146.1   42438.75   16245.622    39281.53   34698.62
6          17134.589            2368   776679.7   272688.34   1762225.6    774510.0   684149.3    873938.3    5140376.8   77667.97   27268.834    77451.00   68414.93
7          10100.025            1575   448223.3   143671.63   1059047.7    515140.8   455040.2    581272.3    3030007.6   44822.33   14367.163    51514.08   45504.02
8          34743.706            3125   890309.6   270800.55   2164728.2   1022104.7   902857.5   1153318.1   10423111.7   89030.96   27080.055   102210.47   90285.75
9           2860.709             520   161407.7    53686.27    377971.5    170078.2   150235.5    191912.1     858212.7   16140.77    5368.627    17007.82   15023.55

Plot showing the estimated quotas (mean and LCI) against the granted quotas. 


```r
dat <- bind_cols(read_csv("results/quotas.csv"), 
                 tibble(granted_quota = c(6750, 6750, 6250, 6750, 6250, 47088, 27390, 24315, 12158))) %>% 
  mutate(Estate = factor(Estate)) %>% 
  mutate(legend = " Granted Quota")

                 
table1 <- dat %>% 
  select(Estate, total_area, suitable_area, mean, lci, uci, nullmean, nulllci, nulluci) %>% 
  write_csv("results/table1.csv")

fig3_dat <- dat %>% 
  select(Estate, spatial_mean = mean, spatial_lci = lci, nonspatial_mean = nullmean, nonspatial_lci = nulllci) %>% 
  mutate_at(vars(spatial_mean, spatial_lci, nonspatial_mean, nonspatial_lci), .funs = function(x) x*0.1) %>% 
  gather(key, value, -Estate) %>%
  separate(key, into = c("model", "measure")) %>% 
  spread(measure, value) %>% 
  mutate(Model = factor(model, levels = c("spatial", "nonspatial"), labels = c(" Spatial model (-95% CI)", " Non-spatial model (-95% CI)")))

ggplot() + 
  geom_bar(data = dat, aes(x = Estate, y = granted_quota, fill = legend), stat = "identity", colour = "black") + 
  scale_fill_manual(values = "white") + 
  geom_point(data = fig3_dat, aes(x = Estate, y = mean, shape = Model), position = position_dodge(1), size = 3) + 
  geom_point(data = fig3_dat, aes(x = Estate, y = lci, group = Model), shape = 95, size = 2, position = position_dodge(1)) + 
  geom_linerange(data = fig3_dat, aes(x = Estate, ymin = lci, ymax = mean, group = Model), position = position_dodge(1)) + 
  ylab("Quota (in number of pieces of sea cucumber)") + 
  theme(legend.position = c(0.15, 0.8),
        legend.title = element_blank())
```

![](C:/Users/lg1u16/GIT_PROJECTS/SIDEPROJ/sea_cucumber/doc/sea_cucumber_range500_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
ggsave(filename = "~/Google Drive/SIDEPROJ/Papers/sea cucumber/3_quotas.jpg", dpi = 300, width = 140, height = 80, units = "mm")
```

