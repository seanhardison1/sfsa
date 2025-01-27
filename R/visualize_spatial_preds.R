library(sdmTMB)
library(sdmTMBextra)
library(tidyverse)
library(sf)
library(ggrepel)
library(rnaturalearth)
library(patchwork)
library(glmmTMB)
library(DHARMa)
library(tsibble)
library(sp)
library(lubridate)
library(readxl)
library(mgcv)
library(terra)
library(smoothr)
library(ggtext)

# data processing-----
process <- F

## herring per hour-----
load(here::here("data/processed_herr_bird5_simple.rdata"))

## island locations-----
ncrs <- "+proj=utm +zone=19 +datum=NAD83 +units=km +no_defs"

locs <- 
  tibble(island_long = c("Matinicus Rock", "Seal Island NWR", "Machias Seal Island",
                         "Eastern Egg Rock", "Petit Manan","Stratton Island",
                         "Outer Green Island","Jenny Island", "Pond Island NWR",
                         "Metinic Island", "Ship Island",
                         "Eastern Brothers Island"),
         island_short = c("MRNWR", "SINWR", "MSI", "EER", "PMNWR", "STI",
                          "OGI", "JI", "PINWR", "METIN", "SHIP","EBI"),
         latitude = c(43.784213,43.884715,44.501395, 
                      43.860507, 44.367173, 43.504799, 
                      43.650106, 43.764960, 43.739154, 
                      43.883776,44.233802,44.557551), 
         longitude = c(-68.854548,-68.744878,-67.101474,
                       -69.382377,-67.865695,-70.311804,
                       -70.123831,-69.908692,-69.770649, 
                       -69.126051, -68.441193, -67.435477)) %>% 
  st_as_sf(., coords = c("longitude","latitude"), crs = 4326) %>% 
  mutate(island_gam = case_when(island_long == "Matinicus Rock" ~ "MR",
                                island_long == "Seal Island NWR" ~ "SINWR",
                                island_long == "Machias Seal Island" ~ "MSI",
                                island_long == "Eastern Egg Rock" ~ "EER",
                                island_long == "Petit Manan" ~ "PMI",
                                island_long == "Stratton Island" ~ "STI",
                                island_long == "Outer Green Island" ~ "OGI",
                                island_long == "Jenny Island" ~ "JI",
                                island_long == "Pond Island NWR" ~ "PINWR",
                                island_long == "Metinic Island" ~ "METIN",
                                island_long == "Ship Island" ~ "SHIP",
                                island_long == "Eastern Brothers Island" ~ "EB")) %>% 
  st_transform(., ncrs)

## bind in island lon/lat----
df <- tern_feeding_index_raw_simple %>% 
  filter(time != 0,
         month != 8) %>% # drop august observations (there are very few overall)
  left_join(.,locs, by = c("island" = "island_gam")) %>% 
  st_as_sf() %>% 
  dream::sfc_as_cols(.,names = c("lon", "lat")) %>% #devtools::install_github("seanhardison1/dream")
  st_set_geometry(NULL) %>% 
  as.data.frame() %>% 
  mutate(ym_fac = factor(ym),
         m_fac = factor(month),
         island = factor(island),
         loc_nest_year = factor(paste(loc_nest, year)),
         ym_num = as.numeric(ym),
         island_long = factor(island_long,
                              levels = unique(.$island_long[order(.$lon)])))

## load fitted model----
load(here::here("data/fitted_dl_model_1014.rdata"))

# bathy details (ETOPO1)-----
bathy <- raster::raster(here::here("data/exportImage.tiff")) 
bathy_df <- 
  bathy %>% 
  dream::rst_to_tib(., var_name = "depth") %>% 
  filter(depth < 200)

# make coastline-----
topo_rast <- raster::raster(here::here("data/exportImage.tiff")) %>% 
  disaggregate(., fact = 5)
topo_rast[topo_rast < 0] <- NA  
topo_rast[topo_rast > 0] <- 1
topo_rast2 <- terra::rast(topo_rast)
coast_sf <- as.polygons(topo_rast2, values=TRUE, dissolve=TRUE) %>% 
  st_as_sf() %>% 
  smooth(method = "densify") %>%
  smooth(method = "chaikin")


# data for prediction-----
ndf <- df %>%
  dplyr::select(lat, lon, island, island_long) %>%
  distinct() %>%
  mutate(loc_nest_year = factor(NA),
         ym_fac = factor(NA),
         year_fac = factor(NA),
         year = NA)

# visualize spatial REs----
s_rfs <- 
  predict(m1_simple, re_form = NULL) %>% 
  dplyr::select(island_long,
                lon, lat, 
                omega_s1, 
                omega_s2) %>% 
  distinct() %>% 
  st_as_sf(., coords = c("lon","lat"), crs = ncrs) %>% 
  st_transform(., crs = 4326)

# visualize
x_min <- -70.5
x_max <- -66.0
y_min <- 43.0
y_max <- 45.25
depth_pal <- pals::kovesi.diverging_linear_bjr_30_55_c53(100)
prob_res <- 
  s_rfs %>% 
  ggplot() +
  geom_raster(data = bathy_df,
              aes(x = longitude,
                  y = latitude, fill = depth),
              show.legend = F) +
  geom_sf(data = coast_sf,
          show.legend = F,
          fill = "#d3f8e2ff",
          color = "black",
          linewidth = 0.01) +
  geom_sf(aes(color = omega_s1),
             size = 4) + 
  geom_text_repel(aes(label = island_long,
                        geometry = geometry),
                      stat = "sf_coordinates",
                  size = 4,
                  nudge_x = 0.2,
                  nudge_y = -0.3) + 
  scale_color_gradientn(colors=as.vector(depth_pal)) +
  coord_sf(xlim = c(x_min, x_max),
           ylim = c(y_min, y_max)) +
  labs(color = "Spatial RE",
       title = "Probability model") +
  theme(panel.background = 
          element_rect(fill = "#90daeeff"),
        text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(
          colour = "black", linewidth = 0.5, fill = "transparent"),
        axis.title = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.background = element_rect(fill = "transparent"))


catch_res <- 
  s_rfs %>% 
  ggplot() +
  geom_raster(data = bathy_df,
              aes(x = longitude,
                  y = latitude, fill = depth),
              show.legend = F) +
  geom_sf(data = coast_sf,
          show.legend = F,
          fill = "darkgreen",
          color = "black",
          linewidth = 0.01) +
  geom_sf(aes(color = omega_s2),
          size = 4) + 
  geom_text_repel(aes(label = island_long,
                      geometry = geometry),
                  stat = "sf_coordinates",
                  size = 4,
                  nudge_x = 0.2,
                  nudge_y = -0.3) + 
  scale_color_gradientn(colors=as.vector(depth_pal)) +
  coord_sf(xlim = c(x_min, x_max),
           ylim = c(y_min, y_max)) +
  labs(color = "Spatial RE",
       title = "Positive catch model") +
  theme(panel.background = 
          element_rect(fill = "#90daeeff"),
        text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(
          colour = "black", linewidth = 0.5, fill = "transparent"),
        axis.title = element_blank(),
        legend.position = c(0.1, 0.8),
        legend.background = element_rect(fill = "transparent"))

spatial_res <- 
  prob_res + catch_res + 
  plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")",
                  tag_prefix = "(")


ggsave(spatial_res,
       filename = here::here("figs/spatial_REs.png"),
       dpi = 300,
       height = 8.91,
       width = 7.39)

# spatiotemporal predictions-----
## data for prediction-----
ndf <- df %>%
  dplyr::select(year_fac, year) %>%
  distinct() %>%
  expand_grid(.,
              df %>%
                dplyr::select(lat, lon, island, island_long) %>%
                distinct()) %>%
  mutate(loc_nest_year = factor(NA),
         ym_fac = factor(NA))

## island-level predictions-----
load(here::here("data/st_model_pred_sims.rdata"))

## process predictions-----
isl_preds <-
  p %>%
  bind_cols(., ndf %>%
              dplyr::select(-loc_nest_year,-ym_fac)) %>%
  gather(sim, value, -year_fac:-island_long) %>%
  group_by_at(vars(year_fac:island_long)) %>%
  dplyr::summarise(est_link = mean(value),
                   sd_link = sd(value)) %>%
  mutate(fit_resp = exp(est_link),
         sd_resp = exp(sd_link),
         upr_resp = exp(est_link + (2 * sd_link)),
         lwr_resp = exp(est_link - (2 * sd_link))) %>% 
  st_as_sf(., coords = c("lon","lat"), crs = ncrs) %>% 
  st_transform(., crs = 4326)

## coastline----
coast_simple <- 
  rnaturalearthhires::states10 %>% 
  filter(name_en %in% c("Maine",
                        "New Brunswick",
                        "Quebec"))

## visualize 1988-2004----
facets <- c(1988:2004)
h_pal <- pals::kovesi.linear_bgy_10_95_c74(100)
plot_list <- lapply(facets, function(facet) {
  proj_subset <- isl_preds %>% filter(year == facet)
  ggplot() +
    geom_sf(data = coast_simple,
            show.legend = F,
            fill = "#d3f8e2ff",
            color = "black",
            linewidth = 0.02) +
    geom_sf(data = proj_subset,
            aes(color = fit_resp),
            size = 3) + 
    coord_sf(xlim = c(x_min, x_max),
             ylim = c(y_min, y_max)) +
    labs(color = "Herring provisioning\nindex") +
    scale_color_gradientn(colors=as.vector(h_pal)) +
    facet_wrap(~year) +
    dream::theme_fade() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(angle = 45,
                                     size = 6),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.key.height = unit(0.1, "cm"),
          plot.title = element_markdown(size = 10,
                                        hjust = 0.5),
          panel.grid.major = element_line(linewidth = 0.1,
                                          color = "grey70"))
})
combined_plot1 <- wrap_plots(plotlist = plot_list, ncol = 5)
ggsave(combined_plot1,
       filename = here::here("figs/1988_2004_predictions.png"),
       dpi = 300,
       height = 10,
       width = 8)


## visualize 2005-2021----
facets <- c(2005:2021)
plot_list <- lapply(facets, function(facet) {
  proj_subset <- isl_preds %>% filter(year == facet)
  ggplot() +
    geom_sf(data = coast_simple,
            show.legend = F,
            fill = "#d3f8e2ff",
            color = "black",
            linewidth = 0.02) +
    geom_sf(data = proj_subset,
            aes(color = fit_resp),
            size = 3) + 
    coord_sf(xlim = c(x_min, x_max),
             ylim = c(y_min, y_max)) +
    labs(color = "Herring provisioning\nindex") +
    scale_color_gradientn(colors=as.vector(h_pal)) +
    facet_wrap(~year) +
    dream::theme_fade() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(angle = 45,
                                     size = 6),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.key.height = unit(0.1, "cm"),
          plot.title = element_markdown(size = 10,
                                        hjust = 0.5),
          panel.grid.major = element_line(linewidth = 0.1,
                                          color = "grey70"))
})
combined_plot2 <- wrap_plots(plotlist = plot_list, ncol = 5)
ggsave(combined_plot2,
       filename = here::here("figs/2005_2021_predictions.png"),
       dpi = 300,
       height = 10,
       width = 8)
