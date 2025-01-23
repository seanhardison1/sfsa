library(tidyverse)
library(raster)
library(sf)
library(terra)
library(pals)
library(gganimate)
library(stars)
library(ggspatial)
library(ggnewscale)
library(patchwork)
library(ggtext)
library(ggrepel)
library(smoothr)

# US/Canada border----
ca_brd <- st_read(here::here("data/kx-canada-and-us-border-SHP"))

# inset map of US and Canada---- 
us_ca <- 
  rnaturalearthhires::countries10 %>% 
  filter(NAME_EN %in% c("United States of America",
                        "Canada")) %>% 
  group_by(NAME_EN) %>% 
  dplyr::summarise() 

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

# island locations-----
ncrs <- "+proj=utm +zone=19 +datum=NAD83 +units=km +no_defs"

locs <- 
  tibble(island_long = c("Matinicus Rock", "Seal Island", "Machias Seal Island",
                         "Eastern Egg Rock", "Petit Manan Island","Stratton Island",
                         "Outer Green Island","Jenny Island", "Pond Island",
                         "Metinic Island", "Ship Island",
                         "Eastern Brothers Island"),
         island_short = c("MRNWR", "SINWR", "MSI", "EER", "PMNWR", "STI",
                          "OGI", "JI", "PINWR", "METIN", "SHIP","EBI"),
         latitude = c(43.784213,
                      43.884715,
                      44.501395, 
                      43.860507, 44.367173, 43.504799, 
                      43.650106, 43.764960, 43.739154, 
                      43.883776,44.233802,44.557551), 
         longitude = c(-68.854548,-68.744878,-67.101474,
                       -69.382377,-67.865695,-70.311804,
                       -70.123831,-69.908692,-69.770649, 
                       -69.126051, -68.441193, -67.435477)) %>% 
  st_as_sf(., coords = c("longitude","latitude"), crs = 4326)

# write out locations----
write_sf(locs, here::here("data/island_locs.kml"))

# GOM map making-----
depth_pal <- pals::kovesi.linear_blue_95_50_c20(100)
x_min <- -70.5
x_max <- -66.0
y_min <- 43.0
y_max <- 45.25

map <- 
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
  geom_sf(data = ca_brd,
          linewidth = 0.5) +
  geom_sf(data = locs) +
  ggrepel::geom_text_repel(data = locs,
                           aes(label = island_long,
                               geometry = geometry),
                           stat = "sf_coordinates",
                           size = 3,
                           max.overlaps = 10,
                           segment.color = "black",
                           nudge_y = -0.5)+
  coord_sf(xlim = c(x_min, x_max),
           ylim = c(y_min, y_max)) +
  scale_fill_gradientn(colors=as.vector(depth_pal)) +
  theme(panel.background = 
          element_rect(fill = "#90daeeff"),
        text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(
          colour = "black", linewidth = 0.5, fill = "transparent"),
        axis.title = element_blank()) +
  annotation_scale(location = "br") +
  annotation_north_arrow(which_north = "true",
                         style = north_arrow_minimal(),
                         location = "tr",
                         pad_y = unit(10, "mm"))

map
ggsave(map, filename = here::here("figs/map.svg"),
       width = 7,
       height = 5.5)

# island inset map-----
bbox <- st_as_sfc(st_bbox(c(xmin = x_min, xmax = x_max, 
                            ymin = y_min, ymax = y_max),
                          crs = 4326))
inset <- 
  ggplot() +
  geom_sf(data = us_ca,
          fill = "#d3f8e2ff") + 
  geom_sf(data = bbox,
          fill = 'transparent',
          color= "red",
          linewidth = 1) +
  coord_sf(xlim = c(-85.5, -60.5),
           ylim = c(25, 55)) + 
  theme_void() + 
  theme(panel.border = element_rect(
          colour = "black", linewidth = 0.5, fill = "transparent"),
        panel.background = 
          element_rect(fill = "#6994bfff"))
inset
ggsave(inset, filename = here::here("figs/inset.svg"),
       width = 4,
       height = 7)
