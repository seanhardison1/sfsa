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



# inset map of US and Canada---- 
us_ca <- 
  rnaturalearthhires::countries10 %>% 
  filter(NAME_EN %in% c("United States of America",
                        "Canada")) %>% 
  group_by(NAME_EN) %>% 
  dplyr::summarise() 

# bathy details (ETOPO1)-----
# Amante, C. and B.W. Eakins, 2009. ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. 
# NOAA Technical Memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA. doi:10.7289/V5C8276M

bathy <- raster::raster(here::here("data/exportImage.tiff")) 
bathy_df <- 
  bathy %>% 
  dream::rst_to_tib(., var_name = "depth") %>% 
  filter(depth < 200) %>% 
  mutate(depth = ifelse(depth >= 0, 0, depth),
         depth = ifelse(depth < -500, -500, depth))

# contours for map
cont_out <- NULL
conts <-  seq(-300, -50, by = 50)
for (cont in conts){
  cont_int <- bathy %>%
    stars::st_as_stars() %>% 
    st_contour(x = .,
               breaks = cont,
               na.rm = T,
               contour_lines = T) %>%
    # smoothr::drop_crumbs(., threshold = 25) %>%
    # st_intersection(., sf_depth) %>%
    st_transform(., 4326) %>% 
    dplyr::rename(contour = exportImage) %>% 
    mutate(contour = factor(contour, levels = rev(conts)))
  assign("cont_out", rbind(cont_int, cont_out))
}


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
  st_as_sf(., coords = c("longitude","latitude"), crs = 4326) %>% 
  st_transform(., ncrs)

# write out locations----
# write_sf(locs, here::here("data/island_locs.kml"))

# Coastline and US-Canada maritime border------

# US-Canada maritime border and coastline citation
# Flanders Marine Institute (2023). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), 
# version 12. Available online at https://www.marineregions.org/. https://doi.org/10.14284/632


# some tricks to make sure we can fill in the land polygon
bounds  <- st_as_sfc(st_bbox(locs)) %>% st_buffer(., dist = 500) %>% st_transform(., crs = 4326)

coastline <- st_read(here::here("data/coastlines")) %>% 
  filter(SOVEREIGN1 %in% c("United States", "Canada")) %>% 
  st_make_valid() %>% 
  st_intersection(., bounds) %>% 
  st_make_valid()

offshore <- st_union(st_geometry(coastline))
land <- st_difference(st_geometry(bounds), 
                           offshore) %>% 
  st_collection_extract(., "POLYGON") %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_as_sf()

# coastline
us_can <- st_read(here::here("data/World_EEZ_v12_20231025/bounds"))

maritime_border <- us_can %>% 
  st_intersection(., bounds) %>% 
  filter(SOVEREIGN2 %in% c("United States", "Canada")) %>% 
  st_transform(., ncrs) %>% 
  summarise()

# US-Canada land border----
ca_brd <- st_read(here::here("data/kx-canada-and-us-border-SHP")) %>% 
  st_transform(., ncrs) %>%
  st_intersection(., land %>% st_transform(., crs = ncrs)) %>% 
  summarise() 

# GOM map making-----
n <- 100
depth_pal <- rev(pals::kovesi.linear_blue_95_50_c20(n))
depth_pal <- pals::ocean.ice(n)
x_min <- -70.5
x_max <- -66.0
y_min <- 43.0
y_max <- 45.25

# positions of island labels
loc_lab_df <- tibble(island_long =  c("Stratton Island",
                                      "Outer Green<br>Island",
                                      "Jenny Island",
                                      "Pond Island",
                                      "Eastern Egg Rock",
                                      "Metinic Island",
                                      "Matinicus Rock",
                                      "Seal Island",
                                      "Ship Island",
                                      "Petit Manan<br>Island",
                                      "Eastern Brothers<br>Island",
                                      "Machias Seal<br>Island"),
                     longitude = c(-70.311,
                                   -70.4,
                                   -70.5,
                                   -70,
                                   -69.75,
                                   -69.75,
                                   -69.1,
                                   -68.75,
                                   -68.6,
                                   -67.9,
                                   -68.2,
                                   -67.3),
                     latitude = c(43.1,
                                  43.6,
                                  44.1,
                                  43.5,
                                  43.6,
                                  44.1,
                                  43.55,
                                  43.7,
                                  43.9,
                                  44.1,
                                  44.75,
                                  44.25))

# This script gets us 75% of the way there. The rest we do in Inkscape.

map <-
  ggplot() +
  geom_raster(data = bathy_df,
              aes(x = longitude,
                  y = latitude, fill = depth),
              show.legend = T) +
  geom_sf(data = land,# ,
          # fill = "#d3f8e2ff",
          fill = "#b4d4c1",
          linewidth = 0.1) +
  geom_sf(data = ca_brd,
          linewidth = 0.25) +
  geom_sf(data = maritime_border,
          linewidth = 0.25,
          linetype = 2) +
  geom_point(data = locs %>% 
               st_transform(., crs = 4326) %>% 
            dream::sfc_as_cols(., names = c("longitude",
                                            "latitude")),
          aes(x = longitude, y = latitude)) +
  ggrepel::geom_text_repel(data = locs,
                aes(label = island_long,
                    geometry = geometry),
                   stat = "sf_coordinates",
                   size = 3,
                   max.overlaps = 10,
                   segment.color = "black",
                   fill = NA,
                   label.color = NA,
                nudge_x = -0.2)+
  labs(fill = "Depth (m)") +
  coord_sf(xlim = c(x_min, x_max),
           ylim = c(y_min, y_max)) +
  scale_fill_gradientn(colors=as.vector(depth_pal),
                       # transform = "pseudo_log",
                       breaks =c(-400, -300,  -200, -100,  0),
                       labels = c("<400", "300", "200", "100", "0")) +
  theme(panel.background = element_rect(fill = "#6892beff"),
        text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    linewidth = 0.5, 
                                    fill = "transparent"),
        axis.title = element_blank(),
        legend.position = c(.25, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.direction = "horizontal") +
  annotation_scale(location = "br") +
  annotation_north_arrow(which_north = "true",
                         style = north_arrow_minimal(),
                         location = "br",
                         pad_y = unit(10, "mm"))
map

ggsave(map, filename = here::here("figs/map2.svg"),
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
          element_rect(fill = "#6892beff"))
inset
ggsave(inset, filename = here::here("figs/inset.svg"),
       width = 4,
       height = 7)
