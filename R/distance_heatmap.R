library(sdmTMB)
library(reshape2)
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
  arrange(., latitude) %>% 
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

dist_mat <- apply(sf::st_distance(locs), 1, as.numeric)
row.names(dist_mat) <- locs$island_long
colnames(dist_mat) <- locs$island_long


# Mask upper triangle
dist_mat_lower <- dist_mat
dist_mat_lower[lower.tri(dist_mat_lower, diag = FALSE)] <- NA

dist_map_viz <- melt(dist_mat_lower, na.rm = TRUE)

# dist_map_viz <- melt(dist_mat)
dist_heatmap <- 
  ggplot() + 
  geom_tile(data = dist_map_viz %>% 
              mutate(value = ifelse(value == 0, NA, value)),
            aes(y = Var1, x = Var2, fill = value)) +
  geom_text(data = dist_map_viz %>% 
              filter((value < 33.96 & value != 0)) %>% 
              mutate(text = "X"),
            aes(y = Var1, x = Var2, label = text)) +
  geom_text(data = dist_map_viz %>% 
              filter((value < 33.96 & value != 0 & value < 27.26)) %>% 
              mutate(text = "X"),
            aes(y = Var1, x = Var2, label = text),
            size = 5,
            fontface = "bold") +
  labs(fill = "Distance (km)") +
  theme_void() +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_discrete(expand = c(0.05, 0.05)) +
  ggsci::scale_fill_material("light-blue") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, color = "black"),
    axis.text.y.right = element_text(size = 10, color = "black"),  # move labels to right
    axis.ticks.y.right = element_line(),
    axis.text.y.left = element_blank(),  # remove left labels
    axis.ticks.y.left = element_blank(),
    axis.line.y.right = element_blank()
  ) +
  scale_y_discrete(position = "right")  # actually position labels on right
  
dist_heatmap

ggsave(dist_heatmap, 
       file = here::here("figs/dist_heatmap.png"),
       width = 9.5,
       height = 7.5,
       dpi = 300)  
