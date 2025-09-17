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
library(officer)
library(lubridate)
library(gt)


# ChatGPT function to help with data table creation-------
# This function takes a numeric vector of years, sorts it, 
# identifies consecutive sequences, and returns a string like "2012-2014, 2016"
compress_years <- function(x) {
  x <- sort(unique(x))  # sort and ensure unique
  if (length(x) == 0) return(NA_character_)
  
  # Identify groups of consecutive numbers
  # `diff(x) != 1` will be TRUE when a "break" in consecutive years occurs
  group_breaks <- c(TRUE, diff(x) != 1)
  # Convert TRUE/FALSE breaks into group IDs via cumulative sum
  groups <- cumsum(group_breaks)
  
  # Split the years into sub-vectors by these group IDs
  split_years <- split(x, groups)
  
  # For each group of consecutive years, convert to "start-end" if length > 1
  format_range <- function(yrs) {
    if (length(yrs) == 1) {
      as.character(yrs)
    } else {
      paste0(min(yrs), "-", max(yrs))
    }
  }
  
  # Apply the formatting to each group and then combine with commas
  compressed <- vapply(split_years, format_range, character(1))
  paste(compressed, collapse = ", ")
}



# herring per hour-----
load(here::here("data/processed_herr_bird5_simple.rdata"))

# island locations-----
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


island_subset <- st_read(here::here("data/island_subset.kml")) %>% 
  st_zm() %>% 
  st_transform(.,ncrs)

# bind in island lon/lat for spatial models
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


# N years by study site
wrd_tbl <- df %>% 
  dplyr::select(Island = island_long, year) %>% 
  distinct() %>% 
  group_by(Island) %>% 
  dplyr::summarise(`N unique years` = n(),
                   `Data years` = compress_years(year)) %>% 
  mutate(Island = as.character(Island)) %>% 
  arrange(Island) %>% 
  gt() %>% 
  as_word()


doc <- read_docx() %>%
  body_add_xml(str = wrd_tbl, pos = "after")

print(doc, target = here::here("data/seabird_data_collection_table.docx"))
