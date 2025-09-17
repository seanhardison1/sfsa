library(tidyverse)
library(sf)
library(readxl)
library(patchwork)
library(lubridate)
library(hms)
library(tsibble)

# Machias Seal Island----
msi_ct <-
  read_excel(here::here("data/seabird_data/MSI_feeding/MSI Seabird Diet.xlsx"),
             sheet = 'Common Tern', skip = 5)

msi_all_disagg <-
  msi_ct %>% 
  dplyr::select(stint = StintID,
                species = Species, 
                location = Observer_Location,
                nFish = Number_of_Items, 
                start = Date_Time_Start,
                end = Date_Time_End, 
                island = Island,
                prey = Prey_Item,
                time_arrive = Time_Arrive,
                time_depart = Time_Depart,
                nest = Nest,
                provider = Provider) %>% 
  mutate(month = month(end),
         year = year(end),
         ym = yearmonth(paste(year, month, sep = "-")),
         nFish = ifelse(is.na(nFish), 1, nFish),
         stint = str_remove_all(stint, "-"),
         loc_nest_stint = paste(location, nest, stint),
         time = as.numeric(difftime(end, start,
                                    units = "hours")),
         loc_nest = paste(location, nest)) %>%
  filter(nest != "Unknown",
         time < 10)

## Visuals-----
a <- ggplot(msi_all_disagg) + 
  geom_histogram(aes(time)) +
  labs(title = "stint length",
       x = "Time (hrs)") +
  facet_wrap(~island) +
  dream::theme_fade() 

b <-
  ggplot(msi_all_disagg) + 
  geom_histogram(aes(month)) +
  labs(title = "Observation month",
       x = "Month") +
  facet_wrap(~island) +
  dream::theme_fade()

## Summarise----
msi_summ <-
  msi_all_disagg %>%
    group_by(year, month, ym, island, loc_nest, 
             stint, time, loc_nest_stint) %>% 
    dplyr::summarise(n_herring = sum(nFish[prey %in% c("R","r")]),
                     time = unique(time),
                     herring_per_hour = sum(nFish[prey %in% c("R","r")])/time)%>% 
  mutate(ds = 1)

# length(unique(msi_all_disagg$loc_nest_stint))
# length(unique(msi_summ$loc_nest_stint))

c <-
  ggplot(msi_summ) + 
  geom_point(aes(y = herring_per_hour, x = ym),
             size = 0.1) +
  labs(title = "Herring per hour (n Herring/stint length)",
       x = "Year-month") +
  facet_wrap(~island) +
  dream::theme_fade() +
  theme(axis.text.x = element_text(size = 7))

ggsave(c, filename = here::here("figs/hph_msi.png"),
       dpi = 300,
       width = 7,
       height = 5)

# Common tern (EER, JI, MR, OGI, PINWR, SINWR, STI)----
cote_diet_aud_disagg <-
  read_excel(here::here("data/seabird_data/audubon_feeding/Audubon_All_COTE_feedings_from_SP.xlsx")) %>% 
  dplyr::select(stint = StintID,
                species = Species, 
                nFish = Number_of_Items, 
                location = Observer_Location,
                start = Date_Time_Start,
                end = Date_Time_End, 
                island = Island, 
                prey = Prey_Item,
                time_arrive = Time_Arrive,
                time_depart = Time_Depart,
                nest = Nest,
                provider = Provider) %>% 
  mutate(month = month(end),
         year = year(end),
         ym = yearmonth(paste(year, month, sep = "-")),
         nFish = ifelse(is.na(nFish), 1, nFish),
         stint = str_remove_all(stint, "-"),
         loc_nest_stint = paste(location, nest, stint),
         time = as.numeric(difftime(end, start,
                                    units = "hours")),
         loc_nest = paste(location, nest)) %>%
  filter(time < 10,
         month != 1)

## Visuals----
# start/end date by island
cote_diet_aud_disagg %>% 
  group_by(island) %>% 
  dplyr::summarise(year_min = min(year),
                   year_max = max(year))

a <- ggplot(cote_diet_aud_disagg) + 
  geom_histogram(aes(time)) +
  labs(title = "stint length",
       x = "Time (hrs)") +
  facet_wrap(~island) +
  dream::theme_fade() 
  
b <- ggplot(cote_diet_aud_disagg) + 
  geom_histogram(aes(month)) +
  labs(title = "Observation month",
       x = "Month") +
  facet_wrap(~island) +
  dream::theme_fade()


## Summarise----
cote_diet_aud_summ <-
  cote_diet_aud_disagg %>%
  group_by(year, month, ym, island, loc_nest, 
           stint, loc_nest_stint, time) %>% 
  dplyr::summarise(n_herring = sum(nFish[prey %in% c("R","r")]),
                   time = unique(time),
                   herring_per_hour = sum(nFish[prey %in% c("R","r")])/time)%>% 
  mutate(ds = 2)

# length(unique(cote_diet_aud_summ$loc_nest_stint))
# length(unique(cote_diet_aud_disagg$loc_nest_stint))


c <- ggplot(cote_diet_aud_summ) + 
  geom_point(aes(y = herring_per_hour, x = ym),
             size = 0.1) +
  labs(title = "Herring per hour (n Herring/stint length)",
       x = "Year-month") +
  facet_wrap(~island) +
  dream::theme_fade() +
  theme(axis.text.x = element_text(size = 7))

ggsave(c, filename = here::here("figs/hph_aud_modern.png"),
       dpi = 300,
       width = 7,
       height = 5)

# Historical data (ARTE, COTE, ROST)----
hist_all_aud_disagg <-
  read_excel(here::here("data/seabird_data/audubon_feeding/Audubon_historic_ARTE_COTE_ROST_feedings.xlsx")) %>%
  filter(Species == "COTE") %>% 
  dplyr::select(species = Species, 
                nFish = Number_of_Items, 
                location = Observer_Location,
                start = Date_Time_Start,
                end = Date_Time_End, 
                island = Island, 
                prey = Prey_Item,
                time_arrive = Time_Arrive,
                time_depart = Time_Depart,
                nest = Nest,
                provider = Provider) %>% 
  mutate(stint = paste(start, end, island), #no stint ID, so using start-end combination as identifier
         year = ifelse(!is.na(start), year(start), year(time_arrive)),
         month = ifelse(!is.na(start), month(start), month(time_arrive)),
         ym = yearmonth(paste(year, month, sep = "-")),
         nFish = ifelse(is.na(nFish), 1, nFish),
         nFish = ifelse(!is.numeric(nFish), 1, nFish),
         stint = str_remove_all(stint, "-"),
         loc_nest_stint = paste(location, nest, stint),
         time = as.numeric(difftime(end, start,
                  units = "hours")),
         loc_nest = paste(location, nest)) %>% 
  filter(between(time, 0, 10))


## Visuals-----
# start/end date by island
hist_all_aud_disagg %>% 
  group_by(island) %>% 
  dplyr::summarise(year_min = min(year),
                   year_max = max(year))

a <- ggplot(hist_all_aud_disagg) + 
  geom_histogram(aes(time)) +
  labs(title = "stint length",
       x = "Time (hrs)") +
  facet_wrap(~island) +
  dream::theme_fade() 

b <-
  ggplot(hist_all_aud_disagg) + 
  geom_histogram(aes(month)) +
  labs(title = "Observation month",
       x = "Month") +
  facet_wrap(~island) +
  dream::theme_fade()

hist_all_aud_summ <-
  hist_all_aud_disagg %>%
  group_by(year, month, ym, island, loc_nest, 
           stint, loc_nest_stint, time) %>% 
  dplyr::summarise(n_herring = sum(nFish[prey %in% c("R","r")]),
                   time = unique(time),
                   herring_per_hour = sum(nFish[prey %in% c("R","r")])/time)%>% 
  mutate(ds = 3)

# length(unique(hist_all_aud_disagg$loc_nest_stint))
# length(unique(hist_all_aud_summ$loc_nest_stint))

c <-
  ggplot(hist_all_aud_summ) + 
    geom_point(aes(y = herring_per_hour, x = ym),
               size = 0.1) +
    labs(title = "Herring per hour (n Herring/stint length)",
         x = "Year-month") +
    facet_wrap(~island) +
    dream::theme_fade() +
    theme(axis.text.x = element_text(size = 7))

ggsave(c, filename = here::here("figs/hph_aud_historic.png"),
       dpi = 300,
       width = 7,
       height = 5)

# Common terns (EB, METIN, PWNWR, SHIP) ----
cote_diet_mcinwr_disagg <- 
  read_excel(here::here("data/seabird_data/MCINWR_feeding/MCINWR_cote_feedings_from_sp.xlsx")) %>% 
  filter(Species == "COTE") %>% 
  dplyr::select(stint = StintID,
                species = Species, 
                nFish = Number_of_Items, 
                location = Observer_Location,
                start = Date_Time_Start,
                end = Date_Time_End, 
                island = Island,
                prey = Prey_Item,
                time_arrive = Time_Arrive,
                time_depart = Time_Depart,
                nest = Nest,
                provider = Provider) %>% 
  mutate(month = month(end),
         year = year(end),
         ym = yearmonth(paste(year, month, sep = "-")),
         nFish = ifelse(is.na(nFish), 1, nFish),
         stint = str_remove_all(stint, "-"),
         loc_nest_stint = paste(location, nest, stint),
         time = as.numeric(difftime(end, start,
                                    units = "hours")),
         loc_nest = paste(location, nest))

## Visuals----
# start/end date by island
cote_diet_mcinwr_disagg %>% 
  group_by(island) %>% 
  dplyr::summarise(year_min = min(year),
                   year_max = max(year))


a <-
  ggplot(cote_diet_mcinwr_disagg) + 
    geom_histogram(aes(time)) +
    labs(title = "stint length",
         x = "Time (hrs)") +
    facet_wrap(~island) +
    dream::theme_fade() 

b <- ggplot(cote_diet_mcinwr_disagg) + 
      geom_histogram(aes(month)) +
      labs(title = "Observation month",
           x = "Month") +
      facet_wrap(~island) +
      dream::theme_fade()

## Summarise----
cote_diet_mcinwr_summ <-
  cote_diet_mcinwr_disagg %>%
  group_by(year, month, ym, island, loc_nest, 
           stint, loc_nest_stint, time) %>% 
  dplyr::summarise(n_herring = sum(nFish[prey %in% c("R","r")]),
                   time = unique(time),
                   herring_per_hour = sum(nFish[prey %in% c("R","r")])/time)%>% 
  mutate(ds = 4)

# length(unique(cote_diet_mcinwr_disagg$loc_nest_stint))
# length(unique(cote_diet_mcinwr_summ$loc_nest_stint))

c <-
  ggplot(cote_diet_mcinwr_summ) + 
  geom_point(aes(y = herring_per_hour, x = ym),
             size = 0.1) +
  labs(title = "Herring per hour (n Herring/stint length)",
       x = "Year-month") +
  facet_wrap(~island) +
  dream::theme_fade() +
  theme(axis.text.x = element_text(size = 7))

ggsave(c, filename = here::here("figs/hph_mcinwr_modern.png"),
       dpi = 300,
       width = 7,
       height = 5)

# Historic MCINWR----
mcinwr_disagg <-
  read_csv(here::here("data/seabird_data/mcinwr_feeding/MCINWR_historic_ARTE_COTE_ROST_feedings.csv")) %>%
  filter(Species == "COTE") %>% 
  dplyr::rename(start = Date_Time_Start,
                end = Date_Time_End) %>% 
  mutate(start = as_datetime(start,
                             format = c("%m/%d/%y %H:%M",
                                        "%m/%d/%Y %H:%M")),
         end = as_datetime(end,
                           format = c("%m/%d/%y %H:%M",
                                      "%m/%d/%Y %H:%M")),
         year = year(start),
         month = month(start),
         time = as.numeric(end - start)/60/60,
         stint = paste(start, end)) %>% 
  filter(year != 1900) %>% 
  mutate(year = ifelse(year == 8, 2008, 
                       ifelse(year == 9, 2009,
                              ifelse(year == 10, 2010,
                                     year)))) %>% 
  dplyr::select(stint,
                species = Species, 
                nFish = Number_of_Items, 
                location = Observer_Location,
                start,
                end,
                island = Island,
                year,
                month,
                prey = Prey_Item_Code,
                time,
                time_arrive = Time_Arrive,
                time_depart = Time_Depart,
                nest = Nest,
                provider = Provider) %>% 
    mutate(month = month(end),
           year = year(end),
           ym = yearmonth(paste(year, month, sep = "-")),
           nFish = ifelse(is.na(nFish), 1, nFish),
           stint = str_remove_all(stint, "-"),
           loc_nest_stint = paste(location, nest, stint),
           time = as.numeric(difftime(end, start,
                                      units = "hours")),
           loc_nest = paste(location, nest))

## Visuals----
mcinwr_disagg %>% 
  group_by(island) %>% 
  dplyr::summarise(year_min = min(year),
                   year_max = max(year))

a <-
  ggplot(mcinwr_disagg) + 
    geom_histogram(aes(time)) +
    labs(title = "stint length",
         x = "Time (hrs)") +
    facet_wrap(~island) +
    dream::theme_fade() 

b <- 
  ggplot(mcinwr_disagg) + 
    geom_histogram(aes(month)) +
    labs(title = "Observation month",
         x = "Month") +
    facet_wrap(~island) +
    dream::theme_fade()


## Summarise----
mcinwr_summ <-
  mcinwr_disagg %>%
  group_by(year, month, ym, island, loc_nest, 
           stint, loc_nest_stint, time) %>% 
  dplyr::summarise(n_herring = sum(nFish[prey %in% c("R","r")]),
                   time = unique(time),
                   herring_per_hour = sum(nFish[prey %in% c("R","r")])/time) %>% 
  mutate(ds = 5)

# length(unique(mcinwr_disagg$loc_nest_stint))
# length(unique(mcinwr_summ$loc_nest_stint))


c <-
  ggplot(mcinwr_summ) + 
  geom_point(aes(y = herring_per_hour, x = ym),
             size = 0.1) +
  labs(title = "Herring per hour (n Herring/stint length)",
       x = "Year-month") +
  facet_wrap(~island) +
  dream::theme_fade() +
  theme(axis.text.x = element_text(size = 7))

ggsave(c, filename = here::here("figs/hph_mcinwr_historic.png"),
       dpi = 300,
       width = 6,
       height = 2.75)

tern_feeding_index_raw_simple <- 
  bind_rows(msi_summ,
            cote_diet_aud_summ,
            hist_all_aud_summ,
            cote_diet_mcinwr_summ,
            mcinwr_summ)  %>% 
  mutate(year_fac = factor(year),
         month_fac = factor(month),
         loc_nest_stint = factor(loc_nest_stint),
         island = factor(ifelse(island == "MRNWR", "MR",
                                ifelse(island == "PMNWR", "PMI",
                                       island))),
         ds = factor(ds)) 

save(tern_feeding_index_raw_simple,
     file = here::here("data/processed_herr_bird5_simple.rdata"))
