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
library(gratia)

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

## load fitted model
load(here::here("data/fitted_dl_model_1014.rdata"))

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
p <- predict(m1_simple,
             newdata = ndf,
             model = NA,
             re_form = NULL,
             re_form_iid = NA,
             nsim = 1000,
             type = "link")

## annual index
x <- get_index_sims(p,
                    agg_function = function(x) median(exp(x)))

# read ASAP outputs
run5 <- dget(here::here("data/RUN5.RDAT"))

# extract numbers at age and join with seabird index
n_age <-
  run5$N.age %>% 
  as.data.frame() %>% 
  mutate(year = as.numeric(row.names(.))) %>%
  rename_at(.vars = 1:8, function(x){paste0("age",x)}) %>% 
  as_tibble() %>% 
  left_join(., x, by = "year") %>% 
  filter(year >= 1988) %>% 
  mutate(age2_next = lead(age2, 1),
         age3_next = lead(age3, 2),
         age4_next = lead(age4, 3),
         lest = log(est),
         lage1 = log(age1),
         sage1 = as.numeric(scale(age1))) 

# age1----
m1 <- lm(age1 ~ est, 
         data = n_age)
m1_sims <- simulateResiduals(m1, n = 1000)
plot(m1_sims)
summary(m1)
m1_pred <- predict(m1, se.fit = T)


m1_l <- lm(age1 ~ lest, 
         data = n_age)
m1l_sims <- simulateResiduals(m1_l, n = 1000)
plot(m1l_sims)
summary(m1_l)
# log transforming the bird index resolves the 
# residual patterning but doesn't markedly alter the outcome

# age2----
m2 <- lm(age2_next ~ est, 
         data = n_age)
m2_sims <- simulateResiduals(m2, n = 1000)
plot(m2_sims)
summary(m2)
m2_pred <- predict(m2, se.fit = T)


m2_l <- lm(age2_next ~ lest, 
           data = n_age)
m2l_sims <- simulateResiduals(m2_l, n = 1000)
plot(m2l_sims)
summary(m2_l)
# log transforming the bird index resolves the 
# residual patterning but doesn't markedly alter the outcome

# age3----
m3 <- lm(age3_next ~ est, 
         data = n_age)
m3_sims <- simulateResiduals(m3, n = 1000)
plot(m3_sims)
summary(m3)
m3_pred <- predict(m3, se.fit = T)

# age4----
m4 <- lm(age4_next ~ est, 
         data = n_age)
m4_sims <- simulateResiduals(m4, n = 1000)
plot(m4_sims)
summary(m4)
m4_pred <- predict(m4, se.fit = T)


pred_df <- tibble(m1_fit = m1_pred$fit,
                  m1_se = m1_pred$se.fit,
                  
                  m2_fit = c(m2_pred$fit, NA),
                  m2_se = c(m2_pred$se.fit, NA),
                  
                  m3_fit = c(m3_pred$fit, NA, NA),
                  m3_se = c(m3_pred$se.fit, NA, NA),
                  
                  m4_fit = c(m4_pred$fit, NA, NA, NA),
                  m4_se = c(m4_pred$se.fit, NA, NA, NA),
                  
                  est = n_age$est) %>% 
  left_join(., 
            n_age, by = "est") %>% 
  mutate(
    m1_low  = m1_fit - 1.96 * m1_se,
    m1_high = m1_fit + 1.96 * m1_se,
    
    m2_low  = m2_fit - 1.96 * m2_se,
    m2_high = m2_fit + 1.96 * m2_se,
    
    m3_low  = m3_fit - 1.96 * m3_se,
    m3_high = m3_fit + 1.96 * m3_se,
    
    m4_low  = m4_fit - 1.96 * m4_se,
    m4_high = m4_fit + 1.96 * m4_se
  )


# visualize-----
a <- ggplot(pred_df) +
  geom_point(aes(y = age1, x = est)) +
  geom_line(aes(y = m1_fit, x = est)) + 
  geom_ribbon(aes(ymin = m1_low, ymax = m1_high,
                  x = est), alpha = 0.25) + 
  labs(y = "Age-1 herring abundance",
       x = "Seabird provisioning index (herring per hour)") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade()

b <- ggplot(pred_df %>% 
              filter(!is.na(age2_next))) +
  geom_point(aes(y = age2_next, x = est)) +
  geom_line(aes(y = m2_fit, x = est)) + 
  geom_ribbon(aes(ymin = m2_low, ymax = m2_high,
                  x = est), alpha = 0.25) + 
  labs(y = "Age-2 herring abundance (1-yr lead)",
       x = "Seabird provisioning index (herring per hour)") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade()


c <- ggplot(pred_df %>% 
              filter(!is.na(age3_next))) +
  geom_point(aes(y = age3_next, x = est)) +
  geom_line(aes(y = m3_fit, x = est)) + 
  geom_ribbon(aes(ymin = m3_low, ymax = m3_high,
                  x = est), alpha = 0.25) + 
  labs(y = "Age-3 herring abundance (2-yr lead)",
       x = "Seabird provisioning index (herring per hour)") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade()


d <- ggplot(pred_df %>% 
              filter(!is.na(age4_next))) +
  geom_point(aes(y = age4_next, x = est)) +
  geom_line(aes(y = m4_fit, x = est)) + 
  geom_ribbon(aes(ymin = m4_low, ymax = m4_high,
                  x = est), alpha = 0.25) + 
  labs(y = "Age-4 herring abundance (3-yr lead)",
       x = "Seabird provisioning index (herring per hour)") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade() 



plt_out <- 
  a + b + 
  c + d + 
  plot_layout(nrow = 2) + 
  plot_annotation(
    tag_level = "a",
    tag_prefix = "(",
    tag_suffix = ")") & 
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

ggsave(plt_out,
       filename = here::here("figs/asap_herring_comps.png"),
       dpi = 300,
       width = 8,
       height = 7)
                  