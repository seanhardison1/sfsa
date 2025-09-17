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
library(ggsci)
library(lubridate)
library(readxl)
library(ggtext)
library(here)

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

# calculate mean inter-island distance-----
dist_mat <- sf::st_distance(locs)
row.names(dist_mat) <- locs$island_long
colnames(dist_mat) <- locs$island_long
dist_mat
dist_vec <- dist_mat[lower.tri(dist_mat, diag = FALSE)]
(mean_dist <- mean(as.numeric(dist_vec)))
(range_dist <- range(as.numeric(dist_vec)))

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

## make mesh (10 knots)----
mesh <- make_mesh(df, c("lon", "lat"), n_knots = 10)
plot(mesh)

# fit model using 10 knots-----
# you will need to run the scripts where we fit the spatiotemporal models rather than load them
# because the output files "fitted_dl_model_1014.rdata" and "fitted_dl_model_12k.rdata" 
# are too large to push to github. Email seanhardison@gmail.com if I can help troubleshoot
process <- F
if (process){
  # fit model
  m1_simple <- 
    sdmTMB(herring_per_hour ~ 0 + 
             year_fac + 
             (1|ym_fac) +
             (1|loc_nest_year),
           spatial = "on",
           spatiotemporal = "ar1",
           time = "year",
           mesh = mesh,
           data = df,
           family = delta_lognormal()
    )
  
  save(m1_simple, file = here::here("data/fitted_dl_model_1014.rdata"))
  
  ## MCMC residuals----
  # samp1 <- predict_mle_mcmc(m1_simple, model = 1, mcmc_iter = 5000)
  samp1 <- predict_mle_mcmc(m1_simple, model = 1, mcmc_iter = 501, mcmc_warmup = 500)
  r1 <- residuals(m1_simple, type = "mle-mcmc", mcmc_samples = samp1)
  png(filename = here::here("figs/d_ln_prob_mod_resids_simple.png"))
  qqnorm(r1, main = "Herring model predation residuals:\ndelta_lognormal probability model")
  qqline(r1)
  dev.off()
  
  samp2 <- predict_mle_mcmc(m1_simple, model = 2, mcmc_iter = 501, mcmc_warmup = 500)
  r2 <- residuals(m1_simple, type = "mle-mcmc", mcmc_samples = samp2)
  png(filename = here::here("figs/d_ln_catch_mod_resids_simple.png"))
  qqnorm(r2, main = "Herring model predation residuals:\ndelta_lognormal positive catch model")
  qqline(r2)
  dev.off()
  save(r1, r2, file = here::here("data/catch_mod_mcmc_resids_simple.rdata"))
} else {
  load(here::here("data/fitted_dl_model_1014.rdata"))
  load(here::here("data/catch_mod_mcmc_resids_simple.rdata"))
  qqnorm(r1, main = "Herring model predation residuals:\ndelta_lognormal probability model")
  qqline(r1)
  qqnorm(r2, main = "Herring model predation residuals:\ndelta_lognormal positive catch model")
  qqline(r2)
}

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

if (process){
  ## island-level predictions-----
  p <- predict(m1_simple,
               newdata = ndf,
               model = NA,
               re_form = NULL,
               re_form_iid = NA,
               nsim = 1000,
               type = "link")
  save(p, file = here::here("data/st_model_pred_sims.rdata"))
} else {
  load(here::here("data/st_model_pred_sims.rdata"))
}


## process predictions-----
isl_preds <-
  p %>%
  bind_cols(., ndf %>%
              dplyr::select(-loc_nest_year,-ym_fac)) %>%
  gather(sim, value, -year_fac:-island_long) %>%
  group_by_at(vars(year_fac:island_long)) %>%
  dplyr::summarise(
    est_link = median(value, na.rm = TRUE),
    q025_link = stats::quantile(value, 0.025, na.rm = TRUE),
    q975_link = stats::quantile(value, 0.975, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    fit_resp = exp(est_link),
    lwr_resp = exp(q025_link),
    upr_resp = exp(q975_link)
  )

## visualize----
isl_preds_plt <-
  isl_preds %>%
  ggplot() +
    geom_ribbon(aes(ymax = upr_resp, ymin = lwr_resp, x = year), alpha = 0.1) +
    geom_point(data = df, aes(y = herring_per_hour,
                              x = year),
               size = 0.01) +
    geom_line(aes(y = fit_resp, x = year),
              color = "darkorange", linewidth = 0.75) +
    facet_wrap(~island_long, scales =
                 "free_y") +
    labs(x = "Year",
         y = "Herring provisioning index, <i>H</i> (herring hour<sup>-1</sup>)") +
    dream::theme_fade() +
    theme(axis.title.y = element_markdown(size = 13),
          axis.title.x = element_markdown(size = 13))

ggsave(isl_preds_plt,
       filename = here::here("figs/isl_preds_plt.png"),
       dpi = 300,
       height = 6,
       width = 8)

## across-island predictions----
x <- get_index_sims(p,
                    agg_function = function(x) median(exp(x)))

## write out index-----
x2 <- x %>% 
  dplyr::mutate(log_se = log(se)) %>% 
  dplyr::rename(estimate = est,
                log_estimate = log_est
                )
write_csv(x2, here::here("data/herring_provisioning_index1_27.csv"))

## visualize-----
reg_pred_plt <-
  ggplot(x) +
    geom_errorbar(aes(x = year, ymin = lwr, ymax = upr), alpha = 0.5, 
                  width = 0) +
    geom_line(aes(y = est, x = year), linewidth = 0.75) +
    geom_point(aes(y = est, x = year)) +
    labs(x = "Year",
         y = "Herring provisioning index, <i>H</i> (herring hour<sup>-1</sup>)") +
      dream::theme_fade() +
      theme(axis.title = element_markdown(size = 12))

reg_pred_plt

ggsave(reg_pred_plt,
       filename = here::here("figs/reg_pred_plt.png"),
       dpi = 300,
       height = 3,
       width = 5)

# fit model using 12 and 8 knots-----
mesh_12k <- make_mesh(df, c("lon", "lat"), n_knots = 12)
mesh_8k <- make_mesh(df, c("lon", "lat"), n_knots = 8)

if (process){
  # 8 knot model----
  m1_8k <- 
    sdmTMB(herring_per_hour ~ 0 + 
             year_fac + 
             (1|ym_fac) +
             (1|loc_nest_year),
           spatial = "on",
           spatiotemporal = "ar1",
           time = "year",
           mesh = mesh_8k,
           data = df,
           family = delta_lognormal()
    )
  
  
  # 12 knot model-----
  m1_12k <- 
    sdmTMB(herring_per_hour ~ 0 + 
             year_fac + 
             (1|ym_fac) +
             (1|loc_nest_year),
           spatial = "on",
           spatiotemporal = "ar1",
           time = "year",
           mesh = mesh_12k,
           data = df,
           family = delta_lognormal()
    )
  
  ## island-level predictions-----
  p_8k <- predict(m1_8k,
                   newdata = ndf,
                   model = NA,
                   re_form = NULL,
                   re_form_iid = NA,
                   nsim = 1000,
                   type = "link")
  
  p_12k <- predict(m1_12k,
               newdata = ndf,
               model = NA,
               re_form = NULL,
               re_form_iid = NA,
               nsim = 1000,
               type = "link")
  
  # process predictions-----
  isl_preds_8k <-
    p_8k %>%
    bind_cols(., ndf %>%
                dplyr::select(-loc_nest_year,-ym_fac)) %>%
    gather(sim, value, -year_fac:-island_long) %>%
    group_by_at(vars(year_fac:island_long)) %>%
    dplyr::summarise(
      est_link = median(value, na.rm = TRUE),
      q025_link = stats::quantile(value, 0.025, na.rm = TRUE),
      q975_link = stats::quantile(value, 0.975, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      fit_resp = exp(est_link),
      lwr_resp = exp(q025_link),
      upr_resp = exp(q975_link)
    )
  
  
  isl_preds_12k <-
    p_12k %>%
    bind_cols(., ndf %>%
                dplyr::select(-loc_nest_year,-ym_fac)) %>%
    gather(sim, value, -year_fac:-island_long) %>%
    group_by_at(vars(year_fac:island_long)) %>%
    dplyr::summarise(
      est_link = median(value, na.rm = TRUE),
      q025_link = stats::quantile(value, 0.025, na.rm = TRUE),
      q975_link = stats::quantile(value, 0.975, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      fit_resp = exp(est_link),
      lwr_resp = exp(q025_link),
      upr_resp = exp(q975_link)
    )
  
  # across-island predictions----
  x_12k <- get_index_sims(p_12k,
                      agg_function = function(x) median(exp(x)))
  
  x_8k <- get_index_sims(p_8k,
                          agg_function = function(x) median(exp(x)))
  

  save(m1_12k,
       m1_8k,
       isl_preds_12k,
       isl_preds_8k,
       x_12k,
       x_8k,
       file = here::here("data/fitted_dl_model_12k.rdata"))
} else {
  load(here::here("data/fitted_dl_model_12k.rdata"))
}

# compare preds with 10 vs 12 knots----
x_knots <- 
  bind_rows(
    x %>% 
      mutate(`N knots` = "10"),
    x_12k %>% 
      mutate(`N knots` = "12"),
    x_8k %>% 
      mutate(`N knots` = "8")
  )

knot_compare <-
  ggplot() + 
  geom_ribbon(data = x_knots,
              aes(x = year, ymin = lwr, ymax = upr,
                  fill = `N knots`), 
              alpha = 0.1) +
  geom_line(data = x_knots,
            aes(y = est, x = year,
                color = `N knots`)) +
  geom_point(data = x_knots,
             aes(y = est, x = year,
                 color = `N knots`)) +

  scale_x_continuous(expand = c(0.01, 0.01)) +
  ggsci::scale_color_d3() + 
  ggsci::scale_fill_d3() + 
  labs(x = "Year",
       y = "Herring provisioning index, <i>H</i> (herring hour<sup>-1</sup>)") +
  dream::theme_fade() +
  theme(axis.title = element_markdown(size = 12),
        legend.position = c(0.8,0.8))
knot_compare

ggsave(knot_compare,
       filename = here::here("figs/knot_compare.png"),
       dpi = 300,
       width = 6,
       height = 3.5)

# predict at SW islands only----
ndf_sw_subset <- ndf %>% 
  filter(lat < 4860)

## island-level predictions-----
p_sw_subset <- predict(m1_simple,
                 newdata = ndf_sw_subset,
                 model = NA,
                 re_form = NULL,
                 re_form_iid = NA,
                 nsim = 1000,
                 type = "link")

## process predictions-----
isl_preds_sw_subset <-
  p_sw_subset %>%
  bind_cols(., ndf_sw_subset %>%
              dplyr::select(-loc_nest_year,-ym_fac)) %>%
  gather(sim, value, -year_fac:-island_long) %>%
  group_by_at(vars(year_fac:island_long)) %>%
  dplyr::summarise(
    est_link = median(value, na.rm = TRUE),
    q025_link = stats::quantile(value, 0.025, na.rm = TRUE),
    q975_link = stats::quantile(value, 0.975, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    fit_resp = exp(est_link),
    lwr_resp = exp(q025_link),
    upr_resp = exp(q975_link)
  )

## across-island predictions----
x_sw_subset <- get_index_sims(p_sw_subset,
                        agg_function = function(x) median(exp(x)))

# predict at ne islands only----
ndf_ne_subset <- ndf %>% 
  filter(lat > 4860)

## island-level predictions-----
p_ne_subset <- predict(m1_simple,
                       newdata = ndf_ne_subset,
                       model = NA,
                       re_form = NULL,
                       re_form_iid = NA,
                       nsim = 1000,
                       type = "link")

## process predictions-----
isl_preds_ne_subset <-
  p_ne_subset %>%
  bind_cols(., ndf_ne_subset %>%
              dplyr::select(-loc_nest_year,-ym_fac)) %>%
  gather(sim, value, -year_fac:-island_long) %>%
  group_by_at(vars(year_fac:island_long)) %>%
  dplyr::summarise(
    est_link = median(value, na.rm = TRUE),
    q025_link = stats::quantile(value, 0.025, na.rm = TRUE),
    q975_link = stats::quantile(value, 0.975, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    fit_resp = exp(est_link),
    lwr_resp = exp(q025_link),
    upr_resp = exp(q975_link)
  )

## across-island predictions----
x_ne_subset <- get_index_sims(p_ne_subset,
                              agg_function = function(x) median(exp(x)))


# compare spatial predictions------
x_spat_sub <- 
  bind_rows(
    x %>% 
      mutate(`Pred. region` = "All islands"),
    x_sw_subset %>% 
      mutate(`Pred. region` = "SW islands"),
    x_ne_subset %>% 
      mutate(`Pred. region` = "NE islands")
  ) %>% 
  mutate(`Pred. region` = factor(`Pred. region`, levels = c( "All islands","NE islands","SW islands")))

pal <- c("SW islands"= "#FFC20A", 
         "NE islands" = "#0C7BDC", 
         "All islands" = "black")

spat_compare <- 
  ggplot() + 
  geom_point(data = x_spat_sub %>% filter(`Pred. region` %in% "NE islands"),
             aes(y = est, x = year - 0.15,
                 color = `Pred. region`),
             alpha = 0.5) +
  geom_errorbar(data = x_spat_sub %>% filter(`Pred. region` %in% "NE islands"),
              aes(x = year - 0.15, ymin = lwr, ymax = upr,
                  color = `Pred. region`, width = 0),
              alpha = 0.5) +
  
  geom_point(data = x_spat_sub %>% filter(`Pred. region` %in% "SW islands"),
             aes(y = est, x = year + 0.15,
                 color = `Pred. region`),
             alpha = 0.5) +
  geom_errorbar(data = x_spat_sub %>% filter(`Pred. region` %in% "SW islands"),
                aes(x = year + 0.15, ymin = lwr, ymax = upr,
                    color = `Pred. region`, width = 0),
                alpha = 0.5) +
  
  
  geom_point(data = x_spat_sub %>% filter(`Pred. region` %in% "All islands"),
             aes(y = est, x = year,
                 color = `Pred. region`)) +
  geom_line(data = x_spat_sub %>% filter(`Pred. region` %in% "All islands"),
            aes(y = est, x = year,
                color = `Pred. region`)) +
  geom_errorbar(data = x_spat_sub %>% filter(`Pred. region` %in% "All islands"),
                aes(x = year, ymin = lwr, ymax = upr,
                    color = `Pred. region`, width = 0)) +
  
  scale_x_continuous(expand = c(0.01, 0.01),
                     labels = c(seq(1990, 2020, by  = 5)),
                     breaks = c(seq(1990, 2020, by  = 5))) +
  scale_color_manual(values = pal) +
  # scale_fill_manual(values = pal) +
  # ggsci::scale_color_d3() + 
  # ggsci::scale_fill_d3() + 
  labs(x = "Year",
       y = "Herring provisioning index, <i>H</i> (herring hour<sup>-1</sup>)") +
  theme_bw() +
  theme(axis.title = element_markdown(size = 12),
        legend.position = c(0.8,0.8))


ggsave(spat_compare,
       filename = here::here("figs/spat_compare.png"),
       dpi = 300,
       width = 7.5,
       height = 4)

# trend in herring provisioning index-----
xt_mod <- lm(est ~ year,data = x)
summary(xt_mod)
acf(residuals(xt_mod))
xt_mod_sims <- simulateResiduals(xt_mod, n = 1000)
plot(xt_mod_sims)


# save sensitivity outputs-----
save(x_spat_sub,
     x_knots,
     file = here::here("data/herring_index_iterations.rdata"))
