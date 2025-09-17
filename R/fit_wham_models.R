# Use this version of WHAM
# pak::pkg_install("timjmiller/wham@e1b6a70")
library(wham)
library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)

## function to extract results from fitted models-----
get_res <- function(mod=NULL){
  if(ncol(mod$rep$log_index_resid)==7){
    res=data.frame("Year"=mod$input$years_full,
                   "LogRecruits"=log(mod$rep$NAA[1,1,,1]),
                   "LogResidual"=mod$rep$log_index_resid[,7],
                   "Catchability"=mod$rep$q[,7])
  } else {
    res=data.frame("Year"=mod$input$years_full,
                   "LogRecruits"=log(mod$rep$NAA[1,1,,1]),
                   "LogResidual"=NA,
                   "Catchability"=NA)  
  }
  return(res)
}


## model fitting function-----
fit_model <- 
  function(hpi = NULL,
           mod = mods[[2]], 
           output_dir = here::here("data","output"),
           plot_output = F,
           add_bird = T){
    if (add_bird){
      hpi$CV=hpi$se/hpi$estimate
      hpi$CV=ifelse(hpi$CV>0.9, 0.9, hpi$CV)
      input <- mod$input
      input$data$agg_indices[2:35,7] <- hpi$estimate
      input$data$agg_index_sigma[2:35,7] <- hpi$se
    } else {
      input <- mod$input
    }

    
    # re-fit the model with the selected seabird index
    mod <- fit_wham(input, 
                    do.proj=FALSE,
                    do.osa = TRUE,
                    do.retro = TRUE,
                    do.check = T)
    if(!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    }
    saveRDS(mod,
            file=file.path(output_dir,
                           paste0(mod$model_name,".rds")))
    if (plot_output){
      plot_wham_output(mod,
                       out.type="png",
                       dir.main = output_dir)
    }
    
    return(mod)
  }


# load everything if not running the models
process <- T
if (process){
  
  # import data-----
  ## iterations on herring provisioning index---
  load(here::here("data/herring_index_iterations.rdata"))
  
  ## import herring provisioning index (10 knots)----
  birdindex_file <- here::here("data",  "herring_provisioning_index1_27.csv")
  birdindex <- read.csv(birdindex_file)
  
  ## fitted models (main results)----
  runs <- c("Run0_nonaa", # base
            "nonaa_nonecov_Run0", # linear
            "nonaa_birdecov_Run0") # AR!
  mod.list <- file.path(
    here::here("data", "birdrunsforsh"),
    paste0(runs, ".rds")
  )
  mods <- lapply(mod.list, readRDS)
  names(mods)=runs
  
  ## Sensitivity runs (not run here)-----
  wham_sens_runs=c("mm192_nonaa_nonecov", # same as Base
                   "mm211_nonaa_noecov_btsbad", # linear
                   "mm211_nonaa_ecov_btsbad") #AR1
  wham_sens_mod.list <- file.path(
    here::here("data", "birdrunsforsh"),
    paste0(wham_sens_runs, ".rds")
  )
  wham_sens_mods <- lapply(wham_sens_mod.list, readRDS)
  names(wham_sens_mods)=wham_sens_runs
  
  ## Initial scenario runs------
  
  # Fit models (i.e. re-fit models for reproducibility). The fit_model function
  # drops in the bird index of choice and re-fits the model. The default index is the
  # one estimated with 10 knots.
  
  ## base model (no herring provisioning index)----
  base_mod <- 
    fit_model(mod = mods[[1]], #"Run0_nonaa"
              output_dir = here::here("data","Run0_nonaa"),
              add_bird = FALSE)
  
  ## assumes linear relationship between recruitment and herring provisioning index----
  lin_mod <- 
    fit_model(hpi = birdindex,
              mod = mods[[2]], #"nonaa_nonecov_Run0"
              output_dir = here::here("data","nonaa_nonecov_Run0"))
  
  ## Now allows for AR1 catchability-----
  ar1_mod <- 
    fit_model(hpi = birdindex,
              mod = mods[[3]], #"nonaa_birdecov_Run0"
              output_dir = here::here("data","nonaa_birdecov_Run0"))

  # collect recruitment and residual outputs----
  base_recruitment <- get_res(mod = base_mod)
  lin_recruitment <- get_res(mod = lin_mod)
  ar1_recruitment <- get_res(mod = ar1_mod)
  base_sens <- get_res(mod = wham_sens_mods[[1]])
  lin_sens <- get_res(mod = wham_sens_mods[[2]])
  ar1_sens <- get_res(mod = wham_sens_mods[[3]])
  
  df <- bind_rows(base_recruitment %>% mutate(`WHAM Scenario` = "Base"),
                  lin_recruitment %>% mutate(`WHAM Scenario` = "Linear"),
                  ar1_recruitment %>% mutate(`WHAM Scenario` = "AR1"),
                  lin_sens %>% mutate(`WHAM Scenario` = "Linear sensitivity"),
                  ar1_sens %>% mutate(`WHAM Scenario` = "AR1 sensitivity")) %>% 
    mutate(`WHAM Scenario` = factor(`WHAM Scenario`, levels = c("Base",
                                                                "Linear",
                                                                "AR1",
                                                                "Linear sensitivity",
                                                                "AR1 sensitivity")))
  
  # eight knots
  k8_wham_mod <- 
    fit_model(hpi = x_knots %>% dplyr::rename(estimate = est) %>% filter(`N knots` == 8),
              mod = mods[[2]],
              output_dir = here::here("data","nonaa_nonecov_Run0_8knots"),
              plot_output = FALSE)
  
  # ten knots
  k10_wham_mod <- 
    fit_model(hpi = x_knots %>% dplyr::rename(estimate = est) %>% filter(`N knots` == 10),
              mod = mods[[2]], 
              output_dir = here::here("data","nonaa_nonecov_Run0_10knots"),
              plot_output = FALSE)
  
  # 12 knots
  k12_wham_mod <- 
    fit_model(hpi = x_knots %>% dplyr::rename(estimate = est) %>% filter(`N knots` == 12),
              mod = mods[[2]], 
              output_dir = here::here("data","nonaa_nonecov_Run0_12knots"),
              plot_output = FALSE)
  
  # NE island grouping
  ne_wham_mod <- 
    fit_model(hpi = x_spat_sub %>% dplyr::rename(estimate = est) %>% filter(`Pred. region` == "NE islands"),
              mod = mods[[2]], 
              output_dir = here::here("data","nonaa_nonecov_Run0_NE_islands"),
              plot_output = FALSE)
  
  # SW island grouping
  sw_wham_mod <- 
    fit_model(hpi = x_spat_sub %>% dplyr::rename(estimate = est) %>% filter(`Pred. region` == "SW islands"),
              mod = mods[[2]], 
              output_dir = here::here("data","nonaa_nonecov_Run0_SW_islands"),
              plot_output = FALSE)
  
  
  k8_wham_mod_recruitment <- get_res(mod = k8_wham_mod)
  k10_wham_mod_recruitment <- get_res(mod = k10_wham_mod)
  k12_wham_mod_recruitment <- get_res(mod = k12_wham_mod)
  ne_wham_mod_recruitment <- get_res(mod = ne_wham_mod)
  sw_wham_mod_recruitment <- get_res(mod = sw_wham_mod)
  
  
  # collect recruitment and residual outputs----
  herring_ind_sens_results <- 
    bind_rows(k8_wham_mod_recruitment %>% mutate(`delta Model` = "8 knots"),
              k10_wham_mod_recruitment %>% mutate(`delta Model` = "10 knots"),
              k12_wham_mod_recruitment %>% mutate(`delta Model` = "12 knots"),
              ne_wham_mod_recruitment %>% mutate(`delta Model` = "NE islands"),
              sw_wham_mod_recruitment %>% mutate(`delta Model` = "SW islands"),
              base_recruitment  %>% mutate(`delta Model` = "None")) %>% 
    mutate(LogResidual = ifelse(LogResidual == 0, NA, LogResidual))
  
  # fit regression models to evaluate sensitivity-----
  herring_ind_sens_results %>% 
    filter(`delta Model` != "None") %>% 
    split(.$`delta Model`) %>% 
    map(\(df) lm(LogResidual ~ LogRecruits, data = df)) |>
    map(summary)
  
  # save outputs----
  save(base_recruitment,
       lin_recruitment,
       ar1_recruitment,
       base_sens,
       lin_sens,
       ar1_sens,
       df,
       herring_ind_sens_results,
       k8_wham_mod_recruitment,
       k10_wham_mod_recruitment,
       k12_wham_mod_recruitment,
       ne_wham_mod_recruitment,
       sw_wham_mod_recruitment,
    file = here::here("data/index_sensitivity_analyses.rdata"))
} else {
  load(here::here("data/index_sensitivity_analyses.rdata"))
}


# Visualize recruitment sensivity -----
p <- c("AR1" = "purple", "Base" = "darkorange", "Linear" = "darkblue",
       "Linear sensitivity" = "darkblue",
       "AR1 sensitivity" = "purple")

lt <- c("AR1" = 1, "Base" = 1, "Linear" = 1,
       "Linear sensitivity" = 2,
       "AR1 sensitivity" = 2)

a <- 
  ggplot(df %>% 
           filter(`WHAM Scenario` %in% c("Base", "Linear", "AR1"))) +
    geom_line(aes(y = exp(LogRecruits), 
                  x = Year, 
                  color = `WHAM Scenario`,
                  linetype = `WHAM Scenario`)) +
    scale_color_manual(values = p) +
    scale_linetype_manual(values = lt) +
    scale_y_continuous(labels = scales::label_scientific(),
                       limits = c(0, 1.25e7)) +
    labs(y = "Recruitment (age-1 abundance)",
         x = "Year") +
    theme_bw()

b <- 
  ggplot(df %>% 
           filter(`WHAM Scenario` %in% c("Base", 
                                         "Linear sensitivity",
                                         "AR1 sensitivity"))) +
    geom_line(aes(y = exp(LogRecruits), 
                  x = Year, 
                  color = `WHAM Scenario`,
                  linetype = `WHAM Scenario`)) +
    scale_color_manual(values = p) +
    scale_linetype_manual(values = lt) +
    scale_y_continuous(labels = scales::label_scientific(),
                       limits = c(0, 1.25e7)) +
    labs(y = "Recruitment (age-1 abundance)",
         x = "Year") +
    theme_bw() +
    theme(axis.title.y = element_blank())

rec_plot <- 
  (a + b) + 
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(tag_levels = "a",
                    tag_prefix = "(",
                    tag_suffix = ")") &
    theme(legend.position = "bottom")

ggsave(rec_plot, #post-processing in inkscape to make pretty
       filename=here::here("figs/rec_plot.svg"),
       width = 9,
       height = 4)

# Visualize residual sensitivity to index specifications------ 
index_sens_plot <- 
  ggplot(herring_ind_sens_results %>% 
           filter(`delta Model` != "None")) +
    geom_point(aes(y = LogResidual, x = LogRecruits, color = `delta Model`, group = `delta Model`)) +
    geom_smooth(aes(y = LogResidual, x = LogRecruits, color = `delta Model`),
                method = "lm", se = F) +
    ggsci::scale_color_uchicago() +
    theme_bw() +
    xlab("Recruitment (age-1 abundance, log)")+
    ylab("Herring provisioning index residual (log)")  +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.text = element_text(size = 12))


ggsave(index_sens_plot,
       filename=here::here("figs/index_sens_plot.png"),
       dpi = 300,
       width = 6,
       height = 4)


df2 <- bind_rows(ar1_recruitment %>% mutate(Model = "AR1"),
                 base_recruitment %>% mutate(Model = "Base"),
                 lin_recruitment %>% mutate(Model = "Linear")) %>% 
  mutate(recruits = exp(LogRecruits),
         residuals = exp(LogResidual),
         Model = factor(Model, levels = c("Linear", "AR1", "Base")))

# Visualize log-residual ~ log-recruitment-----
residplot=ggplot(df2[df2$Model != "Base",])+
  geom_point(aes(x=LogRecruits,y=LogResidual),size=1)+
  geom_smooth(aes(x=LogRecruits,y=LogResidual),method="lm",se=FALSE, color = "black")+
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Model)+
  theme_bw() +
  xlab("Recruitment (age-1 abundance, log)")+
  ylab("Herring provisioning index residual (log)") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
residplot

ggsave(residplot,
       filename=here::here("figs/residual_plots.png"),
       dpi = 300,
       width = 7,
       height = 4)

#Visualize q ~ recruits from AR1 fit-----
q_rec=ggplot(df2[df2$Model=="AR1",])+
  geom_point(aes(x=exp(LogRecruits),y=Catchability),size=1)+
  geom_smooth(aes(x=exp(LogRecruits),y=Catchability),method=loess,se=FALSE, color = "black")+
  theme_bw() +
  xlab("Recruitment (age-1 abundance)")+
  ylab("Catchability") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))

ggsave(q_rec,
       filename=here::here("figs/catchability.png"),
       dpi = 300,
       width = 5,
       height = 4)
