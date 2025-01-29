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
library(ggtext)
library(readxl)

## load fitted model
load(here::here("data/fitted_dl_model_1014.rdata"))

## data for prediction-----
load(here::here("data/st_model_pred_sims.rdata"))

## annual index
x <- get_index_sims(p,
                    agg_function = function(x) median(exp(x)))

# read WHAM base model outputs
load(here::here("data/wham_outputs.rdata"))

# extract numbers at age and join with seabird index
n_age <- wham_out %>% 
  filter(Model == "Base") %>% 
  mutate(age1 = exp(LogRecruits)) %>% 
  dplyr::rename(year = Year) %>% 
  left_join(., x, by = "year") %>% 
  filter(!is.na(est))

# age1----
m1 <- lm(age1 ~ est, 
         data = n_age)
m1_sims <- simulateResiduals(m1)
plot(m1_sims)
summary(m1)
m1_pred <- predict(m1, se.fit = T)

# tidy summaries
bind_rows(
  broom::tidy(m1) %>% 
    mutate(model = "age1",
           r2 = summary(m1)$r.square)
) %>% 
  filter(term == "est")

pred_df <- tibble(m1_fit = m1_pred$fit,
                  m1_se = m1_pred$se.fit,
                  est = n_age$est) %>% 
  left_join(., 
            n_age, by = "est") %>% 
  mutate(
    m1_low  = m1_fit - 1.96 * m1_se,
    m1_high = m1_fit + 1.96 * m1_se,
  )

r2_df <- tibble(x = 0.2,
                y = pred_df %>% filter(age1 == max(age1)) %>% pull(age1),
                r2 = paste0("R<sup>2</sup> = ",round(summary(m1)$r.square,3)),
                model = 1)

# visualize-----
plt_out <-
  ggplot(pred_df) +
  geom_point(aes(y = age1, x = est)) +
  geom_line(aes(y = m1_fit, x = est)) + 
  geom_ribbon(aes(ymin = m1_low, ymax = m1_high,
                  x = est), alpha = 0.25) +
  geom_richtext(data = r2_df %>% 
              filter(model == 1),
            aes(x = x, y = y,
                label = r2),
            label.colour = "transparent") +
  scale_y_continuous(labels = scales::scientific) +
  labs(y = "Age-1 Recruitment",
       x = "Herring provisioning index (herring per hour)") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

ggsave(plt_out,
       filename = here::here("figs/wham_herring_comps.png"),
       dpi = 300,
       width = 6,
       height = 5)
                  