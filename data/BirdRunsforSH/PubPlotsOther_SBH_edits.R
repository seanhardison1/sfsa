#------------------------------------------------------------
# run_wham_code.R
#
# Example script using the 'here' package for relative paths.
#------------------------------------------------------------

# Load necessary packages
library(here)
library(wham)       # from devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel") or CRAN/pak
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

#--------------------------------------------------------------------
# Source your custom WHAM functions
#--------------------------------------------------------------------
# This file was originally at:
#  "C:\\Herring\\2025 Assessment RT\\Assessments\\WHAM\\GetWHAMrepofxns.R"
# Adjust the path inside here::here(...) to match your folder structure:
source(
  here::here("Assessments", "WHAM", "GetWHAMrepofxns.R")
)

#--------------------------------------------------------------------
# Point to your local WHAM repo directory (GitHub clone/fork)
#--------------------------------------------------------------------
# Originally: "C:\\Users\\jonathan.deroba\\Documents\\GitHub\\READ-PDB-COLLABORATIONS\\WHAM"
WHAMrepodirect <- here::here("READ-PDB-COLLABORATIONS", "WHAM")

# Use your custom function to source repository code
sourceWHAMrepo(direct = WHAMrepodirect)

#--------------------------------------------------------------------
# Load your base models
#--------------------------------------------------------------------
# runs = c("mm192_nonaa_nonecov","mm196_nonecov","mm196_birdecov")
# Adjust to match whichever runs you need. For illustration:
runs <- c("mm192_nonaa_nonecov", "mm211_nonaa_noecov_linebird", "mm211_nonaa_birdecov")

# Build paths to .rds files using here::here
mod.list <- file.path(
  here::here("data", "birdrunsforsh"),
  paste0(runs, ".rds")
)

# Read them in
mods <- lapply(mod.list, readRDS)
names(mods) <- runs

#--------------------------------------------------------------------
# Read the seabird index data
#--------------------------------------------------------------------
birdindex_file <- here::here("data",  "herring_provisioning_index1_27.csv")
birdindex <- read.csv(birdindex_file)

# Create capped CV (example as in your code)
birdindex$CV <- birdindex$se / birdindex$estimate
birdindex$CV <- ifelse(birdindex$CV > 0.9, 0.9, birdindex$CV)

#--------------------------------------------------------------------
# Example: Fit a model with new bird index data
#--------------------------------------------------------------------
# This replicates your modifications, e.g., for "mm211_nonaa_noecov_linebird"
oldinput <- mods[["mm192_nonaa_nonecov"]]$input
input <- oldinput

# Update indices with new bird index data
input$data$agg_indices[2:35, 7] <- birdindex$estimate
input$data$agg_index_sigma[2:35, 7] <- birdindex$se

mod.dir <- "mm211_nonaa_noecov_linebird"
write.dir <- here::here("Assessments", "WHAM", mod.dir)
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = TRUE)

# Fit the model
setwd(write.dir) # (Optional) set working directory
mod <- fit_wham(input, do.proj = FALSE, do.osa = TRUE, do.retro = TRUE, do.check = TRUE)

# Save and plot results
saveRDS(mod, file = file.path(write.dir, paste0(mod.dir, ".rds")))
plot_wham_output(mod, out.type = "png")

#--------------------------------------------------------------------
# Another example model with AR1 catchability, e.g. "mm211_nonaa_birdecov"
#--------------------------------------------------------------------
oldinput <- mods[["mm196_birdecov"]]$input
input <- oldinput
input$data$agg_indices[2:35, 7] <- birdindex$estimate
input$data$agg_index_sigma[2:35, 7] <- birdindex$se

mod.dir <- "mm211_nonaa_birdecov"
write.dir <- here::here("Assessments", "WHAM", mod.dir)
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = TRUE)

setwd(write.dir) # (Optional) set working directory
mod <- fit_wham(input, do.proj = FALSE, do.osa = TRUE, do.retro = TRUE, do.check = TRUE)
saveRDS(mod, file = file.path(write.dir, paste0(mod.dir, ".rds")))
plot_wham_output(mod, out.type = "png")

#--------------------------------------------------------------------
# Compare models
#--------------------------------------------------------------------
# Re-load or read again with updated runs:
runs <- c("mm192_nonaa_nonecov", "mm211_nonaa_noecov_linebird", "mm211_nonaa_birdecov")
mod.list <- file.path(
  here::here("data", "birdrunsforsh"),
  paste0(runs, ".rds")
)
mods <- lapply(mod.list, readRDS)
names(mods) <- runs

# Compare
mods.p <- mods
names(mods.p) <- c("base","linear","ar1")

# Store comparison output in 'Survey/seabird/compare_png' or similar
comparison_dir <- here::here("Survey", "seabird")
res <- compare_wham_models(
  mods.p, 
  fdir = comparison_dir,
  do.table = TRUE,
  calc.rho = FALSE,
  plot.opts = list(which = 1:2, years = 1987:2023)
)

#--------------------------------------------------------------------
# Extract and plot results needed
#--------------------------------------------------------------------
getres <- function(mod = NULL){
  # Check if log_index_resid has 7 columns
  if(ncol(mod$rep$log_index_resid) == 7){
    data.frame(
      "Year"         = mod$input$years_full,
      "LogRecruits"  = log(mod$rep$NAA[1,1,,1]),
      "LogResidual"  = mod$rep$log_index_resid[,7],
      "Catchability" = mod$rep$q[,7]
    )
  } else {
    data.frame(
      "Year"         = mod$input$years_full,
      "LogRecruits"  = log(mod$rep$NAA[1,1,,1]),
      "LogResidual"  = NA,
      "Catchability" = NA
    )
  }
}

base   <- getres(mods$mm192_nonaa_nonecov);          base$Model   <- "Base"
linear <- getres(mods$mm211_nonaa_noecov_linebird);  linear$Model <- "Linear"
ar1    <- getres(mods$mm211_nonaa_birdecov);         ar1$Model    <- "AR1"

df <- rbind(base, linear, ar1)
df$Model <- factor(df$Model, levels = c("Base","Linear","AR1"))
wham_out <- df
save(wham_out, file = here::here("data/wham_outputs.rdata"))

# Recruitment time series
recplot <-
  ggplot(df, aes(x = Year, y = exp(LogRecruits))) +
  geom_line(aes(color = Model, linewidth = Model)) +
  ylab("Age-1 Recruitment") +
  xlab("Year") +
  theme_bw() + 
  # ggsci::scale_color_d3() +
  scale_color_manual(values = c("Base" = "black", "Linear" = "purple", "AR1" = "darkorange")) +
  scale_linewidth_manual(values = c("Base" = 1, "Linear" = 0.5, "AR1" = 0.5)) +
  scale_y_continuous(labels = scales::scientific) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade() +
  theme(axis.title = element_text(size = 11))

ggsave(recplot,
       filename = here::here("figs/recplot.png"),
       dpi = 300,
       width = 6,
       height = 3.5)

# Bird index residuals vs. recruitment
residplot <-
  ggplot(df[df$Model != "Base", ], 
                    aes(x = LogRecruits, y = LogResidual)) +
  geom_hline(aes(yintercept = 0), linewidth = 0.5, alpha = 0.25, linetype = 2) +
  geom_point(size = 1) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  facet_wrap(~Model) +
  xlab("Log Recruitment") +
  ylab("Log Residual") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  dream::theme_fade() +
  theme(axis.title = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(residplot,
       filename = here::here("figs/residplot.png"),
       dpi = 300,
       width = 6,
       height = 3)

# Catchability vs. recruits from ecov fit
q_rec <-
  ggplot(df[df$Model == "AR1", ], 
                aes(x = exp(LogRecruits), y = Catchability)) +
  geom_point(size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  xlab("Age-1 Recruitment") +
  ylab("Catchability") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  dream::theme_fade() +
  theme(axis.title = element_text(size = 11))

ggsave(q_rec,
       filename = here::here("figs/q_rec.png"),
       dpi = 300,
       width = 6,
       height = 3.5)

#--------------------------------------------------------------------
# End of script
#--------------------------------------------------------------------

