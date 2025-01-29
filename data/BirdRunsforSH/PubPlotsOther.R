source("C:\\Herring\\2025 Assessment RT\\Assessments\\WHAM\\GetWHAMrepofxns.R")

WHAMrepodirect="C:\\Users\\jonathan.deroba\\Documents\\GitHub\\READ-PDB-COLLABORATIONS\\WHAM"
sourceWHAMrepo(direct=WHAMrepodirect)

##############
devtools::install_github("timjmiller/wham", dependencies=TRUE)
devtools::install_github("timjmiller/wham", dependencies=TRUE,ref="devel")
pak::pkg_install("timjmiller/wham@devel")

library(wham)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
##############
##***These are obsolete as Sean provided new index Jan 27, 2025. These re re-run below.
runs=c("mm192_nonaa_nonecov","mm196_nonecov","mm196_birdecov") #mm196 is mm192 but with linear seabird index included and NAA off so model can't tradeoff fits to data.
mod.list <- file.path(paste(paste("C:/Herring/2025 Assessment RT/Assessments/WHAM",runs,sep="/"),paste0(runs,".rds"),sep="/"))
mods <- lapply(mod.list, readRDS)
names(mods)=runs

##############
#re-do with herring provisioning index Sean provided January 27.
birdindex=read.csv(file.path("C:/Herring/2025 Assessment RT/Survey/seabird/herring_provisioning_index1_27.csv"))
birdindex$CV=birdindex$se/birdindex$estimate
birdindex$CV=ifelse(birdindex$CV>0.9, 0.9, birdindex$CV)

#This assumes linear relationship between index and recruitment; No NAA RE and
#the dummy covariate does not effect the catchability in this run
oldinput=mods$mm196_nonecov$input
input=oldinput
input$data$agg_indices[2:35,7]=birdindex$estimate
input$data$agg_index_sigma[2:35,7]=birdindex$se

mod.dir="mm211_nonaa_noecov_linebird"
write.dir <- paste("C:/Herring/2025 Assessment RT/Assessments/WHAM",mod.dir,sep="/")
dir.create(write.dir)
setwd(write.dir)
mod <- fit_wham(input, do.proj=FALSE,do.osa = TRUE,do.retro = TRUE,do.check = T)
saveRDS(mod,file=file.path(write.dir,paste0(mod.dir,".rds")))
plot_wham_output(mod,out.type="png")

#This run has AR1 catchability informed by dummy covariate
oldinput=mods$mm196_birdecov$input
input=oldinput
input$data$agg_indices[2:35,7]=birdindex$estimate
input$data$agg_index_sigma[2:35,7]=birdindex$se

mod.dir="mm211_nonaa_birdecov"
write.dir <- paste("C:/Herring/2025 Assessment RT/Assessments/WHAM",mod.dir,sep="/")
dir.create(write.dir)
setwd(write.dir)
mod <- fit_wham(input, do.proj=FALSE,do.osa = TRUE,do.retro = TRUE,do.check = T)
saveRDS(mod,file=file.path(write.dir,paste0(mod.dir,".rds")))
plot_wham_output(mod,out.type="png")


##############
#read in new fits
runs=c("mm192_nonaa_nonecov","mm211_nonaa_noecov_linebird","mm211_nonaa_birdecov") #mm211 is mm192 but with linear seabird index included and NAA off so model can't tradeoff fits to data.
mod.list <- paste0(here::here("data/birdrunsforsh",runs),".rds")
mods <- lapply(mod.list, readRDS)
names(mods)=runs

##############
#plot and compare base, linear, and ar1 (AIC etc. should be comparable because ecov=None in base and linear)
mods.p <-mods  #c(list(mods$mm196_birdecov_rw,mod))
names(mods.p)= c("base","linear","ar1")
#names(mods.p) <- c("ASAP",df.mods$Model)
res <- compare_wham_models(mods.p, fdir="C:\\Herring\\2025 Assessment RT\\Survey\\seabird", 
                           do.table=T,calc.rho=F, do.retro = T,plot.opts = list(which=c(1:2),years=1987:2023))

#####
#function to extract results needed for plotting from a fit
getres=function(mod=NULL){
  if(ncol(mod$rep$log_index_resid)==7){
  res=data.frame("Year"=mod$input$years_full,"LogRecruits"=log(mod$rep$NAA[1,1,,1]),"LogResidual"=mod$rep$log_index_resid[,7],"Catchability"=mod$rep$q[,7])
  } else {
    res=data.frame("Year"=mod$input$years_full,"LogRecruits"=log(mod$rep$NAA[1,1,,1]),"LogResidual"=NA,"Catchability"=NA)  
  }
}

base=getres(mods$mm192_nonaa_nonecov)
base$Model="Base"
linear=getres(mods$mm211_nonaa_noecov_linebird)
linear$Model="Linear"
ar1=getres(mods$mm211_nonaa_birdecov)
ar1$Model="AR1"

df=rbind(base,linear,ar1)
df$Model=factor(df$Model,levels=c("Base","Linear","AR1"))

#recruitment time series
recplot=ggplot(df, aes(x=Year, y=exp(LogRecruits))) +
  geom_line(aes(linetype=Model),linewidth=1)+
  ylab("Age-1 Recruitment") +
  xlab("Year")+

  theme_bw()
  
recplot
ggsave(recplot,filename="recplot.jpeg",path="C:\\Herring\\2025 Assessment RT\\Survey\\seabird")

#bird index residuals versus recruitment
residplot=ggplot(df[df$Model != "Base",])+
  geom_point(aes(x=LogRecruits,y=LogResidual),size=1)+
  geom_smooth(aes(x=LogRecruits,y=LogResidual),method=lm,se=FALSE)+
  facet_wrap(~Model)+
  xlab("Log Recruitment")+
  ylab("Log Residual")
  ggsave(residplot,filename="residplot.jpeg",path="C:\\Herring\\2025 Assessment RT\\Survey\\seabird")

#q versus recruits from ecov fit
q_rec=ggplot(df[df$Model=="AR1",])+
  geom_point(aes(x=exp(LogRecruits),y=Catchability),size=1)+
  geom_smooth(aes(x=exp(LogRecruits),y=Catchability),method=loess,se=FALSE)+
  xlab("Age-1 Recruitment")+
  ylab("Catchability")
  ggsave(q_rec,filename="q_rec.jpeg",path="C:\\Herring\\2025 Assessment RT\\Survey\\seabird")
  

