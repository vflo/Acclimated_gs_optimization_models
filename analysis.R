# Improving stomatal optimization models for accurate prediction of photosynthesis under drought conditions.
# Victor Flo1,2, Jaideep Joshi3, Manon Sabot4,  David Sandoval1, Iain Colin Prentice1
# 1 Imperial College London, Department of Life Sciences, Silwood Park Campus, Ascot SL5 7PY, UK
# 2 Univ Autònoma de Barcelona, Cerdanyola del Vallès 08193, Spain
# 3Advancing Systems Analysis Program, International Institute for Applied Systems Analysis, 2361 Laxenburg, Austria
# 4ARC Centre of Excellence for Climate Extremes and Climate Change Research Centre, University of New South Wales, Sydney, NSW 2052, Australia
# Year 2023
#### RESULTS ####
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggpubr)
library(scales)
library(ggalt)
library(grid)
# devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)
library(rphydro)
library(ggplot2)
library(ggpattern)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(DEoptim)
library(MASS)
library(tidyverse)
source("stomatal_optimization_functions.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")
source("themes_and_results_functions.R")


# LOAD DATA --------------------------------------------------------------------

#### TRAITS DATA ####
traits <- read.csv(file="DATA/imputation_df.csv") %>% rename(species = 'Species', genus = 'Genus',Species = "Binomial")
traits[which(traits$Species == "Olea europaea var. meski"), "Species"] <- "Olea europaea var. Meski"
traits[which(traits$Species == "Olea europaea var. chemlali"), "Species"] <- "Olea europaea var. Chemlali"
traits[which(traits$Species == "Beta maritima subsp.marcosii"), "Species"] <- "Beta maritima subsp. marcosii"
traits[which(traits$Species == "Beta maritima subsp. maritima"), "Species"] <- "Beta maritima subsp. maritima"


#### TABLE LIFE FORM ####
lifeform <- read.csv(file="DATA/table_lifeform.csv")

#### SIMULATIONS ALPHA = 0.1 ####
# load(file = "DATA/simulations_kmax.RData")
# df_a_fix <- df %>% 
#   mutate(chi = Ciest/ca,
#          acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")) %>% 
#   filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera")))### Broussonetia papyrifera did not converged in SOX !!!!!!
# 
# 
# df_a_fix_n <- df_a_fix %>% 
#   group_by(scheme,acclimation,Species,source) %>% 
#   filter(!is.na(A)) %>% 
#   summarise(n_dist = n())
# 
# 
# df_a_fix %>% 
#   select(Species, scheme, acclimation, K.scale, gamma) %>% 
#   group_by(Species, scheme, acclimation) %>% 
#   mutate(FLAG = case_when(any(K.scale>49,gamma>49)~1,
#                           TRUE~2),
#          FLAG2 = case_when(any(K.scale>5.9&K.scale<6,gamma>5.9&gamma<6)~1,
#                            TRUE~2)
#   ) %>% filter(FLAG == 1)

## PARAMETERS DATA
# path_par <- "DATA/parameters_kmax/"
# par_data <- list.files(path_par) %>% 
#   purrr::map_df(function(x){
#     readr::read_csv(paste0(path_par,x))
#   })
# 
# df_a_fix_param <- par_data %>% 
#   mutate(acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")
#   )

## CALCULATE WW JMAX AND VCMAX
# vcmax_jmax_ww <- df_a_fix %>%
#   group_by(scheme,source,Species) %>%
#   do(get_vcmax_jmax_ww(.))


#### SIMULATIONS KMAXWW ALPHA ####
load(file = "DATA/simulations_kmaxww_alpha.RData")
df_kmaxww_a <- df  %>% 
  mutate(chi = Ciest/ca_pa,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated")),
         scheme = as_factor(scheme),
         scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2"),
         calibration_type = as_factor(calibration_type),
         calibration_type = fct_recode(calibration_type, 
                                       "Calibrated \u03B1" = "alpha",
                                       "Average \u03B1"="alpha_fix",
                                       "Not acclimated"="no_alpha"))#%>% 
  # filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera")))%>% ### Broussonetia papyrifera did not converged in SOX !!!!!!
df_kmaxww_a_n <- df_kmaxww_a %>%
    group_by(scheme,acclimation,calibration_type,Species,source) %>%
    filter(!is.na(A)) %>%
    summarise(n_dist = n())

df_kmaxww_a <- df_kmaxww_a %>%   
  left_join(df_kmaxww_a_n)






#### SIMULATIONS KMAXWW ALPHA P50ox ####
load(file = "DATA/simulations_kmaxww_alpha_1_3.RData")
df_kmaxww_a_1_3 <- df  %>% 
  mutate(chi = Ciest/ca_pa,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated")),
         scheme = as_factor(scheme),
         scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2"),
         calibration_type = as_factor(calibration_type),
         calibration_type = fct_recode(calibration_type, 
                                       "Calibrated \u03B1" = "alpha",
                                       "Average \u03B1"="alpha_fix",
                                       "Not acclimated"="no_alpha"))#%>% 
# filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera")))%>% ### Broussonetia papyrifera did not converged in SOX !!!!!!
df_kmaxww_a_n_1_3 <- df_kmaxww_a_1_3 %>%
  group_by(scheme,acclimation,calibration_type,Species,source) %>%
  filter(!is.na(A)) %>%
  summarise(n_dist = n())

df_kmaxww_a_1_3 <- df_kmaxww_a_1_3 %>%   
  left_join(df_kmaxww_a_n_1_3)
# df_kmaxww_a %>%
#   dplyr::select(Species, scheme, acclimation,
#                 calibration_type, K.scale, gamma,alpha) %>%
#   group_by(Species, scheme, acclimation) %>%
#   mutate(FLAG = case_when(any(K.scale>14,gamma>14)~1,
#                           TRUE~2),
#          FLAG2 = case_when(any(K.scale>5.9&K.scale<6,gamma>5.9&gamma<6)~1,
#                            TRUE~2)
#   ) %>% filter(FLAG == 1) %>% View()

## PARAMETERS DATA
# path_par <- "DATA/parameters_kmaxww_alpha/"
# par_data <- list.files(path_par) %>% 
#   purrr::map_df(function(x){
#     readr::read_csv(paste0(path_par,x))
#   })
# 
# df_kmaxww_a_param <- par_data %>% 
#   mutate(acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")
#   )


# #### KMAX NO ALPHA DATASET ####
# path_par_kmax_no_alpha <- "DATA/parameters_kmax_no_alpha/"
# par_data_kmax_no_alpha <- list.files(path_par_kmax_no_alpha) %>% 
#   purrr::map_df(function(x){
#     print(x)
#     readr::read_csv(paste0(path_par_kmax_no_alpha,x))
#   })
# 
# df_param_kmax_no_alpha <- par_data_kmax_no_alpha %>% 
#   mutate(acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")
#   )%>%
#   filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera")))### Broussonetia papyrifera did not converged in SOX !!!!!!%>%
# 
# 
# #### KMAX ALPHA DATASET ####
# path_par_kmax_alpha <- "DATA/parameters_kmax_alpha/"
# par_data_kmax_alpha <- list.files(path_par_kmax_alpha) %>% 
#   purrr::map_df(function(x){
#     readr::read_csv(paste0(path_par_kmax_alpha,x))
#   })
# 
# df_param_kmax_alpha <- par_data_kmax_alpha %>% 
#   mutate(stomatal_model = scheme,
#          acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")
#          )%>% 
#   filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera"))) %>% ### Broussonetia papyrifera did not converged in SOX !!!!!!%>%
#   left_join(vcmax_jmax_ww) %>% 
#   filter(K.scale > 1e-04)
# 
# 
# df_param_kmax_alpha <- df_param_kmax_alpha %>% 
#   left_join(df_a_fix %>%
#               dplyr::select(scheme,source,Species,Iabs_growth,D,T,ca) %>%
#               group_by(scheme,source,Species) %>%  
#               summarise_all(mean,na.rm=TRUE) 
#             ) %>% 
#   left_join(df_a_fix_n) %>%
#   left_join(lifeform) %>% 
#   mutate(var = 0,
#          psi0 = purrr::pmap(list(var, T,Iabs_growth,D,ca,K.scale,p50_opt,b_opt,alpha,gamma,stomatal_model),
#                             ~model_numerical(tc =..2, ppfd = ..3,
#                                           vpd = ..4*101325, co2 = ..5,
#                                           elv = 0, fapar = .99, kphio = 0.087,
#                                           psi_soil = ..1, rdark = 0.02,
#                                           par_plant=list(conductivity = ..6*1e-16,
#                                                          psi50 = ..7%>% unique(),
#                                                          b = ..8%>% unique()),
#                                           par_cost = list(
#                                             alpha = ..9,
#                                             gamma = ..10
#                                           ), stomatal_model = ..11 %>% unique())))%>%
#   unnest_wider(psi0,names_sep = ".") %>%
#   mutate(var = LWP_ww,
#          psi_ww = purrr::pmap(list(var, T,Iabs_growth,D,ca,K.scale,p50_opt,b_opt,alpha,gamma,stomatal_model), 
#                             ~model_numerical(tc =..2, ppfd = ..3, 
#                                              vpd = ..4*101325, co2 = ..5, 
#                                              elv = 0, fapar = .99, kphio = 0.087, 
#                                              psi_soil = ..1, rdark = 0.02, 
#                                              par_plant=list(conductivity = ..6*1e-16,
#                                                             psi50 = ..7%>% unique(),
#                                                             b = ..8%>% unique()), 
#                                              par_cost = list(
#                                                alpha = ..9,
#                                                gamma = ..10
#                                              ), stomatal_model = ..11 %>% unique())))%>% 
#   unnest_wider(psi_ww,names_sep = ".")
# 

# #### KMAX ALPHA PMODEL DATASET ####
# path_par_kmax_alpha_pmodel <- "DATA/parameters_kmax_alpha_pmodel/"
# par_data_kmax_alpha_pmodel <- list.files(path_par_kmax_alpha_pmodel) %>% 
#   purrr::map_df(function(x){
#     readr::read_csv(paste0(path_par_kmax_alpha_pmodel,x))
#   })
# 
# df_param_kmax_alpha_pmodel <- par_data_kmax_alpha_pmodel %>% 
#   mutate(stomatal_model = scheme,
#          acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")
#   )%>% 
#   filter(!(scheme %in% c("SOX")&Species %in%c("Broussonetia papyrifera"))) %>% ### Broussonetia papyrifera did not converged in SOX !!!!!!%>%
#   left_join(vcmax_jmax_ww) %>% 
#   filter(K.scale > 1e-04)%>% 
#   left_join(df_a_fix %>%
#               dplyr::select(scheme,source,Species,Iabs_growth,D,T,ca) %>%
#               group_by(scheme,source,Species) %>%  
#               summarise_all(mean,na.rm=TRUE) 
#   ) %>% 
#   left_join(df_a_fix_n) %>%
#   left_join(lifeform)

#### SIMULATIONS KMAX ALPHA ####
load(file = "DATA/simulations_kmax_alpha.RData")
df_kmax_alpha_n <- df %>%
  group_by(scheme,acclimation,Species,source) %>%
  filter(!is.na(A)) %>%
  summarise(acclimation = factor(acclimation, 
                                 levels = c('TRUE','FALSE'),
                                 labels = c("Acclimated", "Not acclimated")),
            scheme = as_factor(scheme),
            scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2"),
            n_dist = n()) %>% 
  unique()
df_kmax_alpha <- df %>% 
  mutate(chi = Ciest/ca_pa,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated")),
         scheme = as_factor(scheme),
         scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")) %>% 
  left_join(df_kmax_alpha_n) 

vcmaxww25_jmaxww25_pred <- df_kmax_alpha %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(LWP>=LWP_q_90) %>%
  summarise(vcmaxww25_pred =  mean(vcmax25_pred, na.rm = TRUE),
            jmaxww25_pred =  mean(jmax25_pred, na.rm = TRUE))

df_param_kmax_alpha <- df_kmax_alpha %>% 
  filter(acclimation == "Acclimated") %>% 
  filter(!is.na(Iabs_growth)) %>% 
  group_by(acclimation,scheme,Species, source) %>% 
  dplyr::select(n_dist,K.scale,b_opt,p50_opt,
                gamma,alpha, SLA..cm2.g.1., Huber.value, 
                KL..kg.m.1.MPa.1.s.1.,Iabs_growth,T,
                ca_pa, P50, b, Height.max..m., D,
                vcmaxww25,jmaxww25) %>% 
  summarise_all(mean) %>% 
  left_join(vcmaxww25_jmaxww25_pred)


df_param_kmax_no_alpha <- df_kmax_alpha %>% 
  filter(acclimation == "Not acclimated") %>% 
  filter(!is.na(Iabs_growth)) %>% 
  group_by(acclimation,scheme,Species, source) %>% 
  dplyr::select(n_dist,K.scale,b_opt,p50_opt,
                gamma,alpha, SLA..cm2.g.1., Huber.value, 
                KL..kg.m.1.MPa.1.s.1.,Iabs_growth,T,
                ca_pa, P50, b, Height.max..m., D,
                vcmaxww25,jmaxww25) %>% 
  summarise_all(mean) %>% 
  left_join(vcmaxww25_jmaxww25_pred)

# #### SIMULATIONS KMAXWW ALPHA ####
# load(file = "DATA/simulations_kmax_alpha.RData")
# df_kmax_alpha <- df %>% 
#   mutate(chi = Ciest/ca,
#          acclimation = factor(acclimation, 
#                               levels = c('TRUE','FALSE'),
#                               labels = c("Acclimated", "Not acclimated")),
#          scheme = fct_recode(scheme, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2"))
# 







#### Jmax vs alpha and Jmax vs Jmax ####

## jmax0
# jmax0_alpha_mod <- lmerTest::lmer(alpha~log(psi0.jmax)+  
#                                     (1|scheme) + (1|source) + (1|Species), 
#                                   data = df_param_kmax_alpha, weights = log(n_dist))
# anova(jmax0_alpha_mod)
# summary(jmax0_alpha_mod)
# r2 <- MuMIn::r.squaredGLMM(jmax0_alpha_mod)
# r2 <- round(r2[1],2)
# r2
# jmax0_alpha <- summary(jmax0_alpha_mod)$coefficients
# jmax0_sim <- seq(min(log(df_param_kmax_alpha$psi0.jmax),na.rm = TRUE),
#                 max(log(df_param_kmax_alpha$psi0.jmax),na.rm = TRUE),
#                 0.01)
# alpha_sim_jmax0 <- jmax0_alpha[1,1]+jmax0_alpha[2,1]*jmax0_sim
# alpha_sd_max_jmax0 <- jmax0_alpha[1,1] + jmax0_alpha[1,2] + (jmax0_alpha[2,1]+jmax0_alpha[2,2])*jmax0_sim
# alpha_sd_min_jmax0 <- jmax0_alpha[1,1] - jmax0_alpha[1,2] + (jmax0_alpha[2,1]-jmax0_alpha[2,2])*jmax0_sim
# resume_jmax0_alpha <- df_param_kmax_alpha %>% 
#   group_by(Species,source) %>% 
#   mutate(ln_psi0.jmax = log(psi0.jmax)) %>% 
#   dplyr::select(ln_psi0.jmax,alpha) %>% 
#   summarise_all(tibble::lst(mean,sd))
# (p1 <- ggplot()+
#   geom_point(data = df_param_kmax_alpha %>% mutate(ln_psi0.jmax = log(psi0.jmax)),
#                mapping=aes(x=ln_psi0.jmax,alpha,group=interaction(Species,source),size=log(n_dist)),
#                color="grey60", show.legend = FALSE,alpha=0.5)+
#   geom_point(resume_jmax0_alpha,mapping = aes(x=ln_psi0.jmax_mean, y = alpha_mean),color="grey20")+
#   geom_errorbar(resume_jmax0_alpha,
#                 mapping = aes(x=ln_psi0.jmax_mean, ymin = alpha_mean-alpha_sd,
#                               ymax = alpha_mean+alpha_sd),
#                 color="grey20")+
#   geom_errorbar(resume_jmax0_alpha,
#                 mapping = aes(y=alpha_mean, xmin = ln_psi0.jmax_mean-ln_psi0.jmax_sd,
#                               xmax = ln_psi0.jmax_mean+ln_psi0.jmax_sd),
#                 color="grey20")+
#   geom_line(aes(x = jmax0_sim, y = alpha_sim_jmax0),color = "grey20")+
#   geom_line(aes(x = jmax0_sim, y = alpha_sd_max_jmax0),color = "grey20",linetype=2)+
#   geom_line(aes(x = jmax0_sim, y = alpha_sd_min_jmax0),color = "grey20",linetype=2)+
#   annotate(geom = "text",x=5,y=0.155,hjust=0,vjust=0, parse = TRUE,
#            label=as.character(as.expression("***"~~italic(R)^2~"="~0.41)),size=4)+
#   # stat_summary(data = df_param_kmax_alpha,
#   #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
#   # geom="pointrange", color="blue")+
#   mytheme()+
#   # scale_radius(range=c(0.2,1))+
#   xlab(expression(J[max0]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
#   ylab(expression(alpha))
# )
# 
# 
# ##jmaxww_modelled
# df_params_kmax_alpha_filtered <- df_param_kmax_alpha %>% 
#   mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
#   filter(ln_psi_ww.jmax>0)
# jmax_ww_alpha_mod <- lmerTest::lmer(alpha~ln_psi_ww.jmax+I(ln_psi_ww.jmax^2)+ 
#                                     (1|scheme) + (1|source) + (1|Species), 
#                                   data = df_params_kmax_alpha_filtered,weights = log(n_dist))
# anova(jmax_ww_alpha_mod)
# summary(jmax_ww_alpha_mod)
# r2 <- MuMIn::r.squaredGLMM(jmax_ww_alpha_mod)
# r2 <- round(r2[1],2)
# r2
# jmax_ww_alpha <- summary(jmax_ww_alpha_mod)$coefficients
# jmax_ww_sim <- seq(min(log(df_params_kmax_alpha_filtered$psi_ww.jmax),na.rm = TRUE),
#                  max(log(df_params_kmax_alpha_filtered$psi_ww.jmax),na.rm = TRUE),
#                  0.01)
# alpha_sim_jmax_ww <- jmax_ww_alpha[1,1]+jmax_ww_alpha[2,1]*jmax_ww_sim+
#   jmax_ww_alpha[3,1]*jmax_ww_sim^2
# 
# alpha_sd_max_jmax_ww <- jmax_ww_alpha[1,1] + jmax_ww_alpha[1,2] + 
#   (jmax_ww_alpha[2,1]-jmax_ww_alpha[2,2])*jmax_ww_sim+ 
#   (jmax_ww_alpha[3,1]+jmax_ww_alpha[3,2])*jmax_ww_sim^2
# 
# alpha_sd_min_jmax_ww <- jmax_ww_alpha[1,1] - jmax_ww_alpha[1,2] + 
#   (jmax_ww_alpha[2,1]+jmax_ww_alpha[2,2])*jmax_ww_sim+ 
#   (jmax_ww_alpha[3,1]-jmax_ww_alpha[3,2])*jmax_ww_sim^2
# 
# resume_jmax_ww_alpha <- df_param_kmax_alpha %>% 
#   group_by(Species,source) %>% 
#   mutate(ln_psi_ww.jmax = log(psi_ww.jmax)) %>% 
#   filter(ln_psi_ww.jmax>0)%>% 
#   dplyr::select(ln_psi_ww.jmax,alpha) %>% 
#   summarise_all(tibble::lst(mean,sd))
# (p2 <- ggplot()+
#   geom_point(data = df_param_kmax_alpha %>% 
#                mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
#                filter(ln_psi_ww.jmax>0),
#              mapping=aes(x=ln_psi_ww.jmax,alpha,group=interaction(Species,source),size=log(n_dist)),
#              color="grey60", show.legend = FALSE,alpha=0.5)+
#   geom_point(resume_jmax_ww_alpha,mapping = aes(x=ln_psi_ww.jmax_mean, y = alpha_mean),color="grey20")+
#   # geom_smooth(data = df_param_kmax_alpha %>% 
#   #               mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
#   #               filter(ln_psi_ww.jmax>0),
#   #             mapping=aes(x=ln_psi_ww.jmax,alpha),
#   #             method = "lm",formula = y~poly(x,2))+
#   geom_errorbar(resume_jmax_ww_alpha,
#                 mapping = aes(x=ln_psi_ww.jmax_mean, ymin = alpha_mean-alpha_sd,
#                               ymax = alpha_mean+alpha_sd),
#                 color="grey20")+
#   geom_errorbar(resume_jmax_ww_alpha,
#                 mapping = aes(y=alpha_mean, xmin = ln_psi_ww.jmax_mean-ln_psi_ww.jmax_sd,
#                               xmax = ln_psi_ww.jmax_mean+ln_psi_ww.jmax_sd),
#                 color="grey20")+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sim_jmax_ww),color = "grey20")+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sd_max_jmax_ww),color = "grey20",linetype=2)+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sd_min_jmax_ww),color = "grey20",linetype=2)+
#   annotate(geom = "text",x=5,y=0.155,hjust=0,vjust=0, parse = TRUE,
#            label=as.character(as.expression("***"~~italic(R)^2~"="~0.33)),size=4)+
#   # stat_summary(data = df_param_kmax_alpha,
#   #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
#   # geom="pointrange", color="blue")+
#   mytheme3()+
#   # scale_radius(range=c(0.2,1))+
#   xlab(expression(J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
#   ylab(expression(alpha))
# )
# 
# 
# ##jmaxww observed
# # df_params_kmax_alpha_filtered <- df_param_kmax_alpha %>% 
# #   mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
# #   filter(ln_psi_ww.jmax>0)
# jmax_ww_alpha_mod <- lmerTest::lmer(alpha~log(jmax_ww)+ 
#                                       (1|scheme) + (1|source) + (1|Species), 
#                                     data = df_param_kmax_alpha, weights = log(n_dist))
# anova(jmax_ww_alpha_mod)
# summary(jmax_ww_alpha_mod)
# r2 <- MuMIn::r.squaredGLMM(jmax_ww_alpha_mod)
# r2 <- round(r2[1],2)
# r2
# jmax_ww_alpha <- summary(jmax_ww_alpha_mod)$coefficients
# jmax_ww_sim <- seq(min(log(df_param_kmax_alpha$jmax_ww),na.rm = TRUE),
#                    max(log(df_param_kmax_alpha$jmax_ww),na.rm = TRUE),
#                    0.01)
# alpha_sim_jmax_ww <- jmax_ww_alpha[1,1]+jmax_ww_alpha[2,1]*jmax_ww_sim
# alpha_sd_max_jmax_ww <- jmax_ww_alpha[1,1] + jmax_ww_alpha[1,2] + (jmax_ww_alpha[2,1]+jmax_ww_alpha[2,2])*jmax_ww_sim
# alpha_sd_min_jmax_ww <- jmax_ww_alpha[1,1] - jmax_ww_alpha[1,2] + (jmax_ww_alpha[2,1]-jmax_ww_alpha[2,2])*jmax_ww_sim
# resume_jmax_ww_alpha <- df_param_kmax_alpha %>% 
#   group_by(Species,source) %>% 
#   mutate(ln_jmax_ww = log(jmax_ww)) %>% 
#   filter(ln_jmax_ww>0)%>% 
#   dplyr::select(ln_jmax_ww,alpha) %>% 
#   summarise_all(tibble::lst(mean,sd))
# (p3 <- ggplot()+
#   geom_point(data = df_param_kmax_alpha %>% 
#                mutate(ln_jmax_ww = log(jmax_ww))%>% 
#                filter(ln_jmax_ww>0),
#              mapping=aes(x=ln_jmax_ww,alpha,group=interaction(Species,source),size=log(n_dist)),
#              color="grey60", show.legend = FALSE,alpha=0.5)+
#   geom_point(resume_jmax_ww_alpha,mapping = aes(x=ln_jmax_ww_mean, y = alpha_mean),color="grey20", size = 3)+
#   geom_errorbar(resume_jmax_ww_alpha,
#                 mapping = aes(x=ln_jmax_ww_mean, ymin = alpha_mean-alpha_sd,
#                               ymax = alpha_mean+alpha_sd),
#                 color="grey20")+
#   geom_errorbar(resume_jmax_ww_alpha,
#                 mapping = aes(y=alpha_mean, xmin = ln_jmax_ww_mean-ln_jmax_ww_sd,
#                               xmax = ln_jmax_ww_mean+ln_jmax_ww_sd),
#                 color="grey20")+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sim_jmax_ww),color = "grey20")+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sd_max_jmax_ww),color = "grey20",linetype=2)+
#   geom_line(aes(x = jmax_ww_sim, y = alpha_sd_min_jmax_ww),color = "grey20",linetype=2)+
#   annotate(geom = "text",x=4.5,y=0.155,hjust=0,vjust=0, parse = TRUE,
#            label=as.character(as.expression("***"~~italic(R)^2~"="~0.38)),size=5)+
#   # stat_summary(data = df_param_kmax_alpha,
#   #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
#   # geom="pointrange", color="blue")+
#   mytheme5()+
#   # scale_radius(range=c(0.2,1))+
#   xlab(expression("Observed"~J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
#   ylab(expression("Calibrated"~alpha))
# )
# 
# 
# (p4 <- ggplot(data = df_param_kmax_alpha,
#        mapping=aes(log(psi0.jmax),log(jmax_ww),color=scheme))+
#        geom_point()+
#        geom_abline(slope=1,intercept=0,color="grey20", linetype=2)+
#        geom_smooth(method="lm",se=FALSE)+
#        mytheme5()+
#        scale_colour_manual(breaks = col_df$scheme, 
#                           values = unique(as.character(col_df$col)))+
#        theme(legend.title = element_blank())+
#        xlab(expression("Modelled"~J[max0]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
#        ylab(expression("Observed"~J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))
#   )
#   
# (p5 <- ggplot(data = df_param_kmax_alpha,
#        mapping=aes(psi_ww.jmax,jmax_ww,color=scheme, weight=log(n_dist)))+
#        geom_point(aes(size =log(n_dist)),alpha = 0.5,show.legend = FALSE)+
#        geom_abline(slope=1,intercept=0,color="grey20", linetype=3,show.legend = FALSE)+
#        geom_smooth(method="lm",se = FALSE)+
#        mytheme5()+
#        scale_colour_manual(breaks = col_df$scheme, 
#                            values = unique(as.character(col_df$col)))+
#        theme(legend.title = element_blank())+
#        xlab(expression("Modelled"~J[maxWW]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#        ylab(expression("Observed"~J[maxWW]~"["*mu*"mol"~m^-2~s^-1*"]"))
#   )
# 
# jmax_jmax_plot <- ggarrange(p3,p5,
#                             align='h', labels=c('a', 'b'),
#                             common.legend = T,legend="bottom",
#                             ncol=1, nrow = 2)
# 
# ggsave("PLOTS/alpha_jmax_jmax_model.png",jmax_jmax_plot, width = 14, height = 28, units = "cm")
# 


#### Alpha multivariate analysis ####
# Customize upper panel
# upper.panel<-function(x, y){
#   points(x,y, pch=19, col=c("red", "green3", "blue")[iris$Species])
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   text(0.5, 0.9, txt)
# }
# pairs(df_param_kmax_alpha %>%
#         mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#                hv = Huber.value*1e4,
#                lk.scale =log(K.scale),
#                Height.max..m. = log(Height.max..m.),
#                l_b = log(b),
#                l_b_opt = log(b_opt),
#                l_p50_opt = log(-p50_opt),
#                D = D*101325/1000,
#                jmax_ww = log(jmaxww25),
#                Ta = T) %>%
#         ungroup() %>% 
#         dplyr::select(#alpha,SLA..cm2.g.1.,jmax_ww, KL,ca_pa,
#           Ta,
#           Iabs_growth,D #, l_b , lk.scale,
#           # P50,
#           # hv, #l_b_opt, l_p50_opt,
#           # Height.max..m.
#         ) %>% 
#         as.data.frame(), lower.panel = NULL, 
#       upper.panel = upper.panel)
library(psych)
pairs.panels(df_param_kmax_alpha %>%
               mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                      hv = Huber.value*1e4,
                      lk.scale =log(K.scale),
                      Height.max..m. = log(Height.max..m.),
                      l_b = log(b),
                      l_b_opt = log(b_opt),
                      l_p50_opt = log(-p50_opt),
                      D = D*101325/1000,
                      jmax_ww = log(jmaxww25),
                      Ta = T) %>%
               ungroup() %>% 
        dplyr::select(#alpha,SLA..cm2.g.1.,jmax_ww, KL,ca_pa,
                      Ta,
                      Iabs_growth,D , l_b ,# lk.scale,
                      P50,
                      # hv, #l_b_opt, l_p50_opt,
                      # Height.max..m.
                      ) %>% 
          rename(`Tair` = Ta,
                 Iabs = Iabs_growth,
                 VPD = D) %>% 
          as.data.frame(),
        ellipses = FALSE,
        hist.col = "grey40",
        smooth = FALSE
        )

# pairs.panels(df_param_kmax_alpha %>%
#                mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#                       # hv = Huber.value*1e4,
#                       lk.scale =log(K.scale),
#                       Height.max..m. = log(Height.max..m.),
#                       l_b = log(b),
#                       l_b_opt = log(b_opt),
#                       l_p50_opt = log(-p50_opt),
#                       D = D*101325/1000,
#                       jmax_ww = log(jmaxww25),
#                       # Ta = T
#                       ) %>%
#                ungroup() %>% 
#                filter(scheme=="PMAX") %>%
#                dplyr::select(#alpha,
#                       # SLA..cm2.g.1.,
#                       # jmax_ww,
#                       # KL,
#                       # ca_pa,
#                       # T,
#                       # Iabs_growth,
#                       # D,
#                       alpha,
#                       l_b ,
#                       P50,
#                       # lk.scale,P50,
#                       # hv, 
#                       l_b_opt, l_p50_opt,
#                       Height.max..m.
#                       ))
## FINAL alpha-traits ANALYSIS ##
#l_b and P50 are moderately-strongly correlated (we use only P50)
# df_param_kmax_alpha[is.na(df_param_kmax_alpha$lifeform),"lifeform"] <- "Evergreen needleleaf gymnosperm"
# df_param_kmax_alpha[is.na(df_param_kmax_alpha$lifeform_comp),"lifeform_comp"] <- "Evergreen"
df_param_kmax_alpha_lifeform <- df_param_kmax_alpha %>%left_join(lifeform)
alpha_mod <- lmerTest::lmer(alpha ~ SLA..cm2.g.1.+ KL+
                              Iabs_growth +
                              # D +
                              jmax_ww +# p50_opt + l_b_opt +
                              # Ta +
                              ca_pa+ 
                              P50+  
                              # l_b+
                              hv+
                              Height.max..m.+
                              # (1|source)+
                              (1|scheme),# + (1|lifeform),# + (1|Species),
                            data = df_param_kmax_alpha_lifeform %>% 
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                                     hv = Huber.value*1e4,
                                     lk.scale =log(K.scale),
                                     Height.max..m. =  log(Height.max..m.),
                                     l_b = log(b),
                                     l_b_opt = log(b_opt),
                                     l_p50_opt = log(-p50_opt),
                                     D = D*101325/1000,
                                     jmax_ww = log(jmaxww25),
                                     Ta = T), 
                            weights = log(n_dist)
                            )
summary(alpha_mod)
step(alpha_mod)
alpha_mod <- lmerTest::lmer(
  alpha ~ #SLA..cm2.g.1. + 
    # D + 
    Iabs_growth +
    jmax_ww + 
    # P50 + 
    # l_b + 
    Height.max..m. + (1 | scheme),
  data = df_param_kmax_alpha_lifeform %>%
    mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
           hv = log(Huber.value*1e4),
           lk.scale =log(K.scale),
           Height.max..m. = log(Height.max..m.),
           l_b = log(b),
           l_b_opt = log(b_opt),
           l_p50_opt = log(-p50_opt),
           D = D*101325/1000,
           jmax_ww = log(jmaxww25),
           Ta = T),
  weights = log(n_dist)
)
# alpha_mod <- lmerTest::lmer(
#   # alpha ~ SLA..cm2.g.1. + Iabs_growth + jmax_ww + hv + MATbest + P50 + (1 | source) + (1 | scheme),
#   # alpha ~ SLA..cm2.g.1. + Iabs_growth + jmax_ww + hv + Ta + P50 + (1 | source) + (1 | scheme),
#   alpha ~ Iabs_growth + jmax_ww + Ta + P50 + Height.max..m. + (1 | source) + (1 | scheme) + (1 | Species),
#   data = df_param_kmax_alpha %>% 
#     mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#            hv = log(Huber.value*1e4),
#            lk.scale =log(K.scale),
#            Height.max..m. = log( Height.max..m.),
#            l_b = log(b),
#            l_b_opt = log(b_opt),
#            l_p50_opt = log(-p50_opt),
#            D = D*101325/1000,
#            jmax_ww = log(jmax_ww),
#            Ta = T),
#   weights = log(n_dist)
# )
# anova(alpha_mod,alpha_mod2)
summary(alpha_mod)
car::vif(alpha_mod)
performance::check_collinearity(alpha_mod)
# confint(alpha_mod)
MuMIn::r.squaredGLMM(alpha_mod)


library(effects)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)

# 
# eff<-effect("SLA..cm2.g.1.", partial.residuals=T, alpha_mod)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, SLA..cm2.g.1. = eff$x$SLA..cm2.g.1.)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$SLA..cm2.g.1.)] + eff$residuals)
# 
# g_eff1 <- ggplot(x, aes(x = SLA..cm2.g.1., y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
#              col = "grey60", shape = 1, show.legend = FALSE) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression("SLA  [cm"^2~"g"^-1*"]"))

eff<-effect("Height.max..m.", partial.residuals=T, alpha_mod)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Height.max..m. = eff$x$Height.max..m.)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Height.max..m.)] + eff$residuals)

g_eff10 <- ggplot(x, aes(x = Height.max..m., y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression("Max height ["*log[e]*"(m)]"))

eff<-effect("Iabs_growth", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Iabs_growth = eff$x$Iabs_growth)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Iabs_growth)] + eff$residuals)

g_eff2 <- ggplot(x, aes(x = Iabs_growth, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression(I[abs]~"["*mu*"mol"~m^-2~s^-1*"]"))


# eff<-effect("hv", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, hv = eff$x$hv)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$hv)] + eff$residuals)
# 
# g_eff3 <- ggplot(x, aes(x = hv, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y, shape = df_param_kmax_alpha_lifeform$gym_ang), col = "grey60", size = 2, show.legend = FALSE) +
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression("Huber value [ln("*cm[sw]^2~m[leaf]^-2*")]"))
# 
# eff<-effect("D", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, D = eff$x$D)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$D)] + eff$residuals)
# 
# g_eff4 <- ggplot(x, aes(x = D, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
#              col = "grey60", shape = 1, show.legend = FALSE) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression("VPD [kPa]"))

# eff<-effect("Ta", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Ta = eff$x$Ta)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Ta)] + eff$residuals)
# 
# g_eff8 <- ggplot(x, aes(x = Ta, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(T[a]~"[ºC]"))
# 
# 
# eff<-effect("P50", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, P50 = eff$x$P50)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$P50)] + eff$residuals)
# 
# g_eff7 <- ggplot(x, aes(x = P50, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
#              col = "grey60", shape = 1, show.legend = FALSE) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(psi[50]~"[MPa]"))

# eff<-effect("l_b", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, l_b = eff$x$l_b)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$l_b)] + eff$residuals)
# 
# g_eff11 <- ggplot(x, aes(x = l_b, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
#              col = "grey60", shape = 1, show.legend = FALSE) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(log[e]*"(b)"))

# eff<-effect("KL", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, KL = eff$x$KL)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$KL)] + eff$residuals)
# 
# g_eff5 <- ggplot(x, aes(x = KL, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(K[L]~"[ln(mol "~m^{-1}~MPa^{-1}~s^{-1}~")]"))
# 
eff<-effect("jmax_ww", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, jmax_ww = eff$x$jmax_ww)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$jmax_ww)] + eff$residuals)

g_eff6 <- ggplot(x, aes(x = jmax_ww, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y), 
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression(J["maxWW,25"]~"["*log[e]*"("*mu*"mol"~m^-2~s^-1*")]"))

alpha_mod_plot <- ggarrange(g_eff2,g_eff6,
                            g_eff10,
          align='hv', labels=c('a', 'b','c'),
          ncol=3, nrow = 1)+
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

# ggsave("PLOTS/alpha_model.png", plot = alpha_mod_plot, width = 30, height = 20, units = "cm")
ggsave("PLOTS/alpha_model_5.pdf", plot = alpha_mod_plot, width = 32, height = 11, units = "cm")




#### Alpha - scheme ####
alpha_scheme <- lmerTest::lmer(alpha~ scheme + (1|Species)+ (1|source), 
                               data = df_param_kmax_alpha %>% 
                                mutate(scheme = factor(scheme, 
                                                       levels = c("SOX2",
                                                                  "SOX",
                                                                  "PMAX3",
                                                                  "PMAX2",
                                                                  "PMAX",
                                                                  "CGAIN2",
                                                                  "CGAIN",
                                                                  "PHYDRO"))), 
                              weights = log(n_dist))
alpha_scheme_mean <- lmerTest::lmer(alpha~ (1|scheme)+ (1|source)+(1|Species), 
                               data = df_param_kmax_alpha , weights = log(n_dist))
summary(alpha_scheme)
model_means <- emmeans(alpha_scheme, "scheme")
  
  model_means_cld <- cld(object = model_means,
                         adjust = "sidak",
                         Letters = letters,
                         alpha = 0.05)

  test(model_means, null = alpha_scheme_mean@beta )
  test(model_means, null = 0.1)

ggplot(data = model_means_cld) +
  geom_pointrange(mapping=aes(scheme, emmean, ymin=emmean - SE, ymax= emmean + SE)) +
  geom_abline(slope=0,intercept=fixef(alpha_scheme_mean),linetype=2)+
  geom_text(aes(y = emmean - SE -0.001, x = scheme, label = str_trim(.group)
            ))+
  ylab(expression(alpha)) +
  mytheme3()+
  xlab("")+
  coord_flip()

ggsave("PLOTS/alpha_scheme.png", width = 14, height = 10, units = "cm")







#### DARK RESPIRATION ####

df_param_kmax_alpha_lifeform <- df_param_kmax_alpha %>%left_join(lifeform)
alpha_mod <- lmerTest::lmer(alpha*jmaxww25 ~ SLA..cm2.g.1.+ KL+
                              Iabs_growth +
                              # D +
                              # jmax_ww +# p50_opt + l_b_opt +
                              # Ta +
                              ca_pa+ 
                              P50+  
                              # l_b+
                              hv+
                              Height.max..m.+
                              (1|source)+
                              (1|scheme),# + (1|lifeform),# + (1|Species),
                            data = df_param_kmax_alpha_lifeform %>% 
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                                     hv = Huber.value*1e4,
                                     lk.scale =log(K.scale),
                                     Height.max..m. =  log(Height.max..m.),
                                     l_b = log(b),
                                     l_b_opt = log(b_opt),
                                     l_p50_opt = log(-p50_opt),
                                     D = D*101325/1000,
                                     jmax_ww = log(jmaxww25),
                                     Ta = T), 
                            weights = log(n_dist)
)
summary(alpha_mod)
step(alpha_mod,reduce.random = FALSE)
alpha_mod <- lmerTest::lmer(
  # alpha * jmaxww25 ~ SLA..cm2.g.1. + KL + P50 + #l_b + hv + 
  #   Height.max..m. + (1 | source) + (1 | scheme),
  alpha * jmaxww25 ~ SLA..cm2.g.1. + KL + Iabs_growth + P50 + #l_b + 
    hv + Height.max..m. + (1 | source) + (1 | scheme),
  data = df_param_kmax_alpha_lifeform %>%
    mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
           hv = log(Huber.value*1e4),
           lk.scale =log(K.scale),
           Height.max..m. = log(Height.max..m.),
           l_b = log(b),
           l_b_opt = log(b_opt),
           l_p50_opt = log(-p50_opt),
           D = D*101325/1000,
           jmax_ww = log(jmaxww25),
           Ta = T),
  weights = log(n_dist)
)

# anova(alpha_mod,alpha_mod2)
summary(alpha_mod)
car::vif(alpha_mod)
performance::check_collinearity(alpha_mod)
# confint(alpha_mod)
MuMIn::r.squaredGLMM(alpha_mod)


library(effects)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)

eff<-effect("SLA..cm2.g.1.", partial.residuals=T, alpha_mod)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, SLA..cm2.g.1. = eff$x$SLA..cm2.g.1.)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$SLA..cm2.g.1.)] + eff$residuals)

g_eff1 <- ggplot(x, aes(x = SLA..cm2.g.1., y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression("SLA  [cm"^2~"g"^-1*"]"))

eff<-effect("Height.max..m.", partial.residuals=T, alpha_mod)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Height.max..m. = eff$x$Height.max..m.)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Height.max..m.)] + eff$residuals)

g_eff10 <- ggplot(x, aes(x = Height.max..m., y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression("Max height ["*log[e]*"(m)]"))


eff<-effect("P50", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, P50 = eff$x$P50)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$P50)] + eff$residuals)

g_eff7 <- ggplot(x, aes(x = P50, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression(psi[50]~"[MPa]"))


eff<-effect("KL", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, KL = eff$x$KL)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$KL)] + eff$residuals)

g_eff5 <- ggplot(x, aes(x = KL, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression(K[L]~"["*log[e]*"(mol "~m^{-1}~MPa^{-1}~s^{-1}~")]"))


eff<-effect("Iabs_growth", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Iabs_growth = eff$x$Iabs_growth)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Iabs_growth)] + eff$residuals)

g_eff2 <- ggplot(x, aes(x = Iabs_growth, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression(I[abs]~"["*mu*"mol"~m^-2~s^-1*"]"))


eff<-effect("hv", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, hv = eff$x$hv)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$hv)] + eff$residuals)

g_eff3 <- ggplot(x, aes(x = hv, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y, size = alpha_mod@frame$`(weights)`), 
             col = "grey60", shape = 1, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(R[photo]~"["*mu*"mol"~m^-2~s^-1*"]"))+
  xlab(expression("Huber value [ln("*cm[sw]^2~m[leaf]^-2*")]"))


alpha_mod_plot <- ggarrange(g_eff2,g_eff5,g_eff7,g_eff3,g_eff1,g_eff10,
                            align='hv', labels=c('a', 'b','c','d','e','f'),
                            ncol=3, nrow = 2)

# ggsave("PLOTS/alpha_model.png", plot = alpha_mod_plot, width = 30, height = 20, units = "cm")
ggsave("PLOTS/rdark_model.pdf", plot = alpha_mod_plot, width = 32, height = 20, units = "cm")



# 
# #### Alpha - scheme ####
# alpha_scheme <- lmerTest::lmer(alpha~ scheme + (1|Species)+ (1|source), 
#                                data = df_param_kmax_alpha %>% 
#                                  mutate(scheme = factor(scheme, 
#                                                         levels = c("SOX2",
#                                                                    "SOX",
#                                                                    "PMAX3",
#                                                                    "PMAX2",
#                                                                    "PMAX",
#                                                                    "CGAIN2",
#                                                                    "CGAIN",
#                                                                    "PHYDRO"))), 
#                                weights = log(n_dist))
# alpha_scheme_mean <- lmerTest::lmer(alpha~ (1|scheme)+ (1|source)+(1|Species), 
#                                     data = df_param_kmax_alpha , weights = log(n_dist))
# summary(alpha_scheme)
# model_means <- emmeans(alpha_scheme, "scheme")
# 
# model_means_cld <- cld(object = model_means,
#                        adjust = "sidak",
#                        Letters = letters,
#                        alpha = 0.05)
# 
# test(model_means, null = alpha_scheme_mean@beta )
# 
# ggplot(data = model_means_cld) +
#   geom_pointrange(mapping=aes(scheme, emmean, ymin=emmean - SE, ymax= emmean + SE)) +
#   geom_abline(slope=0,intercept=fixef(alpha_scheme_mean),linetype=2)+
#   geom_text(aes(y = emmean - SE -0.001, x = scheme, label = str_trim(.group)
#   ))+
#   ylab(expression(alpha)) +
#   mytheme3()+
#   xlab("")+
#   coord_flip()
# 
# ggsave("PLOTS/alpha_scheme.png", width = 14, height = 10, units = "cm")


# 
# library("factoextra")
# df_k_a_pca <- df_param_kmax_alpha%>% filter(scheme=="PROFITMAX") %>% 
#   mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#          hv = log(Huber.value),
#          lk.scale =log(K.scale),
#          l_b = log(b),
#          l_b_opt = log(b_opt),
#          D = D*101325/1000) %>% 
#   dplyr::select(alpha,T, ca,Height.max..m.,SLA..cm2.g.1.,
#                 Iabs_growth,
#                 D, l_b , lk.scale,P50,hv,p50_opt,l_b_opt
#                 )%>%
#   drop_na
# res.pca <- prcomp(df_k_a_pca, scale = TRUE)
# fviz_eig(res.pca)
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# fviz_pca_var(res.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )





#### A #####
df_a <- df_kmaxww_a %>%
# df_a <- df_kmaxww_a_1_3 %>%
  group_by(scheme,Species,source,calibration_type) %>% 
  filter(!is.na(A))  %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type))%>% 
  mutate(diff_a = a_pred - A) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>% 
  filter(beta> -10)

df_a_1_3 <- df_kmaxww_a_1_3 %>%
  group_by(scheme,Species,source,calibration_type) %>% 
  filter(!is.na(A))  %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type))%>% 
  mutate(diff_a = a_pred - A) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>% 
  filter(beta> -10)

df_a <- df_a %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_a_1_3 %>% filter(scheme == "PHYDRO"))

df_a %>% 
  group_by(scheme,calibration_type) %>% 
  dplyr::select(r,bias,rmse,beta) %>% 
  summarise(r_median = median(r),
            bias_median = median(bias),
            rmse_median = median(rmse),
            beta_median = median(beta),
            r_mean = mean(r),
            bias_mean = mean(bias),
            rmse_mean = mean(rmse),
            beta_mean = mean(beta)) ->foo

r_a <- lmerTest::lmer(r~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                      )

r_a_p <- pairs(emmeans(r_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value),
         calibration_type = as_factor(calibration_type)) %>% 
  dplyr::select(scheme, calibration_type,adj.p.value)
r_a <- emmeans(r_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>%
  mutate(calibration_type = as_factor(calibration_type)) %>% 
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,r,fill = calibration_type, color = calibration_type, group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" r Pearson's correlation"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1, 19))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


##BETA
beta_a <- lmerTest::lmer(beta~scheme*calibration_type + (1|Species), data = df_a%>% 
                           filter(beta<= 40), weights = log(n_dist)
                         )
test(emmeans::emmeans(beta_a, "calibration_type",by='scheme'),1)

beta_a_p <- pairs(emmeans(beta_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>%
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
beta_a <- emmeans(beta_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_a%>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated")) %>% 
  filter(beta<= 15)%>% 
  ggplot(aes(scheme,beta, fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = beta_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(-10,30)+
  ylab(expression(A[net]~italic(m)))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()+
  # scale_y_log10()+
  NULL


##RMSE
rmse_a <- lmerTest::lmer(rmse~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                         )
# emmeans(rmse_a, "calibration_type",by='scheme')
model_rmse <- emmeans::emmeans(rmse_a, "scheme",by='calibration_type')
cld(object = model_rmse,
    adjust = "sidak",
    Letters = letters,
    alpha = 0.05) %>% broom::tidy()
rmse_a_p <- pairs(emmeans(rmse_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
rmse_a <- emmeans(rmse_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

# rmse_a_accl <- lmerTest::lmer(rmse~scheme + (1|Species), data = df_a %>% filter(acclimation=="Acclimated"), weights = log(n_dist))
# summary(rmse_a_accl)
# rmse_a_p_accl <- emmeans(rmse_a_accl, "scheme")

p3 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,rmse,fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = rmse_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                         )
test(emmeans(bias_a, "calibration_type",by='scheme'))
bias_a_p <- pairs(emmeans(bias_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
bias_a <- emmeans(bias_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,bias,fill = calibration_type, color = calibration_type, 
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = bias_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" Bias"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


ggarrange(p1,
          p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'
                               ),
          common.legend = T,ncol=2, nrow = 2)

# ggarrange(p2,p3,
#           align='hv', labels=c("",""
#           ),
#           common.legend = T,ncol=2, nrow = 1)

ggsave("PLOTS/A_metrics_combined.png", width = 20, height = 14, units = "cm")


# 
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("PROFITMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX2") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("PROFITMAX2 Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# df_kmaxww_a_1_3 %>%
#   filter(scheme == "SOX") %>%
#   ggplot(aes(a_pred,A,color = calibration_type))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("SOX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "PHYDRO") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("PHYDRO Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# 
# df_a_fix %>% 
#   filter(scheme == "CMAX") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("CMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# 
# df_a_fix %>% 
#   filter(scheme == "CGAIN") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("CMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 




#### G ####
df_g <- df_kmaxww_a %>% 
  group_by(scheme,calibration_type,Species,source) %>% 
  filter(!is.na(gC)) %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type)) %>% 
  mutate(diff_g =  g_pred-gC) %>% 
  summarise(n_dist = n(),
            r = cor(gC, g_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_g,na.rm = TRUE)/mean(gC,na.rm = TRUE),
            rmse = Metrics::rmse(gC,g_pred),
            beta = lm(gC~g_pred)$coefficients[2])

df_g_1_3 <- df_kmaxww_a_1_3 %>% 
  group_by(scheme,calibration_type,Species,source) %>% 
  filter(!is.na(gC)) %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type)) %>% 
  mutate(diff_g =  g_pred-gC) %>% 
  summarise(n_dist = n(),
            r = cor(gC, g_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_g,na.rm = TRUE)/mean(gC,na.rm = TRUE),
            rmse = Metrics::rmse(gC,g_pred),
            beta = lm(gC~g_pred)$coefficients[2]) 
  
df_g <- df_g %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_g_1_3 %>% filter(scheme == "PHYDRO"))

df_g %>% 
  group_by(scheme,calibration_type) %>% 
  dplyr::select(r,bias,rmse,beta) %>% 
  summarise(r_median = median(r),
            bias_median = median(bias),
            rmse_median = median(rmse),
            beta_median = median(beta),
            r_mean = mean(r),
            bias_mean = mean(bias),
            rmse_mean = mean(rmse),
            beta_mean = mean(beta)) ->foo

r_a <- lmerTest::lmer(r~scheme*calibration_type + (1|source)+ (1|Species), data = df_g, weights = log(n_dist)
)

r_a_p <- pairs(emmeans(r_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value),
         calibration_type = as_factor(calibration_type)) %>% 
  dplyr::select(scheme, calibration_type,adj.p.value)
r_a <- emmeans(r_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_g %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,r,fill = calibration_type, color = calibration_type, group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(g[s]~" r Pearson's correlation"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1, 19))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


##BETA
beta_a <- lmerTest::lmer(beta~scheme*calibration_type + (1|Species), data = df_g%>% 
                           filter(beta<= 40), weights = log(n_dist)
)
test(emmeans::emmeans(beta_a, "calibration_type",by='scheme'),1)

beta_a_p <- pairs(emmeans(beta_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
beta_a <- emmeans(beta_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>%
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_g%>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>% 
  filter(beta<= 15)%>% 
  ggplot(aes(scheme,beta, fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = beta_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(-10,30)+
  ylab(expression(g[s]~italic(m)))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()+
  # scale_y_log10()+
  NULL


##RMSE
rmse_a <- lmerTest::lmer(rmse~scheme*calibration_type + (1|source)+ (1|Species), data = df_g, weights = log(n_dist)
)
# emmeans(rmse_a, "acclimation",by='scheme')
model_rmse <- emmeans::emmeans(rmse_a, "scheme",by='calibration_type')
# cld(object = model_rmse,
#     adjust = "sidak",
#     Letters = letters,
#     alpha = 0.05)
rmse_a_p <- pairs(emmeans(rmse_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
rmse_a <- emmeans(rmse_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

# rmse_a_accl <- lmerTest::lmer(rmse~scheme + (1|Species), data = df_g %>% filter(acclimation=="Acclimated"), weights = log(n_dist))
# summary(rmse_a_accl)
# rmse_a_p_accl <- emmeans(rmse_a_accl, "scheme")

p3 <- df_g %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,rmse,fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = rmse_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(g[s]~" RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*calibration_type + (1|source)+ (1|Species), data = df_g, weights = log(n_dist)
)
bias_a_p <- pairs(emmeans(bias_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
bias_a <- emmeans(bias_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_g %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,bias,fill = calibration_type, color = calibration_type, 
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5))+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = bias_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(g[s]~" Bias"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


ggarrange(p1,
          p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)

ggsave("PLOTS/gs_metrics.png", width = 20, height = 14, units = "cm")
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX2") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop(" PROFITMAX2 Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "SOX") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("SOX Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "PHYDRO") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("PHYDRO Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 
# 
# df_a_fix %>% 
#   filter(scheme == "CGAIN") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("CGAIN Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "CMAX") %>% 
#   ggplot(aes(g_pred,gC,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("CMAX Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 


#### A P50 1/3 #####
df_a <- df_kmaxww_a %>%
  # df_a <- df_kmaxww_a_1_3 %>%
  group_by(scheme,Species,source,calibration_type) %>% 
  filter(!is.na(A))  %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type))%>% 
  mutate(diff_a = a_pred - A) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>% 
  filter(beta> -10)

df_a_1_3 <- df_kmaxww_a_1_3 %>%
  group_by(scheme,Species,source,calibration_type) %>% 
  filter(!is.na(A))  %>% 
  mutate(scheme = factor(scheme, 
                         levels = c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO")),
         calibration_type = as_factor(calibration_type))%>% 
  mutate(diff_a = a_pred - A) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) %>% 
  filter(beta> -10)

df_a <- df_a %>% filter(scheme == "PHYDRO") %>% 
  bind_rows(df_a_1_3 %>% filter(scheme != "PHYDRO"))

df_a %>% 
  group_by(scheme,calibration_type) %>% 
  dplyr::select(r,bias,rmse,beta) %>% 
  summarise(r_median = median(r),
            bias_median = median(bias),
            rmse_median = median(rmse),
            beta_median = median(beta),
            r_mean = mean(r),
            bias_mean = mean(bias),
            rmse_mean = mean(rmse),
            beta_mean = mean(beta)) ->foo

r_a <- lmerTest::lmer(r~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
)

r_a_p <- pairs(emmeans(r_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value),
         calibration_type = as_factor(calibration_type)) %>% 
  dplyr::select(scheme, calibration_type,adj.p.value)
r_a <- emmeans(r_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>%
  mutate(calibration_type = as_factor(calibration_type)) %>% 
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,r,fill = calibration_type, color = calibration_type, group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5), alpha = 0.5)+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" r Pearson's correlation"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1, 19))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


##BETA
beta_a <- lmerTest::lmer(beta~scheme*calibration_type + (1|Species), data = df_a%>% 
                           filter(beta<= 40), weights = log(n_dist)
)
test(emmeans::emmeans(beta_a, "calibration_type",by='scheme'),1)

beta_a_p <- pairs(emmeans(beta_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>%
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
beta_a <- emmeans(beta_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_a%>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated")) %>% 
  filter(beta<= 15)%>% 
  ggplot(aes(scheme,beta, fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5), alpha = 0.5)+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = beta_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(-10,30)+
  ylab(expression(A[net]~italic(m)))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()+
  # scale_y_log10()+
  NULL


##RMSE
rmse_a <- lmerTest::lmer(rmse~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
)
# emmeans(rmse_a, "calibration_type",by='scheme')
model_rmse <- emmeans::emmeans(rmse_a, "scheme",by='calibration_type')
cld(object = model_rmse,
    adjust = "sidak",
    Letters = letters,
    alpha = 0.05) %>% broom::tidy()
rmse_a_p <- pairs(emmeans(rmse_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
rmse_a <- emmeans(rmse_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

# rmse_a_accl <- lmerTest::lmer(rmse~scheme + (1|Species), data = df_a %>% filter(acclimation=="Acclimated"), weights = log(n_dist))
# summary(rmse_a_accl)
# rmse_a_p_accl <- emmeans(rmse_a_accl, "scheme")

p3 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,rmse,fill = calibration_type, color = calibration_type,
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5), alpha = 0.5)+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = rmse_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" RMSE  ("*mu*"mol m"^{-2}~"s"^{-1}*")"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*calibration_type + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
)
test(emmeans(bias_a, "calibration_type",by='scheme'))
bias_a_p <- pairs(emmeans(bias_a, "calibration_type",by='scheme'), simple = "each")$`simple contrasts for calibration_type`%>% 
  broom::tidy() %>% 
  separate(contrast,into = c('calibration_type',"contrast"),sep = " - ") %>% 
  mutate(calibration_type = case_when(calibration_type == "Calibrated α" & contrast == "Average α"~ "Not acclimated",
                                      TRUE~calibration_type),
         adj.p.value = case_when(calibration_type == "Not acclimated"~1,
                                 TRUE~adj.p.value)) %>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  dplyr::select(scheme, calibration_type,adj.p.value)
bias_a <- emmeans(bias_a,~calibration_type*scheme) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  mutate(calibration_type = as_factor(calibration_type)) %>%
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_a %>% 
  mutate(calibration_type = fct_recode(calibration_type,
                                       "Calibrated \u03B1  "="Calibrated α",
                                       "Average \u03B1  " = "Average α", 
                                       "Not acclimated  " = "Not acclimated"))%>%
  ggplot(aes(scheme,bias,fill = calibration_type, color = calibration_type, 
             group = calibration_type))+
  geom_point(shape= 21,position=position_dodge(width = 0.5), alpha = 0.5)+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = bias_a%>% 
                    mutate(calibration_type = fct_recode(calibration_type,
                                                         "Calibrated \u03B1  "="Calibrated α",
                                                         "Average \u03B1  " = "Average α", 
                                                         "Not acclimated  " = "Not acclimated")), 
                  aes(scheme,estimate,color = calibration_type,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.5),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(A[net]~" Bias"))+
  scale_color_manual(values = c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_fill_manual(values =c("#FFC20A","#0C7BDC","#A2B56F"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  guides(colour = guide_legend(override.aes = list(size=4)))+
  coord_flip()


ggarrange(p1,
          p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'
          ),
          common.legend = T,ncol=2, nrow = 2)



ggsave("PLOTS/A_metrics_combined_1_3.png", width = 20, height = 14, units = "cm")




# STOMATAL NON-STOMATAL PARTITION ----------------------------------------------


partition_full <- df_kmaxww_a %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_kmaxww_a_1_3 %>% filter(scheme == "PHYDRO")) %>%
  filter(!is.na(A)) %>% 
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))),
         calibration_type = as_factor(calibration_type)) %>% 
  # group_by(scheme,species,source,calibration_type) %>% 
  # mutate(n_dist = n())%>% 
  # ungroup() %>% 
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>% 
  do(get_partition_a(.)) %>% 
  pivot_longer(2:5) %>% 
  mutate(Partition = factor(name,
                            levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
                            labels = c("Residuals  ","Species  ", "Non-stomatal  ","Stomatal  "))) %>%
  ggplot()+
  geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
  ylab(expression(R^2))+
  mytheme7()+
  ylim(0,1)+
  # scale_fill_manual(values = c("#d2d2d2","#ffae49","#44b7c2","#024b7a"))+
  scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
  theme(legend.title = element_blank())

partition_dry <- df_kmaxww_a %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_kmaxww_a_1_3 %>% filter(scheme == "PHYDRO")) %>% 
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))),
         calibration_type = as_factor(calibration_type)) %>% 
  group_by(scheme,calibration_type,Species,source) %>% 
  mutate(LWP_q50 = quantile(LWP, 0.5,na.rm = TRUE)) %>%
  filter(!is.na(A),LWP<=LWP_q50) %>% 
  mutate(n_dist = n())%>%
  ungroup() %>% 
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>% 
  do(get_partition_a(.)) %>%
  pivot_longer(2:5) %>% 
  mutate(Partition = factor(name,
                            levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
                            labels = c("Residuals  ","Species  ", "Non-stomatal  ","Stomatal  "))) %>% 
  ggplot()+
  geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
  ylab(expression(R^2~"dry conditions"))+
  mytheme7()+
  ylim(0,1)+
  scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
  theme(legend.title = element_blank())


partition_wet <- df_kmaxww_a %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_kmaxww_a_1_3 %>% filter(scheme == "PHYDRO")) %>% 
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                    "SOX",
                                    "PMAX3",
                                    "PMAX2",
                                    "PMAX",
                                    "CGAIN2",
                                    "CGAIN",
                                    "PHYDRO"))),
         calibration_type = as_factor(calibration_type)) %>% 
  group_by(scheme,calibration_type,Species,source) %>% 
  mutate(LWP_q50 = quantile(LWP, 0.5,na.rm = TRUE)) %>%
  filter(!is.na(A),LWP>=LWP_q50) %>% 
  mutate(n_dist = n())%>%
  ungroup() %>%
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>%
  do(get_partition_a(.)) %>%
  mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
                                 TRUE ~stomatal_r2),
         species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
                               TRUE ~species_r2)
  ) %>%
  pivot_longer(2:5) %>%
  mutate(pattern_1 = case_when(value<0~"stripe",
                             TRUE ~ "none"),
         value = case_when(value<0~ -1*value,
                           TRUE ~ value)) %>%
  mutate(Partition = factor(name,
                            levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
                            labels = c("Residuals  ","Species  ", "Non-stomatal  ","Stomatal  "))) %>%
  ggplot()+
  geom_col_pattern(aes(x = `Stomatal model`,y=value,fill=Partition,pattern = pattern_1),
                   pattern_color = NA,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   show.legend = FALSE)+
  scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
  ylab(expression(R^2~"wet conditions"))+
  mytheme7()+
  ylim(0,1)+
  scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
  theme(legend.title = element_blank())

ggarrange(partition_full, partition_dry,
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=2, nrow = 1)

ggsave("PLOTS/partition.png", width = 32, height = 14, units = "cm")

# 
# library(ggpattern)
# partition_dry_80 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q80 = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q80) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>% 
#   do(get_partition_a(.)) %>%
#   pivot_longer(2:5) %>% 
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>% 
#   ggplot()+
#   geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
#   ylab(expression(R^2~"80th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_dry_60 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q60 = quantile(LWP, 0.6,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q60) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>% 
#   do(get_partition_a(.)) %>%
#   pivot_longer(2:5) %>% 
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>% 
#   ggplot()+
#   geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
#   ylab(expression(R^2~"60th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_dry_70 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q70 = quantile(LWP, 0.7,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q70) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>% 
#   do(get_partition_a(.)) %>%
#   pivot_longer(2:5) %>% 
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>% 
#   ggplot()+
#   geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
#   ylab(expression(R^2~"70th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_dry_40 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q40 = quantile(LWP, 0.4,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q40) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>% 
#   do(get_partition_a(.)) %>%
#   pivot_longer(2:5) %>% 
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>% 
#   ggplot()+
#   geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
#   ylab(expression(R^2~"40th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_dry_30 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q30 = quantile(LWP, 0.3,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q30) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>%
#   do(get_partition_a(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                              TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Stomatal model`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"30th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_dry_20 <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q20 = quantile(LWP, 0.2,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP<=LWP_q20) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>%
#   do(get_partition_a(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                               TRUE ~species_r2)
#          ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>% 
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Stomatal model`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"20th percentile"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())

# library(ggpattern)
# partition_wet <- df_kmaxww_a %>% 
#   group_by(scheme,calibration_type,Species,source) %>% 
#   mutate(LWP_q50 = quantile(LWP, 0.5,na.rm = TRUE)) %>%
#   filter(!is.na(A),LWP>=LWP_q50) %>% 
#   mutate(n_dist = n())%>%
#   ungroup() %>% 
#   rename(`Stomatal model`=scheme) %>%
#   group_by(`Stomatal model`) %>% 
#   do(get_partition_a(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          Residuals= case_when(non_stomatal_r2<0~(Residuals+non_stomatal_r2),
#                               TRUE ~Residuals)) %>% 
#   pivot_longer(2:5) %>% 
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                              TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>% 
#   ggplot()+
#   geom_col_pattern(aes(x = `Stomatal model`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"wet conditions"))+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())





# ggarrange(partition_dry_80,partition_dry_70,partition_dry_60,partition_dry_40,partition_dry_30,partition_dry_20,
#           align='hv', labels=c('a', 'b','c','d','e','f','g'),
#           common.legend = T,ncol=3, nrow = 2)
# 
# ggsave("PLOTS/partition_multiple.png", width = 40, height = 24, units = "cm")






# 
# partition_PMAX <- df_kmaxww_a %>% 
#   filter(scheme == "PMAX") %>% 
#   group_by(calibration_type,Species,source) %>% 
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   # mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#   #                                TRUE ~stomatal_r2),
#   #        species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#   #                              TRUE ~species_r2)
#   # ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col(aes(x = `Drought quantile`,y=value,fill=Partition))+
#   # geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#   #                  pattern_color = NA,
#   #                  pattern_fill = "black",
#   #                  pattern_angle = 45,
#   #                  pattern_density = 0.5,
#   #                  pattern_spacing = 0.025,
#   #                  pattern_key_scale_factor = 1,
#   #                  show.legend = TRUE)+
#   # scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"PMAX"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# 
# partition_PMAX2 <- df_kmaxww_a %>% 
#   filter(scheme == "PMAX2") %>% 
#   group_by(calibration_type,Species,source) %>%
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   # mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#   #                                TRUE ~stomatal_r2),
#   #        species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#   #                              TRUE ~species_r2)
#   # ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"PMAX2"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_SOX <- df_kmaxww_a %>% 
#   filter(scheme == "SOX") %>% 
#   group_by(calibration_type,Species,source) %>% 
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"SOX"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_SOX2 <- df_kmaxww_a %>% 
#   filter(scheme == "SOX2") %>% 
#   group_by(calibration_type,Species,source) %>% 
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"SOX2"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_CGAIN <- df_kmaxww_a %>% 
#   filter(scheme == "CGAIN") %>% 
#   group_by(calibration_type,Species,source) %>% 
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   # mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#   #                                TRUE ~stomatal_r2),
#   #        species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#   #                              TRUE ~species_r2)
#   # ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"CGAIN"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_PHYDRO <- df_kmaxww_a %>% 
#   filter(scheme == "PHYDRO") %>% 
#   group_by(calibration_type,Species,source) %>%
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>%
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   # mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#   #                                TRUE ~stomatal_r2),
#   #        species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#   #                              TRUE ~species_r2)
#   # ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Drought quantile`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"PHYDRO"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# partition_CMAX <- df_kmaxww_a %>% 
#   filter(scheme == "CMAX") %>% 
#   group_by(calibration_type,Species,source) %>% 
#   mutate(`20th` = quantile(LWP, 0.2,na.rm = TRUE),
#          `30th` = quantile(LWP, 0.3,na.rm = TRUE),
#          `40th` = quantile(LWP, 0.4,na.rm = TRUE),
#          `50th` = quantile(LWP, 0.5,na.rm = TRUE),
#          `60th` = quantile(LWP, 0.6,na.rm = TRUE),
#          `70th` = quantile(LWP, 0.7,na.rm = TRUE),
#          `80th` = quantile(LWP, 0.8,na.rm = TRUE)) %>%
#   ungroup() %>% 
#   pivot_longer(71:77, names_to = 'Drought', values_to = 'value_quantile') %>% 
#   rename(`Drought quantile`=Drought) %>%
#   group_by(`Drought quantile`) %>%
#   do(get_partition_a_multiple(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Drought quantile` = factor(`Drought quantile`,
#                                      levels = c("80th","70th","60th","50th","40th","30th","20th"))) %>%
#   ggplot()+
#   geom_col_pattern(aes(x = `Drought quantile`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"CMAX"))+
#   xlab(psi[s]~" quantile")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# leg <- get_legend(partition_PMAX)
# 
# partition_PMAX <- partition_PMAX+ theme(legend.position = "none")
# 
# ggarrange(partition_PHYDRO,partition_CMAX,partition_CGAIN,
#           partition_PMAX,partition_PMAX2,partition_SOX,partition_SOX2,leg,
#           align='hv', labels=c('a', 'b','c','d','e','f',"g",""),
#           common.legend = FALSE,ncol=4, nrow = 2)
# 
# ggsave("PLOTS/partition_multiple_model.png", width = 40, height = 24, units = "cm")
# 


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# partition_PMAX <- df_kmaxww_a %>% 
#   filter(scheme == "PMAX", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Quantile bins` = factor(`Quantile bins`,
#                                   levels = c('q5','q4','q3', 'q2', 'q1'),
#                                   labels = c("80th-100th","60th-80th","40th-60th","20th-40th","0th-20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Quantile bins`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Quantile bins`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"PMAX"))+
#   xlab(psi[s]~" quantile bins")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# 
# partition_PMAX2 <- df_kmaxww_a %>% 
#   filter(scheme == "PMAX2", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          quite a lot of  scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_SOX <- df_kmaxww_a %>% 
#   filter(scheme == "SOX", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Quantile bins` = factor(`Quantile bins`,
#                                   levels = c('q5','q4','q3', 'q2', 'q1'),
#                                   labels = c("80th-100th","60th-80th","40th-60th","20th-40th","0th-20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Quantile bins`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Quantile bins`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"SOX"))+
#   xlab(psi[s]~" quantile bins")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_SOX2 <- df_kmaxww_a %>% 
#   filter(scheme == "SOX2", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Quantile bins` = factor(`Quantile bins`,
#                                   levels = c('q5','q4','q3', 'q2', 'q1'),
#                                   labels = c("80th-100th","60th-80th","40th-60th","20th-40th","0th-20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Quantile bins`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Quantile bins`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"SOX2"))+
#   xlab(psi[s]~" quantile bins")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# partition_CGAIN <- df_kmaxww_a %>% 
#   filter(scheme == "CGAIN", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Quantile bins` = factor(`Quantile bins`,
#                                   levels = c('q5','q4','q3', 'q2', 'q1'),
#                                   labels = c("80th-100th","60th-80th","40th-60th","20th-40th","0th-20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Quantile bins`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Quantile bins`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"CGAIN"))+
#   xlab(psi[s]~" quantile bins")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# 
# 
# partition_PHYDRO <- df_kmaxww_a %>% 
#   filter(scheme == "PHYDRO", !is.na(A)) %>% 
#   group_by(calibration_type,Species,source) %>% 
#   dplyr::select(A,a_pred,LWP) %>% 
#   # mutate(percentiles = binr::bins.quantiles(LWP,5,6)) %>% View()
#   mutate(`1` = quantile(LWP, 0,na.rm = TRUE),
#          `2` = quantile(LWP, 0.2,na.rm = TRUE),
#          `3` = quantile(LWP, 0.4,na.rm = TRUE),
#          `4` = quantile(LWP, 0.6,na.rm = TRUE),
#          `5` = quantile(LWP, 1,na.rm = TRUE),
#          `6` = quantile(LWP, 1,na.rm = TRUE),
#          q1 = "0th-20th",
#          q2 = "20th-40th",
#          q3 = "40th-60th",
#          q4 = "60th-100th",
#          q5 = "80th-100th"
#   ) %>%
#   ungroup() %>% 
#   pivot_longer(13:17, names_to = 'bins', values_to = 'Bin') %>%
#   rename(`Quantile bins`=bins) %>%
#   group_by(`Quantile bins`) %>%
#   do(get_partition_a_multiple_bin(.)) %>%
#   mutate(stomatal_r2 = case_when(non_stomatal_r2<0~(stomatal_r2+non_stomatal_r2),
#                                  TRUE ~stomatal_r2),
#          species_r2= case_when(non_stomatal_r2<0~(species_r2+non_stomatal_r2),
#                                TRUE ~species_r2)
#   ) %>%
#   pivot_longer(2:5) %>%
#   mutate(pattern_1 = case_when(value<0~"stripe",
#                                TRUE ~ "none"),
#          value = case_when(value<0~ -1*value,
#                            TRUE ~ value)) %>%
#   mutate(Partition = factor(name,
#                             levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
#                             labels = c("Residuals","Species", "Non-stomatal","Stomatal")),
#          `Quantile bins` = factor(`Quantile bins`,
#                                   levels = c('q5','q4','q3', 'q2', 'q1'),
#                                   labels = c("80th-100th","60th-80th","40th-60th","20th-40th","0th-20th"))) %>%
#   ggplot()+
#   # geom_col(aes(x = `Quantile bins`,y=value,fill=Partition))+
#   geom_col_pattern(aes(x = `Quantile bins`,y=value,fill=Partition,pattern = pattern_1),
#                    pattern_color = NA,
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.5,
#                    pattern_spacing = 0.025,
#                    pattern_key_scale_factor = 1,
#                    show.legend = FALSE)+
#   scale_pattern_manual(values = c(stripe = "stripe", none = "none"))+
#   ylab(expression(R^2~"PHYDRO"))+
#   xlab(psi[s]~" quantile bins")+
#   mytheme7()+
#   ylim(0,1)+
#   scale_fill_manual(values = c("#7DA5BB","#E7934C","#6448AC","#329587"))+
#   theme(legend.title = element_blank())
# 
# 
# ggarrange(partition_PHYDRO,partition_CGAIN,partition_PMAX,partition_PMAX2,partition_SOX,partition_SOX2,
#           align='hv', labels=c('a', 'b','c','d','e','f'),
#           common.legend = T,ncol=3, nrow = 2)
# 
# ggsave("PLOTS/partition_multiple_model_bin.png", width = 40, height = 24, units = "cm")
# 
# 
# 
# 





#### Species-scheme A r pearson's ####

library(gridExtra)
library(grid)
df_a <- df_kmaxww_a  %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_kmaxww_a_1_3 %>% filter(scheme == "PHYDRO")) %>%  
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))),
         calibration_type = as_factor(calibration_type)) %>% 
  group_by(scheme,calibration_type,Species,source) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2])

df_a %>% 
  ungroup() %>% 
  filter(calibration_type == "Calibrated α") %>% 
  dplyr::select(Species,scheme,r,source) %>% 
  mutate(Species = factor(Species)) %>% 
  # pivot_wider(values_from = 'r', names_from = 'scheme') %>% 
  # grid.table(theme = ttheme_minimal(
  #   core=list(bg_params = list(fill = blues9[1:4], col=NA),
  #             fg_params=list(fontface=3))))
ggpubr::ggballoonplot( 
  # aes(x = scheme, y = Species, fill = r)) +
  # geom_tile(color = "grey20") +
  # scale_fill_viridis_b(name = "r Pearson's") +
  # scale_x_discrete(name = "") +
  # scale_y_discrete(name = "")+
  # mytheme2()
x = "scheme", y = "Species",
size = "r", fill = "r",
size.range = c(1, 10)) +
  # scale_fill_viridis_c(name = "r Pearson's") +
  scale_fill_gradient(low = '#FFC20A', high = '#0C7BDC', name = "r Pearson's")+
  guides(size = "none")+scale_y_discrete(limits=rev)+
  theme(axis.text.y = element_text(face = "italic"))+
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))
    
ggsave("PLOTS/species_scheme_A.png", width = 16, height = 26, units = "cm")



#### Estimated vs Actual ####
df_kmax_alpha %>% 
  ggplot()+
  geom_point(aes(gC,A),shape = 1)+
  # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
  #             se = FALSE,method="lm",linetype = 1, size=0.5)+
  geom_smooth(aes(g_pred,a_pred,color=scheme),se = FALSE,method="lm", formula = y~poly(x,2))+
  facet_wrap(~acclimation)+
  mytheme2()+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
        legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("A ("*mu*"mol m"^-2~"s"^-1~")")))+
  xlab(expression(atop("g"[s]~"(mol m"^-2~"s"^-1*")")))+
  guides(colour = guide_legend(nrow = 1))+
  NULL

df_kmax_alpha %>% 
  ggplot()+
  geom_point(aes(gC,g_pred),shape = 1)+
  geom_abline(intercept=0,slope=1,color="grey20")+
  # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
  #             se = FALSE,method="lm",linetype = 1, size=0.5)+
  geom_smooth(aes(gC,g_pred,color=scheme),se = FALSE,method="lm", formula = y~poly(x,2))+
  facet_wrap(~acclimation)+
  mytheme2()+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("g"[s]~"estimated (mol m"^-2~"s"^-1~")")))+
  xlab(expression(atop("g"[s]~" actual (mol m"^-2~"s"^-1*")")))+
  guides(colour = guide_legend(nrow = 1))+
  NULL

df_kmax_alpha %>% 
  group_by(scheme,acclimation) %>% 
  filter(!is.na(A),!is.na(a_pred)) %>% 
  mutate(dens = get_density(A,a_pred,n = 50)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_point(aes(A,a_pred, fill=sqrt(dens)),
             size = 3, shape=21, alpha = 0.4, colour="transparent",show.legend = FALSE)+
  geom_abline(intercept=0,slope=1,color="grey20")+
  # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
  #             se = FALSE,method="lm",linetype = 1, size=0.5)+
  geom_smooth(aes(A,a_pred,color=scheme),se = FALSE,method="lm", formula = y~poly(x,2))+
  facet_grid(scheme~acclimation)+
  mytheme4()+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  scale_fill_gradient(low=c("#d2d2d2"),high = ("#ffae49"))+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("A estimated ("*mu*"mol m"^-2~"s"^-1~")")))+
  xlab(expression(atop("A actual ("*mu*"mol m"^-2~"s"^-1*")")))+
  guides(colour = guide_legend(nrow = 1))+
  NULL





#### Supplementary raw data ####
df_kmaxww_a %>%
  dplyr::filter(scheme == "PMAX", calibration_type == "Calibrated α") %>% 
  dplyr::group_by(Source,Species) %>% 
  dplyr::mutate(LWP_q90 = quantile(LWP, 0.8, na.rm = TRUE),
                is_90 = case_when(LWP >= LWP_q90 ~ ">80",
                                  TRUE ~ "<80")
           ) %>%
  ungroup() %>% 
  ggplot(aes(LWP,A,color = is_90)) +
  geom_point(show.legend = FALSE, alpha=0.8)+
  facet_wrap(~Species+Source,ncol=5, labeller = label_value, scales = "free_y")+
  # mytheme6()+
  theme_bw()+
  theme(strip.text = element_text(face = "italic",size=10),
        strip.background = element_rect(fill="transparent"))+
  ylab(expression(atop(A[net]~"("*mu*"mol m"^-2~"s"^-1*")")))+
  xlab(expression(italic(psi)[s]~"(MPa)"))+
  scale_colour_manual(values = c("grey20","#ffae49"))


ggsave("PLOTS/A_raw.png", width = 32, height = 40, units = "cm")




df_kmaxww_a %>%
  dplyr::filter(scheme == "PMAX", calibration_type == "Calibrated α") %>% 
  dplyr::group_by(Source,Species) %>% 
  dplyr::mutate(LWP_q90 = quantile(LWP, 0.8, na.rm = TRUE),
                is_90 = case_when(LWP >= LWP_q90 ~ ">80",
                                  TRUE ~ "<80")
  ) %>%
  ungroup() %>% 
  ggplot(aes(LWP,gC,color = is_90)) +
  geom_point(show.legend = FALSE, alpha=0.8)+
  facet_wrap(~Species+Source,ncol=5, labeller = label_value, scales = "free_y")+
  # mytheme6()+
  theme_bw()+
  theme(strip.text = element_text(face = "italic",size=10),
        strip.background = element_rect(fill="transparent"))+
  ylab(expression(g[s]~"("*mol[CO[2]]*m^-2~"s"^-1*")"))+
  xlab(expression(italic(psi)[s]~"(MPa)"))+
  scale_colour_manual(values = c("grey20","#ffae49"))

ggsave("PLOTS/gs_raw.png", width = 32, height = 40, units = "cm")




# 
# 
# #### Optimal parameters distributions ####
df_kmax_alpha %>%
  group_by(scheme,acclimation,Species,source) %>%
  # mutate(n_dist = n())%>%
  ungroup() %>%
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>%
  do(get_partition_a_accl(.)) %>%
  mutate(total_A_r2 = stomatal_r2+non_stomatal_r2)
# 
df_kmax_alpha %>%
  group_by(scheme,acclimation,Species,source) %>%
  # mutate(n_dist = n())%>%
  ungroup() %>%
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>%
  do(get_partition_g_accl(.)) %>%
  mutate(total_g_r2 = stomatal_r2+non_stomatal_r2)
# 
# 
# library(ggh4x)
# actual_dist <- df_param_kmax_alpha %>% 
#   bind_rows(df_param_kmax_no_alpha) %>%
#   mutate(Kmax1 = KL..kg.m.1.MPa.1.s.1.*55.5,
#          Kmax = case_when(Kmax1>0.5~NA_real_,
#                           TRUE~Kmax1)) %>%
#   group_by(scheme,acclimation) %>% 
#   dplyr::select(P50,Kmax,b) %>% 
#   pivot_longer(c(P50,Kmax,b)) %>% 
#   mutate(name = factor(name, levels = c("Kmax","P50","b","alpha","gamma")),
#          name = fct_recode(name,`italic(b)`="b",
#                            `gamma*"/"*omega1*"/"*beta[1]`='gamma',
#                            `K[max]~"(mol"~m^{-2}~s^{-1}~MPa^{-1}*")"`="Kmax",
#                            `psi[50]~"(MPa)"`='P50')) 
# 
# data_opt <- df_param_kmax_alpha %>%
#   bind_rows(df_param_kmax_no_alpha %>% left_join(df_param_kmax_alpha %>% dplyr::select(source,Species, scheme,`T`)))%>%
#   dplyr::select(scheme,acclimation,K.scale,gamma,alpha,p50_opt,b_opt,T) %>% 
#   rowwise() %>% 
#   mutate(Ta = T,
#          viscosity_water = calc_viscosity_h2o(Ta, 101325),
#          density_water = calc_density_h2o(Ta,101325),
#          Kmax = K.scale*1e-16/viscosity_water,
#          Kmax = Kmax * density_water * 55.5,
#          Kmax = Kmax * 1e6,
#          b_opt = case_when(b_opt>10~NA_real_,
#                           TRUE~b_opt)
#          )%>% 
#   group_by(scheme,acclimation) %>%
#   pivot_longer(c(Kmax,gamma,alpha,p50_opt,b_opt)) %>% 
#   mutate(name = factor(name, levels = c("Kmax","p50_opt","b_opt","alpha","gamma")),
#          name = fct_recode(name,`italic(b)`="b_opt",
#                            `gamma*"/"*omega1*"/"*beta[1]`='gamma',
#                            `K[max]~"(mol"~m^{-2}~s^{-1}~MPa^{-1}*")"`="Kmax",
#                            `psi[50]~"(MPa)"`='p50_opt')) 
# 
# ggplot()+
#   geom_density(data=actual_dist,aes(x=value), linewidth = 1,alpha=0.5)+
#   geom_density(data=data_opt,aes(x =value, color = acclimation),alpha=0.5,show.legend = TRUE, linewidth = 1)+
#   ggh4x::facet_grid2(c("scheme","name"), 
#               labeller = label_parsed, 
#               scales = "free", 
#               independent = "all")+
#   scale_color_manual(values = c("#FFC20A","#018571"))+
#   mytheme5()+
#   theme(legend.position="bottom",
#         axis.title.x = element_blank(),
#         legend.title = element_blank())
#   
#   
# ggsave("PLOTS/param_distri.png", width = 30, height = 30, units = "cm")
#   
#   

#### Optimal parameters boxplots ####
df_param_kmaxww <- df_kmaxww_a  %>% filter(scheme != "PHYDRO") %>% 
  bind_rows(df_kmaxww_a_1_3 %>% filter(scheme == "PHYDRO")) %>%  
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))),
         calibration_type = as_factor(calibration_type)) %>% 
    filter(calibration_type == "Calibrated α",
           !is.na(Iabs_growth)) %>% 
    rowwise() %>%
    mutate(Ta = mean(T),
           viscosity_water = calc_viscosity_h2o(Ta, 101325),
           density_water = calc_density_h2o(Ta,101325),
           Kmax = K.scale*1e-16/viscosity_water,
           Kmax = Kmax * density_water * 55.5,
           Kmax = Kmax * 1e6,
           KL = KL..kg.m.1.MPa.1.s.1.*55.5) %>%
    group_by(scheme,Species,source) %>%
    dplyr::select(n_dist,Kmax,
                  gamma,alpha, KL,
                  P50, b) %>% 
    summarise_all(mean)

df_param_kmax_no_alpha_comp <- df_param_kmax_no_alpha %>% 
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO")))) %>% 
  filter(!is.na(Iabs_growth)) %>% 
  rowwise() %>% 
  mutate(Ta = mean(T),
         viscosity_water = calc_viscosity_h2o(Ta, 101325),
         density_water = calc_density_h2o(Ta,101325),
         Kmax = K.scale*1e-16/viscosity_water,
         Kmax = Kmax * density_water * 55.5,
         Kmax = Kmax * 1e6,
         KL = KL..kg.m.1.MPa.1.s.1.*55.5) %>%  
  dplyr::select(n_dist,Kmax, p50_opt, b_opt,
                gamma,alpha, KL,
                P50, b,source,
                scheme,
                Species) %>% 
  rename(Kmax.cal3 = Kmax,
         p50_opt.cal3 = p50_opt,
         b_opt.cal3 = b_opt)

df_param_compar <- df_param_kmax_alpha %>% 
  mutate(scheme = factor(scheme, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO")))) %>% 
  filter(!is.na(Iabs_growth)) %>% 
  rowwise() %>% 
  mutate(Ta = mean(T),
         viscosity_water = calc_viscosity_h2o(Ta, 101325),
         density_water = calc_density_h2o(Ta,101325),
         Kmax = K.scale*1e-16/viscosity_water,
         Kmax = Kmax * density_water * 55.5,
         Kmax = Kmax * 1e6,
         KL = KL..kg.m.1.MPa.1.s.1.*55.5) %>%  
  dplyr::select(n_dist,Kmax, p50_opt, b_opt,
                gamma,alpha, KL,
                P50, b,source,scheme,Species) %>% 
  left_join(df_param_kmaxww, 
            by= c('scheme','Species','source', 'n_dist','P50','b','KL'),
            suffix = c(".cal1",".cal2")) %>% 
  left_join(df_param_kmax_no_alpha_comp, 
            by= c('scheme','Species','source', 'n_dist','P50','b','KL'))

  

# df_param_compar %>% 
#   ggplot(aes(log(Kmax.cal1),log(KL), color = scheme))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method="lm")
library(ggpubr)
library(rstatix)
#KMAX
df_foo <- df_param_compar %>%
  pivot_longer(cols = c(Kmax.cal1,Kmax.cal3,Kmax.cal2,KL), names_to = "cal_type") %>%
  dplyr::select(scheme, value, cal_type,Species) %>%
  rowwise() %>% 
  mutate(scheme = as.character(scheme),
         scheme = case_when(cal_type %in% "KL" ~ "Species",
                            TRUE ~ scheme),
         # scheme = as.factor(scheme),
         cal_type = as.factor(cal_type),
         cal_type = fct_recode(cal_type,
             "Species"="KL",
             "Calibration 1"='Kmax.cal1',
             "Calibration 1 not acclimated"='Kmax.cal3',
             "Calibration 2 and 3"='Kmax.cal2'
             ),
         log_value = log(value)) %>%
  dplyr::distinct()
  
stat.test <- df_foo %>%
  filter(scheme != "Species") %>% 
  group_by(scheme) %>%
  t_test(log_value ~ cal_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "scheme", dodge = 0.9)
stat.test <- stat.test %>% rbind(data.frame('scheme' = "Species",
                                 '.y.' = 'log_value',
                                 'group1' = "Species",
                                 'group2' = "Species",
                                 'n1' = 39,
                                 'n2' = 39,
                                 'statistic' = NA,
                                 'df' = NA,
                                 'p' = 1,
                                 'p.adj' = 1,
                                 'p.adj.signif' = "ns",
                                 'y.position' = 5.5,
                                 'groups' = "<chr>",
                                 'x' = 7,
                                 'xmin' = NA,
                                 'xmax' = NA))

p1 <- df_foo %>% 
  ggplot(aes(scheme,log(value), color = cal_type))+
  geom_violin(width = 1,position=position_dodge(.9))+
  geom_boxplot(width = 0.3,position=position_dodge(.9))+
  # geom_signif(comparisons=list(c("CGAIN", "PHYDRO","PMAX","PMAX2","SOX","SOX2")), 
  #             annotations="***",
  #             y_position = 0, tip_length = 0, vjust=0.4)+
  ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                              tip.length = 0, hide.ns = FALSE)+
  mytheme5()+
  ylab(expression(log[e]*"("*K[max]*")"~"["*log[e]*"(mol"~m^{-1}~s^{-1}~MPa^{-1}*")]"))+
  xlab("Stomatal model")+
  scale_color_manual(values = c("#FFC20A","#A2B56F","#0C7BDC","#AD9E77"))+
  theme(legend.position="top")+
  guides(color=guide_legend(title=""))+
  NULL

  

  
#P50
df_foo <- df_param_compar %>%
    mutate(P50 = P50) %>% 
    pivot_longer(cols = c(P50,p50_opt,p50_opt.cal3), names_to = "cal_type") %>%
    dplyr::select(scheme, value, cal_type,Species,source) %>%
    rowwise() %>% 
    mutate(scheme = as.character(scheme),
           scheme = case_when(cal_type %in% "P50" ~ "Species",
                              TRUE ~ scheme),
           # scheme = as.factor(scheme),
           cal_type = as.factor(cal_type),
           cal_type = fct_recode(cal_type,
                                 "Species"="P50",
                                 "Calibration 1"='p50_opt',
                                 "Calibration 1 not acclimated"='p50_opt.cal3'),
           value = -value
           ) %>%
    dplyr::distinct() %>% 
  ungroup()
  
stat.test <- df_foo %>%
  filter(scheme != "Species") %>% 
  group_by(scheme) %>%
  t_test(value ~ cal_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "scheme", dodge = 0.9)
stat.test <- stat.test %>% rbind(data.frame('scheme' = "Species",
                                            '.y.' = 'log_value',
                                            'group1' = "Species",
                                            'group2' = "Species",
                                            'n1' = 39,
                                            'n2' = 39,
                                            'statistic' = NA,
                                            'df' = NA,
                                            'p' = 1,
                                            'p.adj' = 1,
                                            'p.adj.signif' = "ns",
                                            'y.position' = 5.5,
                                            'groups' = "<chr>",
                                            'x' = 7,
                                            'xmin' = NA,
                                            'xmax' = NA))

  
p2 <- df_foo %>% 
    ggplot(aes(scheme, value, color = cal_type))+
    geom_violin(width = 1,position=position_dodge(.9))+
    geom_boxplot(width = 0.3,position=position_dodge(.9))+
    # ggpubr::stat_pvalue_manual(stat.test,  label = "{p.adj.signif}",
    #                             tip.length = 0, hide.ns = TRUE)+
    ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                              tip.length = 0, hide.ns = FALSE)+
    mytheme5()+
    ylab(expression(psi[50]~"[-"*MPa*"]"))+
    xlab("Stomatal model")+
    scale_color_manual(values = c("#AD9E77","#FFC20A","#A2B56F"))+
    theme(legend.position="top")+
    guides(color=guide_legend(title=""))+
  NULL

  
#b
df_foo <- df_param_compar %>%
  pivot_longer(cols = c(b,b_opt,b_opt.cal3), names_to = "cal_type") %>%
  dplyr::select(scheme, value, cal_type,Species,source) %>%
  rowwise() %>% 
  mutate(scheme = as.character(scheme),
         scheme = case_when(cal_type %in% "b" ~ "Species",
                            TRUE ~ scheme),
         # scheme = as.factor(scheme),
         cal_type = as.factor(cal_type),
         cal_type = fct_recode(cal_type,
                               "Species"="b",
                               "Calibration 1"='b_opt',
                               "Calibration 1 not acclimated"='b_opt.cal3')
  ) %>%
  dplyr::distinct() %>% 
  ungroup()
  
stat.test <- df_foo %>%
  filter(scheme != "Species") %>% 
  group_by(scheme) %>%
  t_test(value ~ cal_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "scheme", dodge = 0.9)
stat.test <- stat.test %>% rbind(data.frame('scheme' = "Species",
                                            '.y.' = 'log_value',
                                            'group1' = "Species",
                                            'group2' = "Species",
                                            'n1' = 39,
                                            'n2' = 39,
                                            'statistic' = NA,
                                            'df' = NA,
                                            'p' = 1,
                                            'p.adj' = 1,
                                            'p.adj.signif' = "ns",
                                            'y.position' = 5.5,
                                            'groups' = "<chr>",
                                            'x' = 7,
                                            'xmin' = NA,
                                            'xmax' = NA))
  
  
p3 <- df_foo %>% 
    ggplot(aes(scheme, value, color = cal_type))+
    geom_violin(width = 1,position=position_dodge(.9))+
    geom_boxplot(width = 0.3,position=position_dodge(.9))+
    ggpubr::stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", y.position = 0,
                              tip.length = 0, hide.ns = FALSE)+
    mytheme5()+
    ylab(expression("b"))+
    xlab("Stomatal model")+
    scale_color_manual(values = c("#AD9E77","#FFC20A","#A2B56F"))+
    theme(legend.position="top")+
    guides(color=guide_legend(title=""))+
    NULL
 
  
#alpha
df_foo <- df_param_compar %>%
    pivot_longer(cols = c(alpha.cal1,alpha.cal2), names_to = "cal_type") %>%
    dplyr::select(scheme, value, cal_type,Species) %>%
    rowwise() %>% 
    mutate(scheme = as.character(scheme),
           # scheme = as.factor(scheme),
           cal_type = as.factor(cal_type),
           cal_type = fct_recode(cal_type,
                                 "Calibration 1"='alpha.cal1',
                                 "Calibration 2 and 3"='alpha.cal2'),
           value = value) %>%
    dplyr::distinct()
  
  stat.test <- df_foo %>%
    group_by(scheme) %>%
    t_test(value ~ cal_type) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")%>%
    add_xy_position(x = "scheme", dodge = 0.9)

  
p4 <- df_foo %>% 
    ggplot(aes(scheme,value, color = cal_type))+
    geom_violin(width = 1,position=position_dodge(.9))+
    geom_boxplot(width = 0.3,position=position_dodge(.9))+
    # geom_signif(comparisons=list(c("CGAIN", "PHYDRO","PMAX","PMAX2","SOX","SOX2")), 
    #             annotations="***",
    #             y_position = 0, tip_length = 0, vjust=0.4)+
    ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                                tip.length = 0, hide.ns = FALSE)+
    mytheme5()+
    ylab(expression(alpha))+
    xlab("Stomatal model")+
    scale_color_manual(values = c("#FFC20A","#0C7BDC"))+
    theme(legend.position="top")+
    guides(color=guide_legend(title=""))+
    NULL

ggarrange(p1,p2,p3,p4,
          align='hv', labels=c('a', 'b','c','d'),
          common.legend = T,ncol=2, nrow = 2)
  

ggsave("PLOTS/par_comp.pdf", width = 36, height = 25, units = "cm")







#P50 only observed values
df_foo <- df_param_compar %>% 
  filter(Species %in%c('Broussonetia papyrifera',
                       'Platycarya longipes',
                       'Pteroceltis tatarinowii',
                       'Cedrus atlantica',
                       'Pseudotzuga menziesii',
                       'Eucalyptus populnea',
                       'Helianthus annuus',
                       'Olea europaea var. Chemlali',
                       'Olea europaea var. Meski',
                       'Picea abies',
                       'Pinus sylvestris',
                       'Betula pendula',
                       'Populus tremula',
                       'Quercus suber',
                       'Quercus ilex',
                       'Quercus petraea',
                       'Quercus pubescens'
  ))%>% 
  pivot_longer(cols = c(P50,p50_opt,p50_opt.cal3), names_to = "cal_type") %>%
  dplyr::select(scheme, value, cal_type,Species,source) %>%
  rowwise() %>% 
  mutate(scheme = as.character(scheme),
         scheme = case_when(cal_type %in% "P50" ~ "Species",
                            TRUE ~ scheme),
         # scheme = as.factor(scheme),
         cal_type = as.factor(cal_type),
         cal_type = fct_recode(cal_type,
                               "Species"="P50",
                               "Calibration 1"='p50_opt',
                               "Calibration 1 not acclimated"='p50_opt.cal3'),
         value = -value
  ) %>%
  dplyr::distinct() %>% 
  ungroup()

stat.test <- df_foo %>%
  filter(scheme != "Species") %>% 
  group_by(scheme) %>%
  t_test(value ~ cal_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "scheme", dodge = 0.9)
stat.test <- stat.test %>% rbind(data.frame('scheme' = "Species",
                                            '.y.' = 'log_value',
                                            'group1' = "Species",
                                            'group2' = "Species",
                                            'n1' = 39,
                                            'n2' = 39,
                                            'statistic' = NA,
                                            'df' = NA,
                                            'p' = 1,
                                            'p.adj' = 1,
                                            'p.adj.signif' = "ns",
                                            'y.position' = 5.5,
                                            'groups' = "<chr>",
                                            'x' = 7,
                                            'xmin' = NA,
                                            'xmax' = NA))


df_foo %>% 
  ggplot(aes(scheme, value, color = cal_type))+
  geom_violin(width = 1,position=position_dodge(.9))+
  geom_boxplot(width = 0.3,position=position_dodge(.9))+
  # ggpubr::stat_pvalue_manual(stat.test,  label = "{p.adj.signif}",
  #                             tip.length = 0, hide.ns = TRUE)+
  ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                              tip.length = 0, hide.ns = FALSE)+
  mytheme5()+
  ylab(expression(psi[50]~"[-"*MPa*"]"))+
  xlab("Stomatal model")+
  scale_color_manual(values = c("#AD9E77","#FFC20A","#A2B56F"))+
  theme(legend.position="top")+
  guides(color=guide_legend(title=""))+
  NULL


df_param_compar %>% 
  filter(Species %in%c('Broussonetia papyrifera',
                       'Platycarya longipes',
                       'Pteroceltis tatarinowii',
                       'Cedrus atlantica',
                       'Pseudotzuga menziesii',
                       'Eucalyptus populnea',
                       'Helianthus annuus',
                       'Olea europaea var. Chemlali',
                       'Olea europaea var. Meski',
                       'Picea abies',
                       'Pinus sylvestris',
                       'Betula pendula',
                       'Populus tremula',
                       'Quercus suber',
                       'Quercus ilex',
                       'Quercus petraea',
                       'Quercus pubescens'
  )) %>% 
  dplyr::select(scheme, Species,source,P50,p50_opt,p50_opt.cal3) %>% 
  # filter(scheme=="PMAX") %>% 
  mutate(p50_diff_accl = P50 - p50_opt,
         p50_diff_no_accl = P50 - p50_opt.cal3,
         p50_diff_no_accl_1_3 = P50/3 - p50_opt.cal3) %>% 
  pivot_longer(cols = c(p50_diff_accl, p50_diff_no_accl,p50_diff_no_accl_1_3)) %>% 
  mutate(name = fct_recode(name,
                           "Acclimated"='p50_diff_accl',
                           "Not acclimated"='p50_diff_no_accl',
                           "Not Acclimated P50/3"='p50_diff_no_accl_1_3')) %>% 
  ggplot(aes(Species, value, color = name))+
  geom_abline(slope=0,intercept=0,linetype=2)+
  # geom_violin(width = 1,position=position_dodge(.9))+
  stat_summary(,geom="pointrange",width = 0.3,position=position_dodge(.4))+
  stat_summary(aes(Species,P50),geom="point",fun = "mean", shape=3, color="black")+
  # # ggpubr::stat_pvalue_manual(stat.test,  label = "{p.adj.signif}",
  # #                             tip.length = 0, hide.ns = TRUE)+
  # ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
  #                             tip.length = 0, hide.ns = FALSE)+
  # mytheme5()+
  theme_bw()+
  ylab(expression(psi[50]*" - "*psi[`50opt`]~"[MPa] "))+
  xlab("Species")+
  scale_color_manual(values = c("#AD9E77","#FFC20A","#A2B56F"))+
  theme(legend.position="top")+
  guides(color=guide_legend(title=""))+
  coord_flip()+
  NULL


df_param_compar %>% 
  filter(Species %in%c('Broussonetia papyrifera',
                       'Platycarya longipes',
                       'Pteroceltis tatarinowii',
                       'Cedrus atlantica',
                       'Pseudotzuga menziesii',
                       'Eucalyptus populnea',
                       'Helianthus annuus',
                       'Olea europaea var. Chemlali',
                       'Olea europaea var. Meski',
                       'Picea abies',
                       'Pinus sylvestris',
                       'Betula pendula',
                       'Populus tremula',
                       'Quercus suber',
                       'Quercus ilex',
                       'Quercus petraea',
                       'Quercus pubescens'
  )) %>% 
  dplyr::select(scheme, Species,source,b,b_opt,b_opt.cal3) %>% 
  # filter(scheme=="PMAX") %>% 
  mutate(b_diff_accl = b - b_opt,
         b_diff_no_accl = b - b_opt.cal3,
         b_diff_no_accl_1_3 = 1 - b_opt.cal3
         ) %>% 
  pivot_longer(cols = c(b_diff_accl, b_diff_no_accl,b_diff_no_accl_1_3)) %>% 
  mutate(name = fct_recode(name,
                           "Acclimated"='b_diff_accl',
                           "Not acclimated b=1" = "b_diff_no_accl_1_3",
                           "Not acclimated"='b_diff_no_accl')) %>% 
  ggplot(aes(Species, value, color = name))+
  geom_abline(slope=0,intercept=0,linetype=2)+
  # geom_violin(width = 1,position=position_dodge(.9))+
  stat_summary(,geom="pointrange",width = 0.3,position=position_dodge(.4))+
  stat_summary(aes(Species,b),geom="point",fun = "mean", shape=3, color="black")+
  # # ggpubr::stat_pvalue_manual(stat.test,  label = "{p.adj.signif}",
  # #                             tip.length = 0, hide.ns = TRUE)+
  # ggpubr:: stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
  #                             tip.length = 0, hide.ns = FALSE)+
  # mytheme5()+
  theme_bw()+
  ylab(expression("b - "*b[opt]))+
  xlab("Species")+
  scale_color_manual(values = c("#AD9E77","#FFC20A","#A2B56F"))+
  theme(legend.position="top")+
  guides(color=guide_legend(title=""))+
  coord_flip()+
  NULL
















#### Example plots ####
  
# par_plant1 = list(
#     conductivity = 1e-16,
#     psi50 = -1,                   #
#     b=2
#   )
#   
# par_plant2 = list(
#     conductivity = 1e-16,
#     psi50 = -1.5,                   #
#     b=2
#   )
  
  

    
calc_examples <- function(scheme,par_cost,par_plant,par_plant_1_3){
  df1 <- tibble(var=seq(-3,0,length.out=20)) %>% 
    mutate(out_hydraulics_1 = purrr::map(var, ~model_numerical(tc = 25, ppfd = 600,
                                                               vpd = 1000, co2 = 400, elv = 0,
                                                               fapar = 0.99, kphio = 0.087, 
                                                               psi_soil = ., par_plant = par_plant, rdark = 0.015,
                                                               par_cost = list(alpha=par_cost[[1]]+0.02, gamma = par_cost[[2]]),
                                                               stomatal_model = scheme
                                                               ))) %>% 
    unnest_wider(out_hydraulics_1)%>% 
    bind_cols(stomatal_model = scheme, alpha_mod=1, p50_mod=1)
  df2 <- tibble(var=seq(-3,0,length.out=20)) %>% 
    mutate(out_hydraulics_1 = purrr::map(var, ~model_numerical(tc = 25, ppfd = 600,
                                                               vpd = 1000, co2 = 400, elv = 0,
                                                               fapar = 0.99, kphio = 0.087, 
                                                               psi_soil = ., par_plant = par_plant, rdark = 0.015,
                                                               par_cost = list(alpha=par_cost[[1]]-0.02, gamma = par_cost[[2]]),
                                                               stomatal_model = scheme
    ))) %>% 
    unnest_wider(out_hydraulics_1)%>% 
    bind_cols(stomatal_model = scheme, alpha_mod=2, p50_mod=1)
  df3 <- tibble(var=seq(-3,0,length.out=20)) %>% 
    mutate(out_hydraulics_1 = purrr::map(var, ~model_numerical(tc = 25, ppfd = 600,
                                                               vpd = 1000, co2 = 400, elv = 0,
                                                               fapar = 0.99, kphio = 0.087, 
                                                               psi_soil = ., par_plant = par_plant_1_3, rdark = 0.015,
                                                               par_cost = list(alpha=par_cost[[1]]+0.02, gamma = par_cost[[2]]),
                                                               stomatal_model = scheme
    ))) %>% 
    unnest_wider(out_hydraulics_1)%>% 
    bind_cols(stomatal_model = scheme, alpha_mod=1, p50_mod=2)
  df4 <- tibble(var=seq(-3,0,length.out=20)) %>% 
    mutate(out_hydraulics_1 = purrr::map(var, ~model_numerical(tc = 25, ppfd = 600,
                                                               vpd = 1000, co2 = 400, elv = 0,
                                                               fapar = 0.99, kphio = 0.087, 
                                                               psi_soil = ., par_plant = par_plant_1_3, rdark = 0.015,
                                                               par_cost = list(alpha=par_cost[[1]]-0.02, gamma = par_cost[[2]]),
                                                               stomatal_model = scheme
    ))) %>% 
    unnest_wider(out_hydraulics_1)%>% 
    bind_cols(stomatal_model = scheme, alpha_mod=2, p50_mod=2)

  bind_rows(df1,df2,df3,df4)
}
  
  
PMAX <- calc_examples("PROFITMAX", list(0.0945,NA),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                      list(conductivity =1e-16,psi50 =-1.58/3,b=1))
PMAX2 <- calc_examples("PROFITMAX2", list(0.0945,NA),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                       list(conductivity =1e-16,psi50 =-1.58/3,b=1))
PMAX3 <- calc_examples("PMAX3", list(0.0945,NA),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                       list(conductivity =1e-16,psi50 =-1.58/3,b=1))
SOX <- calc_examples("SOX", list(0.0945,NA),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                     list(conductivity =1e-16,psi50 =-1.58/3,b=1))
SOX2 <- calc_examples("SOX2", list(0.0945,2),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                      list(conductivity =1e-16,psi50 =-1.58/3,b=1))
PHYDRO <- calc_examples("PHYDRO", list(0.0945,1),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                        list(conductivity =1e-16,psi50 =-1.58/3,b=1))
# CMAX <- calc_examples("CMAX", list(0.1,1.2))
CGAIN <- calc_examples("CGAIN", list(0.0945,5),list(conductivity =1e-16,psi50 =-1.58,b=3.12),
                       list(conductivity =1e-16,psi50 =-1.58/3,b=1))
  
DF_example <- bind_rows(PMAX,PMAX2,PMAX3,SOX,SOX2,PHYDRO,CGAIN) %>% 
  mutate(scheme = factor(stomatal_model, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))))

DF_example <- DF_example %>% 
  dplyr::select(var, jmax25, vcmax25, gs, a, stomatal_model, alpha_mod, p50_mod) %>% 
  pivot_longer(cols = c('jmax25','vcmax25','gs','a'))%>% 
  mutate(alpha_mod = as.factor(alpha_mod),
         p50_mod = as.factor(p50_mod)) %>% 
  mutate(name = factor(name, levels = c('gs','a','jmax25','vcmax25')),
         name = fct_recode(name,
                           `g[s]~"(mol"~m^{-2}~s^{-1}*")"`="gs",
                           `"A ("*mu*"mol"~m^{-2}~s^{-1}~MPa^{-1}*")"`="a",
                           `J["max,25"]~"("*mu*"mol"~m^{-2}~s^{-1}~MPa^{-1}*")"`="jmax25",
                           `V["max,25"]~"("*mu*"mol"~m^{-2}~s^{-1}~MPa^{-1}*")"`="vcmax25"),
         stomatal_model = as.factor(stomatal_model),
         stomatal_model = fct_recode(stomatal_model, "PMAX" = "PROFITMAX","PMAX2"="PROFITMAX2")) %>% 
  mutate(stomatal_model = factor(stomatal_model, 
                         levels = rev(c("SOX2",
                                        "SOX",
                                        "PMAX3",
                                        "PMAX2",
                                        "PMAX",
                                        "CGAIN2",
                                        "CGAIN",
                                        "PHYDRO"))))

# viscosity_water = calc_viscosity_h2o(25, 101325)
# density_water = calc_density_h2o(25,101325)
# Kmax = 1e-16/viscosity_water
# Kmax = Kmax * density_water * 55.5
# Kmax = Kmax * 1e6
  
DF_example %>% 
  ggplot()+
  geom_line(aes(var,value,color=alpha_mod,linetype =p50_mod),show.legend = FALSE)+
  ggh4x::facet_grid2(c("stomatal_model","name"), 
                     labeller = label_parsed, 
                     scales = "free", 
                     independent = "all")+
  # facet_wrap(stomatal_model~name, scales = "free",ncol = 4)+
  mytheme5()+
  scale_color_grey()+
  theme(legend.position="bottom",
        legend.title = element_blank())+
  labs(x=expression(psi[s]~"(MPa)"))

ggsave("PLOTS/example_plot.pdf", width = 24, height = 30, units = "cm")






#### Example plot temperature rise #####
# par_plant1 = list(
#   conductivity = 1e-16,
#   psi50 = -1,                   #
#   b=2
# )
# df_T <- tibble(var=seq(20,30,length.out=10) ) %>% 
#   mutate(
#     VPD = bigleaf::rH.to.VPD(0.8,var)*1000,
#     out_hydraulics_1 = purrr::pmap(list(var, VPD),
#                                    ~model_numerical(tc = ..1, ppfd = 800,
#                                                     vpd = ..2, co2 = 400, elv = 0,
#                                                     fapar = 0.99, kphio = 0.087, 
#                                                     psi_soil = 0, par_plant = par_plant1, rdark = 0.02,
#                                                     par_cost = list(alpha=0.1, gamma = NA),
#                                                     stomatal_model = "PROFITMAX"
#   ))) %>% 
#   unnest_wider(out_hydraulics_1)
# 
# 
# df_T %>% 
#   ggplot(aes(var,jmax/vcmax))+
#   geom_line()
# 
# df_T %>% 
#   ggplot()+
#   geom_line(aes(var,vcmax))+
#   geom_line(aes(var,jmax), color = "red")+
#   geom_line(aes(var,jmax+vcmax), color = "blue")
# 
# 
# df_T %>%
#   rowwise() %>% 
#   mutate(visco_25 = calc_viscosity_h2o(25, 101325),
#          viscosity_water = calc_viscosity_h2o(var, 101325),
#          visco_star = viscosity_water/visco_25) %>%
#   ggplot()+
#   geom_line(aes(var,146*0.0005/visco_star*vcmax + 0.41/4*jmax))+
#   geom_line(aes(var,0.1*jmax), color = "red")+
#   geom_line(aes(var,seq(0.11,0.09, length.out=10)*jmax), color = "blue")
# 
# 
# df_T %>%
#   rowwise() %>% 
#   mutate(visco_25 = calc_viscosity_h2o(25, 101325),
#          viscosity_water = calc_viscosity_h2o(var, 101325),
#          visco_star = viscosity_water/visco_25) %>%
#   ggplot()+
#   geom_point(aes(seq(0.11,0.09, length.out=10)*jmax,146*0.001/visco_star*vcmax + 0.41/4*jmax, color =var))+
#   geom_abline(slope=1, intercept=0)









#### alpha-traits pmodel ANALYSIS ####
#l_b and P50 are moderately-strongly correlated (we use only P50)
# df_param_kmax_alpha[is.na(df_param_kmax_alpha$lifeform),"lifeform"] <- "Evergreen needleleaf gymnosperm"
# alpha_mod <- lm(b_cost~ SLA..cm2.g.1.+ #KL+ #lifeform_comp +
#                               Iabs_growth+vcmax_ww+
#                               jmax_ww+hv+ #p50_opt + l_b_opt +
#                               Ta + ca+ P50+  Height.max..m.,
#                             data = df_param_kmax_alpha_pmodel %>%
#                               rowwise() %>%
#                               mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#                                      hv = log(Huber.value*1e4),
#                                      lk.scale =log(K.scale),
#                                      Height.max..m. = log( Height.max..m.),
#                                      l_b = log(b),
#                                      l_b_opt = log(b_opt),
#                                      l_p50_opt = log(-p50_opt),
#                                      D = D*101325/1000,
#                                      jmax_ww = log(jmax_ww),
#                                      Ta = T,
#                                      visco_25 = calc_viscosity_h2o(25, 101325),
#                                      viscosity_water = calc_viscosity_h2o(Ta, 101325),
#                                      visco_star = viscosity_water/visco_25,
#                                      b_cost = 146*alpha/visco_star),
#                             weights = log(n_dist)
# )
# summary(alpha_mod)
# step(alpha_mod)
# alpha_mod <- lm(formula = b_cost ~ Iabs_growth + vcmax_ww + hv + Ta + ca + 
#                   Height.max..m.,
#                 data = df_param_kmax_alpha_pmodel %>%
#                   rowwise() %>%
#                   mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
#                          hv = log(Huber.value*1e4),
#                          lk.scale =log(K.scale),
#                          Height.max..m. = log( Height.max..m.),
#                          l_b = log(b),
#                          l_b_opt = log(b_opt),
#                          l_p50_opt = log(-p50_opt),
#                          D = D*101325/1000,
#                          jmax_ww = log(jmax_ww),
#                          Ta = T,
#                          visco_25 = calc_viscosity_h2o(25, 101325),
#                          viscosity_water = calc_viscosity_h2o(Ta, 101325),
#                          visco_star = viscosity_water/visco_25,
#                          b_cost = 146*alpha/visco_star),
#                 weights = log(n_dist))
# # anova(alpha_mod,alpha_mod2)
# summary(alpha_mod)
# car::vif(alpha_mod)
# performance::check_collinearity(alpha_mod)
# confint(alpha_mod)
# MuMIn::r.squaredGLMM(alpha_mod)
# 
# 
# library(effects)
# closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
# 
# 
# eff<-effect("SLA..cm2.g.1.", partial.residuals=T, alpha_mod)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, SLA..cm2.g.1. = eff$x$SLA..cm2.g.1.)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$SLA..cm2.g.1.)] + eff$residuals)
# 
# g_eff1 <- ggplot(x, aes(x = SLA..cm2.g.1., y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression("SLA  [cm"^2~"g"^-1*"]"))
# 
# eff<-effect("Iabs_growth", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Iabs_growth = eff$x$Iabs_growth)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Iabs_growth)] + eff$residuals)
# 
# g_eff2 <- ggplot(x, aes(x = Iabs_growth, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(I[growth]~"["*mu*"mol"~m^-2~s^-1*"]"))
# 
# 
# # eff<-effect("hv", partial.residuals=T, alpha_mod)
# # # plot(eff)
# # x.fit <- unlist(eff$x.all)
# # trans <- I
# # x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, hv = eff$x$hv)
# # xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$hv)] + eff$residuals)
# #
# # g_eff3 <- ggplot(x, aes(x = hv, y = fit)) +
# #   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
# #   geom_line(size = 1) +
# #   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
# #   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
# #   geom_smooth(data = xy, aes(x = trans(x), y = y),
# #               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
# #   mytheme3()+
# #   ylab(expression(alpha))+
# #   xlab(expression("Huber value [ln("*cm[sw]^2~m[leaf]^-2*")]"))
# 
# 
# eff<-effect("Ta", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Ta = eff$x$Ta)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Ta)] + eff$residuals)
# 
# g_eff4 <- ggplot(x, aes(x = Ta, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(T[a]~"[ºC]"))
# 
# 
# # eff<-effect("P50", partial.residuals=T, alpha_mod)
# # # plot(eff)
# # x.fit <- unlist(eff$x.all)
# # trans <- I
# # x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, P50 = eff$x$P50)
# # xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$P50)] + eff$residuals)
# #
# # g_eff5 <- ggplot(x, aes(x = P50, y = fit)) +
# #   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
# #   geom_line(size = 1) +
# #   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
# #   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
# #   geom_smooth(data = xy, aes(x = trans(x), y = y),
# #               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
# #   mytheme3()+
# #   ylab(expression(alpha))+
# #   xlab(expression(psi[50]~"[MPa]"))
# #
# 
# eff<-effect("jmax_ww", partial.residuals=T, alpha_mod)
# # plot(eff)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, jmax_ww = eff$x$jmax_ww)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$jmax_ww)] + eff$residuals)
# 
# g_eff6 <- ggplot(x, aes(x = jmax_ww, y = fit)) +
#   geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
#   geom_line(size = 1) +
#   geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
#   geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
#   geom_smooth(data = xy, aes(x = trans(x), y = y),
#               method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   mytheme3()+
#   ylab(expression(alpha))+
#   xlab(expression(J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))
# 
# alpha_mod_plot <- ggarrange(g_eff2,g_eff4,g_eff6,g_eff1,#g_eff3,g_eff5,
#                             align='hv', labels=c('a', 'b','c','d'),
#                             ncol=2, nrow = 2)
# 
# ggsave("PLOTS/alpha_model_2.png", plot = alpha_mod_plot, width = 25, height = 20, units = "cm")
# 
# 
# 
# 






#### Actual Jmax estimated Jmax ####

# mod_vcmax <- lmer(vcmax_obs~vcmax25_pred*scheme + (vcmax25_pred|Species)+(1|source), data=df_kmax_alpha %>%
#                     filter(!is.na(vcmax_obs), acclimation == "Acclimated"),
#                   weights = log(n_dist)# source !="Galmes et al. (2007)")
#                            )
#        
# summary(mod_vcmax)
# emmeans::emtrends(mod_vcmax, var = "vcmax_pred",specs = "scheme")
# MuMIn::r.squaredGLMM(mod_vcmax)
# 
# df_kmax_alpha %>%
#   filter(!is.na(vcmax_obs), acclimation == "Acclimated") %>% 
#   ggplot(aes(vcmax25_pred, vcmax_obs,  color = Species))+
#   geom_point(aes(size=log(n_dist)), alpha=0.7, show.legend = FALSE, shape = 21)+
#   geom_smooth(aes(group=interaction(Species,scheme)), method = "lm", se=FALSE)+
#   geom_abline(slope = 1, intercept = 0, linetype = 3)+
#   facet_wrap(~scheme,nrow = 2)+
#   scale_color_viridis_d()+
#   # ylim(0,80)+
#   mytheme4()+
#   xlab(expression("Predicted"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   ylab(expression("Observed"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))
# 
# 
# # mod_jmax <- lmer(jmax_obs~jmax_pred*scheme + (jmax_pred|Species), data=df_kmax_alpha %>%
# #                     filter(!is.na(jmax_obs), acclimation == "Acclimated", source !="Galmes et al. (2007)"))
# # 
# # summary(mod_jmax)
# # emmeans::emtrends(mod_jmax, var = "jmax_pred",specs = "scheme")
# # MuMIn::r.squaredGLMM(mod_jmax)
# # 
# # df_kmax_alpha %>%
# #   filter(!is.na(jmax_obs), acclimation == "Acclimated", source !="Galmes et al. (2007)") %>% 
# #   ggplot(aes(jmax_obs, jmax_pred, color = Species))+
# #   geom_point()+
# #   geom_smooth(aes(group=interaction(Species,scheme)), method = "lm")+
# #   geom_abline(slope = 1, intercept = 0)+
# #   facet_grid(~scheme)
# 
# 
# 
# 
# 
# mod_vcmax <- lmer(vcmax_obs~vcmax_pred*scheme + (vcmax_pred|Species)+(1|source), data=df_kmaxww_a %>%
#                     filter(!is.na(vcmax_obs), acclimation == "Acclimated"),
#                   weights = log(n_dist)# source !="Galmes et al. (2007)")
# )
# 
# summary(mod_vcmax)
# emmeans::emtrends(mod_vcmax, var = "vcmax_pred",specs = "scheme")
# MuMIn::r.squaredGLMM(mod_vcmax)
# 
# df_kmaxww_a %>%
#   filter(!is.na(vcmax_obs), acclimation == "Acclimated") %>% 
#   ggplot(aes(vcmax_pred, vcmax_obs,  color = Species))+
#   geom_point(aes(size=log(n_dist)), alpha=0.7, show.legend = FALSE, shape = 21)+
#   geom_smooth(aes(group=interaction(Species,scheme)), method = "lm", se=FALSE)+
#   geom_abline(slope = 1, intercept = 0, linetype = 3)+
#   facet_wrap(~scheme,nrow = 2)+
#   scale_color_viridis_d()+
#   # ylim(0,80)+
#   mytheme4()+
#   xlab(expression("Predicted"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   ylab(expression("Observed"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))
# 
# 
# 
# 
# mod_vcmax <- lmer(vcmax_obs~vcmax_pred*scheme + (vcmax_pred|Species)+(1|source), data=df_a_fix %>% left_join(df_a_fix_n) %>% 
#                     filter(!is.na(vcmax_obs), acclimation == "Acclimated"),
#                   weights = log(n_dist)# source !="Galmes et al. (2007)")
# )
# 
# summary(mod_vcmax)
# emmeans::emtrends(mod_vcmax, var = "vcmax_pred",specs = "scheme")
# MuMIn::r.squaredGLMM(mod_vcmax)
# 
# df_a_fix %>% left_join(df_a_fix_n) %>%
#   filter(!is.na(vcmax_obs), acclimation == "Acclimated") %>% 
#   ggplot(aes(vcmax_pred, vcmax_obs,  color = Species))+
#   geom_point(aes(size=log(n_dist)), alpha=0.7, show.legend = FALSE, shape = 21)+
#   geom_smooth(aes(group=interaction(Species,scheme)), method = "lm", se=FALSE)+
#   geom_abline(slope = 1, intercept = 0, linetype = 3)+
#   facet_wrap(~scheme,nrow = 2)+
#   scale_color_viridis_d()+
#   # ylim(0,80)+
#   mytheme4()+
#   xlab(expression("Predicted"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   ylab(expression("Observed"~V[cmax]~"["*mu*"mol"~m^-2~s^-1*"]"))+
#   theme(legend.title = element_blank(),
#         legend.text = element_text(face = "italic"))
