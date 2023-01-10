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
library(tidyverse)
library(scales)
library(ggalt)
library(grid)
library(ggConvexHull)
library(rphydro)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(DEoptim)
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




#### SIMULATIONS ALPHA = 0.1 ####
load(file = "DATA/simulations_kmax.RData")
df_a_fix <- df %>% 
  mutate(chi = Ciest/ca,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated")))

df_a_fix_n <- df_a_fix %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(A)) %>% 
  summarise(n_dist = n())


## PARAMETERS DATA
path_par <- "DATA/parameters_kmax/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })

df_a_fix_param <- par_data %>% 
  mutate(acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated"))
  )

## CALCULATE WW JMAX AND VCMAX
vcmax_jmax_ww <- df_a_fix %>% 
  group_by(scheme,source,Species) %>% 
  do(get_vcmax_jmax_ww(.))


#### SIMULATIONS KMAXWW ALPHA ####
load(file = "DATA/simulations_kmaxww_alpha.RData")
df_kmaxww_a <- df %>% 
  mutate(chi = Ciest/ca,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated")))

## PARAMETERS DATA
path_par <- "DATA/parameters_kmaxww_alpha/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })

df_kmaxww_a_param <- par_data %>% 
  mutate(acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated"))
  )


#### KMAX ALPHA DATASET ####

path_par_kmax_alpha <- "DATA/parameters_kmax_alpha/"
par_data_kmax_alpha <- list.files(path_par_kmax_alpha) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par_kmax_alpha,x))
  })

df_param_kmax_alpha <- par_data_kmax_alpha %>% 
  mutate(stomatal_model = scheme,
         acclimation = factor(acclimation, 
                              levels = c('TRUE','FALSE'),
                              labels = c("Acclimated", "Not acclimated"))
         )%>% 
  # filter(
  #   !Species %in% c('Broussonetia papyrifera'))%>%
  left_join(vcmax_jmax_ww) %>% 
  filter(K.scale > 1e-04)

# df_param_kmax_alpha[which(df_param_kmax_alpha$Species %in% 
#                  c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
#          "Huber.value"] <- olea_hv$Huber.value
# df_param_kmax_alpha[which(df_param_kmax_alpha$Species %in% 
#                  c("Olea europaea var. Chemlali","Olea europaea var. Meski")),
#          "KL..kg.m.1.MPa.1.s.1."] <- olea_kl$KL..kg.m.1.MPa.1.s.1.

df_param_kmax_alpha <- df_param_kmax_alpha %>% 
  left_join(df_a_fix %>% 
              dplyr::select(scheme,source,Species,Iabs_growth,D,T,ca) %>%
              group_by(scheme,source,Species) %>%  
              summarise_all(mean,na.rm=TRUE)
            ) %>% 
  left_join(df_a_fix_n) %>%
  mutate(var = 0,
         psi0 = purrr::pmap(list(var, T,Iabs_growth,D,ca,K.scale,P50,b,alpha,gamma,stomatal_model), 
                            ~model_numerical(tc =..2, ppfd = ..3, 
                                          vpd = ..4*101325, co2 = ..5, 
                                          elv = 0, fapar = .99, kphio = 0.087, 
                                          psi_soil = ..1, rdark = 0.02, 
                                          par_plant=list(conductivity = ..6*1e-16,
                                                         psi50 = ..7%>% unique(),
                                                         b = ..8%>% unique()), 
                                          par_cost = list(
                                            alpha = ..9,
                                            gamma = ..10
                                          ), stomatal_model = ..11 %>% unique())))%>% 
  unnest_wider(psi0,names_sep = ".") %>% 
  mutate(var = LWP_ww,
         psi_ww = purrr::pmap(list(var, T,Iabs_growth,D,ca,K.scale,P50,b,alpha,gamma,stomatal_model), 
                            ~model_numerical(tc =..2, ppfd = ..3, 
                                             vpd = ..4*101325, co2 = ..5, 
                                             elv = 0, fapar = .99, kphio = 0.087, 
                                             psi_soil = ..1, rdark = 0.02, 
                                             par_plant=list(conductivity = ..6*1e-16,
                                                            psi50 = ..7%>% unique(),
                                                            b = ..8%>% unique()), 
                                             par_cost = list(
                                               alpha = ..9,
                                               gamma = ..10
                                             ), stomatal_model = ..11 %>% unique())))%>% 
  unnest_wider(psi_ww,names_sep = ".")


#### Jmax vs alpha and Jmax vs Jmax ####

## jmax0
jmax0_alpha_mod <- lmerTest::lmer(alpha~log(psi0.jmax)+  
                                    (1|scheme) + (1|source/Species), 
                                  data = df_param_kmax_alpha, weights = log(n_dist))
anova(jmax0_alpha_mod)
summary(jmax0_alpha_mod)
r2 <- MuMIn::r.squaredGLMM(jmax0_alpha_mod)
r2 <- round(r2[1],2)
r2
jmax0_alpha <- summary(jmax0_alpha_mod)$coefficients
jmax0_sim <- seq(min(log(df_param_kmax_alpha$psi0.jmax),na.rm = TRUE),
                max(log(df_param_kmax_alpha$psi0.jmax),na.rm = TRUE),
                0.01)
alpha_sim_jmax0 <- jmax0_alpha[1,1]+jmax0_alpha[2,1]*jmax0_sim
alpha_sd_max_jmax0 <- jmax0_alpha[1,1] + jmax0_alpha[1,2] + (jmax0_alpha[2,1]+jmax0_alpha[2,2])*jmax0_sim
alpha_sd_min_jmax0 <- jmax0_alpha[1,1] - jmax0_alpha[1,2] + (jmax0_alpha[2,1]-jmax0_alpha[2,2])*jmax0_sim
resume_jmax0_alpha <- df_param_kmax_alpha %>% 
  group_by(Species,source) %>% 
  mutate(ln_psi0.jmax = log(psi0.jmax)) %>% 
  dplyr::select(ln_psi0.jmax,alpha) %>% 
  summarise_all(tibble::lst(mean,sd))
(p1 <- ggplot()+
  geom_point(data = df_param_kmax_alpha %>% mutate(ln_psi0.jmax = log(psi0.jmax)),
               mapping=aes(x=ln_psi0.jmax,alpha,group=interaction(Species,source),size=log(n_dist)),
               color="grey60", show.legend = FALSE,alpha=0.5)+
  geom_point(resume_jmax0_alpha,mapping = aes(x=ln_psi0.jmax_mean, y = alpha_mean),color="grey20")+
  geom_errorbar(resume_jmax0_alpha,
                mapping = aes(x=ln_psi0.jmax_mean, ymin = alpha_mean-alpha_sd,
                              ymax = alpha_mean+alpha_sd),
                color="grey20")+
  geom_errorbar(resume_jmax0_alpha,
                mapping = aes(y=alpha_mean, xmin = ln_psi0.jmax_mean-ln_psi0.jmax_sd,
                              xmax = ln_psi0.jmax_mean+ln_psi0.jmax_sd),
                color="grey20")+
  geom_line(aes(x = jmax0_sim, y = alpha_sim_jmax0),color = "grey20")+
  geom_line(aes(x = jmax0_sim, y = alpha_sd_max_jmax0),color = "grey20",linetype=2)+
  geom_line(aes(x = jmax0_sim, y = alpha_sd_min_jmax0),color = "grey20",linetype=2)+
  annotate(geom = "text",x=5,y=0.155,hjust=0,vjust=0, parse = TRUE,
           label=as.character(as.expression("***"~~italic(R)^2~"="~0.33)),size=4)+
  # stat_summary(data = df_param_kmax_alpha,
  #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
  # geom="pointrange", color="blue")+
  mytheme()+
  # scale_radius(range=c(0.2,1))+
  xlab(expression(J[max0]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
  ylab(expression(alpha))
)


##jmaxww_modelled
df_params_kmax_alpha_filtered <- df_param_kmax_alpha %>% 
  mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
  filter(ln_psi_ww.jmax>0)
jmax_ww_alpha_mod <- lmerTest::lmer(alpha~ln_psi_ww.jmax+I(ln_psi_ww.jmax^2)+ 
                                    (1|scheme) + (1|source/Species), 
                                  data = df_params_kmax_alpha_filtered,weights = log(n_dist))
anova(jmax_ww_alpha_mod)
summary(jmax_ww_alpha_mod)
r2 <- MuMIn::r.squaredGLMM(jmax_ww_alpha_mod)
r2 <- round(r2[1],2)
r2
jmax_ww_alpha <- summary(jmax_ww_alpha_mod)$coefficients
jmax_ww_sim <- seq(min(log(df_params_kmax_alpha_filtered$psi_ww.jmax),na.rm = TRUE),
                 max(log(df_params_kmax_alpha_filtered$psi_ww.jmax),na.rm = TRUE),
                 0.01)
alpha_sim_jmax_ww <- jmax_ww_alpha[1,1]+jmax_ww_alpha[2,1]*jmax_ww_sim+
  jmax_ww_alpha[3,1]*jmax_ww_sim^2

alpha_sd_max_jmax_ww <- jmax_ww_alpha[1,1] + jmax_ww_alpha[1,2] + 
  (jmax_ww_alpha[2,1]-jmax_ww_alpha[2,2])*jmax_ww_sim+ 
  (jmax_ww_alpha[3,1]+jmax_ww_alpha[3,2])*jmax_ww_sim^2

alpha_sd_min_jmax_ww <- jmax_ww_alpha[1,1] - jmax_ww_alpha[1,2] + 
  (jmax_ww_alpha[2,1]+jmax_ww_alpha[2,2])*jmax_ww_sim+ 
  (jmax_ww_alpha[3,1]-jmax_ww_alpha[3,2])*jmax_ww_sim^2

resume_jmax_ww_alpha <- df_param_kmax_alpha %>% 
  group_by(Species,source) %>% 
  mutate(ln_psi_ww.jmax = log(psi_ww.jmax)) %>% 
  filter(ln_psi_ww.jmax>0)%>% 
  dplyr::select(ln_psi_ww.jmax,alpha) %>% 
  summarise_all(tibble::lst(mean,sd))
(p2 <- ggplot()+
  geom_point(data = df_param_kmax_alpha %>% 
               mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
               filter(ln_psi_ww.jmax>0),
             mapping=aes(x=ln_psi_ww.jmax,alpha,group=interaction(Species,source),size=log(n_dist)),
             color="grey60", show.legend = FALSE,alpha=0.5)+
  geom_point(resume_jmax_ww_alpha,mapping = aes(x=ln_psi_ww.jmax_mean, y = alpha_mean),color="grey20")+
  # geom_smooth(data = df_param_kmax_alpha %>% 
  #               mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
  #               filter(ln_psi_ww.jmax>0),
  #             mapping=aes(x=ln_psi_ww.jmax,alpha),
  #             method = "lm",formula = y~poly(x,2))+
  geom_errorbar(resume_jmax_ww_alpha,
                mapping = aes(x=ln_psi_ww.jmax_mean, ymin = alpha_mean-alpha_sd,
                              ymax = alpha_mean+alpha_sd),
                color="grey20")+
  geom_errorbar(resume_jmax_ww_alpha,
                mapping = aes(y=alpha_mean, xmin = ln_psi_ww.jmax_mean-ln_psi_ww.jmax_sd,
                              xmax = ln_psi_ww.jmax_mean+ln_psi_ww.jmax_sd),
                color="grey20")+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sim_jmax_ww),color = "grey20")+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sd_max_jmax_ww),color = "grey20",linetype=2)+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sd_min_jmax_ww),color = "grey20",linetype=2)+
  annotate(geom = "text",x=5,y=0.155,hjust=0,vjust=0, parse = TRUE,
           label=as.character(as.expression("***"~~italic(R)^2~"="~0.33)),size=4)+
  # stat_summary(data = df_param_kmax_alpha,
  #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
  # geom="pointrange", color="blue")+
  mytheme3()+
  # scale_radius(range=c(0.2,1))+
  xlab(expression(J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
  ylab(expression(alpha))
)


##jmaxww observed
# df_params_kmax_alpha_filtered <- df_param_kmax_alpha %>% 
#   mutate(ln_psi_ww.jmax = log(psi_ww.jmax))%>% 
#   filter(ln_psi_ww.jmax>0)
jmax_ww_alpha_mod <- lmerTest::lmer(alpha~log(jmax_ww)+ 
                                      (1|scheme) + (1|source/Species), 
                                    data = df_param_kmax_alpha, weights = log(n_dist))
anova(jmax_ww_alpha_mod)
summary(jmax_ww_alpha_mod)
r2 <- MuMIn::r.squaredGLMM(jmax_ww_alpha_mod)
r2 <- round(r2[1],2)
r2
jmax_ww_alpha <- summary(jmax_ww_alpha_mod)$coefficients
jmax_ww_sim <- seq(min(log(df_param_kmax_alpha$jmax_ww),na.rm = TRUE),
                   max(log(df_param_kmax_alpha$jmax_ww),na.rm = TRUE),
                   0.01)
alpha_sim_jmax_ww <- jmax_ww_alpha[1,1]+jmax_ww_alpha[2,1]*jmax_ww_sim
alpha_sd_max_jmax_ww <- jmax_ww_alpha[1,1] + jmax_ww_alpha[1,2] + (jmax_ww_alpha[2,1]+jmax_ww_alpha[2,2])*jmax_ww_sim
alpha_sd_min_jmax_ww <- jmax_ww_alpha[1,1] - jmax_ww_alpha[1,2] + (jmax_ww_alpha[2,1]-jmax_ww_alpha[2,2])*jmax_ww_sim
resume_jmax_ww_alpha <- df_param_kmax_alpha %>% 
  group_by(Species,source) %>% 
  mutate(ln_jmax_ww = log(jmax_ww)) %>% 
  filter(ln_jmax_ww>0)%>% 
  dplyr::select(ln_jmax_ww,alpha) %>% 
  summarise_all(tibble::lst(mean,sd))
(p3 <- ggplot()+
  geom_point(data = df_param_kmax_alpha %>% 
               mutate(ln_jmax_ww = log(jmax_ww))%>% 
               filter(ln_jmax_ww>0),
             mapping=aes(x=ln_jmax_ww,alpha,group=interaction(Species,source),size=log(n_dist)),
             color="grey60", show.legend = FALSE,alpha=0.5)+
  geom_point(resume_jmax_ww_alpha,mapping = aes(x=ln_jmax_ww_mean, y = alpha_mean),color="grey20", size = 3)+
  geom_errorbar(resume_jmax_ww_alpha,
                mapping = aes(x=ln_jmax_ww_mean, ymin = alpha_mean-alpha_sd,
                              ymax = alpha_mean+alpha_sd),
                color="grey20")+
  geom_errorbar(resume_jmax_ww_alpha,
                mapping = aes(y=alpha_mean, xmin = ln_jmax_ww_mean-ln_jmax_ww_sd,
                              xmax = ln_jmax_ww_mean+ln_jmax_ww_sd),
                color="grey20")+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sim_jmax_ww),color = "grey20")+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sd_max_jmax_ww),color = "grey20",linetype=2)+
  geom_line(aes(x = jmax_ww_sim, y = alpha_sd_min_jmax_ww),color = "grey20",linetype=2)+
  annotate(geom = "text",x=4.5,y=0.155,hjust=0,vjust=0, parse = TRUE,
           label=as.character(as.expression("***"~~italic(R)^2~"="~0.41)),size=5)+
  # stat_summary(data = df_param_kmax_alpha,
  #              mapping=aes(log(jmax_ww),alpha,group=interaction(species,source)),fun.data=mean_sdl, fun.args = list(mult=1),
  # geom="pointrange", color="blue")+
  mytheme5()+
  # scale_radius(range=c(0.2,1))+
  xlab(expression("Observed"~J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
  ylab(expression("Calibrated"~alpha))
)


(p4 <- ggplot(data = df_param_kmax_alpha,
       mapping=aes(log(psi0.jmax),log(jmax_ww),color=scheme))+
       geom_point()+
       geom_abline(slope=1,intercept=0,color="grey20", linetype=2)+
       geom_smooth(method="lm",se=FALSE)+
       mytheme5()+
       scale_colour_manual(breaks = col_df$scheme, 
                          values = unique(as.character(col_df$col)))+
       theme(legend.title = element_blank())+
       xlab(expression("Modelled"~J[max0]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))+
       ylab(expression("Observed"~J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))
  )
  
(p5 <- ggplot(data = df_param_kmax_alpha,
       mapping=aes(psi_ww.jmax,jmax_ww,color=scheme, weight=log(n_dist)))+
       geom_point(aes(size =log(n_dist)),alpha = 0.5,show.legend = FALSE)+
       geom_abline(slope=1,intercept=0,color="grey20", linetype=3,show.legend = FALSE)+
       geom_smooth(method="lm",se = FALSE)+
       mytheme5()+
       scale_colour_manual(breaks = col_df$scheme, 
                           values = unique(as.character(col_df$col)))+
       theme(legend.title = element_blank())+
       xlab(expression("Modelled"~J[maxWW]~"["*mu*"mol"~m^-2~s^-1*"]"))+
       ylab(expression("Observed"~J[maxWW]~"["*mu*"mol"~m^-2~s^-1*"]"))
  )

jmax_jmax_plot <- ggarrange(p3,p5,
                            align='h', labels=c('a', 'b'),
                            common.legend = T,legend="bottom",
                            ncol=1, nrow = 2)

ggsave("PLOTS/alpha_jmax_jmax_model.png",jmax_jmax_plot, width = 14, height = 28, units = "cm")



#### Alpha multivariate analysis ####

library(psych)
pairs.panels(df_param_kmax_alpha %>%
        mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
               hv = log(Huber.value),
               lk.scale =log(K.scale),
               l_b = log(b),
               l_b_opt = log(b_opt),
               l_p50_opt = log(-p50_opt),
               D = D*101325/1000) %>% 
        dplyr::select(SLA..cm2.g.1.,
                      Iabs_growth,D, l_b , lk.scale,P50,
                      hv, l_b_opt, l_p50_opt,
                      Height.max..m.))

alpha_mod <- lmerTest::lmer(alpha~ SLA..cm2.g.1.+
                              Iabs_growth+jmax_ww+hv+ 
                              l_b +P50+  Height.max..m.+(1|Species)+ (1|scheme), 
                            data = df_param_kmax_alpha %>%
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                                     hv = log(Huber.value),
                                     lk.scale =log(K.scale),
                                     l_b = log(b),
                                     D = D*101325/1000,
                                     jmax_ww = jmax_ww), 
                            weights = log(n_dist)
                            )
summary(alpha_mod)
step(alpha_mod)
alpha_mod <- lmerTest::lmer(alpha~ Iabs_growth + jmax_ww + hv + l_b + 
                              P50 + Height.max..m. + 
                              (1 | Species) + (1 | scheme),
                            # alpha ~  jmax_ww + l_b + P50 + (1 |source/Species) + (1 | scheme),
                            data = df_param_kmax_alpha %>%
                              mutate(KL = log(KL..kg.m.1.MPa.1.s.1.*55.5),
                                     hv = log(Huber.value),
                                     lk.scale =log(K.scale),
                                     l_b = log(b),
                                     D = D*101325/1000,
                                     jmax_ww = jmax_ww), 
                            weights = log(n_dist)
                            )
summary(alpha_mod)
car::vif(alpha_mod)
confint(alpha_mod)
MuMIn::r.squaredGLMM(alpha_mod)


library(effects)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)


eff<-effect("Iabs_growth", partial.residuals=T, alpha_mod)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Iabs_growth = eff$x$Iabs_growth)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Iabs_growth)] + eff$residuals)

g_eff1 <- ggplot(x, aes(x = Iabs_growth, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y), 
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression(I[growth]~"["*mu*"mol"~m^-2~s^-1*"]"))

eff<-effect("Height.max..m.", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Height.max..m. = eff$x$Height.max..m.)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Height.max..m.)] + eff$residuals)

g_eff2 <- ggplot(x, aes(x = Height.max..m., y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y),
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab("Maximum species height [m]")


eff<-effect("hv", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, hv = eff$x$hv)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$hv)] + eff$residuals)

g_eff3 <- ggplot(x, aes(x = hv, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y), 
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression("Huber value [ln("*cm[sw]^2~m[leaf]^-2*")]"))


eff<-effect("l_b", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, lb = eff$x$l_b)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$lb)] + eff$residuals)

g_eff4 <- ggplot(x, aes(x = lb, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y), 
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression("Parameter b [ln(b)]"))


eff<-effect("P50", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, P50 = eff$x$P50)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$P50)] + eff$residuals)

g_eff5 <- ggplot(x, aes(x = P50, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y), 
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression(psi[50]~"[MPa]"))


eff<-effect("jmax_ww", partial.residuals=T, alpha_mod)
# plot(eff)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, jmax_ww = eff$x$jmax_ww)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$jmax_ww)] + eff$residuals)

g_eff6 <- ggplot(x, aes(x = jmax_ww, y = fit)) +
  geom_point(data = xy, aes(x = x, y = y), col = "grey60", size = 2) +
  geom_line(size = 1) +
  geom_line(aes(y= lower), alpha = 0.5,linetype=2) +
  geom_line(aes(y= upper), alpha = 0.5,linetype=2) +
  geom_smooth(data = xy, aes(x = trans(x), y = y), 
              method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  mytheme3()+
  ylab(expression(alpha))+
  xlab(expression(J[maxWW]~"[ln("*mu*"mol"~m^-2~s^-1*")]"))

alpha_mod_plot <- ggarrange(g_eff1,g_eff6,g_eff2,g_eff3,g_eff4,g_eff5,
          align='hv', labels=c('a', 'b','c','d','e','f'),
          ncol=3, nrow = 2)

ggsave("PLOTS/alpha_model.png", plot = alpha_mod_plot, width = 28, height = 14, units = "cm")


alpha_scheme <- lmerTest::lmer(alpha~ scheme + (1|Species)+ (1|source), 
                              data = df_param_kmax_alpha %>% 
                                mutate(scheme = factor(scheme, 
                                                       levels = c("SOX",
                                                                  "PROFITMAX",
                                                                  "PROFITMAX2",
                                                                  "CGAIN",
                                                                  "CMAX",
                                                                  "PHYDRO"))))
alpha_scheme_mean <- lmerTest::lmer(alpha~ (1|Species)+ (1|source), 
                               data = df_param_kmax_alpha )
summary(alpha_scheme)
model_means <- emmeans(alpha_scheme, "scheme")
  
  model_means_cld <- cld(object = model_means,
                         adjust = "sidak",
                         Letters = letters,
                         alpha = 0.05)


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

df_a <- df_a_fix %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2]) 



r_a <- lmerTest::lmer(r~scheme*acclimation + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                      )
r_a_p <- emmeans::contrast(emmeans(r_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_a <- emmeans(r_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_a %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 0, slope = 0, color = "grey20")+
  geom_pointrange(data = r_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(italic(A)~" r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()


##BETA
beta_a <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_a%>% 
                           filter(beta<= 40), weights = log(n_dist)
                         )
beta_a_p <- emmeans::contrast(emmeans(beta_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_a <- emmeans(beta_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
# beta_a_accl <- lmerTest::lmer(beta~scheme + (1|Species), data = df_a %>% 
#                                 filter(beta<= 40)%>% filter(acclimation=="Acclimated"),
#                               weights = log(n_dist))
# summary(beta_a_accl)
# beta_a_p_accl <- emmeans(beta_a_accl, "scheme")
# plot(beta_a_p_accl)
p2 <- df_a%>% 
  filter(beta<= 40)%>% 
  ggplot(aes(scheme,beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 1, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = beta_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(-10,30)+
  ylab(expression(italic(A)~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(19))+
  coord_flip()+
  # scale_y_log10()+
  NULL


##RMSE
rmse_a <- lmerTest::lmer(rmse~scheme*acclimation + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                         )
rmse_a_p <- emmeans::contrast(emmeans(rmse_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_a <- emmeans(rmse_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

rmse_a_accl <- lmerTest::lmer(rmse~scheme + (1|Species), data = df_a %>% filter(acclimation=="Acclimated"), weights = log(n_dist))
summary(rmse_a_accl)
rmse_a_p_accl <- emmeans(rmse_a_accl, "scheme")

p3 <- df_a %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_pointrange(data = rmse_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(italic(A)~" RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_a <- lmerTest::lmer(bias~scheme*acclimation + (1|source)+ (1|Species), data = df_a, weights = log(n_dist)
                         )
bias_a_p <- emmeans::contrast(emmeans(bias_a, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_a <- emmeans(bias_a,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_a_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_a %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 0, slope = 0, color = "grey20", linetype=3)+
  geom_pointrange(data = bias_a, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(italic(A)~" Bias"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,
          p2,p3,p4, 
          align='hv', labels=c('a', 'b','c','d'
                               ),
          common.legend = T,ncol=2, nrow = 2)

ggsave("PLOTS/A_metrics.png", width = 20, height = 14, units = "cm")


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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("PROFITMAX2 Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 
# df_a_fix %>% 
#   filter(scheme == "SOX") %>% 
#   ggplot(aes(a_pred,A,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression(textstyle("Predicted A")))+
#   ylab(expression(textstyle("Observed A")))+
#   ggtitle(expression(atop("CMAX Assimilation rate,", italic(A)*" ("*mu*"mol m"^-2~"s"^-1*")")))
# 




#### G ####

df_g <- df_a_fix %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(gC)) %>% 
  mutate(diff_g = gC - g_pred) %>% 
  summarise(n_dist = n(),
            r = cor(gC, g_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_g,na.rm = TRUE)/mean(gC,na.rm = TRUE),
            rmse = Metrics::rmse(gC,g_pred),
            beta = lm(gC~g_pred)$coefficients[2]) 



r_g <- lmerTest::lmer(r~scheme*acclimation + (1|source)+ (1|Species), data = df_g, weights = log(n_dist))
r_g_p <- emmeans::contrast(emmeans(r_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
r_g <- emmeans(r_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE)%>% 
  left_join(r_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p1 <- df_g %>%
  ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 1, slope = 0, color = "grey20",linetype =3)+
  geom_pointrange(data = r_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange", position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(italic(g)[s]~" r Pearson's correlation"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



beta_g <- lmerTest::lmer(beta~scheme*acclimation + (1|source)+ (1|Species), data = df_g%>% 
                           filter(beta>-10, beta<20), weights = log(n_dist))
beta_g_p <- emmeans::contrast(emmeans(beta_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
beta_g <- emmeans(beta_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(beta_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 

p2 <- df_g%>% 
  filter(beta>-10, beta<20)%>% 
  ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 1, slope = 0, color = "grey20",linetype =3)+
  geom_pointrange(data = beta_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  # ylim(0,2)+
  ylab(expression(italic(g)[s]~beta))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1, 19))+
  coord_flip()



rmse_g <- lmerTest::lmer(rmse~scheme*acclimation + (1|source)+ (1|Species), data = df_g, weights = log(n_dist))
rmse_g_p <- emmeans::contrast(emmeans(rmse_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
rmse_g <- emmeans(rmse_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(rmse_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p3 <- df_g %>%
  ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_pointrange(data = rmse_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig),size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylab(expression(italic(g)[s]~"RMSE  ("*mu*"mol m"^-2~"s"^-1*")"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()

bias_g <- lmerTest::lmer(bias~scheme*acclimation + (1|source)+ (1|Species), data = df_g, weights = log(n_dist))
bias_g_p <- emmeans::contrast(emmeans(bias_g, "acclimation",by='scheme'))%>% 
  broom::tidy() %>% 
  dplyr::select(scheme,adj.p.value) %>% 
  summarise_all(unique)
bias_g <- emmeans(bias_g,~scheme*acclimation) %>% 
  broom::tidy(conf.int = TRUE) %>% 
  left_join(bias_g_p) %>% 
  mutate(sig = case_when(adj.p.value >= 0.05~"NO",
                         TRUE~"YES")) 
p4 <- df_g %>%
  ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
  geom_point(shape= 21,position=position_dodge(width = 0.4))+
  geom_abline(intercept = 0, slope = 0, color = "grey20",linetype =3)+
  geom_pointrange(data = bias_g, 
                  aes(scheme,estimate,color = acclimation,
                      ymin = conf.low,ymax = conf.high,
                      shape = sig), size = 0.8,
                  position=position_dodge(width = 0.4),
                  show.legend = FALSE)+
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="pointrange",position=position_dodge(width = 0.4))+
  mytheme()+
  theme(legend.title = element_blank())+
  xlab("")+
  ylim(-0.45,0.45)+
  ylab(expression(italic(g)[s]~"Bias"))+
  scale_color_manual(values = c("#A6611A","#018571"))+
  scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
  scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
  coord_flip()


ggarrange(p1,p2,p3,p4, 
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
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
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted g"[s]))+
#   ylab(expression("Observed g"[s]))+
#   ggtitle(expression(atop("CMAX Stomatal conductance,","g"[s]~"(mol m"^-2~"s"^-1*")")))
# 




# #### CHI ####
# df_c <- df_a_fix %>% 
#   group_by(scheme,acclimation,Species,source) %>% 
#   filter(!is.na(chi)) %>% 
#   mutate(diff_c = chi - c_pred) %>% 
#   summarise(n_dist = n(),
#             r = cor(chi, c_pred, use = "pairwise.complete.obs"),
#             bias = mean(diff_c,na.rm = TRUE)/mean(chi,na.rm = TRUE),
#             rmse = Metrics::rmse(chi,c_pred),
#             beta = lm(chi~c_pred)$coefficients[2]) %>% 
#   filter(beta<10)
# 
# 
# 
# r_c <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_c, weights = log(n_dist))
# r_c_p <- emmeans::contrast(emmeans(r_c, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# r_c <- emmeans(r_c,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE)%>% 
#   left_join(r_c_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p1 <- df_c %>%
#   ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 0, slope = 0, color = "grey20")+
#   geom_pointrange(data = r_c, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange", position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio r Pearson's correlation"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1, 19))+
#   coord_flip()
# 
# 
# 
# beta_c <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_c, weights = log(n_dist))
# beta_c_p <- emmeans::contrast(emmeans(beta_c, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# beta_c <- emmeans(beta_c,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(beta_c_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# 
# p2 <- df_c%>% 
#   ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 1, slope = 0, color = "grey20")+
#   geom_pointrange(data = beta_c, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   # ylim(0,2)+
#   ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio"~beta))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1, 19))+
#   coord_flip()
# 
# 
# 
# rmse_c <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_c, weights = log(n_dist))
# rmse_c_p <- emmeans::contrast(emmeans(rmse_c, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# rmse_c <- emmeans(rmse_c,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(rmse_c_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p3 <- df_c %>%
#   ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_pointrange(data = rmse_c, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio RMSE"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
#   coord_flip()
# 
# bias_c <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_c, weights = log(n_dist))
# bias_c_p <- emmeans::contrast(emmeans(bias_c, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# bias_c <- emmeans(bias_c,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(bias_c_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p4 <- df_c %>%
#   ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 0, slope = 0, color = "grey20")+
#   geom_pointrange(data = bias_c, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig), size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Leaf internal-to-ambient CO"[2]~"ratio BIAS"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
#   coord_flip()
# 
# 
# ggarrange(p1,p2,p3,p4, 
#           align='hv', labels=c('a', 'b','c','d'),
#           common.legend = T,ncol=2, nrow = 2)

# 
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX") %>% 
#   ggplot(aes(c_pred,chi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~chi))+
#   ylab(expression("Observed "~chi))+
#   ggtitle(expression(atop("PROFITMAX Leaf internal-to-ambient CO"[2]~"ratio")))
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX2") %>% 
#   ggplot(aes(c_pred,chi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~chi))+
#   ylab(expression("Observed "~chi))+
#   ggtitle(expression(atop("PROFITMAX2 Leaf internal-to-ambient CO"[2]~"ratio")))
# 
# df_a_fix %>% 
#   filter(scheme == "SOX") %>% 
#   ggplot(aes(c_pred,chi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~chi))+
#   ylab(expression("Observed "~chi))+
#   ggtitle(expression(atop("SOX Leaf internal-to-ambient CO"[2]~"ratio")))
# 
# df_a_fix %>% 
#   filter(scheme == "PHYDRO") %>% 
#   ggplot(aes(c_pred,chi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~chi))+
#   ylab(expression("Observed "~chi))+
#   ggtitle(expression(atop("PHYDRO Leaf internal-to-ambient CO"[2]~"ratio")))
# 
# 
# 
# 
# #### DPSI ####
# df_d <- df_a_fix %>% 
#   group_by(scheme,acclimation,Species,source) %>% 
#   filter(!is.na(Dpsi)) %>% select(Dpsi,d_pred)%>% 
#   mutate(diff_d = Dpsi - d_pred) %>% 
#   summarise(n_dist = n(),
#             r = cor(Dpsi, d_pred, use = "pairwise.complete.obs"),
#             bias = mean(diff_d,na.rm = TRUE)/mean(Dpsi,na.rm = TRUE),
#             rmse = Metrics::rmse(Dpsi,d_pred),
#             beta = lm(Dpsi~d_pred)$coefficients[2]) %>% 
#   filter(beta> (-50),beta <50, rmse<100)
# 
# 
# 
# r_d <- lmerTest::lmer(r~scheme*acclimation + (1|Species), data = df_d, weights = log(n_dist))
# r_d_p <- emmeans::contrast(emmeans(r_d, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# r_d <- emmeans(r_d,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE)%>% 
#   left_join(r_d_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p1 <- df_d %>%
#   ggplot(aes(scheme,r,fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 0, slope = 0, color = "grey20")+
#   geom_pointrange(data = r_d, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange", position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~"ratio r Pearson's correlation"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1, 19))+
#   coord_flip()
# 
# 
# 
# beta_d <- lmerTest::lmer(beta~scheme*acclimation + (1|Species), data = df_d, weights = log(n_dist))
# beta_d_p <- emmeans::contrast(emmeans(beta_d, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# beta_d <- emmeans(beta_d,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(beta_d_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# 
# p2 <- df_d%>% 
#   ggplot(aes(scheme, beta, fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 1, slope = 0, color = "grey20")+
#   geom_pointrange(data = beta_d, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   # ylim(0,2)+
#   ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~beta))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1, 19))+
#   coord_flip()
# 
# 
# 
# rmse_d <- lmerTest::lmer(rmse~scheme*acclimation + (1|Species), data = df_d, weights = log(n_dist))
# rmse_d_p <- emmeans::contrast(emmeans(rmse_d, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# rmse_d <- emmeans(rmse_d,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(rmse_d_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p3 <- df_d %>%
#   ggplot(aes(scheme,rmse, fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_pointrange(data = rmse_d, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig),size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~"(MPa) RMSE"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
#   coord_flip()
# 
# bias_d <- lmerTest::lmer(bias~scheme*acclimation + (1|Species), data = df_d, weights = log(n_dist))
# bias_d_p <- emmeans::contrast(emmeans(bias_d, "acclimation",by='scheme'))%>% 
#   broom::tidy() %>% 
#   dplyr::select(scheme,adj.p.value) %>% 
#   summarise_all(unique)
# bias_d <- emmeans(bias_d,~scheme*acclimation) %>% 
#   broom::tidy(conf.int = TRUE) %>% 
#   left_join(bias_d_p) %>% 
#   mutate(sig = case_when(adj.p.value >= 0.05~"NO",
#                          TRUE~"YES")) 
# p4 <- df_d %>%
#   ggplot(aes(scheme,bias,fill = acclimation, color = acclimation, group = acclimation))+
#   geom_point(shape= 21,position=position_dodge(width = 0.4))+
#   geom_abline(intercept = 0, slope = 0, color = "grey20")+
#   geom_pointrange(data = bias_d, 
#                   aes(scheme,estimate,color = acclimation,
#                       ymin = conf.low,ymax = conf.high,
#                       shape = sig), size = 0.8,
#                   position=position_dodge(width = 0.4),
#                   show.legend = FALSE)+
#   # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#   #              geom="pointrange",position=position_dodge(width = 0.4))+
#   mytheme2()+
#   theme(legend.title = element_blank())+
#   xlab("")+
#   ylab(expression("Soil-leaf water-potential difference"~italic(Delta*psi)~" BIAS"))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   scale_fill_manual(values =c("#DFC27D","#80CDC1"))+
#   scale_shape_manual(values = c(1,19))+ #since all the pairs have significant difference, select only sig shape
#   coord_flip()
# 
# 
# ggarrange(p1,p2,p3,p4, 
#           align='hv', labels=c('a', 'b','c','d'),
#           common.legend = T,ncol=2, nrow = 2)
# 
# 
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX") %>% 
#   ggplot(aes(d_pred,Dpsi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~Delta*psi))+
#   ylab(expression("Observed "~Delta*psi))+
#   ggtitle(expression(atop("PROFITMAX Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))
# 
# df_a_fix %>% 
#   filter(scheme == "PROFITMAX2") %>% 
#   ggplot(aes(d_pred,Dpsi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~Delta*psi))+
#   ylab(expression("Observed "~Delta*psi))+
#   ggtitle(expression(atop("PROFITMAX2 Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))
# 
# df_a_fix %>% 
#   filter(scheme == "SOX") %>% 
#   ggplot(aes(d_pred,Dpsi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~Delta*psi))+
#   ylab(expression("Observed "~Delta*psi))+
#   ggtitle(expression(atop("SOX Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))
# 
# df_a_fix %>% 
#   filter(scheme == "PHYDRO") %>% 
#   ggplot(aes(d_pred,Dpsi,color = acclimation))+
#   geom_point()+
#   geom_abline(intercept = 0,slope = 1, color = "grey20")+
#   geom_smooth(method = "lm", se = FALSE)+
#   facet_wrap(~Species)+
#   mytheme2()+
#   theme(legend.title = element_blank(),
#         legend.position="top",
#         plot.title = element_text(vjust = -10))+
#   scale_color_manual(values = c("#A6611A","#018571"))+
#   xlab(expression("Predicted "~Delta*psi))+
#   ylab(expression("Observed "~Delta*psi))+
#   ggtitle(expression(atop("PHYDRO Soil-leaf water-potential","difference,"~italic(Delta*psi)~"(MPa)")))
# 
# 
# 

# STOMATAL NON-STOMATAL PARTITION ----------------------------------------------
get_partition <- function(x){
  mod_not <- lmer(A~a_pred + (1|Species),data = x %>% filter(acclimation == "Not acclimated"))
  r2_not <- MuMIn::r.squaredGLMM(mod_not)
  mod <- lmer(A~a_pred + (1|Species),data = x %>% filter(acclimation == "Acclimated"))
  r2 <- MuMIn::r.squaredGLMM(mod)
  
  return(tibble(stomatal_r2 = r2_not[1],
             non_stomatal_r2 = r2[1]-r2_not[1],
             species_r2 = r2[2]-r2[1],
             Residuals = 1-r2[2]))
}

partition_full <- df_kmaxww_a %>% 
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>% 
  do(get_partition(.)) %>% 
  pivot_longer(2:5) %>% 
  mutate(Partition = factor(name,
                            levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
                            labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>%
  ggplot()+
  geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
  ylab(expression(R^2))+
  mytheme4()+
  ylim(0,1)+
  scale_fill_manual(values = c("#d2d2d2","#ffae49","#44b7c2","#024b7a"))+
  theme(legend.title = element_blank())

partition_dry <- df_kmaxww_a %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  mutate(LWP_q50 = quantile(LWP, 0.5,na.rm = TRUE)) %>%
  filter(!is.na(A),LWP<=LWP_q50) %>% 
  ungroup() %>% 
  rename(`Stomatal model`=scheme) %>%
  group_by(`Stomatal model`) %>% 
  do(get_partition(.)) %>% 
  pivot_longer(2:5) %>% 
  mutate(Partition = factor(name,
                            levels = c("Residuals","species_r2","non_stomatal_r2","stomatal_r2" ),
                            labels = c("Residuals","Species", "Non-stomatal","Stomatal"))) %>%
  ggplot()+
  geom_col(aes(x = `Stomatal model`,y=value,fill=Partition))+
  ylab(expression(R^2~"DRY CONDITIONS"))+
  mytheme4()+
  ylim(0,1)+
  scale_fill_manual(values = c("#d2d2d2","#ffae49","#44b7c2","#024b7a"))+
  theme(legend.title = element_blank())

ggarrange(partition_full, partition_dry,
          align='hv', labels=c('a', 'b'),
          common.legend = T,ncol=1, nrow = 2)

ggsave("PLOTS/partition.png", width = 16, height = 22, units = "cm")




#### Species-scheme A r pearson's ####

library(gridExtra)
library(grid)
df_a <- df_a_fix %>% 
  group_by(scheme,acclimation,Species,source) %>% 
  filter(!is.na(A)) %>% 
  mutate(diff_a = A - a_pred) %>% 
  summarise(n_dist = n(),
            r = cor(A, a_pred, use = "pairwise.complete.obs"),
            bias = mean(diff_a,na.rm = TRUE)/mean(A,na.rm = TRUE),
            rmse = Metrics::rmse(A,a_pred),
            beta = lm(A~a_pred)$coefficients[2])

df_a %>% 
  ungroup() %>% 
  filter(acclimation == "Acclimated") %>% 
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
size = "r", fill = "r") +
  scale_fill_viridis_b(name = "r Pearson's") +
  guides(size = "none")+scale_y_discrete(limits=rev)
    
ggsave("PLOTS/species_scheme_A.png", width = 16, height = 22, units = "cm")



#### Estimated vs Actual ####
df %>% 
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

df %>% 
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

df %>% 
  ggplot()+
  geom_point(aes(A,a_pred),shape = 1)+
  geom_abline(intercept=0,slope=1,color="grey20")+
  # geom_smooth(aes(g_pred,a_pred,color=scheme, group=interaction(Species,scheme)),
  #             se = FALSE,method="lm",linetype = 1, size=0.5)+
  geom_smooth(aes(A,a_pred,color=scheme),se = FALSE,method="lm", formula = y~poly(x,2))+
  facet_wrap(~acclimation)+
  mytheme2()+
  scale_colour_manual(breaks = col_df$scheme, 
                      values = unique(as.character(col_df$col)))+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.background = element_rect(fill="transparent"))+
  ylab(expression(atop("A estimated ("*mu*"mol m"^-2~"s"^-1~")")))+
  xlab(expression(atop("A actual ("*mu*"mol m"^-2~"s"^-1*")")))+
  guides(colour = guide_legend(nrow = 1))+
  NULL
