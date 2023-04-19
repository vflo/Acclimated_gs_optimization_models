################## optimization code ###################################

rm(list=ls())

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(rphydro)
library(DEoptim)
source("stomatal_optimization_functions.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")

#LOAD DATA
dat <-  read.csv(file="DATA/drying_experiments_meta-analysis_extended.csv")
names = colnames(dat)
# [5] "Predawn.LWP..MPa."                                 
# [6] "A..umol.m.2.s.1."                                  
# [7] "gC..mol.m.2.s.1."
# [8] "vcmax_obs"
# [9] "jmax_obs"
# [10] "ca..ppm."                                          
# [11] "D..unitless...Pa...Pa.atm.."                       
# [12] "T..deg.C."                                                          
names[5:12] = c("LWP", "A", "gC", "vcmax_obs","jmax_obs","ca", "D", "T")
names[16] = "Iabs_growth"
names[15] = "Iabs_used"
colnames(dat) = names
dat <- dat %>% filter(!is.na(D),!is.na(LWP),!is.na(T),!is.na(ca),!is.na(Iabs_used))


dpsi_df = read.csv(file = "DATA/drying_experiments_dpsi_extended.csv")

# template = read.csv("DATA/fitted_params_template.csv")
# path_par <- "DATA/parameters/"
path_par <- "DATA/parameters_kmaxww_alpha/"
par_data <- list.files(path_par) %>% 
  purrr::map_df(function(x){
    readr::read_csv(paste0(path_par,x))
  })



################################################################################

fun_no_accl = function(data, dpsi_calib=T,  k=7,
                       vcmax25 = vcmax25,
                       jmax25 = jmax25,
                       stomatal_model = stomatal_model_now,
                       par_plant_acclimation = par_plant_acclimation,
                       par_cost_acclimation = par_cost_acclimation,
                       par_plant_now = par_plant,
                       par_cost_now = par_cost,
                       Species_now = species){
  
  data = data %>% 
    mutate(  patm = calc_patm(0,T),
             ca_pa = ca*1e-6 * patm,
             Ciest = ca_pa-(A*1e-6)/(gC/patm))
  dpsi_data = dpsi_df %>% filter(Species == Species_now)
  
  lwp_actual = data$LWP
  dat1 = tibble(var = lwp_actual, jmax25_a=jmax25, vcmax25_a=vcmax25) %>% 
    cbind(data %>% select(t=T,Iabs_used, D,ca)) %>%
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           p = purrr::pmap(list(var, jmax25_a, vcmax25_a,t,Iabs_used,D,ca),
                           ~model_numerical_instantaneous(tc = ..4,
                                                          ppfd = ..5,
                                                          vpd = ..6*101325,
                                                          co2 = ..7, elv = 0,
                                                          fapar = .99, kphio = 0.087,
                                                          psi_soil = ..1, rdark = 0.015,
                                                          par_plant=par_plant_now,
                                                          par_cost = par_cost_now,
                                                          jmax25 = ..2, vcmax25 = ..3,
                                                          stomatal_model = stomatal_model)) ) %>%
           # p = pmap(list(var, jmax25_a, vcmax25_a), 
           #          ~model_numerical_instantaneous(tc = mean(data$T,na.rm = TRUE), 
           #                                         ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
           #                                         vpd = mean(data$D*101325,na.rm = TRUE), 
           #                                         co2 = mean(data$ca,na.rm = TRUE), 
           #                                         elv = 0, fapar = .99, kphio = 0.087, 
           #                                         psi_soil = ..1, rdark = 0.015, 
           #                                         par_plant=par_plant_now, 
           #                                         par_cost = par_cost_now, 
           #                                         jmax25 = ..2, vcmax25 = ..3,
           #                                         stomatal_model = stomatal_model))) %>%  
    unnest_wider(p)
  
  ndays = mean(data$Drydown.days)
  psi_crit = par_plant_now$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant_now$b)
  if(min(data$LWP,na.rm = TRUE)<psi_crit){
    psi_min = psi_crit
  }else{
    psi_min = min(data$LWP,na.rm = TRUE) #-6
  }
  psi_max = 0
  
  lwp = seq(psi_min,0, length.out=20)
  dat2 = tibble(var = lwp, jmax25_a=jmax25, vcmax25_a=vcmax25) %>%
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           p = purrr::pmap(list(var, jmax25_a, vcmax25_a), 
                           ~model_numerical_instantaneous(tc = mean(data$T,na.rm = TRUE), 
                                                          ppfd = mean(data$Iabs_used,na.rm = TRUE), 
                                                          vpd = mean(data$D*101325,na.rm = TRUE), 
                                                          co2 = mean(data$ca,na.rm = TRUE), elv = 0, 
                                                          fapar = .99, kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.015, 
                                                          par_plant = par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax25 = ..2, vcmax25 = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
  
  dat2 <- dat2 %>% filter(gs>=1e-40)
  gx = log(dat2$gs)
  gy = dat2$var
  fpsi = splinefun(x = gx, y=gy, method = "natural")
  gs0 = dat2$gs[which(dat2$var==0)]
  psi88S = fpsi(log(gs0*0.12))
  dpx = dat2$var
  dpy = dat2$dpsi
  f1 = splinefun(dpy~dpx)
  dp88S = f1(psi88S)
  psiL88S = psi88S-dp88S
  
  # A, G, CHI
  
  data_f <- data
  
  a_pred = dat1$a
  g_pred = dat1$gs
  c_pred = dat1$chi
  jmax25_pred = dat1$jmax25
  vcmax25_pred = dat1$vcmax25
  
  #DPSI
  d_spl = splinefun(lwp, y=dat2$dpsi)
  d_pred = d_spl(dpsi_data$SWP)
  d_pred2 = d_spl(data_f$LWP)
  
  data_f <- data_f %>% 
    cbind(a_pred = a_pred,
          g_pred = g_pred,
          c_pred = c_pred,
          d_pred = d_pred2,
          jmax25_pred = jmax25_pred,
          vcmax25_pred = vcmax25_pred,
          psi88S = psi88S,
          psiL88S = psiL88S,
          dp88S = dp88S,
          gs0 = gs0) 
  
  dpsi_data_f <- dpsi_data %>% 
    dplyr::select(Species,DAY,MDWP,SWP,Dpsi) %>%
    dplyr::rename(LWP = SWP) %>% 
    cbind(d_pred = d_pred) %>% 
    dplyr::bind_cols(data_f %>% 
                       dplyr::select(Source,Drydown.days) %>%
                       dplyr::summarise(Source = unique(Source),
                                        Drydown.days = unique(Drydown.days))) 
  
  # JOIN EVERYTHING
  res <- full_join(data_f, dpsi_data_f)%>% dplyr::distinct()
  
  return(res)
  
}


################################################################################

fun_accl = function(data, dpsi_calib=T, inst=F, k=7, 
                    stomatal_model = stomatal_model_now,
                    par_plant_now = par_plant,
                    par_cost_now = par_cost,
                    Species_now = species){
  
  data = data %>% 
    mutate(  patm = calc_patm(0,T),
             ca_pa = ca*1e-6 * patm,
             Ciest = ca_pa-(A*1e-6)/(gC/patm))
  dpsi_data = dpsi_df %>% filter(Species == Species_now)
  
  
  psi_crit = par_plant_now$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant_now$b)
  if(min(data$LWP,na.rm = TRUE)<psi_crit){
    psi_min = psi_crit
  }else{
    psi_min = min(data$LWP,na.rm = TRUE) #-6
  }
  psi_max = 0
  
  lwp = seq(psi_min,0, length.out=20)
  lwp_day = function(day_num){psi_max + day_num/ndays * (psi_min-psi_max)}
  ndays = mean(data$Drydown.days)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  actual_day = ndays * (data$LWP-psi_max)/(psi_min-psi_max)
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
  spl = splinefun(x = max(day):0, y=lwp_week)
  
  #ACCLIMATED
  dat_acc = tibble(var = spl(actual_day)) %>% 
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           pmod = map(var, ~model_numerical(tc = mean(data$T,na.rm = TRUE), ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
                                            vpd = mean(data$D*101325,na.rm = TRUE), co2 = mean(data$ca,na.rm = TRUE), 
                                            elv = 0, fapar = .99, kphio = 0.087, 
                                            psi_soil = ., rdark = 0.015, par_plant=par_plant_now, 
                                            par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  lwp_actual = data$LWP
  
  #INSTANTANEOUS
  dat1 = tibble(var = lwp_actual, jmax25_a=dat_acc$jmax25, vcmax25_a=dat_acc$vcmax25) %>% 
    cbind(data %>% select(t=T,Iabs_used, D,ca)) %>%
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           p = purrr::pmap(list(var, jmax25_a, vcmax25_a,t,Iabs_used,D,ca),
                           ~model_numerical_instantaneous(tc = ..4,
                                                          ppfd = ..5,
                                                          vpd = ..6*101325,
                                                          co2 = ..7, elv = 0,
                                                          fapar = .99, kphio = 0.087,
                                                          psi_soil = ..1, rdark = 0.015,
                                                          par_plant=par_plant_now,
                                                          par_cost = par_cost_now,
                                                          jmax25 = ..2, vcmax25 = ..3,
                                                          stomatal_model = stomatal_model)) ) %>%
           # p = pmap(list(var, jmax25_a, vcmax25_a), 
           #          ~model_numerical_instantaneous(tc = mean(data$T,na.rm = TRUE), 
           #                                         ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
           #                                         vpd = mean(data$D*101325,na.rm = TRUE), 
           #                                         co2 = mean(data$ca,na.rm = TRUE), 
           #                                         elv = 0, fapar = .99, kphio = 0.087, 
           #                                         psi_soil = ..1, rdark = 0.015, 
           #                                         par_plant=par_plant_now, 
           #                                         par_cost = par_cost_now, 
           #                                         jmax25 = ..2, vcmax25 = ..3,
           #                                         stomatal_model = stomatal_model))) %>% 
    unnest_wider(p)
  
  #Calculate all the dry-down to estimate psi88S and DPSI
  dat_acc = tibble(var = spl(day)) %>% 
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           pmod = map(var, ~model_numerical(tc = mean(data$T,na.rm = TRUE), 
                                            ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
                                            vpd = mean(data$D*101325,na.rm = TRUE), 
                                            co2 = mean(data$ca,na.rm = TRUE), 
                                            elv = 0, fapar = .99, kphio = 0.087, 
                                            psi_soil = ., rdark = 0.015, par_plant=par_plant_now, 
                                            par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  dat2 = tibble(var = lwp, jmax25_a=dat_acc$jmax25, vcmax25_a=dat_acc$vcmax25) %>%
    mutate(p = purrr::pmap(list(var, jmax25_a, vcmax25_a), 
                           ~model_numerical_instantaneous(tc = mean(data$T,na.rm = TRUE), 
                                                          ppfd = mean(data$Iabs_used,na.rm = TRUE), 
                                                          vpd = mean(data$D*101325,na.rm = TRUE), 
                                                          co2 = mean(data$ca,na.rm = TRUE), elv = 0, 
                                                          fapar = .99, kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.015, 
                                                          par_plant=par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax25 = ..2, vcmax25 = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
  
  gx = log(dat2$gs)
  gy = dat2$var
  fpsi = splinefun(x = gx, y=gy, method = "natural")
  gs0 = dat2$gs[which(dat2$var==0)]
  psi88S = fpsi(log(gs0*0.12))
  dpx = dat2$var
  dpy = dat2$dpsi
  f1 = splinefun(dpy~dpx)
  dp88S = f1(psi88S)
  psiL88S = psi88S-dp88S
  
  # A, G, CHI
  data_f <- data
  
  a_pred = dat1$a
  g_pred = dat1$gs
  c_pred = dat1$chi
  jmax25_pred = dat1$jmax25
  vcmax25_pred = dat1$vcmax25
  
  #DPSI
  d_spl = splinefun(lwp, y=dat2$dpsi)
  d_pred = d_spl(dpsi_data$SWP)
  d_pred2 = d_spl(data_f$LWP)
  
  data_f <- data_f %>% 
    cbind(a_pred = a_pred,
          g_pred = g_pred,
          c_pred = c_pred,
          d_pred = d_pred2,
          jmax25_pred = jmax25_pred,
          vcmax25_pred = vcmax25_pred,
          psi88S = psi88S,
          psiL88S = psiL88S,
          dp88S = dp88S,
          gs0 = gs0) 
  
  dpsi_data_f <- dpsi_data %>% 
    dplyr::select(Species,DAY,MDWP,SWP,Dpsi) %>%
    dplyr::rename(LWP = SWP) %>% 
    cbind(d_pred = d_pred) %>% 
    dplyr::bind_cols(data_f %>% 
                       dplyr::select(Source,Drydown.days) %>%
                       dplyr::summarise(Source = unique(Source),
                                        Drydown.days = unique(Drydown.days)))
  
  # JOIN EVERYTHING
  res <- full_join(data_f, dpsi_data_f) %>% dplyr::distinct()
  
  return(res)
}


##### SIMULATION #####
get_simulations <- function(x){
  species = x$Species %>% unique()
  stomatal_model_now = x$scheme %>% unique()
  print(species)
  print(stomatal_model_now)
  dpsi_calib = x$dpsi %>% unique()
  inst = x$inst %>% unique()
  data1=filter(dat, Species==species, Source == x$source %>% unique())
  
  
  ##### SIMULATION WITHOUT ACCLIMATION ####
  par_plant = list(
    conductivity= x[2,"K.scale"][[1]]*1e-16,
    psi50 = x[2,"P50"][[1]], 
    b= x[2,"b"][[1]]
  )
  par_cost = list(
    alpha  = x[2,"alpha"][[1]], 
    gamma = x[2,"gamma"][[1]]
  )
  
  data_ww <- data1 %>% 
    filter(!is.na(gC)) %>% 
    mutate(LWP_q90 = quantile(LWP, 0.8, na.rm = TRUE),
           patm = calc_patm(0,T),
           ca_pa = ca*1e-6 * patm,
           ci = ca_pa-(A*1e-6)/(gC/patm)) %>% 
    filter(LWP >= LWP_q90) 
  
  vcmax <- calc_vcmax_no_acclimated_ww(A = data_ww$A,
                                       ci = data_ww$ci,
                                       tc = data_ww$T,
                                       patm = calc_patm(0,data_ww$T),
                                       rdark = 0.0150
  )
  
  vcmax25 = calc_vcmax_arrhenius(vcmaxT1 = vcmax, T1 = (273.15+data_ww$T),
                                 T2 = 298.15)
  vcmax = mean(vcmax)
  vcmax25 = mean(vcmax25)
  jmax <- calc_jmax_no_acclimated_ww(A = data_ww$A,
                                     vcmax = vcmax,
                                     ci = data_ww$ci,
                                     I = data_ww$Iabs_growth,
                                     tc = data_ww$T,
                                     patm = calc_patm(0,data_ww$T),
                                     kphio = 0.087
  )
  jmax25 = calc_jmax_arrhenius(jmaxT1 = jmax, T1 = (273.15+data_ww$T),
                               T2 = 298.15)
  jmax25 = mean(jmax25)
  if(is.na(jmax25)){jmax25 = 1.6*vcmax25} 
  
  no_accl <- fun_no_accl(data1,
                         dpsi_calib = dpsi_calib,
                         vcmax25 = vcmax25,
                         jmax25 = jmax25,
                         stomatal_model = stomatal_model_now,
                         par_plant_acclimation = par_plant_acclimation,
                         par_cost_acclimation = par_cost_acclimation,
                         par_plant_now = par_plant,
                         par_cost_now = par_cost,
                         Species_now = species)
  
  no_accl <- no_accl %>% 
    left_join(x %>% filter(acclimation == FALSE)) %>% 
    cbind(calibration_type = 'no_alpha',
          jmaxww25 = jmax25,
          vcmaxww25 = vcmax)
  
  
  ##### SIMULATION WITH ACCLIMATION #####
  par_plant_acclimation = as.list(tibble(
    conductivity = x[1,"K.scale"][[1]]*1e-16 ,
    psi50 = x[1,"P50"][[1]], 
    b= x[1,"b"][[1]]
    ))
  par_cost_acclimation = as.list(data.frame(
    alpha  = x[1,"alpha"][[1]], 
    gamma = x[1,"gamma"][[1]]
    ))

  accl <- fun_accl(data1,
           Species_now=species,
           dpsi_calib = dpsi_calib, 
           inst = inst,
           par_plant_now = par_plant_acclimation,
           par_cost_now = par_cost_acclimation,
           stomatal_model = stomatal_model_now)
  
  accl <- accl %>% 
    left_join(x %>% filter(acclimation == TRUE)) %>% 
    cbind(calibration_type = 'alpha',
          jmaxww25 = jmax25,
          vcmaxww25 = vcmax)
  
  ##### SIMULATION WITH ACCLIMATION ALPHA FIX #####
  par_plant_acclimation = as.list(tibble(
    conductivity = x[1,"K.scale"][[1]]*1e-16 ,
    psi50 = x[1,"P50"][[1]], 
    b= x[1,"b"][[1]]
  ))
  par_cost_acclimation = as.list(data.frame(
    alpha  = 0.0930, 
    gamma = x[1,"gamma"][[1]]
  ))
  
  accl_alpha_fix <- fun_accl(data1,
                   Species_now=species,
                   dpsi_calib = dpsi_calib, 
                   inst = inst,
                   par_plant_now = par_plant_acclimation,
                   par_cost_now = par_cost_acclimation,
                   stomatal_model = stomatal_model_now)
  
  accl_alpha_fix <- accl_alpha_fix %>% 
    left_join(x %>% filter(acclimation == TRUE)) %>% 
    cbind(calibration_type = 'alpha_fix',
          jmaxww25 = jmax25,
          vcmaxww25 = vcmax25) %>% 
    mutate(alpha = 0.0930)


  
  ##### RETURN #####
  
  return(bind_rows(accl,accl_alpha_fix,no_accl))
  
}


df <- par_data %>% 
  # filter(!scheme %in% c("CMAX")) %>%
  # rbind(par_data_extra) %>%
  group_split(Species,scheme,dpsi,source) %>%
  purrr::map(get_simulations) %>%
  bind_rows()
# 
save(df, file = "DATA/simulations_kmaxww_alpha.RData")



