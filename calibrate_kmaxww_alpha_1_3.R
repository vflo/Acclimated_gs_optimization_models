################## optimization code ###################################

rm(list=ls())

# library(dplyr)
library(purrr)
# library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(rphydro)
library(DEoptim)
library(tidyverse)
library(grid)
# library(rvalues)
# library(furrr)
# plan('multisession', workers = 6)
# options('future.global.maxsize'=2*1024*1024^2)
source("stomatal_optimization_functions.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")
select = dplyr::select

#Convenience functions
opt_curve_param_vc <- function(par,p88,p50,p12){
  p <- par[1]
  b <- par[2]
  # if(p12>=-1){
  psi <- c(p88,p50)
  res <- (1/2)^((psi/p)^b)
  vulne_curve <- c(0.12,0.5)
  rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
  #   
  # }else{
  #   psi <- c(p88,p50,p12)
  #   res <- (1/2)^((psi/p)^b)
  #   vulne_curve <- c(0.12,0.5,0.88)
  #   rmse <- sqrt(sum((vulne_curve-res)^2)/length(vulne_curve))
  # }
  
  if(b<1){rmse <- 1e6}
  return(rmse)
}

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

dpsi_df <-  read.csv(file = "DATA/drying_experiments_dpsi_extended.csv")

traits <- read.csv(file="DATA/imputation_df.csv") %>% rename(species = 'Species', genus = 'Genus',Species = "Binomial")
traits[which(traits$Species == "Olea europaea var. meski"), "Species"] <- "Olea europaea var. Meski"
traits[which(traits$Species == "Olea europaea var. chemlali"), "Species"] <- "Olea europaea var. Chemlali"
traits[which(traits$Species == "Beta maritima subsp.marcosii"), "Species"] <- "Beta maritima subsp. marcosii"
traits[which(traits$Species == "Beta maritima subsp. maritima"), "Species"] <- "Beta maritima subsp. maritima"
# traits <- traits %>% separate(taxon,c("genus","species", "subsp"),sep="_")

template <-  read.csv("DATA/fitted_params_template.csv")
template <- template %>% dplyr::select(-c(species,genus))%>% 
  left_join(traits, by=c('Species','ID')) %>%
  filter(!is.na(P50..MPa.)) %>%
  rowwise() %>%
  mutate(P50 = optimr::optimr(c(p=-2,b=2),
                              fn=opt_curve_param_vc,
                              p88=P88..MPa.,
                              p50=P50..MPa.,
                              p12=P12..MPa.)$par[1],
         b = optimr::optimr(c(p=-2,b=2),
                            fn=opt_curve_param_vc,
                            p88=P88..MPa.,
                            p50=P50..MPa.,
                            p12=P12..MPa.)$par[2]) %>% 
  ungroup() 



plot_all = function(df_w_vol, varname, species, 
                    data, dpsi_data=NULL,
                    stomatal_model, param,
                    param_cost,analytical=F){
  # df_w_vol = df_w_vol[complete.cases(df_w_vol),]
  # View(df_w_vol)
  
  gx = log(df_w_vol$gs+1e-20)
  gy = df_w_vol$var
  # gx1 = seq(max(min(gx), -20), max(gx), length.out=100)
  f = splinefun(x = gx, y=gy, method = "monoH.FC")
  gs0 = df_w_vol$gs[which(df_w_vol$var==0)]
  psi88S = f(log(gs0*0.12))
  psi50S = f(log(gs0*0.50))
  psi12S = f(log(gs0*0.88))
  cat("psi88S = ", psi88S, "\n")
  cat("psi50S = ", psi50S, "\n")
  cat("psi12S = ", psi50S, "\n")
  
  dpx = df_w_vol$var
  dpy = df_w_vol$dpsi
  f1 = splinefun(dpy~dpx)
  dpx1 = seq(min(dpx), max(dpx), length.out=100)
  dp88S = f1(psi88S)
  dp50S = f1(psi50S)
  dp12S = f1(psi12S)
  cat("psiL88S = ", psi88S-dp88S, "\n")
  cat("psiL50S = ", psi50S-dp50S, "\n")
  cat("psiL12S = ", psi12S-dp12S, "\n")
  
  cat(psi88S, "\t", psi50S, "\t", psi12S, "\t", psi88S-dp88S, "\t", psi50S-dp50S, psi12S-dp12S, "\n")
  
  subdata = data %>% filter(LWP > -6)
  
  p1 <- df_w_vol %>% 
    ggplot() +
    geom_line(aes(x = var, y = vcmax25), col="green3", linewidth=1) +
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Vcmax))+
    # geom_vline(xintercept = par_plant_std$psi50, col="orange") +
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p1 = p1 + geom_line(aes(x = var, y = vcmax ), col="grey", linewidth=1) 
  
  p2 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = dpsi), col="blue", linewidth=1)+
    # geom_vline(xintercept = psi88S, col="grey", linewidth=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (!is.null(dpsi_data)){
    cat("Adding points...\n")
    p2 <- p2 + 
      geom_point(data=dpsi_data, aes(x=SWP, y=Dpsi))
  }
  
  p3 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = gs), col="cyan2", linewidth=1)+
    geom_point(data=subdata, aes(x=LWP, y=gC))+
    # geom_vline(xintercept = psi88S, col="grey", linewidth=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p3 = p3 + geom_line(aes(x = var, y = gs), col="grey", linewidth=1)
  
  p4 <- df_w_vol %>%
    # mutate(chi = ci/out_hydraulics_ca) %>% 
    ggplot() +
    geom_line(aes(x = var, y = chi), col="magenta", linewidth=1)+
    geom_point(data=subdata, aes(x=LWP, y=1-A/gC/ca))+
    # geom_vline(xintercept = psi88S, col="grey", linewidth=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p4 = p4 + geom_line(aes(x = var, y = chi), col="grey", linewidth=1)
  
  p5 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = jmax25), col="goldenrod1", linewidth=1) +
    # geom_vline(xintercept = par_plant_std$psi50, col="orange")+
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  
  
  p6 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = a), col="green4", linewidth=1) +
    geom_point(data=subdata, aes(x=LWP, y=A))+
    # geom_vline(xintercept = psi88S, col="grey", linewidth=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p6 = p6 + geom_line(aes(x = var, y = gpp), col="grey",linewidth=1) 
  
  if(stomatal_model %in% par_scheme_gamma){
    grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2,
                 top = textGrob(paste("Model:",stomatal_model,
                                      "/ Species:",species,
                                      "/ gamma:",param_cost$gamma),
                                gp=gpar(fontsize=12,font=3)))
  }else{
    grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2,
                 top = textGrob(paste("Model:",stomatal_model,
                                      "/ Species:",species,
                                      "/ K:",param$conductivity*1e16),
                                gp=gpar(fontsize=12,font=3)))
  }
  
}

################################################################################

error_fun_no_accl = function(x, data, data_template, plot=F, 
                             dpsi_data = dpsi_data,
                             stomatal_model = stomatal_model_now,
                             vcmax25 = vcmax25,
                             jmax25 = jmax25, 
                             Species_now = species,
                             K_PROFITMAX = K_PROFITMAX_no_acclimate){
  
  data = data %>% 
    mutate(  patm = calc_patm(0),
             ca_pa = ca*1e-6 * patm,
             Ciest = ca_pa-(A*1e-6)/(gC/patm))
  
  if(stomatal_model %in% par_scheme_gamma){
    
    par_plant_now = list(
      conductivity = K_PROFITMAX$K_PROFITMAX*1e-16,
      psi50 = data_template$P50 %>% unique()/3,
      b = 1
    )
    
    par_cost_now = list(
      alpha = 0.1,
      gamma = x[1]
    )
    
  }else{
    
    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = data_template$P50%>% unique()/3,
      b = 1
    )
    
    par_cost_now = list(
      alpha  = 0.1
    )
    
  }
  
  lwp = data$LWP
  dat1 = tibble(var = lwp, jmax25_a=jmax25, vcmax25_a=vcmax25) %>% 
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
                                                          stomatal_model = stomatal_model)) 
    ) %>% 
    unnest_wider(p)
  
  if(plot==T) dat1 %>% plot_all(varname = "psi_soil", 
                                species=Species_now, 
                                data = data, 
                                dpsi_data=dpsi_data,
                                param = par_plant_now,
                                param_cost = par_cost_now,
                                stomatal_model=stomatal_model)
  
  # dat2 <- dat1 %>% filter(gs>=1e-40)
  # gx = log(dat2$gs)
  # gy = dat2$var
  # fpsi = splinefun(x = gx, y=gy, method = "natural")
  # gs0 = dat2$gs[which(dat2$var==0)]
  # psi88S = fpsi(log(gs0*0.12))
  # dpx = dat2$var
  # dpy = dat2$dpsi
  # f1 = splinefun(dpy~dpx)
  # dp88S = f1(psi88S)
  # psiL88S = psi88S-dp88S
  
  data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
  y2 = mean((dat1$gs - data_f$gC)^2,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2
  y1 = mean((dat1$a - data_f$A)^2,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
  y4 = mean((dat1$chi - (data_f$Ciest/data_f$ca_pa))^2,na.rm  = TRUE)/mean((data_f$Ciest/data_f$ca_pa),na.rm  = TRUE)^2
  
  # if (!is.null(dpsi_data)){
  #   d_spl = splinefun(lwp, y=dat1$dpsi)
  #   dpsi_data_f <- dpsi_data #%>% filter(SWP >= psi88S) #use only values over Psi88S
  #   y3 = mean((d_spl(dpsi_data_f$SWP) - dpsi_data_f$Dpsi)^2,na.rm  = TRUE)/mean(dpsi_data_f$Dpsi,na.rm  = TRUE)^2 #*40
  #   cat("d_spl:", d_spl(dpsi_data_f$SWP), "\n")
  # }else{
  #   y3=0
  # }
  
  y=y2+y1#+y4
  
  cat(x, "|", y2, " / ", y1, " / ", y4, " / ",y, "\n")
  cat(x, "|", y, "\n")
  
  y
  
}

################################################################################

error_fun_kmax_alpha = function(x, data_all, data_ww, data_template,  plot=F, 
                                k=7, stomatal_model = stomatal_model_now, 
                                dpsi_data = dpsi_data, Species_now = species,
                                res_ww =res_ww, K_PROFITMAX = K_PROFITMAX_acclimate){
  # parameter_max <- c(0.2)
  # 
  # if(x[1]<=0 | x[1]>parameter_max[1] ){over <- TRUE}else{over <- FALSE} #set boundaries
  # 
  # if(over){
  #   print(paste("Over-limits"))
  #   return(1e6)
  # }else{
  
  data = data_ww %>% 
    mutate(  patm = calc_patm(0),
             ca_pa = ca*1e-6 * patm,
             Ciest = ca_pa-(A*1e-6)/(gC/patm))
  
  if(stomatal_model %in% par_scheme_gamma){
    par_plant_now = list(
      conductivity = res_ww$K.scale*1e-16,
      psi50 = data_template$P50%>% unique()/3,
      b = 1
    )
    par_cost_now = list(
      alpha = x[1],
      gamma = res_ww$gamma
    )
  }else{
    par_plant_now = list(
      conductivity = res_ww$K.scale*1e-16,
      psi50 = data_template$P50 %>% unique()/3,
      b = 1
    )
    par_cost_now = list(
      alpha  = x[1]
    )
  }
  
  ndays = mean(data_all$Drydown.days)
  psi_crit = data_template$P50 * (log(1000)/log(2)) ^ ( 1/data_template$b)
  if(min(data$LWP,na.rm = TRUE)<psi_crit){
    psi_min = psi_crit
  }else{
    psi_min = min(data_all$LWP,na.rm = TRUE) #-6
  }
  psi_max = 0 #max(data$LWP)
  # cat(ndays,"\n")
  
  lwp = seq(psi_min,0, length.out=20)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  actual_day = ndays * (data$LWP-psi_max)/(psi_min-psi_max)
  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }
  
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
  spl = splinefun(x = max(day):0, y=lwp_week)
  dat_acc = try(
    tibble(var = spl(actual_day)) %>% 
      mutate(var = case_when(var>0~0,
                             TRUE~var),
             pmod = map(var, ~model_numerical(tc = mean(data$T,na.rm = TRUE), ppfd = mean(data$Iabs_growth,na.rm = TRUE), 
                                              vpd = mean(data$D*101325,na.rm = TRUE), co2 = mean(data$ca,na.rm = TRUE), 
                                              elv = 0, fapar = .99, kphio = 0.087, 
                                              psi_soil = ., rdark = 0.015, par_plant=par_plant_now, 
                                              par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
      unnest_wider(pmod),
    silent = TRUE)
  
  if(any(class(dat_acc) %in% "try-error")){
    print(paste("error try"))
    return(1e6)
  }else{
    
    lwp = data$LWP
    dat1 = try(
      tibble(var = lwp, jmax25_a=dat_acc$jmax25, vcmax25_a=dat_acc$vcmax25) %>% 
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
        unnest_wider(p),
      silent = TRUE)
    
    if(any(class(dat1) %in% "try-error")){
      print(paste("error try"))
      return(1e6)
    }else{
      if(plot==T) dat1 %>% plot_all(varname = "psi_soil", 
                                    species=Species_now, 
                                    data = data, 
                                    dpsi_data=dpsi_data,
                                    param = par_plant_now,
                                    param_cost = par_cost_now,
                                    stomatal_model=stomatal_model)
      
      
      data_f <- data #%>% filter(LWP >= psi88S) #use only values over Psi88S
      y2 = mean((dat1$gs - data_f$gC)^2,na.rm  = TRUE)/mean(data_f$gC,na.rm  = TRUE)^2
      y1 = mean((dat1$a - data_f$A)^2,na.rm  = TRUE)/mean(data_f$A,na.rm  = TRUE)^2
      y4 = mean((dat1$chi - (data_f$Ciest/data_f$ca_pa))^2,na.rm  = TRUE)/mean((data_f$Ciest/data_f$ca_pa),na.rm  = TRUE)^2
      
      y=y2+y1#+y4
      
      cat(x, "|", y2, " / ", y1, " / ", y4, " / ",y, "\n")
      cat(x, "|", y, "\n")
      
      y
    }
  }
}


par_scheme_gamma <- list("PHYDRO","CGAIN", "CMAX", "SOX2")
par_scheme_no_gamma <- list("PROFITMAX2","SOX","PROFITMAX","PMAX3")

##### PARAMETERIZATION #####
get_parameters_kmaxww_alpha <- function(x){
    species = x$Species %>% unique()
    stomatal_model_now = x$scheme %>% unique()
    data_template_now = x
    
    data1 = filter(dat, Species==species, Source == unique(x$source))
    
    if(species %in% dpsi_df$Species){
      dpsi_data = dpsi_df %>% filter(Species == species)
    }else{
      dpsi_data = NULL
    }
    
    if(!is.null(K_PROFITMAX)){
    K_PROFITMAX_no_acclimate = K_PROFITMAX %>% 
      filter(Species == species,
             acclimation == FALSE,
             source == x$source) %>% 
      dplyr::select(K_PROFITMAX) %>% 
      unique()
    # K_PROFITMAX_acclimate = K_PROFITMAX %>% 
    #   filter(Species == species,
    #          acclimation == TRUE,
    #          source == x$source) %>% 
    #   dplyr::select(K_PROFITMAX) %>% 
    #   unique()
    }
    
    ##### PARAMETERIZATION WITHOUT ACCLIMATION WW #####
    data_ww <- data1 %>% 
        filter(!is.na(gC)) %>% 
        mutate(LWP_q90 = quantile(LWP, 0.8, na.rm = TRUE),
               patm = calc_patm(0),
               ca_pa = ca*1e-6 * patm,
               ci = ca_pa-(A*1e-6)/(gC/patm)) %>% 
        filter(LWP >= LWP_q90) 

    vcmax <- calc_vcmax_no_acclimated_ww(A = data_ww$A,
                                         ci = data_ww$ci,
                                         tc = data_ww$T,
                                         patm = calc_patm(0),
                                         rdark = 0.015
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
                                       patm = calc_patm(0),
                                       kphio = 0.087
    )
    jmax25 = calc_jmax_arrhenius(jmaxT1 = jmax, T1 = (273.15+data_ww$T),
                                 T2 = 298.15)
    jmax25 = mean(jmax25)
    if(is.na(jmax25)){jmax25 = 1.6*vcmax25} 
    
    ##### PARAMETERIZATION WITHOUT ACCLIMATION #####
    print(stomatal_model_now)
    print(species)
    parameter_max <- 100
    if(stomatal_model_now %in% c("CGAIN")&
       !species %in%c("Broussonetia papyrifera")){parameter_max <- 60}
    if(stomatal_model_now %in% c("PROFITMAX")&
       species %in%c("Broussonetia papyrifera","Malva subovata")){parameter_max <- 300}
    # if(stomatal_model_now %in% c("CMAX") &
    #    species %in% c("Pyracantha fortuneana")){parameter_max <- 6}
    if(stomatal_model_now %in% c("CMAX")){parameter_max <- 50}
    # if(!stomatal_model_now %in% c("CMAX")&
    #    species %in%c("Broussonetia papyrifera")){parameter_max <- 200}
    if(stomatal_model_now %in% c("PHYDRO")&
       species %in%c("Broussonetia papyrifera",
                     "Pteroceltis tatarinowii",
                     "Helianthus annuus",
                     "Malva subovata",
                     "Pseudotsuga menziesii")){parameter_max <- 300}
    if(stomatal_model_now %in% c("SOX")&
       species %in%c("Broussonetia papyrifera",
                     "Pteroceltis tatarinowii",
                     "Malva subovata")){parameter_max <- 300}
    if(stomatal_model_now %in% c("PROFITMAX2")&
       species %in%c("Broussonetia papyrifera",
                     "Helianthus annuus",
                     "Malva subovata")){parameter_max <- 200}
    if(stomatal_model_now %in% c("SOX2")){parameter_max <- 5}
    # if(stomatal_model_now %in% c("CGAIN2")){parameter_max <- 3}
    
    optimise(error_fun_no_accl,
             interval = c(0,parameter_max),
             data=data_ww,
             data_template = data_template_now,
             dpsi_data = dpsi_data,
             stomatal_model = stomatal_model_now,
             vcmax25 = vcmax25,
             jmax25 = jmax25, 
             Species_now = species,
             K_PROFITMAX = K_PROFITMAX_no_acclimate
    ) -> opt_no_accl
    x_no_accl <- opt_no_accl$minimum

    if(stomatal_model_now %in% par_scheme_gamma){
      res_ww <- tibble(K.scale=K_PROFITMAX_no_acclimate$K_PROFITMAX,
                       gamma=x_no_accl[1])
    }else{
      res_ww <- tibble(K.scale=x_no_accl[1],
                       gamma=NA)
    }

    error_fun_no_accl(x_no_accl, data1, 
                      data_template = data_template_now, plot=T,
                      dpsi_data = dpsi_data,
                      stomatal_model = stomatal_model_now,
                      vcmax25 = vcmax25,
                      jmax25 = jmax25,
                      Species_now = species,
                      K_PROFITMAX = K_PROFITMAX_no_acclimate)
    
    if(stomatal_model_now %in% par_scheme_gamma){
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=K_PROFITMAX_no_acclimate$K_PROFITMAX,
                            alpha=NA,
                            gamma=x_no_accl[1])
    }else{
      res_no_accl <- tibble(x,
                            acclimation = FALSE,
                            K.scale=x_no_accl[1],
                            alpha=NA,
                            gamma=NA)
    }
    
    ##### PARAMETERIZATION WITH ACCLIMATION #####
    print(stomatal_model_now)
    print(species)
    # parameter_ini <- c(0.1) #hydraulic parameter and alpha
    optimise(error_fun_kmax_alpha,
             interval = c(0.000001,0.2),
             data_all = data1,
             data_ww = data_ww,
             data_template = data_template_now,
             dpsi_data = dpsi_data,
             stomatal_model = stomatal_model_now,
             Species_now = species,
             K_PROFITMAX = K_PROFITMAX_no_acclimate,
             res_ww = res_ww#,
                     # control = list(maxit = 500, maximize = TRUE,
                     #                REPORT=0, trace=0, reltol=1e-4)
             ) -> opt_accl
      # x_accl <- opt_accl$par
      x_accl <- opt_accl$minimum
    
    error_fun_kmax_alpha(x_accl, data_all = data1,
                         data1, data_template = data_template_now,
                         dpsi_data = dpsi_data, plot=T, 
                         stomatal_model = stomatal_model_now, 
                         Species_now = species, res_ww =res_ww,
                         K_PROFITMAX = K_PROFITMAX_no_acclimate)
    
    if(stomatal_model_now %in% par_scheme_gamma){
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=res_ww$K.scale,
                         alpha=x_accl[1],
                         gamma=res_ww$gamma)
    }else{
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=res_ww$K.scale,
                         alpha=x_accl[1],
                         gamma=res_ww$gamma)
    }
    
  df <- bind_rows(res_accl,res_no_accl)
  
  readr::write_csv(df,file=paste0("DATA/parameters_kmaxww_alpha_1_3_no_chi/",stomatal_model_now,"_",species,"_",x$source,".csv"))
  
  return(bind_rows(res_accl,res_no_accl))

}

##### COMPUTE PARAMETERS #####
#First compute PROFITMAX model to obtain Kmax for CMAX. CGAIN, WUE and PHYDRO models
K_PROFITMAX <- NULL
template %>% filter(scheme == "PROFITMAX"#,Species %in% c("Diplotaxis ibicensis")
                    ) %>%
  group_split(scheme, dpsi, Species,source) %>%
  purrr::map_df(get_parameters_kmaxww_alpha)->res
#
save(res,file = "DATA/Kmax_PROFITMAX_kmaxww_alpha_1_3_no_chi.RData")
# 
load(file = "DATA/Kmax_PROFITMAX_kmaxww_alpha_1_3_no_chi.RData")

K_PROFITMAX <- res %>% 
  dplyr::select(Species,K_PROFITMAX = K.scale,dpsi,acclimation,source) %>% 
  group_by(Species,dpsi, acclimation,source) %>% 
  summarise_all(unique)

#Compute the other models
template %>% 
  filter(!scheme %in% c("PROFITMAX","CMAX")#,scheme %in% c("PMAX3")
         # Species %in% c("Pteroceltis tatarinowii")
         #   "Ficus tikoua"
           # "Malva subovata"
         #   # "Rosa cymosa",
           # "Broussonetia papyrifera"#,
           # "Helianthus annuus"
         #   # "Cinnamomum bodinieri",
         #   # "Platycarya longipes",
           # "Pteroceltis tatarinowii"
         #   # "Picea abies",
         #   # "Betula pendula",
         #   # "Pinus sylvestris",
         #   # "Populus tremula"
         # )
         ) %>%
  # filter(scheme %in% c("CGAIN")) %>%
  # filter(Species == "Quercus ilex", source == "Epron and Dreyer (1990)") %>%
  group_split(scheme, dpsi, Species,source) %>% 
  purrr::map_df(get_parameters_kmaxww_alpha)
