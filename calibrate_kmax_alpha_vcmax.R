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
# library(furrr)
# plan('multisession', workers = 6)
# options('future.global.maxsize'=2*1024*1024^2)
source("stomatal_optimization_functions.R")
source('hydraulic_functions.R')
source('photosynthetic_functions.R')
source("QUADP.R")

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

template <-  read.csv("DATA/fitted_params_template_vcmax.csv")
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
  

plot_all <- function(df_w_vol, varname, species, data, dpsi_data=NULL, analytical=F){
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
  f1 = splinefun(dpy~dpx, method = "monoH.FC")
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
    geom_point(data=data %>% 
                 mutate(vcmax25_obs = calc_vcmax_arrhenius(vcmax_obs,T+273.15)),
               aes(x=LWP, y=vcmax25_obs))+
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
    # geom_point(data=filter(dat, Species==species), aes(x=LWP, y=Jmax))+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if(any(!is.na(data$jmax_obs))){
    p5 = p5 + geom_point(data=data %>%
                         mutate(jmax25_obs = calc_jmax_arrhenius(jmax_obs,T+273.15)),
                        aes(x=LWP, y=jmax25_obs))
  }
  
  
  p6 <- df_w_vol %>%
    ggplot() +
    geom_line(aes(x = var, y = a), col="green4", linewidth=1) +
    geom_point(data=subdata, aes(x=LWP, y=A))+
    # geom_vline(xintercept = psi88S, col="grey", linewidth=0.8)+
    expand_limits(y=0)+
    xlab(varname)+
    theme_bw()
  if (analytical) p6 = p6 + geom_line(aes(x = var, y = gpp), col="grey",linewidth=1) 
  
  grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}


################################################################################

error_fun_kmax_alpha = function(x, 
                                data, 
                                data_template, 
                                plot=F, 
                                k=7, 
                                stomatal_model = stomatal_model_now, 
                                Species_now = species){

  if(x[1]<=0| x[2]<=0 |x[3]>-0.5| x[4]<1){over <- TRUE}else{over <- FALSE} #set boundaries
  
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
    
    data = data %>% 
      mutate(patm = calc_patm(0,T),
             ca_pa = ca*1e-6 * patm,
             Ciest = ca_pa-(A*1e-6)/(gC/patm))

    par_plant_now = list(
      conductivity = x[1]*1e-16,
      psi50 = x[3],#data_template$P50 %>% unique(),
      b = x[4]#data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha  = x[2]
    )

    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now)
    }
    else{
      dpsi_data=NULL
    }

  ndays = mean(data$Drydown.days)
  psi_crit = data_template$P50 * (log(1000)/log(2)) ^ ( 1/data_template$b)
  if(min(data$LWP,na.rm = TRUE)<psi_crit){
    psi_min = psi_crit
  }else{
    psi_min = min(data$LWP,na.rm = TRUE) #-6
  }
  psi_max = 0 #max(data$LWP)
  
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
                                                  psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
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
                                                                  psi_soil = ..1, rdark = 0.02, 
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
          if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, 
                                        data = data, dpsi_data=dpsi_data)

          y2 = mean((dat1$gs - data$gC)^2,na.rm  = TRUE)/mean(data$gC,na.rm  = TRUE)^2
          y1 = mean((dat1$a - data$A)^2,na.rm  = TRUE)/mean(data$A,na.rm  = TRUE)^2
          y4 = mean((dat1$chi - (data$Ciest/data$ca_pa))^2,na.rm  = TRUE)/
            mean((data$Ciest/data$ca_pa),na.rm  = TRUE)^2
          y3 = mean((calc_vcmax_arrhenius(dat1$vcmax25,298.15,dat1$t+273.15) - data$vcmax_obs)^2,na.rm  = TRUE)/
            mean(data$vcmax_obs,na.rm  = TRUE)^2
          if(all(is.na(data$jmax_obs))){
            y5=0
          }else{
            y5 = mean((calc_jmax_arrhenius(dat1$jmax25,298.15,dat1$t+273.15) - data$jmax_obs)^2,na.rm  = TRUE)/
              mean(data$jmax_obs,na.rm  = TRUE)^2
          }

          
          y=y2+y1+y3+y4+y5
            
          cat(x, "|", y, "\n")
            
          y
            
        }
      }
  }
}



error_fun_kmax_alpha_gamma = function(x, 
                                      data, 
                                      data_template, 
                                      plot=F, 
                                      k=7, 
                                      stomatal_model = stomatal_model_now, 
                                      Species_now = species){

  if(x[1]<=0| x[2]<=0| x[3]<=0 |x[4]>-0.5| x[5]<1
  ){over <- TRUE}else{over <- FALSE} #set boundaries
  
  if(over){
    print(paste("Over-limits"))
    return(1e6)
  }else{
    
    data = data %>% 
      mutate(  patm = calc_patm(0,T),
               ca_pa = ca*1e-6 * patm,
               Ciest = ca_pa-(A*1e-6)/(gC/patm))

      par_plant_now = list(
        conductivity =x[1]*1e-16,
        psi50 = x[4],#data_template$P50%>% unique(),
        b = x[5]#data_template$b%>% unique()
      )
      
      par_cost_now = list(
        alpha = x[3],
        gamma = x[2]
      )

    if (Species_now != ""){
      dpsi_data = dpsi_df %>% filter(Species == Species_now)
    }
    else{
      dpsi_data=NULL
    }
    
    ndays = mean(data$Drydown.days)
    psi_crit = data_template$P50 * (log(1000)/log(2)) ^ ( 1/data_template$b)
    if(min(data$LWP,na.rm = TRUE)<psi_crit){
      psi_min = psi_crit
    }else{
      psi_min = min(data$LWP,na.rm = TRUE) #-6
    }
    psi_max = 0 #max(data$LWP)

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
                                                psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
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
                                                                psi_soil = ..1, rdark = 0.02, 
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
        if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=species, data = data, dpsi_data=dpsi_data)

        y2 = mean((dat1$gs - data$gC)^2,na.rm  = TRUE)/mean(data$gC,na.rm  = TRUE)^2
        y1 = mean((dat1$a - data$A)^2,na.rm  = TRUE)/mean(data$A,na.rm  = TRUE)^2
        y4 = mean((dat1$chi - (data$Ciest/data$ca_pa))^2,na.rm  = TRUE)/
          mean((data$Ciest/data$ca_pa),na.rm  = TRUE)^2
        y3 = mean((calc_vcmax_arrhenius(dat1$vcmax25,298.15,dat1$t+273.15) - data$vcmax_obs)^2,na.rm  = TRUE)/
          mean(data$vcmax_obs,na.rm  = TRUE)^2
        if(all(is.na(data$jmax_obs))){
          y5=0
        }else{
          y5 = mean((calc_jmax_arrhenius(dat1$jmax25,298.15,dat1$t+273.15) - data$jmax_obs)^2,na.rm  = TRUE)/
            mean(data$jmax_obs,na.rm  = TRUE)^2
        }
        
        y=y2+y1+y3+y4+y5

        cat(x, "|", y, "\n")
        
        y
        
      }
    }
  }
}




par_scheme_gamma <- list("PHYDRO","CGAIN", "CMAX", "SOX2")
par_scheme_no_gamma <- list("PROFITMAX2","SOX","PROFITMAX")

##### PARAMETERIZATION #####
get_parameters_kmax_alpha_vcmax <- function(x){
    species = x$Species %>% unique()
    stomatal_model_now = x$scheme %>% unique()
    inst = x$inst
    data_template_now = x
    
    data1 = filter(dat, Species==species, Source == unique(x$source))

    
    ##### PARAMETERIZATION WITH ACCLIMATION #####
    print(stomatal_model_now)
    print(species)
    
    if(stomatal_model_now %in% par_scheme_no_gamma){
      parameter_ini <- c(1,0.1,-1,2) #hydraulic parameters and alpha
      optim(fn = error_fun_kmax_alpha,
            par = parameter_ini,
            data=data1,
            data_template = data_template_now,
            stomatal_model = stomatal_model_now,
            Species_now = species,
            control = list(maxit = 500,#, maximize = TRUE, 
                           parscale = c(1,0.01,0.1,0.1),
                           REPORT=0, trace=0, reltol=1e-3)
      ) -> opt_accl         
      x_accl <- opt_accl$par
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         alpha=x_accl[2],
                         gamma=NA,
                         p50_opt=x_accl[3],
                         b_opt=x_accl[4])
      error_fun_kmax_alpha(x_accl, data1, 
                           data_template = data_template_now,
                           plot=T, 
                           stomatal_model = stomatal_model_now, 
                           Species_now = species)
    }
    
    if(stomatal_model_now %in% par_scheme_gamma){
      parameter_ini <- c(1,1,0.1,-1,2) #hydraulic parameters and alpha
      optim(fn = error_fun_kmax_alpha_gamma,
            par = parameter_ini,
            data=data1,
            data_template = data_template_now,
            stomatal_model = stomatal_model_now,
            Species_now = species,
            control = list(maxit = 500,#, maximize = TRUE, 
                           parscale = c(1,1,0.01,0.1,0.1),
                           REPORT=0, trace=0, reltol=1e-3)
      ) -> opt_accl         
      x_accl <- opt_accl$par
      res_accl <- tibble(x,
                         acclimation = TRUE,
                         K.scale=x_accl[1],
                         alpha=x_accl[3],
                         gamma=x_accl[2],
                         p50_opt=x_accl[4],
                         b_opt=x_accl[5])
      error_fun_kmax_alpha_gamma(x_accl, data1, 
                           data_template = data_template_now,
                           plot=T, 
                           stomatal_model = stomatal_model_now, 
                           Species_now = species)
    }
    
    
  df <- res_accl
  
  readr::write_csv(df,file=paste0("DATA/parameters_kmax_alpha_vcmax/",stomatal_model_now,"_",species,"_",x$source,".csv"))
  
  return(res_accl)

}

##### CALIBRATE PARAMETERS #####

template %>% 
  # filter(scheme %in% c("CMAX")
  # # #        # Species %in% c(
  # # #        #   # "Rosa cymosa",
  # # #        #   # "Broussonetia papyrifera",
  # # #        #   # "Cinnamomum bodinieri",
  # # #        #   # "Platycarya longipes",
  # # #        #   # "Pteroceltis tatarinowii"
  # # #        #   # "Picea abies",
  # # #        #   # "Betula pendula",
  # # #        #   # "Pinus sylvestris",
  # # #        #   # "Populus tremula"
  # # #        # )
  #        ) %>%
  # slice(32:39) %>% 
  # filter(scheme %in% c("CGAIN")) %>%
  # filter(Species == "Quercus ilex", source == "Epron and Dreyer (1990)") %>%
  group_split(scheme, Species,source) %>% 
  purrr::map_df(get_parameters_kmax_alpha_vcmax)
