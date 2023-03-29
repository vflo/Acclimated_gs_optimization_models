################## optimization code ###################################

rm(list=ls())


library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(zoo)
library(stringr)
library(rphydro)
library(DEoptim)
library(dplyr)
select = dplyr::select 
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
  
  grid.arrange(p3,p6,p1,p2,p5,p4, ncol=2)
}


################################################################################


error_fun_kmax_no_alpha = function(x, 
                                data, 
                                data_template, 
                                plot=F, 
                                k=7, 
                                stomatal_model = stomatal_model_now, 
                                vcmax25 = vcmax25,
                                jmax25 = jmax25, 
                                Species_now = species){
  
  if(x[1]<=0| x[2]>-0.5| x[3]<1
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
      conductivity = x[1]*1e-16,
      psi50 = x[2],#data_template$P50 %>% unique(),
      b = x[3]#data_template$b%>% unique()
    )
    par_cost_now = list(
      alpha  = NA
    )
    
    lwp = data$LWP
    dat1 = try(
      tibble(var = lwp, jmax25_a=jmax25, vcmax25_a=vcmax25) %>% 
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
                                                              stomatal_model = stomatal_model)) 
        ) %>% 
        unnest_wider(p),
      silent = TRUE)
    if(any(class(dat1) %in% "try-error")){
      print(paste("error try"))
      return(1e2)
    }else{
      if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=Species_now, 
                                    data = data)

      
      y2 = mean((dat1$gs - data$gC)^2,na.rm  = TRUE)/mean(data$gC,na.rm  = TRUE)^2
      y1 = mean((dat1$a - data$A)^2,na.rm  = TRUE)/mean(data$A,na.rm  = TRUE)^2
      y4 = mean((dat1$chi - (data$Ciest/data$ca_pa))^2,na.rm  = TRUE)/
        mean((data$Ciest/data$ca_pa),na.rm  = TRUE)^2

      
      y=y2+y1+y4

      cat(x, "|", y, "\n")
      
      y
      
    }
  }
}



error_fun_kmax_no_alpha_gamma = function(x, 
                                      data, 
                                      data_template, 
                                      plot=F, 
                                      k=7, 
                                      stomatal_model = stomatal_model_now,
                                      vcmax25 = vcmax25,
                                      jmax25 = jmax25, 
                                      Species_now = species){
  
  if(x[1]<=0| x[2]<=0| x[3]>-0.5| x[4]<1
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
      psi50 = x[3],#data_template$P50%>% unique(),
      b = x[4]#data_template$b%>% unique()
    )
    
    par_cost_now = list(
      alpha = NA,
      gamma = x[2]
    )
    
    lwp = data$LWP
    dat1 = try(
      tibble(var = lwp, jmax25_a=jmax25, vcmax25_a=vcmax25) %>% 
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
                                                              stomatal_model = stomatal_model)) 
        ) %>% 
        unnest_wider(p),
      silent = TRUE)
    if(any(class(dat1) %in% "try-error")){
      print(paste("error try"))
      return(1e2)
    }else{
      if(plot==T) dat1 %>% plot_all(varname = "psi_soil", species=Species_now, 
                                    data = data)

      y2 = mean((dat1$gs - data$gC)^2,na.rm  = TRUE)/mean(data$gC,na.rm  = TRUE)^2
      y1 = mean((dat1$a - data$A)^2,na.rm  = TRUE)/mean(data$A,na.rm  = TRUE)^2
      y4 = mean((dat1$chi - (data$Ciest/data$ca_pa))^2,na.rm  = TRUE)/
        mean((data$Ciest/data$ca_pa),na.rm  = TRUE)^2

      y=y2+y1+y4
      
      cat(x, "|", y, "\n")
      
      y
      
    }
  }
}

par_scheme_gamma <- list("PHYDRO","CGAIN", "CMAX", "SOX2")
par_scheme_no_gamma <- list("PROFITMAX2","SOX","PROFITMAX")

##### PARAMETERIZATION #####
get_parameters_kmax_no_alpha <- function(x){
  
    species = x$Species %>% unique()
    stomatal_model_now = x$scheme %>% unique()
    inst = x$inst
    data_template_now = x
    
    data1 = filter(dat, Species==species, Source == unique(x$source))

    ##### PARAMETERIZATION WITHOUT ACCLIMATION WW #####
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
                                         rdark = 0.020
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
    if(is.na(jmax25)){jmax25 = 1.6*vcmax25} #some Jmax calculation may fail using 1.6 relation
    
    ##### PARAMETERIZATION WITHOUT ACCLIMATION #####
    print(stomatal_model_now)
    print(species)
    if(stomatal_model_now %in% par_scheme_no_gamma){
      parameter_ini <- c(1,-1,2) #hydraulic parameters and alpha
      optim(fn = error_fun_kmax_no_alpha,
            par = parameter_ini,
            data=data1,
            data_template = data_template_now,
            stomatal_model = stomatal_model_now,
            vcmax25 = vcmax25,
            jmax25 = jmax25, 
            Species_now = species,
            control = list(maxit = 500,#, maximize = TRUE, 
                           parscale = c(1,0.1,0.1),
                           REPORT=0, trace=0, reltol=1e-3)
      ) -> opt_accl         
      x_accl <- opt_accl$par
      res <- tibble(x,
                    acclimation = TRUE,
                    K.scale=x_accl[1],
                    alpha=NA,
                    gamma=NA,
                    p50_opt=x_accl[2],
                    b_opt=x_accl[3])
      error_fun_kmax_no_alpha(x_accl, data1, 
                           data_template = data_template_now,
                           plot=T, 
                           stomatal_model = stomatal_model_now,
                           vcmax25 = vcmax25,
                           jmax25 = jmax25, 
                           Species_now = species)
    }
    
    if(stomatal_model_now %in% par_scheme_gamma){
      parameter_ini <- c(1,1,-1,2) #hydraulic parameters and alpha
      optim(fn = error_fun_kmax_no_alpha_gamma,
            par = parameter_ini,
            data=data1,
            data_template = data_template_now,
            stomatal_model = stomatal_model_now,
            vcmax25 = vcmax25,
            jmax25 = jmax25, 
            Species_now = species,
            control = list(maxit = 500,#, maximize = TRUE, 
                           parscale = c(1,1,0.1,0.1),
                           REPORT=0, trace=0, reltol=1e-3)
      ) -> opt_accl         
      x_accl <- opt_accl$par
      res <- tibble(x,
                    acclimation = TRUE,
                    K.scale=x_accl[1],
                    alpha=NA,
                    gamma=x_accl[2],
                    p50_opt=x_accl[3],
                    b_opt=x_accl[4])
      error_fun_kmax_no_alpha_gamma(x_accl, data1, 
                                 data_template = data_template_now,
                                 plot=T, 
                                 stomatal_model = stomatal_model_now,
                                 vcmax25 = vcmax25,
                                 jmax25 = jmax25, 
                                 Species_now = species)
    }
    
  df <- res
  
  readr::write_csv(df,file=paste0("DATA/parameters_kmax_no_alpha/",stomatal_model_now,"_",species,"_",x$source,".csv"))
  # 
  # return(res)

}

##### CALIBRATE PARAMETERS #####
template %>% 
  filter(scheme %in% c("CMAX")#scheme %in% c("CGAIN2") #scheme == "CGAIN",Species == "Diplotaxis ibicensis"
  #        Species %in% c( "Diplotaxis ibicensis",
  #                        'Helianthus annuus',
  #                        "Malva subovata"
  # #        #   # "Rosa cymosa",
  # #        #   # "Broussonetia papyrifera",
  # #        #   # "Cinnamomum bodinieri",
  # #        #   # "Platycarya longipes",
  # #        #   # "Pteroceltis tatarinowii"
  # #        #   # "Picea abies",
  # #        #   # "Betula pendula",
  # #        #   # "Pinus sylvestris",
  # #        #   # "Populus tremula"
  #        )
         ) %>%
  # filter(scheme %in% c("CGAIN")) %>%
  # filter(Species == "Quercus ilex", source == "Epron and Dreyer (1990)") %>%
  group_split(scheme, Species, source) %>% 
  purrr::map_df(get_parameters_kmax_no_alpha)
