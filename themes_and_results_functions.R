# PLOTTING THEMES --------------------------------------------------------------
makeTransparent = function(col, alpha=0.7){
  rgb(t(col2rgb(col)/255), alpha = alpha)
}

# show_col(brewer_pal(palette = "BrBG")(5))

mytheme = function(){
  theme_classic()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          plot.tag.position = "topleft") 
  
}

mytheme2 = function(){
  theme_classic()+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=10),
          legend.text=element_text(size=10),
          plot.tag.position = "topleft") 
  
}

mytheme3 = function(){
  theme_classic()+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.tag.position = "topleft") 
  
}

mytheme4 = function(){
  theme_classic()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=16),
          legend.text=element_text(size=16),
          plot.tag.position = "topleft") 
  
}

mytheme5 = function(){
  theme_classic()+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          legend.text=element_text(size=16),
          plot.tag.position = "topleft") 
  
}


col_df <- tibble(scheme = factor(c('CMAX','CGAIN','PHYDRO','PROFITMAX2','SOX','PROFITMAX')),
                 col = brewer_pal(palette = "Dark2")(6)
)


# WORK FUNCTIONS ---------------------------------------------------------------
get_vcmax_jmax_ww <- function(x){
  data_ww <- x %>% 
    filter(!is.na(gC)) %>% 
    mutate(LWP_q90 = case_when(source == "Galmes et al. (2007)"~quantile(LWP, 0.5,na.rm = TRUE), # we consider quantile 0.5 because there is only 4 data points
                               TRUE~quantile(LWP, 0.9,na.rm = TRUE)),
           ci = ca-A/gC) %>% 
    filter(LWP>=LWP_q90) %>% 
    dplyr::select(LWP,A,gC,T,ci,Iabs_growth) %>% 
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  
  vcmax <- calc_vcmax_no_acclimated_ww(A = data_ww$A,
                                       ci = data_ww$ci,
                                       tc = data_ww$T,
                                       patm = calc_patm(0,data_ww$T),
                                       rdark = 0.02
  )
  jmax <- calc_jmax_no_acclimated_ww(A = data_ww$A,
                                     vcmax = vcmax,
                                     ci = data_ww$ci,
                                     I = data_ww$Iabs_growth,
                                     tc = data_ww$T,
                                     patm = calc_patm(0,data_ww$T),
                                     kphio = 0.087)
  return(tibble(vcmax_ww = vcmax,jmax_ww = jmax,LWP_ww = data_ww$LWP))
}

get_no_acclimated_response <- function(x,tc_now,ppfd_now,vpd_now,co2_now){
  stomatal_model = x$stomatal_model
  par_plant_now = list(
    conductivity = x$K.scale*1e-16,
    psi50 = x$P50%>% unique(),
    b = x$b%>% unique()
  )
  par_cost_now = list(
    alpha = 0.1,
    gamma = x$gamma
  )
  ndays = 30
  psi_max = 0
  psi_min = -6
  lwp = seq(psi_min, 0, length.out=50)
  
  dat_a <- rpmodel::rpmodel(tc =tc_now, ppfd = ppfd_now, 
                            vpd = vpd_now, co2 = co2_now, 
                            elv = 0, fapar = .99)
  
  dat1 = tibble(var = lwp, jmax_a=dat_a$jmax, vcmax_a=dat_a$vcmax) %>%
    mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                           ~model_numerical_instantaneous(tc =tc_now, ppfd = ppfd_now, 
                                                          vpd = vpd_now, co2 = co2_now, 
                                                          elv = 0, fapar = .99, 
                                                          kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.02, 
                                                          par_plant=par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax = ..2, vcmax = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
  dat2 <- dat1 %>% filter(gs>=1e-40)
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
  return(dat1 %>% cbind(scheme = x$scheme, Species = x$Species, psi88S=psi88S,psiL88S=psiL88S,gs0 = gs0))
  
}



get_acclimated_response <- function(x,tc_now,ppfd_now,vpd_now,co2_now){
  stomatal_model = x$stomatal_model
  par_plant_now = list(
    conductivity = x$K.scale*1e-16,
    psi50 = x$P50 %>% unique(),
    b = x$b %>% unique()
  )
  par_cost_now = list(
    alpha = x$alpha,
    gamma = x$gamma
  )
  ndays = 30
  psi_max = 0
  psi_min = -7
  lwp = seq(psi_min, 0, length.out=100)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  
  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }
  
  k = 7
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
  
  spl = splinefun(x = max(day):0, y=lwp_week)
  
  dat_acc = tibble(var = spl(day)) %>% 
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           pmod = map(var, ~model_numerical(tc =tc_now, ppfd = ppfd_now, 
                                            vpd = vpd_now, co2 = co2_now, 
                                            elv = 0, fapar = .99, kphio = 0.087, 
                                            psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                            par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  dat1 = tibble(var = lwp, jmax_a=dat_acc$jmax, vcmax_a=dat_acc$vcmax) %>%
    mutate(p = purrr::pmap(list(var, jmax_a, vcmax_a), 
                           ~model_numerical_instantaneous(tc =tc_now, ppfd = ppfd_now, 
                                                          vpd = vpd_now, co2 = co2_now, 
                                                          elv = 0, fapar = .99, 
                                                          kphio = 0.087, 
                                                          psi_soil = ..1, rdark = 0.02, 
                                                          par_plant=par_plant_now, 
                                                          par_cost = par_cost_now, 
                                                          jmax = ..2, vcmax = ..3, 
                                                          stomatal_model = stomatal_model)) ) %>% 
    unnest_wider(p)
  dat2 <- dat1 %>% filter(gs>=1e-40)
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
  return(dat1 %>% cbind(scheme = x$scheme, Species = x$Species, psi88S=psi88S,psiL88S=psiL88S,gs0 = gs0))
  
}


get_cost <- function(x,tc_now,ppfd_now,vpd_now,co2_now){
  stomatal_model = x$stomatal_model
  par_plant_now = list(
    conductivity = x$K.scale*1e-16,
    psi50 = x$P50%>% unique(),
    b = x$b%>% unique()
  )
  par_cost_now = list(
    alpha = x$alpha,
    gamma = x$gamma
  )
  ndays = x$Drydown.days
  psi_max = 0
  psi_min = -7
  lwp = seq(psi_min, 0, length.out=20)
  day = ndays * (lwp-psi_max)/(psi_min-psi_max)
  
  lwp_day = function(day_num){
    psi_max + day_num/ndays * (psi_min-psi_max)
  }
  
  k = 7
  lwp_week = rollmean(x = lwp_day(c(max(day):0, rep(0,k-1))), k = k, align = "right")
  
  spl = splinefun(x = max(day):0, y=lwp_week)
  
  dat_acc = tibble(var = spl(day)) %>% 
    mutate(var = case_when(var>0~0,
                           TRUE~var),
           pmod = map(var, ~model_numerical(tc =tc_now, ppfd = ppfd_now, 
                                            vpd = vpd_now, co2 = co2_now, 
                                            elv = 0, fapar = .99, kphio = 0.087, 
                                            psi_soil = ., rdark = 0.02, par_plant=par_plant_now, 
                                            par_cost = par_cost_now, stomatal_model = stomatal_model))) %>% 
    unnest_wider(pmod)
  
  
  return(dat_acc %>% cbind(scheme = x$scheme, Species = x$Species, source = x$source))
  
}