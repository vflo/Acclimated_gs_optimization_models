###########################################
## ACCLIMATION OPTIMIZATION SCHEMES
###########################################

fn_profit <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (log(jmax) is supplied by the optimizer)
  dpsi = exp(par[2])#      # delta Psi in MPa (log(dpsi) is supplied)
  psi_leaf = psi_soil-dpsi #MPa
  
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s (include residual gs to avoid A calculation error)
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j = calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a = a_j$a
  ci = a_j$ci
  # vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth) # vcmax based on coordination theory

  ## PHYDRO
  if(stomatal_model == "PHYDRO"){
    out = a - par_cost$alpha * jmax - par_cost$gamma * dpsi^2 
  }

  ## PROFITMAX2
  if(stomatal_model == "PROFITMAX2"){
    cost = (par_cost$alpha*jmax) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b) #Calculate e over e_crit
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = a * e_e_crit
    out = exp(a - out_e - cost) #exp is used to optimize only in the positive space
  }
  
  # LEAST_COST
  if(stomatal_model == "LEAST_COST"){
    out = exp((a)/(par_cost$gamma*e*1e3+vcmax)-(par_cost$alpha*jmax)) #umolco2 umolh2o-1 m-2 s-1
  }
  
  #CGAIN
  if(stomatal_model == "CGAIN"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    out = exp(a - par_cost$alpha*jmax - par_cost$gamma*(K-kl)/K)
  }

  # WUE
  if(stomatal_model == "WUE"){
    out = exp(a - par_cost$alpha*jmax - par_cost$gamma*e*1e3)
  }

  ## CMAX
  if(stomatal_model == "CMAX"){
    bb = 1 # as for Sabot et al. 2022
    out = ((a - par_cost$gamma*psi_leaf^2 - bb*psi_leaf) - (par_cost$alpha * jmax))
  }
  
  ## PROFITMAX
  if(stomatal_model == "PROFITMAX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    g = e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a = calc_assimilation_limiting(vcmax, jmax, g, par_photosynth) #calculate maximum assimilation at e_crit
    amax =  c_a$a
    ks = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    out = exp(a - amax*(1-kl/ks) - par_cost$alpha * jmax)
  }

  ## SOX
  if(stomatal_model == "SOX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out = exp((a)*reduction_factor - (par_cost$alpha * jmax))
  }

  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}


##############################################
# INSTANTANEOUS OPTIMIZATION SCHEMES
##############################################
#Jmax is no longer included as a cost in the instantaneous optimizations since it is fixed
fn_profit_inst_schemes <- function(par, jmax, vcmax, psi_soil, e_crit, p_crit, par_cost, 
                                   par_photosynth, par_plant, par_env, 
                                   stomatal_model, do_optim){
  dpsi = exp(par[1])#      # delta Psi in MPa (log(dpsi) is supplied)
  psi_leaf = psi_soil-dpsi #MPa

  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270  # gs in mol_co2/m2/s/Mpa
  e = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol_h2o/m2_ground/s
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a
  
  ## PHYDRO
  if(stomatal_model == "PHYDRO"){
    profit = A - par_cost$gamma * dpsi^2 
  }

  ## PROFITMAX2
  if(stomatal_model == "PROFITMAX2"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b)
    if(is.na(e_e_crit)){e_e_crit = 1} #avoid errors when e is close to e_crit
    profit = exp(A*(1-e_e_crit))
  }
  
  ## LEAST_COST 
  if(stomatal_model == "LEAST_COST"){
    profit = exp(A/(par_cost$gamma*e*1e3+vcmax))
  }
  
  ## CGAIN
  if(stomatal_model == "CGAIN"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    profit = exp(A - par_cost$gamma*(K-kl)/K)
  }
  
  ## WUE
  if(stomatal_model == "WUE"){
    profit = exp(A - (par_cost$gamma*e*1e3))
  }

  ## CMAX
  if(stomatal_model == "CMAX"){
    bb = 1
    profit = (A - par_cost$gamma*si_leaf^2 - bb*si_leaf)
  }
  
  ## PROFITMAX
  if(stomatal_model == "PROFITMAX"){
    K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    g = e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a = calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax =c_a$a
    ks = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    profit = exp(A - amax*(1-kl/ks))
  }
  
  ## SOX
  if(stomatal_model == "SOX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    profit = exp(A*reduction_factor)
  }
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


##############################################################
# INSTANTANEOUS OPTIMIZATION SCHEMES SEPARATED COST CALCULATOR
##############################################################
fn_profit_schemes_cost_calc <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                        par_plant, par_env, stomatal_model){
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = exp(par[2])     # delta Psi in MPa (log(dpsi) is supplied)
  psi_leaf = psi_soil-dpsi #MPa

  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j = calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a = a_j$a
  ci = a_j$ci
  # vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth) # vcmax based on coordination theory
  
  ## phydro
  if(stomatal_model == "PHYDRO"){
    out <- tibble(jmax_cost = par_cost$alpha * jmax, h_cost = par_cost$gamma * dpsi^2)
  }
  
  ## PROFITMAX2
  if(stomatal_model == "PROFITMAX2"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b) #Calculate e over e_crit
    if(is.na(e_e_crit)){e_e_crit = 1}
    out = tibble(jmax_cost = par_cost$alpha * jmax, h_cost = a * e_e_crit)
  }
  
  ## CGAIN
  if(stomatal_model == "CGAIN"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    out = tibble(jmax_cost = par_cost$alpha * jmax, h_cost = par_cost$gamma*(K-kl)/K)
  }
  
  ## CMAX
  if(stomatal_model == "CMAX"){
    bb = 1
    out = tibble(jmax_cost = par_cost$alpha * jmax, 
                 h_cost = par_cost$gamma*psi_leaf^2 - bb*psi_leaf)
  }
  
  ## PROFITMAX
  if(stomatal_model == "PROFITMAX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    g = e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a = calc_assimilation_limiting(vcmax, jmax, g, par_photosynth)
    amax = c_a$a
    ks = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    out = tibble(jmax_cost = par_cost$alpha * jmax, h_cost = amax*(1-kl/ks)k)
  }

  ## SOX
  if(stomatal_model == "SOX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- tibble(jmax_cost = par_cost$alpha * jmax, h_cost = a*reduction_factor)
  }
    return(out)
}




##########################################
# ACCLIMATED OPTIMIZATION PROCESS
##########################################
optimise_stomata_phydro_schemes <- function(fn_profit, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  if(is.null(jmax_lim)){jmax_lim = 7}

  jmax_ini = c(0) # logarithm of Jmax
  dpsi_ini = c(0) # logarithm of 1 MPa
  lower <- c(-10, -10)
  
  ## This three models ('SOX','CMAX','CGAIN') have some problems when psi decreases and 
  ## the optimr algorithm can not find the initial gradient. We use 
  ## evolutionary global optimization via the Differential Evolution algorithm
  ## This algorithm is really slow... But performs better
  if(!stomatal_model %in% c('SOX','CMAX','CGAIN')){ 
    out_optim <- optimr::optimr(
      par       = c(jmax_ini, dpsi_ini),
      lower     = lower,
      upper     = c(jmax_lim, 20),
      fn             = fn_profit,
      psi_soil       = psi_soil,
      e_crit         = e_crit,
      p_crit         = p_crit,
      par_cost       = par_cost,
      par_photosynth = par_photosynth,
      par_plant      = par_plant,
      par_env        = par_env,
      do_optim       = do_optim, 
      stomatal_model = stomatal_model,
      method         = "L-BFGS-B",
      control        = list(maxit = 500, maximize = TRUE, fnscale = 1e10,
                            REPORT=0, trace=0)
    )
    out_optim$value <- -out_optim$value
  }else{
    out_optim <- DEoptim::DEoptim(
      lower     = lower,
      upper     = c(jmax_lim, 20),
      fn             = fn_profit,
      psi_soil       = psi_soil,
      e_crit         = e_crit,
      p_crit         = p_crit,
      par_cost       = par_cost,
      par_photosynth = par_photosynth,
      par_plant      = par_plant,
      par_env        = par_env,
      do_optim       = do_optim, 
      stomatal_model = stomatal_model,
      control = DEoptim.control(NP = 30,trace = FALSE,itermax = 200)
    )
    }

  
  if (return_all){
    out_optim
  } else {
    if(!stomatal_model %in% c('SOX','CMAX','CGAIN')){
      return(out_optim$par)
    }else{
      return(out_optim$optim$bestmem)
    }
  }
}


#############################################
# INSTANTANEOUS OPTIMIZATION PROCESS
#############################################
optimise_shortterm_schemes <- function(fn_profit_inst, jmax, vcmax, psi_soil, e_crit, p_crit,
                               par_cost, par_photosynth, par_plant, par_env, 
                               stomatal_model, return_all = FALSE, do_optim){
  
  dpsi_ini = 0 # logarithm of 1 MPa
  
  ## This three models ('SOX','CMAX','CGAIN') have some problems when psi decreases and 
  ## the optimr algorithm can not find the initial gradient. We use 
  ## evolutionary global optimization via the Differential Evolution algorithm
  ## This algorithm is really slow... But performs better
  if(!stomatal_model %in% c('SOX','CMAX','CGAIN')){
    out_optim <- optimr::optimr(
      par       = dpsi_ini,
      lower     = c(-10),
      upper     = c(20),
      fn        = fn_profit_inst,
      psi_soil  = psi_soil,
      jmax      = jmax,
      vcmax     = vcmax,
      e_crit    = e_crit,
      p_crit    = p_crit,
      par_cost  = par_cost,
      par_photosynth = par_photosynth,
      par_plant = par_plant,
      par_env   = par_env,
      do_optim  = do_optim, 
      stomatal_model = stomatal_model,
      method    = "L-BFGS-B",
      control   = list(maxit = 500, maximize = TRUE, fnscale = 1e1000)
    )
    out_optim$value <- -out_optim$value
  }else{
    out_optim <- DEoptim::DEoptim(
      lower     = c(-10),
      upper     = c(20),
      fn        = fn_profit_inst,
      psi_soil  = psi_soil,
      jmax      = jmax,
      vcmax     = vcmax,
      e_crit    = e_crit,
      p_crit    = p_crit,
      par_cost  = par_cost,
      par_photosynth = par_photosynth,
      par_plant = par_plant,
      par_env   = par_env,
      do_optim  = do_optim, 
      stomatal_model = stomatal_model,
      control = DEoptim.control(NP = 15,trace = FALSE,itermax = 100)
    )
  }
  
  if (return_all){
    out_optim
  } else {
    if(!stomatal_model %in% c('SOX','CMAX','CGAIN')){
      return(out_optim$par)
    }else{
        return(out_optim$optim$bestmem)
      }
  }
}


###########################################################
## MODEL OPTIMIZATION AND DATA PREPARATION WITH ACCLIMATION
###########################################################
# This function has as input the environmental variables, 
# as well as the hydraulic and cost parameters and the stomatal model to be calculated. 
# The optimisation process of jmax and dpsi is then carried out. 
# The function returns jmax, dpsi, psi_leaf, gs, E, A, ci, chi and vcmax.
model_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, 
                            rdark, par_plant, par_cost, stomatal_model){
  # 1. input preparation
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
  )
  
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )
  
    # Calculate psi and critical E which are based on the species parameters and 
    # for critical E also on the water potential in the soil.
    K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (ground) s-1 MPa-1
    p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
    e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b)+1e-270 #mol m-2 (ground) s-1

    # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit 
    # low_swp <- FALSE
    # if(psi_soil <= p_crit){
    #   psi_soil = p_crit*0.95
    #   low_swp  = TRUE
    #   print("WARNING: soil water potential is lower than critical plant water potential.")
    # }

    # 3. Optimizer
    lj_dps = optimise_stomata_phydro_schemes(fn_profit, 
                                             psi_soil = psi_soil,
                                             e_crit = e_crit,
                                             p_crit = p_crit,
                                             par_cost  = par_cost, 
                                             par_photosynth = par_photosynth, 
                                             par_plant = par_plant, 
                                             par_env = par_env,
                                             jmax_lim = 7, 
                                             return_all = FALSE, 
                                             do_optim=TRUE, 
                                             stomatal_model = stomatal_model)
  
    # 4. Calculate output variables
    jmax  = exp(lj_dps[1]) %>% unname()
    dpsi  = exp(lj_dps[2]) %>% unname()
    psi_l = psi_soil-dpsi
    gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
    a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
    a     = a_j$a
    ci    = a_j$ci
    vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
    E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
    gs    = gs #transform to umol m-2(ground) s-1 Pa-1
    
    cost <- fn_profit_schemes_cost_calc(par = lj_dps, psi_soil = psi_soil,e_crit = e_crit,
                                        p_crit = p_crit,par_cost  = par_cost, 
                                        par_photosynth = par_photosynth, 
                                        par_plant = par_plant, par_env = par_env,
                                        stomatal_model = stomatal_model)
    jmax_cost = cost$jmax_cost
    h_cost = cost$h_cost 

    # 5. Prepare output list and return
    return(list(
      jmax         = jmax,
      dpsi         = dpsi,
      p_leaf       = psi_l,
      e_crit       = e_crit,
      gs           = gs,
      E            = E,
      a            = a,
      ci           = ci,
      chi          = ci/par_photosynth$ca,
      vcmax        = vcmax,
      jmax_cost    = jmax_cost,
      h_cost       = h_cost
    ))
}


########################################################
## INSTANTANEOUS MODEL OPTIMIZATION AND DATA PREPARATION
########################################################
# This function has as input the environmental variables, 
# as well as the hydraulic and cost parameters and the stomatal model to be calculated. 
# The optimisation process of jmax and dpsi is then carried out. 
# The function returns jmax, dpsi, psi_leaf, gs, E, A, ci, chi and vcmax.
model_numerical_instantaneous <- function(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, 
                                          kphio, psi_soil,  rdark, par_plant, 
                                          par_cost, stomatal_model){
  
  # 1. input preparation
  patm = rpmodel::calc_patm(elv)
  par_photosynth <- list(
    kmm       = rpmodel::calc_kmm(tc, patm),
    gammastar = rpmodel::calc_gammastar(tc, patm),
    phi0      = kphio * rpmodel::calc_ftemp_kphio(tc),
    Iabs      = ppfd * fapar,
    ca        = co2 * patm * 1e-6,  # Convert to partial pressure
    patm      = patm,
    delta     = rdark
  )
  par_env = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    density_water   = rpmodel::calc_density_h2o(tc, patm),  # Needs to be imported from rpmodel.R
    patm            = patm,
    tc              = tc,
    vpd             = vpd
  )

  # Calculate psi and critical E which are based on the species parameters and 
  # for critical E also on the water potential in the soil.
  K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
  p_crit = par_plant$psi50 * (log(1000)/log(2)) ^ ( 1/par_plant$b)
  e_crit = K * -integral_P_ecrit(psi_soil, par_plant$psi50, par_plant$b)+1e-270 #mol m-2 (ground) s-1
    
  
  # 2. if soil psi is lower than p_crit, set soil psi to 95% p_crit
  # low_swp <- FALSE
  # if(psi_soil <= p_crit){
  #   psi_soil = p_crit*0.95
  #   low_swp  = TRUE
  #   print("WARNING: soil water potential is lower than critical plant water potential.")
  # }

  # 3. Optimizer
  lj_dps = optimise_shortterm_schemes(fn_profit_inst_schemes,
                                      jmax = jmax, 
                                      vcmax = vcmax,
                                      psi_soil = psi_soil,
                                      e_crit = e_crit,
                                      p_crit= p_crit,
                                      par_cost  = par_cost, 
                                      par_photosynth = par_photosynth, 
                                      par_plant = par_plant, 
                                      par_env = par_env,
                                      return_all = FALSE, 
                                      do_optim = TRUE, 
                                      stomatal_model = stomatal_model)
  
  # 4. Calculate output variables
  dpsi  = exp(lj_dps[1]) %>% unname() # delta Psi in MPa
  psi_l = psi_soil-dpsi  # leaf water potential
  gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
  a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  a     = a_j$a
  ci    = a_j$ci
  E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
  gs    = gs #transform to umol m-2(ground) s-1 Pa-1
  
  # 5. Prepare output list and return
  return(list(
      jmax         = jmax,
      dpsi         = dpsi,
      p_leaf       = psi_l,
      gs           = gs,
      E            = E,
      a            = a,
      ci           = ci,
      chi          = ci/par_photosynth$ca,
      vcmax        = vcmax
  ))

}
