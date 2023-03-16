###########################################
## ACCLIMATION OPTIMIZATION SCHEMES
###########################################

fn_profit <- function(par, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                      par_plant, par_env, do_optim = FALSE, stomatal_model){
  jmax25 = exp(par[1])  # Jmax25 in umol/m2/s (log(jmax25) is supplied by the optimizer)
  dpsi = exp(par[2])#      # delta Psi in MPa (log(dpsi) is supplied)
  psi_leaf = psi_soil-dpsi #MPa
  
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s (include residual gs to avoid A calculation error)
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  jmax <- calc_jmax_arrhenius(jmaxT1 = jmax25, T1 = 298.15, T2 = (par_env$tc + 273.15))
  a_j = calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a = a_j$a
  ci = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth) # vcmax based on coordination theory

  ## PHYDRO
  if(stomatal_model == "PHYDRO"){
    out = a - par_cost$alpha * jmax25 - par_cost$gamma * dpsi^2 
  }

  ## PROFITMAX2
  if(stomatal_model == "PROFITMAX2"){
    cost = (par_cost$alpha*jmax25) #umolco2 umolh2o m-2 s-1
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b) #Calculate e over e_crit
    if(is.na(e_e_crit)){e_e_crit = 1}
    out_e = a * e_e_crit
    out = exp(a - out_e - cost) #exp is used to optimize only in the positive space
  }
  
  # # LEAST_COST
  # if(stomatal_model == "LEAST_COST"){
  #   out = exp((a)/(par_cost$gamma*e*1e3+vcmax)-(par_cost$alpha*jmax25)) #umolco2 umolh2o-1 m-2 s-1
  # }
  
  #CGAIN
  if(stomatal_model == "CGAIN"){
    K      = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl     = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    out = exp(a - par_cost$alpha*jmax25 - par_cost$gamma*(K-kl)/K)
  }
  
  ## SOX2
  if(stomatal_model == "SOX2"){
    out = (a - a*(abs(psi_leaf)^par_plant$b/
                    (abs(par_plant$psi50*par_cost$gamma)^par_plant$b+
                       abs(psi_leaf)^par_plant$b)) - (par_cost$alpha * jmax25))
  }

  # WUE
  if(stomatal_model == "WUE"){
    out = exp(a - par_cost$alpha*jmax25 - par_cost$gamma*e*1e3)
  }

  # CMAX
  # if(stomatal_model == "CMAX"){
  #   foo <- pracma::grad(f=A_grad,x0=c(dpsi,jmax25),
  #                      heps = .Machine$double.eps^(1/3),
  #                      psi_soil, par_cost, par_plant, par_env, par_photosynth)
  #   # bb = 1
  #   # dpsi_grad <- (pracma::grad(f=A_grad_dpsi,x0=c(dpsi),
  #   #               heps = .Machine$double.eps^(1/3),
  #   #               jmax,psi_soil, par_cost, par_plant,
  #   #               par_env, par_photosynth)- par_cost$gamma*psi_leaf - bb)
  #   # jmax_grad <- (pracma::grad(f=A_grad_jmax,x0=c(jmax25),
  #   #               heps = .Machine$double.eps^(1/3),
  #   #               dpsi,psi_soil, par_cost, par_plant,
  #   #               par_env, par_photosynth) -par_cost$alpha)
  #   dpsi_grad <- foo[1]
  #   jmax_grad <- foo[2]
  #   
  #   
  #   out <- exp(-1*(abs(dpsi_grad) * abs(jmax_grad)))
  #   # print(paste("dpsi:",dpsi,'jmax25:',jmax25,"dpsi_grad:",dpsi_grad,"jmax_grad:",jmax_grad,'out:',out))
  # }

  
  ## PROFITMAX
  if(stomatal_model == "PROFITMAX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    g = e_crit/(1.6*(par_env$vpd/par_env$patm)) #mol m-2 (ground) s-1
    c_a = calc_assimilation_limiting(vcmax, jmax, g, par_photosynth) #calculate maximum assimilation at e_crit
    amax =  c_a$a
    ks = K *(1/2)^((psi_soil/par_plant$psi50)^par_plant$b)
    kl = K *(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    k_crit = 0
    out = exp(a - amax*(1-kl/ks) - par_cost$alpha * jmax25)
  }

  ## SOX
  if(stomatal_model == "SOX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out = exp((a)*reduction_factor - (par_cost$alpha * jmax25))
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
#jmax25 is no longer included as a cost in the instantaneous optimizations since it is fixed
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
  
  # ## LEAST_COST 
  # if(stomatal_model == "LEAST_COST"){
  #   profit = exp(A/(par_cost$gamma*e*1e3+vcmax))
  # }
  # 
  ## CGAIN
  if(stomatal_model == "CGAIN"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    profit = exp(A - par_cost$gamma*(K-kl)/K)
  }
  
  ## SOX2
  if(stomatal_model == "SOX2"){
    profit = (A-A*(abs(psi_leaf)^par_plant$b/(abs(par_plant$psi50*par_cost$gamma)^par_plant$b+abs(psi_leaf)^par_plant$b)))
  }
  
  ## WUE
  if(stomatal_model == "WUE"){
    profit = exp(A - (par_cost$gamma*e*1e3))
  }

  # CMAX
  if(stomatal_model == "CMAX"){
    # dpsi_dp = exp(par[1])+0.0001   #      # delta Psi in MPa (log(dpsi) is supplied)
    # gs_dp = calc_gs_phydro(dpsi_dp, psi_soil, par_plant, par_env)+1e-270  # gs in mol_co2/m2/s/Mpa
    # A_dp = calc_assimilation_limiting(vcmax, jmax, gs_dp, par_photosynth)$a
    bb = 0
    # profit = -1*((A_dp-A)/(0.0001) - par_cost$gamma*psi_leaf - bb) #calculated using gradient and multiply by -1 to do maximization instead of minimization
    profit =-1*abs(pracma::grad(f=A_grad_dpsi,x0=c(dpsi),
                     heps = .Machine$double.eps^(1/3),
                     jmax,psi_soil, par_cost, par_plant, 
                     par_env, par_photosynth)- par_cost$gamma*abs(psi_leaf) - bb)
    
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
  jmax25 = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = exp(par[2])     # delta Psi in MPa (log(dpsi) is supplied)
  psi_leaf = psi_soil-dpsi #MPa

  jmax <- calc_jmax_arrhenius(jmaxT1 = jmax25, T1 = 298.15, T2 = (par_env$tc + 273.15))
  
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s
  e  = 1.6*gs*(par_env$vpd/par_env$patm)         # E in mol/m2ground/s
  
  ## light-limited assimilation
  a_j = calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2ground/s
  a = a_j$a
  ci = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth) # vcmax based on coordination theory
  
  ## phydro
  if(stomatal_model == "PHYDRO"){
    out <- tibble(jmax25_cost = par_cost$alpha * jmax25, h_cost = par_cost$gamma * dpsi^2)
  }
  
  ## PROFITMAX2
  if(stomatal_model == "PROFITMAX2"){
    e_e_crit = -integral_P_e_ecrit(dpsi, psi_soil, par_plant$psi50,par_plant$b) #Calculate e over e_crit
    if(is.na(e_e_crit)){e_e_crit = 1}
    out = tibble(jmax25_cost = par_cost$alpha * jmax25, h_cost = a * e_e_crit)
  }
  
  ## CGAIN
  if(stomatal_model == "CGAIN"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    kl = K*(1/2)^((psi_leaf/par_plant$psi50)^par_plant$b)
    out = tibble(jmax25_cost = par_cost$alpha * jmax25, h_cost = par_cost$gamma*(K-kl)/K)
  }
  
  ## SOX2
  if(stomatal_model == "SOX2"){
    out = tibble(jmax25_cost = par_cost$alpha * jmax25, 
                 h_cost =  a*(abs(psi_leaf)^par_plant$b/(abs(par_plant$psi50*par_cost$gamma)^par_plant$b+abs(psi_leaf)^par_plant$b)))
  }
  
  ## CMAX
  if(stomatal_model == "CMAX"){
    bb = 1# as for Sabot et al. 2022
    out = tibble(jmax25_cost = par_cost$alpha * jmax25,
                 h_cost = 1/2*par_cost$gamma*psi_leaf^2 - bb*psi_leaf)
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
    out = tibble(jmax25_cost = par_cost$alpha * jmax25, h_cost = amax*(1-kl/ks))
  }

  ## SOX
  if(stomatal_model == "SOX"){
    K = scale_conductivity(par_plant$conductivity, par_env)  #mol m-2 (ground) s-1 MPa-1
    ks = K*(1/2)^(psi_soil/par_plant$psi50)^par_plant$b
    kl = K*(1/2)^(psi_leaf/par_plant$psi50)^par_plant$b
    k_crit = 0
    reduction_factor = (kl - k_crit)/(K - k_crit)
    out <- tibble(jmax25_cost = par_cost$alpha * jmax25, h_cost = a*reduction_factor)
  }
    return(out)
}




##########################################
# ACCLIMATED OPTIMIZATION PROCESS
##########################################
optimise_stomata_phydro_schemes <- function(fn_profit, psi_soil, par_cost, e_crit, p_crit, par_photosynth, 
                                            par_plant, par_env, jmax_lim, return_all = FALSE, do_optim=TRUE, stomatal_model){
  if(is.null(jmax_lim)){jmax_lim = 7}
  set.seed(7)
  jmax_ini = c(0) # logarithm of Jmax
  dpsi_ini = c(0) # logarithm of 1 MPa
  lower <- c(-10, -10)
  upper <- c(jmax_lim, 4)
  NP = 30
  itermax = 200
  

    out_optim <- DEoptim::DEoptim(
      lower          = lower,
      upper          = upper,
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
      control = DEoptim.control(NP = NP, trace = FALSE, itermax = itermax)
    )
  
  
  
  if (return_all){
    out_optim
  } else {
      return(out_optim$optim$bestmem)
  }
}


#############################################
# INSTANTANEOUS OPTIMIZATION PROCESS
#############################################
optimise_shortterm_schemes <- function(fn_profit_inst, jmax, vcmax, psi_soil, e_crit, p_crit,
                               par_cost, par_photosynth, par_plant, par_env, 
                               stomatal_model, return_all = FALSE, do_optim){
  
  dpsi_ini = 0 # logarithm of 1 MPa
  

    out_optim <- DEoptim::DEoptim(
      lower     = c(-10),
      upper     = c(4),
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
      control = DEoptim.control(NP = 15, trace = FALSE, itermax = 100)
    )
  # }
  
  if (return_all){
    out_optim
  } else {
      return(out_optim$optim$bestmem)
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
    visco_star      = rpmodel::calc_viscosity_h2o(tc, patm)/rpmodel::calc_viscosity_h2o(25, patm),
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

    if(stomatal_model == "CMAX"){
      dpsi_bounds = calc_dpsi_bound(psi_soil = psi_soil, par_plant = par_plant, 
                                    par_env = par_env, par_photosynth = par_photosynth, par_cost = par_cost)
      dpsi_max = dpsi_bounds$Iabs_bound
      u = try(uniroot(f = function(dpsi){dFdx(dpsi, psi_soil = psi_soil, par_plant = par_plant, 
                                          par_env = par_env, par_photosynth = par_photosynth, 
                                          par_cost = par_cost)$dP_dx}, 
                  interval = c(dpsi_max*0.01, dpsi_max*0.99)),
          silent = TRUE)
      # if(any(class(u) %in% "try-error")){
      #   # print(paste("error try"))
      #   dpsi=dpsi_max*0.01
      # }else{
      dpsi=u$root
      # }
      x = calc_x_from_dpsi(dpsi, psi_soil = psi_soil, par_plant = par_plant, 
                           par_env = par_env, par_photosynth = par_photosynth, 
                           par_cost = par_cost)
      psi_l = psi_soil-dpsi
      gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # In mol/m2/s
      E  = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
      J = calc_J(gs, x, par_photosynth)
      jmax = calc_jmax_from_J(J, par_photosynth)
      jmax25 = calc_jmax_arrhenius(jmaxT1 = jmax, T1 = (tc + 273.15), T2 = 298.15)
      a_j = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
      a = a_j$a
      ci = a_j$ci
      vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
      vcmax25 = calc_vcmax_arrhenius(vcmaxT1 = vcmax, T1 = (tc + 273.15), T2 = 298.15)
      jmax25_cost = NA
      h_cost = NA
    }else{
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
      jmax25  = exp(lj_dps[1]) %>% unname()
      dpsi  = exp(lj_dps[2]) %>% unname()
      psi_l = psi_soil-dpsi
      gs    = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env) # gs in mol m-2 (ground) s-1
      jmax = calc_jmax_arrhenius(jmaxT1 = jmax25, T1 = 298.15, T2 = (tc + 273.15))
      a_j   = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
      a     = a_j$a
      ci    = a_j$ci
      vcmax = calc_vcmax_coordinated_numerical(a, ci, par_photosynth)
      vcmax25 = calc_vcmax_arrhenius(vcmaxT1 = vcmax, T1 = (tc + 273.15), T2 = 298.15)
      E     = 1.6 * gs * par_env$vpd/patm  # E in mol m-2 (ground) s-1
      gs    = gs #transform to umol m-2(ground) s-1 Pa-1
      
      cost <- fn_profit_schemes_cost_calc(par = lj_dps, psi_soil = psi_soil,e_crit = e_crit,
                                          p_crit = p_crit,par_cost  = par_cost, 
                                          par_photosynth = par_photosynth, 
                                          par_plant = par_plant, par_env = par_env,
                                          stomatal_model = stomatal_model)
      jmax25_cost = cost$jmax25_cost
      h_cost = cost$h_cost 
    
    }


    # 5. Prepare output list and return
    return(list(
      jmax25         = jmax25,
      vcmax25        = vcmax25,
      dpsi           = dpsi,
      p_leaf         = psi_l,
      e_crit         = e_crit,
      gs             = gs,
      E              = E,
      a              = a,
      ci             = ci,
      chi            = ci/par_photosynth$ca,
      jmax25_cost    = jmax25_cost,
      h_cost         = h_cost
    ))
}
 
########################################################
## INSTANTANEOUS MODEL OPTIMIZATION AND DATA PREPARATION
########################################################
# This function has as input the environmental variables, 
# as well as the hydraulic and cost parameters and the stomatal model to be calculated. 
# The optimisation process of jmax and dpsi is then carried out. 
# The function returns jmax, dpsi, psi_leaf, gs, E, A, ci, chi and vcmax.
model_numerical_instantaneous <- function(vcmax25, jmax25, tc, ppfd, vpd, co2, elv, fapar, 
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

  jmax = calc_jmax_arrhenius(jmaxT1 = jmax25, T1 = 298.15, T2 =  (tc + 273.15))
  vcmax = calc_vcmax_arrhenius(vcmaxT1 = vcmax25,  T1 = 298.15, T2 = (tc + 273.15))
  
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
      jmax25         = jmax25,
      vcmax25        = vcmax25,
      dpsi           = dpsi,
      p_leaf         = psi_l,
      gs             = gs,
      E              = E,
      a              = a,
      ci             = ci,
      chi            = ci/par_photosynth$ca
  ))

}


# A_grad <- function(u, psi_soil, par_cost, par_plant, par_env, par_photosynth){
#   x = u[1]
#   y = u[2]
#   psi_leaf = psi_soil-x
#   gs = calc_gs_phydro(x, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s (include residual gs to avoid A calculation error)
#   jmax = calc_jmax_arrhenius(jmaxT1 = y, T1 = 298.15, T2 = (par_env$tc + 273.15))
#   a = calc_assim_light_limited(gs, jmax, par_photosynth)$a
#   b = 1 # as for Sabot et al. 2022
#   C = 0
#   gamma = par_cost$gamma
#   alpha = par_cost$alpha
#   a  - alpha * y - 1/2*gamma*abs(psi_leaf)^2 - b*abs(psi_leaf)
# }

A_grad_dpsi <- function(u, jmax, psi_soil, par_cost, par_plant, par_env, par_photosynth){
  x = u[1]
  # psi_leaf = psi_soil-x
  gs = calc_gs_phydro(x, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s (include residual gs to avoid A calculation error)
  a = calc_assim_light_limited(gs, jmax, par_photosynth)$a
  # b = 1 # as for Sabot et al. 2022
  # C = 0
  # gamma = par_cost$gamma
  # alpha = par_cost$alpha
  a
}

# A_grad_jmax <- function(u, dpsi, psi_soil, par_cost, par_plant, par_env, par_photosynth){
#   y = u[1]
#   jmax = calc_jmax_arrhenius(jmaxT1 = y, T1 = 298.15, T2 = (par_env$tc + 273.15))
#   # psi_leaf = psi_soil-dpsi
#   gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)+1e-270   # gs in mol/m2ground/s (include residual gs to avoid A calculation error)
#   a = calc_assim_light_limited(gs, jmax, par_photosynth)$a
#   # b = 1 # as for Sabot et al. 2022
#   # C = 0
#   # gamma = par_cost$gamma
#   # alpha = par_cost$alpha
#   a
# }


calc_dpsi_bound <- function(psi_soil, par_plant, par_env, par_photosynth, par_cost){
  gstar = par_photosynth$gammastar/par_photosynth$patm*1e6
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  y = par_cost$gamma
  K = scale_conductivity(par_plant$conductivity, par_env)
  K = K/(1.6*par_env$vpd/par_env$patm)
  Pox = P(psi_soil, par_plant$psi50, par_plant$b)
  Ppox = Pprime(psi_soil, par_plant$psi50, par_plant$b)
  Pppox = Pprimeprime(psi_soil, par_plant$psi50, par_plant$b)
  
  f1 = function(dpsi){
    gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)
    x = calc_x_from_dpsi(dpsi,psi_soil,par_plant, par_env, par_photosynth, par_cost)
    # if(x == 1){
    #   J = 0
    # }else{
      J = calc_J(gs, x, par_photosynth)-4*par_photosynth$phi0*par_photosynth$Iabs
    # }
    
    J
  }
  
  a = (ca + 2*gstar)*K*Pppox*4/8
  b = -(y + (ca + 2*gstar)*K*Ppox)
  c = (ca + 2*gstar)*K*Pox
  del = b^2-4*a*c
  
  # approx_O2 = (-b-sqrt(b^2-4*a*c))/2/a
  # exact = try(uniroot(f = function(dpsi){(-(y*abs(psi_soil - dpsi)+0) + (ca + 2*gstar)*
  #                                       calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, 
  #                 interval = c(0,10))$root,silent = TRUE)
  # 
  # if(any(class(exact) %in% "try-error")){
  #   # print(paste("error try"))
  #   use_bound = approx_O2
  #   exact = "error"
  # }else{
  #   use_bound = exact
  # }
  
  # 
  if (del>=0){
    approx_O2 = (-b-sqrt(b^2-4*a*c))/2/a

    del_x = (-2*approx_O2*y + (ca + 2*gstar)*calc_gsprime(approx_O2, psi_soil, par_plant, par_env))
    ajmax_m_phi0iabs = f1(approx_O2)*f1(approx_O2*0.0001)
    cat("signs of ajmax-phi0Iabs = ", ajmax_m_phi0iabs,"\n")
    if (del_x < 0 | ajmax_m_phi0iabs < 0){
      cat("Warning: Need exact calc: delx = ", del_x, ", ajm_by_phi0iabs = ", ajmax_m_phi0iabs,"\n")
      exact = uniroot(f = function(dpsi){(-(y*abs(psi_soil - dpsi)+0) + (ca + 2*gstar)*
                                            calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
      use_bound = exact
    }
    else{
      use_bound = approx_O2
      ## v--- Only for debug sake
      exact = uniroot(f = function(dpsi){(-(y*abs(psi_soil - dpsi)+0) + (ca + 2*gstar)*
                                             calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, interval = c(0,10))$root
    }

  }
  else{
    approx_O2 = NA
    exact = try(uniroot(f = function(dpsi){(-(y*abs(psi_soil - dpsi)+0) + (ca + 2*gstar)*
                                 calc_gsprime(dpsi, psi_soil, par_plant, par_env))}, 
                        interval = c(0,10))$root,
                silent= TRUE)
    
    if(any(class(exact) %in% "try-error")){
      cat("Warning: Need approx \n")
      exact = "error"
      approx = (ca*K*Pox + 2*gstar*K*Pox)/(ca*K*Ppox +
                                             2*gstar*K*Ppox + y)
      use_bound = approx
    }else{
      cat("Warning: Need exact calc: approx_O2 = NA \n")
      use_bound = exact
    }
  }
  # 
  
  # cat(psi_soil, ":", exact, " ", approx_O2, " ", use_bound, "\n")
  Iabs_bound = uniroot(f = f1, interval = c(use_bound*0.01,use_bound*0.99))$root
  
  # dpsi=seq(exact*0.001,exact*0.99, length.out=200)
  # plot(y=sapply(X = dpsi, FUN = f1), x=dpsi, type="l")
  
  list(exact=exact, Iabs_bound=Iabs_bound, approx_O2 = approx_O2)
}


P <- function(psi, psi50, b){
  0.5^((psi/psi50)^b)
}


Pprime <- function(psi, psi50, b){
  log(0.5)*P(psi, psi50, b)*b*(psi/psi50)^(b-1)/psi50
}


Pprimeprime <- function(psi, psi50, b){
  log(0.5)*b*(psi/psi50)^(b-1)/psi50*Pprime(psi, psi50, b) + 
    log(0.5)*P(psi, psi50, b)/psi50^2*b*(b-1)*(psi/psi50)^(b-2)
}


calc_x_from_dpsi <- function(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost){
  gstar = par_photosynth$gammastar/par_photosynth$patm*1e6
  Km = par_photosynth$kmm/par_photosynth$patm*1e6
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  br = par_photosynth$delta
  y = par_cost$gamma
  
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)
  
  x = (ca^2*((3 - 2*br)*gstar + br*Km)*gsprime
       -ca*(y*abs(psi_soil - dpsi)+0)*(gstar + br*Km) -
         sqrt(ca^2*(y*abs(psi_soil - dpsi)+0)*((-3 + 2*br)*gstar - br*Km)*
                ((-1 + br)*ca + gstar + br*Km)*
             ((ca + 2*gstar)*gsprime-(y*abs(psi_soil - dpsi)+0))))/
    (ca^2*((-1 + br)*(y*abs(psi_soil - dpsi)+0) + ((3 - 2*br)*gstar + br*Km)*
             gsprime)) 
  
  # if(is.na(x)){
  #   x <- 1
  # }else{
    x[x<(gstar + br*Km)/(ca - br*ca)]=(gstar + br*Km)/(ca - br*ca)+1e-12
  # }
  
  x
}


calc_J <- function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  d = par_photosynth$delta
  4*gs*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k)) 
}

calc_gsprime <- function(dpsi, psi_soil, par_plant, par_env){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D*P(psi_soil-dpsi, par_plant$psi50, par_plant$b)
}


dFdx <- function(dpsi, psi_soil, par_photosynth, par_plant, par_env, par_cost){
  gs = calc_gs_phydro(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)#* 1e6/par_photosynth$patm
  
  X =  calc_x_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost)
  
  J = calc_J(gs, X, par_photosynth)
  
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  g = par_photosynth$gammastar/par_photosynth$ca
  
  djmax_dJ = calc_djmax_dJ(J, par_photosynth)
  dJ_dchi = calc_dJ_dchi(gs, X, par_photosynth)
  
  dP_dx = -gs*ca - par_cost$alpha * djmax_dJ * dJ_dchi
  
  # dP_ddpsi = gsprime*ca*(1-X) - par_cost$alpha * calc_djmax_dAjmax(ajmax, par_photosynth) * calc_dAjmax_ddpsi(gsprime, X, par_photosynth) - 2*par_cost$gamma*dpsi #/par_plant$psi50^2
  # # cat(c(dP_dx, dP_ddpsi), "\n")
  # c(dP_dx, dP_ddpsi)
  list(dP_dx=dP_dx, J=J, djmax_dJ=djmax_dJ, dJ_dchi=dJ_dchi)
}


calc_jmax_from_J <- function(J, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  4*p/((4*p/J)^2-1)^(1/2)
}


calc_djmax_dJ <- function(J, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  (4*p)^3/((4*p)^2-J^2)^(3/2)
}


calc_dJ_dchi <- function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca/par_photosynth$patm*1e6
  d = par_photosynth$delta
  # gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
  4*gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) - ((x-g)^2+3*g*(1-g)))/(d*(k + x) + g - x)^2)
  #gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
}


chi_jmax_limited <- function(par_photosynth, par_cost){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  b = par_photosynth$delta
  a = par_cost$alpha
  
  #(g*(1-4*a) + 2*sqrt(3)*sqrt(a*(1-4*a)*(1-g)*g))/(1-4*a)
  # (2*sqrt(-a*(4*a + b - 1)*(2*b^2*g*k + 2*b^2*g + b^2*(-k^2) - b^2*k + 2*b*g^2 - 4*b*g*k - 5*b*g + b*k - 3*g^2 + 3*g)) - 4*a*b*k - 4*a*g + b^2*(-k) - b*g + b*k + g)/((b - 1)*(4*a + b - 1))
  (2*sqrt(-a*(4*a + b - 1)*(-3*g + 2*b*g - b*k)*(-1 + b + g + b*k)) - (4*a + b - 1)*(b*k + g))/((b - 1)*(4*a + b - 1))
}

