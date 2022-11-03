#Some code comes from Jaideep Joshi PHYDRO repository (https://github.com/jaideep777/phydro)

###########################################
## HYDRAULIC FUNCTIONS
###########################################
## Returns Huber value for SAPFLUXNET data. LAI, stand basal area, plant basal area and plant sapwood area must be provided
huber_value <- function(par_env, par_plant){
  par_plant$pl_sapw_area / par_plant$pl_ba * par_plant$st_basal_area * 1e-4 / par_env$LAI
}


## Returns conductivity in mol/m2/s/Mpa ksmax is provided
scale_conductivity_ks = function(K, par_env, par_plant, do_backtransform = TRUE){
  if(do_backtransform){
    #from Kg m-1 s-1 MPa-1 to m3 m-2
    K = (K*par_env$viscosity_water*huber_value(par_env, par_plant))/(par_env$density_water*1e6*par_plant$height) #back transform to m3 m-2
  }
  
  # Flow rate in m3/m2/s/Pa
  K2 = K / par_env$viscosity_water
  
  # Flow rate in mol/m2/s/Pa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2 * par_env$density_water * mol_h20_per_kg_h20
  
  # Flow rate in mol/m2/s/Mpa
  K4 = K3*1e6
  
  return(K4)  
}


## Returns conductivity from sapwood linear permeability (m-1) to leaf surface (m-2)
calc_conductivity_m = function(sapwood_perm, hv, height){
  sapwood_perm*hv/height
}


## Returns conductivity in mol/m2/s/Mpa
scale_conductivity = function(K, par_env){
  # Flow rate in m3/m2/s/Pa
  K2 = K / par_env$viscosity_water
  
  # Flow rate in mol/m2/s/Pa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2 * par_env$density_water * mol_h20_per_kg_h20
  
  # Flow rate in mol/m2/s/Mpa
  K4 = K3*1e6
  
  return(K4)  
}


## Calculate dpsi if gs is provided
calc_dpsi_phydro <- function (gs, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  P50 = par_plant$psi50
  b = par_plant$b
  
  #calculate dpsi using the inverse function of the incomplete gamma function
  l2 = log(2)
  int_P = -gs*1.6*D/K
  int_s = integral_P_soil(psi_soil,P50, b)
  int_l = (int_s + int_P)
  if(int_l<0){int_l = 0}
  inst_l_gamma = int_l/(-(P50/b)*(l2^(-1/b)))
  
  pl_ = zipfR::Igamma.inv(a = 1/b, y = inst_l_gamma,lower = FALSE)
  pl = (pl_/l2)^(1/b)
  dpsi = psi_soil - pl*P50 # dpsi in MPa
  dpsi
}

## Calculate gs if dpsi is provided. 
calc_gs_phydro <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity(par_plant$conductivity, par_env) #mol m-2 (leaf) s-1 MPa-1

  D = (par_env$vpd/par_env$patm)

  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}

## Calculate gs if dpsi  and ksmax Kg m-1 s-1 MPa-1 in is provided. 
calc_gs <- function (dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity_ks(par_plant$conductivity, par_env, par_plant, do_backtransform = TRUE) #mol m-2 (leaf) s-1 MPa-1
  # K = K * par_plant$pl_sapw_area/par_plant$pl_ba*par_plant$st_basal_area*1e-4 #kg m-2 (ground) s-1 MPa-1
  # K = K * par_env$LAI #kg m-2 (leaf) s-1 MPa-1
  D = (par_env$vpd/par_env$patm)
  # D = par_env$vpd
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b) #mol m-2 (leaf) s-1
}

## Return the integral defined by dpsi of the vulnerability curve
integral_P <- function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  pl_ = l2*pl^b
  ps_ = l2*ps^b
  -(psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = pl_)-expint::gammainc(a = 1/b, x = ps_))
}

## Return the integral defined by dpsi of the vulnerability curve divided 
## by the integral from soil psi to -infinite (e_crit). Equivalent to PROFITMAX2 impairment
integral_P_e_ecrit <- function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  pl_ = l2*pl^b
  ps_ = l2*ps^b
  -1*((expint::gammainc(a = 1/b, x = pl_)-expint::gammainc(a = 1/b, x = ps_)))/
    (-1*expint::gammainc(a = 1/b, x = ps_))
}

## Return the integral from soil psi to -infinite (e_crit)
integral_P_ecrit <- function(psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  l2 = log(2)
  ps_ = l2*ps^b
  (psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = ps_))
}


# integral_P_soil <- function(psi_soil, psi50, b, ...){
#   ps = psi_soil/psi50
#   l2 = log(2)
#   ps_ = l2*ps^b
#   
#   -(psi50/b)*(l2^(-1/b))*expint::gammainc(a = 1/b, x = ps_)
# }