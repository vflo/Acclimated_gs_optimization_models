#Some code comes from Jaideep Joshi PHYDRO repository (https://github.com/jaideep777/phydro) or from RPMODEL

###########################################
## PHOTOSYNTHETIC FUNCTIONS
###########################################

# calculate jmax from temperature
GetPhotosyntheticJmax <- function(jmax25, tem){
  ha=50300.0
  hd=152044.0
  sv=495.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  jmax = jmax25 * factor
  return(jmax)
  }

# calculate vcmax from temperature
GetPhotosyntheticVcmax <- function(vcmax25, tem){
  ha=73637.0
  hd=149252.0
  sv=486.0
  t0=298.15
  r=8.315
  c = 1.0 + exp((sv*t0 -hd)/(r*t0))
  t1 = tem + 273.15
  factor = c * exp(ha/r/t0*(1.0-t0/t1)) / (1.0 + exp((sv*t1-hd)/(r*t1)))
  vcmax = vcmax25 * factor
  return(vcmax)
  }


calc_assimilation_limiting = function(vcmax, jmax, gs, par_photosynth){
  # gs = calc_gs(dpsi, psi_soil, par_plant = par_plant, par_env = par_env)
  
  # We need not employ numerical root-finding. calculate chi independently assuming Ac and Aj, and bigger of the two will be the limiting one. Accordingly return Ac or Aj
  Ac = calc_assim_rubisco_limited(gs = gs, vcmax = vcmax, par_photosynth = par_photosynth)
  Aj = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  
  A = list(ac = Ac$a, aj=Aj$a)
  
  if (Ac$ci > Aj$ci ){
    A$ci = Ac$ci
    A$a  = Ac$a
  } else {
    A$ci = Aj$ci
    A$a  = Aj$a
  }
  
  A
}


calc_assim_rubisco_limited <- function(gs, vcmax, par_photosynth){
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  d = par_photosynth$delta
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * par_photosynth$kmm - vcmax*(1-d)
  C <- gs * ca * par_photosynth$kmm + vcmax * (par_photosynth$gammastar + par_photosynth$kmm*d)
  
  ci <- QUADM(A, B, C)
  # a_c <- vcmax * (ci - par$gammastar) / (ci + par$kmm)
  a_c <- gs*(ca-ci) 
  
  return(list(a=a_c, ci=ci))
}


calc_assim_light_limited <- function(gs, jmax, par_photosynth){
  
  ## Only light is limiting
  ## Solve Eq. system
  ## A = gs (ca- ci)
  ## A = phi0 * Iabs * jlim * (ci - gammastar)/(ci + 2*gamma_star)
  
  ## This leads to a quadratic equation:
  ## A * ci^2 + B * ci + C  = 0
  ## 0 = a + b*x + c*x^2
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  phi0iabs = par_photosynth$phi0 * par_photosynth$Iabs
  jlim = phi0iabs / sqrt(1+ (4*phi0iabs/jmax)^2)
  
  d = par_photosynth$delta 
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * 2 * par_photosynth$gammastar - jlim*(1-d)
  C <- gs * ca * 2*par_photosynth$gammastar + jlim * (par_photosynth$gammastar + d*par_photosynth$kmm)
  
  ci <- QUADM(A, B, C)
  aj <- gs*(ca-ci)
  # vcmax_pot <- a*(ci + par$kmm)/(ci - par$gammastar)
  
  return(list(a=aj, ci=ci))
}


calc_vcmax_coordinated_numerical <-  function(aj, ci, par_photosynth){
  d = par_photosynth$delta
  vcmax_coord = aj*(ci + par_photosynth$kmm)/(ci*(1-d) - (par_photosynth$gammastar+par_photosynth$kmm*d))
  return(vcmax_coord)
}


calc_vcmax_no_acclimated_ww <- function(A, ci, tc, patm, rdark = 0.02){
  gammastar <- calc_gammastar(tc,patm)
  rd <- calc_ftemp_inst_rd(tc)*rdark
  kmm <- calc_kmm(tc,patm)
  vcmax <- (A+rd) * (ci+kmm)/(ci-gammastar) # from eqn no 2.20 de von Caemmerer - Biochemical models of leaf photosynthesis
  # vcmax <- A*(ci + kmm)/(ci*(1-rd) - (gammastar+kmm*rd)) #calculate vcmax using Ac formulation
  return(vcmax)
}

calc_jmax_no_acclimated_ww <- function(A, vcmax, ci, I, tc, patm, kphio = 0.087){
  gammastar <- calc_gammastar(tc,patm)
  # rd <- calc_ftemp_inst_rd(tc)*rdark
  phi0 <- calc_ftemp_kphio(tc)*kphio
  kmm <- calc_kmm(tc,patm)
  # j <- 4*vcmax*(ci+2*gammastar)/(ci+kmm)#calculate j using Aj calculation based in photosynthetic-coordination hypothesis
  # jmax <- 4*phi0*I/(sqrt((4*phi0*I/j)^2-1)) #calculate jmax as the inversion of the j saturation with I
  m <- (ci-gammastar)/( ci+2*gammastar)
  jmax <- 4*phi0*I* (1/sqrt(  (phi0 * I * m/A)^2 -1 )) #eq. 13  inversion from https://doi.org/10.5194/gmd-13-1545-2020
  
  return(jmax)
}

calc_ftemp_kphio <- function( tc ){
  
  ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
  
  return(ftemp)
}

calc_gammastar <- function( tc, patm ) {
  #-----------------------------------------------------------------------
  # Input:    float, air temperature, degrees C (tc)
  # Output:   float, gamma-star, Pa (gammastar)
  # Features: Returns the temperature-dependent photorespiratory
  #           compensation point, Gamma star (Pascals), based on constants
  #           derived from Bernacchi et al. (2001) study.
  # Ref:      Bernacchi et al. (2001), Improved temperature response
  #           functions for models of Rubisco-limited photosynthesis,
  #           Plant, Cell and Environment, 24, 253--259.
  #-----------------------------------------------------------------------
  dha    <- 37830       # (J/mol) Activation energy, Bernacchi et al. (2001)
  gs25_0 <- 4.332       # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  
  gammastar <- gs25_0 * patm / calc_patm(0.0) * calc_ftemp_arrh( (tc + 273.15), dha=dha )
  
  return( gammastar )
}

calc_kmm <- function( tc, patm ) {
  
  dhac   <- 79430      # (J/mol) Activation energy, Bernacchi et al. (2001)
  dhao   <- 36380      # (J/mol) Activation energy, Bernacchi et al. (2001)
  kco    <- 2.09476e5  # (ppm) O2 partial pressure, Standard Atmosphere
  
  ## k25 parameters are not dependent on atmospheric pressure
  kc25 <- 39.97   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  ko25 <- 27480   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  
  ## conversion to Kelvin
  tk <- tc + 273.15
  
  kc <- kc25 * calc_ftemp_arrh( tk, dha=dhac )
  ko <- ko25 * calc_ftemp_arrh( tk, dha=dhao )
  
  po  <- kco * (1e-6) * patm         # O2 partial pressure
  kmm <- kc * (1.0 + po/ko)
  
  return(kmm)
}

calc_ftemp_inst_rd <- function( tc ){
  
  # loal parameters
  apar <- 0.1012
  bpar <- 0.0005
  
  fr <- exp( apar * (tc - 25.0) - bpar * (tc^2 - 25.0^2) )
  
  return(fr)
}

calc_ftemp_arrh <- function( tk, dha, tkref = 298.15 ){
  
  # Note that the following forms are equivalent:
  # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
  # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
  # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
  #-----------------------------------------------------------------------
  kR   <- 8.3145     # Universal gas constant, J/mol/K
  ftemp <- exp( dha * (tk - tkref) / (tkref * kR * tk) )
  
  return(ftemp)
}
