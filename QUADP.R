# solves quadratic equation given by y = a*x^2 + bx + c
# Based on MAESTRA equivalent (B. Medlyn)

#' Quadratic Solver
#'
#' 
#'
#' @param A numeric, Stomatal conductance in mol/m2/s, for e.g. as calculated from \code{\link{calc_gs}} 
#' @param B numeric, Carboxylation capacity (umol/m2/s)
#' @param C numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#'
#' @return root
#'
#' @export
#'
# Calculate Assimilation rate (umol/m2/s), given
# - gs (mol/m2/s)
# - vcmax (umol/m2/s)
# - Photosynthesis parameters with concentrations converted to partial pressures (Pa)

# - larger root
QUADP <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
    }
  }
  
}

#' @export
# - smaller root
QUADM <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B - sqrt(B^2 - 4*A*C)) / (2*A))
    } 
  }
  
}
