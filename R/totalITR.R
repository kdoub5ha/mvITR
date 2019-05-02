#' @title Total value function used for treatment assignment value calculation.  
#' 
#' @description This function is used inside the tree growing functions.  
#' 
#' @param dat dataset being assessed 
#' @param z new (alternative) treatment assignment in the splitting procedure.  
#' @param n0 minimum number of observations needed to make a split. 
#' @param aug logical indicator for AIPWE estimation (TRUE) or IPWE estimation (FALSE).
#' @return ITR value from the new treatment assignments
#' @examples 
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # Assign every other patient to treatment
#' z <- rep(c(0,1), nrow(dat)/2)
#' # IPWE value
#' itrtest(dat, z, 5, FALSE)
#' # AIPWE value
#' itrtest(dat, z, 5, TRUE)
#' @export



totalITR <- function(dat, z, n0){
  y <- dat$y.time - min(dat$y.time)
  ae <- dat$y.ae - min(dat$y.ae)
  trt <- dat$trt
  prtx <- dat$prtx
  status <- dat$status
  KM.cens <- dat$KM.cens
  
  itr <- NA
  n <- nrow(dat)
  
  if (length(y)!=length(z)) stop("the vector z must have the same length as data.")
  
    if(length(trt == 1) > 0 & length(trt == 0) > 0){
    
      ae.levels <- unique(sort(dat$ae))
      
      itr.out <- NULL
      for(i in 1:length(ae.levels)){
        index <- which(dat$ae == ae.levels[i])
        itr.ae.num <- sum((trt*ae*z/prtx+(1-trt)*ae*(1-z)/(1-prtx))[index])
        itr.ae.den <- sum((trt*z/prtx+(1-trt)*(1-z)/(1-prtx))[index])
        itr.ae <- itr.ae.num/itr.ae.den
        
        itr.time.num <- sum((trt*y*z/prtx+(1-trt)*y*(1-z)/(1-prtx))[index])
        itr.time.den <- sum((trt*z/prtx+(1-trt)*(1-z)/(1-prtx))[index])
        itr.time <- itr.time.num/itr.time.den
        
        itr.out <- c(itr.out, itr.ae*itr.time)
      }
      
      itr <- sum(itr.out, na.rm = TRUE)
    
    # # ===========================================================
    # # Augmentation not supported yet
    # # ===========================================================
    # 
    # if(aug){
    #   stop("Augmentation not supported yet")
    #   #   first <- (trt*z+(1-trt)*(1-z))*y/(prtx*trt+(1-prtx)*(1-trt))
    #   #   second <- (((trt*z+(1-trt)*(1-z)) - (prtx*trt+(1-prtx)*(1-trt)))/((prtx*trt+(1-prtx)*(1-trt))))*(m1*z+m0*(1-z))
    #   #   itr <- mean(first-second)
    # }
  }
  return(round(itr,4))
}