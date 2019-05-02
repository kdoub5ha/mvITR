#' @title Simulates data from an RCT design
#' @description Simulates data from an RCT according to the following model:
#' 2 + 2*sign(x1<cut2) + beta1*trt*subgrp + beta2*(1-trt)*subgrp + N(0,sigma)
#' If depth=1, then subgrp=(x1<=cut1)
#' If depth!=1 then subgrp=(x1>=0.3 & x3>=0.1)
#' 
#' @param n size of the dataset to be generate.  Defaults to 500.
#' @param depth1,depth2 gives the number of interacting covariates. If set to 1, then 
#'  then covariate X1 interacts with treatment. If set to another value, then 
#'  covariates X1 and X3 both interact with treatment effect (one-way interactions). Defaults to 1. 
#' @param design Indicates the design of the study (rct or emr)
#' @param beta11 controls the strength of the treatment effect for survival outcome in treatment group. Defaults to 0.5. 
#' @param beta12 controls the strength of the control effect for survival outcome in control group. Defaults to 0.5. 
#' @param beta21 controls the strength of the treatment effect for AE outcome in treatment group. Defaults to 0.5. 
#' @param beta22 controls the strength of the control effect for AE outcome in control group. Defaults to 0.5.  
#' @param K internal variable used to generate the fineness of the unit interval on which covariates are generated. Defaults to 50.
#' @param K.E internal variable used to generate the fineness of the unit interval on which noise covariates are generated. Defaults to 10.
#' @param sigma controls standard deviation of random variation.  Defaults to 1. 
#' @param rateC censoring rate. Defaults to 0.1. 
#' @param StudyPeriod maximum observation period. Defaults to 36. Not currently in use. 
#' @param n.noise.cov number of noise covariates to generate. Defaults to 10. 
#' @param cut1,cut2,cut3,cut4 controls where the cutpoints are to define subgroups. 
#' @param offset1 The multiplier of the time scale in survival estimates
#' @param lambda Scale parameter for the survival times
#' @param rho Shape parameter for the survival times
#' @param dispersion Dispersion parameter (size = dispersion) for NB dist
#' @return data frame containing y (outcome), x1-x4 (covariates), trt (treatment), prtx (probability of being in treatment group), and ID variable
#' @export
#' @examples
#' 
#' # This generates a dataframe with 500 observations, X1 as the only variable interacting with 
#' # the treatment, and a signal to noise ratio of 2/2=1.
#' data<-rdat(n=500)

mvdat3 <- function(n = 500, ## The sample size
                   K = 50,  ## The thinness of the uniform generator (1/K)
                   K.E = 10, ## The thinness of the uniform generator for excess covariates
                   design = c("rct", "emr"), ## Design of the study
                   beta1 = 1, ## Effect for T == 1 and in subgroup
                   beta2 = 0.25, ## Scaling effect for T == 1 and in subgroup
                   rateC = 0.1,  ## Censoring rate
                   StudyPeriod = 36,  ## Maximum follow up time (not currently in use)
                   n.noise.cov = 10,  ## Number of noise covariates to generate
                   cut1 = 0.5,  ## Controls the subgroup
                   cut2 = 0.3,  ## Controls the subgroup
                   cut3 = 0.1,  ## Controls the subgroup
                   cut4 = 0.7,  ## Controls the subgroup
                   depth1 = 1,  ## Controls the number of interacting covariates in subgroup
                   depth2 = 1,  ## Controls the number of interacting covariates in subgroup
                   offset1 = 50, ## The multiplier of the time scale in survival estimates
                   lambda = 10, ## Scale parameter for the survival times
                   rho = 3, ## Shape parameter for the survival times
                   dispersion = 20, ## Dispersion parameter (size = dispersion) for NB dist
                   sigma = 0.5) ## Standard deviation for the errors in log survival time
{
  
  require(survival)
  options(stringsAsFactors = FALSE)
  design <- match.arg(design)
  #### Generate Covariates
  
  X <- matrix(sample(1:K, 4*n, replace = TRUE) / K, 
              n, 4)
  colnames(X) <- paste0("x", 1:4)
  
  E <- matrix(sample(1:K.E, n*n.noise.cov, replace = TRUE)/K.E, 
              n, n.noise.cov)
  colnames(E) <- paste0("E", 1:n.noise.cov)
  
  if(design == "rct"){
    trt <- sample(c(0,1), n, replace = TRUE)
    prtx <- rep(0.5, n)
  } else{
    stopifnot(design == "emr")
    expLogit  <- exp(-2 + 2*X[,1] + 2*X[,3])
    prtx  <- round(expLogit/(1+expLogit), 3)
    trt  <- rbinom(n, 1, prtx)
  }
  
  ### 
  if(depth1 == 1){
    subgroup1 <- (X[,1] - cut1)
  }else if(depth1 == 2){ 
    subgroup1 <- (X[,1] + X[,3] - 1)
  }
  
  if(depth2 == 1){
    if(depth1 == 1){
      subgroup21 <- (cut2 - X[,4])*(X[,1] <= cut1)
      subgroup22 <- (cut4 - X[,3])*(X[,1] > cut1)
    } else if(depth1 == 2){
      subgroup21 <- (X[,4] - cut1)*(X[,1] + X[,3] <= 1)
      subgroup22 <- (cut2 - X[,2])*(X[,1] + X[,3] > 1)
    }
  }else if(depth2 == 2){ 
    if(depth1 == 1){
      subgroup21 <- (X[,2] + X[,4] - 1)*(X[,1] > cut1)
      subgroup22 <- (1 - X[,2] - X[,4])*(X[,1] <= cut1)
    } else if(depth1 == 2){
      subgroup21 <- (X[,2] + X[,4] - cut4)*(subgroup1 <= 1)
      subgroup22 <- (1 - X[,2] - X[,4])*(subgroup1 > 1)
    }
  }
  
  # ========================================
  # Generate Survival Times
  # ========================================
  
  if(depth1 == 1){
    design.mat1 <- as.matrix(cbind(1, X, E, subgroup1*(2*trt-1)))
    betas1 <- c(2, 0.25, -0.25, rep(0, 2 + n.noise.cov), beta1)
  } else if(depth1 == 2){
    design.mat1 <- as.matrix(cbind(1, X, E, subgroup1*(2*trt-1)))
    betas1 <- c(2, 0.25, -0.25, rep(0, 2 + n.noise.cov), beta1)
  }
  status <- rbinom(n, 1, 1-rateC)
  time <- as.numeric((-log(runif(n, 0, 1)) / (lambda*exp(betas1 %*% t(design.mat1))))^(1 / rho))
  
  # Obtain KM estimate of survival
  df.tmp <- data.frame(time, status, X, E)
  surv <- survfit(coxph(Surv(time, status)~ ., data = df.tmp), data = df.tmp)
  s.estimate <- surv$surv
  s.estimate <- round(ifelse(s.estimate == 1, 1-1/n, 
                             ifelse(s.estimate == 0, 1/n, s.estimate)), 4)
  
  # Obtain KM estimate of censoring
  dat.tmp <- data.frame(time, status, trt, X, E)
  surv.c <- survfit(coxph(Surv(time, 1-status) ~ ., data = dat.tmp))
  c.estimate <- surv.c$surv
  c.estimate <- round(ifelse(c.estimate == 1, 1-1/n, 
                             ifelse(c.estimate == 0, 1/n, c.estimate)), 4)
  
  # Re-order survival estimates to be in proper order
  index <- match(time, surv$time)
  s.estimate <- s.estimate[index]
  
  index <- match(time, surv.c$time)
  c.estimate <- c.estimate[index]
  
  
  # ===============================================
  # Generate Adverse Event Outcomes
  #  from Poisson process (lambda = f(surv.time))
  # ===============================================
  
  design.mat2 <- as.matrix(cbind(1, log(time), X, E, 
                                 subgroup21*(2*trt-1), 
                                 subgroup22*(2*trt-1)))
  betas2 <- c(0.01, 1, -0.01, 0.01, rep(0, 2 + n.noise.cov), 
              beta2, beta2)
  means2 <- exp(design.mat2 %*% betas2)
  AE <- rnbinom(n, size = dispersion, mu = means2)
  
  output <- data.frame(id = 1:n, 
                       X, 
                       E, 
                       time = time, 
                       status = status, 
                       ae = pmin(AE, 10),
                       KM.surv = s.estimate, 
                       KM.cens = c.estimate,
                       trt = trt, 
                       prtx = prtx
  )
  return(output)
}






