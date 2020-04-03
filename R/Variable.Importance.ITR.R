#' @title Calcuates variable importance measures for a random forest object.  
#' 
#' @description This function accepts a forest object from the `Build.RF.ITR` function and estimates the importance of each
#' predictor. This is accomplished by considering each tree in the forest, obtaining the out-of-bag value for each predictor in that tree, 
#' obtaining the permuted out-of-bag value for each predictor in the tree, and comparing the values. A larger discrepancy between 
#' the original value and permuted value indicates the predictor is more important in predicting treatment. The function
#' returns the variables in order of importance along with the importance measure, scaled to be out of 1.  
#' 
#' @param RF.fit forest object from Build.RF.ITR. Required input. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 2. 
#' @param sort sort the variable importance measure? Defaults to TRUE. 
#' @param N0 minimum number of observations needed in a split to call a node terminal. Defaults to 20. 
#' @param details print details of each tree as the function progresses. Defaults to FALSE.
#' @param depth internal variable.
#' @param AIPWE indicator for AIPWE estimation.
#' @return Returns list of total (VI), efficacy (VI.Efficacy), and risk (VI.Risk) 
#'         ordered variable importance measure calculated for each splitting variable. 
#' @import randomForest
#' @examples 
#' set.seed(1)
#' dat <- generateData(n = 500)
#' # Build a forest with 100 trees
#' fit <- Build.RF.ITR(dat = dat, split.var = 1:10, ntree = 200,
#'                     risk.control = TRUE, risk.threshold = 2.75, 
#'                     lambda = 1)
#' VI <- Variable.Importance.ITR(fit)
#' @export


Variable.Importance.ITR <- function(RF.fit, 
                                    n0 = 5, 
                                    N0 = 20, 
                                    sort = TRUE, 
                                    details = FALSE,
                                    depth = 1, 
                                    AIPWE = FALSE){

  trees <- RF.fit$TREES
  id.boots <- RF.fit$ID.Boots.Samples
  # ARGUMENTS FOR MODEL SPECIFICATION 
  Model.Specification <- RF.fit$Model.Specification
  dat0 <- Model.Specification$data
  col.y <- Model.Specification$efficacy
  col.r <- Model.Specification$risk
  col.trt <- Model.Specification$col.trt
  col.prtx <- Model.Specification$col.prtx
  split.var <- Model.Specification$split.var
  risk.threshold <- Model.Specification$risk.threshold
  lambda <- Model.Specification$lambda
  ctg <- Model.Specification$ctg
  vnames <- colnames(dat0)[split.var]
  # 
  ntree <- length(trees)
  p <- length(split.var)
  VI <- rep(0, p)
  VI.Efficacy <- rep(0, p)
  VI.Risk <- rep(0, p)
  
  for (b in 1:ntree){
    id.b <- id.boots[[b]]
    dat.oob <- dat0[-sort(unique(id.b)), ] 
    n.oob <- nrow(dat.oob)	
    tre.b <- trees[[b]]
    ########## NOTE THAT revise.tree=T HERE! ##########
    out0.b.Total <- send.down.VI.ITR(dat.new = dat.oob, 
                                     tre = tre.b, 
                                     col.y = col.y, 
                                     col.r = col.r, 
                                     col.trt = col.trt,
                                     col.prtx = col.prtx,
                                     lambda = lambda,
                                     risk.threshold = risk.threshold,
                                     ctg = ctg, 
                                     n0 = n0, 
                                     N0 = N0,
                                     revise.tree = TRUE,
                                     depth = depth, 
                                     AIPWE = AIPWE)  
    
    out0.b.Efficacy <- send.down.VI.ITR.Efficacy(dat.new = dat.oob, 
                                                 tre = tre.b, 
                                                 col.y = col.y, 
                                                 col.r = col.r, 
                                                 col.trt = col.trt,
                                                 col.prtx = col.prtx,
                                                 lambda = lambda,
                                                 risk.threshold = risk.threshold,
                                                 ctg = ctg, 
                                                 n0 = n0, 
                                                 N0 = N0,
                                                 revise.tree = TRUE,
                                                 depth = depth, 
                                                 AIPWE = AIPWE)  
    
    out0.b.Risk <- send.down.VI.ITR.Risk(dat.new = dat.oob, 
                                         tre = tre.b, 
                                         col.y = col.y, 
                                         col.r = col.r, 
                                         col.trt = col.trt,
                                         col.prtx = col.prtx,
                                         lambda = lambda,
                                         risk.threshold = risk.threshold,
                                         ctg = ctg, 
                                         n0 = n0, 
                                         N0 = N0,
                                         revise.tree = TRUE,
                                         depth = depth, 
                                         AIPWE = AIPWE)  
    
    tre0.b.Total <- out0.b.Total$tre0	
    tre0.b.Efficacy <- out0.b.Efficacy$tre0	
    tre0.b.Risk <- out0.b.Risk$tre0	
    if (nrow(tre0.b.Total) > 0) {					### AVOID NULL TREES	
      Xs.b.Total <- sort(unique(na.omit(tre0.b.Total$var))) 
      Xs.b.Efficacy <- sort(unique(na.omit(tre0.b.Efficacy$var))) 
      Xs.b.Risk <- sort(unique(na.omit(tre0.b.Risk$var))) 
      G.oob.Total <- out0.b.Total$score  # Overall score from original bootstrap tree
      G.oob.Efficacy <- out0.b.Efficacy$score  # Overall efficacy score from original bootstrap tree
      G.oob.Risk <- out0.b.Risk$score  # Overall risk score from original bootstrap tree
      for (j in 1:p) {
        if (details) print(j)
        G.j.Total <- G.oob.Total   # Initialize score for covariate j to be to bootstrap score value
        G.j.Efficacy <- G.oob.Efficacy   # Initialize score for covariate j to be to bootstrap score value
        G.j.Risk <- G.oob.Risk   # Initialize score for covariate j to be to bootstrap score value
        col.xj <- split.var[j] 
        if (is.element(col.xj, Xs.b.Total)){
          x.j <- dat.oob[, col.xj]
          dat.permuted <- dat.oob
          dat.permuted[ , col.xj] <- x.j[sample(1:n.oob,n.oob, replace=FALSE)]
          ########## NOTE THAT revise.tree=F HERE! ##########
          out0.bj.Total <- send.down.VI.ITR(dat.new = dat.permuted, 
                                            tre = tre0.b.Total, 
                                            col.y = col.y,
                                            col.r = col.r, 
                                            col.trt = col.trt, 
                                            col.prtx = col.prtx,
                                            lambda = lambda,
                                            risk.threshold = risk.threshold,
                                            ctg = ctg, 
                                            n0 = n0,
                                            N0 = N0, 
                                            revise.tree = FALSE,
                                            depth = 1,
                                            AIPWE = AIPWE)
        
          out0.bj.Efficacy <- send.down.VI.ITR.Efficacy(dat.new = dat.permuted, 
                                                        tre = tre0.b.Efficacy, 
                                                        col.y = col.y,
                                                        col.r = col.r, 
                                                        col.trt = col.trt, 
                                                        col.prtx = col.prtx, 
                                                        lambda = lambda,
                                                        risk.threshold = risk.threshold,
                                                        ctg = ctg, 
                                                        n0 = n0,
                                                        N0 = N0, 
                                                        revise.tree = FALSE,
                                                        depth = 1,
                                                        AIPWE = AIPWE)
          
          out0.bj.Risk <- send.down.VI.ITR.Risk(dat.new = dat.permuted, 
                                                tre = tre0.b.Risk, 
                                                col.y = col.y,
                                                col.r = col.r, 
                                                col.trt = col.trt, 
                                                col.prtx = col.prtx, 
                                                lambda = lambda,
                                                risk.threshold = risk.threshold,
                                                ctg = ctg, 
                                                n0 = n0,
                                                N0 = N0, 
                                                revise.tree = FALSE,
                                                depth = 1,
                                                AIPWE = AIPWE)
          
          tre0.bj.Total <- out0.bj.Total$tre0		
          tre0.bj.Efficacy <- out0.bj.Efficacy$tre0		
          tre0.bj.Risk <- out0.bj.Risk$tre0		
          G.j.Total <- ifelse(nrow(tre0.bj.Total) == 1, G.oob.Total, out0.bj.Total$score)
          G.j.Efficacy <- ifelse(nrow(tre0.bj.Efficacy) == 1, G.oob.Efficacy, out0.bj.Efficacy$score)
          G.j.Risk <- ifelse(nrow(tre0.bj.Risk) == 1, G.oob.Risk, out0.bj.Risk$score)
        }
        if (G.j.Total > G.oob.Total) G.j.Total <- G.oob.Total
        if (G.j.Efficacy > G.oob.Efficacy) G.j.Efficacy <- G.oob.Efficacy
        if (G.j.Risk < G.oob.Risk) G.j.Risk <- G.oob.Risk
        ##################### PREVENTS NEGATIVE IMPORTANCE VALUES 
        VI[j] <- VI[j] + (G.oob.Total - G.j.Total)/G.oob.Total
        VI.Efficacy[j] <- VI.Efficacy[j] + (G.oob.Efficacy - G.j.Efficacy)/G.oob.Efficacy
        VI.Risk[j] <- VI.Risk[j] + (G.j.Risk - G.oob.Risk)/G.oob.Risk
      }
    }	
  }

  names(VI) <- vnames
  names(VI.Efficacy) <- vnames
  names(VI.Risk) <- vnames
  if (sort){
    VI <- sort(VI, decreasing = TRUE)
    VI.Efficacy <- sort(VI.Efficacy, decreasing = TRUE)
    VI.Risk <- sort(VI.Risk, decreasing = TRUE)
  }  
  VI <- VI/sum(VI)
  VI.Efficacy <- VI.Efficacy/sum(VI.Efficacy)
  VI.Risk <- VI.Risk/sum(VI.Risk)
  
  out <- list(VI = VI, 
              VI.Efficacy = VI.Efficacy,
              VI.Risk = VI.Risk)
  return(out)
}