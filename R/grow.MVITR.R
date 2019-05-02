#' @title Grows a large interaction tree
#' 
#' @description This function grows an interaction tree using either the IPWE or AIPWE method (AIPWE=F, T). 
#' 
#' 
#' @param data.in data set from which the tree is to be grown.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param test.in testing data
#' @param split.var.in columns of potential spliting variables. Required input.
#' @param min.ndsz.in minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg.in identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param n0.in minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth.in controls the maximum depth of the tree. Defaults to 15. 
#' @param AIPWE.in logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. 
#' @param stabilize.type.in gives the method used for calculating residuals. Current options are 'rf' for random forest and 'linear' for linear model. 
#' @return Summary of a single interaction tree. Each `node` begins with "0" indicating the root node, 
#' followed by a "1" or "2" indicating the less than (or left) child node or greater than (or right) child node. 
#' Additionally, the number of observations `size`, number treated `n.1`, number on control `n.0`, and treatment effect `trt.effect`
#' summaries are provided.  The splitting information includes the column of the chosen splitting variable `var`, the variable name 'vname',
#' the direction the treatment is sent `cut.1` ("r" for right child node, and "l" for left), the chosen split value `cut.2`, 
#' and the estimated value function `score`.
#' @import randomForest
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # Generates tree using simualated EMR data with splitting variables located in columns 1-4.
#' tree <- grow.ITR(data=dat, split.var=1:4)


grow.MVITR <- function(data.in, 
                       split.var.in, 
                       test.in = NULL, 
                       min.ndsz.in = 50,
                       n0.in = 5, 
                       ctg.in = NULL, 
                       max.depth.in = 15, 
                       mtry.in = length(split.var.in), 
                       AIPWE.in = FALSE, 
                       standardize.ae = FALSE,
                       standardize.surv = FALSE, 
                       stabilize.ae = TRUE,
                       stabilize.surv = TRUE,
                       n.folds = 10,
                       stabilize.type.in = c('linear', 'rf'), 
                       use.other.nodes = TRUE, 
                       reverse.ae.scale = FALSE)
{
  
  stabilize.type.in <- match.arg(stabilize.type.in)
  
  # =====================================================
  # Generate an initial split using each outcome
  #  and evaluate the total value
  # This is used to determine the optimal splitting order
  # =====================================================
  
  tre.base.s <- grow.ITR(data = data.in, outcome = "time",
                         split.var = split.var.in, 
                         min.ndsz = min.ndsz.in,
                         n0 = n0.in, ctg = ctg.in, max.depth = 1, 
                         mtry = mtry.in, stabilize = stabilize.surv, 
                         standardize = TRUE, reverse.ae.scale = reverse.ae.scale,
                         stabilize.type = stabilize.type.in)
  
  dat.temp <- data.in
  tre.base.ae <- grow.ITR(data = data.in, outcome = "ae",
                          split.var = split.var.in, 
                          min.ndsz = min.ndsz.in,
                          n0 = n0.in, ctg = ctg.in, max.depth = 1, 
                          mtry = mtry.in, stabilize = stabilize.ae, 
                          standardize = TRUE, reverse.ae.scale = reverse.ae.scale,
                          stabilize.type = stabilize.type.in)
  
  preds.tmp.s <-  predict.ITR(tre.base.s$tree, data.in, ctgs = ctg.in)$trt.pred
  preds.tmp.ae <- predict.ITR(tre.base.ae$tree, dat.temp, ctgs = ctg.in)$trt.pred
  
  # p.change.s <- s.itrtest(data.frame(data.in, y = data.in$time), preds.tmp.s, 5, FALSE) / 
  #   s.itrtest(data.frame(data.in, y = data.in$time), preds.tmp.ae, 5, FALSE)
  # p.change.ae <- itrtest(data.frame(data.in, y = data.in$ae), preds.tmp.ae, 5, FALSE) / 
  #   itrtest(data.frame(data.in, y = data.in$ae), preds.tmp.s, 5, FALSE)
  
  if(nrow(tre.base.ae$tree) > 1 & nrow(tre.base.s$tree) > 1){
    if(tre.base.ae$tree$score[1] > tre.base.s$tree$score[1]){
      first.var <- "ae"
    } else{
      first.var <- "survival"
    }
  } else if(nrow(tre.base.ae$tree) > 1){
    first.var <- "ae"
  } else if(nrow(tre.base.s$tree) > 1){
    first.var <- "survival"
  } else{
    first.var <- "neither"
  }
  
  # =====================================================
  # Generate tree with survival outcome first
  # =====================================================
  # browser()
  if(first.var == "survival"){
    
    tre.base <- grow.ITR(data = data.in, outcome = "time",
                         split.var = split.var.in, 
                         min.ndsz = min.ndsz.in, reverse.ae.scale = reverse.ae.scale,
                         n0 = n0.in, ctg = ctg.in, max.depth = max.depth.in, 
                         mtry = mtry.in, stabilize = stabilize.surv, 
                         standardize = standardize.surv, use.other.nodes = use.other.nodes,
                         stabilize.type = stabilize.type.in, scale.ae = FALSE)

    if(sum(!is.na(tre.base$tree$score)) == 1){
      params <- seq(0, 1.25*tre.base$tree$score[!is.na(tre.base$tree$score)], length.out = 100)
    } else{
      params <- seq(0, 1.25*diff(range(tre.base$tree$score, na.rm = T)), length.out = 100)
    }
    prune.base <- treeCV(tre.base, sp.var = split.var.in, data.in, 
                         param = params, N0 = min.ndsz.in, 
                         outcome = "time", n0 = n0.in, ctgs = ctg.in,
                         stabilize.type = stabilize.type.in, 
                         stabilize = stabilize.surv, sort = TRUE,
                         standardize = standardize.surv, 
                         nfolds = n.folds, scale.ae = FALSE, 
                         use.other.nodes = use.other.nodes,
                         reverse.ae.scale = reverse.ae.scale)
    
    rm(params)
    best.tre.base <- prune.base$best.tree
    best.tre.base <- data.frame(best.tre.base, stringsAsFactors = FALSE)
    if(nrow(best.tre.base) == 1){
      stop("No model identified")
    }

    preds1 <- predict.ITR(best.tre.base, data.in, ctgs = ctg.in)$trt.pred
    dat.tmp <- data.in
    dat.tmp$y <- data.in$time
    value.surv.1 <- s.itrtest(dat.tmp, preds1, 5, FALSE)
    
    nodes <- send.down(data.in, prune.base$best.tree, ctgs = ctg.in)$data$node
    u.nodes <- sort(unique(nodes))
    
    outcome1.nodes <- best.tre.base$node[!best.tre.base$node %in% u.nodes] 

    value.ae.1 <- NULL
    for(t in 1:length(u.nodes)){
      dat.temp <- data.in[which(nodes == u.nodes[t]),]
      # dat.temp$ae <- max(dat.temp$ae) - dat.temp$ae + min(dat.temp$ae)
      temp.tre <- grow.ITR(dat.temp, split.var = split.var.in, 
                           stabilize = stabilize.ae,
                           stabilize.type = stabilize.type.in, outcome = "ae", 
                           standardize = standardize.ae, ctg = ctg.in, 
                           min.ndsz = min.ndsz.in/2, n0 = n0.in/2,  
                           reverse.ae.scale = reverse.ae.scale,
                           use.other.nodes = use.other.nodes, scale.ae = FALSE)
      
      if(nrow(temp.tre$tree) > 1){
        if(sum(!is.na(temp.tre$tree$score)) == 1){
          params <- seq(0, 2*temp.tre$tree$score[!is.na(temp.tre$tree$score)], length.out = 100)
        } else{
          params <- seq(0, 2*diff(range(temp.tre$tree$score, na.rm = T)), length.out = 100)
        }
        # browser()
        prune.temp <- treeCV(temp.tre, dat.temp, 
                             param = params, 
                             outcome = "ae", scale.ae = FALSE, 
                             sp.var = split.var.in, stabilize.type = stabilize.type.in, 
                             stabilize = stabilize.ae, ctgs = ctg.in, 
                             standardize = standardize.ae, sort = TRUE,
                             N0 = max(min.ndsz.in/10, 20), n0 = max(n0.in/20, 5), nfolds = n.folds, 
                             reverse.ae.scale = reverse.ae.scale, use.other.nodes = use.other.nodes)
        temp.best.tre <- prune.temp$best.tree
        
        # Reset node values for merging with parent structure
        temp.best.tre$node <- gsub("0", u.nodes[t], temp.best.tre$node)
        best.tre.base[which(best.tre.base$node == u.nodes[t]), ] <- temp.best.tre[which(temp.best.tre$node == u.nodes[t]), ]
        best.tre.base <- rbind(best.tre.base, 
                               temp.best.tre[-which(temp.best.tre$node == u.nodes[t]),])
      }
      
      preds <- predict.ITR(best.tre.base, dat.temp, ctgs = ctg.in)$trt.pred
      dat.temp$y <- dat.temp$ae
      value.ae.1 <- c(value.ae.1, itrtest(dat.temp, preds, 5, FALSE)*nrow(dat.temp))
    }
    
    out1 <- best.tre.base
    preds.out1 <- predict.ITR(out1, data.in, ctgs = ctg.in)$trt.pred
    out1$outcome <- ifelse(out1$node %in% outcome1.nodes, "Survival", "AdverseEvent")
    
    data.in$y.time <- tre.base.s$y
    data.in$y.ae <- tre.base.ae$y
    value.out1 <- totalITR(data.in, preds.out1, 5)
    
  }
  # =====================================================
  # Generate tree with adverse event outcome first
  # =====================================================
  
  if(first.var == "ae"){
    
    tre.base <- grow.ITR(data = data.in, outcome = "ae",
                         split.var = split.var.in, 
                         min.ndsz = min.ndsz.in, scale.ae = FALSE,
                         n0 = n0.in, ctg = ctg.in, max.depth = max.depth.in, 
                         mtry = mtry.in, stabilize = stabilize.ae, 
                         standardize = standardize.ae, reverse.ae.scale = reverse.ae.scale,
                         stabilize.type = stabilize.type.in, use.other.nodes = use.other.nodes)
    
    if(sum(!is.na(tre.base$tree$score)) == 1){
      params <- seq(0, 1.25*tre.base$tree$score[!is.na(tre.base$tree$score)], length.out = 100)
    } else{
      params <- seq(0, 1.25*diff(range(tre.base$tree$score, na.rm = T)), length.out = 100)
    }
    prune.base <- treeCV(tre.base, data.in,                        
                         param = params, scale.ae = FALSE,
                         stabilize = stabilize.ae,
                         sp.var = split.var.in, N0 = min.ndsz.in/10, 
                         outcome = "ae", n0 = n0.in/5, ctgs = ctg.in,
                         stabilize.type = stabilize.type.in, 
                         standardize = standardize.ae, nfolds = n.folds, 
                         reverse.ae.scale = reverse.ae.scale, use.other.nodes = use.other.nodes)
    best.tre.base <- prune.base$best.tree
    best.tre.base <- data.frame(best.tre.base, stringsAsFactors = FALSE)
    
    if(nrow(best.tre.base) == 1){
      stop("No model identified")
    }
    
    preds1 <- predict.ITR(best.tre.base, data.in, ctgs = ctg.in)$trt.pred
    dat.tmp <- data.in
    dat.tmp$y <- dat.tmp$ae
    value.ae.1 <- itrtest(dat.tmp, preds1, 5, FALSE)
    
    nodes <- send.down(data.in, prune.base$best.tree, ctgs = ctg.in)$data$node
    u.nodes <- sort(unique(nodes))
    
    outcome1.nodes <- best.tre.base$node[!best.tre.base$node %in% u.nodes] 
    
    value.survival.1 <- NULL
    for(t in 1:length(u.nodes)){
      dat.temp <- data.in[which(nodes == u.nodes[t]),]
      temp.tre <- grow.ITR(dat.temp, split.var = split.var.in, stabilize = stabilize.surv, 
                           stabilize.type = stabilize.type.in, outcome = "time", 
                           ctg = ctg.in, standardize = standardize.surv, 
                           min.ndsz = min.ndsz.in/2, n0 = n0.in/2,  
                           reverse.ae.scale = reverse.ae.scale,
                           use.other.nodes = use.other.nodes)
      if(nrow(temp.tre$tree) > 1){
        if(sum(!is.na(temp.tre$tree$score)) == 1){
          params <- seq(0, 2*temp.tre$tree$score[!is.na(temp.tre$tree$score)], length.out = 100)
        } else{
          params <- seq(0, 2*diff(range(temp.tre$tree$score, na.rm = T)), length.out = 100)
        }
        prune.temp <- treeCV(temp.tre, dat.temp, 
                             param = params, 
                             outcome = "time", scale.ae = FALSE, 
                             sp.var = split.var.in, stabilize.type = stabilize.type.in, 
                             stabilize = stabilize.ae, ctgs = ctg.in, 
                             standardize = standardize.ae, sort = TRUE,
                             N0 = min.ndsz.in/5, n0 = n0.in/5, nfolds = n.folds, 
                             reverse.ae.scale = reverse.ae.scale, use.other.nodes = use.other.nodes)
        temp.best.tre <- prune.temp$best.tree
        
        # Reset node values for merging with parent structure
        temp.best.tre$node <- gsub("0", u.nodes[t], temp.best.tre$node)
        best.tre.base[which(best.tre.base$node == u.nodes[t]), ] <- temp.best.tre[which(temp.best.tre$node == u.nodes[t]), ]
        best.tre.base <- rbind(best.tre.base, 
                               temp.best.tre[-which(temp.best.tre$node == u.nodes[t]),])
      }
      preds <- predict.ITR(best.tre.base, dat.temp, ctgs = ctg.in)$trt.pred
      dat.temp$y <- dat.temp$time
      value.survival.1 <- c(value.survival.1, s.itrtest(dat.temp, preds, 5, FALSE)*nrow(dat.temp))
    }
    
    out1 <- best.tre.base
    preds.out1 <- predict.ITR(out1, data.in, ctgs = ctg.in)$trt.pred
    out1$outcome <- ifelse(out1$node %in% outcome1.nodes, "AdverseEvents", "Survival")
    
    data.in$y.time <- tre.base.s$y
    data.in$y.ae <- tre.base.ae$y
    value.out1 <- totalITR(data.in, preds.out1, 5)
    
  }
  
  if(first.var %in% c("ae", "survival")){
    best.model <- out1
    best.value <- value.out1
    
    out <- list()
    out$best.model <- best.model
    out$best.model.value <- best.value
  } else{
    out <- list()
    out$best.model <- NA
    out$best.model.value <- NA
  }
  
  out$best.model$outcome[is.na(out$best.model$score)] <- NA
  
  return(out)
}
