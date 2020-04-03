#' @title Grows rcDT model
#' @description Main function in mvITR package. Constructs an rcDT model which maximizes an efficacy outcome (y) subject to #'              a risk constraint (r). 
#' @param data data.frame. Data used to construct rcDT model.  
#' Must contain efficacy variable (y), 
#' risk variable (r), 
#' binary treatment indicator coded as 0 / 1 (trt), 
#' propensity score (prtx),
#' candidate splitting covariates.
#' @param split.var numeric vector. Columns of spliting variables.
#' @param efficacy char. Efficacy outcome column. Defaults to 'y'.
#' @param risk char. Risk outcome column. Defaults to 'r'.
#' @param col.trt char. Treatment indicator column name. Should be of form 0/1 or -1/+1.
#' @param col.prtx char. Propensity score column name. 
#' @param risk.control logical. Should risk be controlled? Defaults to TRUE.
#' @param risk.threshold numeric. Desired level of risk control. 
#' @param lambda numeric. Penalty parameter for risk scores. Defaults to 0, i.e. no constraint.
#' 
#' Optional arguments
#' @param test data.frame of testing observations. Should be formatted the same as 'data'.
#' @param min.ndsz numeric specifying minimum number of observations required to call a node terminal. Defaults to 20.
#' @param n0 numeric specifying minimum number of treatment/control observations needed in a split to declare a node terminal. Defaults to 5. 
#' @param max.depth numeric specifying maximum depth of the tree. Defaults to 15 levels. 
#' @param mtry numeric specifying the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param stabilize logical indicating if efficacy should be modeled using residuals. Defaults to TRUE. 
#' @param stabilize.type character specifying method used for estimating residuals. Current options are 'linear' for linear model (default) and 'rf' for random forest. 
#' @param use.other.nodes logical. Should global estimator of objective function be used. Defaults to TRUE. 
#' @param ctg numeric vector corresponding to the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. Not available yet. 
#' @param extremeRandomized logical. Experimental for randomly selecting cutpoints in a random forest model. Defaults to FALSE and users should change this at their own peril. 
#' @param print.summary logical. Should a summary of the tree building be printed? Defaults to TRUE for single trees.
#' @return Summary of a single interaction tree. Each `node` begins with "0" indicating the root node, 
#' followed by a "1" or "2" indicating the less than (or left) child node or greater than (or right) child node. 
#' Additionally, the number of observations `size`, number treated `n.1`, number on control `n.0`, and treatment effect `trt.effect`
#' summaries are provided.  The splitting information includes the column of the chosen splitting variable `var`, the variable name 'vname',
#' the direction the treatment is sent `cut.1` ("r" for right child node, and "l" for left), the chosen split value `cut.2`, 
#' and the estimated value function `score`.
#' @import randomForest
#' @export
#' @examples
#' dat <- generateData()
#' # Generates tree using simualated EMR data with splitting variables located in columns 1-4.
#' tree <- grow.ITR(data = dat, split.var = 1:10, 
#'                  risk.control = TRUE, risk.threshold = 2.75, 
#'                  lambda = 1)


grow.ITR <- function(data, 
                     split.var, 
                     test = NULL, 
                     ctg = NULL, 
                     efficacy = "y",
                     risk = "r",
                     col.trt = "trt",
                     col.prtx = "prtx",
                     risk.control = TRUE, 
                     risk.threshold = NA,
                     lambda = 0, 
                     min.ndsz = 20,
                     n0 = 5, 
                     stabilize = TRUE, 
                     stabilize.type = c("linear", "rf"), 
                     use.other.nodes = TRUE, 
                     mtry = length(split.var), 
                     max.depth = 15, 
                     AIPWE = FALSE, 
                     extremeRandomized = FALSE,
                     print.summary = TRUE)
{
  
  # input checks
  if(!is.data.frame(data)) stop("data argument must be dataframe")
  if(!is.numeric(split.var)) stop("split.var must be numeric vector")

  if(!is.character(efficacy)) stop("efficacy argument must be character")
  if(!efficacy %in% colnames(data)) stop("efficacy argument is not in data")
  if(!is.character(risk)) stop("risk argument must be character")
  if(!risk %in% colnames(data)) stop("risk argument is not in data")

  if(risk.control & is.na(risk.threshold)) warning("risk.contrl is TRUE, but risk.threshold not specified")
  if(lambda < 0) stop("lambda must be > 0")
  
  
  if(!any(stabilize.type %in% c("linear", "rf"))) stop("linear and rf values supported for stabilize.type")
  if(mtry > length(split.var)){
    warning("mtry is larger than split.var length -- setting mtry to length(split.var)")
    mtry <- length(split.var)
  } 
  
  if(print.summary) 
    sprintf("rcDT model risk control specified as %s; Penalty specified as %s", risk.threshold, lambda)
  
  # initialize variables and libraries.
  out <- NULL
  list.nd <- NULL
  list.test <- NULL
  temp.list <- NULL
  temp.test <- NULL
  temp.name <- NULL
  name <- "0"
  full.set <- data
  max.score <- NULL

  # Fill in extra variables if not provided (just as placeholders)
  # These variables will be utilized later as survival endpoints 
  #   are included in the available analyses
  if(is.null(data$KM.cens)) data$KM.cens <- rep(1, nrow(data))
  if(is.null(data$status)) data$status <- rep(1, nrow(data))
  if(is.null(data$id)) data$id <- 1:nrow(data)
  if(is.null(.subset2(data, col.trt))) stop("treatment argument col.trt not found in data")
  if(is.null(.subset2(data, col.prtx))) stop("propensity argument col.prtx not found in data")
  if(col.trt != "trt") data$trt <- data[,col.trt]
  if(col.prtx != "prtx") data$prtx <- data[,col.prtx]

  
  if(sum(data$trt %in% c(0,1)) != nrow(data)){
    data$trt <- ifelse(data$trt == -1, 0, 1)
    warning("Assuming trt indicator is of form -1/+1 and changed values to 0/1")
  }
  
  if(!is.null(test)){
    if(is.null(test$KM.cens)) test$KM.cens <- rep(1, nrow(test))
    if(is.null(test$status)) test$status <- rep(1, nrow(test))
    if(is.null(test$id)) test$id <- 1:nrow(test)
  }  
  
  # Model residuals if requested
  stabilize.type <- match.arg(stabilize.type)
  if(stabilize.type == "rf") require(randomForest)

  if(stabilize){
    dat.tmp <- data.frame(data[ ,split.var], 
                          y = .subset2(data ,efficacy))
    if(stabilize.type == "rf"){
      fit <- randomForest(y ~ ., data = dat.tmp)
      data$y <- fit$y - fit$predicted
    } else if(stabilize.type == "linear"){
      fit <- lm(y ~ ., data = dat.tmp)
      data$y <- fit$residuals
    }
    
    rm(dat.tmp)
    if(!is.null(test)){
      colnames(test) <- gsub(":x", ".x", colnames(test))
      test$y <- test[,efficacy] - predict(fit, test)
      colnames(test) <- gsub("[.]x", ":x", colnames(test))
    }
    fit.y <- fit
    rm(fit)
  } else{
    data$y <- data[ ,efficacy]
    if(!is.null(test)){
      test$y <- test[, efficacy]
    }
  }
  
  # Define risk variables
  data$r <- .subset2(data ,risk)
  if(!is.null(test)){
    test$r <- .subset2(test, risk)
  }
  
  # record total dataset for spliting 
  list.nd <- list(data)
  if (!is.null(test)) list.test <- list(test)
  
  # loop over dataset for spliting 
  while(length(list.nd) != 0) {
    for(ii in sample(1:length(list.nd), size = length(list.nd))){    # select node
      if(!is.null(dim(list.nd[[ii]])) && nrow(list.nd[[ii]]) > 1){   # check that dataset is not degenerate
        
        if (!is.null(test)){         # define test set if requested
          test0 <- list.test[[ii]]
        } else{
          test0 <- NULL
        }
        
        # define temp tree structure for lower level splits
        if(length(list.nd) <= 1) temp.tree <- NULL
        
        # define temporary list of splitting information
        tmp.list <- list(y = .subset2(data, 'y'), 
                         prtx = .subset2(data, 'prtx'), 
                         ae = .subset2(data, 'r'),
                         trt = .subset2(data, 'trt'), 
                         KM.cens = .subset2(data, 'KM.cens'), 
                         maxRisk = risk.threshold,
                         status = .subset2(data, 'status'), 
                         n0 = n0, 
                         lambda = lambda)
        
        # Get current value of the objective
        if(name[ii] == "0"){
          max.score <- max(sapply(0:1, function(iii)
            estITR(append(list(z = rep(iii, nrow(data))), tmp.list))))
        } else{
          trt.pred <- predict.ITR(temp.tree, data, split.var)$trt.pred
          max.score <- estITR(append(list(z = trt.pred), tmp.list))
          dat.rest <- data.frame(data[,colnames(data) != "trt.new"], 
                                 trt.new = trt.pred)
          dat.rest <- dat.rest[!dat.rest$id %in% list.nd[[ii]]$id,]
        }
        rm(tmp.list)
        
        # Determine best split across all covariates
        if(!risk.control){
          
          split <- risk.partition.ITR(dat = list.nd[[ii]], 
                                      test = test0, 
                                      name = name[ii], 
                                      min.ndsz = min.ndsz, 
                                      n0 = n0, 
                                      split.var = split.var, 
                                      ctg = ctg,
                                      max.depth = max.depth, 
                                      mtry = mtry, 
                                      dat.rest = dat.rest, 
                                      max.score = max.score, 
                                      AIPWE = AIPWE, 
                                      risk.threshold = risk.threshold, 
                                      use.other.nodes = use.other.nodes, 
                                      extremeRandomized = extremeRandomized,
                                      lambda = 0)
        } else{
          if(is.na(risk.control)) stop("Level for risk must be specified numeric")
          split <- risk.partition.ITR(dat = list.nd[[ii]], 
                                      test = test0, 
                                      name = name[ii], 
                                      min.ndsz = min.ndsz, 
                                      n0 = n0, 
                                      split.var = split.var, 
                                      ctg = ctg,
                                      max.depth = max.depth, 
                                      mtry = mtry, 
                                      dat.rest = dat.rest, 
                                      max.score = max.score, 
                                      AIPWE = AIPWE, 
                                      risk.threshold = risk.threshold, 
                                      use.other.nodes = use.other.nodes, 
                                      extremeRandomized = extremeRandomized,
                                      lambda = lambda)
        }
        
        out <- rbind(out, split$info)
        if(!is.null(nrow(split$left))&&!is.null(nrow(split$right))){
          min.n <- min(nrow(split$left),nrow(split$right))
        }
        if (!is.null(split$left) && min.n >= min.ndsz && is.null(test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } else if (!is.null(split$left) && min.n >= min.ndsz && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }
      }
      temp.tree <- out
    }
    list.nd <- temp.list
    list.test <- temp.test
    name <- temp.name
    temp.tree<-out
    temp.list <- NULL
    temp.test <- NULL
    temp.name <- NULL
  }
  
  out <- out[order(out$node), ]
  out$var[apply(out, 1, function(j) is.na(de(j['node'], out)[1]))] <- NA
  out[is.na(out$var),c('var','vname','cut.1','cut.2','score')] <- NA
  
  if(!is.null(test)){
    out[is.na(out$var), c('score.test')] <- NA
  }
  
  if(!stabilize) fit.y <- NULL
  
  if(!is.null(test)){
    out.tmp <- list(tree = out, y = data$y, y.test = test$y)
    i.nodes <- out.tmp$tree$node[!is.na(out.tmp$tree$var)]
    names(i.nodes) <- i.nodes
    test.value <- sapply(i.nodes, function(jjj){
      keep <- c(setdiff(out.tmp$tree$node, de(jjj, out.tmp$tree)), paste0(jjj, 1:2))
      tmp.tre <- out.tmp$tree[out.tmp$tree$node %in% keep, ]
      tmp.tre[tmp.tre$node %in% paste0(jjj, 1:2), c('var','vname','cut.1','cut.2','score', 'score.test')] <- NA
      preds <- predict.ITR(tmp.tre, test, split.var, ctgs = ctg)$trt.pred
      score.test <- estITR(list(y = .subset2(test, 'y'), 
                                trt = .subset2(test, 'trt'), 
                                ae = .subset2(test, 'r'),
                                prtx = .subset2(test, 'prtx'), 
                                status = .subset2(test, 'status'), 
                                KM.cens = .subset2(test, 'KM.cens'), 
                                n0 = 5, z = preds, 
                                lambda = lambda, maxRisk = risk.threshold))
    })
    out$score.test <- test.value[match(out$node, names(test.value))]
    setNames(list(out, data$y, test$y, risk.threshold, data, test, fit.y, split.var),
             c("tree", "y", "y.test", "risk.threshold", "data", "test", "fit.y", "split.var"))
  } else{
    setNames(list(out, data$y, risk.threshold, data, fit.y, split.var),
             c("tree", "y", "risk.threshold", "data", "fit.y", "split.var"))
  }
}
