#' @title Grows a large interaction tree
#' @description This function grows an interaction tree using either the IPWE or AIPWE method (AIPWE=F, T). 
#' @param data data set from which the tree is to be grown.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param test testing data
#' @param split.var columns of potential spliting variables. Required input.
#' @param min.ndsz minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 15. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. 
#' @param in.forest logical for if the tree is being constructed in a forest. Should not be changed from defaults.
#' @param stabilize.type gives the method used for calculating residuals. Current options are 'rf' for random forest and 'linear' for linear model. 
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


grow.ITR <- function(data, 
                     split.var, 
                     test = NULL, 
                     ctg = NULL, 
                     outcome = c("time", "ae"),
                     scale.ae = FALSE, 
                     stabilize.type = c("linear", "rf"), 
                     min.ndsz = 20,
                     n0 = 5, 
                     mtry = length(split.var), 
                     max.depth = 15, 
                     standardize = FALSE,
                     AIPWE = FALSE, 
                     stabilize = TRUE, 
                     haoda.method = FALSE, 
                     haoda.ae.level = NA, 
                     reverse.ae.scale = FALSE, 
                     use.other.nodes = TRUE)
{
  # initialize variables and libraries.
  library(randomForest)
  out <- NULL
  list.nd <- NULL
  list.test <- NULL
  temp.list <- NULL
  temp.test <- NULL
  temp.name <- NULL
  name <- "0"
  full.set <- data
  max.score <- NULL
  outcome <- match.arg(outcome)

  # Scale AE if requested
  if(outcome == "ae" | haoda.method){
    if(scale.ae){
      data$ae <- data$ae / data$time
      if(!is.null(test)) test$ae <- test$ae / test$time
    }
    if(outcome == "ae"){
      if(reverse.ae.scale){
        data$ae <- max(data$ae) - data$ae
        if(!is.null(test)) test$ae <- max(data$ae) - test$ae
      }
    }
  }
  
  # Model residuals if requested
  stabilize.type <- match.arg(stabilize.type)
  if(stabilize){
    dat.tmp <- data.frame(data[,split.var], y = data[,outcome])
    switch(stabilize.type, 
           rf = {fit <- randomForest(y ~ ., data = dat.tmp, maxnodes = 10)
           data$y <- fit$y - predict(fit, dat.tmp)}, 
           linear = {dat.tmp$trt <- data$trt
           fit <- lm(y ~ ., data = dat.tmp)
           data$y <- fit$residuals})
    rm(dat.tmp)
    if(!is.null(test)){
      colnames(test) <- gsub(":x", ".x", colnames(test))
      test$y <- test[,outcome] - predict(fit, test)
      colnames(test) <- gsub("[.]x", ":x", colnames(test))
    }
  }

  # Standardize if requested
  if(standardize){
    data$y <- if(stabilize){
      scale(data$y)
    } else{
      scale(data[,outcome])
    }
    if(!is.null(test)){
      test$y <- if(stabilize){
        scale(test$y, mean(data$y), sd(data$y))
      } else{
        scale(test[,outcome], mean(data[,outcome]), sd(data[,outcome]))
      }
    }
  } else{
    if(!stabilize){
      data$y <- data[,outcome]
      test$y <- test[,outcome]
    }
  }

  # record total dataset for spliting 
  list.nd <- list(data)
  if (!is.null(test)) list.test <- list(test)

  # loop over dataset for spliting 
  while(length(list.nd) != 0) {
    for(ii in sample(1:length(list.nd), size = length(list.nd))){
      if(!is.null(dim(list.nd[[ii]])) && nrow(list.nd[[ii]]) > 1){
        if (!is.null(test)){
          test0 <- list.test[[ii]]
        } else{
          test0 <- NULL
        }
        
        if(length(list.nd) <= 1) temp.tree <- NULL
        
        if(name[ii] == "0"){
          if(!haoda.method){
            max.score <- switch(outcome, 
                                time = max(sapply(0:1, function(t){
                                  s.itrtest(data, rep(t, nrow(data)), n0, AIPWE)})), 
                                ae = max(sapply(0:1, function(t){
                                  itrtest(data, rep(t, nrow(data)), n0, AIPWE)})))
          } else{
            max.score <- 1E-5
          }
        }else{
          
          send<-send.down(data, temp.tree, ctgs = ctg)
          node<-substr(send$data$node,1,nchar(as.character(send$data$node))-1)
          direction<-substr(send$data$node,
                            nchar(as.character(send$data$node)),
                            nchar(as.character(send$data$node)))
          trt.dir <- temp.tree[match(node,temp.tree$node),]$cut.1
          
          trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                           ifelse(trt.dir=="r" & direction=="2",1,
                                  ifelse(trt.dir=="l" & direction=="1",1,0)))
          
          if(!haoda.method){
            max.score <- switch(outcome, 
                                time = s.itrtest(dat = data, z = trt.pred, n0, AIPWE), 
                                ae = itrtest(dat = data, z = trt.pred, n0, AIPWE))
          } else{
            max.score <- itrtest(dat = data, z = trt.pred, n0, AIPWE)
          }
          
          #Obtain treatments for those not included in the node
          dat.rest <- data[!data$id %in% list.nd[[ii]]$id,]
          dat.rest$trt.new <- trt.pred[!data$id %in% list.nd[[ii]]$id]
        }

        # Determine best split across all covariates
        if(!haoda.method){
          split <- partition.ITR(dat = list.nd[[ii]], test = test0, 
                                 name = name[ii], min.ndsz = min.ndsz, 
                                 n0 = n0, split.var = split.var, ctg = ctg,
                                 max.depth = max.depth, mtry = mtry, 
                                 dat.rest = dat.rest, max.score = max.score, 
                                 AIPWE = AIPWE, outcome = outcome, 
                                 use.other.nodes = use.other.nodes)
        } else{
          if(is.na(haoda.ae.level)) stop("Level for risk must be specified numeric")
          split <- haoda.partition.ITR(dat = list.nd[[ii]], test = test0, 
                                       name = name[ii], min.ndsz = min.ndsz, 
                                       n0 = n0, split.var = split.var, ctg = ctg,
                                       max.depth = max.depth, mtry = mtry, 
                                       dat.rest = dat.rest, max.score = max.score, 
                                       AIPWE = AIPWE, outcome = outcome, 
                                       haoda.ae.level = haoda.ae.level, 
                                       use.other.nodes = use.other.nodes)
        }

        out <- rbind(out, split$info)
        if(!is.null(nrow(split$left))&&!is.null(nrow(split$right))){
          min.n <- min(nrow(split$left),nrow(split$right))
        }
        if (!is.null(split$left) && min.n>min.ndsz && is.null(test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } else if (!is.null(split$left) && min.n>min.ndsz && !is.null(test) && !is.null(split$left.test)) {
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
  
  if(!is.null(test)){
    out.tmp <- list(tree = out, y = data$y, y.test = test$y)
    i.nodes <- out.tmp$tree$node[!is.na(out.tmp$tree$var)]
    names(i.nodes) <- i.nodes
    test.value <- sapply(i.nodes, function(jjj){
      keep <- c(setdiff(out.tmp$tree$node, de(jjj, out.tmp$tree)), paste0(jjj, 1:2))
      tmp.tre <- out.tmp$tree[out.tmp$tree$node %in% keep, ]
      tmp.tre[tmp.tre$node %in% paste0(jjj, 1:2), c('var','vname','cut.1','cut.2','score', 'score.test')] <- NA
      preds <- predict.ITR(tmp.tre, test, ctgs = ctg)$trt.pred
      itrtest(test, preds, 5, F)
    })
    out$score.test <- test.value[match(out$node, names(test.value))]
    setNames(list(out, data$y, test$y), c("tree", "y", "y.test"))
  } else{
    setNames(list(out, data$y), c("tree", "y"))
  }
}
