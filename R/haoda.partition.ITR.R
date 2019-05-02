#' @title Generates partition summary based on itr value. 
#' 
#' @description The partitioning function is used inside of the tree growing functions.  It 
#' selects the best split among a set of predictors based on the IPWE or AIPWE value. 
#' 
#' @param dat data set from which the partition is to be made.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param test testing data set 
#' @param name internal variable the keeps track of the current node(s). 
#' @param split.var columns of potential spliting variables. Required input.
#' @param min.ndsz minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 15. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param max.score controls the minimum score required to make an additional split (internally controlled).
#' @param dat.rest data outside the current node being split
#' @param AIPWE indicator for AIPWE estimation
#' @param haoda.ae.level cutoff for acceptable risk in a given split (numeric)
#' @return summary of the best split for a given data frame. 
#' @export



haoda.partition.ITR<-function(dat, 
                              test = NULL, 
                              name = "0", 
                              min.ndsz = 20, 
                              n0 = 5, 
                              split.var, 
                              outcome = c('time'),
                              ctg = ctg,
                              max.depth = 15, 
                              mtry = length(split.var), 
                              dat.rest = NULL, 
                              max.score = NULL, 
                              AIPWE = AIPWE,
                              haoda.ae.level = NA, 
                              use.other.nodes = TRUE)
{   
  # inialize various variable
  call <- match.call()
  out <- match.call(expand = F)
  out$info <- NULL
  out$name.l <- NULL
  out$name.r <- NULL
  out$left <- NULL
  out$right <- NULL
  out$... <- NULL
  outcome <- match.arg(outcome)
  # label the binary tree by 1 (left) and 2 (right).
  name.l <- paste(name, 1, sep="")
  name.r <- paste(name, 2, sep="")
  # sample size
  n <- nrow(dat)
  # check whether testing data is provided
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
  }
  # prepare for the first cut these variable were used to store final cut information
  var <- NA
  vname <- NA
  cut <- NA
  
  if(name=="0") {
    dat.comb <- dat[ ,c('y', 'trt', 'prtx', 'ae')]
  } else{
    dat.comb <- rbind(dat.rest[, c('y', 'trt', 'prtx', 'ae')], 
                      dat[,c('y', 'trt', 'prtx', 'ae')])
  }
  
  # extract value from data
  trt <- dat$trt
  y <- dat$y
  vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  trt.effect <- NA
  n.0 = length(y[trt==0])
  n.1 = n - n.0
  if (min(n.1, n.0) >0) {
    trt.effect <- mean(y[trt==1]) - mean(y[trt==0])
  }
  # CONTROL THE MAX TREE DEPTH
  # name is the current tree label.
  # only when currently depth < max.depth and n > min terminal node proceed.
  depth <- nchar(as.character(name)) 
  
  if(depth <= max.depth && n >= min.ndsz) {
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)
    
    # Fulfill splitting conditions
    splitVars <- sample(split.var, size = m.try, replace = FALSE)
    names(splitVars) <- colnames(dat)[splitVars]
    
    # apply search algorithm to each potential split
    outSplit <- lapply(splitVars, function(i){
      x <- .subset2(dat, i)
      v.name <- vnames[i]
      temp <- sort(unique(x))
      is.ctg <- is.element(i, ctg)

      if(length(temp) > 1) {
        # handle categorial variable first, otherwise take out the final value as no cut after it.
        if(is.ctg){
          zcutCat <- power.set(temp)
          zcutCat <- do.call(cbind, lapply(zcutCat, 
                                           function(zzz) as.numeric(x %in% zzz)))
          zcut <- 1:10
        } else{
          zcut <- temp[-length(temp)]
          zcutCat <- matrix(1)
        }
        
        if(name == "0"){
          if(!is.ctg){
            x.tmp <- x
          } else{
            x.tmp <- as.numeric(factor(x))
          }
          datMatrix <- list(y = dat.comb$y, 
                            ae = dat.comb$ae, 
                            x = x.tmp, 
                            prtx = dat.comb$prtx, 
                            trt = dat.comb$trt, 
                            trtNew = rep(-1000, length(x)),
                            inNode = rep(1, length(x)))
        } else{
          if(!is.ctg){
            x.tmp <- c(dat.rest[,i], x)
          } else{
            x.tmp <- as.numeric(factor(c(dat.rest[,i], x)))
          }
          
          datMatrix <- list(y = dat.comb$y, 
                            ae = dat.comb$ae, 
                            x = x.tmp, 
                            prtx = dat.comb$prtx, 
                            trt = dat.comb$trt,
                            trtNew = c(dat.rest$trt.new, 
                                       rep(-1000, length(x))),
                            inNode = c(rep(0, length(dat.rest$y)), 
                                       rep(1, length(x))))
        }
        # browser(expr = (v.name == "x1"))
        
        # set up split parameters
        splitParams <- list(isCtg = is.ctg, 
                            nodeSize = min.ndsz, 
                            trtSize = n0, 
                            maxRisk = haoda.ae.level, 
                            useOtherNodes = use.other.nodes, 
                            maxScore = max.score)
        tmp <- splitConditional(zcut, zcutCat, datMatrix, splitParams)
      return(tmp)
    } else{
      return(list(output = NA, direction = NA))
    }
  })
    
    # Extract best split info
    out.idx <- which.max(lapply(outSplit, function(xx) max(xx$output)))
    x <- .subset2(dat, names(out.idx))
    temp <- sort(unique(x))
    tmp.idx <- which.max(outSplit[[out.idx]]$output)
    if(outSplit[[out.idx]]$direction[tmp.idx] %in% c("l", "r")){
      vname <- names(out.idx)
      var <- which(colnames(dat) == vname)
      max.score <- max(outSplit[[out.idx]]$output)
      if(!is.null(ctg) & names(out.idx) %in% colnames(dat)[ctg]){
        best.cut <- paste(unlist(power.set(temp)[tmp.idx]), collapse=",")
      } else{
        best.cut <- temp[-length(temp)][tmp.idx]
      }
      
      cut <- cbind(outSplit[[out.idx]]$direction[tmp.idx], 
                   as.character(best.cut))
    }
  }
  
  # when testing data is provided, assess new treatment assignment 
  # using testing sample and the rule calculated from training sample
  # var is the covariates calcualted before where spliting adopts. 
  # best.cut is the cutting point.
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
    if (!is.na(var)) {
      if (is.element(var,ctg)) {
        if(cut[1] == "l"){
          grp.test <- sign(is.element(test[,var], unlist(strsplit(best.cut, ","))))
        } else{
          grp.test <- sign(!is.element(test[,var], unlist(strsplit(best.cut, ","))))
        }
      } else  {
        if(cut[1]=="l"){
          grp.test <- sign(test[,var] <= best.cut)
        } else{
          grp.test <- sign(test[,var] > best.cut)
        }
      }
      
      if(is.null(test$status)) test$status <- 1
      if(is.null(test$KM.cens)) test$KM.cens <- 1
      score.test <- switch(outcome, 
                           time = s.itrtest(test, z=grp.test, n0=n0, aug = AIPWE), 
                           ae = itrtest(test, z=grp.test, n0=n0, aug = AIPWE))
      if (!is.na(score.test)){
        out$name.l <- name.l
        out$name.r <- name.r
        if(cut[1]=="l"){
          out$left.test <- test[grp.test==1,  ]
          out$right.test <- test[grp.test==0,  ]
        } else{
          out$left.test <- test[grp.test==0,  ]
          out$right.test <- test[grp.test==1,  ]
        }
        if (is.element(var,ctg)) {
          if(cut[1] == "l"){
            out$left  <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
            out$right <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))), ]
          } else{
            out$left  <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
            out$right <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))), ]
          }
        } else {
          if(cut[1]=='l'){
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(1,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(0,n=sum(dat[,var]> best.cut)))
          }else{
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(0,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(1,n=sum(dat[,var]> best.cut)))
          }  
        }
      } else {
        var <- NA
        vname <- NA
        cut <- NA
        max.score <- NA
      }
      # output results from both testing and training data.
      if(!is.na(var)){  
        out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, 
                               trt.effect=trt.effect,var = var, 
                               vname=vname, cut.1 = cut[1], cut.2 = cut[2], 
                               score=ifelse(max.score==-1e20, NA, max.score), 
                               score.test=score.test, size.test=n.test)
      } else{
        out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, 
                               trt.effect=trt.effect, var = NA, 
                               vname=NA, cut.1 = NA, cut.2 = NA, 
                               score=NA,score.test=NA, size.test=n.test)
      }
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = NA, 
                             vname=NA, cut.1=NA, cut.2=NA, score=NA,score.test=NA, size.test=n.test)
    }
  } else {
    # if no testing data output results from training data only.
    if (!is.na(var)) {
      out$name.l <- name.l
      out$name.r <- name.r
      if (is.element(var,ctg)){
        out$left  <- dat[is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
        out$right <- dat[!is.element(dat[,var], unlist(strsplit(best.cut, ","))),]
      } else {
        if(cut[1]=='l'){
          out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(1,n=sum(dat[,var]<= best.cut)))
          out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(0,n=sum(dat[,var]> best.cut)))
        }else{
          out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(0,n=sum(dat[,var]<= best.cut)))
          out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(1,n=sum(dat[,var]> best.cut)))
        }  
      }
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = var, 
                             vname=vname, cut.1 = unique(cut[,1]), cut.2 = paste(cut[,2], collapse = ','),
                             score=ifelse(max.score==-1e20, NA, max.score))
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,var=NA, 
                             vname=NA, cut.1= NA,cut.2=NA, score=NA)
    }
  }
  out 
}