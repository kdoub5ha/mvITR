## ----setup, include=FALSE------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE)
library(mvITR)


## ----Generate data set, results='markup', echo=TRUE----------------------
set.seed(123)
dat <- generateData(n = 1000)
str(dat)


## ----Summary plots, results='markup', echo=FALSE, results='show', fig.cap="Figure 1. Risk score distribution for simulated data", fig.width=7----

boxplot(dat$r ~ dat$trt, boxwex = 0.25,
        xlab = "Original Treatment Group",
        main = "Risk Distribution by Treatment Group",
        ylab = "Risk Score", axes = FALSE)
axis(1, at = 1:2, labels = c("Control", "Treated"),
     col = "white"); axis(2, las = 2)


## ----Grow a tree tau 2.75 - lambda 1, results='markup', echo=TRUE--------
tre1 <- grow.ITR(data = dat, 
                 split.var = 1:10,
                 risk.threshold = 2.75,
                 lambda = 1,
                 efficacy = "y",
                 risk = "r",
                 col.trt = "trt",
                 col.prtx = "prtx")
tre1$tree

## ----Grow a tree tau 2.75 - lambda 2, results='markup', echo=TRUE--------
tre2 <- grow.ITR(data = dat, 
                 split.var = 1:10,
                 risk.threshold = 2.75,
                 lambda = 2,
                 efficacy = "y",
                 risk = "r",
                 col.trt = "trt",
                 col.prtx = "prtx")
tre2$tree

## ----Code Tree Pruning, results='markup'---------------------------------
pruned1 <- prune(tre1, a = 0, risk.threshold = 2.75, lambda = 1)
pruned.display <- pruned1$result[,c(1:6,10,11)]
pruned.display$alpha <- sprintf("%.3f", as.numeric(pruned.display$alpha))
pruned.display$V <- sprintf("%.3f", as.numeric(pruned.display$V))
pruned.display$Benefit <- sprintf("%.3f", as.numeric(pruned.display$Benefit))
pruned.display$Risk <- sprintf("%.3f", as.numeric(pruned.display$Risk))
pruned.display


## ---- Cross Validated Pruning Model, results='hide', echo=TRUE-----------
rcDT.fit <- treeCV(dat = dat, 
                   split.var = 1:10,
                   lambda = 1,
                   risk.threshold = 2.75, 
                   efficacy = "y",
                   risk = "r",
                   col.trt = "trt",
                   col.prtx = "prtx",
                   nfolds = 5)


## ---- Code Forest Growth, results='markup'-------------------------------
set.seed(2)
rcRF.fit <- Build.RF.ITR(dat, 
                         split.var = 1:10, 
                         efficacy = "y", 
                         risk = "r",
                         col.trt = "trt",
                         col.prtx = "prtx",
                         risk.threshold = 2.75,
                         ntree = 100, 
                         lambda = 0.5)

## ---- Treatment Prediction, results='markup', echo=TRUE------------------
preds.rcDT <- predict.ITR(rcDT.fit$best.tree.alpha, 
                          new.data = dat, 
                          split.var = 1:10)$trt.pred
preds.rcRF <- predict.ITR(rcRF.fit, 
                          new.data = dat, 
                          split.var = 1:10)$trt.pred


## ---- Treatment Prediction Plot, results='markup', echo=FALSE, fig.cap="Figure 2. Prediction Comparisons for rcDT and rcRF Models.", fig.width=10, fig.height=6----

par(mfrow = c(1,2))
par(mar = c(5,5,5,1))

plot(dat$x1, dat$x2, pch = 16, 
     cex = 100, col = "lightgray", 
     xlab = expression(X[1]),
     ylab = expression(X[2]),
     main = paste0("Predictions from rcDT (Tree) Model\n", 
                   "Efficacy = ", 
                   sprintf("%.2f", mean(dat$y * (dat$trt == preds.rcDT) / 0.5)),
                   "\nRisk = ", 
                   sprintf("%.2f", mean(dat$r * (dat$trt == preds.rcDT) / 0.5))),
     axes = FALSE)
points(dat$x1, dat$x2, pch = 16, 
       cex = 0.8, 
       col = ifelse(preds.rcDT == 1, "forestgreen", "hotpink"))
axis(1); axis(2, las = 2)

plot(dat$x1, dat$x2, pch = 16, 
     cex = 100, col = "lightgray", 
     xlab = expression(X[1]),
     ylab = expression(X[2]),
     main = paste0("Predictions from rcRF (Forest) Model\n", 
                   "Efficacy = ", 
                   sprintf("%.2f", mean(dat$y * (dat$trt == preds.rcRF) / 0.5)),
                   "\nRisk = ", 
                   sprintf("%.2f", mean(dat$r * (dat$trt == preds.rcRF) / 0.5))),
     axes = FALSE)
points(dat$x1, dat$x2, pch = 16, 
       cex = 0.8, 
       col = ifelse(preds.rcRF == 1, "forestgreen", "hotpink"))
axis(1); axis(2, las = 2)


## ---- Code Variable Importance, results='markup'-------------------------
VI <- Variable.Importance.ITR(rcRF.fit, sort = FALSE)
do.call(cbind, VI)

