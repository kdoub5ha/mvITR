\name{mvITR-package}
\alias{mvITR-package}
\alias{mvITR}
\docType{package}
\title{
  Discovery of risk controlled individualized treatment rules
}
\description{
  This package can be used for discovery of risk controlled individualized treatment rules (ITRs) via risk controlled decision trees (rcDT) or risk controlled random forests (rcRF).
}
\details{
  Use the grow.ITR() and Build.RF.ITR() functions for risk controlled ITR discovery. 
}
\author{
Kevin Doubleday, kdoub5ha@gmail.com.

Maintainer: Kevin Doubleday <kdoub5ha@gmail.com>
}
\references{
  [In Submission] Doubleday, K., Zhou, J., Fu, H. (2020), "Risk Controlled Decision Trees and Random Forests for Individualized Treatment Recommendations".
}
\keyword{ package }
#\seealso{
#  Optional links to other man pages
#}
\examples{
  \dontrun{
    set.seed(123)
    # Generate example data
    dat <- generateData(n = 1000)
    
    # Construct rcDT
    tre1 <- grow.ITR(data = dat, 
                 split.var = 1:10,
                 risk.threshold = 2.75,
                 lambda = 1,
                 efficacy = "y",
                 risk = "r",
                 col.trt = "trt",
                 col.prtx = "prtx")
    
    # Prune and select rcDT
    pruned1 <- prune(tre1, a = 0, risk.threshold = 2.75, lambda = 1)
    rcDT.fit <- treeCV(dat = dat, 
                       split.var = 1:10,
                       lambda = 1,
                       risk.threshold = 2.75, 
                       efficacy = "y",
                       risk = "r",
                       col.trt = "trt",
                       col.prtx = "prtx",
                       nfolds = 5)
    
    # Construct rcRF
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
    
    # Variable importances
    VI <- Variable.Importance.ITR(rcRF.fit, sort = FALSE)
    
    # Predictions
    preds.rcDT <- predict.ITR(rcDT.fit$best.tree.alpha, 
                              new.data = dat, 
                              split.var = 1:10)
    preds.rcRF <- predict.ITR(rcRF.fit, 
                              new.data = dat, 
                              split.var = 1:10)
  }
}
