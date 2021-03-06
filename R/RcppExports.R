# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

SendDown <- function(cutPoint, splitVar, Data, treNodes, direction) {
    .Call('_mvITR_SendDown', PACKAGE = 'mvITR', cutPoint, splitVar, Data, treNodes, direction)
}

estITR <- function(input) {
    .Call('_mvITR_estITR', PACKAGE = 'mvITR', input)
}

estOpt <- function(y, r, trt, prtx, rule, tau, lambda) {
    .Call('_mvITR_estOpt', PACKAGE = 'mvITR', y, r, trt, prtx, rule, tau, lambda)
}

splitConditional <- function(zcut, zcutCat, datMatrix, parameters) {
    .Call('_mvITR_splitConditional', PACKAGE = 'mvITR', zcut, zcutCat, datMatrix, parameters)
}

splitContinuous <- function(zcut, zcutCat, datMatrix, parameters) {
    .Call('_mvITR_splitContinuous', PACKAGE = 'mvITR', zcut, zcutCat, datMatrix, parameters)
}

splitSurvival <- function(zcut, zcutCat, datMatrix, parameters) {
    .Call('_mvITR_splitSurvival', PACKAGE = 'mvITR', zcut, zcutCat, datMatrix, parameters)
}

