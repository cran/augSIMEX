# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

repeatgenerate <- function(W, lambda) {
    .Call('_augSIMEX_repeatgenerate', PACKAGE = 'augSIMEX', W, lambda)
}

scorecloglogcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scorecloglogcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

scoregaussiancpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scoregaussiancpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

scorelogitcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scorelogitcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

scoremofifiedglmcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, linkinv, var, mueta) {
    .Call('_augSIMEX_scoremofifiedglmcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, linkinv, var, mueta)
}

scorepoilogcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scorepoilogcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

scorepoisqrtcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scorepoisqrtcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

scoreprobitcpp <- function(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset) {
    .Call('_augSIMEX_scoreprobitcpp', PACKAGE = 'augSIMEX', beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset)
}

