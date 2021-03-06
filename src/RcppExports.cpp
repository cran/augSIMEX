// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// repeatgenerate
NumericVector repeatgenerate(NumericMatrix W, double lambda);
RcppExport SEXP _augSIMEX_repeatgenerate(SEXP WSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(repeatgenerate(W, lambda));
    return rcpp_result_gen;
END_RCPP
}
// scorecloglogcpp
NumericVector scorecloglogcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scorecloglogcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scorecloglogcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}
// scoregaussiancpp
NumericVector scoregaussiancpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scoregaussiancpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scoregaussiancpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}
// scorelogitcpp
NumericVector scorelogitcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scorelogitcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scorelogitcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}
// scoremofifiedglmcpp
NumericVector scoremofifiedglmcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, Function linkinv, Function var, Function mueta);
RcppExport SEXP _augSIMEX_scoremofifiedglmcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP linkinvSEXP, SEXP varSEXP, SEXP muetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< Function >::type linkinv(linkinvSEXP);
    Rcpp::traits::input_parameter< Function >::type var(varSEXP);
    Rcpp::traits::input_parameter< Function >::type mueta(muetaSEXP);
    rcpp_result_gen = Rcpp::wrap(scoremofifiedglmcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, linkinv, var, mueta));
    return rcpp_result_gen;
END_RCPP
}
// scorepoilogcpp
NumericVector scorepoilogcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scorepoilogcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scorepoilogcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}
// scorepoisqrtcpp
NumericVector scorepoisqrtcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scorepoisqrtcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scorepoisqrtcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}
// scoreprobitcpp
NumericVector scoreprobitcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, NumericVector offset);
RcppExport SEXP _augSIMEX_scoreprobitcpp(SEXP betaSEXP, SEXP YSEXP, SEXP DataMSEXP, SEXP DataM0SEXP, SEXP DataM1SEXP, SEXP phatSEXP, SEXP qhatSEXP, SEXP weightSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM(DataMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM0(DataM0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataM1(DataM1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phat(phatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qhat(qhatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreprobitcpp(beta, Y, DataM, DataM0, DataM1, phat, qhat, weight, offset));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_augSIMEX_repeatgenerate", (DL_FUNC) &_augSIMEX_repeatgenerate, 2},
    {"_augSIMEX_scorecloglogcpp", (DL_FUNC) &_augSIMEX_scorecloglogcpp, 9},
    {"_augSIMEX_scoregaussiancpp", (DL_FUNC) &_augSIMEX_scoregaussiancpp, 9},
    {"_augSIMEX_scorelogitcpp", (DL_FUNC) &_augSIMEX_scorelogitcpp, 9},
    {"_augSIMEX_scoremofifiedglmcpp", (DL_FUNC) &_augSIMEX_scoremofifiedglmcpp, 11},
    {"_augSIMEX_scorepoilogcpp", (DL_FUNC) &_augSIMEX_scorepoilogcpp, 9},
    {"_augSIMEX_scorepoisqrtcpp", (DL_FUNC) &_augSIMEX_scorepoisqrtcpp, 9},
    {"_augSIMEX_scoreprobitcpp", (DL_FUNC) &_augSIMEX_scoreprobitcpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_augSIMEX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
