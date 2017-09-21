#include <Rcpp.h>
using namespace Rcpp;



NumericMatrix scoreglm(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericVector weight, Function linkinv, Function var, Function mueta) {
    int i,j;
    int nrow=DataM.nrow();
    int ncol=DataM.ncol();
    NumericMatrix out(nrow,ncol);
    double eta=0;
    double mu;

    for(i = 0; i < nrow; ++i)
     { eta=0;
     for (j = 0; j < ncol; ++j){
      eta += beta(j)*DataM(i,j);
    }
      mu = Rcpp::as<double>(linkinv(eta));
     for (j = 0; j < ncol; ++j){
      out(i,j) = weight(i) * Rcpp::as<double>(mueta(eta))/Rcpp::as<double>(var(mu))*(Y(i)-mu)* DataM(i,j);
     }
    }
    return out;
  }


  // [[Rcpp::export]]
NumericVector scoremofifiedglmcpp(NumericVector beta, NumericVector Y, NumericMatrix DataM, NumericMatrix DataM0, NumericMatrix DataM1, NumericVector phat, NumericVector qhat, NumericVector weight, Function linkinv, Function var, Function mueta) {
  int i,j;
  int nrow=DataM.nrow();
  int ncol=DataM.ncol();
  NumericMatrix scorez0(nrow,ncol),scorez1(nrow,ncol);
  NumericVector value(ncol);

  scorez0=scoreglm(beta,Y,DataM0,weight,linkinv,var,mueta);
  scorez1=scoreglm(beta,Y,DataM1,weight,linkinv,var,mueta);

  for (j = 0; j < ncol; ++j){
    value(j)=0;
    for(i = 0; i < nrow; ++i){
      value(j)+=((1-DataM(i,ncol-1))*(scorez0(i,j)*(1-phat(i))-scorez1(i,j)*qhat(i))-DataM(i,ncol-1)*(scorez0(i,j)*phat(i)-scorez1(i,j)*(1-qhat(i))))/(1-phat(i)-qhat(i));
    }
  }

  return value;
}
