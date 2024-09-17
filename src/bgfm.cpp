// This script implement generalized factor model based on variational inference.
// Date: 2023-02-11

// Revised log:
// 2023-01-16: replace the for loop with matrix operation in updating variational parameters of Z of E-step!
// 2023-07-31 Regard the factor matrix as an infinite parameter and estimate it.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//'
// diag(W0^t* Cki * W0)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}




double calELBO(const field<mat>& Xf, const vec& typeID, const sp_mat& A,  const mat& Mu_y,
               const mat& S_y,
               const vec& invLambda, const mat& B, const rowvec& mu, const mat& H){

  int i, p, d = Xf.n_elem, n = Xf(0).n_rows, q= B.n_cols;
  double  logS, entropy=0.0, ELBO=0.0;
  double term1=0.0, term2=0.0;
  vec p_vec(d+1, fill::zeros);
  for(i=0; i<d; ++i){
    p = Xf(i).n_cols;
    p_vec(i+1) = p_vec(i) + p;

    if(typeID(i)==2){
      term1 += accu(Xf(i) % Mu_y.cols(p_vec(i),p_vec(i+1)-1));
      term2 += -accu(exp(Mu_y.cols(p_vec(i),p_vec(i+1)-1) + S_y.cols(p_vec(i),p_vec(i+1)-1)/2));
    }
    if(typeID(i)==3){
      rowvec nvec = max(Xf(i));
      term1 += accu((Xf(i)-repmat(nvec, n, 1)) % Mu_y.cols(p_vec(i),p_vec(i+1)-1));  //
      term2 += -accu(repmat(nvec, n, 1) %  exp(-Mu_y.cols(p_vec(i),p_vec(i+1)-1) + S_y.cols(p_vec(i),p_vec(i+1)-1)/2)); // term2 is not written now!
    }
  }



  // if typeID = 1: Mu_y submatrix == x_1
  mat dX = (Mu_y - H * B.t() - repmat(mu, n, 1) - A) % repmat(sqrt(invLambda.t()), n, 1);
  mat LamB = B % repmat(sqrt(invLambda), 1, q);
  double dimreduction_term1 = -0.5*(accu(dX % dX)+
                                    accu(S_y % repmat(invLambda.t(), n, 1)) - n* accu(log(invLambda)) );


  entropy = 0.5*accu(log(S_y+ (S_y==0)));

  ELBO = term1 + term2 + dimreduction_term1 + entropy;
  return ELBO;
}


void VB_Estep(const field<mat>& Xf, const sp_mat& A, const vec& typeID, mat& Mu_y,  mat& S_y,
              const vec& invLambda, const mat& B, const rowvec& mu, mat& H){

  // Basic information:
  int i,  d = Xf.n_elem, n=Xf(0).n_rows, q= B.n_cols, p;
  vec p_vec(d+1, fill::zeros);
  // Mu_y is a n*p matrix
  mat tZ = H * B.t() + repmat(mu, n, 1) + A;
  mat invLamMat = repmat(invLambda.t(), n, 1);

  //  ## VB E-step
  // update posterior variance of y: S_y
  // update posterior mean of y: Mu_y
  // double elbo0 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H, S_h, Sigma_h);

  for(i=0; i<d; ++i){ // update Mu_y and S_y by the variable types

    p = Xf(i).n_cols;
    p_vec(i+1) = p_vec(i) + p;
    //Rprintf("good 2, i= %d!\n", i);
    // if(typeID(i) == 1){ // Gaussian:
    //   tX.cols(p_vec(i), p_vec(i+1)-1) = Xf(i) - repmat(mu.subvec(p_vec(i), p_vec(i+1)-1), n, 1);
    // }
    //Rprintf("good 3, i= %d!\n", i);
    if(typeID(i)==2){ // Poisson
      mat Mu_y1 = Mu_y.cols(p_vec(i), p_vec(i+1)-1);
      Mu_y1 = (Xf(i) -  exp(Mu_y1) % (1-Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1)% tZ.cols(p_vec(i), p_vec(i+1)-1)) /
        (exp(Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
      Mu_y.cols(p_vec(i), p_vec(i+1)-1) = Mu_y1;
      S_y.cols(p_vec(i), p_vec(i+1)-1) = 1.0 / (exp(Mu_y1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );

      //tX.cols(p_vec(i), p_vec(i+1)-1) = Mu_y1 - repmat(mu.subvec(p_vec(i), p_vec(i+1)-1), n, 1);
    }

    //Rprintf("good 4, i= %d!\n", i);

    if(typeID(i)==3){ // Binomial
      rowvec nvec = max(Xf(i));
      mat Mu_y2 = Mu_y.cols(p_vec(i), p_vec(i+1)-1); // take out this submatrix
      Mu_y2 = (Xf(i)- (1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1)% tZ.cols(p_vec(i), p_vec(i+1)-1)) /
        ((1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );
      Mu_y.cols(p_vec(i), p_vec(i+1)-1) = Mu_y2;
      S_y.cols(p_vec(i), p_vec(i+1)-1) = 1.0 / ( (1/(1+exp(-Mu_y2))) %(1- 1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + invLamMat.cols(p_vec(i), p_vec(i+1)-1) );

      //tX.cols(p_vec(i), p_vec(i+1)-1) = Mu_y2 - repmat(mu.subvec(p_vec(i), p_vec(i+1)-1), n, 1);
    }
    // Rprintf("good 5, i= %d!\n", i);

  }

  // double elbo1 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H, S_h, Sigma_h);
  // Rprintf("VB dY= %4f\n", elbo1 - elbo0);


}




// [[Rcpp::export]]
Rcpp::List VB_GFMcpp(const Rcpp::List& XList, const arma::vec& typeID, const arma::sp_mat& A,
                      const arma::mat& Mu_y_int,
                      const arma::mat& S_y_int,
                      const arma::vec& invLambda_int, const arma::mat& B_int,
                      const arma::rowvec& mu_int, const arma::mat& H_int,
                      const double& epsELBO, const int& maxIter, const bool& verbose){
  // typeID: 1 means Gaussian; 2 means Poisson; 3 means Binomial;
  // whether add nvec for Binomial variables???? Now assume n_j = 1 for all.
  int i, d = XList.length();
  int d2 = typeID.n_elem;
  if(d != d2){
    stop("The length of XList must be equal to the length of typeID!");
  }

  field<mat> Xf(d);
  for(i=0; i<d; ++i){ // put the component of list into a field.
    mat Xtmp = XList[i];
    Xf(i) = Xtmp;
  }
  int n = Xf(0).n_rows, q = B_int.n_cols;


  // Initialize
  mat Mu_y(Mu_y_int), S_y(S_y_int),  H(H_int), B(B_int);
  vec invLambda(invLambda_int);
  rowvec mu(mu_int);

  vec ELBO_vec(maxIter);
  ELBO_vec(0) = -1e20;
  mat dX, BL;
  int iter;


  for(iter = 1; iter < maxIter; ++iter){


    // Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(Xf, A, typeID, Mu_y, S_y, invLambda, B, mu, H);

    // Rprintf("Finish E step!\n");
    //VB M-step

    // double elbo1 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);
    //update B
    // Rprintf("update B\n");
    B = trans(Mu_y-repmat(mu, n, 1)-A)  * H * inv_sympd(H.t() * H);
    // double elbo2 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);
    // Rprintf("dB= %4f \n", elbo2 - elbo1);

    // update H
    //BL = B % repmat(sqrt(invLambda), 1, q);
    //H = (Mu_y-repmat(mu, n, 1)-A) * diagmat(invLambda) * B * inv_sympd(BL.t() * BL);
    H = (Mu_y-repmat(mu, n, 1)-A) * (B % repmat(invLambda, 1, q)) * inv(B.t()*  (B % repmat(invLambda, 1, q)));
    // double elbo3 =  calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);
    // Rprintf("dH= %4f\n", elbo3 - elbo2);

    // update mu
    // Rprintf("update mu\n");
    mu = mean(Mu_y - H * B.t()-A); // Mu_y for gaussian: x_1
    // double elbo5 = calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);
    // Rprintf("dmu= %4f \n", elbo5 - elbo3);


    // update Lambda
    // Rprintf("update Lambda\n");
    dX = Mu_y -  H * B.t() - repmat(mu, n, 1)-A;
    vec Lambda =  trans(mean(dX % dX + S_y)); // typeID=1: S_y submatix must be zeros.
    Lambda.elem( find(Lambda < 1e-4) ).fill(1e-4); // increase the numerical stability!
    invLambda = 1.0 / Lambda;
    // double elbo6 =  calELBO( Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);
    // Rprintf("dLambda= %4f\n", elbo6 - elbo5);


    ELBO_vec(iter) =  calELBO(Xf, typeID, A, Mu_y, S_y, invLambda, B, mu, H);

    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n",
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }

  // output return value
  List resList = List::create(
    Rcpp::Named("H") = H,
    Rcpp::Named("B") = B,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("invLambda") = invLambda,
    Rcpp::Named("Mu_y") = Mu_y,
    Rcpp::Named("S_y") = S_y,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);

}
