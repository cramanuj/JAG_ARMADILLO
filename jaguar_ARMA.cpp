/*
 Author: Chaitanya Acharya
 Date: Aug 27, 2015
 ***RCPP ARMADILLO version of JAGUAR***

 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: Y -> Residuals from the model
			k -> vector indicating number of tissues per observation
			snp -> vector of mean centered genotypes 
			X -> Stacked matrix of genotypes of dimension k*nobs by kmax where kmax is max(k)
			Eps -> \hat{\Epsilon} from the model
			Tau -> \hat{\Tau} from the model
 
 Returns: P-value
 
*/

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdint.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

SEXP jagARMA(const arma::mat Yhat, arma::vec k, arma::vec snp, const arma::mat X, double Eps, double Tau){
	
	// Initialize the parameters
	int ksum = sum(k); int kmax = max(k); int nobs = snp.size();
	arma::mat iV(ksum,ksum);		// V^{-1} variance-cov matrix of Y
	arma::mat G(ksum,1);			// G matrix (function of mean-centered genotypes)
	iV.fill(0);  G.fill(0);

	int row_st=0; int col_st=0; int row_end=0; int col_end=0; int nt = 0;
	int row_st_x=0; int row_end_x=0;
	
	//Loop over every observation
	for(int i=0; i<nobs; i++){
		double V1 = (Eps+(k[i]-1)*Tau) / ( (Eps*Eps)+k[i]*Eps*Tau );
		double V2 = - (Tau)/(Eps*Eps+k[i]*Eps*Tau);
		arma::mat id_iV(k[i],k[i]); id_iV.fill(V2);
		arma::vec V1_vec = rep(V1,k[i]);
		arma::vec iG_rep = rep(snp[i],k[i]);
		id_iV.diag() = V1_vec;
		nt+=k[i];
		row_end = col_end = nt;
		row_end_x = nt-1;
		G.submat(row_st_x,0,row_end_x,0) = iG_rep;
		row_st_x = row_end_x+1;
		iV.submat(row_st,col_st,row_end-1,col_end-1)=id_iV;
		row_st = col_st = col_end;
	}
	
	const arma::mat GGt = G*trans(G);										// Matrix crossproduct
	const arma::mat XXt = X*trans(X);										// Matrix crossproduct
	
	// Compute variance-cov matrix and the optimal weights

	double meanB = trace(trans(G)*iV*G);
	double meanG = 0.5 * trace(trans(X)*iV*X);
	double varB = 2 * trace(trans(G)*iV*G*trans(G)*iV*G);
	double varG = 0.5 * trace(trans(X)*iV*X*trans(X)*iV*X);
	double cov = trace(iV*GGt*iV*XXt);
	double a_gam = (varB-cov)/(varB+varG-2*cov); 
	double a_beta = (varG-cov)/(varB+varG-2*cov);

	// Compute the score test statistic and the corresponding p-value using
	// Satterthwaite method
	
	arma::mat U_psi = trans(Yhat)*iV*( (a_beta*GGt) + (0.5*a_gam*XXt) )*iV*Yhat;
	double meanPSI = (a_beta*meanB)+ (a_gam*meanG);
	double varPSI = (a_beta*varB*a_beta) + (a_gam*a_gam*varG);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = U_psi(0,0) / a1;
	return wrap(1 - R::pchisq(scaled_upsi,a2,1,0));
//	return R_NilValue;
	
}
