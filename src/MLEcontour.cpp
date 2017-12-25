

#include "contour.h"
#include "MLEmodel.h"
#include "MLEcontour.h"
#include <math.h>

    using namespace Rcpp ;
	
// class implementation

MLEcontour::MLEcontour( SEXP arg1, arma::colvec arg2, int arg3, double arg4, double arg5) {
	model = MLEmodel(arg1);
	par_hat = arg2;
	dist_num = arg3;
	MLLx =  arg4;
	RatioLL = arg5;
}

int MLEcontour::compareMLLx()  {
	int value = 0;
	double calcMLLx = model.tryLL(par_hat, dist_num);
	if(abs(MLLx-calcMLLx) < 1e-10) {
		value = 1;
	}
	return value;
}



	// Exported Functions
SEXP testMLLx(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5) {
	arma::colvec par=Rcpp::as<arma::colvec>(arg2);
	int dist_num=Rcpp::as<int>(arg3);
	double MLLx = Rcpp::as<double>(arg4);
	double RatioLL = Rcpp::as<double>(arg5);
	MLEcontour mycontour(arg1, par, dist_num, MLLx, RatioLL);

	return wrap(mycontour.compareMLLx());
}