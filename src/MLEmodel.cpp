// traditional header file content - class declaration
#include "contour.h"
#include <math.h>

    using namespace Rcpp ;

class MLEmodel {

arma::colvec time;
arma::colvec qty;
Rcpp::NumericVector N;

double failcomp;
double suscomp;
double discomp;
double intcomp;
int endf;
int ends;
int endd;
int endil;
int endir;
arma::colvec fail;
arma::colvec nf;
arma::colvec susp;
arma::colvec ns;
arma::colvec disc;
arma::colvec nd;
arma::colvec left;
arma::colvec right;
arma::colvec ni;

public:
MLEmodel(SEXP);
double LogLike(arma::colvec, int, int, double);
SEXP MLE_Simplex(SEXP, arma::colvec, double, int);
SEXP dMaxLLdx( SEXP, arma::colvec, double);
};
// end of class declaration


// class implementation

MLEmodel::MLEmodel( SEXP arg1) {
	Rcpp::List L(arg1);
	time=Rcpp::as<arma::colvec>(L["fsdi"]);
	qty=Rcpp::as<arma::colvec>(L["q"]);
	N=L["N"];
// Provide a first non-sense element to front of time vector
// so that position math works when Nf is zero.
	time.insert_rows(0,1);
	qty.insert_rows(0,1);

	endf=N[0];
	ends=endf+N[1];
	endd=ends+N[2];
	endil=endd+N[3];
	endir=endil+N[3];

	if(N[0]>0)  {
		fail=time.rows(1,endf);
		nf=qty.rows(1,endf);
	}

	if(N[1]>0)  {
		susp=time.rows(endf+1,ends);
		ns=qty.rows(endf+1,ends);
	}

	if(N[2]>0)  {
		disc=time.rows(ends+1,endd);
		nd=qty.rows(ends+1,endd);
	}
	if(N[3]>0)  {
		left=time.rows(endd+1,endil);
		right=time.rows(endil+1,endir);
		ni=qty.rows(endd+1,endil);
	}

}


////************* Method tryLL ***************/
double MLEmodel::tryLL(arma::colvec par, int sign, int dist_num, double tz)  {

		double failcomp =0.0;
		double suscomp =0.0;
		double discomp =0.0;
		double intcomp =0.0;

		if(dist_num==1) {
			if(N[0]>0)  {
				for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dweibull(fail(i)-tz,par(0),par(1),1);
				}
			}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
//  Need to exempt any zero or negative values
					if(susp(i)-tz>0)  {
						suscomp=suscomp+ns(i)*R::pweibull(susp(i)-tz,par(0),par(1),0,1);
					}
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::pweibull(disc(i)-tz,par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
// if left(i) should become zero or less, then this data needs to be handled as if left censored on right(i)
					if(left(i)-tz>0)  {
						intcomp=intcomp+ni(i)*log(
						R::pweibull(left(i)-tz,par(0),par(1),0,0) -
						R::pweibull(right(i)-tz,par(0),par(1),0,0)
						);
					}else{
						intcomp=intcomp+ni(i)*log(1-R::pweibull(right(i)-tz,par(0),par(1),0,0));
					}
				}
			}

		}
		else if(dist_num==2)  {
				if(N[0]>0)  {
					for(int i=0; i<N[0]; i++)  {
					failcomp=failcomp+nf(i)*R::dlnorm(fail(i)-tz,par(0),par(1),1);
					}
				}

			if(N[1]>0)  {
				for(int i=0; i<N[1]; i++)  {
//  Need to exempt any zero or negative values
					if(susp(i)-tz>0)  {
						suscomp=suscomp+ns(i)*R::plnorm(susp(i)-tz,par(0),par(1),0,1);
					}
				}
			}
			if(N[2]>0)  {
				for(int i=0; i<N[2]; i++)  {
					discomp=discomp+nd(i)*log(1-R::plnorm(disc(i)-tz,par(0),par(1),0,0));
				}
			}
			if(N[3]>0)  {
				for(int i=0; i<N[3]; i++)  {
// if left(i) should become zero or less, then this data needs to be handled as if left censored on right(i)
					if(left(i)-tz>0)  {
						intcomp=intcomp+ni(i)*log(
						R::plnorm(left(i)-tz,par(0),par(1),0,0) -
						R::plnorm(right(i)-tz,par(0),par(1),0,0)
						);
					}else{
						intcomp=intcomp+ni(i)*log(1-R::plnorm(right(i)-tz,par(0),par(1),0,0));
					}
				}
			}
		}
	return sign*(failcomp+suscomp+discomp+intcomp);
}


	// Exported Functions

	SEXP MLEtryLL(SEXP arg1, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6)  {
		MLEmodel mymodel(arg1);
		arma::colvec par=Rcpp::as<arma::colvec>(arg3);
		int dist_num=Rcpp::as<int>(arg4);
		int sign=Rcpp::as<int>(arg5);
		double tz=Rcpp::as<double>(arg6);
		return wrap(mymodel.tryLL(par, sign, dist_num, tz));
	}




