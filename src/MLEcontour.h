#ifndef _MLEcontour_H
#define _MLEcontour_H

class MLEcontour {

MLEmodel model;
arma::colvec par_hat;
int dist_num;
double MLLx;
double RatioLL;
double RadLimit;




public:
MLEcontour(SEXP, arma::colvec, int, double, double, double);
double compareMLLx();
arma::rowvec getContourPt( double);
double get_par_hat();

};
// end of class declaration


#endif