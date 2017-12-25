#ifndef _contour_H
#define _contour_H


#ifdef __cplusplus

#include <RcppArmadillo.h>

// MLEloglike included sign and tz arguments
//RcppExport SEXP MLEloglike(SEXP arg1, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6);
RcppExport SEXP MLEtryLL(SEXP arg1, SEXP arg3, SEXP arg4);
RcppExport SEXP testMLLx(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4, SEXP arg5);


RcppExport  void R_init_contour(DllInfo* info);

#endif
#endif
