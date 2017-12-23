#ifndef _contour_H
#define _contour_H


#ifdef __cplusplus

#include <RcppArmadillo.h>

// abremDebias code
RcppExport SEXP MLEtryLL(SEXP arg1, SEXP arg3, SEXP arg4, SEXP arg5, SEXP arg6);




RcppExport  void R_init_contour(DllInfo* info);

#endif
#endif
