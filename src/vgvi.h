#ifndef VGVI
#define VGVI

#include <Rcpp.h>

extern Rcpp::NumericVector VGVI_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
									Rcpp::S4 &greenspace, const Rcpp::NumericVector &greenspace_values,
									const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
									const int radius, const Rcpp::NumericVector & h0,
									const int fun, const double m, const double b);

#endif