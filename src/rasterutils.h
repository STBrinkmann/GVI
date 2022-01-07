#ifndef RASTERUTILS
#define RASTERUTILS

#include <Rcpp.h>

extern void recycle(std::vector<int> &x, std::vector<int> &y);
extern Rcpp::IntegerVector cellFromColRowSensitive(Rcpp::S4 &raster, Rcpp::IntegerVector &rcpp_col, Rcpp::IntegerVector &rcpp_row);
extern Rcpp::IntegerVector cellFromColRow(const Rcpp::IntegerVector &x, const Rcpp::IntegerVector &y, const int ncol);
extern Rcpp::IntegerMatrix colRowFromCell(const Rcpp::IntegerVector &cell, const int ncol);
extern Rcpp::NumericMatrix xyFromCell(Rcpp::S4 &raster, const Rcpp::IntegerVector &cell);
extern Rcpp::NumericMatrix xyFromCell(Rcpp::S4 &raster, int cell);
extern Rcpp::IntegerVector cellFromXY(Rcpp::S4 &raster, Rcpp::NumericMatrix xy);
extern int cellFromXY2(Rcpp::S4 &raster, double x, double y);

#endif