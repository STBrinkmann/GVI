#include <Rcpp.h>
#include "rsinfo.h"
#include "rasterutils.h"
#include "bresenham.h"

using namespace Rcpp;

int Sign2(const int dxy)
{
  if(dxy<0) return -1; 
  else if(dxy>0) return 1; 
  else return 0;
}


Rcpp::IntegerMatrix bresenham_map(const int x0, const int y0, const int radius, const int nc)
{
  int x1;
  int y1 = y0+radius;
  
  int Dx, Dy, Sx, Sy, X, Y, c;
  double R;
  
  Rcpp::IntegerMatrix out(radius+1,radius);
  std::fill( out.begin(),out.end(),NA_INTEGER );
  
  Dy = y1 - y0;
  Sy = Sign2(Dy); 
  
  for(int r=0; r<=radius; r++) {
    x1 = x0+r;
    
    // Segment length
    Dx = x1 - x0;
    
    // Increments
    Sx = Sign2(Dx);
    
    // Initial remainder
    R = radius / 2;
    
    X = x0;
    Y = y0;
    
    // Initial update (X, Y) and R
    Y+= Sy; R+= Dx; // Lateral move
    if(R >= Dy) {    
      X+= Sx; 
      R-= Dy; // Diagonal move
    }
    
    c = 0;
    while( ((x0-X)*(x0-X)+(y0-Y)*(y0-Y)) <= (radius*radius) ) {
      out(r,c) = Y*nc + X;
      
      // Update (X, Y) and R
      Y+= Sy; R+= Dx; // Lateral move
      if(R >= Dy) {    
        X+= Sx; 
        R-= Dy; // Diagonal move
      }
      c++;  
    }
  }
  
  return out;
}


Rcpp::NumericVector tangentsMap(const int x0, const int y0, const double h0, 
                                const int nr, const int nc,
                                const int r, const Rcpp::NumericVector &dsm_vec){
  const int range = (2*r)+1;
  Rcpp::NumericVector out(range*range, NA_REAL);
  int i = 0;
  
  for (int y = y0-r; y <= y0+r; ++y) {
    for (int x = x0-r; x<=x0+r; ++x) {
      const int cell = y*nc + x;
      
      if( !(y<0 || y >= nr || x<0 || x >= nc || Rcpp::NumericVector::is_na(dsm_vec[cell])) ){
        const double distance_traveled = sqrt((x0 - x)*(x0 - x) + (y0 - y)*(y0 - y));
        if(distance_traveled <= r) {
          out[i] = (dsm_vec[cell] - h0) / (distance_traveled);
        }
      }
      i++;
    }
  }
  return out;
}


void LoS(const Rcpp::IntegerMatrix &los_ref_mat,
         const Rcpp::NumericVector &tan_vec,
         Rcpp::LogicalVector &visibility_vec,
         const int &j, const int &r){
  
  double max_tangent = -9999;
  
  for(int i = 0; i < r; ++i) {
    int cell = los_ref_mat( j,i );
    
    if(!Rcpp::IntegerVector::is_na(cell)) {
      double this_tangent = tan_vec[cell];
      
      if(!Rcpp::NumericVector::is_na(this_tangent)) {
        if (this_tangent > max_tangent) {
          max_tangent = this_tangent;
          visibility_vec[cell] = TRUE;
        }
      }
    } else {
      break;
    }
  }
}