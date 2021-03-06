#include <Rcpp.h>
#include "rsinfo.h"
#include "rasterutils.h"
#include "bresenham.h"
#include "integrate.h"
#include "vgvi.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector VGVI_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                             Rcpp::S4 &greenspace, const Rcpp::NumericVector &greenspace_values,
                             const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                             const int radius, const Rcpp::NumericVector & h0,
                             const int fun, const double m, const double b)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  Rcpp::NumericVector output(input_cells.size(), NA_REAL);
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);

  // Main loop
  for(int k=0; k<input_cells.size(); k++)
  {
    // Outpur: x0/y0 is always visible
    Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
    viewshed[c0_ref] = input_cells[k]+1;
    
    // A. Viewshed Analysis
    if(h0[k] > dsm_values[input_cells[k]]) {
      Rcpp::NumericVector max_tan_vec(r);
      int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
      for(int i = 0; i < (r*8); i++){
        // Re-use tangents calculated from prior LoS
        int k_i = 0;
        if(i > 0){
          for(int j = 0; j < r; j++){
            k_i = j;
            if(los_ref_vec[i*r + j] != los_ref_vec[(i-1)*r + j]){
              break;
            }
          }
        }
        double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
        
        
        for(int j = k_i; j < r; j++){
          const int los_ref_cell = los_ref_vec[i*r + j];
          if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
            // Project reference cells to actual cell value
            const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
            const int row = trunc(cell/ras.ncol);
            const int col = cell - (row * ras.ncol);
            const int dcol = abs(col-x0_o[k]);
            
            if(!(cell<0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
              // Compute tangents of LoS_vec and update output
              const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
              const double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
              
              if(this_tan > max_tan){
                max_tan = this_tan;
                viewshed[los_ref_cell] = cell+1;
              }
            }
            max_tan_vec[j] = max_tan;
          } else {
            break;
          }
        }
      }  
    }
    
    Rcpp::IntegerVector viewshed_na = Rcpp::na_omit(viewshed);
    
    
    // B: Greenspace Visibility Index
    const int cs = viewshed_na.size();
    
    // Get XY coordinates of visible cells
    const Rcpp::NumericMatrix dsm_xy = xyFromCell(dsm, viewshed_na);
    const Rcpp::NumericMatrix xy0 = xyFromCell(dsm, input_cells[k]);
    
    // Calculate distance
    Rcpp::IntegerVector dxy(cs);
    for(int i = 0; i < cs; i++){
      double d = sqrt( ((xy0(0,0) - dsm_xy(i,0))*(xy0(0,0) - dsm_xy(i,0))) + ((xy0(0,1) - dsm_xy(i,1))*(xy0(0,1) - dsm_xy(i,1))) );
      int di = round(d);
      if(di < 1) {
        dxy[i] = 1;
      } else {
        dxy[i] = di;
      }
    }
    
    
    // Intersect XY with greenspace mask
    const Rcpp::IntegerVector greenspace_cells = cellFromXY(greenspace, dsm_xy);
    Rcpp::NumericVector greenspace_cell_values(cs, 0.0);
    for(int i = 0; i < cs; i++){
      int gs_cell = greenspace_cells[i]; 
      if(!Rcpp::IntegerVector::is_na(gs_cell)){
        double gs = greenspace_values[gs_cell];
        greenspace_cell_values[i] =  Rcpp::NumericVector::is_na(gs) ? 0.0 : gs;
      }
    }
    
    // Get number of green visible cells and total visible cells per distance
    const double max_d = max(dxy);
    const Rcpp::IntegerVector dxy_seq = seq(1, max_d);
    const int n = dxy_seq.size();
    Rcpp::IntegerVector visibleTotal(n);
    Rcpp::IntegerVector visibleGreen(n);
    
    for(int i = 0; i < cs; i++){
      int this_dxy = dxy[i];
      visibleTotal[this_dxy-1] += 1;
      visibleGreen[this_dxy-1] += greenspace_cell_values[i];
    }
    for(int i = 0; i < n; i++){
      if(visibleTotal[i] == 0){
        visibleTotal[i] = 1;
      }
    }
    
    
    // Proportion of visible green cells
    if(max_d == 1){
      output[k] = visibleGreen[0]/visibleTotal[0];
    }
    Rcpp::NumericVector raw_GVI(n);
    for(int i = 0; i < n; i++){
      raw_GVI[i] = (double)visibleGreen[i] / visibleTotal[i];
    }
    
    // Normalize distance
    Rcpp::NumericVector nDxy(n);
    for(int i = 0; i < n; i++){
      nDxy[i] = dxy_seq[i] / (double)radius; //max_d; // TAKE radius instead and save max_d / r for upper limit
    }
    
    // Calculate weights by taking the proportion of the integral of each step from the integral of the whole area
    //double big_integral = integrate(0, 1, 200, fun, m, b);
    const double min_dist = min(nDxy);
    Rcpp::NumericVector decayWeights(n);
    for(int i = 0; i < n; i++){
      double d = nDxy[i];
      decayWeights[i] = integrate(d-min_dist, d, 200, fun, m, b);
    }
    const double big_integral = sum(decayWeights);
    
    // Proportion of visible green
    double vgvi_sum = 0.0; 
    for(int i = 0; i < n; i++){
      vgvi_sum += raw_GVI[i] * (decayWeights[i]/big_integral);
    }
    
    output[k] = vgvi_sum;
  }
  
  return output;
}


// [[Rcpp::export]]
Rcpp::IntegerVector viewshed_cpp_sum1(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                      const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                      const int radius, const Rcpp::NumericVector & h0)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  Rcpp::IntegerVector output(input_cells.size(), NA_INTEGER);
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);

  // Main loop
  for(int k=0; k<input_cells.size(); k++){
    // Outpur: x0/y0 is always visible
    Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
    viewshed[c0_ref] = input_cells[k]+1;
    
    if(h0[k] > dsm_values[input_cells[k]]){
      // 3. Viewshed Analysis
      Rcpp::NumericVector max_tan_vec(r);
      int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
      for(int i = 0; i < (r*8); i++){
        double max_tan = -9999.0;
        
        for(int j = 0; j < r; j++){
          const int los_ref_cell = los_ref_vec[i*r + j];
          if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
            // Project reference cells to actual cell value
            const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
            const int row = trunc(cell/ras.ncol);
            const int col = cell - (row * ras.ncol);
            const int dcol = abs(col-x0_o[k]);
            
            if(!(cell<0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
              // Compute tangents of LoS_vec and update output
              const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
              const double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
              
              if(this_tan > max_tan){
                max_tan = this_tan;
                viewshed[los_ref_cell] = cell+1;
              }
            }
            max_tan_vec[j] = max_tan;
          } else {
            break;
          }
        }
      }
    }
    Rcpp::IntegerVector viewshed_na = Rcpp::na_omit(viewshed);
    output[k] = viewshed_na.size();
  }
  
  
  
  return output;
}


// [[Rcpp::export]]
Rcpp::IntegerVector viewshed_cpp_sum2(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                      const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                      const int radius, const Rcpp::NumericVector & h0)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  Rcpp::IntegerVector output(input_cells.size(), NA_INTEGER);
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  const Rcpp::IntegerVector los_start = shared_LoS(r, los_ref_vec);
  
  // Main loop
  for(int k=0; k<input_cells.size(); k++){
    // Outpur: x0/y0 is always visible
    Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
    viewshed[c0_ref] = input_cells[k]+1;
    
    if(h0[k] > dsm_values[input_cells[k]]){
      // 3. Viewshed Analysis
      Rcpp::NumericVector max_tan_vec(r);
      int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
      for(int i = 0; i < (r*8); i++){
        // Re-use tangents calculated for prior LoS
        int k_i = los_start[i];
        double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
        
        
        for(int j = k_i; j < r; j++){
          const int los_ref_cell = los_ref_vec[i*r + j];
          if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
            // Project reference cells to actual cell value
            const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
            const int row = trunc(cell/ras.ncol);
            const int col = cell - (row * ras.ncol);
            const int dcol = abs(col-x0_o[k]);
            
            if(!(cell<0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
              // Compute tangents of LoS_vec and update output
              const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
              const double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
              
              if(this_tan > max_tan){
                max_tan = this_tan;
                viewshed[los_ref_cell] = cell+1;
              }
            }
            max_tan_vec[j] = max_tan;
          } else {
            break;
          }
        }
      }
    }
    Rcpp::IntegerVector viewshed_na = Rcpp::na_omit(viewshed);
    output[k] = viewshed_na.size();
  }
  
  
  
  return output;
}

// [[Rcpp::export]]
Rcpp::IntegerVector viewshed_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                 const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                 const int radius, const Rcpp::NumericVector & h0)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)round(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  
  int k=0;
  
  // Outpur: x0/y0 is always visible
  Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
  viewshed[c0_ref] = input_cells[k]+1;
  
  if(h0[k] > dsm_values[input_cells[k]]){
    // 3. Viewshed Analysis
    Rcpp::NumericVector max_tan_vec(r);
    int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
    
    for(int i = 0; i < (r*8); i++){
      // Re-use tangents calculated for prior LoS
      int k_i = 0;
      if(i > 0){
        for(int j = 0; j < r; j++){
          k_i = j;
          if(los_ref_vec[i*r + j] != los_ref_vec[(i-1)*r + j]){
            break;
          }
        }
      }
      double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
      
      
      for(int j = k_i; j < r; j++){
        const int los_ref_cell = los_ref_vec[i*r + j];
        if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
          // Project reference cells to actual cell value
          const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
          const int row = trunc(cell/ras.ncol);
          const int col = cell - (row * ras.ncol);
          const int dcol = abs(col-x0_o[k]);
          
          if(!(cell<0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
            // Compute tangents of LoS_vec and update output
            const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
            const double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
            
            if(this_tan > max_tan){
              max_tan = this_tan;
              viewshed[los_ref_cell] = cell+1;
            }
          }
          max_tan_vec[j] = max_tan;
        } else {
          break;
        }
      }
    }
  }
  
  return Rcpp::na_omit(viewshed);
}


#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
Rcpp::IntegerMatrix viewshed_distance_analysis_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                                   const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                                   const int radius, const Rcpp::NumericVector & h0,
                                                   const int ncores=1, const bool display_progress=false)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  const int s = input_cells.size();
  
  // Output
  Rcpp::IntegerMatrix output( r*s,3 );
  for(int k = 0; k < s; k++) {
    for(int d = 0; d<r; d++) {
      output( (d*s + k),0 ) = (d+1) * (int)ras.res;
    }
  }
  //for(int d = 0; d<r; d++) {
  //  output( d,0 ) = d+1;
  //}
  
  
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  
  // Progress bar
  Progress p(s, display_progress);
  
  // Main loop
  for(int k = 0; k < s; k++)
  {
    if(h0[k] > dsm_values[input_cells[k]]){
      // Viewshed Analysis
      const int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
#if defined(_OPENMP)
      omp_set_num_threads(ncores);
#pragma omp parallel for shared(output)
#endif    
      for(int i = 0; i < (r*8); i++){
        if (!p.is_aborted()) {
          Progress::check_abort();
          
          double max_tan = -9999.0;
          
          for(int j = 0; j < r; j++){
            const int los_ref_cell = los_ref_vec[i*r + j];
            if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
              // Project reference cells to actual cell value
              int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
              int row = trunc(cell/ras.ncol);
              int col = cell - (row * ras.ncol);
              int dcol = abs(col-x0_o[k]);
              
              if(!(cell < 0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
                // Compute tangents of LoS_vec and update output
                double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
                double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
                
                int d = round(distance_traveled) - 1;
                d = (d < 0) ? 0 : d;
                d = (d >= r) ? (r-1) : d;
                
                output( (d*s + k),1 ) += 1;
                //output( d,1 ) += 1;
                
                if(this_tan > max_tan){
                  max_tan = this_tan;
                  
                  output( (d*s + k),2 ) += 1;
                  //output( d,2 ) += 1;
                }
              }
            } else {
              break;
            }
          }
        }
      }
    }
    p.increment();
  }
  
  return output;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix viewshed_and_greenness_distance_analysis_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                                                 Rcpp::S4 &greenspace, const Rcpp::NumericVector &greenspace_values,
                                                                 const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                                                 const int radius, const Rcpp::NumericVector & h0,
                                                                 const int ncores=1, const bool display_progress=false)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  const RasterInfo gs_ras(greenspace);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  const int s = input_cells.size();
  
  // Output
  // v1: Distance
  // V2: Total number of cells per distance
  // V3: Number of all visible cells per distance
  // V4: Total number of green cells per distance
  // V5: Number of green visible cells per distance
  Rcpp::IntegerMatrix output( r*s,5 );
  for(int k = 0; k < s; k++) {
    for(int d = 0; d<r; d++) {
      output( (d*s + k),0 ) = (d+1) * (int)ras.res;
    }
  }
  
  
  
  // Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  
  // Progress bar
  Progress p(s, display_progress);
  
  // Main loop
#if defined(_OPENMP)
  omp_set_num_threads(ncores);
#pragma omp parallel for shared(output)
#endif
  for(int k = 0; k < s; k++)
  {
    if(h0[k] > dsm_values[input_cells[k]]){
      // Viewshed Analysis
      const int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
      for(int i = 0; i < (r*8); i++){
        if (!p.is_aborted()) {
          Progress::check_abort();
          
          double max_tan = -9999.0;
          
          for(int j = 0; j < r; j++){
            const int los_ref_cell = los_ref_vec[i*r + j];
            if(los_ref_cell != NA_INTEGER){ //!Rcpp::IntegerVector::is_na(los_ref_cell)
              // Project reference cells to actual cell value
              int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
              
              // row/col
              int row = trunc(cell/ras.ncol);
              int col = cell - (row * ras.ncol);
              int dcol = abs(col-x0_o[k]);
              
              // X/Y coordinates
              double x_j = ras.xmin + (col + 0.5) * ras.res;
              double y_j = ras.ymax - (row + 0.5) * ras.res;
              
              // Cell of greenspace raster
              // OpenMP can't call cellFromXY2 function without causing R to crash.
              // Therefore, I'll calcualte the cell here manually
              int cell_gs = cellFromXY2(gs_ras, x_j, y_j);
              
              
              
              
              if(!(cell < 0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){ 
                // Compute tangents of LoS_vec and update output
                double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
                double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
                
                int d = round(distance_traveled) - 1;
                d = (d < 0) ? 0 : d;
                d = (d >= r) ? (r-1) : d;
                
                output( (d*s + k),1 ) += 1;
                if(cell_gs != NA_INTEGER && greenspace_values[cell_gs] == 1){
                  output( (d*s + k),3 ) += 1;  
                }
                
                if(this_tan > max_tan){
                  max_tan = this_tan;
                  
                  output( (d*s + k),2 ) += 1;
                  if(cell_gs != NA_INTEGER && greenspace_values[cell_gs] == 1){
                    output( (d*s + k),4 ) += 1;  
                  }
                }
              }
            } else {
              break;
            }
          }
        }
      }
    }
    p.increment();
  }
  
  return output;
}

// Alternative version of viewshed. Slightly faster for smaller radius, but more RAM intensive.
Rcpp::IntegerVector viewshed_cpp2(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                                  const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                                  const int radius, const Rcpp::NumericVector & h0)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int l = r+1;
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  int k=0;
  
  
  // 1. Reference Bresenham Lines for the first eigth of circle
  const Rcpp::IntegerMatrix bh_mat = bresenham_map(x0_ref, y0_ref, r, nc_ref);
  
  // 2. Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells 
  Rcpp::IntegerMatrix los_ref_mat = Rcpp::IntegerMatrix(l*8 - 8, r);
  
  
  for(int i=0; i<l; i++){
    const Rcpp::IntegerMatrix bs_xy = colRowFromCell(bh_mat( i,_ ), nc_ref);
    
    // Fill los_ref_mat with projected Bresenham Lines, using Bresenham's Midpoint algorithm
    for(int j=0; j<r; j++){
      
      if ( Rcpp::IntegerVector::is_na(bs_xy( j,0 )) || Rcpp::IntegerVector::is_na(bs_xy( j,1 )) ) {
        los_ref_mat( i,j )             = NA_INTEGER; 
        los_ref_mat( 2*(l-1) + i,j )   = NA_INTEGER;
        los_ref_mat( 4*(l-1) + i,j )   = NA_INTEGER;
        los_ref_mat( 6*(l-1) + i,j )   = NA_INTEGER;
        
        if( i != 0 && i != (l-1) ) {
          los_ref_mat( 2*(l-1) - i,j ) = NA_INTEGER;
          los_ref_mat( 4*(l-1) - i,j ) = NA_INTEGER;
          los_ref_mat( 6*(l-1) - i,j ) = NA_INTEGER;
          los_ref_mat( 8*(l-1) - i,j ) = NA_INTEGER;
        }
      } else {
        const int x = bs_xy( j,0 )-x0_ref;
        const int y = bs_xy( j,1 )-x0_ref;
        los_ref_mat( i,j )             = (y+y0_ref)*nc_ref + (x+x0_ref); 
        los_ref_mat( 2*(l-1) + i,j )   = (-x+y0_ref)*nc_ref + (y+x0_ref);
        los_ref_mat( 4*(l-1) + i,j )   = (-y+y0_ref)*nc_ref + (-x+x0_ref);
        los_ref_mat( 6*(l-1) + i,j )   = (x+y0_ref)*nc_ref + (-y+x0_ref);
        
        if( i != 0 && i != (l-1) ) {
          los_ref_mat( 2*(l-1) - i,j ) = (x+y0_ref)*nc_ref + (y+x0_ref);
          los_ref_mat( 4*(l-1) - i,j ) = (-y+y0_ref)*nc_ref + (x+x0_ref);
          los_ref_mat( 6*(l-1) - i,j ) = (-x+y0_ref)*nc_ref + (-y+x0_ref);
          los_ref_mat( 8*(l-1) - i,j ) = (y+y0_ref)*nc_ref + (-x+x0_ref);
        }
      }
    }
  }
  
  
  // A: Tangents Map
  const Rcpp::NumericVector tan_vec = tangentsMap(x0_o[k], y0_o[k], h0[k], ras.nrow, ras.ncol, r, dsm_values);
  
  // B: Visibility
  Rcpp::LogicalVector visibility_vec(nc_ref*nr_ref); 
  
  for(int j = 0; j < (l*8 - 8); j++){ 
    
    // Compute tangents of LoS_vec and update visibility_vec
    LoS(los_ref_mat, tan_vec, visibility_vec, j, r);
  }
  
  
  // Project reference cells to actual cell value
  Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
  viewshed[c0_ref] = input_cells[k]+1;
  int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
  for(int j = 0; j < visibility_vec.size(); j++){
    if(visibility_vec[j]){
      const int cell = x + j + trunc(j/nc_ref)*(ras.ncol-nc_ref);
      const int row = cell - (trunc(cell/ras.ncol) * ras.ncol);
      const int drow = abs(row-x0_o[k]);
      
      viewshed[j] = (cell<0 || cell > ras.ncell ||
        Rcpp::NumericVector::is_na(dsm_values[cell]) ||
        drow>r) ? NA_INTEGER : (cell+1);
    }
  }
  
  return Rcpp::na_omit(viewshed);
}



Rcpp::NumericVector VGVI_cpp2(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                              Rcpp::S4 &greenspace, const Rcpp::NumericVector &greenspace_values,
                              const Rcpp::IntegerVector x0, const Rcpp::IntegerVector y0,
                              const int radius, const Rcpp::NumericVector & h0,
                              const int fun, const double m, const double b,
                              const int ncores=1)
{
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int l = r+1;
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  Rcpp::NumericVector output(input_cells.size(), NA_REAL);
  
  // 1. Reference Bresenham Lines for the first eigth of circle
  const Rcpp::IntegerMatrix bh_mat = bresenham_map(x0_ref, y0_ref, r, nc_ref);
  
  // 2. Reference matrix (los_mat): Project reference Bresenham Lines to all perimeter cells
  // Will be used as a reference for all input points
  Rcpp::IntegerVector los_ref_vec = Rcpp::IntegerVector(r*8 * r);
  
  for(int i=0; i<l; i++){
    const Rcpp::IntegerMatrix bs_xy = colRowFromCell(bh_mat( i,_ ), nc_ref);
    
    for(int j=0; j<r; j++){
      if ( Rcpp::IntegerVector::is_na(bs_xy( j,0 )) || Rcpp::IntegerVector::is_na(bs_xy( j,1 )) ) {
        los_ref_vec[(0*r+i)*r + j] = NA_INTEGER;
        los_ref_vec[(2*r+i)*r + j] = NA_INTEGER;
        los_ref_vec[(4*r+i)*r + j] = NA_INTEGER;
        los_ref_vec[(6*r+i)*r + j] = NA_INTEGER;
        
        if( i != 0 && i != (l-1) ) {
          los_ref_vec[(2*r-i)*r + j] = NA_INTEGER;
          los_ref_vec[(4*r-i)*r + j] = NA_INTEGER;
          los_ref_vec[(6*r-i)*r + j] = NA_INTEGER;
          los_ref_vec[(8*r-i)*r + j] = NA_INTEGER;
        }
      } else {
        const int x = bs_xy( j,0 )-x0_ref;
        const int y = bs_xy( j,1 )-x0_ref;
        los_ref_vec[(0*r+i)*r + j] = (y+y0_ref)*nc_ref + (x+x0_ref);
        los_ref_vec[(2*r+i)*r + j] = (-x+y0_ref)*nc_ref + (y+x0_ref);
        los_ref_vec[(4*r+i)*r + j] = (-y+y0_ref)*nc_ref + (-x+x0_ref);
        los_ref_vec[(6*r+i)*r + j] = (x+y0_ref)*nc_ref + (-y+x0_ref);
        
        if( i != 0 && i != (l-1) ) {
          los_ref_vec[(2*r-i)*r + j] = (x+y0_ref)*nc_ref + (y+x0_ref);
          los_ref_vec[(4*r-i)*r + j] = (-y+y0_ref)*nc_ref + (x+x0_ref);
          los_ref_vec[(6*r-i)*r + j] = (-x+y0_ref)*nc_ref + (-y+x0_ref);
          los_ref_vec[(8*r-i)*r + j] = (y+y0_ref)*nc_ref + (-x+x0_ref);
        }
        
      }
    }
  }
  
  // 3. Main loop
  for(int k=0; k<input_cells.size(); k++)
  {
    // Outpur: x0/y0 is always visible
    Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
    viewshed[c0_ref] = input_cells[k]+1;
    
    // A. Viewshed Analysis
    {
      int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      
#if defined(_OPENMP)
      omp_set_num_threads(ncores);
#pragma omp parallel for shared(viewshed, los_ref_vec)
#endif  
      for(int i = 0; i < (r*8); i++){
        
        double max_tan = -9999.0;
        
        for(int j = 0; j < r; j++){
          const int los_ref_cell = los_ref_vec[i*r + j];
          if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
            // Project reference cells to actual cell value
            const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
            const int row = trunc(cell/ras.ncol);
            const int col = cell - (row * ras.ncol);
            const int dcol = abs(col-x0_o[k]);
            
            if(!(cell<0 || cell > ras.ncell || Rcpp::NumericVector::is_na(dsm_values[cell]) || dcol>r)){
              // Compute tangents of LoS_vec and update output
              const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
              const double this_tan = (dsm_values[cell] - h0[k]) / (distance_traveled);
              
              if(this_tan > max_tan){
                max_tan = this_tan;
                viewshed[los_ref_cell] = cell+1;
              }
            }
          } else {
            break;
          }
        }
      }
    }
    Rcpp::IntegerVector viewshed_na = Rcpp::na_omit(viewshed);
    
    
    // B: Greenspace Visibility Index
    const int cs = viewshed_na.size();
    
    // Get XY coordinates of visible cells
    const Rcpp::NumericMatrix dsm_xy = xyFromCell(dsm, viewshed_na);
    const Rcpp::NumericMatrix xy0 = xyFromCell(dsm, input_cells[k]);
    
    // Calculate distance
    Rcpp::IntegerVector dxy(cs);
    for(int i = 0; i < cs; i++){
      double d = sqrt( ((xy0(0,0) - dsm_xy(i,0))*(xy0(0,0) - dsm_xy(i,0))) + ((xy0(0,1) - dsm_xy(i,1))*(xy0(0,1) - dsm_xy(i,1))) );
      int di = round(d);
      if(di < 1) {
        dxy[i] = 1;
      } else {
        dxy[i] = di;
      }
    }
    
    
    // Intersect XY with greenspace mask
    const Rcpp::IntegerVector greenspace_cells = cellFromXY(greenspace, dsm_xy);
    Rcpp::NumericVector greenspace_cell_values(cs, 0.0);
    for(int i = 0; i < cs; i++){
      int gs_cell = greenspace_cells[i]; 
      if(!Rcpp::IntegerVector::is_na(gs_cell)){
        double gs = greenspace_values[gs_cell];
        greenspace_cell_values[i] =  Rcpp::NumericVector::is_na(gs) ? 0.0 : gs;
      }
    }
    
    // Get number of green visible cells and total visible cells per distance
    const double max_d = max(dxy);
    const Rcpp::IntegerVector dxy_seq = seq(1, max_d);
    const int n = dxy_seq.size();
    Rcpp::IntegerVector visibleTotal(n);
    Rcpp::IntegerVector visibleGreen(n);
    
    for(int i = 0; i < cs; i++){
      int this_dxy = dxy[i];
      visibleTotal[this_dxy-1] += 1;
      visibleGreen[this_dxy-1] += greenspace_cell_values[i];
    }
    for(int i = 0; i < n; i++){
      if(visibleTotal[i] == 0){
        visibleTotal[i] = 1;
      }
    }
    
    
    // Proportion of visible green cells
    if(max_d == 1){
      output[k] = visibleGreen[0]/visibleTotal[0];
    }
    Rcpp::NumericVector raw_GVI(n);
    for(int i = 0; i < n; i++){
      raw_GVI[i] = (double)visibleGreen[i] / visibleTotal[i];
    }
    
    // Normalize distance
    Rcpp::NumericVector nDxy(n);
    for(int i = 0; i < n; i++){
      nDxy[i] = dxy_seq[i] / (double)radius; //max_d; // TAKE radius instead and save max_d / r for upper limit
    }
    
    // Calculate weights by taking the proportion of the integral of each step from the integral of the whole area
    //double big_integral = integrate(0, 1, 200, fun, m, b);
    const double min_dist = min(nDxy);
    Rcpp::NumericVector decayWeights(n);
    for(int i = 0; i < n; i++){
      double d = nDxy[i];
      decayWeights[i] = integrate(d-min_dist, d, 200, fun, m, b);
    }
    const double big_integral = sum(decayWeights);
    
    // Proportion of visible green
    double vgvi_sum = 0.0; 
    for(int i = 0; i < n; i++){
      vgvi_sum += raw_GVI[i] * (decayWeights[i]/big_integral);
    }
    
    output[k] = vgvi_sum;
  }
  return output;
}




