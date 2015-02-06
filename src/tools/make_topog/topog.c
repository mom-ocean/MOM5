#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "create_xgrid.h"
#include "mosaic_util.h"
#include "interp.h"
#include "mpp_io.h"
#include "mpp.h"
#include "mpp_domain.h"
#include "topog.h"

#define D2R (M_PI/180.)

const double deg2metre=111.324e3;
void filter_topo( int nx, int ny, int num_pass, int smooth_topo_allow_deepening, double *depth, domain2D domain);
void set_depth(int nx, int ny, const double *xbnd, const double *ybnd, double alat1, double slon1,
	       double elon1, double alat2, double slon2, double elon2, double depth_in, double *depth);
void show_deepest(int nk, const double *zw, const double *depth, domain2D domain );
void enforce_min_depth(int nx, int ny, int nk, const double *zw, int round_shallow, int deepen_shallow,
		       int fill_shallow, double *depth, int *kmt, int kmt_min, int debug);
void restrict_partial_cells(int ni, int nj, int nk, const double *zw, double *ht, int *kmt,
			    double min_thickness, int open_very_this_cell, double fraction_full_cell, int debug);
void remove_isolated_cells(int ni, int nj, double *ht, int *kmt, int tripolar_grid, int fill_isolated_cells,
			   int dont_change_landmask, int debug);
void process_topo(int nk, double *depth, int *num_levels, const double *zw,
		  int tripolar_grid, int cyclic_x, int cyclic_y, int full_cell,
		  int flat_bottom, int adjust_topo, int fill_isolated_cells,
		  double min_thickness, int open_very_this_cell, double fraction_full_cell, int round_shallow,
		  int deepen_shallow, int fill_shallow, int dont_change_landmask, int kmt_min, 
		  domain2D domain, int debug );
int kmin(int ni, int *kmt, int i, int j);
double hmin(int ni, double *ht, int i, int j);
double max4_double(double d1, double d2, double d3, double d4);
double min4_double(double d1, double d2, double d3, double d4);
int max4_int(int d1, int d2, int d3, int d4);
int min4_int(int d1, int d2, int d3, int d4);
/*********************************************************************
   void create_rectangular_topog()
   Constructing a rectangular basin with a flat bottom, basin_depth can be 0 ( all land).
 ********************************************************************/
void create_rectangular_topog(int nx, int ny, double basin_depth, double *depth)
{
  int i;

  for(i=0; i<nx*ny; i++)  depth[i] = basin_depth;
}; /* create_rectangular_topog  */

/*********************************************************************
   void create_bowl_topog()
   From "Simulation of density-driven frictional downslope flow
   in z-coordinate mocean models"
   Winton et al. JPO, Vol 28, No 11, 2163-2174, November 1998
 ********************************************************************/
void create_bowl_topog(int nx, int ny, const double *x, const double *y, double bottom_depth,
		       double min_depth, double bowl_east, double bowl_south, double bowl_west,
		       double bowl_north, double *depth)
{
  double bottom, xx, yy;
  int i, j, nxp, nyp;

  nxp = nx+1;
  nyp = ny+1;
  
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      xx = (x[j*nxp+i]+x[j*nxp+i+1]+x[(j+1)*nxp+i]+x[(j+1)*nxp+i+1])*0.25;
      yy = (y[j*nxp+i]+y[j*nxp+i+1]+y[(j+1)*nxp+i]+y[(j+1)*nxp+i+1])*0.25;	    
      if (xx <= bowl_west  || xx >= bowl_east || yy <= bowl_south || yy >= bowl_north)            
	bottom = min_depth;
      else
	bottom = min_depth + bottom_depth*(1.0-exp(-pow((yy-bowl_south)/2.0,2)))
	  *(1.0-exp(-pow((yy-bowl_north)/2.0,2)))
	  *(1.0-exp(-pow((xx-bowl_west)/4.0,2)))
	  *(1.0-exp(-pow((xx-bowl_east)/4.0,2)));
      depth[j*nx+i] = bottom;
    }
  }
}; /* create_bowl_topog */


/*********************************************************************
   void create_gaussian_topog()
   Constructing a gaussian bump
 ********************************************************************/
void create_gaussian_topog(int nx, int ny, const double *x, const double *y, double bottom_depth,
			   double min_depth, double gauss_amp, double gauss_scale, double slope_x,
			   double slope_y, double *depth)
{
  double bump_height, bump_scale, xcent, ycent, arg, bottom;
  double xe, xw, ys, yn;
  double *xt, *yt;
  int i, j, nxp, nyp;

  xt = (double *)malloc(nx*ny*sizeof(double));
  yt = (double *)malloc(nx*ny*sizeof(double));
  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++) {
      xt[j*nx+i] = (x[j*nxp+i] + x[j*nxp+i+1] + x[(j+1)*nxp+i] + x[(j+1)*nxp+i+1])*0.25;
      yt[j*nx+i] = (y[j*nxp+i] + y[j*nxp+i+1] + y[(j+1)*nxp+i] + y[(j+1)*nxp+i+1])*0.25;
    }
  }
  
  xw = xt[0];
  ys = yt[0];
  xe = xw;
  yn = ys;

  for(j=1;j<ny;j++){
    for(i=1;i<nx;i++){
      xw = min(xt[j*nx+i],xw); xe = max(xt[j*nx+i],xe);
      ys = min(yt[j*nx+i],ys); yn = max(yt[j*nx+i],yn);
    }
  }

  bump_height = gauss_amp*bottom_depth;
  bump_scale = gauss_scale*min(xe-xw, yn-ys);
  xcent = 0.5*(xe+xw);
  ycent = 0.5*(yn+ys);

  printf("Constructing a gaussian bump of height = %f meters with a scale width of %f degrees.\n",
	 bump_height, bump_scale);
  printf("The bump is centered at (lon,lat) = (%f,%f) deg.\n", xcent, ycent);
  printf("The ocean floor rises with a slope of %f meters/deg towards the east and %f "
	 " meters/deg to the north.", slope_x, slope_y) ;

  for(j=0;i<nx*ny;j++){
      arg = pow(xt[i]-xcent,2) + pow(yt[i]-ycent,2);
      bottom = bottom_depth - bump_height*exp(-arg/pow(bump_scale,2));
      bottom = bottom - slope_x*(xt[i]-xw)- slope_y*(yt[i]-ys);
      depth[i] = max(bottom,min_depth);
  }

  free(xt);
  free(yt);
  
}; /* create_gaussian_topog */


/*********************************************************************

  void create_idealized_topog
   construct a highly "idealized" world ... piece by piece

   note: the purpose of this geometry/topography is to automatically
         map into arbitrary resolution as grid dimensions are
         changed, thereby facilitating the implementation
         and verification of the model on various computer platforms
         without referencing the topographic database.  Although it
         somewhat resembles the real world, it is NOT realistic.
   Note: this routine needs to be re-thought for generalized curvilinear coordinates

 ********************************************************************/
void create_idealized_topog( int nx, int ny, const double *x, const double *y,
			     double bottom_depth, double min_depth, double *depth)
{
  int i, j, nxp, nyp;
  double arg;
  double *xbnd, *ybnd;

  nxp = nx+1;
  nyp = ny+1;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++) depth[j*nx+i]=bottom_depth;
  }
  
  //antarctica
  xbnd = (double *)malloc(nx*sizeof(double));
  ybnd = (double *)malloc(ny*sizeof(double));
  for(i=0; i<nx; i++) xbnd[i] = (x[i]+x[i+1])*0.5;
  for(j=0; j<ny; j++) ybnd[j] = (y[j*nxp]+y[(j+1)*nxp])*0.5;
  
  set_depth (nx, ny, xbnd, ybnd, -90.0, 0.0, 360.0, -80.0, 0.0, 360.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -80.0, 360.0-25.0, 360.0, -70.0, 360.0, 360.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -80.0, 0.0, 360.0, -70.0, 0.0, 170.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -80.0, 360.0-135.0, 360.0-60.0, -68.0, 360.0-75.0, 360.0-60.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -70.0, 0.0, 155.0, -67.0, 50.0, 145.0, 0.0, depth);

  // australia

  set_depth (nx, ny, xbnd, ybnd, -35.0, 116.0, 120.0, -31.0, 114.0, 130.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -38.0, 140.0, 151.0, -31.0, 130.0, 151.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -31.0, 115.0, 153.0, -20.0, 113.0, 149.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -20.0, 113.0, 149.0, -11.0, 131.0, 143.0, 0.0, depth);

  // south america

  set_depth (nx, ny, xbnd, ybnd, -50.0, 360.0-74.0, 360.0-68.0, -40.0, 360.0-73.0, 360.0-62.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -40.0, 360.0-73.0, 360.0-62.0, -20.0, 360.0-70.0, 360.0-40.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -20.0, 360.0-70.0, 360.0-40.0, -16.0, 360.0-81.0, 360.0-35.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, -16.0, 360.0-81.0, 360.0-35.0, 0.0, 360.0-80.0, 360.0-50.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 0.0, 360.0-80.0, 360.0-50.0, 11.0, 360.0-75.0, 360.0-60.0, 0.0, depth);

  // central america

  set_depth (nx, ny, xbnd, ybnd, 6.0, 360.0-78.0, 360.0-75.0, 20.0, 360.0-105.0, 360.0-97.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 20.0, 360.0-105.0, 360.0-97.0, 30.0, 360.0-115.0, 360.0-94.0, 0.0, depth);

  // north america

  set_depth (nx, ny, xbnd, ybnd, 25.0, 360.0-82.0, 360.0-80.0, 30.0, 360.0-85.0, 360.0-81.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 30.0, 360.0-115.0, 360.0-80.0, 40.0, 360.0-124.0, 360.0-74.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 40.0, 360.0-124.0, 360.0-74.0, 50.0, 360.0-124.0, 360.0-57.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 50.0, 360.0-124.0, 360.0-57.0, 60.0, 360.0-140.0, 360.0-64.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 60.0, 360.0-165.0, 360.0-64.0, 65.0, 360.0-140.0, 360.0-64.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 65.0, 360.0-140.0, 360.0-64.0, 70.0, 360.0-162.0, 360.0-72.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 70.0, 360.0-162.0, 360.0-140.0, 72.0, 360.0-157.0, 360.0-157.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 70.0, 360.0-130.0, 360.0-70.0, 75.0, 360.0-120.0, 360.0-80.0, 0.0, depth);

  // greenland

  set_depth (nx, ny, xbnd, ybnd, 60.0, 360.0-45.0, 360.0-45.0, 75.0, 360.0-58.0, 360.0-19.0, 0.0, depth);

  // africa

  set_depth (nx, ny, xbnd, ybnd, -35.0, 19.0, 28.0, 6.0, 8.0, 50.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 6.0, 0.0, 50.0, 18.0, 0.0, 56.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 18.0, 0.0, 56.0, 26.0, 0.0, 59.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 6.0, 360.0-10.0, 360.0, 18.0, 360.0-18.0, 360.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 18.0, 360.0-18.0, 360.0, 26.0,  360.0-15.0, 360.0, 0.0, depth);

  // northern africa and europe and asia

  set_depth (nx, ny, xbnd, ybnd, 26.0, 360.0-15.0, 360.0, 40.0, 360.0-7.0, 360.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 40.0, 360.0-7.0, 360.0, 50.0, 360.0, 360.0, 0.0, depth);

  set_depth (nx, ny, xbnd, ybnd, 8.0, 77.0, 78.0, 26.0, 65.0, 90.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 4.0, 99.0, 100.0, 26.0, 90.0, 115.0, 0.0, depth);

  set_depth (nx, ny, xbnd, ybnd, 26.0, 0.0, 126.0, 40.0, 0.0, 122.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 40.0, 0.0, 130.0, 50.0, 0.0, 140.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 50.0, 0.0, 140.0, 60.0, 8.0, 140.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 60.0, 8.0, 163.0, 65.0, 13.0, 180.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 65.0, 13.0, 188.0, 70.0, 20.0, 180.0, 0.0, depth);
  set_depth (nx, ny, xbnd, ybnd, 70.0, 70.0, 180.0, 75.0, 90.0, 100.0, 0.0, depth);

  // add an "idealized" undulating topography

  nxp = nx+1;
  nyp = ny+1;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      if (depth[j*nx+i] > 0.0 ) {
	arg = bottom_depth*(1-0.4*fabs(cos(((j+1)*M_PI)/nyp)*sin(((i+1)*2*M_PI)/nxp)));
        arg = max(arg, min_depth);      
        depth[j*nx+i] = arg;
      }
    }
  }

  // add "idealized" ridges

  arg = 0.666*bottom_depth;
  
  set_depth (nx, ny, xbnd, ybnd, -20.0, 360.0-20.0, 360.0-10.0, 30.0, 360.0-45.0, 360.0-35.0, arg, depth);
  set_depth (nx, ny, xbnd, ybnd, 30.0, 360.0-45.0, 360.0-35.0, 60.0, 360.0-20.0,  360.0-30.0, arg, depth);
  set_depth (nx, ny, xbnd, ybnd, -60.0,360.0-100.0, 360.0-130.0, 40.0, 360.0-160.0, 180.0, arg, depth);

  arg = 0.5*bottom_depth;
  
  set_depth (nx, ny, xbnd, ybnd, -50.0, 360.0-120.0, 360.0-120.0, 30.0, 190.0, 190.0, arg, depth);

  free(xbnd);
  free(ybnd);
  
} /* create_idealized_topog */

/*********************************************************************
   void set_depth(double alat1, double slon1, double elon1, double alat2,
                  double slon2, double elon2, double depth_in )
        set the topography depth "depth[i][j](i,j)" = "depth_in" within the area of
        the trapezoid bounded by vertices:
        (alat1,slon1), (alat1,elon1), (alat2,slon2), and (alat2,elon2)
   
        inputs:
   
        alat1 = southern latitude of trapezoid (degrees)
        slon1 = starting longitude of southern edge of trapezoid (deg)
        elon1 = ending longitude of southern edge of trapezoid (deg)
        alat2 = northern latitude of trapezoid (degrees)
        slon2 = starting longitude of northern edge of trapezoid (deg)
        elon2 = ending longitude of northern edge of trapezoid (deg)
        depth_in = constant depth value
   
 ********************************************************************/
void set_depth(int nx, int ny, const double *xbnd, const double *ybnd, double alat1, double slon1, double elon1, double alat2,
                      double slon2, double elon2, double depth_in, double *depth)
{
  double rdj, d;
  int i, j, i1, i2, j1, j2, is, ie, js, je, is1, ie1, is2, ie2;
  
  j1 = nearest_index(alat1, ybnd, ny );
  j2 = nearest_index(alat2, ybnd, ny );
  js = min(j1,j2);
  je = max(j1,j2);

  i1  = nearest_index(slon1, xbnd, nx);
  i2  = nearest_index(elon1, xbnd, nx);
  is1 = min(i1,i2);
  ie1 = max(i1,i2);

  i1  = nearest_index(slon2, xbnd, nx );
  i2  = nearest_index(elon2, xbnd, ny );
  is2 = min(i1,i2);
  ie2 = max(i1,i2);

  is = is1;
  ie = ie1;

  // fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
  // the nudging of 1.e-5 is to insure that the test case resolution
  // generates the same topography and geometry on various computers.

  if (js == je) 
    rdj = 1.0;
  else
    rdj = 1.0/(je-js);

  for(j=js;j<=je;j++){
    for(i=is;i<=ie;i++) depth[j*nx+i] = depth_in;
    d = rdj*((j-js)*is2 + (je-j)*is1) + 1.0e-5;
    is = ceil(d);
    if( is - d > 0.5) is = is -1;
    d = rdj*((j-js)*ie2 + (je-j)*ie1) + 1.0e-5;
    ie = ceil(d);
    if( ie - d > 0.5) ie = ie -1;
  }

}; /* set_depth */

/*********************************************************************
   void create_realistic_topog( )
   reading data from source data file topog_file and remap it onto current grid
 ********************************************************************/
void create_realistic_topog(int nx_dst, int ny_dst, const double *x_dst, const double *y_dst, const char *vgrid_file, 
			    const char* topog_file, const char* topog_field, double scale_factor,
			    int tripolar_grid, int cyclic_x, int cyclic_y, 
			    int fill_first_row, int filter_topog, int num_filter_pass,
			    int smooth_topo_allow_deepening, int round_shallow, int fill_shallow,
			    int deepen_shallow, int full_cell, int flat_bottom, int adjust_topo,
			    int fill_isolated_cells, int dont_change_landmask, int kmt_min, double min_thickness,
			    int open_very_this_cell, double fraction_full_cell, double *depth,
			    int *num_levels, domain2D domain, int debug, int use_great_circle_algorithm )
{
  char xname[128], yname[128];
  int nx_src, ny_src, nxp_src, nyp_src, i, j, count, n;
  double *depth_src, *mask_src, *xt_src, *yt_src, *x_src, *y_src;
  double *x_out, *y_out, *xc_src, *yc_src;
  double missing, y_min, y_max, yy;
  int    nzv, nk, k;
  double *zeta=NULL, *zw=NULL;
  int    jstart, jend, jj;
  int    fid, vid;
  int    ny_now;
  size_t start[4], nread[4];


  /* read the vertical grid when vgrid_file is defined */
  if( vgrid_file ) {
    fid = mpp_open(vgrid_file, MPP_READ);
    nzv = mpp_get_dimlen(fid, "nzv");
    if( (nzv-1)%2 ) mpp_error("topog: size of dimension nzv should be 2*nk+1, where nk is the number of model vertical level");
    nk = (nzv-1)/2;
    zw   = (double *)malloc(nk*sizeof(double));
    zeta = (double *)malloc(nzv*sizeof(double));
    vid = mpp_get_varid(fid, "zeta");
    mpp_get_var_value(fid, vid, zeta);
    mpp_close(fid);
    for(k=0; k<nk; k++) zw[k] = zeta[2*(k+1)];
    free(zeta);
  }
    
  /* first read source topography data to get source grid and source topography */
  fid = mpp_open(topog_file, MPP_READ);
  vid = mpp_get_varid(fid, topog_field);
  mpp_get_var_dimname(fid, vid, 1, xname);
  mpp_get_var_dimname(fid, vid, 0, yname);
  nx_src = mpp_get_dimlen( fid, xname );
  ny_src = mpp_get_dimlen( fid, yname );
  nxp_src = nx_src + 1;
  nyp_src = ny_src + 1;

  xt_src    = (double *)malloc(nx_src*sizeof(double));
  yt_src    = (double *)malloc(ny_src*sizeof(double));
  xc_src    = (double *)malloc(nxp_src*sizeof(double));
  yc_src    = (double *)malloc(nyp_src*sizeof(double));

  mpp_get_var_att(fid, vid, "missing_value", &missing);
  vid = mpp_get_varid(fid, xname);
  mpp_get_var_value(fid, vid, xt_src);
  vid = mpp_get_varid(fid, yname);
  mpp_get_var_value(fid, vid, yt_src);

  
  for(i=1; i<nx_src; i++) xc_src[i] = (xt_src[i-1] + xt_src[i])*0.5;
  xc_src[0] = 2*xt_src[0] - xc_src[1];
  xc_src[nx_src] = 2*xt_src[nx_src-1] - xc_src[nx_src-1];

  for(j=1; j<ny_src; j++) yc_src[j] = (yt_src[j-1] + yt_src[j])*0.5;
  yc_src[0] = 2*yt_src[0] - yc_src[1];
  yc_src[ny_src] = 2*yt_src[ny_src-1] - yc_src[ny_src-1];
    

  /* scale grid to radius */
  for(i=0; i<nxp_src; i++) xc_src[i] = xc_src[i]*D2R;
  for(j=0; j<nyp_src; j++) yc_src[j] = yc_src[j]*D2R;

  x_out = (double *)malloc((nx_dst+1)*(ny_dst+1)*sizeof(double));
  y_out = (double *)malloc((nx_dst+1)*(ny_dst+1)*sizeof(double));
  for(i=0; i<(nx_dst+1)*(ny_dst+1); i++) {
    x_out[i] = x_dst[i]*D2R;
    y_out[i] = y_dst[i]*D2R;
  }
  
  /* doing the conservative interpolation */
  y_min = minval_double((nx_dst+1)*(ny_dst+1), y_out);
  y_max = maxval_double((nx_dst+1)*(ny_dst+1), y_out);
  jstart = ny_src; jend = -1;
  for(j=0; j<=ny_src; j++){
    yy = yc_src[j];
    if( yy > y_min && yy < y_max) {
      if(j > jend  ) jend   = j;
      if(j < jstart) jstart = j;
    }
  }
  jstart = max(0, jstart-1);
  jend   = min(ny_src-1, jend+1);
  ny_now = jend-jstart+1;

  x_src    = (double *)malloc(nxp_src*(ny_now+1)*sizeof(double));
  y_src    = (double *)malloc(nxp_src*(ny_now+1)*sizeof(double));
  depth_src = (double *)malloc(nx_src*ny_now*sizeof(double));
  mask_src = (double *)malloc(nx_src*ny_now*sizeof(double));

  for(j=0; j<=ny_now; j++) {
     jj = j+jstart;
     for(i=0; i<nxp_src; i++) x_src[j*nxp_src+i] = xc_src[i];
     for(i=0; i<nxp_src; i++) y_src[j*nxp_src+i] = yc_src[jj];
  }


  for(i=0; i<4; i++) {
     start[i] = 0;
     nread[i] = 1;
  }
  start[0] = jstart;
  nread[0] = ny_now;
  nread[1] = nx_src;
  vid = mpp_get_varid(fid, topog_field);
  {
     nc_type vartype;
     int *depth_int=NULL;
     vartype = mpp_get_var_type(fid, vid);
     switch (vartype) { 
     case NC_DOUBLE:case NC_FLOAT:
        mpp_get_var_value_block(fid, vid, start, nread, depth_src);
        break;
     case NC_INT:
        depth_int = (int *)malloc(nx_src*ny_now*sizeof(int));
        mpp_get_var_value_block(fid, vid, start, nread, depth_int);
        for(i=0; i<nx_src*ny_now; i++) depth_src[i] = depth_int[i];
        free(depth_int);
        break;
     default:
        mpp_error("topog.c: nc_type should be NC_DOUBLE, NC_FLOAT or NC_INT");
     }
  }

  mpp_close(fid);

  for(i=0; i<nx_src*ny_now; i++) {
    if(depth_src[i] == missing)
      mask_src[i] = 0.0;
    else {
      depth_src[i] = depth_src[i]*scale_factor;
      if( depth_src[i] <= 0.0)
	mask_src[i] = 0.0;
      else
	mask_src[i] = 1.0;
    }
  }
  
  if(use_great_circle_algorithm)  
    conserve_interp_great_circle(nx_src, ny_now, nx_dst, ny_dst, x_src, y_src,
		    x_out, y_out, mask_src, depth_src, depth );
  else
    conserve_interp(nx_src, ny_now, nx_dst, ny_dst, x_src, y_src,
		    x_out, y_out, mask_src, depth_src, depth );
  
  if (filter_topog) filter_topo(nx_dst, ny_dst, num_filter_pass, smooth_topo_allow_deepening, depth, domain);
  if(debug) show_deepest(nk, zw, depth, domain);
  
  /* make first row of ocean model all land points for ice model */

  if(fill_first_row) {
    for(i=0;i<nx_dst;i++) {
      depth[i] = 0.0;
    }
  }
  if(vgrid_file) process_topo(nk, depth, num_levels, zw, tripolar_grid, cyclic_x, cyclic_y,
			      full_cell, flat_bottom, adjust_topo, fill_isolated_cells, min_thickness,
			      open_very_this_cell, fraction_full_cell, round_shallow, deepen_shallow, fill_shallow,
			      dont_change_landmask, kmt_min, domain, debug );
  
  free(depth_src);
  free(mask_src);
  free(x_src);
  free(y_src);
  free(xc_src);
  free(yc_src);
  free(xt_src);
  free(yt_src);
  free(x_out);
  free(y_out);
  
}; /* create_realistic_topog */

void process_topo(int nk, double *depth, int *num_levels, const double *zw,
		  int tripolar_grid, int cyclic_x, int cyclic_y, int full_cell,
		  int flat_bottom, int adjust_topo, int fill_isolated_cells,
		  double min_thickness, int open_very_this_cell, double fraction_full_cell, int round_shallow,
		  int deepen_shallow, int fill_shallow, int dont_change_landmask, int kmt_min, 
		  domain2D domain, int debug )
{
  int i,j, isc, iec, jsc, jec, nxc, nyc, nip, njp, ni, nj;
  double *ht=NULL, *tmp=NULL;
  int    *kmt=NULL;
  double ht_prev;

  /* first get the global data */
  mpp_get_compute_domain2d( domain, &isc, &iec, &jsc, &jec);
  mpp_get_global_domain2d( domain, &ni, &nj);
  nxc = iec-isc+1;
  nyc = jec-jsc+1;
  
  nip = ni + 2;
  njp = nj + 2;
  ht  = (double *)malloc(nip*njp*sizeof(double));
  kmt = (int *)malloc(nip*njp*sizeof(double));
  tmp = (double *)malloc(ni*nj*sizeof(double));
  for(i=0; i<nip*njp; i++) {
    kmt[i] = -1;
    ht[i] = 0;
  }
  mpp_global_field_all_double(domain, nxc, nyc, depth, tmp);
  for(j=1; j<=nj; j++) for(i=1; i<=ni; i++) ht[j*nip+i] = tmp[(j-1)*ni+i-1];
  free(tmp);
  
  for(j=1; j<=nj; j++) for(i=1; i<=ni; i++) {
    if( ht[j*nip+i] > 0 ) {
      kmt[j*nip+i] = nearest_index( ht[j*nip+i], zw, nk);
      if(zw[kmt[j*nip+i]] < ht[j*nip+i] ) {
	if(( ht[j*nip+i]-zw[kmt[j*nip+i]]) < min_thickness) {
	  ht_prev = ht[j*nip+i];
	  ht[j*nip+i] = zw[kmt[j*nip+i]];
	  if(debug) {
	    printf("Warning:  very thin partial cell at %d, %d (possibly because of netcdf-accuracy) is changed from %g to %g\n",
		   i, j, ht_prev, ht[j*nip+i]);
	    printf(" If this is not wanted, pass --min_thickness 0 when running make_topog\n");
	  }
	}
        else {
	  if(kmt[j*nip+i] == nk-1 )
	    ht[j*nip+i] = zw[kmt[j*nip+i]];
	  else
	    kmt[j*nip+i]++;
	}
      }
    }
  }

  if(full_cell) {
    printf("Warning from topog: Replacing partial bottom cells with full cell thicknesses.");
    for(j=1; j<=nj; j++) for(i=1; i<=ni; i++) {
      if(ht[j*nip+i] >0 ) ht[j*nip+i] = zw[kmt[j*nip+i]];
    }
  }

  if(flat_bottom) {
    printf("Warning from topog: Replacing the ocean topography with a flat bottom where kmt(i,j)=nk.");
    for(j=1; j<=nj; j++) for(i=1; i<=ni; i++) {
      if(ht[j*nip+i] > 0) {
	ht[j*nip+i] = zw[nk-1];
	kmt[j*nip+i]= nk-1;
      }
    }
  }

  /* fill the boundary condition */
  if(tripolar_grid) {
    for(i=1; i<=ni; i++) {
      kmt[(nj+1)*nip+i] = kmt[nj*nip+(ni-i+1)];
      ht [(nj+1)*nip+i] = ht [nj*nip+(ni-i+1)];
    }
  }
  else if(cyclic_y) {
    for(i=1; i<=ni; i++){
      kmt[i] = kmt[nj*nip+i];
      kmt[(nj+1)*nip+i] = kmt[nip+i];
      ht[i] = ht[nj*nip+i];
      ht[(nj+1)*nip+i] = ht[nip+i];
    }
  }

  if(cyclic_x) {
    for(j=0; j<=nj+1; j++) {
      kmt[j*nip] = kmt[j*nip+ni];
      kmt[j*nip+ni+1] = kmt[j*nip+1];
      ht [j*nip] = ht[j*nip+ni];
      ht [j*nip+ni+1] = ht[j*nip+1];
    }
  }

  if (adjust_topo) {
    remove_isolated_cells(ni, nj, ht, kmt, tripolar_grid, fill_isolated_cells, dont_change_landmask, debug);
  }

  /* copy global data to local data */
  for(j=jsc; j<=jec; j++) for(i=isc; i<=iec; i++) {
    depth[(j-jsc)*nxc+i-isc] = ht[(j+1)*nip+i+1];
    num_levels[(j-jsc)*nxc+i-isc] = kmt[(j+1)*nip+i+1];
  }
  free(ht);
  free(kmt);
  
  if (adjust_topo) {  
    restrict_partial_cells(nxc, nyc, nk, zw, depth, num_levels, min_thickness,
			   open_very_this_cell, fraction_full_cell, debug);
    enforce_min_depth(nxc, nyc, nk, zw, round_shallow, deepen_shallow, fill_shallow,
		      depth, num_levels, kmt_min, debug);
  }

  for(i=0; i<nxc*nyc; i++) num_levels[i]++;
  
  /* do some analysis. This will be added in the future if necessary
    if(debug) analyze_topographic_slopes(dxte, dytn, geolon_t, geolat_t, ht(1:ni,1:nj), kmt(1:ni,1:nj), nk) 
*/

}



/*********************************************************************

    void filter_topo( int nx, int ny, int num_pass, double *depth)

       smooth topographic depth "d" with "num_pass" applications of a 2D
       version of a shapiro filter (weights = 1/4, 1/2, 1/4) . 
       allow filtering to decrease the bottom depth but not increase it.
       do not allow original geometry to change.
       note: depth "d" should be on a grid of uniformly constant spacing

********************************************************************/
void filter_topo( int nx, int ny, int num_pass, int smooth_topo_allow_deepening, double *depth, domain2D domain)
{
  double *rmask, *s;
  double *depth_global;
  double f[9], d_old;
  int n1, n2, i, j, n, m, ip, jp;
  size_t *dims;
  int is, ie, js, je, nxg, nyg, pos;
  
  
  mpp_get_compute_domain2d(domain, &is, &ie, &js, &je);
  if( ie-is+1 != nx && je-js+1 !=ny ) mpp_error("topog.c: nx and ny does not match compute domain");
  mpp_get_global_domain2d(domain, &nxg, &nyg);
  depth_global = (double *)malloc(nxg*nyg*sizeof(double));
  mpp_global_field_all_double(domain, nx, ny, depth, depth_global);
  
  rmask = (double *)malloc(nxg*nyg*sizeof(double));
  s     = (double *)malloc(nxg*nyg*sizeof(double));

  /* 2D symmetric filter weights */

  f[0] = 1.0/16.0;
  f[1] = 1.0/8.0;
  f[2] = 1.0/16.0;
  f[3] = 1.0/8.0;
  f[4] = 1.0/4.0;
  f[5] = 1.0/8.0;
  f[6] = 1.0/16.0;
  f[7] = 1.0/8.0;
  f[8] = 1.0/16.0;
  
  /* geometry mask */
  for(i=0; i<nxg*nyg; i++) {
    s[i] = depth_global[i];
    if(depth_global[i] == 0.0)
      rmask[i] = 0.0;
    else
      rmask[i] = 1.0;
  }

  for(n=1;n<=num_pass;n++) {
    for(j=1;j<nyg-1;j++) {
      for(i=1;i<nxg-1;i++) {
	s[j*nxg+i] = 0.0;
	d_old = depth_global[j*nxg+i];
	for(ip=-1;ip<=1;ip++) {
	  for(jp=-1;jp<=1;jp++) {
	    m = (ip+1)*3 + jp+1;
	    if (rmask[(j+jp)*nxg+i+ip] == 0.0) 
	      s[j*nxg+i] = s[j*nxg+i] + depth_global[j*nxg+i]*f[m];
            else
	      s[j*nxg+i] = s[j*nxg+i] + depth_global[(j+jp)*nxg+i+ip]*f[m];
	  }
	}
	if (! smooth_topo_allow_deepening) {
	  if (s[j*nxg+i] > d_old)  s[j*nxg+i] = d_old;
	}
	s[j*nxg+i] = s[j*nxg+i]*rmask[j*nxg+i];
      }
    }
    for(i=0; i<nxg*nyg; i++) depth_global[i] = s[i];
  }

  /* copy data back to local data */
  pos = 0;
  for(j=js; j<=je; j++) for(i=is; i<=ie; i++) {
    depth[pos] = depth_global[j*nxg+i];
    pos++;
  }

  free(depth_global);
  free(rmask);
  free(s);
  
}; /* filter_topog */

int min4_int(int d1, int d2, int d3, int d4)
{
  int x, x1, x2;

  x1 = min(d1,d2);
  x2 = min(d3,d4);
  x  = min(x1,x2);

  return x;
}


double min4_double(double d1, double d2, double d3, double d4)
{
  double x, x1, x2;

  x1 = min(d1,d2);
  x2 = min(d3,d4);
  x  = min(x1,x2);

  return x;
}

int max4_int(int d1, int d2, int d3, int d4)
{
  int x, x1, x2;

  x1 = max(d1,d2);
  x2 = max(d3,d4);
  x  = max(x1,x2);

  return x;
}

double max4_double(double d1, double d2, double d3, double d4)
{
  double x, x1, x2;

  x1 = max(d1,d2);
  x2 = max(d3,d4);
  x  = max(x1,x2);

  return x;
}

/* get the minimum depth of surrounding cells
   it is assumed the data ht has halo = 1
*/
  double hmin(int ni, double *ht, int i, int j)
{
  int nip;
  nip = ni + 2;
  return min4_double(ht[j*nip+i], ht[j*nip+i+1], ht[(j+1)*nip+i], ht[(j+1)*nip+i+1]);
}  

/* get the minimum levels of surrounding cells
   it is assumed the data ht has halo = 1
*/
  int kmin(int ni, int *kmt, int i, int j)
{
  int nip;
  nip = ni + 2;
  return min4_double(kmt[j*nip+i], kmt[j*nip+i+1], kmt[(j+1)*nip+i], kmt[(j+1)*nip+i+1]);
}  


/*--- remove isolated cells */
void remove_isolated_cells(int ni, int nj, double *ht, int *kmt, int tripolar_grid, int fill_isolated_cells,
			   int dont_change_landmask, int debug)
{
  int i, j, k, n, nj_max, nip;
  double tmp;

  /* when tripolar_grid is true, do not check at j = nj (folded region).
    Will upgrade it if needed.
  */
  if(tripolar_grid)
    nj_max = nj-1;
  else
    nj_max = nj;
  nip = ni+2;
  
  if(fill_isolated_cells) {
    /*  fill isolated cells (potholes and trenches) at all levels in kmt.
        filled kmt array is the maximum of the surrounding kmt levels.
    */
    n=0;
    if(debug) printf(" Searching for isolated T cells...\n");
    for(j=1; j<=nj_max; j++) for(i=1; i<=ni; i++) {
      k = max4_int(kmin(ni, kmt,i,j), kmin(ni, kmt,i-1,j), kmin(ni, kmt,i,j-1), kmin(ni, kmt,i-1,j-1));
      if(k>=0 || !dont_change_landmask) {
	if(kmt[j*nip+i] != k) {
	  n++;
	  tmp = max4_double(hmin(ni, ht,i,j), hmin(ni, ht,i-1,j), hmin(ni, ht,i,j-1),hmin(ni, ht,i-1,j-1));
	  if(debug)printf("Resetting location (%d, %d), kmt from %d to %d, ht from %g to %g\n",
			  i, j, kmt[j*nip+i], k, ht[j*nip+i], tmp);
	  ht[j*nip+i] = tmp;
	  kmt[j*nip+i] = k;
	}
      }
    }

    if(debug) {
      if (n > 0)
	printf("remove_isolated_cells: found %d and filled them in\n", n);
      else
	printf("removed_isolated_cells: none found\n");
    }

  }
  else
    if(debug) printf("remove_isolated_cells: Not filling isolated T cells...\n"); 
}

/* Restricting partial cells, there will be no halo for data ht and kmt */
void restrict_partial_cells(int ni, int nj, int nk, const double *zw, double *ht, int *kmt,
			    double min_thickness, int open_very_this_cell, double fraction_full_cell, int debug)
{
  double *p_cell_min;
  double tmp, tmp1, tmp2;
  int    i, j, n, k, itmp;

  p_cell_min = (double *)malloc(nk*sizeof(double));

  n = 0;
  if(fraction_full_cell==0.0)
    for(k=0; k<nk; k++) p_cell_min[k] = min(50.0, zw[0]);
  else {
    p_cell_min[0] = fraction_full_cell*zw[0];
    for(k=1; k<nk; k++) p_cell_min[k] = max(fraction_full_cell*(zw[k]-zw[k-1]),min_thickness); 
  }
  if(debug) {
    printf("Partial cell restriction\n");
    for(k=0; k<nk; k++)
      printf("Restricting partial cell # %d to a min thickness of %g\n",k, p_cell_min[k]);
  }

  if( open_very_this_cell ) {
    for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
      k = kmt[j*ni+i];
      if(k > 0) {
	if((ht[j*ni+i]-zw[k-1]) < p_cell_min[k]) {
	  tmp = zw[k-1] + p_cell_min[k];
	  if(debug)printf("restrict_partial_cells: Resetting depth at (%d,%d) from %g to %g\n", i,j,ht[j*ni+i],tmp);
	  ht[j*ni+i] = tmp;
	  n++;
	}
      }
    }
  }
  else {
    for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
      k = kmt[j*ni+i];
      if(k > 0) {
	if((ht[j*ni+i]-zw[k-1]) < p_cell_min[k]) {
	  tmp1 = ht[j*ni+i]-zw[k-1];
	  tmp2 = zw[k-1] + p_cell_min[k] - ht[j*ni+i];
	  if (tmp2 > tmp1) {
	    itmp = kmt[j*ni+i] - 1; 
	    tmp = zw[k-1];
	    if(debug)printf("restrict_partial_cells: Resetting kmt at (%d,%d) from %d to %d\n", i,j,kmt[j*ni+i],itmp);
	    kmt[j*ni+i] = itmp;
	  }
	  else
	    tmp = zw[k-1] + p_cell_min[k];
	  if(debug)printf("restrict_partial_cells: Resetting depth at (%d,%d) from %g to %g\n", i,j,ht[j*ni+i],tmp);  
	  ht[j*ni+i] = tmp;
	  n++;
	}
      }
    }
  }

  if(debug) {
    if(n>0)
      printf("restrict_partial_cells: Found %d cells with too thin partical cell thickness, and so reset depth ht for these cells.\n",n );
    else
      printf("restrict_partial_cells: No cells were found whose partial cells were too thin.\n");
  }
}

/*--- print out deepest topography */
void show_deepest(int nk, const double *zw, const double *depth, domain2D domain )
{
  int  i, isc, iec, jsc, jec, nxc, nyc, ni, nj;
  double deepest;
  double *ht;
  
  mpp_get_compute_domain2d( domain, &isc, &iec, &jsc, &jec);
  mpp_get_global_domain2d( domain, &ni, &nj);
  nxc = iec-isc+1;
  nyc = jec-jsc+1;
  ht = (double *)malloc(ni*nj*sizeof(double));
  mpp_global_field_double(domain, nxc, nyc, depth, ht);
  
  if(mpp_pe() == mpp_root_pe()) {
    deepest = 0.0;
    for(i=0; i<ni*nj; i++) {
      if(ht[i] != 0) deepest = max(ht[i], deepest);
    }

    if((deepest - zw[nk-1]) > 1.0) {
      printf("Warning: Topography reaches to a depth of %g m\n", deepest);
      printf("         The deepest model level only reaches %g m\n", zw[nk-1]);
      printf("         Re-think the vertical grid specification if the idea \n");
      printf("         is to accurately capture the specified topography.\n");
    }
    else if( nk > 1 && deepest <= zw[nk-2] ) {
      printf("Warning: Topography reaches to a depth of %g m\n", deepest);
      printf("         The deepest model level reaches %g m\n", zw[nk-1]);
      printf("         Fewer than %d vertical levels are required.\n", nk);
      printf("         The specified number of vertical levels is wasteful.\n");
    }
    else {
      printf("Note: The vertical grid specification contains the correct number of \n");
      printf("      levels to efficiently contain the specified topography.\n");
    }
  }
}
    
/*********************************************************************
  void enforce_min_depth()
  ! limit the minimum value of depth, depth should be >= min_depth. 
 ********************************************************************/
void enforce_min_depth(int nx, int ny, int nk, const double *zw, int round_shallow, int deepen_shallow,
		       int fill_shallow, double *depth, int *kmt, int kmt_min, int debug)
{
  int i, j, n, l, count, kmt_shallow, kmt_min_local;
  char errmsg[256];
  
  
  if(debug) printf(" Enforcing the minimum number of ocean cells in the vertical to be %d\n", kmt_min);
  if(nk < kmt_min) {
    sprintf(errmsg, "topog: number of vertical cells=%d, is less than kmt_min=%d", nk, kmt_min);
    mpp_error(errmsg);
  }

  kmt_min_local = kmt_min-1;
  
  count = 0;
  if(round_shallow) count++;
  if(fill_shallow) count++;
  if(deepen_shallow) count++;
  if(count > 1) mpp_error("topog: at most one of round_shallow/deepen_shallow/fill_shallow can be set to true");
  
  n = 0;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      l = j*nx+i;
      if (depth[l] > 0 && kmt[l] < kmt_min_local) {
	n = n + 1;
	kmt_shallow = kmt[l];
	if( round_shallow || (!deepen_shallow && !fill_shallow) ) {    
	  if (zw[kmt[l]] < 0.5*zw[kmt_min_local]) {
	    if(debug) printf("Making location i,j,kmt= %d,%d,%d to land.\n",i, j, kmt[l]);
	    depth[l] = 0.0;
	    kmt[l] = -1;
	  }
	  else {
	    if(debug) printf("Setting location i,j,kmt= %d,%d,%d to minimum ocean depth.\n",i,j, kmt[l]);
	    depth[l] = zw[kmt_min_local];
	    kmt[l] = kmt_min_local;
	  }
	}
	if(fill_shallow) {
	  if (debug) printf("Setting location i,j,kmt= %d, %d, %d to land\n", i,j, kmt[l]);
	  depth[l] = 0.0;
	  kmt[l] = -1;
	}
	if(deepen_shallow) {
	  if(debug) printf("Setting location i,j,kmt= %d, %d, %d to minimum ocean depth\n", i,j, kmt[l]);
	  depth[l] = zw[kmt_min_local];
	  kmt[l] = kmt_min_local;
	}
      }
    }
  }
  
  if(debug) {
    if(n>0)
      printf("enforce_min_depth: Modified %d shallow cells\n", n);
    else
      printf("enforce_min_depth: No modifications needed\n");
  }
    
}; /* enforce_min_depth */


/********************************************************************************
  void create_box_channel_topog(int nx, int ny, double basin_depth,
				double jwest_south, double jwest_north, double jeast_south,
                                double jeast_north, double *depth)
  construct a box_channel topography.
********************************************************************************/
void create_box_channel_topog(int nx, int ny, double basin_depth,
			      double jwest_south, double jwest_north, double jeast_south,
			      double jeast_north, double *depth)
{
  int i, j;

  for(i=0; i<nx*ny; i++)  depth[i] = basin_depth;
  for(i=0; i<nx; i++) {
    depth[i] = 0;           /* southern boundary */
    depth[(ny-1)*nx+i] = 0; /* northern boundary */
  }
  /* west boundary */
  for(j=0; j<jwest_south; j++) depth[j*nx] = 0;
  for(j=jwest_north; j<ny; j++) depth[j*nx] = 0;
  /* east boundary */
  for(j=0; j<jeast_south; j++) depth[j*nx+nx-1] = 0;
  for(j=jeast_north; j<ny; j++) depth[j*nx+nx-1] = 0;
  
}

/******************************************************************************
void create_dome_topog(int nx, int ny, const double *x, const double *y, double dome_slope,
		       double dome_bottom, double dome_embayment_west, double dome_embayment_east,
		       double dome_embayment_south, double dome_embayment_depth, double *depth)

Construct a domain similar to the DOME configuration used in 
      Legg, Hallberg, and Girton (2005) Ocean Modelling.  
  
  
                         |     |
                         |     |
                        /|     |
                       / |     |
       up             /  |     |
       ^             /   |     |
       |            /    |     |
       |           /     |     |
                  /      |     |
  ---------------/-------|-----|
       --->north 
  
  
*******************************************************************************/
void create_dome_topog(int nx, int ny, const double *x, const double *y, double dome_slope,
		       double dome_bottom, double dome_embayment_west, double dome_embayment_east,
		       double dome_embayment_south, double dome_embayment_depth, double *depth)
{
  int i, j;
  double xw, xe, ys, yn;
  double xw_embay, xe_embay, yn_embay, ys_embay;
  double yn_slope, ys_slope, height;

  /* determine full domain extents */
  xw = x[0];
  ys = y[0];
  xe = xw;
  yn = ys;
  
  for(i=0; i<nx*ny; i++) {
    xw = min(x[i], xw);
    xe = max(x[i], xe);
    ys = min(y[i], ys);
    yn = max(y[i], yn);
  }
  
  /* define embayment location */
  xw_embay = dome_embayment_west;
  xe_embay = dome_embayment_east;
  yn_embay = yn;
  ys_embay = dome_embayment_south;

  /* find grid (lat,lon) corresponding to embayment specifications */
  for(i=0; i<nx-1; i++) {
    if(x[i] <= xw_embay && x[i+1] >= xw_embay) xw_embay = x[i];
    if(x[i] <= xe_embay && x[i+1] >= xe_embay) xe_embay = x[i];
  }
  for(j=0; j<ny-1; j++) {
    if(y[j*nx] <= ys_embay && y[(j+1)*nx] >= ys_embay ) ys_embay = y[j*nx];
  }

  /* define latitudes sloping surface starts and finishes */

  yn_slope = ys_embay;
  ys_slope = yn_slope - (dome_bottom-dome_embayment_depth)/(dome_slope*deg2metre);

  for(j=0; j<ny-1; j++) {
    if(y[j*nx] <= ys_slope && y[(j+1)*nx] >= ys_slope ) ys_slope = y[j*nx];
  }

  printf("Constructing DOME configuration topography\n");
  printf("western boundary (longitude) of full domain = %10.4f\n", xw);
  printf("eastern boundary (longitude) of full domain = %10.4f\n", xe);
  printf("southern boundary (latitude) of full domain = %10.4f\n", ys);
  printf("northern boundary (latitude) of full domain = %10.4f\n", yn);
  printf("depth (m) of the embayment                  = %10.4f\n", dome_embayment_depth);
  printf("depth (m) of the domain bottom              = %10.4f\n", dome_bottom);
  printf("western boundary (longitude) of embayment   = %10.4f\n", xw_embay);
  printf("eastern boundary (longitude) of embayment   = %10.4f\n", xe_embay);
  printf("southern boundary (latitude) of embayment   = %10.4f\n", ys_embay);
  printf("northern boundary (latitude) of embayment   = %10.4f\n", yn_embay);

  printf("slope of the dome configuration             = %10.4f\n", dome_slope);
  printf("northern boundary (latitude) of slope       = %10.4f\n", yn_slope);
  printf("southern boundary (latitude) of slope       = %10.4f\n", ys_slope);
  
  /* define flat bottom */
  for(i=0; i<nx*ny; i++) depth[i] = dome_bottom;

  /* define sloping bottom */
  for(i=0; i<nx*ny; i++) {
    if(y[i] >= ys_slope && y[i] < yn_slope ) {
      height = dome_slope*deg2metre*(y[i]-ys_slope);
      depth[i] = dome_bottom - height;
    }
  }
    
  /* define shelf */
  for(i=0; i<nx*ny; i++) {
    if( y[i] >= yn_slope ) depth[i] = 0.0;
  }

  /* define embayment */
  for(i=0; i<nx*ny; i++) {
    if( x[i] >= xw_embay && x[i] <= xe_embay &&
	y[i] >= ys_embay && y[i] <= yn_embay )
      depth[i] = dome_embayment_depth;
  }
  
}


