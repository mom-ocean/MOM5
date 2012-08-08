/*
 * Modify land grid to match ocean grid at coast and calculate atmos/land,
 * atmos/ocean, and land/ocean overlaps using the Sutherland-Hodgeman polygon
 * clipping algorithm (Sutherland, I. E. and G. W. Hodgeman, 1974:  Reentrant
 * polygon clipping, CACM, 17(1), 32-42).  Code here is non-reentrant for speed.
 *                                                     - Mike Winton (4/01)
 */

/*
 * to compile using libMPI:
 *
 *    cc -O2 -o make_xgrids_parallel make_xgrids.c -Duse_libMPI -I/usr/local/include -L/usr/local/lib -lnetcdf -lmpi -lm
 *
 * to compile without using libMPI:
 *    cc -O2 -o make_xgrids make_xgrids.c -I/usr/local/include -L/usr/local/lib -lnetcdf -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#ifdef use_libMPI
#include <mpi.h>
#endif

/* This constant is used to control the area of exchange grids. */
#define EPS1 (1e-10)
#define EPS3 (1e-4)
#define MAXLOCAL 1e7

int pe, npes, root_pe;
#ifdef use_libMPI
  MPI_Request request=MPI_REQUEST_NULL;
#endif


char *usage[] = {
"",
"             __________________ MAKE_XGRIDS __________________",
"",
"Make_xgrids generates two exchange grids for the FMS coupler, one for fluxes",
"between the atmosphere and surface (sea ice and land), the other for runoff",
"between the land and sea ice.  Make_xgrids expects a NetCDF format input",
"specification for the ice/ocean grid.  Three fields are required to be in",
"this file:",
"",
"  (1)  wet - a 2D array of double precision numbers set to 1.0 where the ice",
"       and ocean models are active and 0.0 else where.  WET has im indices in",
"       the i (pseudo east-west) direction and jm indices in the j (pseudo",
"       north-south) direction.  These correspond to the global arrays of",
"       temperature, salinity and ice thickness in the coupled climate model.",
"",
"  (2)  geolon_vert_t and (3) geolat_vert_t - 2D double precision arrays",
"       (dimensioned im+1 by jm+1) that contain the longitudes and latitudes",
"       (respectively) of the corners of the sea ice and ocean model grid cells",
"       in degrees.",
"",
"Since the atmosphere and land models are required to be latitude/longitude,",
"1D arrays of the latitude and longitude boundaries suffice to specify their",
"grids.  Make_xgrids expects this information to be provided in text files with",
"a single latitude or longitude (in degrees) per line.  In both the 2D and 1D",
"boundary specification arrays, the coordinate values must increase",
"monotonically with index.",
"",
"A sample call to make_xgrids that makes exchange grids for a tripolar grid",
"sea ice/ocean, a T42 atmosphere, and a one-degree land model might look like:",
"",
"If compiled with -Duse_libMPI, the command will be:",
"", 
"  make_xgrids -o ocean_grid_spec.nc -a atmos_grid_spec.nc -l land_grid_spec.nc",
"or",
"  make_xgrids -o tripolar_ocean.nc -a t42_x.dat,t42_y.dat -l 1deg_x.dat,1deg_y.dat",
"",
"",
"If compiled without -Duse_libMPI, the command will be:",
"", 
"  mpirun -np # make_xgrids -o ocean_grid_spec.nc -a atmos_grid_spec.nc -l land_grid_spec.nc",
"or",
"  mpirun -np # make_xgrids -o tripolar_ocean.nc -a t42_x.dat,t42_y.dat -l 1deg_x.dat,1deg_y.dat",
"",
 
"Make_xgrids copies all fields of the ocean grid specification file to its output",
"file, \"grid_spec.nc\", and appends fields that specify the atmosphere and land",
"model grids, and the surface and runoff exchange grids.  The ocean model, sea ice",
"model, land model and coupler expect to find the grid_spec.nc file in the INPUT",
"subdirectory at run time.  The coupler checks that the run time atmosphere model",
"has the same resolution as that in grid_spec.nc.",
"",
"             __________________ MAKE_XGRIDS __________________",
"",
NULL };

#define NC_CALL(F) \
  { int status; \
    if ((status=(F))!=NC_NOERR) { \
      fprintf(stderr, "netCDF error: %s\n", nc_strerror(status)); \
      exit(1); \
    } \
  }

int n_max = 0; /* maximum number of vertices in x-cell */

#define MV 20
typedef struct { double x, y; } VTX;

void error_handler(char *str)
{
  fprintf( stderr, "%s\n", str);
#ifdef use_libMPI      
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(1);
#endif  
}


/*
 * Sutherland-Hodgeman algorithm sequentially clips parts outside 4 boundaries
 */
int clip(VTX v_in[], int n_in, VTX lolf, VTX uprt, VTX v_out[])
{
  VTX v_tmp[MV], v_last, v_isect;
  int i_in, i_out, n_out, inside_last, inside;

  /* clip polygon with LEFT boundary - clip V_IN to V_TMP */
  v_last = v_in[n_in-1];
  inside_last = (v_last.x >= lolf.x);
  for (i_in=0,i_out=0;i_in<n_in;i_in++) {
 
    /* if crossing LEFT boundary - output intersection */
    if ((inside=(v_in[i_in].x >= lolf.x))!=inside_last) {
      v_isect.x = lolf.x;
      v_isect.y = v_last.y + (lolf.x - v_last.x) * (v_in[i_in].y - v_last.y)
                                                 / (v_in[i_in].x - v_last.x);
      v_tmp[i_out++] = v_isect;
    }

    /* if "to" point is right of LEFT boundary, output it */
    if (inside) v_tmp[i_out++] = v_in[i_in];

    v_last = v_in[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with RIGHT boundary - clip V_TMP to V_OUT */
  v_last = v_tmp[n_out-1];
  inside_last = (v_last.x <= uprt.x);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing RIGHT boundary - output intersection */
    if ((inside=(v_tmp[i_in].x <= uprt.x))!=inside_last) {
      v_isect.x = uprt.x;
      v_isect.y = v_last.y + (uprt.x - v_last.x) * (v_tmp[i_in].y - v_last.y)
                                                 / (v_tmp[i_in].x - v_last.x);
      v_out[i_out++] = v_isect;
    }

    /* if "to" point is left of RIGHT boundary, output it */
    if (inside) v_out[i_out++] = v_tmp[i_in];

    v_last = v_tmp[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with BOTTOM boundary - clip V_OUT to V_TMP */
  v_last = v_out[n_out-1];
  inside_last = (v_last.y >= lolf.y);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing BOTTOM boundary - output intersection */
    if ((inside=(v_out[i_in].y >= lolf.y))!=inside_last) {
      v_isect.y = lolf.y;
      v_isect.x = v_last.x + (lolf.y - v_last.y) * (v_out[i_in].x - v_last.x)
                                                 / (v_out[i_in].y - v_last.y);
      v_tmp[i_out++] = v_isect;
    }

    /* if "to" point is above BOTTOM boundary, output it */
    if (inside) v_tmp[i_out++] = v_out[i_in];

    v_last = v_out[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with TOP boundary - clip V_TMP to V_OUT */
  v_last = v_tmp[n_out-1];
  inside_last = (v_last.y <= uprt.y);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing TOP boundary - output intersection */
    if ((inside=(v_tmp[i_in].y <= uprt.y))!=inside_last) {
      v_isect.y = uprt.y;
      v_isect.x = v_last.x + (uprt.y - v_last.y) * (v_tmp[i_in].x - v_last.x)
                                                 / (v_tmp[i_in].y - v_last.y);
      v_out[i_out++] = v_isect;
    }

    /* if "to" point is below TOP boundary, output it */
    if (inside) v_out[i_out++] = v_tmp[i_in];

    v_last = v_tmp[i_in];
    inside_last = inside;
  }
  n_max = n_out>n_max ? n_out : n_max;
  return(i_out);
} /* clip */

#define TRAPEZOID_RULE 0
/*
 * poly_area - obtains area of input polygon by line integrating -sin(lat)d(lon)
 *             * Vertex coordinates must be in degrees.
 *             * Vertices must be listed counter-clockwise around polygon.
 *             * Returns area in fraction of total sphere surface area.
 */
double poly_area(VTX v[], int n)
{
  double area = 0.0;
  double D2R = M_PI/180;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (v[ip].x-v[i].x)*D2R;
    double lat1, lat2;
    lat1 = v[ip].y*D2R;
    lat2 = v[i].y*D2R;
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    
    if (lat1 == lat2) /* cheap area calculation along latitude */
      area -= dx*sin(lat1);
    else
#if TRAPEZOID_RULE
      area -= dx*(sin(lat1)+sin(lat2))/2;
#else
      area += dx*(cos(lat1)-cos(lat2))/(lat1-lat2);
#endif
  }
  return (area/(4*M_PI));
} /* poly_area */

/* This routine is used to calculate the integral which equal the product of
   latitude of the centroid and area of any grid box */

double poly_int3(VTX v[], int n)
{
  double int3 = 0.0;
  double D2R = M_PI/180;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (v[ip].x-v[i].x)*D2R;
    double lat1, lat2;
    lat1 = v[ip].y*D2R;
    lat2 = v[i].y*D2R;
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (lat1 == lat2) /* cheap area calculation along latitude */
      int3 -= dx*(cos(lat1) + lat1*sin(lat1));
    else
#if TRAPEZOID_RULE      
      int3 -= dx*(cos(lat1) + lat1*sin(lat1)+cos(lat2) + lat2*sin(lat2))/2.0;
#else
    int3 -= dx*(2*sin(lat1)-lat1*cos(lat1) - 2*sin(lat2)+lat2*cos(lat2) )/(lat1-lat2);
#endif
  }
  return (int3/(4*M_PI));
} /* poly_int3 */        

/* this routine can be used to calculate the latitude of the centroid of any grid box */
double poly_ctrlat(VTX v[], int n)
{
  double area;
  area = poly_area(v,n);
  if(area < EPS1)
     return (0.0);
  else  
     return (poly_int3(v,n)/area);
} /* poly_ctrlat */     

/* This routine is used to calculate the integral which equal the product of
   lontitude of the centroid and area of any grid box */

double poly_int2(VTX v[], int n, double clon)
{
  double int2 = 0.0;
  double D2R = M_PI/180;
  int    i;

  clon = clon * D2R;
  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
    double f1, f2, fac, fint;
    phi1   = v[ip].x*D2R;
    phi2   = v[i].x*D2R;
    lat1 = v[ip].y*D2R;
    lat2 = v[i].y*D2R;    
    dphi   = phi1 - phi2;
    
    if      (dphi==0.0) continue;

    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    /* this will make sure longitude of centroid is at
       the same interval as the center of any grid   */
    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      int2 -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      int2 -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
	}
    
  }
  return (int2/(4*M_PI));
}   /* poly_int2 */

/* this routine can be used to calculate the lontitude of the centroid of any grid box */
double poly_ctrlon(VTX v[], int n, double clon)
{

  double area;
  area = poly_area(v,n);
  if(area < EPS1)
     return (0.0);
  else
     return (poly_int2(v,n,clon)/area);
} /* poly_ctrlon */

/* end change */

double box_area(VTX ll, VTX ur)
{
  double D2R = M_PI/180;
  double dx = (ur.x-ll.x)*D2R;
  
  if(dx > M_PI)  dx = dx - 2.0*M_PI;
  if(dx < -M_PI) dx = dx + 2.0*M_PI;

  return ( dx*(sin(ur.y*D2R)-sin(ll.y*D2R))/(4*M_PI) );
} /* box_area */

/* This routine is used to calculate the integral which equal the product of
   latitude of the centroid and area of a uniform grid box */

double box_int3(VTX ll, VTX ur)
{
  double D2R = M_PI/180;
  double dphi = (ur.x-ll.x)*D2R;
 
  if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
  if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
  return ( dphi*(cos(ur.y*D2R) + ur.y*D2R*sin(ur.y*D2R)-(cos(ll.y*D2R) +
				      ll.y*D2R*sin(ll.y*D2R)))/(4*M_PI) );
} /* box_int3 */


/* this routine can be used to calculate the latitude of the centroid of a uniform grid box */
double box_ctrlat(VTX ll, VTX ur)
{
  double area;
  area = box_area(ll,ur);
  if(area < EPS1)
     return (0.0);
  else    
     return ( box_int3(ll,ur)/area);
} /* box_ctrlat */

/* This routine is used to calculate the integral which equal the product of
   lontitude of the centroid and area of a uniform grid box */

double box_int2(VTX ll, VTX ur, double clon)
{
  double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
  double f1, f2, fac, fint;  
  double D2R   = M_PI/180;
  double int2  = 0.0;
  int i;
  clon = clon * D2R;  
  for( i =0; i<2; i++) {
    if(i == 0) {
      phi1 = ur.x * D2R;
      phi2 = ll.x * D2R;
      lat1 = lat2 = ll.y * D2R;
    }
    else {
      phi1 = ll.x * D2R;
      phi2 = ur.x * D2R;
      lat1 = lat2 = ur.y * D2R;
    }
    dphi   = phi1 - phi2;
    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    /* make sure the center is in the same grid box. */
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      int2 -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      int2 -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
	}
  }
  return (int2/(4*M_PI));    
	   } /* box_int2 */

/* this routine can be used to calculate the lontitude of the centroid of a uniform grid box */
double box_ctrlon(VTX ll, VTX ur, double clon)
{
  double area;
  area = box_area(ll,ur);
  if(area < EPS1)
     return (0.0);
  else      
     return(box_int2(ll,ur,clon)/area);
}   /* box_ctrlon */



int vtx_delete(VTX v[], int n, int n_del)
{
  for (;n_del<n-1;n_del++) v[n_del] = v[n_del+1];
  return (n-1);
} /* vtx_delete */

int vtx_insert(VTX v[], int n, int n_ins, double x, double y)
{
  int i;

  for (i=n-1;i>=n_ins;i--) v[i+1] = v[i];
  v[n_ins].x = x;
  v[n_ins].y = y;
  return (n+1);
} /* vtx_insert */

void v_print(VTX v[], int n)
{
  int i;

  for (i=0;i<n;i++) printf(" %20g   %20g\n", v[i].x, v[i].y);
} /* v_print */

int lon_fix(VTX v[], int n, double tlon)
{
  double x_sum, dx;
  int i, nn = n, pole = 0;

  for (i=0;i<nn;i++) if (fabs(v[i].y)>=90.0-EPS3) pole = 1;
  if (0&&pole) {
    printf("fixing pole cell\n");
    v_print(v,nn);
    printf("---------");
  }

  /* all pole points must be paired */
  for (i=0;i<nn;i++) if (fabs(v[i].y)>=90.0-EPS3) {
    int im=(i+nn-1)%nn, ip=(i+1)%nn;

    if (v[im].y==v[i].y && v[ip].y==v[i].y) {
      nn = vtx_delete(v, nn, i);
      i--;
    } else if (v[im].y!=v[i].y && v[ip].y!=v[i].y) {
      nn = vtx_insert(v, nn, i, v[i].x, v[i].y);
      i++;
    }
  }
  /* first of pole pair has longitude of previous vertex */
  /* second of pole pair has longitude of subsequent vertex */
  for (i=0;i<nn;i++) if (fabs(v[i].y)>=90.0-EPS3) {
    int im=(i+nn-1)%nn, ip=(i+1)%nn;

    if (v[im].y!=v[i].y) v[i].x = v[im].x;
    if (v[ip].y!=v[i].y) v[i].x = v[ip].x;
  }

  if (nn) x_sum = v[0].x; else return(0);
  for (i=1;i<nn;i++) {
    double dx = v[i].x-v[i-1].x;

    if      (dx < -180) dx = dx + 360;
    else if (dx >  180) dx = dx - 360;
    x_sum += (v[i].x = v[i-1].x + dx);
  }

  dx = (x_sum/nn)-tlon;
  if      (dx < -180) for (i=0;i<nn;i++) v[i].x += 360;
  else if (dx >  180) for (i=0;i<nn;i++) v[i].x -= 360;

  if (0&&pole) {
    printf("area=%g\n", poly_area(v,nn));
    v_print(v,nn);
    printf("---------");
  }

  return (nn);
} /* lon_fix */

/* This routine will routine number of files that holds the grid. This routine assume
   the last field in the grid file is lastfield and file format will be grid.nc, grid1.nc,
   grid2.nc, the return value will be 3.
*/
int get_number_files(char *filename, char *lastfield)
{
  int n, status, varid, ncid;
  char curfile[256], str[24];
  size_t pathlength;
  
  n = 0;
  for(;;) { /* infinite loop */
    if(n == 0) {
      strcpy(curfile, filename);
    }
    else {
      pathlength = strlen(filename) - 3;
      if(strncmp( filename+pathlength, ".nc", 3) != 0)
	error_handler("The suffix of the input file is not .nc");
      strncpy(curfile, filename, pathlength);
      curfile[pathlength] = '\0';
      sprintf(str, "%d", n);
      strcat(curfile, str);
      strcat(curfile, ".nc");
    }
    NC_CALL(nc_open(curfile, NC_NOWRITE, &ncid));
    status = nc_inq_varid(ncid, lastfield, &varid);
    NC_CALL(nc_close(ncid));
    if (status == NC_NOERR) return n+1;
    n++;
  }
} /*get_number_files*/ 


/* get ncid and varid associated with fname and name and return 1 if found the field,
 otherwise return 0
*/
int get_ncid_varid( filename, fieldname, nfiles, ncid, varid )
char *filename;
char *fieldname;
int nfiles;
int *ncid;
int *varid;
{
  int status, n;
  char curfile[256]="", str[24]="";
  size_t pathlength;
  
  n = 0;
  *ncid = 0;
  *varid = 0;
  for(;;) { /* infinite loop */
    if( n >= nfiles ) return 0;
    if(n == 0) {
      strcpy(curfile, filename);
    }
    else {
      pathlength = strlen(filename) - 3;
      if(strncmp( filename+pathlength, ".nc", 3) != 0)
	error_handler("The suffix of the input file is not .nc");
      strncpy(curfile, filename, pathlength);
      curfile[pathlength] = '\0';
      sprintf(str, "%d", n);
      strcat(curfile, str);
      strcat(curfile, ".nc");      
    }
    NC_CALL(nc_open(curfile, NC_NOWRITE, ncid));
    status = nc_inq_varid(*ncid, fieldname, varid);	 
    if (status == NC_NOERR) return 1;
    NC_CALL(nc_close(*ncid));
    n++;
  }

  return 0;

} /* get_ncid_varid */

  
void *get_var2d_field3d(fname, name, nfiles, size1ptr, size2ptr)
char *fname;
char *name;
int nfiles;
int  *size1ptr, *size2ptr;
{
   double *var;
   double *result;
   int  ncid, ncvid, dims[4];
   int status, i,j,m,n;
   size_t im, jm, km;

   status = get_ncid_varid( fname, name, nfiles, &ncid, &ncvid );
   if(status == 0) return 0;
   NC_CALL(nc_inq_vardimid(ncid, ncvid, dims))
   NC_CALL(nc_inq_dimlen(ncid, dims[0], &km)) 
   NC_CALL(nc_inq_dimlen(ncid, dims[1], &jm))
   NC_CALL(nc_inq_dimlen(ncid, dims[2], &im))

   if( km != 4) {
      fprintf(stderr,"each T-cell should have 4 vertices %s\n", "get_nc_var2d");
      exit(1);
   }          
   
   *size1ptr = im;
   *size2ptr = jm;
    var = (double *) malloc ( im*jm*4*sizeof(double) );
   result = (double *) malloc ( (im+1)*(jm+1)*sizeof(double) );
   if (nc_get_var_double(ncid, ncvid, var)==-1) {
      fprintf(stderr,"Couldn't get variable %s\n", name);
      exit(1);
   }
   NC_CALL(nc_close(ncid))
     m=0;
     for(j=0;j<jm;j++) {
       for(i=0;i<im;i++) {	 
	 n = j*im + i;
         result[m] = var[n];
	 m = m + 1;
       }
       n = im*jm + (j+1)*im - 1;
       result[m] = var[n];
       m = m + 1;
     }
   for(i=0;i<im;i++) {
     n = 4*im*jm - im + i;
     result[m] = var[n];
     m = m + 1;
   }
   result[m] = var[3*im*jm-1];
   if(m != (im+1)*(jm+1)-1) {
      fprintf(stderr,"size is not comformable %s\n","get_nc_var2d");
      exit(1);
   }
   return(result);
} /* get_var2d_field3d */

void *get_var2d_field2d(fname, name, nfiles, size1ptr, size2ptr)
char *fname;
char *name;
int nfiles;
int  *size1ptr, *size2ptr;
{
   void *var;
   int  ncid, ncvid, dims[4];
   int status;
   size_t im, jm;

   status = get_ncid_varid( fname, name, nfiles, &ncid, &ncvid );
   if(status == 0) return 0;   
   NC_CALL(nc_inq_vardimid(ncid, ncvid, dims))
   NC_CALL(nc_inq_dimlen(ncid, dims[0], &jm))
   NC_CALL(nc_inq_dimlen(ncid, dims[1], &im))

   *size1ptr = im;
   *size2ptr = jm;

   var = (void *) malloc ( (*size1ptr)*(*size2ptr)*sizeof(double) );

   if (nc_get_var_double(ncid, ncvid, var)==-1) {
      fprintf(stderr,"Couldn't get variable %s\n", name);
      exit(1);
   }
   NC_CALL(nc_close(ncid))
   return(var);
} /* get_var2d_field2d */


void *get_var1d_field2d(fname, name, nfiles, sizeptr, cart)
char *fname;
char *name;
int nfiles;
int  *sizeptr;
char cart;
{
   double *var, *result;
   int  ncid, ncvid, dims[4];
   int status, i, j;
   size_t im, jm;

   status = get_ncid_varid( fname, name, nfiles, &ncid, &ncvid );
   if(status == 0) return 0;   

   NC_CALL(nc_inq_vardimid(ncid, ncvid, dims))
   NC_CALL(nc_inq_dimlen(ncid, dims[0], &jm))
   NC_CALL(nc_inq_dimlen(ncid, dims[1], &im))
   if (cart == 'x')
      *sizeptr = im;
   else if(cart == 'y')
      *sizeptr = jm;
   else {
     fprintf(stderr,"cart should be either x or y %s\n", "get_var1d_field2d");
     exit(1);
   }
   result = (double *) malloc ( (*sizeptr)*sizeof(double) );
   var    = (double *) malloc ( im*jm*sizeof(double) );   
   if (nc_get_var_double(ncid, ncvid, var)==-1) {
      fprintf(stderr,"Couldn't get variable %s\n", name);
      exit(1);
   }
   NC_CALL(nc_close(ncid))
   if (cart == 'x') {
     for(i=0;i<im;i++)
       result[i] = var[i];
   }
   else {
     for(j=0;j<jm;j++)
       result[j] = var[j*im];
   }
   
   return(result);
} /* get_var1d_field2d */

void *get_var1d_field3d(fname, name, nfiles, sizeptr, cart)
char *fname;
char *name;
int nfiles;
int  *sizeptr;
char cart;
{
   double *var;
   double *result;
   int  ncid, ncvid, dims[4];
   int status, i,j,k;
   size_t im, jm, km;

   status = get_ncid_varid( fname, name, nfiles, &ncid, &ncvid );
   if(status == 0) return 0;   
   
   NC_CALL(nc_inq_vardimid(ncid, ncvid, dims))
   NC_CALL(nc_inq_dimlen(ncid, dims[0], &km)) 
   NC_CALL(nc_inq_dimlen(ncid, dims[1], &jm))
   NC_CALL(nc_inq_dimlen(ncid, dims[2], &im))

   if( km != 4) {
      fprintf(stderr,"each T-cell should have 4 vertices %s\n", "get_nc_var2d");
      exit(1);
   }          
   if (cart == 'x')
      *sizeptr = im;
   else if(cart == 'y')
      *sizeptr = jm;
   else {
     fprintf(stderr,"cart should be either x or y%s\n", "get_var1d_field2d");
     exit(1);
   }   

   var = (double *) malloc ( im*jm*4*sizeof(double) );
   result = (double *) malloc ( (*sizeptr+1)*sizeof(double) );
   if (nc_get_var_double(ncid, ncvid, var)==-1) {
      fprintf(stderr,"Couldn't get variable %s\n", name);
      exit(1);
   }
   NC_CALL(nc_close(ncid))

   /* check if the grid is rectangular grid or not */
   if (cart == 'x') {
     for(k=0;k<4;k++)
       for(i=0;i<im;i++) 
	 for(j=1;j<jm;j++)
           if(var[im*jm*k+j*im+i] != var[im*jm*k+i]){
	     fprintf(stderr,"the grid is not rectangular grid.%s\n", name);
	     exit(1);
	   }
   }
   else if(cart == 'y'){
     for(k=0;k<4;k++)
       for(j=0;j<jm;j++)
         for(i=1;i<im;i++) 
           if(var[im*jm*k+j*im+i] != var[im*jm*k+j*im]){
	     fprintf(stderr,"the grid is not rectangular grid.%s\n", name);
	     exit(1);
	   }
   }
   if (cart == 'x') {
     for(i=0;i<im;i++)
       result[i] = var[i];
     result[im] = var[im*jm*2-1];
   }
   else {
     for(j=0;j<jm;j++)
       result[j] = var[j*im];
     result[jm] = var[im*jm*4-1];
   }
     
   return(result);
} /* get_var1d_field3d */

void nc_copy_1d_int(int file1, int var1, int file2, int var2, int n)
{
  size_t i;
  int    x;

  for (i=0;i<n;i++) {
    NC_CALL(nc_get_var1_int(file1, var1, &i, &x))
    NC_CALL(nc_put_var1_int(file2, var2, &i, &x))
  }
} /* nc_copy_1d_int */

void nc_copy_1d_double(int file1, int var1, int file2, int var2, int n)
{
  size_t i;
  double x;

  for (i=0;i<n;i++) {
    NC_CALL(nc_get_var1_double(file1, var1, &i, &x))
    NC_CALL(nc_put_var1_double(file2, var2, &i, &x))
  }
} /* nc_copy_1d_double */

/*
 * read vertices from sdin; clip with given window; output vertices to stdout
 */
double interactive_test(double ll_x, double ll_y, double ur_x, double ur_y)
{
  VTX clip_ll, clip_ur, v_in[MV], v_out[MV];
  double x, y;
  int n_in = 0, n_out;

  clip_ll.x = ll_x; clip_ll.y = ll_y;
  clip_ur.x = ur_x; clip_ur.y = ur_y;
  while (n_in<(MV-1) && scanf("%lf %lf",&x,&y)!=EOF) {
    v_in[n_in].x = x;
    v_in[n_in].y = y;
    n_in++;
  }

  printf("\nINTERACTIVE TEST FIXED INPUT:\n\n");
  n_in = lon_fix(v_in, n_in, (ll_x+ur_x)/2);
  v_print(v_in, n_in);
  printf("\nINTERACTIVE TEST INTERSECTION:\n\n");
  n_out = clip ( v_in, n_in, clip_ll, clip_ur, v_out );
  v_print(v_out, n_out);
  return(poly_area(v_out, n_out));
} /* interactive_test */

/*   send local_count to root_pe and sum up  */
int get_global_count( int local_count ) 
{
  int global_count;
  global_count = local_count;

#ifdef use_libMPI
  {
    int p, recv_count;
    MPI_Status status;
    
    if( pe != root_pe ) { /* send to root pe */
      MPI_Isend( &local_count, 1,  MPI_INT,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    }
    else {  /* receive data */
      for(p = root_pe+1; p<npes; p++) {
	MPI_Recv( &recv_count, 1,  MPI_INT,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	global_count += recv_count;
      }
    }
    if( request != MPI_REQUEST_NULL ) MPI_Wait( &request, &status);
  }
#endif    

  return global_count;
  
} /* get_global_count */ 


void get_global_data_double( double *ldata, double *gdata, int count )
{
  int i;
  
  if( pe == root_pe ) {
    for( i = 0; i< count; i++ ) gdata[i] = ldata[i];
  }

#ifdef use_libMPI
  {
    int p, recv_count;
    MPI_Status status;
    if( pe != root_pe ) { /* send to root pe */
      MPI_Isend( &count, 1,  MPI_INT,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);    
      MPI_Isend( ldata, count,  MPI_DOUBLE,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    }
    else {  /* receive data */
      gdata += count;
      for(p = root_pe+1; p<npes; p++) {
	MPI_Recv( &recv_count, 1,  MPI_INT,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv( gdata, recv_count,  MPI_DOUBLE,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	gdata += recv_count;
      }
    }
    if( request != MPI_REQUEST_NULL ) MPI_Wait( &request, &status);   
  }
#endif    
  
} /*get_global_data_double*/

void get_global_data_int( int *ldata, int *gdata, int count )
{
  int i;
  
  if( pe == root_pe ) {
    for( i = 0; i< count; i++ ) gdata[i] = ldata[i];
  }

#ifdef use_libMPI
  {
    int p, recv_count;
    MPI_Status status;    
    if( pe != root_pe ) { /* send to root pe */
      MPI_Isend( &count, 1,  MPI_INT,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);    
      MPI_Isend( ldata, count,  MPI_INT,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    }
    else {  /* receive data */
      MPI_Status status;
      gdata += count;
      for(p = root_pe+1; p<npes; p++) {
	MPI_Recv( &recv_count, 1,  MPI_INT,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv( gdata, recv_count,  MPI_INT,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	gdata += recv_count;
      }
    }
    if( request != MPI_REQUEST_NULL ) MPI_Wait( &request, &status);    
  }
#endif    
  
} /*get_global_data_int*/

void sum_over_pe( double *data, int count )
{
  
#ifdef use_libMPI
  {
    MPI_Status status;    
    if( pe != root_pe ) { /* send to root pe */
      MPI_Isend( data, count,  MPI_DOUBLE,  root_pe, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    }
    else {  /* receive data */
      int p, i;
      double *ldata;
      ldata = (double *) malloc(count*sizeof(double));
      for(p = root_pe+1; p<npes; p++) {
	MPI_Recv( ldata, count,  MPI_DOUBLE,  p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	for(i=0; i<count; i++) data[i] += ldata[i];
      }
      free(ldata); 
    }
    if( request != MPI_REQUEST_NULL ) MPI_Wait( &request, &status);    
  }
#endif    
  
} /*sum_over_pe*/


/*
 * Example atmosphere and land grids (both regular and lat/lon)
 */
#define MN  10000

int IMA = 144;
int JMA =  90;
int IML = 360;
int JML = 180;

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))

int nums_from_file(fname, nums)
char   *fname;
double *nums[];
{
  FILE   *fp = fopen(fname, "r");
  double d[MN];
  int    i, n;

  for (n=0;n<MN;n++) if (fscanf(fp, " %lg", d+n)<=0) break;
  *nums = (double *) malloc(n*sizeof(double));
  for (i=0;i<n;i++) (*nums)[i] = d[i];
  return (n);
} /* nums_from_file */

void xynums(str, xptr, nxptr, yptr, nyptr)
char *str;
double *xptr[], *yptr[];
int    *nxptr, *nyptr;
{
  char xname[BUFSIZ], yname[BUFSIZ];
  int i;

  if (sscanf(str, "%[^,],%[^,]", xname, yname)!=2) {
    fprintf(stderr, "bad option format: %s\n", str);
    exit(1);
  }
  *nxptr = nums_from_file(xname, xptr);
  for (i=1;i<*nxptr;i++) if ((*xptr)[i]<=(*xptr)[i-1]) {
    fprintf(stderr, "%s numbers not increasing ( ... %g, %g ...)\n",
	       xname, (*xptr)[i-1], (*xptr)[i] );
    exit (1);
  }
  *nyptr = nums_from_file(yname, yptr);
  for (i=1;i<*nyptr;i++) if ((*yptr)[i]<=(*yptr)[i-1]) {
    fprintf(stderr, "%s numbers not increasing ( ... %g, %g ...)\n",
	       yname, (*yptr)[i-1], (*yptr)[i] );
    exit (1);
  }
} /* xynums */

/* the following will define the domain decomposition */
void define_domain(int npts, int ndivs, int *startlist, int *endlist )
{
  int n, npts_left, start, npts_used;
  
  if(ndivs > npts )
    error_handler("define_domain: number of divisions is larger than number of points");

  npts_left = npts;
  start = 0;
  for(n=0; n<ndivs; n++){
    npts_used = ceil(npts_left*1.0/(ndivs-n));
    startlist[n] = start;
    endlist[n] = start + npts_used - 1;
    start = start + npts_used;
    npts_left=npts_left - npts_used;
  }

} /* define_domain */

/* global variable */

main (int argc, char *argv[])
{
  int c;
  extern char *optarg;
  extern int optind;
  char *ofile = NULL;
  char *afile = NULL;
  char *lfile = NULL;
  char *gfile = "grid_spec.nc";
  char *axo_file = "atmXocn.nc";
  char *axl_file = "atmXlnd.nc";
  char *lxo_file = "lndXocn.nc";
  int errflg = (argc == 1);
  int im, jm, n;
  int    IMO=0, JMO=0;
  int    ia, ja, io, jo, il, jl, i, j, na;
  double *wet, *ocnlonb, *ocnlatb, *atmlonb, *atmlatb, *lndlonb, *lndlatb;
  double *atmlont, *atmlatt, *lndlont, *lndlatt ;  
  double *ocn_area, *ocn_area_left, *ocn_minus_runoff, *atm_area, *atm_area_left, *lnd_area, *gatm_area;
  double *lnd_cell_area, *lnd_cell_ctrlon, *lnd_cell_ctrlat;
  double *atm_ctrlon,*atm_ctrlat;
  int ia_max_err=0, ja_max_err=0, io_max_err=0, jo_max_err=0;

  int ncidAO, dim_AO, area_AO_id, ia_AO_id, ja_AO_id, io_AO_id, jo_AO_id;
  int ncidAL, dim_AL, area_AL_id, ia_AL_id, ja_AL_id, il_AL_id, jl_AL_id;
  int dj_AL_id, di_AL_id, dj_AO_id, di_AO_id, dj_LO_id, di_LO_id;
  int ncidLO, dim_LO, area_LO_id, il_LO_id, jl_LO_id, io_LO_id, jo_LO_id;
  size_t iAXO=0, iAXL=0, iLXO=0;
  double tot_area_AXO = 0.0, tot_area_AXL = 0.0;
  int nofile = 0;
  int land_same_as_atmos, isc, iec, natm_local;
  int *startlist, *endlist;

  /* First parallize the code */
#ifdef use_libMPI
  {
    int zero = 0;
    MPI_Init(&zero,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&pe);
    MPI_Comm_size(MPI_COMM_WORLD,&npes);
  }
#else
  pe = 0;
  npes = 1;
#endif  
  root_pe = 0;
  
  /*
   * process command line
   */
  while ((c = getopt(argc, argv, "i:a:l:o:g:")) != -1)
    switch (c) {
    case 'i':
      {
        double ll_x, ll_y, ur_x, ur_y;

        sscanf(optarg, "%lf%*1[,]%lf%*1[,]%lf%*1[,]%lf",
                                                    &ll_x, &ll_y, &ur_x, &ur_y);
        fprintf(stderr, "\nINTERSECTION AREA = %g\n", 
                                      interactive_test(ll_x, ll_y, ur_x, ur_y));
      }
      exit(1);
      break;
    case 'a': 
      afile = optarg;
      break;
    case 'l': 
      lfile = optarg;
      break;
    case 'o':
      ofile = optarg;
      break;
    case 'g':
      gfile = optarg;
      break;
    case '?':
      errflg++;
    }
  if (errflg || !afile ) {
    char **u = usage;

    while (*u) { if(pe == root_pe) fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }

  /*
   * Read ocean grid boundaries and mask (where water is)
   */

  if (ofile) {
    /* The ocean grid may be distributed in multiple files, get number of ocean files */
    nofile = get_number_files(ofile, "wet"); 
    if (!(ocnlonb = get_var2d_field3d(ofile, "x_vert_T", nofile, &io, &jo))) {
      if (!(ocnlonb = get_var2d_field2d(ofile, "geolon_vert_t", nofile, &io, &jo))) {
	fprintf(stderr, "Couldn't get x_vert_T/geolon_vert_t from %s\n", ofile);
	exit(1);
      }
      else {
	io = io-1;
	jo = jo-1;
	if(pe == root_pe)printf("\nGot geolon_vert_t from %s\n",ofile);
      }
    }
    else { if(pe == root_pe) printf("\nGot x_vert_T from %s\n",ofile); }
    
    if (!(ocnlatb = get_var2d_field3d(ofile, "y_vert_T", nofile, &IMO, &JMO))) {
      if (!(ocnlatb = get_var2d_field2d(ofile, "geolat_vert_t", nofile, &IMO, &JMO))) {
	fprintf(stderr, "Couldn't get y_vert_T/geolat_vert_t from %s\n", ofile);
	exit(1);
      }
      else {
	IMO = IMO-1;
	JMO = JMO-1;
	if(pe == root_pe)printf("Got geolat_vert_t from %s\n",ofile);
      }
    } 
    else {if(pe == root_pe)printf("Got y_vert_T from %s\n",ofile);}
    if (IMO != io || JMO != jo) {
      fprintf(stderr, "x_vert_T/geolon_vert_t and y_vert_T/geolat_vert_t must be same size\n");
         exit(1);
      } 

    if (!(wet     = get_var2d_field2d(ofile, "wet", nofile, &IMO, &JMO))) {
      fprintf(stderr, "Couldn't get wet from %s\n"    , ofile);
      exit(1);
    } 
    else {if(pe == root_pe) printf("Got wet from %s\n", ofile);}

    if (IMO != io || JMO != jo) {
      fprintf(stderr, "wet and x_vert_T and y_vert_T must be same size \n");
      exit(1);
    }
  } else {if(pe == root_pe)printf("\nNo \"-o\" - No Ocean Grid\n"); }
  
    /*
   * Read atmosphere grid
   */
  n = strlen(afile);
  if(strcmp(afile+n-3,".nc") == 0) {
    if (!(atmlonb = get_var1d_field3d(afile, "x_vert_T", 1, &im, 'x'))) {
      fprintf(stderr, "Couldn't get x_vert_T from %s\n", afile);
      exit(1);
    } 
    else {if(pe == root_pe)printf("Got x_vert_T from %s\n", afile);}
    if (!(atmlatb = get_var1d_field3d(afile, "y_vert_T", 1, &jm, 'y'))) {
      fprintf(stderr, "Couldn't get y_vert_T from %s\n", afile);
      exit(1);
    } 
    else {if(pe == root_pe)printf("Got y_vert_T from %s\n", afile);}

    IMA = im;
    JMA = jm;

    if (!(atmlont = get_var1d_field2d(afile, "x_T", 1, &im, 'x'))) {
      fprintf(stderr, "Couldn't get x_T from %s\n", afile);
      exit(1);
    } 
    else {if(pe == root_pe) printf("Got x_T from %s\n", afile);}
    if (!(atmlatt = get_var1d_field2d(afile, "y_T", 1, &jm,'y'))) {
      fprintf(stderr, "Couldn't get y_T from %s\n", afile);
      exit(1);
    } 
    else {if(pe == root_pe)printf("Got y_T from %s\n", afile);}
    
  }
  else {
      xynums(afile, &atmlonb, &im, &atmlatb, &jm);
      IMA = im-1;
      JMA = jm-1;     
      atmlont = (double *) malloc (IMA*sizeof(double));
      atmlatt = (double *) malloc (JMA*sizeof(double));
      for(i=0;i<IMA;i++) atmlont[i] = 0.5*(atmlonb[i]+atmlonb[i+1]);
      for(i=0;i<JMA;i++) atmlatt[i] = 0.5*(atmlatb[i]+atmlatb[i+1]);
  }

  if(pe == root_pe ) {
    printf("\nAtmosphere Longitude Boundaries: ");
    for (i=0;i<=IMA;i++) printf(" %g", atmlonb[i]);
    printf("\n");
    printf("\nAtmosphere Latitude Boundaries: ");
    for (i=0;i<=JMA;i++) printf(" %g", atmlatb[i]);
    printf("\n");
  }
    /*
   * Read land grid
   */
  if (lfile) {
    n = strlen(lfile);
    if(strcmp(lfile+n-3,".nc") == 0) {    
      if (!(lndlonb = get_var1d_field3d(lfile, "x_vert_T", 1, &im, 'x'))) {
	fprintf(stderr, "Couldn't get x_vert_T from %s\n", lfile);
	exit(1);
      } 
      else {if(pe == root_pe)printf( "Got x_vert_T from %s\n", lfile);}
      if (!(lndlatb = get_var1d_field3d(lfile, "y_vert_T",  1, &jm, 'y'))) {
	fprintf(stderr, "Couldn't get y_vert_T from %s\n", lfile);
	exit(1);
      } 
      else {if(pe == root_pe)printf( "Got y_vert_T from %s\n", lfile);}
    
      IML = im;
      JML = jm;
      
      if (!(lndlont = get_var1d_field2d(lfile, "x_T", 1, &im,'x'))) {
	fprintf(stderr, "Couldn't get x_T from %s\n", lfile);
	exit(1);
      } 
      else {if(pe == root_pe)printf("Got x_T from %s\n", lfile);}
      if (!(lndlatt = get_var1d_field2d(lfile, "y_T",  1, &jm,'y'))) {
	fprintf(stderr, "Couldn't get y_T from %s\n", lfile);
	exit(1);
      }
      else {if(pe == root_pe)printf("Got y_T from %s\n", lfile);}

    }
    else {
      xynums(lfile, &lndlonb, &im, &lndlatb, &jm);
      IML = im-1;
      JML = jm-1;      
      lndlont = (double *) malloc (IML*sizeof(double));
      lndlatt = (double *) malloc (JML*sizeof(double));
      for(i=0;i<IML;i++) lndlont[i] = 0.5*(lndlonb[i]+lndlonb[i+1]);
      for(i=0;i<JML;i++) lndlatt[i] = 0.5*(lndlatb[i]+lndlatb[i+1]);
    }
  }
  else {   /* if land boundaries not specified, use atmosphere's */
    lndlonb = atmlonb;
    lndlatb = atmlatb;
    lndlont = atmlont;
    lndlatt = atmlatt;
    IML = IMA;
    JML = JMA;
    if(pe == root_pe) printf("\nNo \"-l\" - Land Grid Same As Atmosphere Grid\n");
  }
  land_same_as_atmos = 0;
  if( IML == IMA && JML == JMA ) {
    land_same_as_atmos = 1;
    for(i = 0; i< im; i++ ) {
      if( lndlonb[i] != atmlonb[i] ) {
	land_same_as_atmos = 0;
	break;
      }
    }
    if( land_same_as_atmos == 1 ) {
      for(j = 0; j< jm; j++ ) {
	if( lndlatb[j] != atmlatb[j] ) {
	  land_same_as_atmos = 0;
	  break;
	}
      }
    }
  }

  if( pe == root_pe) {
    for (i=0;i<=IML;i++) printf(" %g", lndlonb[i]);
    printf("\n");
    printf("\nLand Latitude Boundaries: ");
    for (i=0;i<=JML;i++) printf(" %g", lndlatb[i]);
    printf("\n");
  }

  startlist = (int *) malloc( npes*sizeof(int));
  endlist   = (int *) malloc( npes*sizeof(int));  
  define_domain( IMA*JMA, npes, startlist, endlist);
  isc = startlist[pe];
  iec = endlist[pe];  
  natm_local = iec - isc + 1;

  if(natm_local > MAXLOCAL) error_handler(" number of atmos points on this pe is greater MAXLOCAL, increase MAXLOCAL");
  
  atm_area      = (double *) malloc (MAXLOCAL*sizeof(double));
  atm_ctrlon    = (double *) malloc (MAXLOCAL*sizeof(double));
  atm_ctrlat    = (double *) malloc (MAXLOCAL*sizeof(double));  
  atm_area_left = (double *) malloc (MAXLOCAL*sizeof(double));

  for (na = 0; na< natm_local; na ++) {
    VTX ll, ur;
    ia = (na+isc)%IMA;
    ja = (na+isc)/IMA;
    
    ll.x = atmlonb[ia  ]; ll.y = atmlatb[ja  ];
    ur.x = atmlonb[ia+1]; ur.y = atmlatb[ja+1];
    atm_area_left[na]   = atm_area[na] = box_area(ll, ur);
    atm_ctrlon[na]      = box_ctrlon(ll, ur, 0.5*(ll.x+ur.x));
    atm_ctrlat[na]      = box_ctrlat(ll, ur);
  }

  /* get global atm_area */
  gatm_area = (double *)malloc(IMA*JMA*sizeof(double));
  get_global_data_double(atm_area, gatm_area, natm_local);
  
  lnd_area   = (double *) malloc(IML*JML*sizeof(double));

  for (i=0;i<IML*JML;i++) lnd_area[i] = 0.0;

  
  lnd_cell_area   = (double *) malloc(IML*JML*sizeof(double));
  lnd_cell_ctrlon = (double *) malloc(IML*JML*sizeof(double));
  lnd_cell_ctrlat = (double *) malloc(IML*JML*sizeof(double));
  for (jl=0;jl<JML;jl++) for (il=0;il<IML;il++) {
    VTX ll, ur;

    ll.x = lndlonb[il  ]; ll.y = lndlatb[jl  ];
    ur.x = lndlonb[il+1]; ur.y = lndlatb[jl+1];
    lnd_cell_area[jl*IML+il] = box_area(ll, ur);
    lnd_cell_ctrlon[jl*IML+il] = box_ctrlon(ll, ur, 0.5*(ll.x+ur.x));
    lnd_cell_ctrlat[jl*IML+il] = box_ctrlat(ll, ur);
  }

  ocn_area      = (double *) malloc(IMO*JMO*sizeof(double));
  ocn_area_left = (double *) malloc(IMO*JMO*sizeof(double));
  ocn_minus_runoff = (double *) malloc(IMO*JMO*sizeof(double));
  /*find the centroid of each ocean grid */
  
  for (jo=0;jo<JMO;jo++) for (io=0;io<IMO;io++) if (wet[jo*IMO+io]>0.5) {
    VTX v[MV];
    int n_in;

    /* list ocean cell vertices in counter-clockwise order */
    v[0].y = ocnlatb[ jo   *(IMO+1)+io  ]; /* SW - y */
    v[1].y = ocnlatb[ jo   *(IMO+1)+io+1]; /* SE - y */
    v[2].y = ocnlatb[(jo+1)*(IMO+1)+io+1]; /* NE - y */
    v[3].y = ocnlatb[(jo+1)*(IMO+1)+io  ]; /* NW - y */

    v[0].x = ocnlonb[ jo   *(IMO+1)+io  ]; /* SW - x */
    v[1].x = ocnlonb[ jo   *(IMO+1)+io+1]; /* SE - x */
    v[2].x = ocnlonb[(jo+1)*(IMO+1)+io+1]; /* NE - x */
    v[3].x = ocnlonb[(jo+1)*(IMO+1)+io  ]; /* NW - x */
    n_in = lon_fix(v, 4, 180.0);
    if ((ocn_minus_runoff[jo*IMO+io] = ocn_area_left[jo*IMO+io] =
				    ocn_area[jo*IMO+io] = poly_area(v, n_in))<0) {
      fprintf(stderr, "\nOcean cell at i=%d j=%d has negative area!\n", io, jo);
      v_print(v, 4);
      fprintf(stderr,"\n");
      exit(1);
    }

  } else ocn_minus_runoff[jo*IMO+io] = ocn_area_left[jo*IMO+io] = ocn_area[jo*IMO+io] = 0.0;

  /*
   * Define output file variables for atmos/ocean and atmos/land X-grids
   */
  NC_CALL(nc_create(axo_file, NC_NOCLOBBER, &ncidAO))
  NC_CALL(nc_def_dim(ncidAO, "i_atmXocn", NC_UNLIMITED, &dim_AO))
  NC_CALL(nc_def_var(ncidAO, "AREA_ATMxOCN", NC_DOUBLE,1, &dim_AO, &area_AO_id))
  NC_CALL(nc_def_var(ncidAO, "DI_ATMxOCN", NC_DOUBLE,1, &dim_AO, &di_AO_id))    
  NC_CALL(nc_def_var(ncidAO, "DJ_ATMxOCN", NC_DOUBLE,1, &dim_AO, &dj_AO_id))
  NC_CALL(nc_def_var(ncidAO, "I_ATM_ATMxOCN",  NC_INT,  1, &dim_AO, &ia_AO_id))
  NC_CALL(nc_def_var(ncidAO, "J_ATM_ATMxOCN",  NC_INT,  1, &dim_AO, &ja_AO_id))
  NC_CALL(nc_def_var(ncidAO, "I_OCN_ATMxOCN",  NC_INT,  1, &dim_AO, &io_AO_id))
  NC_CALL(nc_def_var(ncidAO, "J_OCN_ATMxOCN",  NC_INT,  1, &dim_AO, &jo_AO_id))
  NC_CALL(nc_enddef(ncidAO))

  NC_CALL(nc_create(axl_file, NC_NOCLOBBER, &ncidAL))
  NC_CALL(nc_def_dim(ncidAL, "I_ATMxLND", NC_UNLIMITED, &dim_AL))
  NC_CALL(nc_def_var(ncidAL, "AREA_ATMxLND", NC_DOUBLE,1, &dim_AL, &area_AL_id))
  NC_CALL(nc_def_var(ncidAL, "DI_ATMxLND", NC_DOUBLE,1, &dim_AL, &di_AL_id))    
  NC_CALL(nc_def_var(ncidAL, "DJ_ATMxLND", NC_DOUBLE,1, &dim_AL, &dj_AL_id))
  NC_CALL(nc_def_var(ncidAL, "I_ATM_ATMxLND",  NC_INT,  1, &dim_AL, &ia_AL_id))
  NC_CALL(nc_def_var(ncidAL, "J_ATM_ATMxLND",  NC_INT,  1, &dim_AL, &ja_AL_id))
  NC_CALL(nc_def_var(ncidAL, "I_LND_ATMxLND",  NC_INT,  1, &dim_AL, &il_AL_id))
  NC_CALL(nc_def_var(ncidAL, "J_LND_ATMxLND",  NC_INT,  1, &dim_AL, &jl_AL_id))
  NC_CALL(nc_enddef(ncidAL))

  /*
   * Define output file variables for land/ocean X-grid
   */
  NC_CALL(nc_create(lxo_file, NC_NOCLOBBER, &ncidLO))
  NC_CALL(nc_def_dim(ncidLO, "i_lndXocn", NC_UNLIMITED, &dim_LO))
  NC_CALL(nc_def_var(ncidLO, "AREA_LNDxOCN", NC_DOUBLE,1, &dim_LO, &area_LO_id))
  NC_CALL(nc_def_var(ncidLO, "DI_LNDxOCN", NC_DOUBLE,1, &dim_LO, &di_LO_id))
  NC_CALL(nc_def_var(ncidLO, "DJ_LNDxOCN", NC_DOUBLE,1, &dim_LO, &dj_LO_id))
  NC_CALL(nc_def_var(ncidLO, "I_LND_LNDxOCN",  NC_INT,  1, &dim_LO, &il_LO_id))
  NC_CALL(nc_def_var(ncidLO, "J_LND_LNDxOCN",  NC_INT,  1, &dim_LO, &jl_LO_id))
  NC_CALL(nc_def_var(ncidLO, "I_OCN_LNDxOCN",  NC_INT,  1, &dim_LO, &io_LO_id))
  NC_CALL(nc_def_var(ncidLO, "J_OCN_LNDxOCN",  NC_INT,  1, &dim_LO, &jo_LO_id))
  NC_CALL(nc_enddef(ncidLO))

  /* 
   * Find overlaps of atmosphere cells (the window) and surface cells
   * Find the information about land.
   */
  {
    int atmXlnd_n = 0;
    int atmXocn_n = 0;
    int atmXlnd_start, atmXlnd_n_prev;
    double *atmXlnd_area, *atmXlnd_int2, *atmXlnd_int3;
    double *atmXlnd_di, *atmXlnd_dj;
    double *atmXocn_area, *atmXocn_di, *atmXocn_dj;
    double *gal_area=0, *gal_di=0, *gal_dj=0;
    int    *gal_il=0, *gal_jl=0, *gal_ia=0, *gal_ja=0;
    double *gao_area=0, *gao_di=0, *gao_dj=0; 
    int    *gao_io=0, *gao_jo=0, *gao_ia=0, *gao_ja=0;    
    
    int    *atmXocn_io, *atmXocn_jo;
    int    *atmXocn_ia, *atmXocn_ja;
    int    *atmXlnd_il, *atmXlnd_jl;
    int    *atmXlnd_ia, *atmXlnd_ja;    
    VTX    *atmXlnd_ll, *atmXlnd_ur;


    atmXlnd_area  = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXlnd_di    = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXlnd_dj    = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXlnd_int2  = (double *) malloc (MAXLOCAL*sizeof(double));  
    atmXlnd_int3  = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXlnd_ll    = (VTX    *) malloc (MAXLOCAL*sizeof(   VTX));
    atmXlnd_ur    = (VTX    *) malloc (MAXLOCAL*sizeof(   VTX));
    atmXlnd_il    = (int    *) malloc (MAXLOCAL*sizeof(   int));
    atmXlnd_jl    = (int    *) malloc (MAXLOCAL*sizeof(   int));
    atmXlnd_ia    = (int    *) malloc (MAXLOCAL*sizeof(   int));
    atmXlnd_ja    = (int    *) malloc (MAXLOCAL*sizeof(   int));    
    atmXocn_area  = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXocn_di    = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXocn_dj    = (double *) malloc (MAXLOCAL*sizeof(double));
    atmXocn_io    = (int    *) malloc (MAXLOCAL*sizeof(   int));
    atmXocn_jo    = (int    *) malloc (MAXLOCAL*sizeof(   int));    
    atmXocn_ia    = (int    *) malloc (MAXLOCAL*sizeof(   int));
    atmXocn_ja    = (int    *) malloc (MAXLOCAL*sizeof(   int));   
    
    for (na = 0; na < natm_local; na ++) {
      VTX atm_ll, atm_ur;

      atmXlnd_start = atmXlnd_n;
      ia = (na+isc)%IMA;
      ja = (na+isc)/IMA;
    
      /* atmosphere cell clipping window */
      atm_ll.x = atmlonb[ia  ]; atm_ll.y = atmlatb[ja  ];
      atm_ur.x = atmlonb[ia+1]; atm_ur.y = atmlatb[ja+1];


      for (jl=0;jl<JML;jl++) for (il=0;il<IML;il++) {
	VTX lnd_ll, lnd_ur;

	lnd_ll.x = lndlonb[il  ]; lnd_ll.y = lndlatb[jl  ];
	lnd_ur.x = lndlonb[il+1]; lnd_ur.y = lndlatb[jl+1];

	/* cyclic condition is considered */
        if(lnd_ll.x > atm_ur.x) {
	  lnd_ll.x -= 360;
	  lnd_ur.x -= 360;
	}
	if(atm_ll.x > lnd_ur.x) {
	  lnd_ll.x += 360;
	  lnd_ur.x += 360;
	}
	  
	atmXlnd_ll[atmXlnd_n].x = MAX(atm_ll.x,lnd_ll.x);
	atmXlnd_ur[atmXlnd_n].x = MIN(atm_ur.x,lnd_ur.x);
	if (atmXlnd_ll[atmXlnd_n].x>atmXlnd_ur[atmXlnd_n].x) continue;

	atmXlnd_ll[atmXlnd_n].y = MAX(atm_ll.y,lnd_ll.y);
	atmXlnd_ur[atmXlnd_n].y = MIN(atm_ur.y,lnd_ur.y);
	if (atmXlnd_ll[atmXlnd_n].y>atmXlnd_ur[atmXlnd_n].y) continue;

	atmXlnd_area[atmXlnd_n] = box_area(atmXlnd_ll[atmXlnd_n],
					   atmXlnd_ur[atmXlnd_n]);
	atmXlnd_int2[atmXlnd_n] = box_int2(atmXlnd_ll[atmXlnd_n],
					   atmXlnd_ur[atmXlnd_n], 0.5*(atm_ll.x+atm_ur.x));
	atmXlnd_int3[atmXlnd_n] = box_int3(atmXlnd_ll[atmXlnd_n],
					   atmXlnd_ur[atmXlnd_n]);      
	atmXlnd_il[atmXlnd_n] = il+1;
	atmXlnd_jl[atmXlnd_n] = jl+1;
	atmXlnd_ia[atmXlnd_n] = ia+1;
	atmXlnd_ja[atmXlnd_n] = ja+1;      
	atmXlnd_n++;
	if(atmXlnd_n > MAXLOCAL) error_handler("atmXlnd_n is greater than MAXLOCAL, increase MAXLOCAL"); 
      };
  
      /* calculate atmos/ocean x-cells */
      for (jo=0;jo<JMO;jo++) for (io=0;io<IMO;io++) if (wet[jo*IMO+io]>0.5) {
	int n_in, n_out;
	VTX v_in[MV], v_out[MV];
	double Xarea, Xctrlon, Xctrlat;

	/* list ocean cell vertices in counter-clockwise order */
	v_in[0].y = ocnlatb[ jo   *(IMO+1)+io  ]; /* SW - y */
	v_in[1].y = ocnlatb[ jo   *(IMO+1)+io+1]; /* SE - y */
	v_in[2].y = ocnlatb[(jo+1)*(IMO+1)+io+1]; /* NE - y */
	v_in[3].y = ocnlatb[(jo+1)*(IMO+1)+io  ]; /* NW - y */

	if (  (v_in[0].y<=atm_ll.y) && (v_in[1].y<=atm_ll.y)
	      && (v_in[2].y<=atm_ll.y) && (v_in[3].y<=atm_ll.y) ) continue;

	if (  (v_in[0].y>=atm_ur.y) && (v_in[1].y>=atm_ur.y)
	      && (v_in[2].y>=atm_ur.y) && (v_in[3].y>=atm_ur.y) ) continue;

	v_in[0].x = ocnlonb[ jo   *(IMO+1)+io  ]; /* SW - x */
	v_in[1].x = ocnlonb[ jo   *(IMO+1)+io+1]; /* SE - x */
	v_in[2].x = ocnlonb[(jo+1)*(IMO+1)+io+1]; /* NE - x */
	v_in[3].x = ocnlonb[(jo+1)*(IMO+1)+io  ]; /* NW - x */

	n_in = lon_fix(v_in, 4, (atm_ll.x+atm_ur.x)/2);
   
	if (  ((n_out = clip ( v_in, n_in, atm_ll, atm_ur, v_out )) > 0)
	      && ((Xarea = poly_area ( v_out, n_out ))              > EPS1) ) {
        
	  if (ocn_area[jo*IMO+io] < EPS1) {
	    fprintf(stderr, "\nOcean cell at i=%d j=%d has invalid area!\n", io, jo);
	    fprintf(stderr,"\n");
	    exit(1);
	  }        
	  if (atm_area[na] < EPS1) {
	    fprintf(stderr, "\nAtmosphter cell at i=%d j=%d has invalid area!\n", ia, ja);
	    fprintf(stderr,"\n");
	    exit(1);
	  }        	
	  Xctrlon = poly_ctrlon ( v_out, n_out, 0.5*(atm_ll.x + atm_ur.x));
	  Xctrlat = poly_ctrlat ( v_out, n_out );
	  atmXocn_di[atmXocn_n] = Xctrlon - atm_ctrlon[na];
	  atmXocn_dj[atmXocn_n] = Xctrlat - atm_ctrlat[na];
	  atmXocn_area[atmXocn_n] = Xarea;
	  atmXocn_io[atmXocn_n] = io+1;
	  atmXocn_jo[atmXocn_n] = jo+1;
	  atmXocn_ia[atmXocn_n] = ia+1;
	  atmXocn_ja[atmXocn_n] = ja+1;
	  atmXocn_n++;	
	  if(atmXocn_n > MAXLOCAL) error_handler("atmXocn_n is greater than MAXLOCAL, increase MAXLOCAL"); 
	  ocn_area_left[jo*IMO+io] -= Xarea;
	  atm_area_left[na] -= Xarea;
	  tot_area_AXO += Xarea;
	}

	/* remove ocean parts from atmos/land x-cells */
	for (i=atmXlnd_start;i<atmXlnd_n;i++) 
	  if (((n_out = clip(v_in, n_in, atmXlnd_ll[i], atmXlnd_ur[i], v_out)) > 0) &&
	      ( (Xarea = poly_area ( v_out, n_out )) > EPS1) ) {
	    atmXlnd_area[i] -= Xarea;
	    atmXlnd_int2[i] -= poly_int2 ( v_out, n_out, 0.5*(atm_ll.x+ atm_ur.x));
	    atmXlnd_int3[i] -= poly_int3 ( v_out, n_out );
	  }

      } /* ocean loop */

      atmXlnd_n_prev = atmXlnd_n;
      atmXlnd_n = atmXlnd_start;
      for (i=atmXlnd_start;i<atmXlnd_n_prev;i++) if (atmXlnd_area[i]>EPS1) {

	/*calculate the centroid of amt/lnd xgrid. */

	atmXlnd_di[atmXlnd_n] = atmXlnd_int2[i] /atmXlnd_area[i] - atm_ctrlon[na];
	atmXlnd_dj[atmXlnd_n] = atmXlnd_int3[i] /atmXlnd_area[i] - atm_ctrlat[na];
      
	atm_area_left[na]                             -= atmXlnd_area[i];
	lnd_area[(atmXlnd_jl[i]-1)*IML+atmXlnd_il[i]-1] += atmXlnd_area[i];
	tot_area_AXL                                  += atmXlnd_area[i];
	atmXlnd_area[atmXlnd_n] = atmXlnd_area[i];
	atmXlnd_il[atmXlnd_n]   = atmXlnd_il[i];
	atmXlnd_jl[atmXlnd_n]   = atmXlnd_jl[i];
	atmXlnd_ia[atmXlnd_n]   = atmXlnd_ia[i];
	atmXlnd_ja[atmXlnd_n]   = atmXlnd_ja[i];
	atmXlnd_n++;
      }
      if( npes == 1 ) { /* print this information in parallel run is complicated */
	if (ia+1==IMA) printf( "to %6.2f: %d AL x-cells (area=%g); %d AO x-cells (area=%g)\n",
			       atmlatb[ja+1], atmXlnd_n, tot_area_AXL,
			       atmXocn_n, tot_area_AXO);
      }
    } /* atmosphere loop */

    iAXL = get_global_count( atmXlnd_n );
    iAXO = get_global_count( atmXocn_n );    
    sum_over_pe(&tot_area_AXL, 1);
    sum_over_pe(&tot_area_AXO, 1);
    if(pe == root_pe) printf( "\n%d AL x-cells (area=%g); %d AO x-cells (area=%g)\n\n",
			      (int)iAXL, tot_area_AXL, (int)iAXO, tot_area_AXO );
    
    /*write out atm/land xgrid, first need to send the data to root pe.
      write out only from root pe*/

    if ( pe == root_pe ) {
      gal_area = (double *) malloc( iAXL * sizeof(double) );
      gal_di   = (double* ) malloc( iAXL * sizeof(double) );
      gal_dj   = (double* ) malloc( iAXL * sizeof(double) );
      gal_ia   = (int*    ) malloc( iAXL * sizeof(int) );  
      gal_ja   = (int*    ) malloc( iAXL * sizeof(int) );      
      gal_il   = (int*    ) malloc( iAXL * sizeof(int) );  
      gal_jl   = (int*    ) malloc( iAXL * sizeof(int) ); 
    }
    get_global_data_double( atmXlnd_area, gal_area, atmXlnd_n );
    get_global_data_double( atmXlnd_di, gal_di, atmXlnd_n );	
    get_global_data_double( atmXlnd_dj, gal_dj, atmXlnd_n );
    get_global_data_int( atmXlnd_il, gal_il, atmXlnd_n );
    get_global_data_int( atmXlnd_jl, gal_jl, atmXlnd_n );
    get_global_data_int( atmXlnd_ia, gal_ia, atmXlnd_n );
    get_global_data_int( atmXlnd_ja, gal_ja, atmXlnd_n );

    /* write out from root pe */
    if( pe == root_pe ) {
      size_t start, nwrite;
      start = 0; nwrite = iAXL;
      nc_put_vara_double( ncidAL, area_AL_id, &start, &nwrite, gal_area );
      nc_put_vara_double( ncidAL, di_AL_id, &start, &nwrite, gal_di );
      nc_put_vara_double( ncidAL, dj_AL_id, &start, &nwrite, gal_dj );
      nc_put_vara_int( ncidAL, ia_AL_id, &start, &nwrite, gal_ia );
      nc_put_vara_int( ncidAL, ja_AL_id, &start, &nwrite, gal_ja );
      nc_put_vara_int( ncidAL, il_AL_id, &start, &nwrite, gal_il );
      nc_put_vara_int( ncidAL, jl_AL_id, &start, &nwrite, gal_jl );
      /* release memory */
      free(gal_area);
      free(gal_di);
      free(gal_dj);
      free(gal_ia);
      free(gal_ja);
      free(gal_il);
      free(gal_jl);
    }

    /*write out atm/ocn xgrid, first need to send the data to root pe.
      write out only from root pe*/    
    if ( pe == root_pe ) {
      gao_area = (double *) malloc( iAXO * sizeof(double) );
      gao_di   = (double* ) malloc( iAXO * sizeof(double) );
      gao_dj   = (double* ) malloc( iAXO * sizeof(double) );      
      gao_io   = (int*    ) malloc( iAXO * sizeof(int) );  
      gao_jo   = (int*    ) malloc( iAXO * sizeof(int) ); 
      gao_ia   = (int*    ) malloc( iAXO * sizeof(int) );  
      gao_ja   = (int*    ) malloc( iAXO * sizeof(int) );
    }
    get_global_data_double( atmXocn_area, gao_area, atmXocn_n );
    get_global_data_double( atmXocn_di, gao_di, atmXocn_n );	
    get_global_data_double( atmXocn_dj, gao_dj, atmXocn_n );
    get_global_data_int( atmXocn_ia, gao_ia, atmXocn_n );
    get_global_data_int( atmXocn_ja, gao_ja, atmXocn_n );
    get_global_data_int( atmXocn_io, gao_io, atmXocn_n );
    get_global_data_int( atmXocn_jo, gao_jo, atmXocn_n );

    /* write out from root pe */
    if( pe == root_pe ) {
      size_t start, nwrite;
      start = 0; nwrite = iAXO;
      nc_put_vara_double( ncidAO, area_AO_id, &start, &nwrite, gao_area );
      nc_put_vara_double( ncidAO, di_AO_id, &start, &nwrite, gao_di );
      nc_put_vara_double( ncidAO, dj_AO_id, &start, &nwrite, gao_dj );
      nc_put_vara_int( ncidAO, ia_AO_id, &start, &nwrite, gao_ia );
      nc_put_vara_int( ncidAO, ja_AO_id, &start, &nwrite, gao_ja );
      nc_put_vara_int( ncidAO, io_AO_id, &start, &nwrite, gao_io );
      nc_put_vara_int( ncidAO, jo_AO_id, &start, &nwrite, gao_jo );
      
      if(land_same_as_atmos == 1) {
	iLXO = iAXO;
	nc_put_vara_double( ncidLO, area_LO_id, &start, &nwrite, gao_area );
	nc_put_vara_double( ncidLO, di_LO_id, &start, &nwrite, gao_di );
	nc_put_vara_double( ncidLO, dj_LO_id, &start, &nwrite, gao_dj );
	nc_put_vara_int( ncidLO, il_LO_id, &start, &nwrite, gao_ia );
	nc_put_vara_int( ncidLO, jl_LO_id, &start, &nwrite, gao_ja );
	nc_put_vara_int( ncidLO, io_LO_id, &start, &nwrite, gao_io );
	nc_put_vara_int( ncidLO, jo_LO_id, &start, &nwrite, gao_jo );	
      }
      
      /* release memory */
      free(gao_area);
      free(gao_di);
      free(gao_dj);
      free(gao_ia);
      free(gao_ja);
      free(gao_io);
      free(gao_jo);
    }
    /* release memory */
    free(atmXocn_area);
    free(atmXocn_di);
    free(atmXocn_dj);
    free(atmXocn_ia);
    free(atmXocn_ja);
    free(atmXocn_io);
    free(atmXocn_jo);
    free(atmXlnd_area);
    free(atmXlnd_di);
    free(atmXlnd_dj);
    free(atmXlnd_ia);
    free(atmXlnd_ja);
    free(atmXlnd_il);
    free(atmXlnd_jl); 
    /* get global lnd_area by summing over all the pe. */
    sum_over_pe(lnd_area, IML*JML);  
  }
    /* 
     * Find overlaps of land cells (the window) and ocean cells for runoff x-grid
     * when land and atmosphere are on different grid.
     */
    if( land_same_as_atmos == 0 ) {
      int lndXocn_n=0;
      int nlnd_local, nl;
      double *lndXocn_area, *lndXocn_di, *lndXocn_dj;
      double *glo_area, *glo_di, *glo_dj;
      int    *lndXocn_il, *lndXocn_jl, *lndXocn_io, *lndXocn_jo;
      int    *glo_il, *glo_jl, *glo_io, *glo_jo;
      
      define_domain( IML*JML, npes, startlist, endlist);
      isc = startlist[pe];
      iec = endlist[pe];
      nlnd_local = iec - isc + 1;

      if(nlnd_local > MAXLOCAL) error_handler("number of local land points is greater than MAXLOCAL, increase MAXLOCAL");
      
      lndXocn_area  = (double *) malloc (MAXLOCAL*sizeof(double));
      lndXocn_di    = (double *) malloc (MAXLOCAL*sizeof(double));
      lndXocn_dj    = (double *) malloc (MAXLOCAL*sizeof(double));
      lndXocn_io    = (int    *) malloc (MAXLOCAL*sizeof(   int));
      lndXocn_jo    = (int    *) malloc (MAXLOCAL*sizeof(   int));    
      lndXocn_il    = (int    *) malloc (MAXLOCAL*sizeof(   int));
      lndXocn_jl    = (int    *) malloc (MAXLOCAL*sizeof(   int));   
      
      for (nl=0;nl<nlnd_local;nl++) {
	VTX lnd_ll, lnd_ur;
	il = (nl+isc)%IML;
	jl = (nl+isc)/IML;
	
	/* land cell clipping window */
	lnd_ll.x = lndlonb[il  ]; lnd_ll.y = lndlatb[jl  ];
	lnd_ur.x = lndlonb[il+1]; lnd_ur.y = lndlatb[jl+1];

	for (jo=0;jo<JMO;jo++) for (io=0;io<IMO;io++) if (wet[jo*IMO+io]>0.5) {
	  int n_in, n_out;
	  VTX v_in[MV], v_out[MV];
	  double Xarea, Xctrlon, Xctrlat;

	  /* list ocean cell vertices in counter-clockwise order */
	  v_in[0].y = ocnlatb[ jo   *(IMO+1)+io  ]; /* SW - y */
	  v_in[1].y = ocnlatb[ jo   *(IMO+1)+io+1]; /* SE - y */
	  v_in[2].y = ocnlatb[(jo+1)*(IMO+1)+io+1]; /* NE - y */
	  v_in[3].y = ocnlatb[(jo+1)*(IMO+1)+io  ]; /* NW - y */

	  if (  (v_in[0].y<=lnd_ll.y) && (v_in[1].y<=lnd_ll.y)
		&& (v_in[2].y<=lnd_ll.y) && (v_in[3].y<=lnd_ll.y) ) continue;

	  if (  (v_in[0].y>=lnd_ur.y) && (v_in[1].y>=lnd_ur.y)
		&& (v_in[2].y>=lnd_ur.y) && (v_in[3].y>=lnd_ur.y) ) continue;

	  v_in[0].x = ocnlonb[ jo   *(IMO+1)+io  ]; /* SW - x */
	  v_in[1].x = ocnlonb[ jo   *(IMO+1)+io+1]; /* SE - x */
	  v_in[2].x = ocnlonb[(jo+1)*(IMO+1)+io+1]; /* NE - x */
	  v_in[3].x = ocnlonb[(jo+1)*(IMO+1)+io  ]; /* NW - x */

	  n_in = lon_fix(v_in, 4, (lnd_ll.x+lnd_ur.x)/2);
   
	  if (  ((n_out = clip ( v_in, n_in, lnd_ll, lnd_ur, v_out )) > 0)
		&& ((Xarea = poly_area ( v_out, n_out ))              > EPS1) ) {

	    /* calculate the centroid of lnd/ocn runoff */
	    Xctrlon = poly_ctrlon ( v_out, n_out, 0.5*(lnd_ll.x+lnd_ur.x ));
	    Xctrlat = poly_ctrlat ( v_out, n_out );
	    lndXocn_area[lndXocn_n] = Xarea;
	    lndXocn_di[lndXocn_n]   = Xctrlon - lnd_cell_ctrlon[jl*IML+il];
	    lndXocn_dj[lndXocn_n]   = Xctrlat - lnd_cell_ctrlat[jl*IML+il];
            lndXocn_il[lndXocn_n]   = il + 1;
	    lndXocn_jl[lndXocn_n]   = jl + 1;
            lndXocn_io[lndXocn_n]   = io + 1;
	    lndXocn_jo[lndXocn_n]   = jo + 1;
	    lndXocn_n++;
	    ocn_minus_runoff[jo*IMO+io] -= Xarea;
	  }
	} /* ocean loop */
      } /* land loop */
      iLXO = get_global_count( lndXocn_n );
      if ( pe == root_pe ) {
	glo_area = (double *) malloc( iLXO * sizeof(double) );
	glo_di   = (double* ) malloc( iLXO * sizeof(double) );
	glo_dj   = (double* ) malloc( iLXO * sizeof(double) );      
	glo_il   = (int*    ) malloc( iLXO * sizeof(int) );  
	glo_jl   = (int*    ) malloc( iLXO * sizeof(int) ); 
	glo_io   = (int*    ) malloc( iLXO * sizeof(int) );  
	glo_jo   = (int*    ) malloc( iLXO * sizeof(int) );
      }
      get_global_data_double( lndXocn_area, glo_area, lndXocn_n );
      get_global_data_double( lndXocn_di, glo_di, lndXocn_n );	
      get_global_data_double( lndXocn_dj, glo_dj, lndXocn_n );
      get_global_data_int( lndXocn_il, glo_il, lndXocn_n );
      get_global_data_int( lndXocn_jl, glo_jl, lndXocn_n );
      get_global_data_int( lndXocn_io, glo_io, lndXocn_n );
      get_global_data_int( lndXocn_jo, glo_jo, lndXocn_n );
      /* write out from root pe */
      if( pe == root_pe ) {
	size_t start, nwrite;
	start = 0; nwrite = iLXO;
	nc_put_vara_double( ncidLO, area_LO_id, &start, &nwrite, glo_area );
	nc_put_vara_double( ncidLO, di_LO_id, &start, &nwrite, glo_di );
	nc_put_vara_double( ncidLO, dj_LO_id, &start, &nwrite, glo_dj );
	nc_put_vara_int( ncidLO, il_LO_id, &start, &nwrite, glo_il );
	nc_put_vara_int( ncidLO, jl_LO_id, &start, &nwrite, glo_jl );
	nc_put_vara_int( ncidLO, io_LO_id, &start, &nwrite, glo_io );
	nc_put_vara_int( ncidLO, jo_LO_id, &start, &nwrite, glo_jo );
	/* release memory */
	free(glo_area);
	free(glo_di);
	free(glo_dj);
	free(glo_io);
	free(glo_jo);
	free(glo_il);
	free(glo_jl);
      }
      free(lndXocn_area);
      free(lndXocn_di);
      free(lndXocn_dj);
      free(lndXocn_il);
      free(lndXocn_jl);
      free(lndXocn_io);
      free(lndXocn_jo);
      printf( "Land/Ocean x-grid complete (%d x-cells)\n", (int) iLXO);
    }
    
  { double a_err = 0.0, o_err = 0.0;
    double a_max_err = 0.0, o_max_err = 0.0;
    int    o_cnt = 0;
    double *ga_area=NULL, *ga_area_left=NULL;

    ga_area = (double *)malloc(IMA*JMA*sizeof(double) ); 
    ga_area_left = (double *)malloc(IMA*JMA*sizeof(double) );
    get_global_data_double( atm_area_left, ga_area_left, natm_local );
    get_global_data_double( atm_area, ga_area, natm_local );

    if( pe == root_pe ) {
      for (i=0;i<IMA*JMA;i++) {
	double err  = ga_area_left[i]/ga_area[i];
	double err2 = err*err;

	if (err2>a_max_err*a_max_err) {
	  a_max_err = -err;
	  ja_max_err = i/IMA;
	  ia_max_err = i-ja_max_err*IMA;
	}

	a_err += err2;
      }
      printf( "\nRMS Atmosphere cell tiling error = %lg\n",
	      sqrt(a_err/(IMA*JMA)));
      printf( "MAX Atmosphere cell tiling error = %lg ", a_max_err);
      printf( "at (i=%d,j=%d), (x=%6.2f,y=%6.2f)\n",
	      ia_max_err+1, ja_max_err+1,
	      (atmlonb[ia_max_err]+atmlonb[ia_max_err+1])/2,
	      (atmlatb[ja_max_err]+atmlatb[ja_max_err+1])/2);
    }
    free(ga_area);
    free(ga_area_left);
    
    if (ofile) {
      if( npes > 1 ) { /* need to sum over all the pe and minus (n-1)*ocn_area */
	sum_over_pe( ocn_area_left, IMO*JMO );			  
        for(i=0; i<IMO*JMO; i++) ocn_area_left[i] -= ocn_area[i]*(npes-1);
      }

      if( pe == root_pe ) {
      for (i=0;i<IMO*JMO;i++) if (wet[i]>0.5) {
	double err  = ocn_area_left[i]/ocn_area[i];
	double err2 = err*err;
	if (err2>o_max_err*o_max_err) {
	  o_max_err = -err;
	  jo_max_err = i/IMO;
	  io_max_err = i-jo_max_err*IMO;
	}
	o_err += err2;
	o_cnt++;
      }      
      if(o_cnt>0) printf( "\nRMS Ocean cell tiling error (atmos/ocean x-cells) = %lg\n", sqrt(o_err/o_cnt));
      printf( "MAX Ocean cell tiling error (atmos/ocean x-cells) = %lg ", o_max_err);
      printf( "at (i=%d,j=%d), (x=%6.2f,y=%6.2f)\n",
                    io_max_err+1, jo_max_err+1,
                    (ocnlonb[ jo_max_err   *(IMO+1)+io_max_err  ]
                    +ocnlonb[ jo_max_err   *(IMO+1)+io_max_err+1]
                    +ocnlonb[(jo_max_err+1)*(IMO+1)+io_max_err+1]
                    +ocnlonb[(jo_max_err+1)*(IMO+1)+io_max_err  ])/4,
                    (ocnlatb[ jo_max_err   *(IMO+1)+io_max_err  ]
                    +ocnlatb[ jo_max_err   *(IMO+1)+io_max_err+1]
                    +ocnlatb[(jo_max_err+1)*(IMO+1)+io_max_err+1]
                    +ocnlatb[(jo_max_err+1)*(IMO+1)+io_max_err  ])/4);
      }
      if(land_same_as_atmos == 0) {
	if( npes > 1 ) { /* need to sum over all the pe and minus (n-1)*ocn_area */
	  sum_over_pe( ocn_minus_runoff, IMO*JMO );			  
	  for(i=0; i<IMO*JMO; i++) ocn_minus_runoff[i] -= ocn_area[i]*(npes-1);
	}

	if( pe == root_pe ) {
	  o_max_err = 0.0;
	  o_err =0.0;
	  o_cnt = 0;
    
	  for (i=0;i<IMO*JMO;i++) if (wet[i]>0.5) {
	    double err  = ocn_minus_runoff[i]/ocn_area[i];
	    double err2 = err*err;

	    if (err2>o_max_err*o_max_err) {
	      o_max_err = -err;
	      jo_max_err = i/IMO;
	      io_max_err = i-jo_max_err*IMO;
	    }
	    o_err += err2;
	    o_cnt++;
	  }
	}
      }
      
      if( pe == root_pe ) {
	printf( "MAX Ocean cell tiling error (land/ocean runoff) = %lg ", o_max_err);
	printf( "at (i=%d,j=%d), (x=%6.2f,y=%6.2f)\n",
		io_max_err+1, jo_max_err+1,
		(ocnlonb[ jo_max_err   *(IMO+1)+io_max_err  ]
		 +ocnlonb[ jo_max_err   *(IMO+1)+io_max_err+1]
		 +ocnlonb[(jo_max_err+1)*(IMO+1)+io_max_err+1]
		 +ocnlonb[(jo_max_err+1)*(IMO+1)+io_max_err  ])/4,
		(ocnlatb[ jo_max_err   *(IMO+1)+io_max_err  ]
		 +ocnlatb[ jo_max_err   *(IMO+1)+io_max_err+1]
		 +ocnlatb[(jo_max_err+1)*(IMO+1)+io_max_err+1]
		 +ocnlatb[(jo_max_err+1)*(IMO+1)+io_max_err  ])/4);
      }
    } 
  }
  
  if( pe == root_pe ) {

    printf("\nGlobe tiling error = %lg\n", tot_area_AXL+tot_area_AXO-1);
    printf("\nMax # x-cell vertices = %d\n", n_max);
  }
  
  if( pe == root_pe ) {
    char cmd[BUFSIZ];
    int  ncid;
    int area_AO_id2, ia_AO_id2, ja_AO_id2, io_AO_id2, jo_AO_id2;
    int di_AO_id2, dj_AO_id2, di_AL_id2, dj_AL_id2, di_LO_id2, dj_LO_id2;
    int area_AL_id2, ia_AL_id2, ja_AL_id2, il_AL_id2, jl_AL_id2;
    int area_LO_id2, il_LO_id2, jl_LO_id2, io_LO_id2, jo_LO_id2;
    int dim_xba, dim_yba, dim_xbl, dim_ybl, xba_id, yba_id, xbl_id, ybl_id;
    int dim_xta, dim_yta, xta_id, yta_id, atm_area_id, a_dims[2];
    int dim_xtl, dim_ytl, xtl_id, ytl_id, lnd_area_id, l_dims[2];
    int lnd_cell_area_id, dim_null, null_id;
    int dim_xto, dim_yto, xto_id, yto_id, ocn_area_id, o_dims[2];
    size_t i;

    if (ofile && nofile == 1) {
      sprintf(cmd, "cp %s %s", ofile, gfile);
      system(cmd);
      NC_CALL(nc_open(gfile, NC_WRITE, &ncid))
      NC_CALL(nc_redef(ncid))
    } else {
      NC_CALL(nc_create(gfile, NC_NOCLOBBER, &ncid))
    }

    dim_null = -1;
    
    if (ofile) {
      if(iAXO ==0) {
	NC_CALL(nc_def_dim(ncid, "null_dimension", iAXO, &dim_null))  
	NC_CALL(nc_def_var(ncid, "null_dimension", NC_DOUBLE,1,&dim_null, &null_id))
	dim_AO = dim_null;
      }
      else {
	NC_CALL(nc_def_dim(ncid, "i_atmXocn", iAXO, &dim_AO)) 
      }
      NC_CALL(nc_def_var(ncid, "AREA_ATMxOCN",NC_DOUBLE,1,&dim_AO,&area_AO_id2))
      NC_CALL(nc_def_var(ncid, "DI_ATMxOCN",NC_DOUBLE,1,&dim_AO,&di_AO_id2))
      NC_CALL(nc_def_var(ncid, "DJ_ATMxOCN",NC_DOUBLE,1,&dim_AO,&dj_AO_id2))
      NC_CALL(nc_def_var(ncid, "I_ATM_ATMxOCN", NC_INT, 1, &dim_AO, &ia_AO_id2))
      NC_CALL(nc_def_var(ncid, "J_ATM_ATMxOCN", NC_INT, 1, &dim_AO, &ja_AO_id2))
      NC_CALL(nc_def_var(ncid, "I_OCN_ATMxOCN", NC_INT, 1, &dim_AO, &io_AO_id2))
      NC_CALL(nc_def_var(ncid, "J_OCN_ATMxOCN", NC_INT, 1, &dim_AO, &jo_AO_id2))
    }

    if(iAXL ==0) {
      if(dim_null == -1) {
	NC_CALL(nc_def_dim(ncid, "null_dimension", iAXL, &dim_null))  
	NC_CALL(nc_def_var(ncid, "null_dimension", NC_DOUBLE,1,&dim_null, &null_id))
	dim_AL = dim_null;	
      }
      else {
        dim_AL = dim_null;
      }
    }
    else {
      NC_CALL(nc_def_dim(ncid, "i_atmXlnd", iAXL, &dim_AL))
    }
     
    NC_CALL(nc_def_var(ncid, "AREA_ATMxLND", NC_DOUBLE,1,&dim_AL, &area_AL_id2))
    NC_CALL(nc_def_var(ncid, "DI_ATMxLND", NC_DOUBLE,1,&dim_AL, &di_AL_id2))
    NC_CALL(nc_def_var(ncid, "DJ_ATMxLND", NC_DOUBLE,1,&dim_AL, &dj_AL_id2))
    NC_CALL(nc_def_var(ncid, "I_ATM_ATMxLND",  NC_INT,  1, &dim_AL, &ia_AL_id2))
    NC_CALL(nc_def_var(ncid, "J_ATM_ATMxLND",  NC_INT,  1, &dim_AL, &ja_AL_id2))
    NC_CALL(nc_def_var(ncid, "I_LND_ATMxLND",  NC_INT,  1, &dim_AL, &il_AL_id2))
    NC_CALL(nc_def_var(ncid, "J_LND_ATMxLND",  NC_INT,  1, &dim_AL, &jl_AL_id2))

    if(iLXO ==0) {
      if(dim_null == -1) {
	NC_CALL(nc_def_dim(ncid, "null_dimension", iLXO, &dim_null))  
	NC_CALL(nc_def_var(ncid, "null_dimension", NC_DOUBLE,1,&dim_null, &null_id))
	dim_LO = dim_null;	
      }
      else {
        dim_LO = dim_null;
      }
    }
    else {
      NC_CALL(nc_def_dim(ncid, "i_lndXocn", iLXO, &dim_LO))
    }

    NC_CALL(nc_def_var(ncid, "AREA_LNDxOCN", NC_DOUBLE,1,&dim_LO, &area_LO_id2))
    NC_CALL(nc_def_var(ncid, "DI_LNDxOCN", NC_DOUBLE,1,&dim_LO, &di_LO_id2))
    NC_CALL(nc_def_var(ncid, "DJ_LNDxOCN", NC_DOUBLE,1,&dim_LO, &dj_LO_id2))
    NC_CALL(nc_def_var(ncid, "I_LND_LNDxOCN",  NC_INT,  1, &dim_LO, &il_LO_id2))
    NC_CALL(nc_def_var(ncid, "J_LND_LNDxOCN",  NC_INT,  1, &dim_LO, &jl_LO_id2))
    NC_CALL(nc_def_var(ncid, "I_OCN_LNDxOCN",  NC_INT,  1, &dim_LO, &io_LO_id2))
    NC_CALL(nc_def_var(ncid, "J_OCN_LNDxOCN",  NC_INT,  1, &dim_LO, &jo_LO_id2))

    NC_CALL(nc_def_dim(ncid, "xba", IMA+1, &dim_xba))
    NC_CALL(nc_def_var(ncid, "xba", NC_DOUBLE, 1, &dim_xba, &xba_id))
    NC_CALL(nc_def_dim(ncid, "yba", JMA+1, &dim_yba))
    NC_CALL(nc_def_var(ncid, "yba", NC_DOUBLE, 1, &dim_yba, &yba_id))

    NC_CALL(nc_def_dim(ncid, "xta", IMA, &dim_xta))
    NC_CALL(nc_def_var(ncid, "xta", NC_DOUBLE, 1, &dim_xta, &xta_id))
    NC_CALL(nc_def_dim(ncid, "yta", JMA, &dim_yta))
    NC_CALL(nc_def_var(ncid, "yta", NC_DOUBLE, 1, &dim_yta, &yta_id))
    a_dims[1] = dim_xta; a_dims[0] = dim_yta;
    NC_CALL(nc_def_var(ncid, "AREA_ATM", NC_DOUBLE, 2, a_dims, &atm_area_id))
      
    NC_CALL(nc_def_dim(ncid, "xbl", IML+1, &dim_xbl))
    NC_CALL(nc_def_var(ncid, "xbl", NC_DOUBLE, 1, &dim_xbl, &xbl_id))
    NC_CALL(nc_def_dim(ncid, "ybl", JML+1, &dim_ybl))
    NC_CALL(nc_def_var(ncid, "ybl", NC_DOUBLE, 1, &dim_ybl, &ybl_id))

    NC_CALL(nc_def_dim(ncid, "xtl", IML, &dim_xtl))
    NC_CALL(nc_def_var(ncid, "xtl", NC_DOUBLE, 1, &dim_xtl, &xtl_id))
    NC_CALL(nc_def_dim(ncid, "ytl", JML, &dim_ytl))
    NC_CALL(nc_def_var(ncid, "ytl", NC_DOUBLE, 1, &dim_ytl, &ytl_id))
    l_dims[1] = dim_xtl; l_dims[0] = dim_ytl;
    NC_CALL(nc_def_var(ncid, "AREA_LND", NC_DOUBLE, 2, l_dims, &lnd_area_id))
    NC_CALL(nc_def_var(ncid, "AREA_LND_CELL", NC_DOUBLE, 2, l_dims, &lnd_cell_area_id))
    if (ofile) {
      NC_CALL(nc_def_dim(ncid, "xto", IMO, &dim_xto))
      NC_CALL(nc_def_var(ncid, "xto", NC_DOUBLE, 1, &dim_xto, &xto_id))
      NC_CALL(nc_def_dim(ncid, "yto", JMO, &dim_yto))
      NC_CALL(nc_def_var(ncid, "yto", NC_DOUBLE, 1, &dim_yto, &yto_id))
      o_dims[1] = dim_xto; o_dims[0] = dim_yto;
      NC_CALL(nc_def_var(ncid, "AREA_OCN", NC_DOUBLE, 2, o_dims, &ocn_area_id))
    }

    NC_CALL(nc_enddef(ncid))

    if (ofile) {
      nc_copy_1d_double(ncidAO, area_AO_id, ncid, area_AO_id2, (int) iAXO);
      nc_copy_1d_double(ncidAO, di_AO_id, ncid, di_AO_id2, (int) iAXO);
      nc_copy_1d_double(ncidAO, dj_AO_id, ncid, dj_AO_id2, (int) iAXO);
      nc_copy_1d_int(ncidAO, ia_AO_id, ncid, ia_AO_id2, (int) iAXO);
      nc_copy_1d_int(ncidAO, ja_AO_id, ncid, ja_AO_id2, (int) iAXO);
      nc_copy_1d_int(ncidAO, io_AO_id, ncid, io_AO_id2, (int) iAXO);
      nc_copy_1d_int(ncidAO, jo_AO_id, ncid, jo_AO_id2, (int) iAXO);
      NC_CALL(nc_close(ncidAO))
    }

    nc_copy_1d_double(ncidAL, area_AL_id, ncid, area_AL_id2, (int) iAXL);
    nc_copy_1d_double(ncidAL, di_AL_id, ncid, di_AL_id2, (int) iAXL);
    nc_copy_1d_double(ncidAL, dj_AL_id, ncid, dj_AL_id2, (int) iAXL);
    nc_copy_1d_int(ncidAL, ia_AL_id, ncid, ia_AL_id2, (int) iAXL);
    nc_copy_1d_int(ncidAL, ja_AL_id, ncid, ja_AL_id2, (int) iAXL);
    nc_copy_1d_int(ncidAL, il_AL_id, ncid, il_AL_id2, (int) iAXL);
    nc_copy_1d_int(ncidAL, jl_AL_id, ncid, jl_AL_id2, (int) iAXL);
    NC_CALL(nc_close(ncidAL))

    if (ofile) {
      nc_copy_1d_double(ncidLO, area_LO_id, ncid, area_LO_id2, (int) iLXO);
      nc_copy_1d_double(ncidLO, di_LO_id, ncid, di_LO_id2, (int) iLXO);
      nc_copy_1d_double(ncidLO, dj_LO_id, ncid, dj_LO_id2, (int) iLXO);
      nc_copy_1d_int(ncidLO, il_LO_id, ncid, il_LO_id2, (int) iLXO);
      nc_copy_1d_int(ncidLO, jl_LO_id, ncid, jl_LO_id2, (int) iLXO);
      nc_copy_1d_int(ncidLO, io_LO_id, ncid, io_LO_id2, (int) iLXO);
      nc_copy_1d_int(ncidLO, jo_LO_id, ncid, jo_LO_id2, (int) iLXO);
      NC_CALL(nc_close(ncidLO))
    }

    for (i=0;i<=IMA;i++) NC_CALL(nc_put_var1_double(ncid, xba_id, &i, &(atmlonb[i])));
    for (i=0;i<=JMA;i++) NC_CALL(nc_put_var1_double(ncid, yba_id, &i, &(atmlatb[i])));

    for (i=0;i<IMA;i++) 
      NC_CALL(nc_put_var1_double(ncid, xta_id, &i, &(atmlont[i])))

    for (i=0;i<JMA;i++)
      NC_CALL(nc_put_var1_double(ncid, yta_id, &i, &(atmlatt[i])))

    for (i=0;i<=IML;i++)
      NC_CALL(nc_put_var1_double(ncid, xbl_id, &i, &(lndlonb[i])))
    for (i=0;i<=JML;i++)
      NC_CALL(nc_put_var1_double(ncid, ybl_id, &i, &(lndlatb[i])))

    for (i=0;i<IML;i++)
      NC_CALL(nc_put_var1_double(ncid, xtl_id, &i, &(lndlont[i])))

    for (i=0;i<JML;i++) 
      NC_CALL(nc_put_var1_double(ncid, ytl_id, &i, &(lndlatt[i])))

    NC_CALL(nc_put_var_double(ncid, atm_area_id, gatm_area)) 
    NC_CALL(nc_put_var_double(ncid, lnd_area_id, lnd_area))        
    NC_CALL(nc_put_var_double(ncid, lnd_cell_area_id, lnd_cell_area))

    if (ofile) {
      NC_CALL(nc_put_var_double(ncid, ocn_area_id, ocn_area))
	}
    NC_CALL(nc_close(ncid))

    sprintf(cmd, "rm -f %s %s %s", axo_file, axl_file, lxo_file);
    system(cmd);
    printf("\n grid_spec.nc is generated successfully\n");
    
  }
#ifdef use_libMPI   
    MPI_Finalize();
#endif  

} /* main */
