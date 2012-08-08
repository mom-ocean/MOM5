#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <netcdf.h>

/*   This subroutine converts several fields in a file from density   */
/* to depth as a vertical coordinate.                                 */

/*   No information is extrapolated below the bottom depth.  Above    */
/* center of the topmost layer with more than EPSILON mass and below  */
/* the center of the bottommost layer with any mass, constant values  */
/* are used, otherwise profiles are linearly interpolated between     */
/* layer centers.                                                     */

/*   All the dimensions used by any variables must be present in the  */
/* first file.                                                        */

/* to compile                                                         */
/*                                                                    */
/* cc -O -o Resample_on_Z Resample_on_Z.c -lm -lnetcdf -I/usr/local/include -L/usr/local/lib */


#define MAX_VARS 20
#define MAX_FILES 200
#define FILE_NAME_SZ 200

void get_depths(char depth_file[], char depth_names[2][NC_MAX_NAME],
                double **lay_depth, double **int_depth, size_t *nz);
int read_field(int fn, int lev, char var_name[], double **array,
               size_t *nlay, size_t *nx, size_t *ny, double st[2]);
int read_axes(char var_name[], double **array, size_t *nz, size_t *nx, size_t *ny);

void write_output(char out_file[], char in_file[][FILE_NAME_SZ], char depth_file[],
     size_t nz, char var_names[MAX_VARS][NC_MAX_NAME], char time_names[2][NC_MAX_NAME], 
     char depth_names[2][NC_MAX_NAME], int *num_var,
     char arglist[]);
void write_field(int lev, int v, double *array);
void copy_time_vals(size_t lev);
void handle_error(int status, char msg[], int id);

// size_t nx, ny, nz;

char zc_name[NC_MAX_NAME], ze_name[NC_MAX_NAME];
char xdim_name[NC_MAX_NAME], ydim_name[NC_MAX_NAME];
int xid, yid;
double *z_in;

double missing_val[MAX_VARS];

int ncInid[MAX_FILES], ncOutid;
int time_out_id[2];
int time_in_id[2];
int var_id[NC_MAX_VARS];
int var_file[NC_MAX_VARS];
int static_var[MAX_VARS];
size_t var_size[MAX_VARS];
size_t var_count[NC_MAX_VARS][NC_MAX_VAR_DIMS];
char master_in_file[FILE_NAME_SZ];

int main(int argc, char *argv[]) {
  double *int_depth, *lay_depth;
  int nvar=0, n_tname = 0, n, m, i, j, k, lev, v, num_files = 0;
  size_t nx, ny, nxy, nz;
  int verbose = 0;
  int status;
  size_t nlay, nint;
  double st_eta[2], st_var[2];
  double *eta = NULL, *var_out = NULL, *var_in = NULL, *ev = NULL;
  double *elc = NULL;
//  double elc[100], depths[100];
  char var_names[MAX_VARS][NC_MAX_NAME], eta_name[NC_MAX_NAME];
  double EPSILON = 1.0e-10;
  double missing = -1.0e34;
  char time_names[2][NC_MAX_NAME], depth_names[2][NC_MAX_NAME];

  char depth_file[FILE_NAME_SZ], in_file[MAX_FILES][FILE_NAME_SZ], out_file[FILE_NAME_SZ];
  char arglist[2000];

  
  // The following block interprets the input fields.
  {
    if (argc < 6) {
      printf("%s requires at least 5 arguments (# = %d):\n\n"
          "useage:\n%s -d'Depth_file' -V:'variable_names'[,'vn2',...] -e'Eta_name' variables\n"
          "\t-T:'Time_varable_name1,Time_varable_name2' -o'Output_file' 'Input_file_name[s]'\n"
          "\nOptions:\n\t-v\t\t Verbose mode"
           "\n\t-z:'Z-center_name,Z-edge_name'"
         "\nvariables:\n\tVariable names from input file.\n",
          argv[0],argc-1,argv[0]);
      return -1;
    }
    
    strcpy(arglist,argv[0]);
    for (n=1;n<argc;n++) {strcat(arglist," "); strcat(arglist,argv[n]);}

    for (n=1;n<argc;n++) if (strcmp(argv[n],"-v")==0) {
      verbose = 1; printf("\nUsing verbose mode.\n\n");
    }

    strcpy(depth_file,""); strcpy(out_file,"");
    for (n=0;n<MAX_FILES;n++) strcpy(in_file[n],"");
    strcpy(zc_name,"zt"); strcpy(ze_name,"zw");
    strcpy(eta_name,""); strcpy(time_names[0],""); strcpy(time_names[1],"");

    for (n=1;n<argc;n++) {
      size_t sz;
      // if (verbose) printf("Arg %d: %s\n",n,argv[n]);
      if (strcmp(argv[n],"-v")==0) continue;
      else if (strncmp(argv[n],"-z:",3)==0) {
        // PARSE FOR depth names.
        m = 3; sz = strlen(argv[n]);
        for (i=m;i<sz;i++) if (argv[n][i] == ',') {break;}
        strncpy(zc_name,&argv[n][m],i-m); zc_name[i-m] = '\0';
        // if (verbose) printf("\tTime name %d is %s.\n",n_tname,time_names[n_tname]);
        m = i+1;
        if (i<sz) {
          strncpy(ze_name,&argv[n][m],sz-m); ze_name[sz-m] = '\0';
        }
      }
      else if (strncmp(argv[n],"-d",2)==0) strcpy(depth_file,&argv[n][2]);
      else if (strncmp(argv[n],"-o",2)==0) strcpy(out_file,&argv[n][2]);
      else if (strncmp(argv[n],"-e",2)==0) strcpy(eta_name,&argv[n][2]);
      else if (strncmp(argv[n],"-V:",3)==0) {
        // PARSE FOR variable names.
        m = 3; sz = strlen(argv[n]);
        for (;;) {
          for (i=m;i<sz;i++) if (argv[n][i] == ',') {break;}
          strncpy(var_names[nvar],&argv[n][m],i-m);
          var_names[nvar][i-m] = '\0';
          // if (verbose) printf("\tVariable %d is %s.\n",nvar,var_names[nvar]);
          nvar++;
          if (i==sz) break;
          m = i+1;
        }
      }
      else if (strncmp(argv[n],"-T:",3)==0) {
        // PARSE FOR time names.
        m = 3; sz = strlen(argv[n]);
        for (;;) {
          for (i=m;i<sz;i++) if (argv[n][i] == ',') {break;}
          strncpy(time_names[n_tname],&argv[n][m],i-m);
          time_names[n_tname][i-m] = '\0';
          // if (verbose) printf("\tTime name %d is %s.\n",n_tname,time_names[n_tname]);
          n_tname++;
          if (i==sz) break;
          m = i+1;
        }
      }
      else if (strncmp(argv[n],"-",1)==0)
        printf("Unrecognized argument %s.\n",argv[n]);
      else {
        // Add input file name.
        if (num_files >= MAX_FILES)
          printf("Unable to add input file %s, because %d files are already in use.\n"
                 "Increase MAX_FILES in %s.c.\n",argv[n],num_files,argv[0]);
        strcpy(in_file[num_files],argv[n]);
        // if (verbose) printf("Adding input file %d with name %s.\n",num_files,in_file[num_files]);
        num_files++;
      }
    }

    if (strlen(depth_file) == 0) {
      printf("Depth file name must be specified as -dDEPTH_FILE_NAME\n");
      exit(-1);
    }
    if (num_files == 0) {
      printf("At least one input file name must be specified.\n");
      exit(-1);
    }
    if (strlen(out_file) == 0) {
      printf("Output file name must be specified as -oOUTPUT_FILE_NAME\n");
      exit(-1);
    }
    if (strlen(eta_name) == 0) {
      printf("Interface height variable name must be specified as -eETA_NAME\n");
      exit(-1);
    }
    if (nvar < 1) {
      printf("At least 1 variable must be specified.\n");
      exit(-1);
    }
  }
      
  strcpy(depth_names[0],zc_name);
  strcpy(depth_names[1],ze_name);
  //  get_vertical_grid(...);
  if (verbose) printf("\tGet depths from %s.\n",depth_file);
  get_depths(depth_file, depth_names, &lay_depth, &int_depth, &nz);

  
//  if (verbose) printf("Read eta as %s from %s.\n",eta_name,in_file);
//  read_eta(in_file,eta_name,&eta,&nlay);
  
  // Create the new NetCDF file.
  {
    if (verbose) printf("\tPrepare %s.\n",out_file);
    write_output(out_file,in_file,depth_file,
                 nz,var_names,time_names,depth_names,&nvar,arglist);
  }

  // Allocate space for the depth space variables.
  // if (verbose) printf("Read eta.\n");
  status = read_field(0,0,eta_name,&eta,&nint,&nx,&ny, st_eta);
  if (status != 0) {
    printf("ERROR: Unsuccessful in reading %s from %s.\n",eta_name,in_file[0]);
    exit(-1);
  }
  if (verbose) printf("\tEta %s starts at %g %g.\n",eta_name,st_eta[0],st_eta[1]);
  nxy = nx*ny;
  nlay = nint-1;

  if (verbose) printf("\tAllocating space for output variables.\n");
  var_out = (double *) calloc((size_t) (nz*nxy), sizeof(double));
  ev = (double *) calloc((size_t) nint, sizeof(double));
  elc = (double *) calloc((size_t) nint, sizeof(double));
  if (verbose) printf("\tDone with allocations.\n");
  

  for (lev=0;;lev++) {
    size_t junk, nlv, nxv, nyv, nxyv;
  
    if (verbose) printf("\tWorking on time level %d.\n",lev);
    
    if (lev>0) {
      // if (verbose) printf("Read %s for time %d.\n",eta_name,lev);
      status = read_field(0, lev, eta_name, &eta, &junk, &junk, &junk, st_eta);
      // if (verbose) printf("Done reading %s for time %d.\n",eta_name,lev);
    }
    if (status != 0) break;
    for (k=nint-1;k>=0;k--) for (i=0;i<nxy;i++) eta[k*nxy+i] -= eta[i];

    for (v=0;v<nvar;v++) {
            // This assumes that the missing value is -1e34.
      double missing_value = -1.0e34;
      int I, grid = 0;
      if ((lev>0) && static_var[v]) break;
     // if (verbose) printf("Read %s for time %d.\n",var_names[v],lev);
      status = read_field(var_file[v], lev, var_names[v], &var_in, &nlv, &nxv, &nyv, st_var);
      nxyv = nxv*nyv;
      if (status != 0) break;
     // if (verbose) printf("Done reading %s for time %d - starts at %g %g.\n",
     //   var_names[v],lev,st_var[0],st_var[1]);
      if (nlv != nlay) printf("Mismatch of layers %d vs expected %d.\n",(int) nlv, (int) nlay);

      if ((st_var[0] == st_eta[0]) && (st_var[1] == st_eta[1])) { grid = 0;
        if (verbose && (lev==0)) printf("\t%s is on the h-grid.\n",var_names[v]);
        if ((nyv != ny) || (nxv != nx)) {
          printf("h-grid variable %s is not same size as %s (%d x %d) vs (%d x %d).\n",
              var_names[v],eta_name,(int)nyv,(int)ny,(int)nxv,(int)nx);
          exit(-1);
        }
      }
      else if ((st_var[0] == st_eta[0]) && (st_var[1] != st_eta[1])) { grid = 1;
        if (verbose && (lev==0)) printf("\t%s is on the u-grid.\n",var_names[v]);
        if ((nyv != ny) || (abs((nxv-nx))>1)) {
          printf("u-grid variable %s is not consistent size with %s (%d x %d) vs (%d x %d).\n",
              var_names[v],eta_name,(int)nyv,(int)nxv,(int)ny,(int)nx);
          exit(-1);
        }
   //     printf("%s is on the u-grid.  No accomodation has yet been made for this grid.\n"
      }
      else if ((st_var[1] == st_eta[1]) && (st_var[0] != st_eta[0])) { grid = 2;
        if (verbose && (lev==0)) printf("\t%s is on the v-grid.\n",var_names[v]);
        if ((nxv != nx) || (abs((nyv-ny))>1)) {
          printf("v-grid variable %s is not consistent size with %s (%d x %d) vs (%d x %d).\n",
              var_names[v],eta_name,(int)nyv,(int)nxv,(int)ny,(int)nx);
          exit(-1);
        }
   //     printf("%s is on the v-grid.  No accomodation has yet been made for this grid.\n"
      } else {
        printf("%s is on the q-grid.  No accomodation has yet been made for this grid.\n"
               "Start locations: %g, %g vs %g, %g. Delta %g, %g\n",var_names[v],st_var[1],st_var[0],
                st_eta[1],st_eta[0],st_var[1]-st_eta[1],st_var[0]-st_eta[0]);
        exit(-1);
      }

/*      { int K;
        for (K=0;K<nz;K++) depths[K] = lay_depth[K]; 
      } */
      
      for (j=0;j<nyv;j++) for (i=0;i<nxv;i++) {
        I = i + nxv*j;
        // Interpolate interface height onto the variable's grid.
        if (grid == 0) {
          for (k=0;k<nint;k++) ev[k] = eta[k*nxy+I];
        } if (grid == 1) {
          if (st_var[1] < st_eta[1]) {
            if (i==0) {
              if (nx!=nxv) for (k=0;k<nint;k++) ev[k] = eta[k*nxy+j*nx];
              else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx]+eta[k*nxy+j*nx+nx-1]);
            }
            else if ((i==nxv-1)&&(nxv>nx)) {
              for (k=0;k<nint;k++) ev[k] = eta[k*nxy+j*nx+i-1];
            }
            else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+j*nx+i-1]);
          } else {
            if (i==nx-1) for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+j*nx]);
            else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+j*nx+i+1]);
          }
        } else if (grid == 2) {
          if (st_var[0] < st_eta[0]) {
            if (j==0) {
              if (ny!=nyv) for (k=0;k<nint;k++) ev[k] = eta[k*nxy+i];
              else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+(ny-1)*nx+i]);
            }
            else if ((j==nyv-1)&&(nyv>ny)) {
              for (k=0;k<nint;k++) ev[k] = eta[k*nxy+(j-1)*nx+i];
            }
            else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+(j-1)*nx+i]);
          } else {
            if (j==ny-1) for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+i]);
            else for (k=0;k<nint;k++) ev[k] = 0.5*(eta[k*nxy+j*nx+i]+eta[k*nxy+(j+1)*nx+i]);
          }
        }
        
        {
          int K, ktop=0, kbot=0;
          double z_sbl = 0.0, z_bbl = 0.0;
//          double e_bot;
          for (k=0;k<nlay;k++) elc[k] = 0.5*(ev[k]+ev[k+1]);
          for (k=0;k<nlay;k++) {
            if (((ev[0]-elc[k]) > 2.0*nlay*EPSILON) && (var_in[k*nxyv+I] != missing_value))
              {ktop = k; z_sbl = elc[k]; break;}
          }
          for (k=nlay-1;k>=0;k--) {
            if (((elc[k]-ev[nint-1]) > 2.0*nlay*EPSILON) && (var_in[k*nxyv+I] != missing_value))
              {kbot = k; z_bbl = elc[k]; break;}
          }

          // Fill in massless interior layers with values from above.
          for (k=ktop+1;k<kbot;k++) {
            if ((ev[k]-ev[k+1]) < 2.0*EPSILON) var_in[k*nxyv+I] = var_in[(k-1)*nxyv+I];
            if (var_in[k*nxyv+I] == missing_value) var_in[k*nxyv+I] = var_in[(k-1)*nxyv+I];
          }
          // Interpolate var into depth space.
          k=ktop;
//          if (z_bbl < 0.0)
//            e_bot = ev[nint-1];
          for (K=0;K<nz;K++) {
            if (lay_depth[K] > z_sbl)
              var_out[K*nxyv+I] = var_in[ktop*nxyv+I];
            else if (lay_depth[K] < ev[nint-1])
              var_out[K*nxyv+I] = missing;
            else if
              (lay_depth[K] < z_bbl) var_out[K*nxyv+I] = var_in[kbot*nxyv+I];
            else {
              while ((lay_depth[K] < elc[k]) && (k <= kbot)) k++;

              if ((lay_depth[K] >= elc[k]) && (lay_depth[K] <= elc[k-1])) {
                var_out[K*nxyv+I] = var_in[k*nxyv+I];
                /**/
                if (fabs(elc[k-1]-elc[k]) > EPSILON)
                  var_out[K*nxyv+I] += (var_in[(k-1)*nxyv+I] - var_in[k*nxyv+I]) *
                    ((lay_depth[K] - elc[k]) / (elc[k-1]-elc[k]));
                   /* */
              } else {
                printf("Unexpected error: k = %d, ktop = %d, kbot = %d, Depth %g is not between %g and %g.\n",
                        k,ktop,kbot,lay_depth[K],elc[k-1],elc[k]);
              }
            }

          }
        }
        
      }
      
      // if (verbose) printf("Done interpolating %s for time %d.\n",var_names[v],lev);
      write_field(lev,v,var_out);
      // if (verbose) printf("Done writing %s for time %d.\n",var_names[v],lev);
    }

    // if (verbose) printf("Start copying time at %d.\n",lev);
    copy_time_vals(lev);
    // if (verbose) printf("Done copying time at %d.\n",lev);
    
    // status = nc_sync(ncOutid);
    // if (status != NC_NOERR) handle_error(status,out_file,status);
  }

  nc_close(ncOutid);
  if (verbose) printf("Done: successfully created %s.\n",out_file);
  for (i=0;i<num_files;i++) if (ncInid[i] >= 0) nc_close(ncInid[i]);

}

void get_depths(char depth_file[], char depth_names[2][NC_MAX_NAME], double **lay_depth,
                double **int_depth, size_t *nz) {
  //  Read the depths of the intended layers from a NetCDF file.
  int ncid, status, layid, layvid, intvid, k;
  size_t nl;
  
  status = nc_open(depth_file,NC_NOWRITE, &ncid);
  if (status != NC_NOERR) handle_error(status,depth_file,status);

  status = nc_inq_dimid(ncid, depth_names[0], &layid);
  if (status != NC_NOERR) handle_error(status,depth_names[0],status);
  status = nc_inq_dimlen(ncid, layid, &nl); 
  if (status != NC_NOERR) handle_error(status,"Getting # Layers",0);
  *nz = nl;

  *lay_depth = (double *) calloc(nl, sizeof(double));
  *int_depth = (double *) calloc(nl+1, sizeof(double));
  if ((lay_depth == NULL) || (int_depth == NULL)) {
    printf("Unable to allocate input densities.\n"); exit(-1);
  }

  status = nc_inq_varid (ncid, depth_names[0], &layvid);
  if (status != NC_NOERR) handle_error(status,depth_names[0],0);
  status = nc_get_var_double(ncid, layvid, *lay_depth);
  if (status != NC_NOERR) handle_error(status,depth_names[0],0);

  status = nc_inq_varid (ncid, depth_names[1], &intvid);
  if (status != NC_NOERR) handle_error(status,depth_names[1],0);
  status = nc_get_var_double(ncid, intvid, *int_depth);
  if (status != NC_NOERR) handle_error(status,depth_names[1],0);

  nc_close(ncid);
  
  /* Check the sign convention and change to the HIM "height" convention. */
  if ((*int_depth)[0] < (*int_depth)[1]) {
    for (k=0;k<nl;k++) (*lay_depth)[k] *= -1.0;
    for (k=0;k<=nl;k++) (*int_depth)[k] *= -1.0;
  }
  
  /* Check for inversions in grid. */
  for (k=0;k<nl;k++) {
    if ((*lay_depth)[k] > (*int_depth)[k])
      printf("Inverted layer/interface %d/%d - %g / %g.\n",
          k,k,(*lay_depth)[k],(*int_depth)[k]);
    if ((*lay_depth)[k] < (*int_depth)[k+1])
      printf("Inverted layer/interface %d/%d - %g / %g.\n",
          k,k+1,(*lay_depth)[k],(*int_depth)[k+1]);
    if ((*int_depth)[k] < (*int_depth)[k+1])
      printf("Inverted interface %d/%d depths - %g / %g.\n",
          k,k+1,(*int_depth)[k],(*int_depth)[k+1]);
  }
  
}

int read_field(int fn, int lev, char var_name[], double **array,
               size_t *nlay, size_t *nx, size_t *ny, double st[2]) {
  static int ncid;
  size_t start[MAX_VAR_DIMS] = {0}, count[MAX_VAR_DIMS];
  int i, varid, ndims, dimids[MAX_VAR_DIMS], recdim, status;
  
  ncid = ncInid[fn];
  if (ncid < 0) {
    printf("ERROR: Required file (%d) for reading %s is not open.\n",fn,var_name);
    exit(-1);
  }

  // printf("Inquire unlimdim of %s - id %d.\n",in_file,ncid);
  status = nc_inq_unlimdim(ncid,&recdim);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  // printf("Done opening %s.\n",in_file);

  status = nc_inq_varid(ncid,var_name,&varid);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  status = nc_inq_varndims(ncid,varid,&ndims);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  status = nc_inq_vardimid(ncid,varid,dimids);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  // printf("Get %d dim lengths %s.\n",ndims,in_file);
  for (i=0;i<ndims;i++) {
    status = nc_inq_dimlen(ncid,dimids[i],&count[i]);
    if (status != NC_NOERR) handle_error(status,var_name,i);
  }
  if (dimids[0] == recdim) {
    if (count[0] < lev+1) return -1;
    start[0] = lev; count[0] = 1;
  } else if (lev > 0) return -1;
  
  for (i=ndims-2;i<ndims;i++) {
    char dimname[NC_MAX_NAME];
    int varid;
    size_t index = 0;
    status = nc_inq_dimname(ncid,dimids[i],dimname);
    if (status != NC_NOERR) handle_error(status,var_name,i);
    status = nc_inq_varid(ncid, dimname, &varid);
    if (status != NC_NOERR) handle_error(status,dimname,status);
    status = nc_get_var1_double(ncid, varid, &index, &st[i-ndims+2]);
    if (status != NC_NOERR) handle_error(status,dimname,status);
  }
  
  // printf("Done getting info about %s.\n",var_name);
  
  if (*array == NULL) {
    size_t sz = 1;
    for (i=0;i<ndims;i++) if (count[i] > 1) sz *= count[i]+1;
    // printf("Attempting to allocate array of size %ld.\n",(long) sz);
    *array = (double *) calloc(sz, sizeof(double));
    if (*array == NULL) printf("Unable to allocate array of size %ld for %s.\n",(long) sz,var_name);
  }
  
  *nlay = count[ndims-3];
  *nx = count[ndims-1];
  *ny = count[ndims-2];

  {
    size_t sz = 1;
    for (i=0;i<ndims;i++) if (count[i] > 1) sz *= count[i];
    // printf("Attempting to read array of size %ld.\n",(long) sz);
  }
  status = nc_get_vara_double(ncid,varid,start,count,*array);
  
  if (status != NC_NOERR) return -1;

  return 0;
}

int read_axes(char var_name[], double **array, size_t *nz, size_t *nx, size_t *ny) {
  size_t count[MAX_VAR_DIMS];
  int ncid, i, varid, ndims, dimids[MAX_VAR_DIMS], status;
  char dimname[NC_MAX_NAME];
  
  ncid = ncInid[0];
  if (ncid < 0) {
    printf("ERROR: Required file (%d) for reading axes of %s is not open.\n",0,var_name);
    exit(-1);
  }

  // printf("Inquire unlimdim of %s - id %d.\n",in_file,ncid);
  // status = nc_inq_unlimdim(ncid,&recdim);
  //if (status != NC_NOERR) handle_error(status,var_name,status);
  // printf("Done opening %s.\n",in_file);

  status = nc_inq_varid(ncid,var_name,&varid);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  status = nc_inq_varndims(ncid,varid,&ndims);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  if (ndims < 3) {
    printf("ERROR: Unable to determine array sizes from %d-D field %s in file 0.\n",
        ndims, var_name); return -1;
  }
  status = nc_inq_vardimid(ncid,varid,dimids);
  if (status != NC_NOERR) handle_error(status,var_name,status);
  // printf("Get %d dim lengths %s.\n",ndims,in_file);
  for (i=0;i<ndims;i++) {
    status = nc_inq_dimlen(ncid,dimids[i],&count[i]);
    if (status != NC_NOERR) handle_error(status,var_name,i);
  }
  *nz = count[ndims-3];
  *ny = count[ndims-2];
  *nx = count[ndims-1];

  // printf("Done getting info about %s.\n",var_name);

  if (*array == NULL) {
    *array = (double *) calloc(*nz, sizeof(double));
    if (*array == NULL) printf("Unable to allocate array of size %ld for %s.\n",(long) *nz,var_name);
  }
  
  status = nc_inq_dimname(ncid,dimids[ndims-3],dimname);
  if (status != NC_NOERR) handle_error(status,var_name,dimids[ndims-3]);
  status = nc_inq_varid(ncid, dimname, &varid);
  if (status != NC_NOERR) handle_error(status,dimname,status);
  status = nc_get_var_double(ncid,varid,*array);
  if (status != NC_NOERR) handle_error(status,dimname,status);
  
  if (status != NC_NOERR) return -1;

  return 0;
}


void write_field(int lev, int v, double *array) {
  int status;
  size_t start[MAX_VAR_DIMS] = {0};
  start[0] = lev;
  
  status = nc_put_vara_double(ncOutid,var_id[v],start,var_count[v],array);
  if (status != NC_NOERR) handle_error(status,"Writing field",v);
}

void copy_time_vals(size_t lev) {
  int status, i;
  double time_io;

  for (i=0;i<2;i++) if ((time_in_id[i] >= 0) && (time_out_id[i] >= 0)) {
    status = nc_get_var1_double(ncInid[0],time_in_id[i],&lev,&time_io);
    if (status != NC_NOERR) handle_error(status,"Reading time",lev);
    status = nc_put_var1_double(ncOutid,time_out_id[i],&lev,&time_io);
    if (status != NC_NOERR) handle_error(status,"Writing time",lev);
  }
}

void write_output(char out_file[], char in_file[][FILE_NAME_SZ], char depth_file[],
     size_t nz, char var_names[MAX_VARS][NC_MAX_NAME], char time_names[2][NC_MAX_NAME], 
     char depth_names[2][NC_MAX_NAME], int *num_var,
     char arglist[]) {
  int i, j, v;
  int ncid, nc_den_id, status, ndims, ndim, nvariables, ngatt, recdim;
  int use_dim[NC_MAX_DIMS];
  int dv_id[NC_MAX_DIMS], dim2id[NC_MAX_DIMS], dim2vid[NC_MAX_DIMS];
  int dimids[NC_MAX_VAR_DIMS], dd2id[2];
  int ddvid[2], dd2vid[2];
  int num_att;
  int use_file[MAX_FILES], max_file;
  char dim_name[10][NC_MAX_NAME];
  char att_name[NC_MAX_NAME];
  nc_type xtype;

  size_t dim_len[NC_MAX_DIMS], ddim_len[2], max_dim_len = 0;
    
  status = nc_open(in_file[0],NC_NOWRITE, &ncid);
  if (status != NC_NOERR) handle_error(status,in_file[0],status);
  ncInid[0] = ncid; use_file[0] = 1;

  status = nc_inq(ncid,&ndims,&nvariables,&ngatt,&recdim);
  if (status != NC_NOERR) handle_error(status,in_file[0],status);
  strcpy(master_in_file,in_file[0]);

  for (i=1;i<MAX_FILES;i++) { ncInid[i] = -1; use_file[i] = 0;}
  for (i=0;i<MAX_FILES;i++) if (strlen(in_file[i]) == 0) break;
  max_file = i;

/* Determine which dimensions need to be created. */
  for (i=0;i<ndims;i++) use_dim[i] = 0.0;
  for (v=0;v<*num_var;v++) {
    int vin_id, fn = 0, id;
    for (fn=0;fn<max_file;fn++) {
      if (ncInid[fn] < 0) {
        status = nc_open(in_file[fn],NC_NOWRITE, &ncInid[fn]);
        if (status != NC_NOERR) handle_error(status,in_file[fn],status);
      }
      status = nc_inq_varid(ncInid[fn], var_names[v], &vin_id);
      if (status == NC_NOERR) break;
    }
    if (fn==max_file) {
      printf("ERROR: Unable to find variable %s in any of %d files.\n",var_names[v],max_file);
      handle_error(status, var_names[v], v);
    }
    id = ncInid[fn];
    status = nc_inq_var(id, vin_id, att_name, &xtype, &ndim, dimids, &num_att);
    if (status != NC_NOERR) handle_error(status, var_names[v], v);

    if (ndim < 2) printf("Variable %s has only 2 dimensions and will be excluded.\n", var_names[v]);
    else {
      use_dim[find_dimid(id,dimids[ndim-1],ncid)] = 1;
      use_dim[find_dimid(id,dimids[ndim-2],ncid)] = 1;
      if (ndim > 4) {
        printf("ERROR: Variable %s has %d dimensions. This program only works with up to 4 dimensions.\n",
                var_names[v],ndim);
        exit(-1);
      }
      if (ndim == 4) {
        int frecdim; status = nc_inq_unlimdim(id,&frecdim);
        if (status != NC_NOERR) handle_error(status,"Finding record dimid",status);
        if (dimids[0] = frecdim) use_dim[recdim] = 1;
        else {
          printf("ERROR: Variable %s has 4 non-record dimensions. This program only works with 3 such dimensions.\n",
                  var_names[v]);
          exit(-1);
        }
      }
      
      var_file[v] = fn;
      use_file[fn] = 1;
    }
  }

  // Close any unneeded files.
  for (i=1;i<max_file;i++) if ((use_file[i] == 0) && (ncInid[i] >= 0)) {
    nc_close(ncInid[i]); ncInid[i] = -1;
  }
  
  status = nc_create(out_file, 0, &ncOutid);
  if (status != NC_NOERR) handle_error(status,out_file,status);
  status = nc_set_fill(ncOutid,NC_NOFILL,&j);
  if (status != NC_NOERR) handle_error(status,out_file,status);

  // printf("Created file %s with id %d.\n",out_file,ncOutid);
    
  // Copy all of the global attributes over.
  for (j=0;j<ngatt;j++) {
    status = nc_inq_attname (ncid, NC_GLOBAL, j, att_name);
    if (status != NC_NOERR) handle_error(status,"Global",j);    
    status = nc_copy_att(ncid, NC_GLOBAL, att_name, ncOutid, NC_GLOBAL);
    if (status != NC_NOERR) handle_error(status,att_name,j);
  }
  {
    char hist[1000];
    status = nc_get_att_text(ncid, NC_GLOBAL,"history",hist);
    if (status == NC_NOERR) {strcat(hist,"\n"); strcat(hist,arglist);}
    else strcpy(hist,arglist);
    status = nc_put_att_text(ncOutid, NC_GLOBAL,"history",strlen(hist),hist);
  }
  
  // Copy all appropriate dimensions over.
  for (i=0;i<ndims;i++) if (use_dim[i]) {
    status = nc_inq_dim(ncid, i, dim_name[i], &dim_len[i]);
    if (status != NC_NOERR) handle_error(status, "", i);
    if (dim_len[i] > max_dim_len) max_dim_len = dim_len[i];

    if (i==recdim)
      status = nc_def_dim(ncOutid, dim_name[i], NC_UNLIMITED, &dim2id[i]);
    else
      status = nc_def_dim(ncOutid, dim_name[i], dim_len[i], &dim2id[i]);
    if (status != NC_NOERR) handle_error(status,dim_name[i],i);

    // Get information about the coordinate variable.
    status = nc_inq_varid (ncid, dim_name[i], &dv_id[i]);
    if (status != NC_NOERR) handle_error(status, dim_name[i], i);
    status = nc_inq_vartype(ncid, dv_id[i], &xtype);
    if (status != NC_NOERR) handle_error(status, dim_name[i], i);
    status = nc_inq_varnatts(ncid, dv_id[i], &num_att);
    if (status != NC_NOERR) handle_error(status,dim_name[i],i);

    // Create the coordinate variables.
    status = nc_def_var (ncOutid, dim_name[i], xtype, 1, &dim2id[i], &dim2vid[i]);
    if (status != NC_NOERR) handle_error(status, dim_name[i],i);    
    // Copy all of the attributes over.
    for (j=0;j<num_att;j++) {
      status = nc_inq_attname (ncid, dv_id[i], j, att_name);
      if (status != NC_NOERR) handle_error(status,att_name,j);
      status = nc_copy_att(ncid, dv_id[i], att_name, ncOutid, dim2vid[i]);
      if (status != NC_NOERR) handle_error(status,att_name,j);
    }
  }
  
  // Copy the vertical dimensions over from depth_file.
  //    if (strlen(depth_file) > 1)
  status = nc_open(depth_file,NC_NOWRITE, &nc_den_id);
  if (status != NC_NOERR) handle_error(status,depth_file,status);

  for (i=0;i<2;i++) {
    int ddid;
    status = nc_inq_dimid (nc_den_id, depth_names[i], &ddid);
    if (status != NC_NOERR) handle_error(status,depth_names[i],0);
    status = nc_inq_dimlen(nc_den_id, ddid, &ddim_len[i]);
    if (status != NC_NOERR) handle_error(status,depth_names[i], i);

    status = nc_def_dim(ncOutid, depth_names[i], ddim_len[i], &dd2id[i]);
    if (status != NC_NOERR) handle_error(status,depth_names[i],i);

    // Get information about the coordinate variable.
    status = nc_inq_varid (nc_den_id, depth_names[i], &ddvid[i]);
    if (status != NC_NOERR) handle_error(status, depth_names[i], i);
    status = nc_inq_vartype(nc_den_id, ddvid[i], &xtype);
    if (status != NC_NOERR) handle_error(status, depth_names[i], i);
    status = nc_inq_varnatts(nc_den_id, ddvid[i], &num_att);
    if (status != NC_NOERR) handle_error(status,depth_names[i],i);

    // Create the coordinate variables.
    status = nc_def_var (ncOutid, depth_names[i], xtype, 1, &dd2id[i], &dd2vid[i]);
    if (status != NC_NOERR) handle_error(status, depth_names[i],i);    
    // Copy all of the attributes over.
    for (j=0;j<num_att;j++) {
      status = nc_inq_attname (nc_den_id, ddvid[i], j, att_name);
      if (status != NC_NOERR) handle_error(status,att_name,j);
      status = nc_copy_att(nc_den_id, ddvid[i], att_name, ncOutid, dd2vid[i]);
      if (status != NC_NOERR) handle_error(status,att_name,j);
    }
  }


  // Create the auxiliary time variable (if it exists) and store the time indices.
  if (recdim != -1) {
    // If there is a record dimension, it must be one of the two named
    // time dimensions.  If none are named, the name of the record
    // dimension is copied into it.
    if ((strcmp(dim_name[recdim],time_names[0]) != 0) &&
        (strcmp(dim_name[recdim],time_names[1]) != 0)) {
      if (strlen(time_names[0]) == 0) strcpy(time_names[0],dim_name[recdim]);
      else if (strlen(time_names[1]) == 0) strcpy(time_names[1],dim_name[recdim]);
      else {
        printf("ERROR: The specified time variables %s and %s do not agree\n"
               "\twith the record variable, %s.\n",time_names[0],time_names[1],
               dim_name[recdim]);
        exit(-1);
      }
    }
  }
  for (v=0;v<2;v++) {
    status = nc_inq_varid(ncOutid, time_names[v], &time_out_id[v]);
    if (status != NC_NOERR) {
      if (strlen(time_names[v]) > 0)  {
        int out_dim;
        status = nc_inq_varid(ncid, time_names[v], &time_in_id[v]);
        if (status != NC_NOERR) handle_error(status, time_names[v], v);
        status = nc_inq_var(ncid, time_in_id[v], att_name, &xtype, &ndim,
                            dimids, &num_att);
        if (status != NC_NOERR) handle_error(status, time_names[v], v);
        if (ndim > 1) {
          printf("ERROR: Time variable %s has %d dimensions in %s.\n",time_names[v],ndim,in_file[0]);
          exit(-1);
        }
        out_dim = dim2id[dimids[0]];

        status = nc_def_var(ncOutid, time_names[v], xtype, ndim, &out_dim, &time_out_id[v]);
        if (status != NC_NOERR) handle_error(status, time_names[v],v);
        // Copy all of the attributes over.
        for (j=0;j<num_att;j++) {
          status = nc_inq_attname(ncid, time_in_id[v], j, att_name);
          if (status != NC_NOERR) handle_error(status,att_name,j);
          status = nc_copy_att(ncid, time_in_id[v], att_name, ncOutid, time_out_id[v]);
          if (status != NC_NOERR) handle_error(status,att_name,j);
        }
      }
    } else {
      status = nc_inq_varid(ncid, time_names[v], &time_in_id[v]);
      if (status != NC_NOERR) handle_error(status, time_names[v], v);
    }
  }

  // Create the output variables, while checking the validity of the list.
  for (v=0;v<*num_var;v++) {
    int id, vin_id, frecdim, valid = 1;
    for (i=0;i<2;i++) if (strcmp(var_names[v],time_names[i])==0) valid = 0;
    for (i=0;i<ndims;i++) if (strcmp(var_names[v],dim_name[i])==0) valid = 0;
    for (i=0;i<2;i++) if (strcmp(var_names[v],depth_names[i])==0) valid = 0;
    
    if (valid) {
      id = ncInid[var_file[v]];
      
      if (var_file[v] == 0) frecdim = recdim;
      else {
        status = nc_inq_unlimdim(id,&frecdim);
        if (status != NC_NOERR) handle_error(status,"Finding record dimid",status);
      }

      status = nc_inq_varid(id, var_names[v], &vin_id);
      if (status != NC_NOERR) handle_error(status, var_names[v], v);
      status = nc_inq_var(id, vin_id, att_name, &xtype, &ndim,
                          dimids, &num_att);
      if (status != NC_NOERR) handle_error(status, var_names[v], v);

      if (ndim <= 2) {
        printf("Variable %s has only 2 dimensions and will be excluded.\n", var_names[v]);
        valid = 0;
      } else if ((ndim == 3) && (dimids[0] == frecdim)) {
        printf("Variable %s uses the record dimension as 3rd dimension and will be excluded.\n", var_names[v]);
        valid = 0;
      }
    }

    if (valid) {
      // Get information about the variable.
      int out_dimids[4];
      if (dimids[0] != frecdim) {
        out_dimids[0] = dd2id[0]; var_count[v][0] = ddim_len[0];
        static_var[v] = 1; i = 1;
      } else {
        out_dimids[0] = dim2id[recdim]; out_dimids[1] = dd2id[0];
        var_count[v][0] = 1; var_count[v][1] = ddim_len[0];
        static_var[v] = 0; i = 2;
      }
      var_size[v] = ddim_len[0];
      for (;i<ndim;i++) {
        int did;
        did = find_dimid(id,dimids[i],ncid);
        out_dimids[i] = dim2id[did];
        var_count[v][i] = dim_len[did];
        var_size[v] *= var_count[v][i];
      }

      status = nc_def_var(ncOutid, var_names[v], xtype, ndim, out_dimids, &var_id[v]);
      if (status != NC_NOERR) handle_error(status, var_names[v],v);
      // Copy all of the attributes over.
      for (j=0;j<num_att;j++) {
        status = nc_inq_attname(id, vin_id, j, att_name);
        if (status != NC_NOERR) handle_error(status,att_name,j);
        status = nc_copy_att(id, vin_id, att_name, ncOutid, var_id[v]);
        if (status != NC_NOERR) handle_error(status,att_name,j);
      }
      status = nc_get_att_double(id, vin_id,"missing_value",&missing_val[v]);
      if (status != NC_NOERR) {
        missing_val[v] = -1.0e34;
        status = nc_put_att_double(ncOutid,var_id[v],"missing_value",xtype,1,&missing_val[v]);
        if (status != NC_NOERR) handle_error(status,"missing_value",v);
      }
    } else {
      for (i=v;i<*num_var-1;i++) strcpy(var_names[i],var_names[i+1]);
      (*num_var)--; v--;
    }
  }

  status = nc_enddef(ncOutid);
  if (status != NC_NOERR) handle_error(status,out_file,status);
  // printf("Finished define mode for %s.\n",out_file);
    
    
    
    
/*
    // Create the vertical coordinates. //
    status = nc_def_dim(ncOutid, zc_name, nz, &layerid);
    if (status != NC_NOERR) handle_error(status,zc_name,0);
    status = nc_def_var (ncOutid, zc_name, NC_DOUBLE, 1, &layerid, &layervid);
    if (status != NC_NOERR) handle_error(status,zc_name,0);

    status = nc_def_dim(ncOutid, ze_name, (size_t) (nz+1), &intid);
    if (status != NC_NOERR) handle_error(status,ze_name,0);
    status = nc_def_var (ncOutid, ze_name, NC_DOUBLE, 1, &intid, &intvid);
    if (status != NC_NOERR) handle_error(status,ze_name,0);
    status = nc_put_att_text(ncOutid, intvid, "units", 2, "m");
    if (status != NC_NOERR) handle_error(status,"Units Interface",0);

    dims[0] = layervid;
    {
      strcpy(att[0][0],"long_name"); strcpy(att[0][1],"Depth of Layer Center");
//      strcpy(att[1][0],"units"); strcpy(att[1][1],"m");
      strcpy(att[1][0],"units"); strcpy(att[1][1],"cm");
      strcpy(att[2][0],"positive"); strcpy(att[2][1],"down");
      strcpy(att[3][0],"edges"); strcpy(att[2][1],ze_name);
      for (j=0;j<=2;j++) {
        status = nc_put_att_text(ncOutid, layervid, att[j][0],
            strlen(att[j][1]), att[j][1]);
        if (status != NC_NOERR) handle_error(status,att[j][0],j);
      }
    
      strcpy(att[0][0],"long_name"); strcpy(att[0][1],"Depth of edges");
//      strcpy(att[1][0],"units"); strcpy(att[1][1],"m");
      strcpy(att[1][0],"units"); strcpy(att[1][1],"cm");
      strcpy(att[2][0],"positive"); strcpy(att[2][1],"down");
      for (j=0;j<=2;j++) {
        status = nc_put_att_text(ncOutid, intvid, att[j][0],
            strlen(att[j][1]), att[j][1]);
        if (status != NC_NOERR) handle_error(status,att[j][0],j);    
      }
    }
*/
     
    // Copy the axis values from the first file.
  {
    double *coord_array;
    coord_array = malloc(sizeof(double)*(max_dim_len));
    for (i=0;i<ndims;i++) if (use_dim[i] && (i!=recdim)) {
      status = nc_get_var_double(ncid,dv_id[i],coord_array);
      if (status != NC_NOERR) handle_error(status,"Read Coordinate Variable",i);
      
      status = nc_put_var_double(ncOutid,dim2vid[i],coord_array);
      if (status != NC_NOERR) handle_error(status,"Write Coordinate Variable",i);
    }
    
    // Copy the vertical axis values from the depth file.
    for (i=0;i<2;i++) {
      status = nc_get_var_double(nc_den_id,ddvid[i],coord_array);
      if (status != NC_NOERR) handle_error(status,"Read Coordinate Variable",i);
      
      status = nc_put_var_double(ncOutid,dd2vid[i],coord_array);
      if (status != NC_NOERR) handle_error(status,"Write Coordinate Variable",i);
    }

    free(coord_array);
  }

    /*
    {
      double *z;
      z = (double *) calloc((size_t) (nz+1), sizeof(double));
      for (i=0;i<=nz;i++) z[i] = (-100.0)*int_depth[i];
      status = nc_put_var_double(ncOutid,intid,z);
      if (status != NC_NOERR) handle_error(status,"Interface Coordinate",0);

      for (i=0;i<nz;i++) z[i] = (-100.0)*lay_depth[i];
      status = nc_put_var_double(ncOutid,layerid,z);
      if (status != NC_NOERR) handle_error(status,"Layer Coordinate",0);
    }
    */

  nc_close(nc_den_id);
}

int find_dimid(int ncid1,int dimf1,int ncid2) {
  // Find the dimension ID in ncid2 with the same name as dimf1 in ncid1.
  char dimname[NC_MAX_NAME], message[200];
  int status, dimid;
  size_t len1, len2;

  if (ncid1==ncid2) return dimf1;
  
  status = nc_inq_dim(ncid1,dimf1,dimname,&len1);
  if (status != NC_NOERR) handle_error(status,"Finding dimension name",status);
  status = nc_inq_dimid(ncid2,dimname,&dimid);
  if (status != NC_NOERR) {
    sprintf(message,"Cannot find dimension %s in first input file.\n",dimname);
    handle_error(status,dimname,status);
  }
  status = nc_inq_dimlen(ncid2,dimid,&len2);

  if (len1 != len2) {
    printf("ERROR: Dimension %s does not agree in size between files; it is %ld and %ld.\n",
           dimname,(long)len2,(long)len1);
    exit(-1);
  }

  // Compare axis data?

  return dimid;
}

void handle_error(int status, char msg[], int id) {
  if (status != NC_NOERR) {
    fprintf(stderr, "%s %s %d\n", nc_strerror(status), msg, id);
    exit(-1);
  }
}
