#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#define  TILE_INDEX_NAME "tile_index"
#define  COHORT_INDEX_NAME "cohort_index"
#define  COHORT_NAME "cohort"
#define  TILE_NAME "tile"
#define  LON_NAME        "lon"
#define  LAT_NAME        "lat"
#define  TIMENAME        "time"
#define  D2R (M_PI/180.)

char *usage[] = {
  "",
  "                                                                                ",
  "                                                                                ",
  "                    Usage of remap_land                                         ",
  "                                                                                ",
  "   remap_land --src_mosaic src_mosaic --src_restart src_restart                 ",
  "              --dst_mosaic dst_msoaic --dst_restart dst_restart                 ",
  "              --dst_cold_restart dst_cold_restart [--remap_file remap_file]     ",             
  "                                                                                ",
	  " remap_land remap land restart file from one mosaic grid to another mosaic grid ",
  " remap_land takes the following flags,                                          ",
  "                                                                                ",
  " REQUIRED:                                                                      ",
  "                                                                                ",
  " --src_mosaic src_mosaic    specify the source mosaic information. This file    ",
  "                            contains list of tile files which specify the grid  ",
  "                            information for each tile.                          ",
  "                                                                                ",
  " --dst_mosaic dst_mosaic    specify the destination mosaic information. This    ",
  "                            file contains list of tile files which specify the  ",
  "                            grid information for each tile.                     ",
  "                                                                                ",
  " --src_restart src_restart  specify the source restart file.                    ",
  "                                                                                ",
  " --dst_restart dst_restart  specify the restart file to be generated on         ",
  "                            destination grid.                                   ",
  "                                                                                ",
  " --dst_cold_restart file    specify the cold restart file destination grid.     ",
  "                            This is the input file. The dst_cold_restart_file   ",
  "                            could be obtained by running the experiment         ",
  "                            for 1 day with --dst_mosaic using cold restart.     ",
  "                                                                                ",
  "                                                                                ",   
  " OPTIONAL FLAGS                                                                 ",
  "                                                                                ",
  "                                                                                ", 
  " --remap_file remap_file    specify the file name that saves remapping          ",
  "                            information. If remap_file is specified and the     ",
  "                            file does not exist, remapping information will be  ",
  "                            calculated and stored in remap_file. If remap_file  ",
  "                            is specified and the file exists, remapping         ",
  "                            information will be read from remap_file.           ",
  "                                                                                ",
  "                                                                                ",
  " Example: remap land restart from C48 grid onto C180 grid.                      ",
  "                                                                                ",
  "   remap_land --src_mosaic /archive/z1l/tools/input/C48/C48_mosaic.nc           ",
  "              --dst_mosaic /archive/z1l/tools/input//C180/C180_mosaic.nc        ",
  "              --src_restart /archive/z1l/tools/input/land_src_restart/soil.res  ",
  "              --dst_cold_restart /archive/z1l/tools/input/land_dst_cold_restart/soil.res ",
  "              --dst_restart soil.res                                            ",
  "                                                                                ",
  NULL };

char grid_version[] = "0.2";
char tagname[] = "$Name: siena_201205_z1l $";


void get_actual_file_name(int nface, int face, const char *file_orig, char *file);
void full_search_nearest(int nface_src, const int *npts_src, const double *lon_src, const double *lat_src,
			 const int *tile_type_src, int npts_dst, const double *lon_dst,
			 const double *lat_dst, const int *tile_type_dst, int *idx_map, int *face_map);


int main(int argc, char* argv[])
{

  char *src_mosaic            = NULL;
  char *dst_mosaic            = NULL;
  char *src_restart_file      = NULL;
  char *dst_restart_file      = NULL;
  char *dst_cold_restart      = NULL;
  char *remap_file            = NULL;
  
  int nface_src = 0, nface_dst = 0;
  int nx_src, ny_src, nx_dst, ny_dst;
  int nidx_tot_src, n_cohort, n_tile;
  int *fid_src=NULL;
  int *cohort_data=NULL, *tile_data=NULL;
  int *nidx_src=NULL, *start_pos=NULL;
  int *idx_src=NULL, *idx_type_src=NULL;
  
  char   history[1280];
  double *x_src=NULL, *y_src=NULL;
  int time_exist, ntime, l;
  int *has_taxis=NULL, *var_type=NULL, *ndim_src=NULL, *nz_src=NULL;
  double *time_data=NULL;
  
  int    npes, nvar_src;
  int    option_index, c;
  int    errflg = (argc == 1);  
  /*
   * process command line
   */

  static struct option long_options[] = {
    {"src_mosaic",       required_argument, NULL, 'a'},
    {"dst_mosaic",       required_argument, NULL, 'b'},
    {"src_restart",      required_argument, NULL, 'c'},
    {"dst_restart",      required_argument, NULL, 'o'},
    {"dst_cold_restart", required_argument, NULL, 'd'},
    {"remap_file",       required_argument, NULL, 'r'},
    {"help",             no_argument,       NULL, 'h'},
    {0, 0, 0, 0},
  };

  /* start parallel */

  mpp_init(&argc, &argv);

  mpp_domain_init();  
   
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
    switch (c) {
    case 'a':
      src_mosaic = optarg;
      break;      
    case 'b':
      dst_mosaic = optarg;
      break;      
    case 'c':
      src_restart_file = optarg;
      break;      
    case 'o':
      dst_restart_file = optarg;
      break;
    case 'd':
      dst_cold_restart = optarg;
      break;      
    case 'r':
      remap_file = optarg;
      break;
    case '?':
      errflg++;
      break;
    }

  if( !src_mosaic ) errflg++;
  if( !dst_mosaic) errflg++;
  if( !src_restart_file) errflg++;
  if( !dst_restart_file ) errflg++;
  if( !dst_cold_restart ) errflg++;

  if (errflg) {
    char **u = usage;
    if(mpp_pe() == mpp_root_pe()) {
      while (*u) { fprintf(stderr, "%s\n", *u); u++; }
      if(!src_mosaic) mpp_error("remap_land: src_mosaic is not specified");
      if(!dst_mosaic) mpp_error("remap_land: dst_mosaic is not specified");
      if(!src_restart_file) mpp_error("remap_land: src_restart_file is not specified");
      if(!dst_restart_file) mpp_error("remap_land: dst_restart_file is not specified");
      if(!dst_cold_restart) mpp_error("remap_land: dst_cold_restart is not specified");
    }
    mpp_error("remap_land: check the command line arguments");
  }

  npes = mpp_npes();
  
  /*----------------------------------------------------------------------------
    get source grid size
    --------------------------------------------------------------------------*/
  {
    int *nx, *ny;
    int n;
    nface_src = read_mosaic_ntiles(src_mosaic);
    nface_dst = read_mosaic_ntiles(dst_mosaic);
    nx        = (int *)malloc(nface_src*sizeof(int));
    ny        = (int *)malloc(nface_src*sizeof(int));
    read_mosaic_grid_sizes(src_mosaic, nx, ny);
    /* nx, ny of source should have the same value on each face */
    for(n=1; n<nface_src; n++) {
      if(nx[n] != nx[0] || ny[n] != ny[0]) mpp_error("remap_land: all the faces of source grid should have the same number of grid points");
    }
    nx_src = nx[0];
    ny_src = ny[0];
    free(nx);
    free(ny);
    nx        = (int *)malloc(nface_dst*sizeof(int));
    ny        = (int *)malloc(nface_dst*sizeof(int));
    read_mosaic_grid_sizes(dst_mosaic, nx, ny);
    /* nx, ny of source should have the same value on each face */
    for(n=1; n<nface_dst; n++) {
      if(nx[n] != nx[0] || ny[n] != ny[0]) mpp_error("remap_land: all the faces of destination grid should have the same number of grid points");
    }
    nx_dst = nx[0];
    ny_dst = ny[0];
    free(nx);
    free(ny);    
  }
    
  /*-----------------------------------------------------------------------------
    read the  tile_index on source grid and source grid.
    ----------------------------------------------------------------------------*/
  {
    double *x_tmp=NULL, *y_tmp=NULL;
    int    n, pos;
    
    fid_src      = (int *)malloc(nface_src*sizeof(int));
    nidx_src     = (int *)malloc(nface_src*sizeof(int));
    nidx_tot_src = 0;
    for(n=0; n<nface_src; n++) {
      int nlon, nlat;
      char file[512];
      get_actual_file_name(nface_src, n, src_restart_file, file);
      fid_src[n] = mpp_open(file, MPP_READ);
      nlon = mpp_get_dimlen(fid_src[n], LON_NAME);
      nlat = mpp_get_dimlen(fid_src[n], LAT_NAME);
      if(nx_src != nlon) mpp_error("remap_land: mismatch on the longitude dimension size "
				   "between source mosaic grid file and src_restart");
      if(ny_src != nlat) mpp_error("remap_land: mismatch on the latitude dimension size "
				   "between source mosaic grid file and src_restart");    
      nidx_src[n] = mpp_get_dimlen(fid_src[n], TILE_INDEX_NAME);
      nidx_tot_src += nidx_src[n];
    }
  
    idx_src      = (int *)malloc(nidx_tot_src*sizeof(int));
    idx_type_src = (int *)malloc(nidx_tot_src*sizeof(int));
    x_src        = (double *)malloc(nidx_tot_src*sizeof(double));
    y_src        = (double *)malloc(nidx_tot_src*sizeof(double));
    for(n=0; n<nidx_tot_src; n++) {
      idx_type_src[n] = -1;
    }
    pos = 0;

    x_tmp = (double *)malloc(nx_src*ny_src*sizeof(double));
    y_tmp = (double *)malloc(nx_src*ny_src*sizeof(double));
    for(n=0; n<nface_src; n++) {
      int vid_src, l, i, j, k;
    
      read_mosaic_grid_data(src_mosaic, "x", nx_src, ny_src, x_tmp, n, 1, 1);
      read_mosaic_grid_data(src_mosaic, "y", nx_src, ny_src, y_tmp, n, 1, 1);
    
      vid_src = mpp_get_varid(fid_src[n], TILE_INDEX_NAME);
      mpp_get_var_value(fid_src[n], vid_src, idx_src+pos);
      for(l=pos; l<pos+nidx_src[n]; l++) {
	if(idx_src[l] < 0) continue;
	i = idx_src[l]%nx_src;
	k = idx_src[l]/nx_src;
	j = k%ny_src;
	k = k/ny_src;
	idx_type_src[l] = k;
	x_src[l]        = x_tmp[j*nx_src+i]*D2R;
	y_src[l]        = y_tmp[j*nx_src+i]*D2R;
      }
      pos += nidx_src[n];
    }
    free(x_tmp);
    free(y_tmp);
  }

  /*-----------------------------------------------------------------------------
    Get the time information
    ---------------------------------------------------------------------------*/
  time_exist = mpp_dim_exist(fid_src[0], TIMENAME);
  ntime = 1;
  if(time_exist) {
    int vid;
    
    ntime = mpp_get_dimlen(fid_src[0], TIMENAME);
    time_data = (double *)malloc(ntime*sizeof(double));
    vid = mpp_get_varid(fid_src[0], TIMENAME);
    mpp_get_var_value(fid_src[0], vid, time_data);
  }
  
  /* loop through each variable of the source data to see get variable information*/  
  nvar_src = mpp_get_nvars(fid_src[0]);
  has_taxis    = (int *)malloc(nvar_src*sizeof(int));
  ndim_src     = (int *)malloc(nvar_src*sizeof(int));
  nz_src       = (int *)malloc(nvar_src*sizeof(int));
  var_type     = (int *)malloc(nvar_src*sizeof(int));

  for(l=0; l<nvar_src; l++) {
    char varname[256];
    int vid, m;
    
    has_taxis[l] = 0;
    nz_src[l] =1;
    mpp_get_varname(fid_src[0], l, varname);
    vid = mpp_get_varid(fid_src[0], varname);
    var_type[l] = mpp_get_var_type(fid_src[0], vid);
    if( var_type[l] != MPP_INT && var_type[l] != MPP_DOUBLE)
	mpp_error("remap_land: field type must be MPP_INT or MPP_DOUBLE");
    
    ndim_src[l] = mpp_get_var_ndim(fid_src[0], vid);
    if(ndim_src[l] > 2) mpp_error("remap_land: number of dimensions for the field in src_restart is greater than 2");
    for(m=0; m<ndim_src[l]; m++) {
      mpp_get_var_dimname(fid_src[0], vid, m, varname);
      if( !strcmp(varname, TIMENAME) ) has_taxis[l] = 1;
    }
    if(ndim_src[l] == 2 && !has_taxis[l] ) {
      mpp_get_var_dimname(fid_src[0], vid, 1, varname);
      nz_src[l] = mpp_get_dimlen(fid_src[0], varname);
    }
  }

  /* get the cohort and tile data when time axis exist*/
  n_cohort = 0;
  n_tile = 0;
  if(time_exist) {
    int dimsize, m, vid, i;
    int *tmp;
    if(!mpp_dim_exist(fid_src[0], COHORT_NAME)) mpp_error("remap_land: dimension cohort should exist when time exist");
    if(!mpp_var_exist(fid_src[0], COHORT_NAME)) mpp_error("remap_land: field cohort should exist when time exist");
    n_cohort = mpp_get_dimlen(fid_src[0], COHORT_NAME);
    cohort_data = (int *)malloc(n_cohort*sizeof(int));
    vid = mpp_get_varid(fid_src[0], COHORT_NAME);
    mpp_get_var_value(fid_src[0], vid, cohort_data);
    for(m=1; m<nface_src; m++) {
      dimsize = mpp_get_dimlen(fid_src[m], COHORT_NAME);
      if(dimsize != n_cohort)mpp_error("remap_land: the dimension size of cohort is different between faces");
      tmp = (int *)malloc(n_cohort*sizeof(int));
      vid = mpp_get_varid(fid_src[m], COHORT_NAME);
      mpp_get_var_value(fid_src[m], vid, tmp);
      for(i=0; i<n_cohort; i++) {
	if(cohort_data[i] != tmp[i]) mpp_error("remap_land: cohort value is different between faces");
      }
      free(tmp);
    }
    
    if(!mpp_dim_exist(fid_src[0], TILE_NAME)) mpp_error("remap_land: dimension tile should exist when time exist");
    if(!mpp_var_exist(fid_src[0], TILE_NAME)) mpp_error("remap_land: field tile should exist when time exist");

    n_tile   = mpp_get_dimlen(fid_src[0], TILE_NAME);
    tile_data   = (int *)malloc(n_tile*sizeof(int));
    vid = mpp_get_varid(fid_src[0], TILE_NAME);
    mpp_get_var_value(fid_src[0], vid, tile_data);  
    for(m=1; m<nface_src; m++) {
      dimsize = mpp_get_dimlen(fid_src[m], TILE_NAME);
      if(dimsize != n_tile)mpp_error("remap_land: the dimension size of tile is different between faces");
      tmp = (int *)malloc(n_tile*sizeof(int));
      vid = mpp_get_varid(fid_src[m], TILE_NAME);
      mpp_get_var_value(fid_src[m], vid, tmp);
      for(i=0; i<n_tile; i++) {
	if(tile_data[i] != tmp[i]) mpp_error("remap_land: tile value is different between faces");
      }
      free(tmp);      
    }
  }

  {
    int n;
    strcpy(history,argv[0]);
    for(n=1;n<argc;n++) {
      strcat(history, " ");
      strcat(history, argv[n]);
    }
  }
  /*------------------------------------------------------------------------------------------
    loop through each face of destination grid, first read the grid, then read the tile_index,
    then find the remapping index, then setup metadata for the destination file,
    last do the remapping and write out the data to dst_restart_file
    ----------------------------------------------------------------------------------------*/
  {
    double *x_tmp, *y_tmp;
    int    n, pos;
    char   cmd[1024];
    start_pos = (int *)malloc(nface_src*sizeof(int));
  
    x_tmp = (double *) malloc(nx_dst*ny_dst*sizeof(double));
    y_tmp = (double *) malloc(nx_dst*ny_dst*sizeof(double));

    for(n=0; n<nface_dst; n++) {
      double *x_dst, *y_dst;
      int *idx_dst, *idx_type_dst, *idx_tmp;
      int *idx_map, *face_map;
      int nlon, nlat, layout[2], l, ll, i, j, k, t;
      int fid_dst, vid_dst, vid_src, fid_read, vid_read;
      int nidx_dst, nvar_dst;
      int isc_dst, iec_dst, jsc_dst, jec_dst, nxc_dst, nyc_dst;
      size_t start[4], nread[4], nwrite[4];
      char file_dst[512], file_cold[512];
      domain2D Dom_dst;
    
      get_actual_file_name(nface_dst, n, dst_restart_file, file_dst);
      get_actual_file_name(nface_dst, n, dst_cold_restart, file_cold);
    
      fid_read = mpp_open(file_cold, MPP_READ);
      nidx_dst = mpp_get_dimlen(fid_read, TILE_INDEX_NAME);

      /* when the time exist ( static vegetation), need to create a new file. otherwise copy
	 dst_cold_restart to dst_restart_file
      */
      if(time_exist) {
	fid_dst = mpp_open(file_dst, MPP_WRITE);
      }
      else {
	if(mpp_pe() == mpp_root_pe() ) {
	  sprintf(cmd, "cp -f %s %s", file_cold, file_dst);
	  system(cmd);
	}
        /* when ndix_dst == 1, no need to do remapping, simply return here. */
        if( nidx_dst == 1) continue;
	fid_dst = mpp_open(file_dst, MPP_APPEND);
      }

      nlon = mpp_get_dimlen(fid_read, LON_NAME);
      nlat = mpp_get_dimlen(fid_read, LAT_NAME);
      if(nx_dst != nlon) mpp_error("remap_land: mismatch on the longitude dimension size "
				   "between destination mosaic grid file and dst_cold_restart");
      if(ny_dst != nlat) mpp_error("remap_land: mismatch on the latitude dimension size "
				   "between destination mosaic grid file and dst_cold_restart");
    
      /* setup domain */
      mpp_define_layout( nidx_dst, 1, npes, layout);
      mpp_define_domain2d( nidx_dst, 1, layout, 0, 0, &Dom_dst);
      mpp_get_compute_domain2d( Dom_dst, &isc_dst, &iec_dst, &jsc_dst, &jec_dst);
      /* for this domain, jsc_dst should equal jec_dst */
      if(jsc_dst != jec_dst) {
	mpp_error("remap_land: This is a 1-D domain decomposition, jsc_dst must equal to jec_dst");
      }
      nxc_dst = iec_dst - isc_dst + 1;
      nyc_dst  = 1;
    
      idx_tmp      = (int *) malloc(nidx_dst*sizeof(int));
      idx_dst      = (int *) malloc(nxc_dst*sizeof(int));
      idx_type_dst = (int *) malloc(nxc_dst*sizeof(int));
      idx_map      = (int *) malloc(nxc_dst*sizeof(int));
      face_map     = (int *) malloc(nxc_dst*sizeof(int));
      x_dst        = (double *) malloc(nxc_dst*sizeof(double));
      y_dst        = (double *) malloc(nxc_dst*sizeof(double));
      read_mosaic_grid_data(dst_mosaic, "x", nx_dst, ny_dst, x_tmp, n, 1, 1);
      read_mosaic_grid_data(dst_mosaic, "y", nx_dst, ny_dst, y_tmp, n, 1, 1);

      vid_read = mpp_get_varid(fid_read, TILE_INDEX_NAME);
      mpp_get_var_value(fid_read, vid_read, idx_tmp);
      for(l=0; l<nxc_dst; l++) {
	idx_type_dst[l] = -1;
      }
    
      for(l=isc_dst; l<=iec_dst; l++) {
	ll = l-isc_dst;
	idx_dst[ll] = idx_tmp[l];
	if(idx_tmp[l] < 0) continue;
	idx_dst[ll] = idx_tmp[l];
	i = idx_dst[ll]%nx_dst;
	k = idx_dst[ll]/nx_dst;
	j = k%ny_dst;
	k = k/ny_dst;
	idx_type_dst[ll] = k;
	x_dst[ll]        = x_tmp[j*nx_dst+i]*D2R;
	y_dst[ll]        = y_tmp[j*nx_dst+i]*D2R;      
      }

      /* define the metadata for dst_restart_file */
      if(time_exist) {
	int dim_time, dim_cohort_index, dim_lat, dim_lon;
	int dim_tile_index, dim_cohort, dim_tile;
      
	dim_time         = mpp_def_dim(fid_dst, TIMENAME, NC_UNLIMITED);
	dim_cohort_index = mpp_def_dim(fid_dst, COHORT_INDEX_NAME, nidx_dst);
	dim_lat          = mpp_def_dim(fid_dst, LAT_NAME, ny_dst);
	dim_lon          = mpp_def_dim(fid_dst, LON_NAME, nx_dst);
	dim_tile_index   = mpp_def_dim(fid_dst, TILE_INDEX_NAME, nidx_dst);
	dim_cohort       = mpp_def_dim(fid_dst, COHORT_NAME, n_cohort);
	dim_tile         = mpp_def_dim(fid_dst, TILE_NAME, n_tile);

	for(l=0; l<nvar_src; l++) {
	  char varname[256], dimname[256];
	  int vid1, vid2, ndim, m, dims[4];
	
	  mpp_get_varname(fid_src[0], l, varname);
	  vid1 = mpp_get_varid(fid_src[0], varname);
	  ndim = mpp_get_var_ndim(fid_src[0], vid1);

	  for(m=0; m<ndim; m++) {
	    mpp_get_var_dimname(fid_src[0], vid1, m, dimname);
	    if( !strcmp(dimname, TIMENAME) )
	      dims[m] = dim_time;
	    else if( !strcmp(dimname, COHORT_INDEX_NAME) )
	      dims[m] = dim_cohort_index;
	    else if( !strcmp(dimname, LAT_NAME ) )
	      dims[m] = dim_lat;
	    else if( !strcmp(dimname, LON_NAME ) )
	      dims[m] = dim_lon;
	    else if( !strcmp(dimname, TILE_INDEX_NAME ))
	      dims[m] = dim_tile_index;
	    else if( !strcmp(dimname, COHORT_NAME) )
	      dims[m] = dim_cohort;
	    else if( !strcmp(dimname, TILE_NAME ))
	      dims[m] = dim_tile;
	    else
	      mpp_error("REMAP_LAND: invalid dimension name");
	  }
	  vid2 = mpp_def_var(fid_dst, varname, var_type[l], ndim, dims, 0);
	  mpp_copy_var_att(fid_src[0], vid1, fid_dst, vid2);  
	}
      }
      else {
	mpp_redef(fid_dst);
      }      
      mpp_def_global_att(fid_dst, "history", history);
      mpp_end_def(fid_dst);
    
    /*--------------------------------------------------------------------------------
      if remap_file exists, read the remap file, otherwise
      Find the remap index
      -------------------------------------------------------------------------------*/
      {
	int remap_file_exist;
	int write_remap_file;
	char file[512];
	remap_file_exist = 0;
	write_remap_file = 0;
	if(remap_file) {
	  get_actual_file_name(nface_dst, n, remap_file, file);
	  remap_file_exist = mpp_file_exist(file);
	  if( !remap_file_exist) write_remap_file = 1;
	}
	if(remap_file_exist) { /* read from remap file */
	  size_t start[4], nread[4];
	  int fid, nidx, vid;

	  if(mpp_pe() == mpp_root_pe()) printf("Read remap information from remap_file \n");
	  for(l=0; l<4; l++) {
	    start[l] = 0; nread[l] =1;
	  }
	  start[0] =isc_dst; nread[0] = nxc_dst;
	  fid = mpp_open(file, MPP_READ);
	  nidx = mpp_get_dimlen(fid, TILE_INDEX_NAME);
	  if(nidx != nidx_dst) mpp_error("remap_land: mismatch of dimension length of tile_index between read_restart and remap_file");
	
	  vid = mpp_get_varid(fid, "remap_index");
	  mpp_get_var_value_block(fid, vid, start, nread, idx_map);
	  vid = mpp_get_varid(fid, "remap_face");
	  mpp_get_var_value_block(fid, vid, start, nread, face_map);
	
	  mpp_close(fid);
	}
	else {
	  full_search_nearest(nface_src, nidx_src, x_src, y_src, idx_type_src, nxc_dst,
			      x_dst, y_dst, idx_type_dst, idx_map, face_map);
	}

	if(write_remap_file) { /* write out restart file */
	  int *gdata;
	  int fid, dim_tile_index, id_remap_index, id_remap_face;
	
	  fid = mpp_open(file, MPP_WRITE);
	  dim_tile_index = mpp_def_dim(fid, TILE_INDEX_NAME, nidx_dst);
	  id_remap_index = mpp_def_var(fid, "remap_index", NC_INT, 1, &dim_tile_index, 1,
				       "standard_name", "remap index");
	  id_remap_face  = mpp_def_var(fid, "remap_face", NC_INT, 1, &dim_tile_index, 1,
				       "standard_name", "remap face");
	  mpp_def_global_att(fid, "history", history);
	  mpp_end_def(fid);
	  gdata = (int *)malloc(nidx_dst*sizeof(int));
	  mpp_gather_field_int_root(nxc_dst, idx_map, gdata);
	  mpp_put_var_value(fid, id_remap_index, gdata);
	  mpp_gather_field_int_root(nxc_dst, face_map, gdata);
	  mpp_put_var_value(fid, id_remap_face, gdata);	
	  mpp_close(fid);
	  free(gdata);
	}
      }
      /*-------------------------------------------------------------------------------
	Remap the data and write out to dst_restart_file
	It is assumed all the fields need to be remapped will have the first dimension
	name "tile_index" or "cohort_index" ( excludes "tile_index" or "cohort_index")
	-----------------------------------------------------------------------------*/

      /* First find number of varaibles in each file */

      nvar_dst = mpp_get_nvars(fid_read);

      if( !time_exist && nvar_dst != nvar_src)
	mpp_error("remap_land: nvar is different in src_restart and dst_read_restart when time does not exist");

      /* loop through each time level */
      for(t=0; t<ntime; t++) {
	for(l=0; l<nvar_src; l++) {
	  double *data_src, *data_dst, *gdata_dst;
	  int    *idata_src, *idata_dst;
	  int    *start_pos;
	  char   varname[128], dimname[128];
	  int    nz_dst;
	  int    m;

	  if( !has_taxis[l] && t>0 ) continue;      
	  mpp_get_varname(fid_src[0], l, varname);
	  vid_dst = mpp_get_varid(fid_dst, varname);
	  vid_src = mpp_get_varid(fid_src[0], varname);

	  if(time_exist) {
	    if(strcmp(varname, TIMENAME) == 0) {
	      /* copy the time data from src_restart_file to dst_restart_file */
	      for(m=0; m<4; m++) {
		start[m]=0; nwrite[m] = 1;
	      }
	      start[0] = t;
	      mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, time_data+t);
	      continue;
	    }
	  }
          else { /* when time axis does not exist, ndim should be the same for src_restart_file and dst_cold_restart */
	    int ndim_dst;
	    vid_read = mpp_get_varid(fid_read, varname);
	    ndim_dst = mpp_get_var_ndim(fid_read, vid_read);
	    if( ndim_src[l] != ndim_dst)
	      mpp_error("remap_land: number of dimensions for the field in dst_read_restart "
			"does not match that in src_restart_file when time does not exist");
	    if(ndim_dst == 0) continue;	
	  }

	  /* when time axis exist, need to copy tile_index, cohort_index, lon and lat from dst_cold_restart
	     to dst_restart_file */
	  if(time_exist) {
	    if(strcmp(varname, TILE_INDEX_NAME) == 0 || strcmp(varname, COHORT_INDEX_NAME) == 0 ||
	       strcmp(varname, LON_NAME) == 0 || strcmp(varname, LAT_NAME) == 0) {
	      int varsize;
	  
	      /* here we are assuming the vanname will have dimension name as its varname. */
	      if( has_taxis[l] ) mpp_error("remap_land: TILE_INDEX, COHORT_INDEX, LON and LAT should not depend on time");
	      varsize = mpp_get_dimlen(fid_read, varname);
	      vid_read = mpp_get_varid(fid_read, varname);
	      data_dst = (double *)malloc(varsize*sizeof(double));
	      mpp_get_var_value( fid_read, vid_read, data_dst );
	      mpp_put_var_value( fid_dst, vid_dst, data_dst );
	      continue;
	    }
	  }  
	  else {
	    if(strcmp(varname, TILE_INDEX_NAME) == 0 || strcmp(varname, COHORT_INDEX_NAME) == 0) continue;
	  }
	  
	  if(has_taxis[l])
	    mpp_get_var_dimname(fid_src[0], vid_src, 1, dimname);
	  else
	    mpp_get_var_dimname(fid_src[0], vid_src, 0, dimname);

	  if( strcmp(dimname, TILE_INDEX_NAME) && strcmp(dimname, COHORT_INDEX_NAME) ) continue;

	  nz_dst = 1;
	  if(!time_exist) {  
	    if(ndim_src[l] == 2) {
	      mpp_get_var_dimname(fid_read, vid_read, 1, dimname);
	      nz_dst = mpp_get_dimlen(fid_read, dimname);
	    }
	    if(nz_dst != nz_src[l])
	      mpp_error("remap_land: the src_restart and dst_cold_restart have the different dimension length.");
	  }
	  
	  if(var_type[l] == MPP_INT) idata_src = (int *)malloc(nidx_tot_src*nz_src[l]*sizeof(int));
	  data_src = (double *)malloc(nidx_tot_src*nz_src[l]*sizeof(double));
	  pos = 0;
	  for(m=0; m<4; m++) {
	    start[m]=0; nread[m] = 1;
	  }

	  for(m=0; m<nface_src; m++) {
	    if(has_taxis[l]) {
	      start[0] = t; nread[1] = nidx_src[m];
	    }
	    else {
	      nread[0] = nidx_src[m];
	      nread[1] = nz_src[l];
	    }		
	    vid_src = mpp_get_varid(fid_src[m], varname);
	    if(var_type[l] == MPP_INT)
	      mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src+pos);
	    else
	      mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src+pos);
	    pos += (nidx_src[m]*nz_src[l]);
	  }
	  if( var_type[l] == MPP_INT )
	    for(m=0; m< nidx_tot_src*nz_src[l]; m++) data_src[m] = idata_src[m];
      
	  start_pos = (int *)malloc(nface_src*sizeof(int));
	  start_pos[0] = 0;
	  for(m=1; m<nface_src; m++) {
	    start_pos[m] = start_pos[m-1] + nz_src[l]*nidx_src[m-1];
	  }
	  data_dst = (double *)malloc(nxc_dst*nz_dst*sizeof(double));
	  for(k=0; k<nz_dst; k++) for(m=0; m<nxc_dst; m++) {
	    int face, lll;
	    if(idx_map[m] < 0) continue;
	    face = face_map[m];
	    lll = start_pos[face]+k*nidx_src[face]+idx_map[m];
	    data_dst[k*nxc_dst+m] = data_src[lll];
	  }
	  /* pack the data and write to output */
	  gdata_dst = (double *) malloc(nidx_dst*nz_dst*sizeof(double));
	  mpp_global_field_double_3D(Dom_dst, nxc_dst, nyc_dst, nz_dst, data_dst, gdata_dst);
	  for(m=0; m<4; m++) {
	    start[m]=0; nwrite[m] = 1;
	  }
	  if(has_taxis[l]){
	    start[0] = t; nwrite[1] = nidx_dst;
	  }
	  else {
	    nwrite[0] = nidx_dst; nwrite[1] = nz_dst;
	  }
	
	  if( var_type[l] == MPP_INT ){
	    idata_dst = (int *) malloc(nidx_dst*nz_dst*sizeof(int));
	    for(m=0; m<nidx_dst*nz_dst; m++) idata_dst[m] = gdata_dst[m];
	    mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_dst);
	  }
	  else
	    mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, gdata_dst);

	  if(var_type[l] == MPP_INT) {
	    free(idata_dst);
	    free(idata_src);
	  }
	  free(data_src);
	  free(start_pos);
	  free(data_dst);
	  free(gdata_dst);
	}
      }
      mpp_close(fid_dst);
      mpp_close(fid_read);
      free(x_dst);
      free(y_dst);
      free(idx_tmp);
      free(idx_dst);
      free(idx_type_dst);    
  
    }
    free(x_tmp);
    free(y_tmp);
  }

  {
    int n;
    for(n=0; n<nface_src; n++) mpp_close(fid_src[n]);
  }
    
  /* release memory */
  free(x_src);
  free(y_src);

  free(fid_src);
  free(nidx_src);
  free(idx_src);

  free(idx_type_src);
  free(start_pos);

  free(has_taxis);
  free(ndim_src);
  free(nz_src);
  free(var_type);
  if(time_exist) free(time_data);
    
  
  if(mpp_pe() == mpp_root_pe() ) printf("\n******** Successfully run remap_land***********\n");

  mpp_end();

  return 0;
}; //main  

  
void get_actual_file_name(int nface, int face, const char *file_orig, char *file)
{
  if(nface == 1)
    strcpy(file, file_orig);
  else
    sprintf(file, "%s.tile%d.nc", file_orig, face+1);

}


  
  
/********************************************************************
 void full_search_nearest
 search the nearest point from the first of source grid to the last.
********************************************************************/
void full_search_nearest(int nface_src, const int *npts_src, const double *lon_src, const double *lat_src,
			 const int *tile_type_src, int npts_dst, const double *lon_dst,
			 const double *lat_dst, const int *tile_type_dst, int *idx_map, int *face_map)
{
  int i_dst, i_src, cur_tile_type;
  int pos, m, idx_cur, face_cur;
  double d_cur,d;
  double p1[2], p2[2];  

  for(i_dst=0; i_dst<npts_dst; i_dst++) {
    cur_tile_type = tile_type_dst[i_dst];
    if( cur_tile_type < 0) {
      idx_map[i_dst] = -1;
      face_map[i_dst] = -1;
      continue;
    }
    d_cur = -1;
    pos = 0;
    p1[0]= lon_dst[i_dst];
    p1[1]= lat_dst[i_dst];
    for(m=0; m<nface_src; m++) {
      for(i_src=0; i_src<npts_src[m]; i_src++) {
	if(tile_type_src[pos+i_src] != cur_tile_type) continue;
        p2[0] = lon_src[pos+i_src];
        p2[1] = lat_src[pos+i_src];
	d = great_circle_distance(p1, p2);
	if( d_cur < 0 || d<d_cur) {
	  idx_cur  = i_src;
	  face_cur = m;
	  d_cur    = d;
	}
      }
      pos += npts_src[m];
    }
    if(d_cur < 0)
      mpp_error("remap_land(full_search_nearest): no nearest point is found, the possible reason is that all the source points does not match the tile value of this point.");
    else {
      idx_map[i_dst]  = idx_cur;
      face_map[i_dst] = face_cur;
    }
  }

} ; /* full_search_nearest */

    
