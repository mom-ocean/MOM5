#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "constant.h"
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#define  TILE_INDEX_NAME "tile_index"
#define  COHORT_INDEX_NAME "cohort_index"
#define  COHORT_NAME "cohort"
#define  TILE_NAME "tile"
#define  LON_NAME  "lon"
#define  LAT_NAME  "lat"
#define  TIMENAME  "time"
#define  FRAC_NAME "frac"
#define  SOIL_NAME "soil"
#define  VEGN_NAME "vegn"
#define  GLAC_NAME "glac"
#define  LAKE_NAME "lake"
#define  LEVEL_NAME "zfull"
#define  N_ACCUM_NAME "n_accum"
#define  NMN_ACM_NAME "nmn_acm"
#define  THETA_AV_NAME "theta_av"
#define  THETA_AV_PHEN_NAME "theta_av_phen"
#define  THETA_AV_FIRE_NAME "theta_av_fire"
#define  PSIST_AV_NAME "psist_av"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)

char *usage[] = {
  "",
  "                                                                                ",
  "                                                                                ",
  "                    Usage of remap_land                                         ",
  "                                                                                ",
  "   remap_land --src_mosaic src_mosaic --src_restart src_restart                 ",
  "              --dst_mosaic dst_msoaic --dst_restart dst_restart                 ",
  "              --dst_cold_restart dst_cold_restart                               ",
  "              --land_src_restart land_src_restart                               ",
  "              --land_dst_cold_restart land_dst_cold_restart                     ",
  "              [--remap_file remap_file] [--print_memory]                        ",             
  "                                                                                ",
	  " remap_land remap land restart file from one mosaic grid to another mosaic grid ",
  " remap_land takes the following flag  s,                                          ",
  "                                                                                  ",
  " REQUIRED:                                                                        ",
  "                                                                                  ",
  " --src_mosaic src_mosaic      specify the source mosaic information. This file    ",
  "                              contains list of tile files which specify the grid  ",
  "                              information for each tile.                          ",
  "                                                                                  ",
  " --dst_mosaic dst_mosaic      specify the destination mosaic information. This    ",
  "                              file contains list of tile files which specify the  ",
  "                              grid information for each tile.                     ",
  "                                                                                  ",
  " --src_restart src_restart    specify the source restart file.                    ",
  "                                                                                  ",
  " --dst_restart dst_restart    specify the restart file to be generated on         ",
  "                              destination grid.                                   ", 
  "                                                                                  ",
  " --dst_cold_restart file      specify the cold restart file destination grid.     ",
  "                              This is the input file. The dst_cold_restart_file   ",
  "                              could be obtained by running the experiment         ",
  "                              for 1 day with --dst_mosaic using cold restart.     ",
  "                                                                                  ",
  " --file_type file_type        spefify file type. Its value could be 'land',       ",
  "                              'cana', 'snow', 'glac', 'lake', 'soil' or 'vegn'.   ",
  "                              when file_type is 'cana' or 'snow', need to         ",
  "                              specify --land_src_restart and --land_cold_restart  ",
  "                                                                                  ",   
  " OPTIONAL FLAGS                                                                   ",
  "                                                                                  ",
  "                                                                                  ",
  " --land_src_restart file      specify the source file of land.res.nc. It is       ",
  "                              required when the restart file is snow.res or       ",
  "                              cana.res                                            ",
  "                                                                                  ",
  " --land_dst_cold_restart file specify the destination file of land.res.nc. It is  ",
  "                              required when the restart file is snow.res or       ",
  "                              cana.res                                            ",
  "                                                                                  ",
  " --remap_file remap_file    specify the file name that saves remapping            ",
  "                            information. If remap_file is specified and the       ",
  "                            file does not exist, remapping information will be    ",
  "                            calculated and stored in remap_file. If remap_file    ",
  "                            is specified and the file exists, remapping           ",
  "                            information will be read from remap_file.             ",
  "                                                                                  ",
  " --print_memory             debug memory usage when it is set                     ",
  "  ",
  "                                                                                  ",
  " Example: remap land restart from C48 grid onto C180 grid.                        ",
  "                                                                                  ",
  "   remap_land --src_mosaic C48/C48_mosaic.nc                                      ",
  "              --dst_mosaic C180/C180_mosaic.nc                                    ",
  "              --src_restart src_restart/soil.res                                  ",
  "              --dst_cold_restart dst_cold_restart/soil.res                        ",
  "              --dst_restart dst_restart/soil.res                                  ",
  "              --land_src_restart src_restart/land.res                             ",
  "              --land_dst_cold_restart dst_cold_restart/land.res                   ",
  "                                                                                  ",
  NULL };

char grid_version[] = "0.2";
char tagname[] = "$Name: tikal $";

double distance(double lon1, double lat1, double lon2, double lat2);
void get_actual_file_name(int nface, int face, const char *file_orig, char *file);
void get_land_tile_info(int fid, const char *name1, const char *name2, int nidx, const int *idx_in, const double *frac_in,
			int nx, int ny, int ntile, int isc, int iec, int *count, double *frac, int *tag1, int *tag2, int *idx, int all_tile);
void full_search_nearest(int nface_src, int npts_src, const double *lon_src, const double *lat_src,
			 const int *mask_src, int npts_dst, const double *lon_dst, const double *lat_dst,
			 const int *mask_dst, int *idx_map, int *face_map);
void compress_int_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
		       const int *data, int *data_global, int all_tile );
void compress_double_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
			  const double *data, double *data_global, int all_tile );
const int LANDTYPE = 1; 
const int SOILTYPE = 2;
const int GLACTYPE = 3; 
const int LAKETYPE = 4; 
const int VEGNTYPE = 5; 
const int CANATYPE = 6; 
const int SNOWTYPE = 7; 

int main(int argc, char* argv[])
{

  char *src_mosaic            = NULL;
  char *dst_mosaic            = NULL;
  char *src_restart_file      = NULL;
  char *dst_restart_file      = NULL;
  char *dst_cold_restart      = NULL;
  char *land_src_restart      = NULL;
  char *land_cold_restart     = NULL;
  char *remap_file            = NULL;
  char *file_type             = NULL;
  
  int face_dst;
  int nface_src = 0, nface_dst = 0;
  int nx_src, ny_src, nx_dst, ny_dst;
  int ntile_src, ntile_cold, ntile_dst;
  int nidx_tot_src, ncohort;
  int *fid_src=NULL;
  int *cohort_data=NULL, *tile_axis_data=NULL;
  int *nidx_src=NULL, *nidx_land_src=NULL;

  int *idx_soil_src=NULL, *idx_glac_src=NULL, *idx_lake_src=NULL;
  int *soil_count_src=NULL, *glac_count_src=NULL, *lake_count_src=NULL;
  int *soil_tag_src=NULL, *glac_tag_src=NULL, *lake_tag_src=NULL, *vegn_tag_src=NULL;
  double *soil_frac_src=NULL, *glac_frac_src=NULL, *lake_frac_src=NULL;  

  int filetype;
  char   history[1280];
  double *x_src=NULL, *y_src=NULL;
  int time_exist, zaxis_exist, ntime, l;
  int *has_taxis=NULL, *var_type=NULL, *ndim_src=NULL, *nz_src=NULL;
  double *time_data=NULL;
  int has_glac=0, has_lake=0;
  int src_has_tile=0, cold_has_tile=0;
  int src_has_cohort=0, cold_has_cohort=0;
  int src_has_theta_av=0, cold_has_theta_av=0;
  int src_has_theta_av_phen=0, cold_has_theta_av_phen=0;
  int src_has_theta_av_fire=0, cold_has_theta_av_fire=0;
  int src_has_psist_av=0, cold_has_psist_av=0;
  int nz;
  double *z_axis_data=NULL;
  int    print_memory=0;
  int    use_all_tile;
  
  int    npes, nvar_src;
  int    option_index, c;
  int    errflg = (argc == 1);  
  /*
   * process command line
   */

  static struct option long_options[] = {
    {"src_mosaic",            required_argument, NULL, 'a'},
    {"dst_mosaic",            required_argument, NULL, 'b'},
    {"src_restart",           required_argument, NULL, 'c'},
    {"dst_restart",           required_argument, NULL, 'o'},
    {"dst_cold_restart",      required_argument, NULL, 'd'},
    {"file_type",             required_argument, NULL, 't'},
    {"land_src_restart",      required_argument, NULL, 'l'},
    {"land_cold_restart",     required_argument, NULL, 'm'},
    {"remap_file",            required_argument, NULL, 'r'},
    {"print_memory",          no_argument,       NULL, 'p'},
    {"help",                  no_argument,       NULL, 'h'},
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
    case 'l':
      land_src_restart = optarg;
      break;
    case 'm':
      land_cold_restart = optarg;
      break;      
    case 't':
      file_type = optarg;
      break;
    case 'r':
      remap_file = optarg;
      break;
    case 'p':
      print_memory = 1;
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
  if( !file_type) errflg++;
  
  if (errflg) {
    char **u = usage;
    if(mpp_pe() == mpp_root_pe()) {
      while (*u) { fprintf(stderr, "%s\n", *u); u++; }
      if(!src_mosaic) mpp_error("remap_land: src_mosaic is not specified");
      if(!dst_mosaic) mpp_error("remap_land: dst_mosaic is not specified");
      if(!src_restart_file) mpp_error("remap_land: src_restart_file is not specified");
      if(!dst_restart_file) mpp_error("remap_land: dst_restart_file is not specified");
      if(!dst_cold_restart) mpp_error("remap_land: dst_cold_restart is not specified");
      if(!file_type) mpp_error("remap_land: file_type is not specified");
    }
    mpp_error("remap_land: check the command line arguments");
  }

  if(print_memory) print_mem_usage("at the begining of remap_land");
  
  /* write out arguments */
  if(mpp_pe() == mpp_root_pe()) {
    printf("src_mosaic       is %s\n", src_mosaic);
    printf("dst_mosaic       is %s\n", dst_mosaic);
    printf("src_restart_file is %s\n", src_restart_file);
    printf("dst_restart_file is %s\n", dst_restart_file);
    printf("dst_cold_restart is %s\n", dst_cold_restart);
    printf("file_type        is %s\n", file_type);
    if(land_src_restart) {
      printf("land_src_restart is %s\n", land_src_restart);
      printf("land_cold_restart is %s\n", land_cold_restart);
    }
    else {
      printf("land_src_restart is not specified\n");
      printf("land_cold_restart is not specified\n");
    }
  }
    
  /* file type must be the land, cana, snow, soil, lake, glac, vegn */
  if(!strcmp(file_type, "land"))
    filetype = LANDTYPE;
  else if(!strcmp(file_type, "soil"))
    filetype = SOILTYPE;
  else if(!strcmp(file_type, "cana") )
    filetype = CANATYPE;
  else if(!strcmp(file_type, "snow"))
    filetype = SNOWTYPE;
  else if(!strcmp(file_type, "lake"))
    filetype = LAKETYPE;
  else if(!strcmp(file_type, "glac"))
    filetype = GLACTYPE;
  else if(!strcmp(file_type, "vegn") )
    filetype = VEGNTYPE;
  else
    mpp_error("remap_land: invalid option in --file_type");

  /* when file_type is cana or snow, land_src_restart and land_cold_restart must be specified */
  if(filetype == LANDTYPE) {
    if(land_src_restart) mpp_error("remap_land: land_src_restart must not be specified "
				   "when file_type is 'land'");
    if(land_cold_restart) mpp_error("remap_land: land_cold_restart must not be specified "
				    "when file_type is 'land'");
    land_src_restart = src_restart_file;
    land_cold_restart = dst_cold_restart;
  }
  else {
    if(!land_src_restart) mpp_error("remap_land: land_src_restart must be specified "
				    "when file_type is not 'land'");
    if(!land_cold_restart) mpp_error("remap_land: land_cold_restart must be specified "
				     "when file_type is not 'land'");
  }
    
  npes = mpp_npes();
  
  /*----------------------------------------------------------------------------
    get source and destination grid size
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
      if(nx[n] != nx[0] || ny[n] != ny[0])
	mpp_error("remap_land: all the faces of source grid should have the same number of grid points");
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
      if(nx[n] != nx[0] || ny[n] != ny[0])
	mpp_error("remap_land: all the faces of destination grid should have the same number of grid points");
    }
    nx_dst = nx[0];
    ny_dst = ny[0];
    free(nx);
    free(ny);    
  }

  
  /*-----------------------------------------------------------------------------
    read the source tile_index size and source grid
    ----------------------------------------------------------------------------*/
  {
    int    n, pos;
    
    fid_src      = (int *)malloc(nface_src*sizeof(int));
    nidx_src     = (int *)malloc(nface_src*sizeof(int));
    nidx_land_src     = (int *)malloc(nface_src*sizeof(int));
    nidx_tot_src = 0;
    for(n=0; n<nface_src; n++) {
      int nlon, nlat;
      char file[512];
      int fid_land;
      
      get_actual_file_name(nface_src, n, src_restart_file, file);
      fid_src[n] = mpp_open(file, MPP_READ);
      nlon = mpp_get_dimlen(fid_src[n], LON_NAME);
      nlat = mpp_get_dimlen(fid_src[n], LAT_NAME);
      if(nx_src != nlon) mpp_error("remap_land: mismatch on the longitude dimension size "
				   "between source mosaic grid file and src_restart");
      if(ny_src != nlat) mpp_error("remap_land: mismatch on the latitude dimension size "
				   "between source mosaic grid file and src_restart");    
      nidx_src[n] = mpp_get_dimlen(fid_src[n], TILE_INDEX_NAME);
      
      /* get the size of tile_index in land_src_restart */
      if(filetype == LANDTYPE)
	nidx_land_src[n] = nidx_src[n];
      else {
	get_actual_file_name(nface_src, n, land_src_restart, file);
	fid_land = mpp_open(file, MPP_READ);;
	nidx_land_src[n] = mpp_get_dimlen(fid_land, TILE_INDEX_NAME);
	mpp_close(fid_land);
      }

      nidx_tot_src += nidx_src[n];
    }
  
    x_src        = (double *)malloc(nface_src*nx_src*ny_src*sizeof(double));
    y_src        = (double *)malloc(nface_src*nx_src*ny_src*sizeof(double));

    for(n=0; n<nface_src; n++) {
      int i;
    
      read_mosaic_grid_data(src_mosaic, "x", nx_src, ny_src, x_src+n*nx_src*ny_src, n, 1, 1);
      read_mosaic_grid_data(src_mosaic, "y", nx_src, ny_src, y_src+n*nx_src*ny_src, n, 1, 1);   
    }
    /* convert to radians */
    for(n=0; n<nface_src*nx_src*ny_src; n++) {
      x_src[n] *= D2R;
      y_src[n] *= D2R;
    }
    
  }

  /*-----------------------------------------------------------------------------
    Get the time information, time axis only exists for static vegetation.
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

  /* check if z-axis exist or not */
  zaxis_exist = mpp_dim_exist(fid_src[0], LEVEL_NAME);
  nz = 1;
  if(zaxis_exist) {
    int vid;

    nz = mpp_get_dimlen(fid_src[0], LEVEL_NAME);
    z_axis_data = (double *)malloc(nz*sizeof(double));
    vid = mpp_get_varid(fid_src[0], LEVEL_NAME);
    mpp_get_var_value(fid_src[0], vid, z_axis_data);
  }
    
  /*-----------------------------------------------------------------------------
    loop through each variable of the source data to see get dimension of each variable
    ----------------------------------------------------------------------------*/
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
      mpp_get_var_dimname(fid_src[0], vid, 0, varname);
      if( strcmp(varname, LEVEL_NAME) ) mpp_error("remap_land: The vertical dimension name should be zfull, contact developer");
      nz_src[l] = mpp_get_dimlen(fid_src[0], varname);
    }
  }

  /*------------------------------------------------------------------------------
    get the cohort data, currently it only contains in vegn data
    -----------------------------------------------------------------------------*/
  {
    
    ncohort = 0;
    src_has_cohort=0;
    if(filetype == VEGNTYPE) {
      int dimsize, m, vid, i;
      int *tmp;
      if(!mpp_dim_exist(fid_src[0], COHORT_NAME)) mpp_error("remap_land: dimension cohort should exist when file_type is vegn");
      ncohort = mpp_get_dimlen(fid_src[0], COHORT_NAME);
      if(ncohort != 1) mpp_error("remap_land: size of cohort should be 1, contact developer");
      src_has_theta_av = mpp_var_exist(fid_src[0], THETA_AV_NAME);
      src_has_theta_av_phen = mpp_var_exist(fid_src[0], THETA_AV_PHEN_NAME);
      src_has_theta_av_fire = mpp_var_exist(fid_src[0], THETA_AV_FIRE_NAME);
      src_has_psist_av = mpp_var_exist(fid_src[0], PSIST_AV_NAME);
      
      if(mpp_var_exist(fid_src[0], COHORT_NAME)) {
	src_has_cohort = 1;
	cohort_data = (int *)malloc(ncohort*sizeof(int));
	vid = mpp_get_varid(fid_src[0], COHORT_NAME);
	mpp_get_var_value(fid_src[0], vid, cohort_data);
	for(m=1; m<nface_src; m++) {
	  dimsize = mpp_get_dimlen(fid_src[m], COHORT_NAME);
	  if(dimsize != ncohort)mpp_error("remap_land: the dimension size of cohort is different between faces");
	  tmp = (int *)malloc(ncohort*sizeof(int));
	  vid = mpp_get_varid(fid_src[m], COHORT_NAME);
	  mpp_get_var_value(fid_src[m], vid, tmp);
	  for(i=0; i<ncohort; i++) {
	    if(cohort_data[i] != tmp[i]) mpp_error("remap_land: cohort value is different between faces");
	  }
	  free(tmp);
	}
      }
    }
  }

  
  /*-----------------------------------------------------------------------------
    Check if the cold_restart has lake or glac. soil is required to be existed
    ---------------------------------------------------------------------------*/
  {
    int max_nidx, face_dst, vid, i;
    int *fid=NULL;
    char file[512];
    int *nidx=NULL;
    int *tmp=NULL;

    cold_has_tile = 0;
    cold_has_cohort = 0;
    has_lake = 0;
    has_glac = 0;
    max_nidx = 0;
    fid = (int *)malloc(nface_dst*sizeof(int));
    nidx = (int *)malloc(nface_dst*sizeof(int));

    /* get ntile for dst_cold_restart */
    {
      int fid_cold;
      
      get_actual_file_name(nface_dst, 0, dst_cold_restart, file);      
      fid_cold = mpp_open(file, MPP_READ);
      cold_has_tile = mpp_var_exist(fid_cold, TILE_NAME);
      if(filetype==VEGNTYPE) {
	cold_has_cohort = mpp_var_exist(fid_cold, COHORT_NAME);
	cold_has_theta_av = mpp_var_exist(fid_cold, THETA_AV_NAME);
	cold_has_theta_av_phen = mpp_var_exist(fid_cold, THETA_AV_PHEN_NAME);
	cold_has_theta_av_fire = mpp_var_exist(fid_cold, THETA_AV_FIRE_NAME);
	cold_has_psist_av = mpp_var_exist(fid_cold, PSIST_AV_NAME);
      }
      mpp_close(fid_cold);
    }


    for(face_dst=0; face_dst<nface_dst; face_dst++) {
      get_actual_file_name(nface_dst, face_dst, land_cold_restart, file);      
      fid[face_dst] = mpp_open(file, MPP_READ);
      nidx[face_dst] = mpp_get_dimlen(fid[face_dst], TILE_INDEX_NAME);
      if(nidx[face_dst] > max_nidx) max_nidx = nidx[face_dst];
    }
      
    ntile_cold   = mpp_get_dimlen(fid[0], TILE_NAME);
    for(face_dst=1; face_dst<nface_dst; face_dst++) {
      int ntile;
      ntile = mpp_get_dimlen(fid[face_dst], TILE_NAME);
      if(ntile != ntile_cold) mpp_error("remap_land: mismatch of tile dimension between different faces");
    }
      
    tmp = (int *)malloc(max_nidx*sizeof(int));
    for(face_dst=0; face_dst<nface_dst; face_dst++) {
      vid = mpp_get_varid(fid[face_dst], GLAC_NAME);
      mpp_get_var_value(fid[face_dst], vid, tmp);
      for(i=0; i<nidx[face_dst]; i++) {
	if(tmp[i] > 0) {
	  has_glac = 1;
	  goto GLAC_CHECK;
	}
      }
    }
  GLAC_CHECK:
    if( mpp_pe() == mpp_root_pe() ) {
      if(has_glac)
         printf("remap_land: there is glac in cold restart file\n");
      else
	printf("remap_land: there is no glac in cold restart file\n");
    }

    for(face_dst=0; face_dst<nface_dst; face_dst++) {
      vid = mpp_get_varid(fid[face_dst], LAKE_NAME);
      mpp_get_var_value(fid[face_dst], vid, tmp);
      for(i=0; i<nidx[face_dst]; i++) {
	if(tmp[i] > 0) {
	  has_lake = 1;
	  goto LAKE_CHECK;
	}
      }
    }
  LAKE_CHECK:
    if( mpp_pe() == mpp_root_pe() ) {
      if(has_lake)
         printf("remap_land: there is lake in cold restart file\n");
      else
	printf("remap_land: there is no lake in cold restart file\n");
    }
    for(face_dst=0; face_dst<nface_dst; face_dst++) mpp_close(fid[face_dst]);
    
    free(tmp);
    free(nidx);
    free(fid);
  }
    
  
  /*-----------------------------------------------------------------------------
    get the tile data
    ----------------------------------------------------------------------------*/  
  if(!mpp_dim_exist(fid_src[0], TILE_NAME)) mpp_error("remap_land: dimension tile should exist");

  ntile_src   = mpp_get_dimlen(fid_src[0], TILE_NAME);
  ntile_dst = ntile_src;
  /*  if(!has_glac) ntile_dst--; */
  /* if(!has_lake) ntile_dst--;  */

  src_has_tile = 0;
  if(mpp_var_exist(fid_src[0], TILE_NAME)){
    int vid, m, i, dimsize;
    int *tmp=NULL;

    src_has_tile=1;
    tile_axis_data   = (int *)malloc(ntile_src*sizeof(int));
    vid = mpp_get_varid(fid_src[0], TILE_NAME);
    mpp_get_var_value(fid_src[0], vid, tile_axis_data);  
    for(m=1; m<nface_src; m++) {
      dimsize = mpp_get_dimlen(fid_src[m], TILE_NAME);
      if(dimsize != ntile_src)mpp_error("remap_land: the dimension size of tile is different between faces");
      tmp = (int *)malloc(ntile_src*sizeof(int));
      vid = mpp_get_varid(fid_src[m], TILE_NAME);
      mpp_get_var_value(fid_src[m], vid, tmp);
      for(i=0; i<ntile_src; i++) {
	if(tile_axis_data[i] != tmp[i]) mpp_error("remap_land: tile value is different between faces");
      }
      free(tmp);      
    }
  }

  /*-----------------------------------------------------------------------
    get source tile information. 
    get the tile information for soil, lake and glac 
    The following will use large memroy for high resolution source restart file,
    May come back to rewrite to decrease memory usage.
    ---------------------------------------------------------------------*/
  {
    soil_count_src = (int *)malloc(nface_src*nx_src*ny_src*sizeof(int));
    soil_frac_src =  (double *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(double));
    soil_tag_src =  (int *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(int));
    vegn_tag_src =  (int *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(int));
    idx_soil_src = (int *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(int));
  
    if(has_glac) {
      glac_count_src = (int *)malloc(nface_src*nx_src*ny_src*sizeof(int));
      glac_frac_src =  (double *)malloc(nface_src*nx_src*ny_src*sizeof(double)); /* at most one tile for glac */
      glac_tag_src =  (int *)malloc(nface_src*nx_src*ny_src*sizeof(int)); /* at most one tile for glac */
      idx_glac_src = (int *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(int));
    }
    if(has_lake) {
      lake_count_src = (int *)malloc(nface_src*nx_src*ny_src*sizeof(int));
      lake_frac_src =  (double *)malloc(nface_src*nx_src*ny_src*sizeof(double)); /* at most one tile for glac */
      lake_tag_src =  (int *)malloc(nface_src*nx_src*ny_src*sizeof(int)); /* at most one tile for glac */
      idx_lake_src = (int *)malloc(nface_src*nx_src*ny_src*ntile_src*sizeof(int));
    }
    
    {
      double *frac_land_src=NULL;
      int *idx_src=NULL;
      int *idx_land_src=NULL;
      char file[512];
      
      int n, max_nidx, pos1, pos2;
      max_nidx = 0;
    
      for(n=0; n<nface_src; n++) {
	if(nidx_land_src[n] >max_nidx) max_nidx = nidx_land_src[n];
      }
      frac_land_src =  (double *)malloc(max_nidx*sizeof(double));
      idx_land_src =  (int *)malloc(max_nidx*sizeof(int));
      
      pos1 = 0;
      pos2 = 0;
      use_all_tile = 0;
      if(filetype==LANDTYPE || filetype==CANATYPE || filetype==SNOWTYPE) use_all_tile = 1;
	
      for(n=0; n<nface_src; n++) {
	int vid, fid;

	if(filetype==LANDTYPE) {
	  fid = fid_src[n];
	}
	else {
	  get_actual_file_name(nface_src, n, land_src_restart, file);
          fid = mpp_open(file, MPP_READ);
	}
	vid = mpp_get_varid(fid, FRAC_NAME);
	mpp_get_var_value(fid, vid, frac_land_src);
	vid = mpp_get_varid(fid, TILE_INDEX_NAME);
	mpp_get_var_value(fid, vid, idx_land_src);
      
	/* soil, vegn */
	get_land_tile_info(fid, SOIL_NAME, VEGN_NAME, nidx_land_src[n], idx_land_src, frac_land_src,
			   nx_src, ny_src, ntile_src, 0, nx_src*ny_src-1, soil_count_src+pos1, soil_frac_src+pos2,
			   soil_tag_src+pos2, vegn_tag_src+pos2, idx_soil_src+pos2, use_all_tile);
    
	/* glac */
	if(has_glac) 
	  get_land_tile_info(fid, GLAC_NAME, NULL, nidx_land_src[n], idx_land_src, frac_land_src,
			     nx_src, ny_src, 1, 0, nx_src*ny_src-1, glac_count_src+pos1, glac_frac_src+pos1,
			     glac_tag_src+pos1, NULL, idx_glac_src+pos1, use_all_tile);
	/* lake */
	if(has_lake)
	  get_land_tile_info(fid, LAKE_NAME, NULL, nidx_land_src[n], idx_land_src, frac_land_src,
			     nx_src, ny_src, 1, 0, nx_src*ny_src-1, lake_count_src+pos1, lake_frac_src+pos1,
			     lake_tag_src+pos1, NULL, idx_lake_src+pos1, use_all_tile);
    

	pos1 += nx_src*ny_src;
	pos2 += nx_src*ny_src*ntile_src;
	if(filetype!=LANDTYPE) mpp_close(fid);

	/* make sure the tile_index consistency between src_restart_file and land_src_restart */
	if(filetype!=LANDTYPE) {
	  if(filetype==CANATYPE || filetype==SNOWTYPE) {
	    if(nidx_land_src[n] != nidx_src[n])
	      mpp_error("remap_land: size of tile_index mismatch between src_restart_file "
			"and land_src_restart for 'cana' or 'snow'");
	  }
	  else {
	    if(nidx_land_src[n] < nidx_src[n])
	      mpp_error("remap_land: size of tile_index mismatch between src_restart_file "
			"and land_src_restart for 'soil', 'vegn', 'glac' or 'lake'");
	  }
	  idx_src =  (int *)malloc(max_nidx*sizeof(int));
	  vid = mpp_get_varid(fid_src[n], TILE_INDEX_NAME);
	  mpp_get_var_value(fid_src[n], vid, idx_src);
	  if(filetype==CANATYPE || filetype==SNOWTYPE) {
	    int i;
	    for(i=0; i<nidx_land_src[n]; i++)
	      if(idx_src[i] != idx_land_src[i]) mpp_error("remap_land: mismatch of tile_index between src_restart_file "
							  "and land_src_restart for 'soil', 'vegn', 'glac' or 'lake'");
	  }
	  else if(filetype==SOILTYPE || filetype==VEGNTYPE) {
	    int i, j, k, l, p, m, idx, p2;
	    for(m=0; m<ntile_src; m++) {
	      for(l=0; l<nx_src*ny_src; l++) {
		p = n*nx_src*ny_src + l;
		idx = idx_soil_src[ntile_src*p+m];
		if(idx==MPP_FILL_INT) continue;
		if(idx> nidx_src[n])mpp_error("remap_land: idx is out of bound for soil consistency check");
		idx = idx_src[idx];
		i = idx%nx_src;
		k = idx/nx_src;
		j = k%ny_src;
		k = k/ny_src;
		p2 = n*nx_src*ny_src + j*nx_src+i;
		if(p!=p2) mpp_error("remap_land: mismatch of tile_index for src soil check");
		if(soil_tag_src[ntile_src*p+k] == MPP_FILL_INT) mpp_error("remap_land: soil_tag_src is not defined for src soil check");
	      }
	    }
	  }
	}
      }
      free(frac_land_src);
      free(idx_land_src);
      if(filetype!=LANDTYPE)free(idx_src);
    }
  }

  /* define history attribute */
  {
    int n;
    strcpy(history,argv[0]);
    for(n=1;n<argc;n++) {
      strcat(history, " ");
      strcat(history, argv[n]);
    }
  }

  if(print_memory) print_mem_usage("After initialization of source information");
  
  /*------------------------------------------------------------------------------------------
    loop through each face of destination grid, first read the grid, then read the tile_index,
    then find the remapping index, then setup metadata for the destination file,
    last do the remapping and write out the data to dst_restart_file
    ----------------------------------------------------------------------------------------*/
  {
    double *data_dst=NULL, *data_src=NULL;
    int    *idata_dst=NULL, *idata_src=NULL;
    int    *start_pos=NULL;    
    double *x_tmp=NULL, *y_tmp=NULL;
    double *lon_axis_dst=NULL, *lat_axis_dst=NULL;
    double *x_dst=NULL, *y_dst=NULL;

    int *idx_dst=NULL;
    int    *glac_tag_dst=NULL, *lake_tag_dst=NULL, *soil_tag_dst=NULL, *vegn_tag_dst=NULL;
    int    *idx_map_soil=NULL, *face_map_soil=NULL;
    int    *idx_map_glac=NULL, *face_map_glac=NULL;
    int    *idx_map_lake=NULL, *face_map_lake=NULL;
    int    *land_idx_map=NULL, *land_face_map=NULL;
    double *glac_frac_cold=NULL, *lake_frac_cold=NULL, *soil_frac_cold=NULL, *tmp_frac_cold=NULL;
    int    *glac_count_cold=NULL, *lake_count_cold=NULL, *soil_count_cold=NULL;
    int    *glac_tag_cold=NULL, *lake_tag_cold=NULL, *soil_tag_cold=NULL, *vegn_tag_cold=NULL;
    int    *land_count_dst=NULL;
    int    *idx_map_land=NULL;
    double *land_frac_dst=NULL;
    
    int      isc_dst, iec_dst, jsc_dst, jec_dst, nxc_dst, nyc_dst;
    int      n, pos, has_int_var;
    int      layout[2];
    domain2D Dom_dst;
    
    lon_axis_dst = (double *)malloc(nx_dst*sizeof(double));
    lat_axis_dst = (double *)malloc(ny_dst*sizeof(double));
    x_tmp = (double *) malloc(nx_dst*ny_dst*sizeof(double));
    y_tmp = (double *) malloc(nx_dst*ny_dst*sizeof(double));
    start_pos = (int *)malloc(nface_src*sizeof(int));

    has_int_var = 0;
    for(l=0; l<nvar_src; l++) {
      if(var_type[l] == MPP_INT) {
	has_int_var = 1;
	break;
      }
    }
    
    idata_src = (int *)malloc(nidx_tot_src*sizeof(int));
    data_src = (double *)malloc(nidx_tot_src*sizeof(double));
    
    /* setup domain */
    layout[0] = npes;
    layout[1] = 1;
    mpp_define_domain2d( nx_dst*ny_dst, 1, layout, 0, 0, &Dom_dst);
    mpp_get_compute_domain2d( Dom_dst, &isc_dst, &iec_dst, &jsc_dst, &jec_dst);
    /* for this domain, jsc_dst should equal jec_dst */
    if(jsc_dst != jec_dst) {
      mpp_error("remap_land: This is a 1-D domain decomposition, jsc_dst must equal to jec_dst");
    }
    nxc_dst = iec_dst - isc_dst + 1;
    nyc_dst  = 1;

    data_dst = (double *)malloc(ntile_dst*nxc_dst*sizeof(double));
    idata_dst = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    x_dst = (double *)malloc(nxc_dst*sizeof(double));
    y_dst = (double *)malloc(nxc_dst*sizeof(double));

    land_idx_map    = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    land_face_map   = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    land_count_dst = (int *)malloc(nxc_dst*sizeof(int));
    idx_dst        = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    
    soil_tag_dst   = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    vegn_tag_dst   = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
    land_frac_dst  = (double *)malloc(ntile_dst*nxc_dst*sizeof(double));
    soil_count_cold = (int *)malloc(nxc_dst*sizeof(int));
    soil_frac_cold  = (double *) malloc(nxc_dst*sizeof(double));
    tmp_frac_cold   = (double *) malloc(ntile_cold*nxc_dst*sizeof(double));
    soil_tag_cold   = (int *)malloc(ntile_cold*nxc_dst*sizeof(int));
    vegn_tag_cold   = (int *)malloc(ntile_cold*nxc_dst*sizeof(int));
    idx_map_soil    = (int *)malloc(nxc_dst*sizeof(int));
    face_map_soil   = (int *)malloc(nxc_dst*sizeof(int));

    if(has_glac){
      glac_tag_dst   = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
      glac_count_cold = (int *)malloc(nxc_dst*sizeof(int));
      glac_frac_cold  = (double *) malloc(nxc_dst*sizeof(double));
      glac_tag_cold   = (int *)malloc(nxc_dst*sizeof(int));
      idx_map_glac  = (int *)malloc(nxc_dst*sizeof(int));
      face_map_glac = (int *)malloc(nxc_dst*sizeof(int));
    }

    if(has_lake) {
      lake_tag_dst   = (int *)malloc(ntile_dst*nxc_dst*sizeof(int));
      lake_count_cold = (int *)malloc(nxc_dst*sizeof(int));
      lake_frac_cold  = (double *) malloc(nxc_dst*sizeof(double));
      lake_tag_cold   = (int *)malloc(nxc_dst*sizeof(int));
      idx_map_lake  = (int *)malloc(nxc_dst*sizeof(int));
      face_map_lake = (int *)malloc(nxc_dst*sizeof(int));
    }
    
    if(print_memory) print_mem_usage("before the loop of face_dst");    
    for(face_dst=0; face_dst<nface_dst; face_dst++) {
      int nlon, nlat, l, ll, i, j, k, t, p;
      int fid_dst, vid_dst, vid_src, vid_cold, fid_cold;
      int nidx_dst, nidx_dst_global;
      int nidx_cold, nvar_cold, vid;
      int fid_land_cold;
      double *frac_cold=NULL;
      int    *idx_cold=NULL;

      double *rdata_global=NULL;
      int    *idata_global=NULL;
      
      size_t start[4], nread[4], nwrite[4];
      char land_cold[512], file_dst[512], file_cold[512];
    
      get_actual_file_name(nface_dst, face_dst, dst_restart_file, file_dst);
      get_actual_file_name(nface_dst, face_dst, dst_cold_restart, file_cold);

      fid_cold = mpp_open(file_cold, MPP_READ);
      if(filetype == LANDTYPE) {
	nidx_cold = mpp_get_dimlen(fid_cold, TILE_INDEX_NAME);
	fid_land_cold = fid_cold;
      }
      else {
        get_actual_file_name(nface_dst, face_dst, land_cold_restart, land_cold);
	fid_land_cold = mpp_open(land_cold, MPP_READ);
	nidx_cold = mpp_get_dimlen(fid_land_cold, TILE_INDEX_NAME);
      }
            
      idx_cold = (int *)malloc(nidx_cold*sizeof(int));
      vid = mpp_get_varid(fid_land_cold, TILE_INDEX_NAME); 
      mpp_get_var_value(fid_land_cold, vid, idx_cold);

      /* get destination grid */
      read_mosaic_grid_data(dst_mosaic, "x", nx_dst, ny_dst, x_tmp, face_dst, 1, 1);
      read_mosaic_grid_data(dst_mosaic, "y", nx_dst, ny_dst, y_tmp, face_dst, 1, 1);

      for(i=0; i<nx_dst; i++) lon_axis_dst[i] = x_tmp[i];
      for(i=0; i<ny_dst; i++) lat_axis_dst[i] = y_tmp[i*nx_dst];
            
      fid_dst = mpp_open(file_dst, MPP_WRITE);

      nlon = mpp_get_dimlen(fid_cold, LON_NAME);
      nlat = mpp_get_dimlen(fid_cold, LAT_NAME);
      if(nx_dst != nlon) mpp_error("remap_land: mismatch on the longitude dimension size "
				   "between destination mosaic grid file and dst_cold_restart");
      if(ny_dst != nlat) mpp_error("remap_land: mismatch on the latitude dimension size "
				   "between destination mosaic grid file and dst_cold_restart");

      if(filetype != LANDTYPE) {
	nlon = mpp_get_dimlen(fid_land_cold, LON_NAME);
	nlat = mpp_get_dimlen(fid_land_cold, LAT_NAME);
	if(nx_dst != nlon) mpp_error("remap_land: mismatch on the longitude dimension size "
				     "between destination mosaic grid file and land_cold_restart");
	if(ny_dst != nlat) mpp_error("remap_land: mismatch on the latitude dimension size "
				     "between destination mosaic grid file and land_cold_restart");
      }

      /* get destination grid */
      for(i=isc_dst; i<=iec_dst; i++) {
	x_dst[i-isc_dst] = x_tmp[i]*D2R;
	y_dst[i-isc_dst] = y_tmp[i]*D2R;
      }
            
      for(i=0; i<nxc_dst; i++) land_count_dst[i] = 0;
      for(i=0; i<ntile_dst*nxc_dst; i++) land_idx_map[i] = -1;	
	
      {
	frac_cold = (double *)malloc(nidx_cold*sizeof(double));
	vid = mpp_get_varid(fid_land_cold, FRAC_NAME);
	mpp_get_var_value(fid_land_cold, vid, frac_cold);
      
	for(i=0; i<nxc_dst; i++) {
	  soil_frac_cold[i] = 0.0;
	}

	for(i=0; i<ntile_dst*nxc_dst; i++) {
	  soil_tag_dst[i] = MPP_FILL_INT;
	  vegn_tag_dst[i] = MPP_FILL_INT;
	  land_frac_dst[i] = MPP_FILL_DOUBLE;
	}
	
	if(has_glac){
	  for(i=0; i<ntile_dst*nxc_dst; i++) glac_tag_dst[i] = MPP_FILL_INT;
	}

	if(has_lake) {
	  for(i=0; i<ntile_dst*nxc_dst; i++) lake_tag_dst[i] = MPP_FILL_INT;
	}
      
	/* soil, vegn */
	get_land_tile_info(fid_land_cold, SOIL_NAME, VEGN_NAME, nidx_cold, idx_cold, frac_cold,
			   nx_dst, ny_dst, ntile_cold, isc_dst, iec_dst, soil_count_cold, tmp_frac_cold,
			   soil_tag_cold, vegn_tag_cold, NULL, 0);

	/* glac */
	if(has_glac)
	  get_land_tile_info(fid_land_cold, GLAC_NAME, NULL, nidx_cold, idx_cold, frac_cold,
			     nx_dst, ny_dst, 1, isc_dst, iec_dst, glac_count_cold, glac_frac_cold,
			     glac_tag_cold, NULL, NULL, 0);

	/* lake */
	if(has_lake)
	  get_land_tile_info(fid_land_cold, LAKE_NAME, NULL, nidx_cold, idx_cold, frac_cold,
			     nx_dst, ny_dst, 1, isc_dst, iec_dst, lake_count_cold, lake_frac_cold,
			     lake_tag_cold, NULL, NULL, 0);

	for(i=0; i<nxc_dst; i++) for(l=0; l<soil_count_cold[i]; l++) soil_frac_cold[i] += tmp_frac_cold[ntile_cold*i+l];
      }
      

      /* find nearest source grid for each destination grid for soil, lake, glac */
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
	  get_actual_file_name(nface_dst, face_dst, remap_file, file);
	  remap_file_exist = mpp_file_exist(file);
	  if( !remap_file_exist) write_remap_file = 1;
	}
	/* The following mpp_sync is necessary. It is possible that the root pe is in the process
	   of writing remap_file and other processors think the remap_file exist
	*/
        mpp_sync();
	if(remap_file_exist) { /* read from remap file */
	  size_t start[4], nread[4];
	  int fid, nidx, vid;

	  if(mpp_pe() == mpp_root_pe()) printf("Read remap information from remap_file \n");
	  for(l=0; l<4; l++) {
	    start[l] = 0; nread[l] =1;
	  }
	  start[0] =isc_dst; nread[0] = nxc_dst;
	  fid = mpp_open(file, MPP_READ);
	  nlon = mpp_get_dimlen(fid, LON_NAME);
	  nlat = mpp_get_dimlen(fid, LAT_NAME);
	  if(nlon != nx_dst) mpp_error("remap_land: mismatch of dimension length of LON between destination grid and remap_file");
	  if(nlat != ny_dst) mpp_error("remap_land: mismatch of dimension length of LAT between destination grid and remap_file");

	  vid = mpp_get_varid(fid, "remap_face_soil");
	  mpp_get_var_value_block(fid, vid, start, nread, face_map_soil);
	  vid = mpp_get_varid(fid, "remap_index_soil");
	  mpp_get_var_value_block(fid, vid, start, nread, idx_map_soil);

	  if(has_glac) {
	    vid = mpp_get_varid(fid, "remap_face_glac");
	    mpp_get_var_value_block(fid, vid, start, nread, face_map_glac);
	    vid = mpp_get_varid(fid, "remap_index_glac");
	    mpp_get_var_value_block(fid, vid, start, nread, idx_map_glac);
	  }

	  if(has_lake) {
	    vid = mpp_get_varid(fid, "remap_face_lake");
	    mpp_get_var_value_block(fid, vid, start, nread, face_map_lake);
	    vid = mpp_get_varid(fid, "remap_index_lake");
	    mpp_get_var_value_block(fid, vid, start, nread, idx_map_lake);
	  }
	  mpp_close(fid);
	}
	else {
	    /* compute remap info for soil */
	    full_search_nearest(nface_src, nx_src*ny_src, x_src, y_src, soil_count_src, nxc_dst,
				x_dst, y_dst, soil_count_cold, idx_map_soil, face_map_soil);
	    /* compute remap info for glac */
	    if(has_glac)
	      full_search_nearest(nface_src, nx_src*ny_src, x_src, y_src, glac_count_src, nxc_dst,
				  x_dst, y_dst, glac_count_cold, idx_map_glac, face_map_glac);
	    /* compute remap info for lake */
	    if(has_lake)
	      full_search_nearest(nface_src, nx_src*ny_src, x_src, y_src, lake_count_src, nxc_dst,
				  x_dst, y_dst, lake_count_cold, idx_map_lake, face_map_lake);
	}

	if(write_remap_file) { /* write out restart file */
	  int *gdata=NULL;
	  int fid;
	  int dim_lon, dim_lat, dim_npts;
	  int id_face_soil, id_index_soil;
	  int id_face_glac, id_index_glac;
	  int id_face_lake, id_index_lake;

	  if(mpp_pe()==mpp_root_pe())printf("write remap information to remap_file\n");
	  fid = mpp_open(file, MPP_WRITE);
      	  dim_lon = mpp_def_dim(fid, LON_NAME, nx_dst);
      	  dim_lat = mpp_def_dim(fid, LAT_NAME, ny_dst);
	  dim_npts = mpp_def_dim(fid, "num_points", nx_dst*ny_dst);
	  gdata = (int *)malloc(nx_dst*ny_dst*sizeof(int));	  
	    
	  id_face_soil =  mpp_def_var(fid, "remap_face_soil", NC_INT, 1, &dim_npts, 1,
				      "standard_name", "soil remap face");
	  id_index_soil =  mpp_def_var(fid, "remap_index_soil", NC_INT, 1, &dim_npts, 1,
				       "standard_name", "soil remap index");
	  if(has_glac) {
	    id_face_glac =  mpp_def_var(fid, "remap_face_glac", NC_INT, 1, &dim_npts, 1,
					"standard_name", "glac remap face");
	    id_index_glac =  mpp_def_var(fid, "remap_index_glac", NC_INT, 1, &dim_npts, 1,
					 "standard_name", "glac remap index");
	  }
	  if(has_lake) {
	    id_face_lake =  mpp_def_var(fid, "remap_face_lake", NC_INT, 1, &dim_npts, 1,
					"standard_name", "lake remap face");
	    id_index_lake =  mpp_def_var(fid, "remap_index_lake", NC_INT, 1, &dim_npts, 1,
					 "standard_name", "lake remap index");	  
	  }
	  mpp_def_global_att(fid, "history", history);
	  mpp_end_def(fid);

	  mpp_gather_field_int_root(nxc_dst, face_map_soil, gdata);
	  mpp_put_var_value(fid, id_face_soil, gdata);
	  mpp_gather_field_int_root(nxc_dst, idx_map_soil, gdata);
	  mpp_put_var_value(fid, id_index_soil, gdata);

	  if(has_glac) {
	    mpp_gather_field_int_root(nxc_dst, face_map_glac, gdata);
	    mpp_put_var_value(fid, id_face_glac, gdata);
	    mpp_gather_field_int_root(nxc_dst, idx_map_glac, gdata);
	    mpp_put_var_value(fid, id_index_glac, gdata);
	  }

	  if(has_lake) {
	    mpp_gather_field_int_root(nxc_dst, face_map_lake, gdata);
	    mpp_put_var_value(fid, id_face_lake, gdata);
	    mpp_gather_field_int_root(nxc_dst, idx_map_lake, gdata);
	    mpp_put_var_value(fid, id_index_lake, gdata);
	  }
	  mpp_close(fid);
	  free(gdata);
	}
      }

      /* get soil, glac, lake and vegn data on the destination grid when file_type is land */
      {
	nidx_dst = 0; 
	for(i=0; i<nxc_dst; i++) {
	  int n, idx, p, count;

	  if( soil_count_cold[i] > 0 ) {
	    double totfrac;
	  
	    n = face_map_soil[i];
	    if(n<0) {
	      printf("soil_count_cold=%d, face_map_soil=%d, pe=%d, i=%d\n", soil_count_cold[i], n, mpp_pe(), i);
	      mpp_error("remap_land: soil_count_cold >0 but face_map_soil<0");
	    }
	    idx = idx_map_soil[i];
	    p = n*nx_src*ny_src+idx;
	    count = soil_count_src[p];
	    if( filetype != GLACTYPE && filetype != LAKETYPE ) {
	      pos = land_count_dst[i];
	      totfrac = 0;
	      for(l=0; l<count; l++) {
		soil_tag_dst[ntile_dst*i+pos] = soil_tag_src[ntile_src*p+l];
		vegn_tag_dst[ntile_dst*i+pos] = vegn_tag_src[ntile_src*p+l];
		land_frac_dst[ntile_dst*i+pos] = soil_frac_src[ntile_src*p+l]*soil_frac_cold[i];
		land_idx_map[ntile_dst*i+pos] = idx_soil_src[ntile_src*p+l];
		land_face_map[ntile_dst*i+pos] = n;
		totfrac += soil_frac_src[ntile_src*p+l];
		pos++;
	      }
	      pos = land_count_dst[i];
	      for(l=pos; l<pos+count; l++) land_frac_dst[ntile_dst*i+l] /= totfrac;
	      nidx_dst += count;
	    }
	    land_count_dst[i] += count;
	  }
	  if(has_glac) {
	    if( glac_count_cold[i] > 0 ) {
	      n = face_map_glac[i];
	      idx = idx_map_glac[i];
	      p = n*nx_src*ny_src+idx;
	      pos = land_count_dst[i];
	      count = glac_count_src[p];
	      if( count != 1) mpp_error("remap_land: glac_count_src should be 1");
	      if( filetype != SOILTYPE && filetype != VEGNTYPE  && filetype != LAKETYPE ){
		glac_tag_dst[ntile_dst*i+pos] = glac_tag_src[p];
		land_frac_dst[ntile_dst*i+pos] = glac_frac_cold[i];  /* preserve fraction in cold restart file */
		land_idx_map[ntile_dst*i+pos] = idx_glac_src[p];
		land_face_map[ntile_dst*i+pos] = n;
		nidx_dst++;
	      }
	      land_count_dst[i]++;
	    }
	  }
	  if(has_lake) {
	    if( lake_count_cold[i] > 0 ) {
	      n = face_map_lake[i];
	      idx = idx_map_lake[i];
	      p = n*nx_src*ny_src+idx;
	      count = lake_count_src[p];
	      pos = land_count_dst[i];
	      if( count != 1) mpp_error("remap_land: lake_count_src should be 1");
	      if( filetype != SOILTYPE && filetype != VEGNTYPE && filetype != GLACTYPE ) {
		lake_tag_dst[ntile_dst*i+pos] = lake_tag_src[p];
		land_frac_dst[ntile_dst*i+pos] = lake_frac_cold[i]; /* preserve fraction in cold restart file */
		land_idx_map[ntile_dst*i+pos] = idx_lake_src[p];;
		land_face_map[ntile_dst*i+pos] = n;
		nidx_dst++;
	      }
	      land_count_dst[i]++;
	    }	
	  }
	}
      }
      nidx_dst_global = nidx_dst;
      mpp_sum_int(1, &nidx_dst_global);
      /* compute tile_index */
      for(i=0; i<ntile_dst*nxc_dst; i++) idx_dst[i] = MPP_FILL_INT;
      
      for(n=0; n<ntile_dst; n++) {
	for(i=0; i<nxc_dst; i++) {
	  if( land_idx_map[ntile_dst*i+n] > -1 ) {
	    idx_dst[ntile_dst*i+n] = n*nx_dst*ny_dst+i+isc_dst;
	  }
	}      
      }
      
      /* define the metadata for dst_restart_file */
      {
	int dim_time, dim_cohort_index, dim_lat, dim_lon;
	int dim_tile_index, dim_cohort, dim_tile, dim_z;
      
	dim_lon          = mpp_def_dim(fid_dst, LON_NAME, nx_dst);
	dim_lat          = mpp_def_dim(fid_dst, LAT_NAME, ny_dst);
	dim_tile         = mpp_def_dim(fid_dst, TILE_NAME, ntile_dst);   
	dim_tile_index   = mpp_def_dim(fid_dst, TILE_INDEX_NAME, max(nidx_dst_global,1) );
	if(zaxis_exist) dim_z            = mpp_def_dim(fid_dst, LEVEL_NAME, nz);
	if(filetype==VEGNTYPE) {
	  dim_cohort       = mpp_def_dim(fid_dst, COHORT_NAME, ncohort);
	  dim_cohort_index = mpp_def_dim(fid_dst, COHORT_INDEX_NAME, max(nidx_dst_global,1)); 
	}
	if(time_exist)dim_time = mpp_def_dim(fid_dst, TIMENAME, NC_UNLIMITED);
	
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
	    else if( !strcmp(dimname, LEVEL_NAME ))
	      dims[m] = dim_z;
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
	if(src_has_tile ==0 && cold_has_tile==1) {
	  int vid2, vid1;
	  vid2 = mpp_def_var(fid_dst, TILE_NAME, MPP_INT, 1, &dim_tile, 0);
	  vid1 = mpp_get_varid(fid_cold, TILE_NAME);
	  mpp_copy_var_att(fid_cold, vid1, fid_dst, vid2); 
	}
	if(src_has_cohort ==0 && cold_has_cohort==1) {
	  int vid2, vid1;
	  vid2 = mpp_def_var(fid_dst, COHORT_NAME, MPP_INT, 1, &dim_cohort, 0);
	  vid1 = mpp_get_varid(fid_cold, COHORT_NAME);
	  mpp_copy_var_att(fid_cold, vid1, fid_dst, vid2); 
	}	
      }

      mpp_def_global_att(fid_dst, "history", history);
      mpp_end_def(fid_dst);

      
      /*-------------------------------------------------------------------------------
	Remap the data and write out to dst_restart_file
	It is assumed all the fields need to be remapped will have the first dimension
	name "tile_index" or "cohort_index" ( excludes "tile_index" or "cohort_index")
	-----------------------------------------------------------------------------*/

      /* First find number of varaibles in each file */

      nvar_cold = mpp_get_nvars(fid_cold);

      /* make sure number of variables match between src_restart_file and dst_cold_restart */
      
      if( !time_exist) { 
	int nvar1, nvar2;
	nvar1 = nvar_src;
	nvar2 = nvar_cold;

        if(src_has_tile==0 && cold_has_tile==1)
	  nvar2--;
	else if(src_has_tile==1 && cold_has_tile==0)
	  nvar1--;

        if(src_has_cohort==0 && cold_has_cohort==1)
	  nvar2--;
	else if(src_has_cohort==1 && cold_has_cohort==0)
	  nvar1--;

        if(src_has_theta_av==0 && cold_has_theta_av==1)
          nvar2--;
        else if(src_has_theta_av==1 && cold_has_theta_av==0)
	  nvar1--;

        if(src_has_theta_av_phen==1 && cold_has_theta_av_phen==0)
          nvar1--;
        else if(src_has_theta_av_phen==0 && cold_has_theta_av_phen==1)
	  nvar2--;

        if(src_has_theta_av_fire==1 && cold_has_theta_av_fire==0)
          nvar1--;
        else if(src_has_theta_av_fire==0 && cold_has_theta_av_fire==1)
	  nvar2--;

        if(src_has_psist_av==1 && cold_has_psist_av==0)
          nvar1--;
        else if(src_has_psist_av==0 && cold_has_psist_av==1)
	  nvar2--;
	
	if(nvar1 != nvar2)
	    mpp_error("remap_land: mismatch of nvar between src_restart_file and dst_cold_restart");
      }
	  
      rdata_global = (double *)malloc(nidx_dst_global*sizeof(double));
      idata_global = (int *)malloc(nidx_dst_global*sizeof(int));

      /* loop through each time level */
      for(t=0; t<ntime; t++) {
	for(l=0; l<nvar_src; l++) {
	  char   varname[128], dimname[128];
	  int    nz_dst;
	  int    m;

	  if( !has_taxis[l] && t>0 ) continue;      
	  mpp_get_varname(fid_src[0], l, varname);
	  if(nidx_dst_global==0) {
	    if(strcmp(varname, LON_NAME) && strcmp(varname, LAT_NAME) &&
	       strcmp(varname, LEVEL_NAME) && strcmp(varname, TILE_NAME) &&
	       strcmp(varname, COHORT_NAME)) continue;
	  }

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
	    int ndim_cold;
	    if( (strcmp(varname, THETA_AV_NAME) || cold_has_theta_av) &&
	        (strcmp(varname, TILE_NAME)     || cold_has_tile    ) &&
	        (strcmp(varname, COHORT_NAME)   || cold_has_cohort  ) &&
                (strcmp(varname, THETA_AV_PHEN_NAME) || cold_has_theta_av_phen) &&
                (strcmp(varname, THETA_AV_FIRE_NAME) || cold_has_theta_av_fire) &&
                (strcmp(varname, PSIST_AV_NAME) || cold_has_psist_av) ) {
	      vid_cold = mpp_get_varid(fid_cold, varname);
	      ndim_cold = mpp_get_var_ndim(fid_cold, vid_cold);
	      if( ndim_src[l] != ndim_cold)
		mpp_error("remap_land: number of dimensions for the field in dst_read_restart "
			  "does not match that in src_restart_file when time does not exist");
	    }
	  }

	  nz_dst = 1;
	  if(!time_exist) {  
	    if(ndim_src[l] == 2) {
	      mpp_get_var_dimname(fid_cold, vid_cold, 0, dimname);
	      nz_dst = mpp_get_dimlen(fid_cold, dimname);
	    }
	    if(nz_dst != nz_src[l])
	      mpp_error("remap_land: the src_restart and dst_cold_restart have the different dimension length.");
	  }
	  
	  pos = 0;
	  for(m=0; m<4; m++) {
	    start[m]=0; nread[m] = 1;
	  }
	  
	  /* write out lon, lat, tile, tile_index */
	  if( !strcmp(varname, N_ACCUM_NAME) || !strcmp(varname, NMN_ACM_NAME) ) { /* for n_accum, nmn_acm */
	    int fid, vid, scalar_data;
	    
	    fid = fid_src[0];
	    if(nface_src == nface_dst) fid = fid_src[face_dst];
	    vid = mpp_get_varid(fid, varname);
	    mpp_get_var_value(fid, vid, &scalar_data);
	    mpp_put_var_value(fid_dst, vid_dst, &scalar_data);
	  }
 	  else if( !strcmp(varname, LON_NAME) )  
	     mpp_put_var_value(fid_dst, vid_dst, lon_axis_dst ); 
 	  else if( !strcmp(varname, LAT_NAME) )
	     mpp_put_var_value(fid_dst, vid_dst, lat_axis_dst );
	  else if( !strcmp(varname, TILE_NAME) )
	     mpp_put_var_value(fid_dst, vid_dst, tile_axis_data );
	  else if( !strcmp(varname, LEVEL_NAME) )
	     mpp_put_var_value(fid_dst, vid_dst, z_axis_data );
	  else if( !strcmp(varname, COHORT_NAME) )
	     mpp_put_var_value(fid_dst, vid_dst, cohort_data );	  
	  else if( !strcmp(varname, TILE_INDEX_NAME) || !strcmp(varname, COHORT_INDEX_NAME)  ) {
	    compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idx_dst, idata_global, use_all_tile);
	    mpp_put_var_value(fid_dst, vid_dst, idata_global);
	  }
	  else if( !strcmp(varname, FRAC_NAME) ) {
	    compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, land_frac_dst, rdata_global, use_all_tile);
 	    mpp_put_var_value(fid_dst, vid_dst, rdata_global);
	  }
 	  else if( !strcmp(varname, GLAC_NAME) ) {
	    if(has_glac) {
	       compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, glac_tag_dst, idata_global, use_all_tile);
 	       mpp_put_var_value(fid_dst, vid_dst, idata_global);
	    }
	  }
 	  else if( !strcmp(varname, LAKE_NAME) ) {
	    if(has_lake) {
	      compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, lake_tag_dst, idata_global, use_all_tile);
 	      mpp_put_var_value(fid_dst, vid_dst, idata_global);
	    }
	  }
 	  else if( !strcmp(varname, SOIL_NAME) ) {
	    compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, soil_tag_dst, idata_global, use_all_tile);
 	    mpp_put_var_value(fid_dst, vid_dst, idata_global);
	  }
	  else if( !strcmp(varname, VEGN_NAME) ) {
	    compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, vegn_tag_dst, idata_global, use_all_tile);
 	    mpp_put_var_value(fid_dst, vid_dst, idata_global);
	  }
	  else { /* other fields, read source data and do remapping */
	    start_pos[0] = 0;
	    for(m=1; m<nface_src; m++) {
	      start_pos[m] = start_pos[m-1] + nidx_src[m-1];
	    }
 
	    for(k=0; k<nz_dst; k++) {
	      pos = 0;
	      /* read the source data */
	      for(m=0; m<nface_src; m++) {
		if(has_taxis[l]) {
		  start[0] = t; nread[1] = nidx_src[m];
		}
		else {
		  start[0] = k;
		  start[ndim_src[l]-1] = 0;
		  nread[ndim_src[l]-1] = nidx_src[m];
		}
		vid_src = mpp_get_varid(fid_src[m], varname);
		if(var_type[l] == MPP_INT)
		  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, idata_src+pos);
		else
		  mpp_get_var_value_block(fid_src[m], vid_src, start, nread, data_src+pos);
		pos += nidx_src[m];
	      }
	      for(m=0; m<4; m++) {
		start[m]=0; nwrite[m] = 1;
	      }
	      if(has_taxis[l]){
		start[0] = t; nwrite[1] = nidx_dst_global;
	      }
	      else {
		start[0] = k; nwrite[ndim_src[l]-1] = nidx_dst_global;
	      }
	      if(var_type[l] == MPP_INT) {
		for(m=0; m<ntile_dst*nxc_dst; m++) {
		  int face, lll;
		  if(land_idx_map[m] < 0) {
		    idata_dst[m] = MPP_FILL_INT;
		    continue;
		  }
		  face = land_face_map[m];
		  lll = start_pos[face]+land_idx_map[m];
		  idata_dst[m] = idata_src[lll];
		}
		compress_int_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, idata_dst, idata_global, use_all_tile);
		mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, idata_global);
	      }
	      else {
		for(m=0; m<ntile_dst*nxc_dst; m++) {
		  int face, lll;
		  if(land_idx_map[m] < 0) {
		    data_dst[m] = MPP_FILL_DOUBLE;
		    continue;
		  }
		  face = land_face_map[m];
		  lll = start_pos[face]+land_idx_map[m];
		  data_dst[m] = data_src[lll];
		}
		compress_double_data(ntile_dst, nxc_dst, nidx_dst, nidx_dst_global, land_count_dst, data_dst, rdata_global, use_all_tile);
		mpp_put_var_value_block(fid_dst, vid_dst, start, nwrite, rdata_global);
	      }
	    }
	  }
	}

	if(t==0 && src_has_tile ==0 && cold_has_tile==1) {
	  int i;
	  int *tile_data=NULL;
	  tile_data = (int *)malloc(ntile_dst*sizeof(int));
	  for(i=0; i<ntile_dst; i++) tile_data[i] = i+1;
	  vid_dst = mpp_get_varid(fid_dst, TILE_NAME);
	  mpp_put_var_value(fid_dst, vid_dst, tile_data);
	  free(tile_data);
	}
	if(t==0 && src_has_cohort ==0 && cold_has_cohort==1) {
	  int i;
	  int *cohort_data=NULL;
	  cohort_data = (int *)malloc(ncohort*sizeof(int));
	  for(i=0; i<ncohort; i++) cohort_data[i] = i+1;
	  vid_dst = mpp_get_varid(fid_dst, COHORT_NAME);
	  mpp_put_var_value(fid_dst, vid_dst, cohort_data);
	  free(cohort_data);
	}	
      }
      mpp_close(fid_dst);
      mpp_close(fid_cold);

      free(frac_cold);
      free(rdata_global);
      free(idata_global);
      
      if(mpp_pe()==mpp_root_pe()) printf("NOTE from remap_land: %s is created\n", file_dst );
      if(print_memory) {
	int n;
	char mesg[128];
	sprintf(mesg, "End of loop face_dst = %d", face_dst);
        print_mem_usage(mesg);
      }
    }


    free(lon_axis_dst);
    free(lat_axis_dst);
    free(x_tmp);
    free(y_tmp);
    free(start_pos);
    free(idata_src);
    free(data_src);
    free(data_dst);
    free(idata_dst);
    free(x_dst);
    free(y_dst);
    free(land_idx_map);
    free(land_face_map);
    free(land_count_dst);
    free(idx_dst);
    free(soil_tag_dst);
    free(vegn_tag_dst);
    free(land_frac_dst);
    free(soil_count_cold);
    free(soil_frac_cold);
    free(tmp_frac_cold);
    free(soil_tag_cold);
    free(vegn_tag_cold);
    free(idx_map_soil);
    free(face_map_soil);
    if(has_glac) {
      free(glac_tag_dst);
      free(glac_count_cold);
      free(glac_frac_cold);
      free(glac_tag_cold);
      free(idx_map_glac);
      free(face_map_glac);
    }
    if(has_lake){
      free(lake_tag_dst);
      free(lake_count_cold);
      free(lake_frac_cold);
      free(lake_tag_cold);
      free(idx_map_lake);
      free(face_map_lake);
    }

    mpp_delete_domain2d(&Dom_dst);
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

/********************************************************************************
void get_land_tile_info( )
********************************************************************************/
 void get_land_tile_info(int fid, const char *name1, const char *name2, int nidx, const int *idx_in, const double *frac_in,
			 int nx, int ny, int ntile, int isc, int iec, int *count, double *frac, int *tag1, int *tag2, int *idx, int all_tile)
 {
  
   int vid, l, i, j, k, p, nxc, pos; 
   int *tmp1=NULL;
   int *tmp2=NULL;

   nxc = iec - isc + 1;
   
   tmp1 = (int *)malloc(nidx*sizeof(int));
   
    vid = mpp_get_varid(fid, name1);
    mpp_get_var_value(fid, vid, tmp1);
    if(tag2) {
      if(!name2) mpp_error("remap_land: name2 can not be NULL when tag2 is not NULL");
      tmp2 = (int *)malloc(nidx*sizeof(int));
      vid = mpp_get_varid(fid, name2);
      mpp_get_var_value(fid, vid, tmp2);
    }
      
    /* set count to 0 */
    for(i=0; i<nxc; i++) count[i] = 0;
    for(i=0; i<ntile*nxc; i++) {
      frac[i] = MPP_FILL_DOUBLE;
      tag1[i] = MPP_FILL_INT;
      if(tag2) tag2[i] = MPP_FILL_INT;
      if(idx) idx[i] = MPP_FILL_INT;
    }
    pos = 0;
    for(l=0; l<nidx; l++) {
      if(tmp1[l] != MPP_FILL_INT ) {
	i = idx_in[l]%nx;
	k = idx_in[l]/nx;
	j = k%ny;
	p = j*nx+i;
	if(p<isc || p >iec) continue;
	p = p - isc;
	if(count[p] > ntile) mpp_error("remap_land: number of tiles is greater than allowed ntiles on one grid cell");
        frac[ntile*p+count[p]] = frac_in[l];
	tag1[ntile*p+count[p]] = tmp1[l];
	if(tag2) tag2[ntile*p+count[p]] = tmp2[l];
	if(idx) {
	  if(all_tile)
	    idx[ntile*p+count[p]] = l;
	  else
	    idx[ntile*p+count[p]] = pos;
	}
	/*	if(p==192) printf("idx=%d,l=%d,pos=%d,count=%d\n", idx[ntile*p+count[p]], l, pos, count[p]); */
	pos++;
	count[p]++;
      }
    }
       
    free(tmp1);
    if(tmp2) free(tmp2);
    

 } 
  

/********************************************************************
 void full_search_nearest
 search the nearest point from the first of source grid to the last.
********************************************************************/
void full_search_nearest(int nface_src, int npts_src, const double *lon_src, const double *lat_src,
			 const int *mask_src, int npts_dst, const double *lon_dst, const double *lat_dst,
			 const int *mask_dst, int *idx_map, int *face_map)
{
  int i_dst, i_src, l, face_cur;
  int pos, m, ind_cur;
  double d_cur,d;
  double p1[2], p2[2];  

  for(i_dst=0; i_dst<npts_dst; i_dst++) {
    if( mask_dst[i_dst] == 0 ){
      face_map[i_dst] = -1;
      idx_map[i_dst] = -1;
      continue;
    }
    d_cur = -1;
    pos = 0;
    p1[0]= lon_dst[i_dst];
    p1[1]= lat_dst[i_dst];
    
    for(m=0; m<nface_src; m++) {
      for(i_src=0; i_src<npts_src; i_src++) {
	l = m*npts_src + i_src;
	/*
	  if(face_cur==2 && i_dst == 16*192+168 && m==2 && (i_src==185 ||i_src==186 ||i_src==187||
							  i_src==233 ||i_src==234 ||i_src==235)) {
	     printf("at i_src=%d, mask_src=%d\n", i_src, mask_src[l]);
	}
	*/
	if( mask_src[l] == 0 ) continue;
        p2[0] = lon_src[l];
        p2[1] = lat_src[l];
	d = great_circle_distance(p1, p2);
	/*
if(face_cur==2 && i_dst == 16*192+168 && m==2 && (i_src==185 ||i_src==186 ||i_src==187||
							  i_src==233 ||i_src==234 ||i_src==235)) {
  printf("grid=%g,%g, d=%g,d_cur=%g,ind_cur=%d\n", p2[0]*R2D, p2[1]*R2D, d, d_cur, ind_cur);
 }
	*/
	if( d_cur < 0 || d<d_cur) {
	  ind_cur  = i_src;
	  face_cur = m;
	  d_cur    = d;
	}
      }
    }
    if(d_cur < 0)
      mpp_error("remap_land(full_search_nearest): no nearest point is found");
    else {
      /*
      if(face_cur==2 && i_dst == 16*192+168) printf("at (169,17), i_dst=%d\n", ind_cur);
      if(face_cur==2 && i_dst == 16*192+169) printf("at (170,17), i_dst=%d\n", ind_cur);
      if(face_cur==2 && i_dst == 16*192+168) printf("at (169,17), grid=%g,%g\n", p1[0]*R2D,p1[1]*R2D);
      */
      face_map[i_dst] = face_cur;
      idx_map[i_dst] = ind_cur;
    }
  }

} ; /* full_search_nearest */

    
/*-------------------------------------------------------------------------
  void compress_double_data ()
  get global compressed data
  ------------------------------------------------------------------------*/
void compress_double_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
			  const double *data, double *data_global, int all_tile )
{
  int pos1, pos2, pos, n, i, count;
  double *data_local=NULL;
  
  if(nidx>0) data_local = (double *)malloc(nidx*sizeof(double));
  for(i=0; i<nidx_global; i++) data_global[i] = MPP_FILL_DOUBLE; 

  pos1 = 0;
  pos2 = 0;
  for(n=0; n<ntile; n++) {
    pos = pos1;
    for(i=0; i<npts; i++) {
      if( n<land_count[i] ) {
	if(data[ntile*i+n] != MPP_FILL_DOUBLE) {
	  data_local[pos] = data[ntile*i+n];
	  if(!all_tile) pos++;
	}
	if(all_tile) pos++;
      }
    }
    count = pos - pos1;
    mpp_gather_field_double_root(count, data_local+pos1, data_global+pos2);
    pos1 = pos;
    mpp_sum_int(1, &count);
    pos2 += count;
  }
  
  if(nidx>0) free(data_local);	    

}

/*-------------------------------------------------------------------------
  void compress_int_data ()
  get global compressed data
  ------------------------------------------------------------------------*/
void compress_int_data(int ntile, int npts, int nidx, int nidx_global, const int *land_count,
		       const int *data, int *data_global, int all_tile )
{
  int pos1, pos2, pos, n, i, count;
  int *data_local=NULL;
  
  if(nidx>0) data_local = (int *)malloc(nidx*sizeof(int));
  for(i=0; i<nidx; i++) data_local[i] = MPP_FILL_INT;
  
  pos1 = 0;
  pos2 = 0;
  for(n=0; n<ntile; n++) {
    pos = pos1;
    for(i=0; i<npts; i++) {
      if( n<land_count[i] ) {
	if(data[ntile*i+n] != MPP_FILL_INT) {
	  data_local[pos] = data[ntile*i+n];
	  if(!all_tile) pos++;
	}
	if(all_tile) pos++;
      }
    }
    count = pos - pos1;
    mpp_gather_field_int_root(count, data_local+pos1, data_global+pos2);
    
    pos1 = pos;
    mpp_sum_int(1, &count);
    pos2 += count;
  }
  if(nidx>0) free(data_local);	    


}

  
