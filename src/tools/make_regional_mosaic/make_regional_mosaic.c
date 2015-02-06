/*
  Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
  This program is distributed under the terms of the GNU General Public
  License. See the file COPYING contained in this directory

  This program generates various types of horizontal grids in netCDF file format

 AUTHOR: Zhi Liang (Zhi.Liang)
          NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <netcdf.h>
#include "mpp.h"
#include "mpp_domain.h"
#include "constant.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "mosaic_util.h"

char *usage[] = {
  "",
  "                                                                                    ",
  "                                                                                    ",
  "                    Usage of make_regional_mosaic                                   ",
  "                                                                                    ",
  "   make_regional_mosaic --global_mosaic global_mosaic --regional_file regional_file ",
  "                                                                                    ",
  "   NOTE: This program can generate horizontal grid and solo mosaic for a regional   ",
  "         output. The created grid and solo mosaic could be used to regrid regional  ",
  "         output data onto regular lat-lon grid.                                     ",
  "                                                                                    ",
  "   make_regional_mosaic take the following flags                                    ",
  "                                                                                    ",
  "   --global_mosaic global_mosaic Specify the mosaic file for the global grid.       ",
  "                                                                                    ",
  "   --regional_file regional_file Specify the regional model output file.            ",
  "                                                                                    ",
  "   Example                                                                          ",
  "                                                                                    ",
  "      make_regional_mosaic --global_mosaic C48_mosaic.nc --regional_file            ",
  "                           rregionatmos_4xdaily.tile1.nc                            ",
  "                                                                                    ",
  NULL };

char tagname[] = "$Name: tikal $";
const char xaxis_name[] = "grid_xt_sub01";
const char yaxis_name[] = "grid_yt_sub01";

int main(int argc, char* argv[])
{
  char *global_mosaic=NULL;
  char *regional_file=NULL;
  char grid_descriptor[128] = "";

  int ntiles, tile;
  int nx, ny, i_min, i_max, j_min, j_max;

  double *x=NULL, *y=NULL;
  char history[512], dir[512];
  char gridfile[512], tilename[32];
  int errflg, c, i;  
  int option_index;

  static struct option long_options[] = {
    {"global_mosaic", required_argument, NULL, 'a'},
    {"regional_file", required_argument, NULL, 'b'},
    {"grid_descriptor", required_argument, NULL, 'c'},
    {0, 0, 0, 0},
  };

  mpp_init(&argc, &argv);

  errflg = argc <3;

  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'a':
      global_mosaic = optarg;
      break;
    case 'b':
      regional_file = optarg;
      break;
    case 'c':
      strcpy(grid_descriptor, optarg);
      break;
    case '?':
      errflg++;
    }
  }

  if( !global_mosaic ) {
    printf("ERROR from make_regional_mosaic: global_mosaic is not specified");
    errflg++;
  }
  if( !regional_file ) {
    printf("ERROR from make_regional_mosaic: regional_file is not specified");
    errflg++;
  }
  
  if (errflg) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }
  
  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }


  /* read the regional_file to get the tile number and regional grid index in global mosaic */
  /* currently we are assuming the regional file is on one file */
  {
    int fid;
    get_file_path(global_mosaic, dir);
    fid = mpp_open(global_mosaic, MPP_READ);
    ntiles = mpp_get_dimlen(fid, "ntiles");
    mpp_close(fid);
  }
  
  /* figure out the tile number from the regional_file name */
  {
    char *pch=NULL;
    char str[32];
    int  n;
    tile = -1;
    for(n=0; n<ntiles; n++) {
      sprintf(str, ".tile%d.nc", n+1);
      pch = strstr (regional_file, str);
      
      if( pch ) {
	tile = n;
	break;
      }
      pch = NULL;
    }
    if(tile == -1) mpp_error("make_regional_mosaic: tile number not found in the regional_file name");
  }
    
  {
    int fid, vid;
    double *tmpx=NULL, *tmpy=NULL;
    int *indx=NULL, *indy=NULL;
    
    fid = mpp_open(regional_file, MPP_READ);
    nx = mpp_get_dimlen(fid, xaxis_name);
    ny = mpp_get_dimlen(fid, yaxis_name);

    indx = (int *)malloc(nx*sizeof(int));
    indy = (int *)malloc(ny*sizeof(int));
    tmpx = (double *)malloc(nx*sizeof(double));
    tmpy = (double *)malloc(ny*sizeof(double));
    vid = mpp_get_varid(fid, xaxis_name);
    mpp_get_var_value(fid, vid, tmpx);
    for(i=0; i<nx; i++) indx[i] = round_to_nearest_int(tmpx[i]); /* use round to avoid truncation issue */
    vid = mpp_get_varid(fid, yaxis_name);
    mpp_get_var_value(fid, vid, tmpy);
    for(i=0; i<ny; i++) indy[i] = round_to_nearest_int(tmpy[i]); /* use round to avoid truncation issue */
    free(tmpx);
    free(tmpy);
    mpp_close(fid);
    i_min = indx[0];
    i_max = indx[nx-1];
    j_min = indy[0];
    j_max = indy[ny-1];

    /* i_max-i_min+1 should equal nx and j_max-j_min+1 should equal ny */
    if( i_max-i_min+1 != nx ) mpp_error("make_regional_mosaic: i_max-i_min+1 != nx ");
    if( j_max-j_min+1 != ny ) mpp_error("make_regional_mosaic: j_max-j_min+1 != ny ");
    
  }
  
  {
    char filename[512], gridfile[512];
    size_t start[4], nread[4];
    int fid, vid, n;
    
    /* read the regional grid from global grid */
    for(n=0; n<4; n++) {
      start[n] = 0;
      nread[n] = 1;
    }

    fid = mpp_open(global_mosaic, MPP_READ);
    start[0] = tile; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    vid = mpp_get_varid(fid, "gridfiles");
    mpp_get_var_value_block(fid, vid, start, nread, filename);
    sprintf(gridfile, "%s/%s", dir, filename);
    mpp_close(fid);
    
    x = (double *)malloc((2*nx+1)*(2*ny+1)*sizeof(double));
    y = (double *)malloc((2*nx+1)*(2*ny+1)*sizeof(double));
  
    fid = mpp_open(gridfile, MPP_READ);
    start[0] = 2*j_min - 2; nread[0] = 2*ny+1;
    start[1] = 2*i_min - 2; nread[1] = 2*nx+1;
    vid = mpp_get_varid(fid, "x");
    mpp_get_var_value_block(fid, vid, start, nread, x);
    vid = mpp_get_varid(fid, "y");
    mpp_get_var_value_block(fid, vid, start, nread, y);
    mpp_close(fid);
  }

  /* write out grid file */
  {
    int dimlist[5], dims[2], id_tile, id_x, id_y;
    int fid;
    int m;
    size_t start[4], nwrite[4];
    
    sprintf(gridfile, "regional_grid.tile%d.nc", tile+1);
    fid = mpp_open(gridfile, MPP_WRITE);
    dimlist[0] = mpp_def_dim(fid, "string", STRING);
    dimlist[1] = mpp_def_dim(fid, "nx", 2*nx);
    dimlist[2] = mpp_def_dim(fid, "ny", 2*ny);
    dimlist[3] = mpp_def_dim(fid, "nxp", 2*nx+1);
    dimlist[4] = mpp_def_dim(fid, "nyp", 2*ny+1);

    id_tile = mpp_def_var(fid, "tile", MPP_CHAR, 1, dimlist, 1, "standard_name", "grid_tile_spec");
    dims[0] = dimlist[4]; dims[1] = dimlist[3];
    id_x = mpp_def_var(fid, "x", MPP_DOUBLE, 2, dims, 2, "standard_name", "geographic_longitude",
		       "units", "degree_east");
    id_y = mpp_def_var(fid, "y", MPP_DOUBLE, 2, dims, 2, "standard_name", "geographic_latitude",
		       "units", "degree_north");
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);

    mpp_end_def(fid);
    sprintf(tilename, "tile%d", tile+1);
    for(m=0; m<4; m++) { start[m] = 0; nwrite[m] = 0; }
    nwrite[0] = strlen(tilename);
    mpp_put_var_value_block(fid, id_tile, start, nwrite, tilename );
    mpp_put_var_value(fid, id_x, x);
    mpp_put_var_value(fid, id_y, y);
    mpp_close(fid);
  }

  /* write out the mosaic file */
  {
    int fid, dim_ntiles, dim_string, id_mosaic, id_gridtiles;
    int dim[2], id_gridfiles;
    size_t start[4], nwrite[4];
    char mosaicfile[] = "regional_mosaic.nc";
    char mosaic_name[] = "regional_mosaic";
    fid = mpp_open(mosaicfile, MPP_WRITE);
    ntiles = 1;
    dim_ntiles = mpp_def_dim(fid, "ntiles", ntiles);
    dim_string = mpp_def_dim(fid, "string", STRING);
    id_mosaic = mpp_def_var(fid, "mosaic", MPP_CHAR, 1, &dim_string, 4, "standard_name",
			    "grid_mosaic_spec", "children", "gridtiles", "contact_regions", "contacts",
			    "grid_descriptor", grid_descriptor);
    dim[0] = dim_ntiles; dim[1] = dim_string;
    id_gridfiles = mpp_def_var(fid, "gridfiles", MPP_CHAR, 2, dim, 0);
    id_gridtiles = mpp_def_var(fid, "gridtiles", MPP_CHAR, 2, dim, 0);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    mpp_end_def(fid);

        /* write out data */
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }

    nwrite[0] = strlen(mosaic_name);
    mpp_put_var_value_block(fid, id_mosaic, start, nwrite, mosaic_name);
    nwrite[0] = 1;
    nwrite[1] = strlen(tilename);
    mpp_put_var_value_block(fid, id_gridtiles, start, nwrite, tilename);
    nwrite[1] = strlen(gridfile);
    mpp_put_var_value_block(fid, id_gridfiles, start, nwrite, gridfile);
    mpp_close(fid);
  }

  printf("congradulation: You have successfully run make_regional_mosaic\n");

  return 0;
  
}; // end of main
