/*
  This program remaps runoff data from a spherical grid onto any grid
  (spherical or tripolar) using conservative scheme.

 AUTHOR: Zhi Liang (Zhi.Liang)
          NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  For the full text of the GNU General Public License,
  write to: Free Software Foundation, Inc.,
            675 Mass Ave, Cambridge, MA 02139, USA.  
-----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "mpp.h"
#include "mpp_io.h"
#include "mpp_domain.h"
#include "create_xgrid.h"
#include "tool_util.h"

#define  EPSLN10 (1.e-10)
#define D2R (M_PI/180)
#define MAXATT   4096

char *usage[] = {
  "",
  "  runoff_regrid --input_file input_file --input_fld_name input_fld_name        ",
  "                --output_mosaic output_mosaic --output_topog output_topog      ",
  "                [--output_file output_file [--output_fid_name output_fid_name] ",
  "                                                                               ",
  "runoff_regrid expects to read runoff data from a netcdf file, which is         ",
  "specfied by option --input_file. The name of the runoff field in input         ",
  "file is specified by option --input_fld_name. The output file is a netcdf      ",
  "file specified by opiton --output_fld_file. The name of the runoff field       ",
  "in output file is specified by option output_fld_name. The output grid         ",
  "is specified by optoin --output_mosaic. The output grid topography is          ",
  "specified by option --output_topog output_topog. This tool only supports       ",
  "mosaic version grid. If there is any land point has runoff data after          ",
  "remapping runoff data onto destination grid, the runoff value of that          ",
  "land point will be moved to the nearest ocean point.                           ",
  "                                                                               ",
  "runoff_regrid takes the following flags:                                       ",
  "                                                                               ",
  "REQUIRED:                                                                      ",
  "                                                                               ",
  " --input_file                      specify the input file name.                ",
  "                                                                               ",
  " --input_fld_name                  specify name of runoff field in input file. ",
  "                                                                               ",
  " --output_mosaic output_mosaic     specify the output mosaic information. This ",
  "                                   file contains list of tile files which      ",
  "                                   specify the grid information for each tile. ",
  "                                                                               ",
  " --output_topog output_topog       specify the file name that contains output  ",
  "                                   grid topography information.                ",
  "                                                                               ",
  "                                                                               ",
  "OPTIONAL FLAGS                                                                 ",
  "                                                                               ",
  " --output_file output_file         specify the output file name. The default   ",
  "                                   value is runff.nc.                          ",
  "                                                                               ",
  " --output_fid_name output_fid_name specify name of runoff field in output file,",
  "                                   The default value is runoff                 ",
  "                                                                               ",
  "                                                                               ",
  "example:                                                                       ",
  "                                                                               ",
  "  runoff_regrid --input_file /archive/z1l/tools/input/runoff.daitren.iaf.10FEB2011.nc ",
  "                --input_fld_name runoff --output_mosaic                        ",
  "                /archive/z1l/tools/input/cm2m_grid/cm2m_ocean_mosaic.nc        ",
  "                --output_topog /archive/z1l/tools/input/cm2m_grid/cm2m_ocean_topog.nc ",
  "                --output_file runoff.cm2m.nc                                   ",
  "                                                                               ",
  NULL};

char tagname[] = "$Name: tikal $";

typedef struct {
  int nx;
  int ny;
  int n_ext;
  double *xt;
  double *yt;
  double *xc;
  double *yc;  
  double *xt1d;
  double *yt1d;
  double *xc1d;
  double *yc1d;
  double *mask;
  double *area;
} Grid_type;

typedef struct {
  int nxgrid;
  int    *i_in;
  int    *j_in;
  int    *i_out;
  int    *j_out;
  double *xgrid_area;
  int    *imap;
  int    *jmap;
} Remap_type;

double distance(double lon1, double lat1, double lon2, double lat2);
void get_input_grid(const char *file, const char *field, Grid_type *grid);
void get_output_grid(const char *mosaic, const char *topog_file, double sea_level, Grid_type *grid);
void setup_remap(const Grid_type *grid_in, const Grid_type *grid_out, Remap_type *remap);
void process_data(const char *infile, const char *fld_name_in, const char *outfile, const char *fld_name_out,
		  const Grid_type *grid_in, const Grid_type *grid_out, const Remap_type *remap, const char *history);
void nearest(int nlon, int nlat, double *mask, const double *lon, const double *lat,
	     double plon, double plat, int *iout, int *jout);  
int main(int argc, char* argv[])
{
  unsigned int opcode = 0;
  int          option_index, c, i;
  char    *input_file=NULL;
  char    *input_fld_name=NULL;
  char    *output_file=NULL;
  char    *output_fld_name=NULL;
  char    *output_mosaic=NULL;           /* output mosaic file name */
  char    *output_topog=NULL;
  double  sea_level = 0;
  char    history[MAXATT];  
  char    default_output_file[] = "runoff.nc";
  char    default_fld_name[]    = "runoff";
  int errflg = (argc == 1);

  Grid_type grid_in;
  Grid_type grid_out;
  Remap_type remap;
  
  static struct option long_options[] = {
    {"input_file",        required_argument, NULL, 'a'},
    {"input_fld_name",    required_argument, NULL, 'b'},
    {"output_mosaic",     required_argument, NULL, 'c'},
    {"output_topog",      required_argument, NULL, 'd'},
    {"output_file",       required_argument, NULL, 'e'},
    {"output_fld_name",   required_argument, NULL, 'f'},
    {"sea_level",         required_argument, NULL, 'g'},
    {"help",              no_argument,       NULL, 'h'},
    {0, 0, 0, 0},
  };  

  mpp_init(NULL,NULL);
  /*  if(mpp_npes() > 1) mpp_error("runoff_regrid: the tool is supposed to run on single processor"); */
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'a':
      input_file  = optarg;
      break;
    case 'b':
      input_fld_name  = optarg;
      break;
    case 'c':
      output_mosaic  = optarg;
      break;
    case 'd':
      output_topog = optarg;
      break;
    case 'e':
      output_file = optarg;
      break;
    case 'f':
      output_fld_name = optarg;
      break;
    case 'g':
      sea_level= atof(optarg);
      break;
    case '?':
      errflg++;
      break;
    }
  }

  if (errflg) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }

  if(!input_file) mpp_error("runoff_regrid: input_file is not specified");
  if(!input_fld_name) mpp_error("runoff_regrid: input_fld_name is not specified");
  if(!output_mosaic) mpp_error("runoff_regrid: output_mosaic is not specified");
  if(!output_topog) mpp_error("runoff_regrid: output_topog is not specified");
  if(!output_file) output_file = default_output_file;
  if(!output_fld_name) output_fld_name = default_fld_name;

  /* define history to be the history in the grid file */
  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }
  
  /* get input grid */
  get_input_grid(input_file, input_fld_name, &grid_in);
  if(mpp_pe() == mpp_root_pe() )printf("\nNOTE from runoff_regrid: complete get_input_grid\n");
  /* get output grid */
  get_output_grid(output_mosaic, output_topog, sea_level, &grid_out);
  if(mpp_pe() == mpp_root_pe() ) printf("\nNOTE from runoff_regrid: complete get_output_grid\n");  
  /* computing remapping information */
  setup_remap(&grid_in, &grid_out, &remap);
  if(mpp_pe() == mpp_root_pe() )printf("\nNOTE from runoff_regrid: complete setup_remap\n");
  
  /* do the remapping and write out data */
  /* process data on the root pe */
  if(mpp_pe() == mpp_root_pe()) {
    process_data(input_file, input_fld_name, output_file, output_fld_name, &grid_in, &grid_out, &remap, history);
  }

  if(mpp_pe() == mpp_root_pe() ) printf("NOTE from runoff_regrid: succefully created runoff data %s\n", output_file);

  mpp_end();

  return 0;
}  
     
  

void get_input_grid(const char *file, const char *field, Grid_type *grid)
{
  int fid, vid, ndim;
  int nx, ny, nxp, nyp, i, j;
  char xname[128], yname[128];
  size_t start[4], nread[4];
  double *data=NULL;
  double missing_value;
  
  fid = mpp_open(file, MPP_READ);
  vid = mpp_get_varid(fid, field);
  ndim = mpp_get_var_ndim(fid, vid);
  /* ndim-1 will be longitude and ndim-2 will be latitude */
  if( ndim !=2 && ndim !=3 ) mpp_error("runoff_regrid: ndim should be 2 or 3");
  mpp_get_var_dimname(fid, vid, ndim-1, xname);
  mpp_get_var_dimname(fid, vid, ndim-2, yname);
  nx = mpp_get_dimlen(fid, xname);
  ny = mpp_get_dimlen(fid, yname);
  nxp = nx+1;
  nyp = ny+1;
  grid->nx = nx;
  grid->ny = ny;
  grid->xt1d = (double *)malloc(nx*sizeof(double));
  grid->yt1d = (double *)malloc(ny*sizeof(double));
  grid->xc1d = (double *)malloc(nxp*sizeof(double));
  grid->yc1d = (double *)malloc(nyp*sizeof(double));
  grid->xc   = (double *)malloc(nxp*nyp*sizeof(double));
  grid->yc   = (double *)malloc(nyp*nyp*sizeof(double));
  grid->mask = (double *)malloc(nx*ny*sizeof(double));
  grid->area = (double *)malloc(nx*ny*sizeof(double));
  
  vid = mpp_get_varid(fid, xname);
  mpp_get_var_value(fid, vid, grid->xt1d);
  vid = mpp_get_varid(fid, yname);
  mpp_get_var_value(fid, vid, grid->yt1d);

  for(i=1; i<nx; i++) grid->xc1d[i] = (grid->xt1d[i-1] + grid->xt1d[i])*0.5;
  grid->xc1d[0]  = 2*grid->xt1d[0]    - grid->xc1d[1];
  grid->xc1d[nx] = 2*grid->xt1d[nx-1] - grid->xc1d[nx-1];
  for(j=1; j<ny; j++) grid->yc1d[j] = (grid->yt1d[j-1] + grid->yt1d[j])*0.5;
  grid->yc1d[0]  = 2*grid->yt1d[0]    - grid->yc1d[1];
  grid->yc1d[ny] = 2*grid->yt1d[ny-1] - grid->yc1d[ny-1];

  /* convert to radians */
  for(i=0; i<nxp; i++) grid->xc1d[i] *= D2R;
  for(j=0; j<nyp; j++) grid->yc1d[j] *= D2R;
    
  for(j=0; j<nyp; j++) for(i=0; i<nxp; i++) {
    grid->xc[j*nxp+i] = grid->xc1d[i];
    grid->yc[j*nxp+i] = grid->yc1d[j];
  }

  get_grid_area(&nx, &ny, grid->xc, grid->yc, grid->area);
  
  /* get the mask */
  vid = mpp_get_varid(fid, field);
  data = (double *)malloc(nx*ny*sizeof(double));
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }  
  nread[ndim-1] = nx;
  nread[ndim-2] = ny;
  mpp_get_var_value_block(fid, vid, start, nread, data);
  missing_value = -1.e+20;
  if(mpp_var_att_exist(fid, vid, "missing_value"))
    mpp_get_var_att_double(fid, vid, "missing_value", &missing_value);
  for(i=0; i<nx*ny; i++){
    if(fabs(data[i]-missing_value) <= EPSLN10)
      grid->mask[i] = 0.0;
    else
      grid->mask[i] = 1.0;
  }

  free(data);
  
  mpp_close(fid);
  
}
   
void get_output_grid(const char *mosaic, const char *topog_file, double sea_level, Grid_type *grid)
{
  int fid, ntile, vid, i, j, ind, n_ext;
  int ni, nj, nip, njp, nx, ny, nxp, nyp, ny_old;
  char gridfile[256], filename[256], dir[256];
  double *tmp=NULL, *depth=NULL;  
  
   get_file_path(mosaic, dir);
   fid = mpp_open(mosaic, MPP_READ);
   ntile  = mpp_get_dimlen(fid, "ntiles");
   if(ntile > 1) mpp_error("runoff_regrid: ntile in mosaic file should be 1");
   vid = mpp_get_varid(fid, "gridfiles");
   mpp_get_var_value(fid, vid, filename);
   sprintf(gridfile, "%s/%s", dir, filename);
   mpp_close(fid);
   fid = mpp_open(gridfile, MPP_READ);
   ni = mpp_get_dimlen(fid, "nx");
   nj = mpp_get_dimlen(fid, "ny");
   if(ni%2 != 0 ) mpp_error("runoff_regrid: ouptut supergrid x-size can not be divided by 2");
   if(nj%2 != 0 ) mpp_error("runoff_regrid: output supergrid y-size can not be divided by 2");
   nip = ni+1;
   njp = nj+1;
   tmp = (double *)malloc(nip*njp*sizeof(double));
   vid = mpp_get_varid(fid, "y");
   mpp_get_var_value(fid, vid, tmp);
   /* we might need to add one extra raw to expand to south pole to ensure conservative */
   n_ext = 0;
   
   if(tmp[0] > -90 + EPSLN10 ) n_ext = 1;
   if(mpp_pe() == mpp_root_pe() ) printf("The south extension is %d\n", n_ext);
   nx  = ni/2;
   ny  = nj/2;
   ny_old = ny;
   ny += n_ext;
   nxp = nx+1;
   nyp = ny+1;   
   grid->nx = nx;
   grid->ny = ny;
   grid->n_ext = n_ext;
   grid->xc = (double *)malloc(nxp*nyp*sizeof(double));
   grid->yc = (double *)malloc(nxp*nyp*sizeof(double));
   grid->xt = (double *)malloc(nx*ny*sizeof(double));
   grid->yt = (double *)malloc(nx*ny*sizeof(double));
   grid->xt1d = (double *)malloc(nx*sizeof(double));
   grid->yt1d = (double *)malloc(ny_old*sizeof(double));
   grid->area = (double *)malloc(nx*ny*sizeof(double));
   
   for(j=n_ext; j<nyp; j++) for(i=0; i<nxp; i++) grid->yc[j*nxp+i] = tmp[2*(j-n_ext)*nip+2*i];
   for(j=n_ext; j<ny ; j++) for(i=0; i<nx ; i++) grid->yt[j*nx +i] = tmp[(2*(j-n_ext)+1)*nip+2*i+1];
   ind = nx/4;
   for(j=0; j<ny_old; j++) grid->yt1d[j] = tmp[(2*j+1)*nip+2*ind];   
   
   if(n_ext >0) {
     for(i=0; i<nxp; i++) grid->yc[i] = -90;
     for(i=0; i<nx; i++) grid->yt[i] = 0.5*(grid->yc[i] + grid->yc[nxp+i]);
   }
        

   vid = mpp_get_varid(fid, "x");
   mpp_get_var_value(fid, vid, tmp);
   for(j=n_ext; j<nyp; j++) for(i=0; i<nxp; i++) grid->xc[j*nxp+i] = tmp[2*(j-n_ext)*nip+2*i];
   for(j=n_ext; j<ny ; j++) for(i=0; i<nx ; i++) grid->xt[j*nx +i] = tmp[(2*(j-n_ext)+1)*nip+2*i+1];
   for(i=0; i<nx; i++) grid->xt1d[i] = tmp[nip+2*i+1];
   if(n_ext >0) {
     for(i=0; i<nxp; i++) grid->xc[i] = grid->xc[nxp+i];
     for(i=0; i<nx; i++) grid->xt[i] =  grid->xt[nx+i];
   }     
   
   free(tmp);
   mpp_close(fid);

   /* convert grid->xt and grid->yt from degree to radian */
   for(i=0; i<nx*ny; i++) {
     grid->xt[i] *= D2R;
     grid->yt[i] *= D2R;
   }
   for(i=0; i<nxp*nyp; i++) {
     grid->xc[i] *= D2R;
     grid->yc[i] *= D2R;
   }

   get_grid_area(&nx, &ny, grid->xc, grid->yc, grid->area);
   
   /* read the topography data to get the land sea mask */
   fid = mpp_open(topog_file, MPP_READ);
   ntile = 1;
   if(mpp_dim_exist(fid, "ntiles"))ntile = mpp_get_dimlen(fid, "ntiles");
   if(ntile > 1) mpp_error("runoff_regrid: ntile in topog file should be 1");
   nx = mpp_get_dimlen(fid, "nx");
   ny = mpp_get_dimlen(fid, "ny");
   if(nx != grid->nx || ny != grid->ny-n_ext)
     mpp_error("runoff_regrid: grid size mismatch between mosaic file and topog file");
   ny += n_ext;
   depth = (double *)malloc(nx*ny_old*sizeof(double));
   grid->mask = (double *)malloc(nx*ny*sizeof(double));
   vid = mpp_get_varid(fid, "depth");
   mpp_get_var_value(fid, vid, depth);
   mpp_close(fid);
   for(i=0; i<nx*ny; i++) grid->mask[i] = 0.0;
   for(j=n_ext; j<ny; j++) for(i=0; i<nx; i++) {
     if(depth[(j-n_ext)*nx+i] >sea_level) grid->mask[j*nx+i] = 1.0;
   }
   free(depth);			    


    
}

void setup_remap(const Grid_type *grid_in, const Grid_type *grid_out, Remap_type *remap)
{
  int nx_in, ny_in, nx_out, ny_out;
  int nxgrid, i, j, ii, jj, l, ll;
  int *i_in, *j_in, *i_out, *j_out;
  double *xgrid_area;
  int npes, layout[2], nxgrid_local;
  int isc, iec, jsc, jec, nxc, nyc;
  domain2D domain;
  double *xc=NULL, *yc=NULL;
  int *imap=NULL, *jmap=NULL;
  
  
  nx_in = grid_in ->nx;
  ny_in = grid_in ->ny;
  nx_out = grid_out ->nx;
  ny_out = grid_out ->ny;
  npes = mpp_npes();
  
  mpp_define_layout(nx_out, ny_out, npes, layout);
  mpp_define_domain2d(nx_out, ny_out, layout, 0, 0, &domain);

  mpp_get_compute_domain2d(domain, &isc, &iec, &jsc, &jec);
  nxc = iec-isc+1;
  nyc = jec-jsc+1;
  xc = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));
  yc = (double *)malloc((nxc+1)*(nyc+1)*sizeof(double));

  /* copy grid to local data */
  for(j=0; j<=nyc; j++) for(i=0; i<=nxc; i++) {
    jj = j+jsc;
    ii = i+isc;
    xc[j*(nxc+1)+i] = grid_out->xc[jj*(nx_out+1)+ii];  
    yc[j*(nxc+1)+i] = grid_out->yc[jj*(nx_out+1)+ii];
  }
    
  i_in       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  j_in       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  i_out      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  j_out      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  xgrid_area = (double *)malloc(MAXXGRID   * sizeof(double));

  nxgrid = create_xgrid_1dx2d_order1(&nx_in, &ny_in, &nxc, &nyc, grid_in->xc1d,
				     grid_in->yc1d,  xc,  yc,
				     grid_in->mask, i_in, j_in, i_out, j_out, xgrid_area);
  /* add isc and jsc to i_out and j_out */
  for(i=0; i<nxgrid; i++) {
     i_out[i] += isc;
     j_out[i] += jsc;
  }

  nxgrid_local = nxgrid;
  
  mpp_sum_int(1, &nxgrid);
  
  remap->nxgrid = nxgrid;
  remap->i_in = (int *)malloc(nxgrid*sizeof(int));
  remap->j_in = (int *)malloc(nxgrid*sizeof(int));  
  remap->i_out = (int *)malloc(nxgrid*sizeof(int));
  remap->j_out = (int *)malloc(nxgrid*sizeof(int));
  remap->xgrid_area = (double *)malloc(nxgrid*sizeof(double));

  mpp_gather_field_int_root( nxgrid_local, i_in, remap->i_in);
  mpp_gather_field_int_root( nxgrid_local, j_in, remap->j_in);
  mpp_gather_field_int_root( nxgrid_local, i_out, remap->i_out);
  mpp_gather_field_int_root( nxgrid_local, j_out, remap->j_out);
  mpp_gather_field_double_root( nxgrid_local, xgrid_area, remap->xgrid_area);

  
  /* compute the nearest ocean points for any land points */
  imap = (int *)malloc(nxc*nyc*sizeof(int));
  jmap = (int *)malloc(nxc*nyc*sizeof(int));
  remap->imap = (int *)malloc(nx_out*ny_out*sizeof(int));
  remap->jmap = (int *)malloc(nx_out*ny_out*sizeof(int));

  for(i=0; i<nxc*nyc; i++) {
    imap[i] = -1;
    jmap[i] = -1;
  }
  
  
  for(j=0; j<nyc; j++) {
    for(i=0; i<nxc; i++) {
      jj = j+jsc;
      ii = i+isc;
      ll = jj*nx_out+ii;
      l = j*nxc+i;
      if( grid_out->mask[ll] == 0) {
	nearest(nx_out, ny_out, grid_out->mask, grid_out->xt, grid_out->yt,
		grid_out->xt[ll], grid_out->yt[ll], &(imap[l]), &(jmap[l]) );
      }
    }
  }
  mpp_global_field_int(domain, nxc, nyc, imap, remap->imap);
  mpp_global_field_int(domain, nxc, nyc, jmap, remap->jmap);
  

  free(imap);
  free(jmap);
    
  free(xc);
  free(yc);
  free(i_in);
  free(j_in);
  free(i_out);
  free(j_out);
  free(xgrid_area);
  mpp_delete_domain2d(&domain);
  
}


void process_data(const char *infile, const char *fld_name_in, const char *outfile, const char *fld_name_out,
		  const Grid_type *grid_in, const Grid_type *grid_out, const Remap_type *remap, const char *history)
{
  int fid_in, vid_in, fid_out, ndim, ntime;
  int id_xt, id_yt, id_time, id_fld, id_time_in;
  int nx_in, ny_in, nx_out, ny_out;
  int i,j,n,l,i1,j1,i2,j2, n_ext;
  int iout, jout;
  int dims[3];
  size_t start[4], nwrite[4], nread[4];
  double xarea;
  double *time_value=NULL;
  double *data_in=NULL, *data_out=NULL;
  double runoff_total_in, runoff_total_out, rate_change;
  char timename[128];
  
  /* get the time information in the file */
  fid_in = mpp_open(infile, MPP_READ);
  vid_in = mpp_get_varid(fid_in, fld_name_in);
  ndim   = mpp_get_var_ndim(fid_in, vid_in);
  ntime = 1;
  if(ndim>2) {
    mpp_get_var_dimname(fid_in, vid_in, 0, timename);
    ntime = mpp_get_dimlen(fid_in, timename);
    time_value = (double *)malloc(ntime*sizeof(double));
    id_time_in = mpp_get_varid(fid_in, timename);
    mpp_get_var_value(fid_in, id_time_in, time_value);
  }
  
  /* set up the output metadata */
  nx_in = grid_in ->nx;
  ny_in = grid_in ->ny;
  nx_out = grid_out ->nx;
  ny_out = grid_out ->ny;
  n_ext = grid_out->n_ext;
  data_in = (double *)malloc(nx_in*ny_in*sizeof(double));
  data_out = (double *)malloc(nx_out*ny_out*sizeof(double));

  fid_out = mpp_open(outfile, MPP_WRITE);
  mpp_def_global_att(fid_out, "history", history);
  dims[ndim-1] = mpp_def_dim(fid_out, "grid_x_T", nx_out);
  dims[ndim-2] = mpp_def_dim(fid_out, "grid_y_T", ny_out-n_ext);
  if(ndim>2) { /* has time axis */
    dims[0] = mpp_def_dim(fid_out, "Time", NC_UNLIMITED);
  }

  id_xt = mpp_def_var(fid_out, "grid_x_T", NC_DOUBLE, 1, dims+(ndim-1), 0);
  mpp_def_var_att(fid_out, id_xt, "long_name", "Nominal Longitude of T-cell center");
  mpp_def_var_att(fid_out, id_xt, "cartesian_axis", "X");
  mpp_def_var_att(fid_out, id_xt, "units", "degree_east");

  id_yt = mpp_def_var(fid_out, "grid_y_T", NC_DOUBLE, 1, dims+(ndim-2), 0);
  mpp_def_var_att(fid_out, id_yt, "long_name", "Nominal Latitude of T-cell center");
  mpp_def_var_att(fid_out, id_yt, "cartesian_axis", "Y");
  mpp_def_var_att(fid_out, id_yt, "units", "degree_north");  

  if(ndim > 2) {
    id_time = mpp_def_var(fid_out, "Time", NC_DOUBLE, 1, dims, 0);
    mpp_copy_var_att(fid_in, id_time_in, fid_out, id_time);
  }
    
  id_fld = mpp_def_var(fid_out, fld_name_out, NC_DOUBLE, ndim, dims, 0);
  mpp_copy_var_att(fid_in, vid_in, fid_out, id_fld);
  mpp_end_def(fid_out);

  /* write out axis data */
  mpp_put_var_value(fid_out, id_xt, grid_out->xt1d);
  mpp_put_var_value(fid_out, id_yt, grid_out->yt1d);
  
  /* loop through ntime */
  for(n=0; n<ntime; n++) {
    
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1; nread[i] = 1;
    }        
    start[0] = n;
    if(ndim>2) {
      mpp_put_var_value_block(fid_out, id_time, start, nwrite, time_value+n);
    }
    nread[ndim-1] = nx_in;
    nread[ndim-2] = ny_in;
    mpp_get_var_value_block(fid_in, vid_in, start, nread, data_in);
    /* remap the data */
    for(i=0; i<nx_out*ny_out; i++) data_out[i] = 0;
    for(l=0; l<remap->nxgrid; l++) {
      i2   = remap->i_out[l];
      j2   = remap->j_out[l];
      i1   = remap->i_in [l];
      j1   = remap->j_in [l];    
      xarea = remap->xgrid_area [l];
      /* since input mask is already considered when computing exchange grid,
	 there is no need to consider missing value
      */
      data_out[j2*nx_out+i2] += data_in[j1*nx_in+i1]*xarea;
    }
    
    /* move the runoff to the nearest ocean points */
    for(j=0; j<ny_out; j++) for(i=0; i<nx_out; i++) {
      l = j*nx_out+i;
      if(data_out[l] > 0 && grid_out->mask[l] == 0) {
	iout = remap->imap[l];
	jout = remap->jmap[l];
        data_out[jout*nx_out+iout] += data_out[l];
	data_out[l] = 0;
      }
    }

    /* divided by ocn area */
    for(l=0; l<nx_out*ny_out; l++) {
      data_out[l] /= grid_out->area[l];
    }
    nwrite[ndim-1]=nx_out;
    nwrite[ndim-2]=ny_out-n_ext;
    if(n_ext>0) {
      double *tmp;
      tmp = (double *)malloc(nx_out*(ny_out-n_ext)*sizeof(double));
      for(j=0; j<ny_out-n_ext; j++) for(i=0; i<nx_out; i++) tmp[j*nx_out+i] = data_out[(j+n_ext)*nx_out+i];
      mpp_put_var_value_block(fid_out, id_fld, start, nwrite, tmp);
      free(tmp);
    }
    else
      mpp_put_var_value_block(fid_out, id_fld, start, nwrite, data_out);
    /* conservative check */
    runoff_total_in = 0;
    runoff_total_out = 0;
    for(i=0; i<nx_in*ny_in; i++) {
      if(grid_in->mask[i]>0)runoff_total_in += (data_in[i]*grid_in->area[i]);
    }
    for(i=0; i<nx_out*ny_out; i++) runoff_total_out += (data_out[i]*grid_out->area[i]);
    if(runoff_total_in > 0 && runoff_total_out > 0) 
      rate_change = (runoff_total_out-runoff_total_in)/runoff_total_in;
    else
      rate_change = 0.0;
    printf("At time step %d, runoff_total_in=%g, runoff_total_out=%g, change rate is =%g%%\n",
	   n, runoff_total_in, runoff_total_out, rate_change);
    
  }
  mpp_close(fid_out);
  mpp_close(fid_in);

  free(time_value);
  free(data_in);
  free(data_out);
  
}

/* calculate a distance in 3-d between points on unit square */
double distance(double lon1, double lat1, double lon2, double lat2)
{
  double dist, z1, z2, x1, x2, y1, y2;

  
  z1 = sin(lat1);
  z2 = sin(lat2);
  x1 = cos(lat1)*cos(lon1);
  x2 = cos(lat2)*cos(lon2);
  y1 = cos(lat1)*sin(lon1);
  y2 = cos(lat2)*sin(lon2);

  dist = sqrt(pow((z1-z2),2)+pow((x1-x2),2)+pow((y1-y2),2));
  return dist;
}
	  
/* find the nearest ocean point of a given land point. */

void nearest(int nlon, int nlat, double *mask, const double *lon, const double *lat,
	     double plon, double plat, int *iout, int *jout)
{
  int i,j, n;
  double r, r1;

  r = distance(0.0,M_PI/2,0.0,-M_PI/2);  /* some big value */
  for(j=0; j<nlat; j++) for(i=0; i<nlon; i++) {
    n = j*nlon+i;
    if(mask[n] == 0.0) continue;
    r1 = distance(plon,plat,lon[n],lat[n]);
    if ( r1 < r ) {
      *iout = i;
      *jout = j;
      r = r1;
    }
  }

}

/* do a radial search starting from (iref, jref). */
/* #define MAX_ITER 1000 */
/* void radial_search(int nlon, int nlat, double *mask, const double *lon, const double *lat, */
/* 		   double plon, double plat, int *iout, int *jout, int iref, int jref) */
/* { */

/*   iter = 0; */
/*   while (iter<MAX_ITER) { */
/*     iter++; */

/*     /\* left boundary *\/ */
/*     i_left = iref-n; */
/*     if (i_left < 0) /\* cyclic condition is assumed here *\/ */
/*       i_left = nlon + i_left; */
/* 	else */
/* 	  i_left = 1 */




/*   } */
  
