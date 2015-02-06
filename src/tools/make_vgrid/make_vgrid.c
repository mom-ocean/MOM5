#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "create_vgrid.h"
#include "mpp_io.h"
#include "mpp.h"
#include "tool_util.h"

#define MAXBOUNDS 100

char *usage[] = {
  "",
  "                                                                                 ",
  "                                                                                 ",
  "                    Usage of make_vgrid                                          ",
  "                                                                                 ",
  "   make_vgrid --nbnds nbnds --bnds z(1),...,z(nbnds) [--nz nz(1),...,nz(nbnds-1)]",
  "              [--dbnds dbnds(1), ...,dbnds(nbnds)] [--stretch stretch]           ",
  "              [--grid_name gridname]    [--center center]                        ",
  "                                                                                 ",
  "   make_vgrid is used to make a vertical grid for an FMS model.  It uses         ",
  "   a piecewise monotonic shape preserving cubic-spline algorithm to              ",
  "   calculate the grid cell location.  The output netcdf file contains            ",
  "   information on the  supergrid with grid size equal to the model grid size     ",
  "   multiplied by a refinement factor of two.                                     ",
  "                                                                                 ",
  "   Two algorithms are provided for creating vertical grid: Monotonic             ",
  "   cubic interpolation and legacy algorithm using cosine function to             ",
  "   configure grid cell location. The monotonic cubic interpolation is            ",
  "   developed by Russell Fiedler and the detail of the alrogithm is               ",
  "   explained in                                                                  ",
  "                                                                                 ",
  "   http://www.mathworks.com/moler/interp.pdf. More details of the                ",
  "   algorithm is available in                                                     ",
  "                                                                                 ",
  "   F. N. Fritsch and R. E. Carlson, Monotone                                     ",
  "   Piecewise Cubic Interpolation, SIAM Journal on Numerical Analysis, 17         ",
  "   (1980), pp. 238-246.                                                          ",
  "   D. Kahaner, C. Moler, and S. Nash, Numerical Methods and Software,            ",
  "   Prentice Hall, Englewood Cliffs, NJ, 1989.                                    ",
  "                                                                                 ",
  "   --nz needs to be specified to use the monotonic cubic algorithm.  The         ",
  "   legacy algorithm is developed by Ron Pacanowski, and this older               ",
  "   algorithm uses cosine function to do the interpolation. --dbnds needs         ",
  "   to be specified to use the legacy algorithm. Since monotonic cubic            ",
  "   interpolation is a higher order interpolation, it will produce                ",
  "   smoother grid distance. It is strongly suggested to use the monotonic         ",
  "   cubic interpolation by specifying argument --nz.                              ",
  "                                                                                 ",
  "   make_vgrid takes the following flags                                          ",
  "                                                                                 ",
  "   Required Flags:                                                               ",
  "                                                                                 ",
  "   --nbnds nbnds             Specify number of vertical regions for varying      ",
  "                             resolution.                                         ",
  "                                                                                 ", 
  "   --bnds z(1),.,z(nbnds)    Specify boundaries for defining vertical regions of ",
  "                             varying resolution.                                 ",
  "                                                                                 ",
  
  "   Optional Flags:                                                               ",
  "                                                                                 ",
  "   --nz nz(1),..,nz(nbnds-1) Number of supergrid points (double number of        ",
  "                             model grid points) for each vertical                ",
  "                             regions of varying resolution.                      ",
  "                                                                                 ",
  " --dbnds dbnds(1),.,dbnds(nbnds) nominal resolution of depth regions             ",
  "                                                                                 ",
  " --stretch stretch           stretching factor for last region.                  ",
  "                             Default is 1 (no stretching)                        ",
  "                                                                                 ",
  " --grid_name gridname        Specify the grid name. The output grid file name    ",
  "                             will be grid_name.nc. The default value is          ",
  "                             vertical_grid.                                      ",
  "                                                                                 ", 
  " --center   center           Specify the center location of grid. The valid      ",
  "                             entry will be 'none', 't_cell' or 'c_cell' with     ",
  "                             default value 'none'. The grid refinement is        ",
  "                             assumed to be 2 in x and y-direction when center    ",
  "                             is not 'none'. 'c_cell' should be used for the grid ",
  "                             used in MOM.                                        ",
  "                                                                                 ",
  "                                                                                 ",
  "   Example 1: Use Monotonic cubic spline interpolation                           ", 
  "                                                                                 ",
  "     make_vgrid --nbnds 3 --bnds 10,200,1000 --nz 10,20                          ",
  "       will make a grid file with 60 supergrid cells, 30 model grid cells        ",
  "       with refinement is 2.                                                     ",
  "                                                                                 ",
  " In order to generate regions of constant resolution the user should specify 3   ",
  " consecutive regions of constant resolution (only 2 for the first and last regions",
  " are required)                                                                   ",
  "                                                                                 ",
  " make_vgrid --nbnds 4 --bnds 0,180,200,5000 --nz 18,2,30 --center c_cell         ",
  " will have constant resolution in the top 200m                                   ",
  "                                                                                 ",
  " make_vgrid --nbnds 6 --bnds 0,100,120,200,220,5000 --nz 20,2,8,2,30  --center c_cell  ",
  " will have constant resolution from 120m-200m                                    ",
  "                                                                                 ",
  "                                                                                 ",
  "   Example 2: Use legacy algorithm                                               ",
  "                                                                                 ",
  "    make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286          ",
  "               --center c_cell                                                   ",
  "      Will make a vertical grid similar to GFDL/CM2/ocn w/ 100 supergrid cells.  ",
  "                                                                                 ",
  "",
  NULL };  

char grid_version[] = "0.2";
char tagname[] = "$Name: tikal $";

int main(int argc, char* argv[])
{
  int nbnds, n1, n2, n3, i, nk;
  int use_legacy;
  double bnds[MAXBOUNDS];
  double dbnds[MAXBOUNDS];
  int    nz[MAXBOUNDS-1];
  char gridname[128]= "vertical_grid";
  char filename[128];
  double *zeta;
  double stretch = 1;
  char center[32] = "none";
  char entry[512];
  char history[256];
  int errflg, c, option_index = 0;
  static struct option long_options[]= {  
    {"nbnds",    required_argument, NULL, 'n'},
    {"bnds",     required_argument, NULL, 'b'},
    {"nz",       required_argument, NULL, 'z'},
    {"dbnds",    required_argument, NULL, 'd'},
    {"stretch",    required_argument, NULL, 's'},
    {"grid_name",required_argument, NULL, 'o'},        
    {"center",   required_argument, NULL, 'c'},    
    {0, 0, 0, 0}
  };
  
    /*
   * process command line
   */
  errflg = argc <4;
  nbnds = 0;
  n1 = 0;
  n2 = 0;
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
    switch (c) {
    case 'n':
      nbnds = atoi(optarg);
      break;
    case 'b':
      strcpy(entry, optarg);
      n1 = get_double_entry(entry, bnds);
      break;
    case 'z':
      strcpy(entry, optarg);
      n2 = get_int_entry(entry, nz);
      break;
    case 'd':
      strcpy(entry, optarg);
      n3 = get_double_entry(entry, dbnds);
      break;
    case 's':
      stretch = atof(optarg);
      break;
    case 'o':
      strcpy(gridname, optarg);
      break;
    case 'c':
      strcpy(center, optarg);
      break;
    case '?':
      errflg++;      
    }

  if (errflg ) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    mpp_error("Wrong usage of this program, check arguments") ;
  }    

  /* check the command-line arguments to make sure the value are suitable */
  if( nbnds < 2 ) mpp_error("number of bounds specified through -nbnd should be an integer greater than 1");
  if( nbnds != n1 ) mpp_error("nbnds does not equal number entry specified through -bnd");
  if( nbnds-1 == n2 )
    use_legacy = 0;
  else if (nbnds == n3 )
    use_legacy = 1;
  else
    mpp_error("nbnds does not match number entry specified through --nz and --dbnds");

  /* generate grid */
  nk = 0;
  if( use_legacy )
    nk = MAX_GRID_LENGTH;
  else
    for(i=0; i<nbnds-1; i++) nk += nz[i];
  zeta = (double *)malloc((nk+1)*sizeof(double));
  {
    int npts;
    create_vgrid(nbnds, bnds, nz, dbnds, stretch, use_legacy, zeta, &npts, center);
    if(use_legacy) nk = npts;
  }
  
  /* define history to be the history in the grid file */
  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }

  sprintf(filename, "%s.nc", gridname);
  
  /* write out vertical grid into a netcdf file */
  {
    int fid, dim, varid;
    
    fid = mpp_open(filename, MPP_WRITE);
    dim  = mpp_def_dim(fid, "nzv", nk+1);
    varid = mpp_def_var(fid, "zeta", NC_DOUBLE, 1, &dim, 2, "standard_name", "vertical_grid_vertex",
			"units", "meters");
    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);    
    mpp_end_def(fid);
    mpp_put_var_value(fid, varid, zeta);
  
    mpp_close(fid);
  }

  /*  if(mpp_pe() == mpp_root_pe()) printf("Successfully generate vertical grid file %s\n", filename); */
  
  mpp_end();

  return 0;
}    
