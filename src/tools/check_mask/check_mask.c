#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "read_mosaic.h"


char *usage[] = {
  "",
  " check_mask --grid_file grid_file [--ocean_topog topog_file] [--min_pe #] [--max_pe #]  ",
  "            [--halo #] [--sea_level #] [--show_valid_only] [--have_obc]                 ",
  "            [--direction d(1)..,d(nobc)] [--is is(1)...,is(nobc)]                       ",
  "            [--ie ie(1)...,ie(nobc)], [--js js(1)...,js(nobc)], [--je je(1)...,je(nobc)]",
  "            [--layout layout(1),layout(2)]                                              ",
  "                                                                                        ",
  " check_mask is used to configure the processors which contains all land points to be    ",
  " masked out. This program is supposed to run on single processor. This tool will print  ",
  " out number of processors to be masked, layout of the domain, masked out processors     ",
  " list. These information will be written into files, with file name is                  ",
  " mask_table.n_mask.layout(1)Xlayout(2). The model is expected to read file mask_table   ",
  " if you want to mask out some processors.                                               ",
  "                                                                                        ",
  " The purpose of using mask domain is to decrease the processor usage without            ",
  " decreasing the performance, without changing bits relative to the case run without the ",
  " mask domain option. The mask domain option acts to mask out regions that have zero     ",
  " ocean points (i.e., points that are all land) in both the ocean and sea ice models.    ", 
  "                                                                                        ",
  " check_mask takes the following flags:                                                  ",
  "                                                                                        ",
  "REQUIRED:                                                                               ",
  "                                                                                        ",
  " --grid_file grid_file    Specify the grid file to be used in the model run. The        ",
  "                          grid_file could be old version of grid (contains field wet)   ",
  "                          or the new version mosaic grid, which contains field gridfiles",
  "                          When grid_file is mosaic grid, it should be the ocean solo    ",
  "                          mosaic and ocean_topog file needs to be specified.            ",
  "                                                                                        ",
  " OPTIONAL FLAGS                                                                         ",
  "                                                                                        ",
  " --ocean_topog topog_file Specify the topography information for ocean mosaic.          ",
  "                          The field name of the topography is depth_tile# or depth      ",
  "                          when ntiles = 1, The topography data is positive down.        ",
  "                                                                                        ",
  " --min_pe #               Specify the smallest processor count to be checked            ",
  "                                                                                        ",
  " --max_pe #               Specify the largest processor count to be checked             ",
  "                                                                                        ",
  " --layout #,#             specify the layout to be checked. When layout is specified,   ",
  "                          --min_pe and --max_pe will be ignored.                        ",
  "                                                                                        ",
  "                                                                                        ",
  " --halo #                 Specify the halo size in the ocean model. When there is no    ",
  "                          ocean points on a processor (including halo data), the        ",
  "                          processor will be masked out. Default value is 1.             ",
  "                                                                                        ",
  " --sea_level #            Specify the sea level ( in meters ) and its value will be used",
  "                          to determine land/sea mask. When topography of a grid cell    ",
  "                          is less than sea level, this grid cell will be land,          ",
  "                          otherwise it will be ocean. Default value is 0                ",
  "                                                                                        ",
  " --show_valid_only        When set, only layouts valid with OBC are shown.              ",
  "                                                                                        ",
  " --nobc #                 number of open boundary condition. Its value should be less   ",
  "                          than MAX_OBC (current is 4 ). default vaule is 0 for nobc.    ",
  "                                                                                        ",
  " --direction d(1)..,d(nobc) open boundary direction. Each element value should be west, ",
  "                            east, south or north.                                       ",
  "                                                                                        ",
  " --is is(1)...,is(nobc)   starting i-index for open boundary condition. The index is    ",
  "                          Fortran index.                                                ",
  "                                                                                        ",
  " --ie ie(1)...,ie(nobc)   ending i-index for open boundary condition. The index is      ",
  "                          Fortran index.                                                ",
  "                                                                                        ",
  " --js js(1)...,js(nobc)   starting j-index for open boundary condition. The index is    ",
  "                          Fortran index.                                                ",
  "                                                                                        ",
  " --je je(1)...,je(nobc)   ending i-index for open boundary condition. The index is      ",
  "                          Fortran index.                                                ",
  "                                                                                        ",
  " Example 1: use mask domain in fully coupled model                                      ",
  "                                                                                        ",
  "  1. check_mask --grid_file /archive/z1l/tools/input/cm2m_grid/ocean_mosaic.nc          ",
  "                --ocean_topog /archive/z1l/tools/input/cm2m_grid/topog.nc               ",
  "                --min_pe 100 --max_pe 200                                               ",
  "     This will create a list of mask_table with filename                                ",
  "     mask_table.n_mask.layout(1)Xlayout(2). For example mask_table.9.18x10 has          ",
  "     layout=(18,10) and 9 domain processor is masked out. For the step 2, 3 and 4 we    ",
  "     use mask_table.9.18x10 as the example                                              ",
  "  2. copy mask_table.9.18x10 into 'INPUT' under working directory.                      ",
  "  3. set ocean_model_nml layout = 18,10                                                 ",
  "         ocean_model_nml mask_table = 'INPUT/mask_table.9.18x10'                        ",
  "  4. For concurrent run, set coupler_nml ocean_npes = 171 (18x10-9)                     ",
  "                                                                                        ",
  "  NOTE: For coupled models using concurrent coupling, you may only wish to set          ",
  "        ocean_model_nml variable mask_table, and not set ice_model_nml mask_table.      ",
  "        The reason for avoiding the ice_model use of mask_table with fully coupled      ",
  "        models is that,                                                                 ", 
  "        1. the mask_table for the sea ice may not be necessary for performance          ",
  "           enhancement in fully coupled models                                          ",
  "        2. most importantly, setting mask_table for sea ice complicates the layout      ",
  "           for the ice model in fully coupled models.                                   ",
  "                                                                                        ",
  " Example 2: use mask domain in sea-ice model ( baltic1 experiment, ocean, ice,          ",
  "            atmosphere and land are all on the same grid )                              ",
  "  1. check_mask --grid_file /archive/z1l/tools/input/baltic1_grid_spec.nc               ",
  "                 --min_pe 60 --max_pe 80                                                ",
  "     This will create a list of mask_table with filename                                ",
  "     mask_table.n_mask.layout(1)Xlayout(2). For example mask_table.26.6x11 has          ",
  "     layout=(6,11) and 26 domain processor is masked out. For the step 2 and 3 we       ",
  "     use mask_table.26.6x11 as the example                                              ",
  "                                                                                        ",
  "  2.  copy mask_table.26.6x11    into 'INPUT' under working directory.                  ",
  "  3.  set ocean_model_nml layout = 6,11                                                 ",
  "          ocean_model_nml mask_table = 'INPUT/mask_table.26.6x11                        ",
  "          ice_model_nml layout = 6,11                                                   ",
  "          ice_model_nml mask_table = 'INPUT/mask_table.26.6x11                          ",
  "          atmos_model_nml layout = 6,11                                                 ",
  "          atmos_model_nml mask_table = 'INPUT/mask_table.26.6x11                        ",
  "          land_model_nml layout = 6,11                                                  ",
  "          land_model_nml mask_table = 'INPUT/mask_table.26.6x11                         ",
  "   This experiment will run on 40 (6x11-26) processors.                                 ", 
  "",
  NULL };

int get_text_entry(char *line, char *value[]);
void get_grid_size( const char *grid_file, int grid_version, int *nx, int *ny );
void get_ocean_mask(const char *grid_file, int grid_version, double *mask, double sea_level, int nx, int ny );
void check_mask(int nx, int ny, const double *wet_in, int cyclic_x, int cyclic_y,
		int is_tripolar, int halo, int min_pe, int max_pe, int layout[], int show_valid_only, int nobc,
		char *direction[], int *is, int *ie, int *js, int *je  );

#define MAX_OBC 4

int main (int argc, char *argv[])
{
  char *grid_file=NULL;
  char *topog_file=NULL;
  int min_pe = 4;
  int max_pe = 128;
  int halo   = 1;
  int show_valid_only = 0;
  double sea_level = 0.0;
  int have_obc = 0;
  int nobc = 0;
  char **direction = NULL;
  int is[] = {-999, -999, -999, -999};
  int ie[] = {-999, -999, -999, -999};
  int js[] = {-999, -999, -999, -999};
  int je[] = {-999, -999, -999, -999};
  int layout[] = {0,0};

  int num_layout_entry;
  int grid_version = 0;
  int cyclic_x=0;
  int cyclic_y=0;
  int is_tripolar=0;
  int nx=0, ny=0;
  double *mask=NULL;
  
  int nobc1=0, nobc2=0, nobc3=0, nobc4=0, nobc5=0;
  char entry[512];
  int errflg=0, c, option_index = 0;

  static struct option long_options[] = {
    {"grid_file",       required_argument, NULL, 'a'},
    {"ocean_topog",     required_argument, NULL, 'b'},
    {"min_pe",          required_argument, NULL, 'c'},
    {"max_pe",          required_argument, NULL, 'd'},
    {"halo",            required_argument, NULL, 'e'},
    {"sea_level",       required_argument, NULL, 'f'},
    {"show_valid_only", no_argument,       NULL, 'g'},
    {"nobc",            required_argument, NULL, 'i'},
    {"direction",       required_argument, NULL, 'j'},
    {"is",              required_argument, NULL, 'k'},
    {"ie",              required_argument, NULL, 'l'},
    {"js",              required_argument, NULL, 'm'},
    {"je",              required_argument, NULL, 'n'},
    {"layout",          required_argument, NULL, 'o'},
    {NULL, 0, NULL, 0}
  };

  mpp_init(&argc, &argv);

  if( mpp_npes() > 1 ) mpp_error("check_mask: This program could only be run on single processor");

  /*
   * process command line
   */
  {
    int n;
    direction = (char **)malloc(MAX_OBC*sizeof(char *));
    for(n=0; n<MAX_OBC; n++) {
      direction[n] = (char *)malloc(10*sizeof(char));
    }
  }
  while ((c = getopt_long(argc, argv, "h:", long_options, &option_index) ) != -1) {
    switch (c) {
    case 'a': 
      grid_file = optarg;
      break;
    case 'b': 
      topog_file = optarg;
      break;
    case 'c':
      min_pe = atoi(optarg);
      break;
    case 'd':
      max_pe = atoi(optarg);
      break;
    case 'e':
      halo = atoi(optarg);
      break;
    case 'f':
      sea_level = atoi(optarg);
      break;
    case 'g':
      show_valid_only = 1;
      break;      
    case 'i':
      nobc = atoi(optarg);
      break;
    case 'j':
      strcpy(entry, optarg);
      nobc1 = get_text_entry(entry, direction);
      break;
    case 'k':
      strcpy(entry, optarg);
      nobc2 = get_int_entry(entry, is);
      break;
    case 'l':
      strcpy(entry, optarg);
      nobc3 = get_int_entry(entry, ie);
      break;
    case 'm':
      strcpy(entry, optarg);
      nobc4 = get_int_entry(entry, js);
      break;
    case 'n':
      strcpy(entry, optarg);
      nobc5 = get_int_entry(entry, je);
      break;
    case 'o':
      strcpy(entry, optarg);
      num_layout_entry = get_int_entry(entry, layout);
      if(num_layout_entry != 2) mpp_error("check_mask: layout should be specified by --layout #,#");
      break;      
    case '?':
      errflg++;
    }
  }
  
  if (errflg || !grid_file) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }

  if( nobc > MAX_OBC) mpp_error("Error from check_mask: nboc > MAX_OBC");
  if( nobc1 != nobc ) mpp_error("Error from check_mask: number of entry specified in --direction is not equal to nobc");
  if( nobc2 != nobc ) mpp_error("Error from check_mask: number of entry specified in --is is not equal to nobc");
  if( nobc3 != nobc ) mpp_error("Error from check_mask: number of entry specified in --ie is not equal to nobc");
  if( nobc4 != nobc ) mpp_error("Error from check_mask: number of entry specified in --js is not equal to nobc");
  if( nobc5 != nobc ) mpp_error("Error from check_mask: number of entry specified in --je is not equal to nobc");

  /* Change the Fortran index to C index */
  if(nobc > 0) {
    int n;
    for(n=0; n<nobc; n++){
      is[n]--;
      ie[n]--;
      js[n]--;
      je[n]--;
    }
  }
      
  /* print out the input arguments */
  {
    int n;
    
    if(layout[0]*layout[1] > 0) {
      min_pe = layout[0]*layout[1];
      max_pe = min_pe;
      printf("\n ===>NOTE from check_mask: when layout is specified, min_pe and max_pe is set to layout(1)*layout(2)=%d\n",
	     layout[0]*layout[1]);
    }
    
    printf("\n ===>NOTE from check_mask: Below is the list of command line arguments.\n\n");
    printf("grid_file = %s\n", grid_file);
    if( topog_file )
      printf("topog_file = %s\n", topog_file);
    else
      printf("topog_file is not specified");
    printf("min_pe = %d\n", min_pe);
    printf("max_pe = %d\n", max_pe);
    printf("layout = %d, %d\n", layout[0], layout[1]);
    printf("halo = %d\n", halo);
    printf("sea_level = %g\n", sea_level);
    if( show_valid_only )
      printf("show_valid_only is set\n");
    else
      printf("show_valid_only is not set\n");

    printf("nobc = %d\n", nobc);
    for(n=0; n<nobc; n++) {
      printf("obc #%d, direction=%s, is=%d, ie=%d, js=%d, je=%d\n", n+1, direction[n], is[n], ie[n], js[n], je[n]);
    }  
    
    printf("\n ===>NOTE from check_mask: End of command line arguments.\n");
  }

  /* find out grid version */
  if(mpp_field_exist(grid_file, "wet") ) {
    grid_version = VERSION_1;
    printf("\n ===>NOTE from check_mask: the grid file is version 1 grid which contains field wet\n");
  }
  else if(mpp_field_exist(grid_file, "gridfiles") ) {
    grid_version = VERSION_2;
    printf("\n ===>NOTE from check_mask: the grid file is version 2 (mosaic grid) grid which contains field gridfiles\n");
  }
  else
    mpp_error("Error from check_mask: both wet and gridfiles do not exist in grid file");

  /* grid_version is VERSION_2, ocean_topog must be specified */
  if( grid_version == VERSION_2 && !topog_file )
    mpp_error("Error from check_mask: when grid is VERSION_2, ocean_topog must be specified");
      
  /* get the boundary condition */
  get_boundary_type( grid_file, grid_version, &cyclic_x, &cyclic_y, &is_tripolar );

  /* get the grid resolution */
  get_grid_size( grid_file, grid_version, &nx, &ny );  

  mask = (double *)malloc(nx*ny*sizeof(double));
  
  /* get the land sea mask */
  if( grid_version == VERSION_1)
    get_ocean_mask( grid_file, grid_version, mask, sea_level, nx, ny );
  else
    get_ocean_mask( topog_file, grid_version, mask, sea_level, nx, ny );

  /* check mask */
  check_mask(nx, ny, mask, cyclic_x, cyclic_y, is_tripolar, halo, min_pe, max_pe, layout, show_valid_only, nobc, direction, is, ie, js, je );
  free(mask);
  
  printf("\n***** Congratulation! You have successfully run check_mask\n");
  
  return 0;
  
}  

int get_text_entry(char *line, char *value[])
{
  char* pch;
  int num;
  
  pch = strtok(line, ", ");
  num = 0;
  while( pch != NULL) {
    strcpy(value[num], pch);
    num++;    
    pch = strtok(NULL, ", ");
  }
  return num;
};

  
void check_mask(int nx, int ny, const double *wet_in, int cyclic_x, int cyclic_y,
		int is_tripolar, int halo, int min_pe, int max_pe, int layout_in[], int show_valid_only, int nobc,
		char *direction[], int *is, int *ie, int *js, int *je  )
{
  int nxd, nyd;
  int isc, iec, jsc, jec;
  int isd, ied, jsd, jed;
  int isg, ieg, jsg, jeg;
  int i, j, ii, jj, n1, n2, np, no, n;
  int on_domain=0;
  int nmask=0, nerror=0;
  int *mask_list=NULL;
  double *wet=NULL;
  int *mask=NULL;
  int *ibegin=NULL, *iend=NULL, *jbegin=NULL, *jend=NULL;
  int *obc_error=NULL;
  
  int layout[2]={0,0};

  nxd = nx+2*halo;
  nyd = ny+2*halo;
  wet = (double *)malloc(nxd*nyd*sizeof(double));
  for(i = 0; i<nxd*nyd ; i++) wet[i] = 0.0;

  for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
    n1 = j*nx+i;
    n2 = (j+halo)*nxd + i+halo;
    wet[n2] = wet_in[n1];
  }
  
  /* fill the global halo */
  if( cyclic_x ) {
    for(j=0; j<ny; j++) {
      for(i=0; i<halo; i++) {  /* west halo */
	n1 = (j+halo)*nxd + nx+halo-1-i;
	n2 = (j+halo)*nxd+i;
	wet[n2] = wet[n1];
      }
      for(i=0; i<halo; i++) { /* east halo */
	n1 = (j+halo)*nxd + halo+i;
	n2 = (j+halo)*nxd + nx+halo+i;
	wet[n2] = wet[n1];
      }
    }
  }

  if( cyclic_y ) {
    for(i=0; i<nx; i++) {
      for(j=0; j<halo; j++) {  /* south halo */
	n1 = (ny+halo-1-j)*nxd + i+halo;
	n2 = j*nxd + i+halo;
	wet[n2] = wet[n1];
      }
      for(j=0; j<halo; j++) { /* north halo */
	n1 = (halo+j)*nxd + i+halo;
	n2 = (ny+halo+j)*nxd + i+halo;
	wet[n2] = wet[n1];
      }
    }
  }

  if(is_tripolar) {
    for(i=0; i<nxd; i++) for(j=0; j<halo; j++) {
	n1 = (ny+halo-1-j)*nxd + nxd-1-i;
	n2 = (ny+halo+j)*nxd + i;
	wet[n2] = wet[n1];
    }
  }

  printf("\n==>NOTE from check_mask: Checking for possible masking:\n");
  printf("==>NOTE from check_mask: Assume %d halo rows\n", halo);
  printf("==>NOTE from check_mask: Total domain size is %d, %d\n", nx, ny);
  if( show_valid_only && nobc>0 )
    printf("==>NOTE from check_mask: Only layouts valid with OBC are shown.");
  else
    if (nobc>0)printf("==>NOTE from check_mask: Also layouts invalid with OBC are shown.");

  isg = 0;
  ieg = nx-1; 
  jsg = 0; 
  jeg = ny-1; 
  
  mask   = (int *)malloc(max_pe*sizeof(int));
  ibegin = (int *)malloc(max_pe*sizeof(int));
  iend   = (int *)malloc(max_pe*sizeof(int));
  jbegin = (int *)malloc(max_pe*sizeof(int));
  jend   = (int *)malloc(max_pe*sizeof(int));
  if(nobc>0) obc_error=(int *)malloc(max_pe*sizeof(int));

  for(np=min_pe; np<=max_pe; np++) {
    if( layout_in[0]*layout_in[1] == np) {
      layout[0] = layout_in[0];
      layout[1] = layout_in[1];
    }
    else {
      mpp_define_layout(nx, ny, np, layout);
    }
    if( layout[0] > nx || layout[1] > ny ) continue;
    mpp_compute_extent(nx,layout[0],ibegin,iend);
    mpp_compute_extent(ny,layout[1],jbegin,jend);
    /* add halo to the compute domain */
    for(i=0; i<layout[0]; i++){
      ibegin[i] += halo;
      iend[i]   += halo;
    }
    for(i=0; i<layout[1]; i++){
      jbegin[i] += halo;
      jend[i]   += halo;
    }    
    for(i=0; i<np; i++) mask[i] = 0;
    if (nobc>0) {
      for(i=0; i<np; i++)obc_error[i] = 0;
    }
    for(j=0; j<layout[1]; j++) {
      jsc = jbegin[j]; 
      jec = jend[j]; 
      jsd = jbegin[j] - halo;
      jed = jend[j] + halo;
      for(i=0; i<layout[0]; i++) {
	isc = ibegin[i];
	iec = iend[i]; 
	isd = ibegin[i] - halo;
	ied = iend[i] + halo;
	for(jj=jsd; jj<=jed; jj++) for(ii=isd; ii<=ied; ii++) {
	  if(wet[jj*nxd+ii] > 0.5) {
	    mask[j*layout[0]+i] = 1;
	    goto found ;
	  }
	}
      found: /* do nothing*/
	if(nobc>0) {
	  for(no=0; no<nobc; no++) {
	    /* check, if the open boundary is at a domain    */
	    on_domain=0;
	    for(jj=jsd; jj<=jed; jj++) for(ii=isd; ii<=ied; ii++) {
	      if(jj <= js[no] && jj <= je[no] && ii >= is[no] && ii <= ie[no]) on_domain = 1;
	    }
	    if( !on_domain ) continue;

	    /* check, if the open boundary is at a domain halo	   */
	    if(!strcmp(direction[no], "west") ) {
	      if (isg != isc) {
		if ((is[no]-isc) < 0)    obc_error[j*layout[0]+i]=no;
	      }
	      if ((ied-is[no]) <= halo)  obc_error[j*layout[0]+i]=no;
	    }
	    else if(!strcmp(direction[no], "east") ) {
	      if (ieg != iec) {
		if ((iec-is[no]) < 0)    obc_error[j*layout[0]+i]=no;
	      }
	      if ((is[no]-isd) <= halo)   obc_error[j*layout[0]+i]=no;
	    }
	    else if(!strcmp(direction[no], "south") ) {
	      if (jsg != jsc) {
		if ((js[no]-jsc) < 0)    obc_error[j*layout[0]+i]=no;
	      }
	      if ((jec-js[no]) <= halo)  obc_error[j*layout[0]+i]=no;
	    }
	    else if(!strcmp(direction[no], "north") ) {
	      if (jeg != jec) {
		if ((jec-js[no]) < 0)    obc_error[j*layout[0]+i]=no;
	      }
	      if ((js[no]-jsc) <= halo)  obc_error[j*layout[0]+i]=no;
	    }
	  }
	}
      }
    }

    nmask=0;
    for(i=0; i<layout[0]*layout[1]; i++) {
      if(mask[i] ==0) nmask++;
    }

    if(nmask == np) printf("(WARNING from program check_mask: all points are land");
    if(nmask > 0) {
      mask_list = (int *)malloc(2*nmask*sizeof(int));

      n = 0;
      for(j=0; j<layout[1]; j++) for(i=0; i<layout[0]; i++) {
	if(mask[j*layout[0]+i] == 0) {
	  mask_list[2*n] = i+1;
	  mask_list[2*n+1] = j+1;
	  n++;
	}
      }

      if(show_valid_only && nobc>0) {
	nerror = 0;
	for(i=0; i<layout[0]*layout[1]; i++) {
	  if(obc_error[i] ==0) nerror++;
	}
	if(nerror>0) continue;
      }
      printf("\n_______________________________________________________________________\n");
      printf("\nNOTE from check_mask: The following is for using model source code with version older than siena_201207,\n");
      printf("Possible setting to mask out all-land points region, for use in coupler_nml");
      printf("Total number of domains = %d\n", np);
      printf("Number of tasks (excluded all-land region) to be used is %d\n", np - nmask);
      printf("Number of regions to be masked out = %d\n", nmask);
      printf("The layout is %d, %d\n", layout[0], layout[1]);
      printf("Masked and used tasks, 1: used, 0: masked\n");
      for(j=layout[1]-1; j>=0; j--) {
	for(i=0; i<layout[0]; i++) printf("%d", mask[j*layout[0]+i]);
	printf("\n");
      }
      if( nobc > 0 && nerror > 0) {
	printf("OBC at halos is not allowed yet. Do not use this configuration.\n");
	for(j=0; j<layout[1]; j++) {
	  for(i=0; i<layout[0]; i++) printf("%d ", obc_error[j*layout[0]+i]);
	  printf("\n");
	}
      }
      
      printf(" domain decomposition\n");
      for(i=0; i<layout[0]; i++) printf("%4d", iend[i]-ibegin[i]+1);
      printf("\n");
      for(j=0; j<layout[1]; j++) printf("%4d", jend[j]-jbegin[j]+1);
      printf("\n");
      printf(" used=%d, masked=%d, layout=%d,%d\n ", np - nmask, nmask, layout[0], layout[1]);
      
      printf("To chose this mask layout please put the following lines in ocean_model_nml and/or ice_model_nml\n ");
      printf("nmask = %d\n", nmask);
      printf("layout = %d, %d\n", layout[0], layout[1]);
      printf("mask_list = ");
      for(n=0; n<nmask; n++) {
	printf("%d,%d", mask_list[2*n], mask_list[2*n+1]);
	if( n < nmask-1 )printf(",");
      }
      printf("\n\n");

      /* Finally write the n_mask, layout_mask and mask_list information to a table file */
      {
	char file[128];
	FILE *fp;

	sprintf(file, "mask_table.%d.%dx%d", nmask, layout[0], layout[1]);
	fp=fopen(file, "w");
	fprintf(fp, "%d\n", nmask);
	fprintf(fp, "%d, %d\n", layout[0], layout[1]);
	for(n=0; n<nmask; n++) {
	  fprintf(fp, "%d,%d\n", mask_list[2*n], mask_list[2*n+1]);
	}
        printf("\n_______________________________________________________________________\n");
        printf("\nNOTE from check_mask: The following is for using model source code with version siena_201207 or newer,\n");
	printf("                      specify ocean_model_nml/ice_model_nml/atmos_model_nml/land_model/nml \n");
	printf("                      variable mask_table with the mask_table created here.\n");
	printf("                      Also specify the layout variable in each namelist using corresponding layout\n"); 
	fclose(fp);
      }
    }
  }
  

}/* check_mask */
  
 
void get_ocean_mask(const char *grid_file, int grid_version, double *mask, double sea_level, int nx, int ny )
{
  int fid, vid;
  double *depth=NULL;
  int i;
  
  if( grid_version == VERSION_1) {
    fid = mpp_open( grid_file, MPP_READ );
    vid = mpp_get_varid(fid, "wet");
    mpp_get_var_value(fid, vid, mask);
    mpp_close(fid);
  }
  else if( grid_version == VERSION_2 ) {
    fid = mpp_open( grid_file, MPP_READ );
    vid = mpp_get_varid(fid, "depth");
    depth = (double *)malloc(nx*ny*sizeof(double));
    mpp_get_var_value(fid, vid, depth);
    mpp_close(fid);
    
    for(i=0; i<nx*ny;i++) {
      if(depth[i] > sea_level)
	mask[i] = 1.0;
      else
	mask[i] = 0.0;
    }
    
    free(depth);
  }
  else {
    mpp_error("==>Error from program check_mask: grid_version should be VERSION_1 or VERSION_2");
  }
    
}


void get_grid_size( const char *grid_file, int grid_version, int *nx, int *ny )
{
  int fid, vid;
  char xaxis_name[32];
  char yaxis_name[32];
  
  if( grid_version == VERSION_1 ) {
    fid = mpp_open(grid_file, MPP_READ);
    vid = mpp_get_varid(fid, "wet");
    mpp_get_var_dimname(fid, vid, 1, xaxis_name);
    mpp_get_var_dimname(fid, vid, 0, yaxis_name);
    *nx = mpp_get_dimlen(fid, xaxis_name );
    *ny = mpp_get_dimlen(fid, yaxis_name );
    mpp_close(fid);
  }
  else if(grid_version == VERSION_2 ) { /* mosaic grid */
    read_mosaic_grid_sizes( grid_file, nx, ny );
  }
  else {
    mpp_error("==>Error from program check_mask: grid_version should be VERSION_1 or VERSION_2");
  }
}


		    

