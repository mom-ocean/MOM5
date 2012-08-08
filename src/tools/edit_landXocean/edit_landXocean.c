
/*******************************************************************************
  contact: Zhi Liang (Zhi.Liang)
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "constant.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"
#include "tool_util.h"

char *usage[] = {
  "",
  "  edit_landXocean --mosaic mosaic_grid                                                ",
  "                                                                                      ",
  "edit_landXocean will adjust the landXocean exchange grid cell area to put the flux    ",
  "onto ocean points that are adjacent to land points.                                   ",
  "                                                                                      ",
  "edit_landXocean takes the following flags:                                            ",
  "                                                                                      ",
  "REQUIRED:                                                                             ",
  "                                                                                      ",
  "--mosaic_grid mosaic_grid   specify the mosaic grid file to be edited. This mosaic    ",
  "                            file should be a coupler mosaic file, which contains link ",
  "                            to ocean topog file and the exchange grid file.           ",
  "                                                                                      ",
  NULL
};  

int main(int argc, char* argv[])
{
  unsigned int opcode = 0;
  int          option_index, c, n;
  char         *mosaic_file     = NULL;
  char         **lXo_file = NULL;
  char         **lXo_base_file = NULL;
  char         history[1024];
  char         mosaic_dir[STRING], ocn_mosaic[STRING], lnd_mosaic[STRING];
  int          **ocn_mask;
  int          *nxgrid;
  double       **area_lXo;
  int          **il_lXo, **jl_lXo, **io_lXo, **jo_lXo;
    
  
  int          *nx_ocn, *ny_ocn, *nx_lnd, *ny_lnd;
  int          nlXo, ntiles_ocn, ntiles_lnd;
  int          x_cyclic, y_cyclic, folded_north;
  
  int errflg = (argc == 1);
  
  static struct option long_options[] = {
    {"mosaic_grid",       required_argument, NULL, 'a'},
    {0, 0, 0, 0},
  };

  mpp_init(&argc, &argv);

  /* This tool will be run on one processor */
  if(mpp_npes() > 1) mpp_error("edit_landXocean: parallel is not supported yet, try running one single processor");
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'a':
      mosaic_file  = optarg;
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
  /* check the arguments */
  if( !mosaic_file    ) mpp_error("fregrid: mosaic_grid is not specified");

  strcpy(history,argv[0]);
  for(n=1;n<argc;n++)  {
    strcat(history, " ");
    strcat(history, argv[n]); 
  }

  /* read the mosaic file to get the landxocean exchange grid file name */
  {
    char tmp_file[STRING];
    char topog_file[STRING];
    size_t start[4], nread[4];
    int fid, vid;
    
    get_file_path(mosaic_file, mosaic_dir);
    fid = mpp_open(mosaic_file, MPP_READ);
    nlXo = mpp_get_dimlen( fid, "nfile_lXo");
    vid = mpp_get_varid(fid, "lXo_file");
    lXo_file = (char **)malloc(nlXo*sizeof(char *));
    lXo_base_file = (char **)malloc(nlXo*sizeof(char *));
    for(n=0; n<4; n++) {start[n]=0; nread[n]=1;}
    nread[1] = STRING;
    for(n=0; n<nlXo; n++) {
      start[0] = n;
      lXo_file[n] = (char *)malloc(STRING*sizeof(char));
      lXo_base_file[n] = (char *)malloc(STRING*sizeof(char));
      mpp_get_var_value_block(fid, vid, start, nread, lXo_base_file[n]);
      sprintf(lXo_file[n],  "%s/%s", mosaic_dir, lXo_base_file[n]);
    }
    
    /* read the ocean mosaic and land mosaic */
    vid = mpp_get_varid(fid, "ocn_mosaic_file");
    mpp_get_var_value(fid, vid, tmp_file);
    sprintf(ocn_mosaic, "%s/%s", mosaic_dir, tmp_file);
    vid = mpp_get_varid(fid, "lnd_mosaic_file");
    mpp_get_var_value(fid, vid, tmp_file);
    sprintf(lnd_mosaic, "%s/%s", mosaic_dir, tmp_file);    
    /* read the ocean topog file name */
    vid = mpp_get_varid(fid, "ocn_topog_file");
    mpp_get_var_value(fid, vid, tmp_file);
    sprintf(topog_file,  "%s/%s", mosaic_dir, tmp_file);
    mpp_close(fid);

  }
	 
  /*  get ntiles for ocean mosaic and land mosaic. currently we are assuming
      ntiles_ocean = 1 and ntiles_land = nlXo
      Also read the grid size for ocn_mosaic and land_mosaic.
  */
  ntiles_ocn = read_mosaic_ntiles(ocn_mosaic);
  ntiles_lnd = read_mosaic_ntiles(lnd_mosaic); 
  nx_ocn = (int *)malloc(ntiles_ocn*sizeof(int));
  ny_ocn = (int *)malloc(ntiles_ocn*sizeof(int));
  read_mosaic_grid_sizes(ocn_mosaic, nx_ocn, ny_ocn);
  nx_lnd = (int *)malloc(ntiles_lnd*sizeof(int));
  ny_lnd = (int *)malloc(ntiles_lnd*sizeof(int));
  read_mosaic_grid_sizes(lnd_mosaic, nx_lnd, ny_lnd);

  {
    char   **ocn_grid, **lnd_grid;
    char   tmp_file[STRING];
    size_t start[4], nread[4];
    int    fid, vid;
    int    ncontact;
    int   *tile1, *istart1, *iend1, *jstart1, *jend1;
    int   *tile2, *istart2, *iend2, *jstart2, *jend2;

    /* get the boundary condition for ocean mosaic. Since we are assuming there is only one
     tile in ocean mosaic, we could assume the boundary will be cyclic in x-direction,
     cyclic or folded-north for y-direction
    */
    x_cyclic     = 0;
    y_cyclic     = 0;
    folded_north = 0;    
    ncontact = read_mosaic_ncontacts(ocn_mosaic);
    tile1   = (int *)malloc(ncontact*sizeof(int));
    tile2   = (int *)malloc(ncontact*sizeof(int));
    istart1 = (int *)malloc(ncontact*sizeof(int));
    iend1   = (int *)malloc(ncontact*sizeof(int));
    jstart1 = (int *)malloc(ncontact*sizeof(int));
    jend1   = (int *)malloc(ncontact*sizeof(int));
    istart2 = (int *)malloc(ncontact*sizeof(int));
    iend2   = (int *)malloc(ncontact*sizeof(int));
    jstart2 = (int *)malloc(ncontact*sizeof(int));
    jend2   = (int *)malloc(ncontact*sizeof(int));
    if(ncontact >0) {
      read_mosaic_contact(ocn_mosaic, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2);
      if(ncontact <= 2) { /* x-cyclic of y-cyclic */
	for(n=0; n<ncontact; n++) {
	  if(istart1[n] == iend1[n] && istart2[n] == iend2[n] ) /* x_cyclic */
	    x_cyclic = 1;
	  else if(jstart1[n] == jend1[n] && jstart2[n] == jend2[n] ) {
	    if( jstart1[n] == jend1[n] )
	      folded_north = 1;
	    else
	      y_cyclic = 1;
	  }
	  else
	    mpp_error("edit_landXocean: for one-tile mosaic, the boundary condition should be either x-cyclic, y-cyclic or folded-north");
	}
      }
      else
	mpp_error("edit_landXocean: the number of contact should be either 0, 1 or 2");    
    }
      
    free(tile1);
    free(tile2);
    free(istart1);
    free(iend1);
    free(jstart1);
    free(jend1);
    free(istart2);
    free(iend2);
    free(jstart2);
    free(jend2);
  }
  
  /* read the landXocean exchange grid information */
  {
    area_lXo = (double **)malloc(nlXo*sizeof(double *));
    il_lXo   = (int    **)malloc(nlXo*sizeof(int    *));
    jl_lXo   = (int    **)malloc(nlXo*sizeof(int    *));
    io_lXo   = (int    **)malloc(nlXo*sizeof(int    *));
    jo_lXo   = (int    **)malloc(nlXo*sizeof(int    *));
    nxgrid   = (int     *)malloc(nlXo*sizeof(int     ));
    for(n=0; n<nlXo; n++) {
      nxgrid[n]   = read_mosaic_xgrid_size(lXo_file[n]);
      area_lXo[n] = (double *)malloc(nxgrid[n]*sizeof(double));
      il_lXo[n]   = (int    *)malloc(nxgrid[n]*sizeof(int   ));  
      jl_lXo[n]   = (int    *)malloc(nxgrid[n]*sizeof(int   ));
      io_lXo[n]   = (int    *)malloc(nxgrid[n]*sizeof(int   ));
      jo_lXo[n]   = (int    *)malloc(nxgrid[n]*sizeof(int   ));
      read_mosaic_xgrid_order1(lXo_file[n], il_lXo[n], jl_lXo[n], io_lXo[n], jo_lXo[n], area_lXo[n]);    
    }
  }

  /*calculate ocean grid mask based on landXocean exchange grid information
    We also assume ocean cell is always full-cell. so the mask value will be
    either 0 (land) or 1 (ocean) */
  {
    int nx, ny, nxp1, nyp1, nxp2, nyp2;
    int n, m, i, j;
    
    ocn_mask = (int **)malloc(ntiles_ocn*sizeof(int *));
    nx   = nx_ocn[0];
    ny   = ny_ocn[0];
    nxp1 = nx + 1;
    nyp1 = ny + 1;
    nxp2 = nx + 2;
    nyp2 = ny + 2;
    for(n=0; n<ntiles_ocn; n++) {
      ocn_mask[n] = (int *)malloc(nxp2*nyp2*sizeof(int) );
      for(m=0; m<nxp2*nyp2; m++) ocn_mask[n][m] = 0;
    }

    for(n=0; n<nlXo; n++) {
      for(m=0; m<nxgrid[n]; m++) {
	i = io_lXo[n][m];   
	j = jo_lXo[n][m];
	ocn_mask[n][j*nxp2+i] = 1;
      }
    }
    /* fill the halo */
    if(x_cyclic) {
      for(j=1; j<nyp1; j++) {
	ocn_mask[0][j*nxp2     ] = ocn_mask[0][j*nxp2+nx];      /* West */
	ocn_mask[0][j*nxp2+nxp1] = ocn_mask[0][j*nxp2+ 1];      /* East */
      }
    }
    if(y_cyclic) {
      for(i=1; i<nxp1; i++) {
	ocn_mask[0][          i] = ocn_mask[0][ny*nxp2+i];      /* south */
	ocn_mask[0][nyp1*nxp2+i] = ocn_mask[0][   nxp2+i];      /* north */
      }
      if(x_cyclic) {
	ocn_mask[0][0]              = ocn_mask[0][ny*nxp2+nx];   /* southwest */
	ocn_mask[0][nxp1]           = ocn_mask[0][ny*nxp2+1];    /* southeast */    
	ocn_mask[0][nyp1*nxp2+nxp1] = ocn_mask[0][nxp2+1];       /* northeast */
	ocn_mask[0][nyp1*nxp2]      = ocn_mask[0][nxp2+nx];      /* northwest */
      }
    }
    else if(folded_north) {
      if( !x_cyclic) mpp_error("Error from program edit_landXocean: when the y-boundary condition "
			       "is folded-north, the x-boundary condition must be cyclic");
      for(i=1; i<nxp1; i++) {
	ocn_mask[0][nyp1*nxp2+i] = ocn_mask[0][ny*nxp2+nxp1-i];
      }
      ocn_mask[0][nyp1*nxp2]   = ocn_mask[0][ny*nxp2+1];
      ocn_mask[0][nyp1*nxp2+nxp1] = ocn_mask[0][ny*nxp2+nx];
    } 
  }
  /* figure out the exchange grid cells that the runoff will put to */
  /* set such point with value 1 and other points with value 0 */
  {
    double *scale, *tot_area, *set_area;
    int nxl, nyl, n, m, io, jo, il, jl, ioff, joff, ii, jj;
    int is_selected, nxp2;
    char orig_file[STRING], cmd[512];
    int fid, id_scale, id_ncell;
    
    nxp2 = nx_ocn[0] + 2;
    for(n=0; n<nlXo; n++) {
      nxl = nx_lnd[n];
      nyl = ny_lnd[n];
      tot_area = (double *)malloc(nxl*nyl*sizeof(double));
      set_area = (double *)malloc(nxl*nyl*sizeof(double));
      scale    = (double *)malloc(nxgrid[n]*sizeof(double));
      for(m=0; m<nxl*nyl; m++) {
	tot_area[m] = 0;
	set_area[m] = 0;
      }
      for(m=0; m<nxgrid[n]; m++) {
	il = il_lXo[n][m] - 1;
	jl = jl_lXo[n][m] - 1;  
	io = io_lXo[n][m];
	jo = jo_lXo[n][m];
	is_selected = 0;
	for(joff=-1; joff<=1; joff++) for(ioff=-1; ioff<=1; ioff++) {
	  if(ioff == 0 && joff == 0) continue;
	  ii = io+ioff;
	  jj = jo+joff;
	  if( ocn_mask[0][jj*nxp2+ii] == 0) {
	    is_selected = 1;
	  }
	}
	tot_area[jl*nxl+il] += area_lXo[n][m];
	if(is_selected) {
	  scale[m] = 1;
	  set_area[jl*nxl+il] += area_lXo[n][m];
	}
	else {
	  scale[m] = 0;
	}
      }

      for(m=0; m<nxgrid[n]; m++) {
	il = il_lXo[n][m] - 1;
	jl = jl_lXo[n][m] - 1;
	if(set_area[jl*nxl+il] == 0)
	  scale[m] = 1;
	else if(set_area[jl*nxl+il] > 0) 
	  scale[m] *= (tot_area[jl*nxl+il]/set_area[jl*nxl+il]);
	else {
	  printf("at point (%d, %d), set_area = %15.12f is negative\n", il, jl, set_area[jl*nxl+il]);
	  mpp_error("program edit_landXocean: set_area is negative for some point");
	}
      }
      free(tot_area);
      free(set_area);
      /* write the scale information into the xgrid file */
      if(!strcmp(lXo_file[n], lXo_base_file[n])) {
	sprintf(orig_file, "%s.orig", lXo_file[n]);
	sprintf(cmd, "cp %s %s", lXo_file[n], orig_file);
	system(cmd);
      }
      else {
	sprintf(cmd, "cp %s %s", lXo_file[n], lXo_base_file[n]);
        system(cmd);
      }
 
      fid = mpp_open(lXo_base_file[n], MPP_APPEND);
      mpp_redef(fid);
      id_ncell = mpp_get_dimid(fid, "ncells");
      id_scale = mpp_def_var(fid, "scale", MPP_DOUBLE, 1, &id_ncell, 0);
      mpp_end_def(fid);
      mpp_put_var_value(fid, id_scale, scale);
      mpp_close(fid);
    }
  }
      
  /* free the memory */
  for(n=0; n<nlXo; n++) {
    free(area_lXo[n]);
    free(il_lXo[n]);
    free(jl_lXo[n]);
    free(io_lXo[n]);
    free(jo_lXo[n]);
  }
  free(area_lXo);
  free(il_lXo);
  free(jl_lXo);
  free(io_lXo);
  free(jo_lXo);

  printf("Successfully running river_regrid and the following output file are generated.\n");

  mpp_end();

  return 0;  
} /* end of main */

			       
    
