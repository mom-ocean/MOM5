/*
   Copyright (C) 2007,2009-2010,2012,2013 Remik Ziemlinski

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

	 20070821 rsz Created.
	 20070824 rsz Works with file test3.nc (0,1,2,3,4)D variables.
	 20121130 rsz Fixes attributed start tile index to 1-based to be compatible with mppnccombine (thanks to Zhi Lang).
   20130303 rsz Uneven partition should use symmetric layout (instead of remainder). Added --io-layout-x/y -j/-i -n -p -w options.
*/
#include "mppncscatter.h"

/*-------------------------------------------------------------------*/
void printsizetarray(size_t *a, int n) {
	int i=0;
	for(;i < n; ++i)
		fprintf(stdout, "%d ", (int)a[i]);

	fprintf(stdout, "\n");
}
/*-------------------------------------------------------------------*/
/* Return number of final x, y divisions depending if io_layout. */
void get_num_divs(mnsopts* opts, int* nx, int* ny) {
  if (opts == NULL) {
    *nx = 0;
    *ny = 0;
  } else if (opts->nxio && opts->nyio) {
    *nx = opts->nxio;
    *ny = opts->nyio;
  } else {
    *nx = opts->nx;
    *ny = opts->ny;  
  }
}
/*-------------------------------------------------------------------*/
int get_num_files(mnsopts* opts) {
  int nx, ny;
  if (opts == NULL) return 0;
  
  get_num_divs(opts, &nx, &ny);
  return nx*ny;
}
/*-------------------------------------------------------------------*/
/* Memory copies subdomain to preallocated output pointer (either t,s,i,f,d datatype pointers).
	 Record variables should not pass in record info in the start,count,ndim,
	 so be sure to pass the pointers offset by 1 entry and ndim-1.
 */
void hyperslabcopy(nc_type type, size_t *dimlen, int *dimids, size_t *start, size_t *count, int ndim, char *ti, short *si, int *ii, float *fi, double *di, char *t, short *s, int *i, float *f, double *d) 
{
	size_t k, j;
	size_t offset=0;
	size_t i0, stridek, stridej;

	switch(ndim) {
	case 0:
		switch(type) {
		case NC_BYTE: case NC_CHAR:
			*t = *ti;
			break;
		case NC_SHORT: 
			*s = *si;
			break;
		case NC_INT:
			*i = *ii;
			break;
		case NC_FLOAT:
			*f = *fi;
			break;
		case NC_DOUBLE:
			*d = *di;
			break;
		}
		break;
	case 1:
		switch(type) {
		case NC_BYTE: case NC_CHAR:
			memcpy(t, ti+start[0], count[0]*sizeof(char));
			break;
		case NC_SHORT: 
			memcpy(s, si+start[0], count[0]*sizeof(short));
			break;
		case NC_INT:
			memcpy(i, ii+start[0], count[0]*sizeof(int));
			break;
		case NC_FLOAT:
			memcpy(f, fi+start[0], count[0]*sizeof(float));
			break;
		case NC_DOUBLE:
			memcpy(d, di+start[0], count[0]*sizeof(double));
			break;
		}
		break;
	case 2:
		i0 = start[0]*dimlen[dimids[1]];
		stridek = dimlen[dimids[1]];
		switch(type) {
		case NC_BYTE: case NC_CHAR:
			for(k=0; k < count[0]; ++k) {
				memcpy(t+offset, ti+i0+k*stridek+start[1], count[1]*sizeof(char));
				offset += count[1];
			}
			break;
		case NC_SHORT:
			for(k=0; k < count[0]; ++k) {
				memcpy(s+offset, si+i0+k*stridek+start[1], count[1]*sizeof(short));
				offset += count[1];
			}
			break;
		case NC_INT:
			for(k=0; k < count[0]; ++k) {
				memcpy(i+offset, ii+i0+k*stridek+start[1], count[1]*sizeof(int));
				offset += count[1];
			}
			break;
		case NC_FLOAT:
			for(k=0; k < count[0]; ++k) {
				memcpy(f+offset, fi+i0+k*stridek+start[1], count[1]*sizeof(float));
				offset += count[1];
			}
			break;
		case NC_DOUBLE:
			for(k=0; k < count[0]; ++k) {
				memcpy(d+offset, di+i0+k*stridek+start[1], count[1]*sizeof(double));
				offset += count[1];
			}
			break;
		}
		break;
	case 3: 
		i0 = start[0]*dimlen[dimids[1]]*dimlen[dimids[2]] + 
			   start[1]*dimlen[dimids[2]] +
			   start[2];
		stridek = dimlen[dimids[1]]*dimlen[dimids[2]];
		stridej = dimlen[dimids[2]];

		switch(type) {
		case NC_BYTE: case NC_CHAR:
			for(k=0; k < count[0]; ++k) {
				for(j=0; j < count[1]; ++j) {
					memcpy(t+offset, ti+i0+k*stridek+j*stridej, count[2]*sizeof(char));
					offset += count[2];
				}
			}
			break;
		case NC_SHORT:
			for(k=0; k < count[0]; ++k) {
				for(j=0; j < count[1]; ++j) {
					memcpy(s+offset, si+i0+k*stridek+j*stridej, count[2]*sizeof(short));
					offset += count[2];
				}
			}
			break;
		case NC_INT:
			for(k=0; k < count[0]; ++k) {
				for(j=0; j < count[1]; ++j) {
					memcpy(i+offset, ii+i0+k*stridek+j*stridej, count[2]*sizeof(int));
					offset += count[2];
				}
			}
			break;
		case NC_FLOAT:
			for(k=0; k < count[0]; ++k) {
				for(j=0; j < count[1]; ++j) {
					memcpy(f+offset, fi+i0+k*stridek+j*stridej, count[2]*sizeof(float));
					offset += count[2];
				}
			}
			break;
		case NC_DOUBLE:
			for(k=0; k < count[0]; ++k) {
				for(j=0; j < count[1]; ++j) {
					memcpy(d+offset, di+i0+k*stridek+j*stridej, count[2]*sizeof(double));
					offset += count[2];
				}
			}
			break;
		}
		break;
	case 4: 
		i0 = start[1]*dimlen[dimids[2]]*dimlen[dimids[3]] + 
			   start[2]*dimlen[dimids[3]] +
			   start[3];
		stridek = dimlen[dimids[2]]*dimlen[dimids[3]];
		stridej = dimlen[dimids[3]];

		switch(type) {
		case NC_BYTE: case NC_CHAR:
			for(k=0; k < count[1]; ++k) {
				for(j=0; j < count[2]; ++j) {
					memcpy(t+offset, ti+i0+k*stridek+j*stridej, count[3]*sizeof(char));
					offset += count[3];
				}
			}
			break;
		case NC_SHORT:
			for(k=0; k < count[1]; ++k) {
				for(j=0; j < count[2]; ++j) {
					memcpy(s+offset, si+i0+k*stridek+j*stridej, count[3]*sizeof(short));
					offset += count[3];
				}
			}
			break;
		case NC_INT:
			for(k=0; k < count[1]; ++k) {
				for(j=0; j < count[2]; ++j) {
					memcpy(i+offset, ii+i0+k*stridek+j*stridej, count[3]*sizeof(int));
					offset += count[3];
				}
			}
			break;
		case NC_FLOAT:
			for(k=0; k < count[1]; ++k) {
				for(j=0; j < count[2]; ++j) {
					memcpy(f+offset, fi+i0+k*stridek+j*stridej, count[3]*sizeof(float));
					offset += count[3];
				}
			}
			break;
		case NC_DOUBLE:
			for(k=0; k < count[1]; ++k) {
				for(j=0; j < count[2]; ++j) {
					memcpy(d+offset, di+i0+k*stridek+j*stridej, count[3]*sizeof(double));
					offset += count[3];
				}
			}
			break;
		}
		break;
	}
}
/*-------------------------------------------------------------------*/
void free_scatter_dims(ScatterDim* dims[], int ndims) {
  int i;
  for(i=0; i < ndims; ++i) {
    ScatterDim_free(dims[i]);
  }
}
/*-------------------------------------------------------------------*/
/* 
Sets each array element to NOSCATTER, SCATTERX, or SCATTERY. 
These tags denote if the dimension should be scattered.
'scatterdims' must be preallocated to size of 'ndims'. 
*/
void get_scatter_dims(int nc, int ndims, int nvars, ScatterDim* scatterdims[], mnsopts *opt)
{
	int dimid, varid, status;
	char name[NC_MAX_NAME];
	char att[256];
	nc_type type;
  size_t dimlen;
	int foundunits = 0;
  scatter_t scatter_type;
  int ndiv;
  int recid; 
  
	status = nc_inq_unlimdim(nc, &recid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error: %s/%d: Failed to query input file for record id.\n", __FILE__, __LINE__);
		exit(-1);
	}

	/* Condition for scatter dimension:
		 Coordvar is x/y based on metadata, and shares dim name. */
	for(dimid=0; dimid < ndims; ++dimid) {
		scatter_type = NOSCATTER;
    dimlen = 0;
    ndiv = 0;
    
		status = nc_inq_dimname(nc, dimid, name);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query number dim name from input file for dimid %d.\n", dimid);
			goto DIMDONE;
		}

    if (dimid == recid) {
			dimlen = NC_UNLIMITED;
      goto DIMDONE;
    } else {
      status = nc_inq_dimlen(nc, dimid, &dimlen);
      if (status != NC_NOERR) {
        fprintf(stderr, "Error. Failed to query dim length from input file for \"%s\".\n", name);
        goto DIMDONE;
      }
    }
    
		/* First check user's x/y dim lists from command-line. */
		if (instringlist(opt->xdims, name, opt->xdims_len)) {
			scatter_type = SCATTERX;
			goto DIMDONE;
		} else if (instringlist(opt->ydims, name, opt->ydims_len)) {
			scatter_type = SCATTERY;
			goto DIMDONE;
		}
		
		/* Try to auto-detect if not in list from command-line. */
		/* Check to see if coordvar is X or Y. */
		status = nc_inq_varid(nc, name, &varid);
		if (status != NC_NOERR) {
			/* No variable that shares dimension name. */
			goto DIMDONE;
		}
		
		/* If not found in user's explicit list of x/y dims, infer via metadata. */
		status = nc_inq_atttype(nc, varid, "units", &type);
		if (status == NC_NOERR) {
			if (type != NC_CHAR) {
				fprintf(stderr, "Error. Coordinate variable units attribute isn't char type for variable %s.  Checking for \"cartesian_axis\" instead.\n", name);
				goto ENDUNITS;
			}
			status = nc_get_att_text(nc, varid, "units", att);
			if (status != NC_NOERR) {
				fprintf(stderr, "Error. Failed to load units attribute for variable %s. Checking for \"cartesian_axis\" instead.\n", name);
				goto ENDUNITS;
			}

			if ( (!strncasecmp(att, "degrees_east",12)) || 
					 (!strncasecmp(att, "degree_east",11)) || 
					 (!strncasecmp(att, "degrees_e",9)) || 
					 (!strncasecmp(att, "degree_e",8)) || 
					 (!strncasecmp(att, "degreee",7)) || 
					 (!strncasecmp(att, "degreese",8)) ) {
				scatter_type = SCATTERX;
				foundunits = 1;
			}	else if	(
								(!strncasecmp(att, "degrees_north",13)) || 
								(!strncasecmp(att, "degree_north",12)) || 
								(!strncasecmp(att, "degrees_n",9)) || 
								(!strncasecmp(att, "degree_n",8)) ||
								(!strncasecmp(att, "degreesn",8)) || 
								(!strncasecmp(att, "degreen",7)) ) {
				scatter_type = SCATTERY;
				foundunits = 1;				
			} else {
				scatter_type = NOSCATTER;
				foundunits = 0;
			}
		}
 
	ENDUNITS:
		if (!foundunits) {
			status = nc_inq_atttype(nc, varid, "cartesian_axis", &type);
			if (status != NC_NOERR) 	goto DIMDONE;
			if (type != NC_CHAR) {
				fprintf(stderr, "Error. Coordinate variable \"cartesian_axis\" attribute isn't char type for variable %s.\n", name);
				exit(-1);
			}
			status = nc_get_att_text(nc, varid, "cartesian_axis", att);
			if (status != NC_NOERR) {
				fprintf(stderr, "Error. Failed to load \"cartesian_axis\" attribute for variable %s.\n", name);
				exit(-1);
			}

			if ( !strncasecmp(att, "x",1) ) 
				scatter_type = SCATTERX;
			else if ( !strncasecmp(att, "y",1) ) 
				scatter_type = SCATTERY;
			else
				scatter_type = NOSCATTER;
		}
		
	DIMDONE:    
    switch(scatter_type) {
      case SCATTERX:
        ndiv = opt->nx;
        break;
      case SCATTERY:
        ndiv = opt->ny;
        break;
      default: 
        ndiv = 0;
        break;
    }
    
    scatterdims[dimid] = ScatterDim_new(dimid, dimlen, name, scatter_type, ndiv);
	}
}
/*-------------------------------------------------------------------*/
/*
Computes the start/end indices per dim if scattered.
*/
void get_scatter_extents(ScatterDim* scatterdims[], int ndims) {
  int i, j;
  ScatterDim * pdim;
  
  for(i=0; i < ndims; ++i) {
    pdim = scatterdims[i];
    if (pdim == NULL) continue;
    if (pdim->scatter_type == NOSCATTER) continue;
    
    mpp_compute_extent(0, pdim->len-1, pdim->scatter_ndiv, pdim->scatter_start, pdim->scatter_end);
    
    /* Compute lengths for scattering. */
    for(j=0; j < pdim->scatter_ndiv; ++j) {
      pdim->scatter_len[j] = pdim->scatter_end[j] - pdim->scatter_start[j] + 1;
    }
  }
}
/*-------------------------------------------------------------------*/
/*
Computes the start/end/len per dim if scattered using optional
io_layout mode, which must equally divide tiling.
*/
void get_scatter_extents_iolayout(ScatterDim* scatterdims[], int ndims, int nxio, int nyio) {
  int d, i, k;
  ScatterDim * pdim;
  size_t *startio = 0;
  size_t *endio = 0;
  size_t *lenio = 0;
  int ndivio;
  int step;
  
  for(d=0; d < ndims; ++d) {
    pdim = scatterdims[d];
    if (pdim == NULL) continue;
    if (pdim->scatter_type == NOSCATTER) continue;

    if (pdim->scatter_type == SCATTERX) {
      ndivio = nxio;
    } else {
      ndivio = nyio;
    }
    
    step = pdim->scatter_ndiv / ndivio;
    
    if (startio) XFREE(startio);
    startio = XMALLOC(size_t, ndivio);

    if (endio) XFREE(endio);
    endio = XMALLOC(size_t, ndivio);

    if (lenio) XFREE(lenio);
    lenio = XMALLOC(size_t, ndivio);
    
    i = k = 0;
    while (i < pdim->scatter_ndiv) {      
      startio[k] = pdim->scatter_start[i];
      endio[k] = pdim->scatter_end[i + step - 1];
      lenio[k] = endio[k] - startio[k] + 1;
      k += 1;
      i += step;
    }
        
    pdim->scatter_ndiv = ndivio;

    if (pdim->scatter_start) XFREE(pdim->scatter_start);
    pdim->scatter_start = startio; /* Own new array. */
    startio = 0;
    
    if (pdim->scatter_end) XFREE(pdim->scatter_end);
    pdim->scatter_end = endio; /* Own new array. */
    endio = 0;
    
    if (pdim->scatter_len) XFREE(pdim->scatter_len);
    pdim->scatter_len = lenio; /* Own new array. */
    lenio = 0;
  }
}
/*-------------------------------------------------------------------*/
void print_scatter_dims(ScatterDim* scatterdims[], int ndims) {
  int d, i;
  ScatterDim * pdim;
  
  for(d=0; d < ndims; ++d) {
    pdim = scatterdims[d];
    if (pdim == NULL) continue;

    fprintf(stdout, "Info: Dimension %d \"%s\":\n", pdim->id, pdim->name);
    fprintf(stdout, "Info:   Will be scattered? %s\n",
          (pdim->scatter_type == NOSCATTER ? "No" : "Yes") );
    
    if (pdim->scatter_type == NOSCATTER) continue;
    
    fprintf(stdout, "Info:   Scatter indices (start, end, length):");
    fflush(stdout);
    
    for(i=0; i < pdim->scatter_ndiv; ++i) {
      fprintf(stdout, " (%zu,%zu,%zu)", pdim->scatter_start[i], pdim->scatter_end[i], pdim->scatter_len[i]);
      fflush(stdout);
    }
    fprintf(stdout, "\n");
  }
}
/*-------------------------------------------------------------------*/
/* Define dimensions with optional scattering in new output files. */
void def_dim(int nc, int *ncids, int ndims, ScatterDim *scatterdims[], mnsopts *opts)
{
	int dimid, xi, yi, i, status, dummy, recid;
	size_t len, newlen;
	char name[NC_MAX_NAME];
  int nx, ny;
  ScatterDim* scatdim = 0;
  
  get_num_divs(opts, &nx, &ny);

	for(dimid=0; dimid < ndims; ++dimid) {
    scatdim = scatterdims[dimid];
    if (scatdim == NULL) continue;
    
		for(yi=0; yi < ny; ++yi) {
			for(xi=0; xi < nx; ++xi) {
				/* Index into complete enumerated list of partitions. Used in new filename. */
				i = yi*nx + xi;
				
				switch(scatdim->scatter_type) {
				case NOSCATTER:
					if (opts->verbose) {
						fprintf(stdout, "Info: Defining dim \"%s\".\n", scatdim->name);
						fflush(stdout);
					}
						
          if (!opts->dryrun) {
            status = nc_def_dim(ncids[i], scatdim->name, scatdim->len, &dummy);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error: %s/%d: Failed to define dimension \"%s\" in output file index %d. Aborting.\n", __FILE__, __LINE__, scatdim->name, i);
              exit(-1);
            }
          }
					break;
					
				case SCATTERX:
          if (!opts->dryrun) {
            status = nc_def_dim(ncids[i], scatdim->name, scatdim->scatter_len[xi], &dummy);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error: %s/%d: Failed to define scattered x dim \"%s\" with length %d in output file index %d. Aborting.\n",  __FILE__, __LINE__, scatdim->name, (int)scatdim->scatter_len[xi], i);
              exit(-1);
            }
          }
					break;
					
				case SCATTERY:
          if (!opts->dryrun) {
            status = nc_def_dim(ncids[i], scatdim->name, scatdim->scatter_len[yi], &dummy);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error: %s/%d: Failed to define scattered y dim \"%s\" with length %d in output file index %d. Aborting.\n",  __FILE__, __LINE__, scatdim->name, (int)scatdim->scatter_len[yi], i);
              exit(-1);
            }
          }				
					break;
				default:
					break;
				}
			}
		}
	}
}
/*-------------------------------------------------------------------*/
void def_var(int nc, int *ncids, int nvars, int ndims, ScatterDim *scatterdims[], mnsopts *opts)
{
	int varid, xi, yi, ifile, j, status, varidnew, natt, ndimvar;
	size_t len, newlen;
	char varname[NC_MAX_NAME], attname[NC_MAX_NAME], dimname[NC_MAX_NAME];
	nc_type type;
	int dimids[NC_MAX_DIMS];
	/* Scatter attribute. */
	int scatatt[4];
  char verbose = opts->verbose;
  char dryrun = opts->dryrun;
  int nx, ny;
  ScatterDim* scatdim = 0;
  
  get_num_divs(opts, &nx, &ny);

  /* 
  From mppnccombine:
     "domain_decomposition = #0, #1, #2, #3 attribute
     #0 starting position of original dimension
     #1 ending position of original dimension
     #2 starting position of decomposed dimension
     #3 ending position of decomposed dimension
  rsz: 
     All values are 1-based.
     #0 is always 1.
     #1 is the original length.
     #2 is 1 based new start.
     #3 is 1 based new length.
  */
	scatatt[0] = 1;

	for(varid=0; varid < nvars; ++varid) {
		status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error: %s/%d: Failed to query var from input file. Aborting.\n", __FILE__, __LINE__);
			exit(-1);
		}

		for(yi=0; yi < ny; ++yi) {
			for(xi=0; xi < nx; ++xi) {
				ifile = yi*nx + xi;

        if (!dryrun) {
          status = nc_def_var(ncids[ifile], varname, type, ndimvar, dimids, &varidnew);
          if (status != NC_NOERR) {
            fprintf(stderr, "Error: %s/%d: Failed to define var \"%s\" in output file index %d. Aborting.\n", __FILE__, __LINE__, varname, ifile+opts->start);
            exit(-1);
          }
        }

				if (verbose) {
					fprintf(stdout, "Info: Defining variable \"%s\" for file %d.\n", varname, ifile+opts->start);
					fflush(stdout);
				}

				/* Copy atts. */
        if (!dryrun) {
          for(j=0; j < natt; ++j) {
            status = nc_inq_attname(nc, varid, j, attname);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error: %s/%d: Failed to define var \"%s\" in output file index %d. Aborting.\n", __FILE__, __LINE__, varname, ifile+opts->start);
              exit(-1);
            }
            
            status = nc_copy_att(nc, varid, attname, ncids[ifile], varidnew);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error: %s/%d: Failed to copy var \"%s\" att \"%s\" to output file index %d. Aborting.\n", __FILE__, __LINE__, varname, attname, ifile+opts->start);
              exit(-1);
            }
          }
        }
				
				/* Need to add new attribute "domain_decomposition" for
				   dimension variables that are 1D. */
				if (ndimvar == 1) {
          scatdim = scatterdims[dimids[0]];
          if (scatdim == NULL) continue;
          
					/* Vars only dimension is scattered and var name is also the dim name. */
					if ( (scatdim->scatter_type == SCATTERX) && 
							 (strcmp(varname, scatdim->name)==0) ) {
						scatatt[1] = (int)scatdim->len;					
  					scatatt[2] = scatdim->scatter_start[xi] + 1; /* 1-based. */
						scatatt[3] = scatdim->scatter_len[xi];

						if (verbose) {
							fprintf(stdout, "Info:   Adding attribute:  %s:domain_decomposition = %d %d %d %d\n", varname, scatatt[0], scatatt[1], scatatt[2], scatatt[3]);
							fflush(stdout);
						}

            if (!dryrun) {
              status = nc_put_att_int(ncids[ifile], varid, "domain_decomposition", NC_INT, 4, scatatt);
              if (status != NC_NOERR) {
                fprintf(stderr, "Error: %s/%d: Failed to set \"domain_decomposition\" attribute data for output file index %d. Aborting.\n", __FILE__, __LINE__, ifile+opts->start);
                exit(-1);
              }
            }
					} else if ( (scatdim->scatter_type == SCATTERY) && 
											(strcmp(varname, scatdim->name)==0) ) {
            scatatt[1] = (int)scatdim->len;					
  					scatatt[2] = scatdim->scatter_start[yi] + 1; 
						scatatt[3] = scatdim->scatter_len[yi];
						
						if (verbose) {
							fprintf(stdout, "Info:   Adding attribute:  %s:domain_decomposition = %d %d %d %d\n", varname, scatatt[0], scatatt[1], scatatt[2], scatatt[3]);
							fflush(stdout);
						}

            if (!dryrun) {
              status = nc_put_att_int(ncids[ifile], varid, "domain_decomposition", NC_INT, 4, scatatt);
              if (status != NC_NOERR) {
                fprintf(stderr, "Error: %s/%d: Failed to set \"domain_decomposition\" attribute data for output file index %d. Aborting.\n", __FILE__, __LINE__, ifile+opts->start);
                exit(-1);
              }
            }
					} 
				} 
			}
		}
	}
}
/*-------------------------------------------------------------------*/
void put_var(int nc, int *ncids, int ndims, int nvars, ScatterDim *scatterdims[], mnsopts *opt)
{
	int varid, xi, yi, i, j, status, varidnew, natt, ndimvar, recid, reci, dimi;
	size_t dimlen[NC_MAX_DIMS], size;
	size_t instart[4];
	size_t outstart[4] = {0,0,0,0};
	size_t count[4], newlen;
	size_t nrec = 0;
  
	size_t maxsize[5] = {0,0,0,0,0}; /* Array of maxsizes for each datatype. */
	enum maxsizeindex {CHAR=0,SHORT,INT,FLOAT,DOUBLE};
  
	char varname[NC_MAX_NAME], attname[NC_MAX_NAME], dimname[NC_MAX_NAME];
	nc_type type;
	int dimids[NC_MAX_DIMS];
	char *tp, *otp;
	short *sp, *osp;
	int *ip, *oip;
	float *fp, *ofp;
	double *dp, *odp;
  ScatterDim *scatdim = 0;
  int nx, ny;

  if (opt->dryrun) return;
  
  get_num_divs(opt, &nx, &ny);
  
	status = nc_inq_unlimdim(nc, &recid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to query input file for record id. Aborting.\n");
		exit(-1);
	}

	if (recid != -1) {
		status = nc_inq_dimlen(nc, recid, &nrec);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to number of records from input file. Aborting.\n");
			exit(-1);
		}
	}

	for(i=0; i < ndims; ++i) {
		status = nc_inq_dimlen(nc, i, &dimlen[i]);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query dim length for dim id %d in input file. Aborting.\n", i);
			exit(-1);
		}
	}

	/* Find largest var to limit alloc to single call because of overhead. */
	for(varid=0; varid < nvars; ++varid) {
		status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query var id %d in input file. Aborting.\n", varid);
			exit(-1);
		}

		size = 1;
		for(i=0; i < ndimvar; ++i) {
			if (dimids[i] == recid)
				continue;

			size *= dimlen[dimids[i]];
		}

		switch(type) {
		case NC_BYTE: case NC_CHAR:
			if (size > maxsize[CHAR])
				maxsize[CHAR] = size;
			break;
		case NC_SHORT:
			if (size > maxsize[SHORT])
				maxsize[SHORT] = size;
			break;
		case NC_INT:
			if (size > maxsize[INT])
				maxsize[INT] = size;
			break;
		case NC_FLOAT:
			if (size > maxsize[FLOAT])
				maxsize[FLOAT] = size;
			break;
		case NC_DOUBLE:
			if (size > maxsize[DOUBLE])
				maxsize[DOUBLE] = size;
			break;
		default:
			fprintf(stderr, "Error. Unknown data type for var %s.\n", varname);
			break;
		}
	}

	/* Allocate all the pointers just once with maximum size. */
	tp = (char*)malloc(sizeof(char)*maxsize[CHAR]);
	otp = (char*)malloc(sizeof(char)*maxsize[CHAR]);
	sp = (short*)malloc(sizeof(short)*maxsize[SHORT]);
	osp = (short*)malloc(sizeof(short)*maxsize[SHORT]);
	ip = (int*)malloc(sizeof(int)*maxsize[INT]);
	oip = (int*)malloc(sizeof(int)*maxsize[INT]);
	fp = (float*)malloc(sizeof(float)*maxsize[FLOAT]);
	ofp = (float*)malloc(sizeof(float)*maxsize[FLOAT]);
	dp = (double*)malloc(sizeof(double)*maxsize[DOUBLE]);
	odp = (double*)malloc(sizeof(double)*maxsize[DOUBLE]);

	/* Process nonrec vars. */
	for(varid=0; varid < nvars; ++varid) {
		status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query var id %d in input file. Aborting.\n", varid);
			exit(-1);
		}
		
		/* Skip record var. */
		if ( (recid != -1) &&
				 (ndimvar > 0) && 
				 (dimids[0] == recid) )
			continue;
		
		if (opt->verbose) {
			fprintf(stdout, "Info: Reading data for static variable \"%s\".\n", varname);
			fflush(stdout);
		}

    switch(type) {
      case NC_BYTE: case NC_CHAR:
        status = nc_get_var_text(nc, varid, tp);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get var %s char data from input file.\n", varname);
          exit(-1);
        }
        break;
      case NC_SHORT:
        status = nc_get_var_short(nc, varid, sp);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get var %s short data from input file.\n", varname);
          exit(-1);
        }
        break;
      case NC_INT:
        status = nc_get_var_int(nc, varid, ip);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get var %s int data from input file.\n", varname);
          exit(-1);
        }
        break;
      case NC_FLOAT:
        status = nc_get_var_float(nc, varid, fp);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get var %s float data from input file.\n", varname);
          exit(-1);
        }
        break;
      case NC_DOUBLE:
        status = nc_get_var_double(nc, varid, dp);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get var %s double data from input file.\n", varname);
          exit(-1);
        }
        break;
      default:
        fprintf(stderr, "Error. Unknown data type for var %s.\n", varname);
        break;
    }

		for(yi=0; yi < ny; ++yi) {
			for(xi=0; xi < nx; ++xi) {
				i = yi*nx + xi;

				if (ndimvar > 0) {
					for(dimi=0; dimi < ndimvar; ++dimi) {
            scatdim = scatterdims[dimids[dimi]];
            if (scatdim == NULL) continue;

						switch(scatdim->scatter_type) {
						case NOSCATTER:
							instart[dimi] = 0;
							count[dimi] = scatdim->len;
							break;
						case SCATTERX:
              instart[dimi] = scatdim->scatter_start[xi];
							count[dimi] = scatdim->scatter_len[xi];
							break;							
						case SCATTERY:
              instart[dimi] = scatdim->scatter_start[yi];
							count[dimi] = scatdim->scatter_len[yi];
							break;
						}
					}
					
					if (opt->verbose) {
            fprintf(stdout, "Info: Performing hyperslab copy into file %d.\n", i+opt->start);
            fprintf(stdout, "Info:   var = \"%s\"\n", varname);
            fprintf(stdout, "Info:   start = ");
            printsizetarray(instart, ndimvar);
            fprintf(stdout, "Info:   count = ");
            printsizetarray(count, ndimvar);
            fprintf(stdout, "Info:   outstart = ");
            printsizetarray(outstart, ndimvar);
            fflush(stdout);
					}

          hyperslabcopy(type, dimlen, dimids, instart, count, ndimvar, tp, sp, ip, fp, dp, otp, osp, oip, ofp, odp); 
					
          switch(type) {
          case NC_BYTE: case NC_CHAR:
            status = nc_put_vara_text(ncids[i], varid, outstart, count, otp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put var \"%s\" char data into output file %d.\n", varname, i);
              handle_error(status);
            }
            break;
          case NC_SHORT:
            status = nc_put_vara_short(ncids[i], varid, outstart, count, osp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put var \"%s\" short data into output file %d.\n", varname, i);
              handle_error(status);
            }
            break;
          case NC_INT:
            status = nc_put_vara_int(ncids[i], varid, outstart, count, oip);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put var \"%s\" int data into output file %d.\n", varname, i);
              handle_error(status);
            }
            break;
          case NC_FLOAT:
            status = nc_put_vara_float(ncids[i], varid, outstart, count, ofp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put var \"%s\" float data into output file %d.\n", varname, i);
              handle_error(status);
            }
            break;
          case NC_DOUBLE:
            status = nc_put_vara_double(ncids[i], varid, outstart, count, odp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put var \"%s\" double data into output file %d.\n", varname, i);
              handle_error(status);
            }
            break;
          }
        } else { /* Scalar variable. */
          switch(type) {
          case NC_BYTE: case NC_CHAR:
            status = nc_put_var_text(ncids[i], varid, tp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put scalar var \"%s\" text data into output file %d.\n", varname, i);
              exit(-1);
            }
            break;
          case NC_SHORT:
            status = nc_put_var_short(ncids[i], varid, sp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put scalar var \"%s\" short data into output file %d.\n", varname, i);
              exit(-1);
            }
            break;
          case NC_INT:
            status = nc_put_var_int(ncids[i], varid, ip);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put scalar var \"%s\" int data into output file %d.\n", varname, i);
              exit(-1);
            }
            break;
          case NC_FLOAT:
            status = nc_put_var_float(ncids[i], varid, fp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put scalar var \"%s\" float data into output file %d.\n", varname, i);
              exit(-1);
            }
            break;
          case NC_DOUBLE:
            status = nc_put_var_double(ncids[i], varid, dp);
            if (status != NC_NOERR) {
              fprintf(stderr, "Error. Failed to put scalar var \"%s\" double data into output file %d.\n", varname, i);
              exit(-1);
            }
            break;
          }
        }
			}
		}
	}

	/* Process recvars. */
	for(reci=0; reci < nrec; ++reci) {
		instart[0] = reci;
		outstart[0] = reci;
		count[0] = 1;
		
		for(varid=0; varid < nvars; ++varid) {
			instart[1] = 0;
			instart[2] = 0;
			instart[3] = 0;
			status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
			if (status != NC_NOERR) {
				fprintf(stderr, "Error. Failed to query var id %d in input file.\n", varid);
				exit(-1);
			}
			
			/* Skip nonrecord var. */
			if ( (dimids[0] != recid) )
				continue;

			if (ndimvar > 1) {
				for(dimi=1; dimi < ndimvar; ++dimi) {
					count[dimi] = dimlen[dimids[dimi]];
				}
			}

			if (opt->verbose) {
				fprintf(stdout, "Info: Reading variable \"%s\", record %d.\n", varname, reci);
			  fflush(stdout);
		  }
			
			/*printsizetarray(instart, ndimvar);
				printsizetarray(count, ndimvar);*/

			switch(type) {
			case NC_BYTE: case NC_CHAR:
				status = nc_get_vara_text(nc, varid, instart, count, tp);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to get var %s char data from input file.\n", varname);
					exit(-1);
				}
				break;
			case NC_SHORT:
				status = nc_get_vara_short(nc, varid, instart, count, sp);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to get var %s short data from input file.\n", varname);
					exit(-1);
				}
				break;
			case NC_INT:
				status = nc_get_vara_int(nc, varid, instart, count, ip);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to get var %s int data from input file.\n", varname);
					exit(-1);
				}
				break;
			case NC_FLOAT:
				status = nc_get_vara_float(nc, varid, instart, count, fp);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to get var %s float data from input file.\n", varname);
					exit(-1);
				}
				break;
			case NC_DOUBLE:
				status = nc_get_vara_double(nc, varid, instart, count, dp);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to get var %s double data from input file.\n", varname);
					exit(-1);
				}
				break;
			default:
				fprintf(stderr, "Error. Unknown data type for var %s.\n", varname);
				break;
			}
			
			for(yi=0; yi < ny; ++yi) {
				for(xi=0; xi < nx; ++xi) {
					i = yi*nx + xi;
					
					for(dimi=1; dimi < ndimvar; ++dimi) {
            scatdim = scatterdims[dimids[dimi]];
            if (scatdim == NULL) continue;

						switch(scatdim->scatter_type) {
						case NOSCATTER:
							instart[dimi] = 0;
							count[dimi] = scatdim->len;
							break;
						case SCATTERX:
              instart[dimi] = scatdim->scatter_start[xi];
							count[dimi] = scatdim->scatter_len[xi];
							break;							
						case SCATTERY:
              instart[dimi] = scatdim->scatter_start[yi];
							count[dimi] = scatdim->scatter_len[yi];
							break;
						}
					}
					
					if (opt->verbose) {
            fprintf(stdout, "Info: Performing hyperslab copy into file %d.\n", i+opt->start);
            fprintf(stdout, "Info:   var = \"%s\"\n", varname);
            fprintf(stdout, "Info:   start = ");
            printsizetarray(instart, ndimvar);
            fprintf(stdout, "Info:   count = ");
            printsizetarray(count, ndimvar);
            fprintf(stdout, "Info:   outstart = ");
            printsizetarray(outstart, ndimvar);
            fflush(stdout);
					}

					hyperslabcopy(type, dimlen, dimids, instart, count, ndimvar == 1?0:ndimvar, tp, sp, ip, fp, dp, otp, osp, oip, ofp, odp); 						

					if (opt->verbose) {
						fprintf(stdout, "Info: Writing variable \"%s\", record %d\n", varname, reci);
						fflush(stdout);
					}
					
					switch(type) {
					case NC_BYTE: case NC_CHAR:
						status = nc_put_vara_text(ncids[i], varid, outstart, count, otp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put var %s data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_SHORT:
						status = nc_put_vara_short(ncids[i], varid, outstart, count, osp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put var %s data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_INT:
						status = nc_put_vara_int(ncids[i], varid, outstart, count, oip);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put var %s data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_FLOAT:
						status = nc_put_vara_float(ncids[i], varid, outstart, count, ofp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put var %s data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_DOUBLE:
						status = nc_put_vara_double(ncids[i], varid, outstart, count, odp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put var %s data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;						
					}
				}
			}
		}
	}
	
	free( tp);
	free(otp);
	free( sp);
	free(osp);
	free( ip);
	free(oip);
	free( fp);
	free(ofp);
	free( dp);
	free(odp);
}
/*-------------------------------------------------------------------*/
void scatter_dims(int nc, int ndims, int nvars, ScatterDim* scatterdims[], mnsopts* opt) {
  /* Populate scatter types and num divisions. */
  get_scatter_dims(nc, ndims, nvars, scatterdims, opt);
  
  /* Compute scatter indices. */
  get_scatter_extents(scatterdims, ndims);
  
  /* Convert to io_layout. */
  if (opt->nxio && opt->nyio) {
    if (opt->nx % opt->nxio) {
      fprintf(stderr, "Error: x divisions are not wholly divisble by io-layout x divisions (%d/%d=%g). Aborting.\n", opt->nx, opt->nxio, ((float)opt->nx)/opt->nxio); 
      exit(1);  
    }

    if (opt->ny % opt->nyio) {
      fprintf(stderr, "Error: y divisions are not wholly divisble by io-layout y divisions (%d/%d=%g). Aborting.\n", opt->ny, opt->nyio, ((float)opt->ny)/opt->nyio); 
      exit(1);  
    }
    
    get_scatter_extents_iolayout(scatterdims, ndims, opt->nxio, opt->nyio);
  }
  
  if (opt->verbose) {
    print_scatter_dims(scatterdims, ndims);
  }
}
/*-------------------------------------------------------------------*/
/*
In:
  isg: Start index for all tiles.
  ieg: End index for all tiles.
  ndivs: Number of tiles.

Out:
  start:  Start indices of symmetric tiling.
          Must be pre-allocated array of size ndivs.
  end:  End indices of symmetric tiling.
        Must be pre-allocated array of size ndivs.
*/
void mpp_compute_extent(size_t isg, size_t ieg, size_t ndivs, size_t* start, size_t* end) {
  size_t n = ieg - isg + 1;
  size_t iss = isg;
  size_t ndiv, imax, ndmax, ie, ndmirror;
  char symmetrize;
  
  for(ndiv=0; ndiv < ndivs; ++ndiv) {
    symmetrize = ( EVEN(ndivs) && EVEN(n) ) ||
      ( ODD(ndivs) && ODD(n) ) ||
      ( ODD(ndivs) && EVEN(n) && (ndivs < (n/2)) );
      
    if (ndiv == 0) {
      imax = ieg;
      ndmax = ndivs;
    }
    
    if ( ndiv < ((ndivs-1)/2+1) ) {
      ie = iss + ceil( ((float)(imax-iss+1.0))/(ndmax-ndiv) ) - 1;
      ndmirror = (ndivs-1) - ndiv;
      if ( (ndmirror > ndiv) && symmetrize) {
        start[ndmirror] = MAX( isg+ieg-ie, ie+1 );
        end[ndmirror]   = MAX( isg+ieg-iss, ie+1 );
        imax = start[ndmirror] - 1;
        ndmax = ndmax - 1;
      }
    } else {
      if (symmetrize) {
        iss = start[ndiv];
        ie = end[ndiv];
      } else {
        ie = iss + ceil( ((float)(imax-iss+1.0))/(ndmax-ndiv) ) - 1;
      }
    }

    start[ndiv] = iss;
    end[ndiv] = ie;

    if (ie < iss) {
      fprintf(stderr, "Error: %s/%d: domain extents must be positive definite. \"ie\"=%zu, \"iss\"=%zu\n", __FILE__, __LINE__, ie, iss);
    }
    if ( (ndiv == (ndivs-1)) && (end[ndiv] != ieg) ) {
      fprintf(stderr, "Error: %s/%d: domain extents do not span space completely.\n", __FILE__, __LINE__);
    }
    
    iss = ie + 1;
  }
}
/*-------------------------------------------------------------------*/
int mppncscatter(mnsopts* opts) 
{
	int status = 0;
	int nfiles = 0; /* number of output files. */
	int nc; /* input file id. */
	int *ncids = NULL; /* store ncid for output files. */
	int i,j,k;
	char output[256]; /* Temporary buffer for creating output filenames. */
	char *prefix, *pchar; /* pointer to root file name (without path). */
	int format = 0; /* Create output files with same format as inputs. */
	int dummy;
	int natts, ngatts, ndims, nvars, unlimdimid;
	nc_type type;
	size_t len;
  ScatterDim * scatterdims[NC_MAX_DIMS];
	char name[NC_MAX_NAME];
  char outnameformat[256]; /* Format string for creating out filenames. */
  char dryrun = opts->dryrun;
  char verbose = opts->verbose;
  
  /*-------------------------------------------------*/
	/* Strip path in file name for output file names. */
	prefix = opts->filein;
	pchar = strstr(prefix, "/");
	if (pchar != NULL) {
		do {
			prefix = pchar + 1;
		} while(pchar = strstr(prefix, "/"));
	}

  /*-------------------------------------------------*/
  /* Get basic input file info. */
	status = nc_open(opts->filein, NC_NOWRITE, &nc);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error: %s/%d: Failed to open input file \"%s\". Aborting.\n", __FILE__, __LINE__, opts->filein);
		return -1;
	}

	status = nc_inq(nc, &ndims, &nvars, &ngatts, &unlimdimid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error: %s/%d: Failed to query general info for input file. Aborting.\n", __FILE__, __LINE__);
		return -1;
	}

	status = nc_inq_format(nc, &format);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error: %s/%d: Failed to query input file format. Aborting.\n", __FILE__, __LINE__);
		return -1;
	}

  switch(format) {
    case NC_FORMAT_64BIT:
      format = NC_64BIT_OFFSET;
      break;
    case NC_FORMAT_NETCDF4:
      format = NC_NETCDF4;
      break;
    case NC_FORMAT_NETCDF4_CLASSIC:
      format = NC_CLASSIC_MODEL | NC_NETCDF4;
      break;
    case NC_FORMAT_CLASSIC:
    default:
      format = 0;
      break;
  }
  
  /*-------------------------------------------------*/
  if (opts->prefix) 
    i = strlen(opts->prefix);
  else
    i = 0;
  
  sprintf(outnameformat, "%s%s%%s.%%0%dd", ((i>0) ? opts->prefix : ""), ((i>0) ? "/" : ""), opts->width);
  
  /*-------------------------------------------------*/
  scatter_dims(nc, ndims, nvars, scatterdims, opts);
  
  nfiles = get_num_files(opts);
	ncids = XMALLOC(int, nfiles);
  
	for(i=0; i < nfiles; ++i) {
		if (sprintf(output, outnameformat, prefix, i+opts->start) < 1) {
			fprintf(stderr, "Error: %s/%d: Failed to create output file name. Aborting.\n", __FILE__, __LINE__);
			return -1;
		}
    
    if (verbose) {
      fprintf(stdout, "Info: Creating file \"%s\".\n", output);
      fflush(stdout);
    }
      
    if (!dryrun) {
      status = nc_create(output, NC_CLOBBER | format, &ncids[i]);
      if (status != NC_NOERR) {
        fprintf(stderr, "Error: %s/%d: Failed to create output file \"%s\". Aborting.\n",  __FILE__, __LINE__, output);
        return -1;
      }
    
      status = nc_set_fill(ncids[i], NC_NOFILL, &dummy);
      if (status != NC_NOERR) {
        fprintf(stderr, "Error: %s/%d: Failed to disable prefill for output file \"%s\". Aborting.\n",  __FILE__, __LINE__, output);
        return -1;
      }

      status = nc_put_att_int(ncids[i], NC_GLOBAL, "NumFilesInSet", NC_INT, 1, &nfiles);
      if (status != NC_NOERR) {
        fprintf(stderr, "Error: %s/%d: Failed to add global attribute \"NumFilesInSet\" to output file \"%s\". Aborting.\n", __FILE__, __LINE__, output);
        return -1;
      }
    }
	}	
 
  /*-------------------------------------------------------*/
	/* Copy all global attributes. */
  if (!dryrun) {
    for(j=0; j < ngatts; ++j) {
      for(i=0; i < nfiles; ++i) {
        status = nc_inq_attname(nc, NC_GLOBAL, j, name);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to get global attribute name from input file.\n");
          return -1;
        }

        status = nc_copy_att(nc, NC_GLOBAL, name, ncids[i], NC_GLOBAL);
        if (status != NC_NOERR) {
          fprintf(stderr, "Error. Failed to copy global attribute to output file %04d.\n", i);
          return -1;
        }
      }
    }
  }
  
  if (verbose) {
		fprintf(stdout, "Info: Defining output dimensions.\n");
		fflush(stdout);
	}
	
	def_dim(nc, ncids, ndims, scatterdims, opts);

	if (verbose) {
		fprintf(stdout, "Info: Defining output variables.\n");
		fflush(stdout);
	}
	
	def_var(nc, ncids, nvars, ndims, scatterdims, opts);

	if (verbose) {
		fprintf(stdout, "Info: Ending define mode.\n");
		fflush(stdout);
	}
	
  if (!dryrun) {
    for(i=0; i < nfiles; ++i)
      nc_enddef(ncids[i]);
  }
  
	put_var(nc, ncids, ndims, nvars, scatterdims, opts);

	if (verbose) {
		fprintf(stdout, "Info: Closing files.\n");
		fflush(stdout);
	}
	
	nc_close(nc);
  
  if (!dryrun) {
    for(i=0; i < nfiles; ++i)
      nc_close(ncids[i]);
  }

	if (verbose) {
		fprintf(stdout, "Info: Freeing memory.\n");
		fflush(stdout);
	}

	if (ncids != NULL)
		free(ncids);

  free_scatter_dims(scatterdims, ndims);

	return 0;
}


