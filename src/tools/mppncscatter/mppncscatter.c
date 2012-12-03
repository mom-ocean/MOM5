/*
   Copyright (C) 2007,2009-2010,2012 Remik Ziemlinski

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
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <netcdf.h>
#include <string.h>
#include <strings.h>
#include "strlist.h"
#include "getopt.h"

#define USAGE "\
mppncscatter -- Decomposes single NetCDF file into many files (converse to mppnccombine).  The output files are created in the current working directory and prefixed with the input file name, i.e. in.nc.0000 ...\n\
\n\
Usage: mppncscatter [OPTION...] in.nc\n\
\n\
  -h, --help         Give this usage message.\n\
      --usage                \n\
  -s, --start N      Start file name suffix numbers from N (default is 0).\n\
  -v, --version      Print program version.\n\
  -V, --verbose      Progressively output messages to stdout.\n\
  -x, --npx N        Try to split domain evenly into N columns.\n\
  -X, --xdims d1,... List of X dimension names to scatter (for those not detectable through metadata).\n\
  -y, --npy N        Try to split domain evenly into N rows.\n\
  -Y, --ydims d1,... List of Y dimension names to scatter (for those not detectable through metadata).\n\
\n\
Report bugs to Remik . Ziemlinski @ noaa gov.\n\
"

typedef struct MNSOPTS {
	char*   filein;		/* Input filename (allocated). */
	int     help;     /* give usage insructions */
	int 		nx;				/* Number of columns to split file into. (Required) */
	int 		ny;				/* Number of rows to split file into. (Required) */
	int     start;    /* Start filename number suffix from this. */
	int     version;  /* give program version */
	int     verbose;  /* Verbose echos to stdout. */
	char**	xdims;		/* List of xdim names to scatter. */
	int			xdims_len;/* Number of names in above list. */
	char**	ydims;		/* List of xdim names to scatter. */
	int			ydims_len;/* Number of names in above list. */
} mnsopts;

#define NOSCATTER 0
#define SCATTERX 1
#define SCATTERY 2

int getmnsopts(int argc, char** argv, mnsopts* popts);
void printusage();
void printversion();
void initmnsopts(mnsopts* popts);
void freemnsopts(mnsopts* popts);
int mppncscatter(mnsopts* popts);

#define handle_error(status) {                      \
	if (status != NC_NOERR) {                         \
		fprintf(stderr, "%s\n", nc_strerror(status));   \
		exit(-1);                                       \
	}                                                 \
}    

void printsizetarray(size_t *a, int n) {
	int i=0;
	for(;i < n; ++i)
		fprintf(stdout, "%d ", (int)a[i]);

	fprintf(stdout, "\n");
}
static struct option const long_options[] =
{ 
	{"help",        no_argument,            0, 'h'},
	{"usage",       no_argument,            0, 'h'},
	{"start" ,      required_argument,      0, 's'},
	{"version",     no_argument,            0, 'v'},
	{"verbose",     no_argument,            0, 'V'},
	{"npx" ,        required_argument,      0, 'x'},
	{"xdims",     	required_argument,      0, 'X'},
	{"npy" ,        required_argument,      0, 'y'},
	{"ydims",     	required_argument,      0, 'Y'},
	{0, 0, 0, 0}
}; 

void initmnsopts(mnsopts* popts)
{
	if (popts == NULL) 
		return;

	popts->filein = NULL;
	popts->help = 0;
	popts->nx = 0;
	popts->ny = 0;
	popts->start = 0;
	popts->verbose = 0;
	popts->version = 0;
	popts->xdims = NULL;
	popts->xdims_len = 0;
	popts->ydims = NULL;
	popts->ydims_len = 0;
}
void freemnsopts(mnsopts* popts)
{
	if (popts == NULL) return;

	if (popts->filein != NULL) {
		free(popts->filein);
		popts->filein = NULL;
	}

	if (popts->xdims != NULL) {
		freestringlist(&popts->xdims, NC_MAX_DIMS);
		popts->xdims = NULL;
	}

	if (popts->ydims != NULL) {
		freestringlist(&popts->ydims, NC_MAX_DIMS);
		popts->ydims = NULL;
	}
}
void printusage()
{
	fprintf(stderr, USAGE);
	fprintf(stderr, "Built with NetCDF %s\n", nc_inq_libvers());
}
void printversion()
{
	fprintf(stderr, "mppncscatter 0.2.0\n");
	fprintf(stderr, "Built with NetCDF %s\n", nc_inq_libvers());
	fprintf(stderr, "Copyright (C) 2007,2009-2010 Remik Ziemlinski\n\
\n\
This program comes with NO WARRANTY, to the extent permitted by law.\n\
You may redistribute copies of this program\n\
under the terms of the GNU General Public License.\n"); 
}
int getmnsopts(int argc, char** argv, mnsopts* popts)
{
	int c;       
	char *token, *cp;
  const char delimiters[] = ",";

	if (popts == NULL) return -1;
  
	if (newstringlist(&popts->xdims, &c, NC_MAX_DIMS)) {
		printf("ERROR: Failed to allocate memory for X dimension list.\n");
		exit(1);
	}

	if (newstringlist(&popts->ydims, &c, NC_MAX_DIMS)) {
		printf("ERROR: Failed to allocate memory for Y dimension list.\n");
		exit(1);
	}
	
	while ( (c = getopt_long(argc, argv, "hs:vVx:y:X:Y:", long_options, 0))
					!= -1
					)
		switch (c)
			{
			case 's':
				popts->start = atoi(optarg);
				break;

			case 'x':
				popts->nx = atoi(optarg);
				break;

			case 'y':
				popts->ny = atoi(optarg);
				break;

			case 'v':
				popts->version = 1;
				break;

			case 'V':
				popts->verbose = 1;
				break;
			
			case ':':
				fprintf(stderr, "Error, -%c without argument.\n\n", optopt); 
				popts->help = 1;
				break;
			case '?':
				fprintf(stderr, "Error, Unknown argument %c.\n\n", optopt);
				popts->help = 1;
				break;
			case 'h':
				popts->help = 1;
				break;

			case 'X':
				getstringlist(optarg, &popts->xdims, &popts->xdims_len);
				break;

			case 'Y':
				getstringlist(optarg, &popts->ydims, &popts->ydims_len);
				break;
			}

	if (popts->help)		{
		printusage();
		return -1;
	}       
	if (popts->version)		{
		printversion();
		return -1;
	}
	if (optind == argc)        {
		fprintf(stderr, "Error, missing operand after `%s'.\n\n", argv[argc - 1]);
		printusage();
		return -1;
	} 

	/* get filename argument */   
	argc -= optind;
	argv += optind;
	if ((argc < 1) || (argv[0] == NULL))        {
		fprintf(stderr, "Error, input file required.\n\n");
		printusage();
		popts->filein = NULL;
		return -1;
	} else        {
		/* store filename */
		popts->filein = (char*)malloc(sizeof(char)*(strlen(argv[0]) + 1));
		strcpy(popts->filein, argv[0]);
	}

	return 0;
}
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
//----------------------------------------------------------------------------
// Returns dimension size for an even partitioning.
// In:
//  i: Partition index (0 <= i < n).
//  len: Entire length of original dimension.
//  n: Number of partitions for the dimension.
size_t dimlen_even(int i, size_t len, int n) {
	size_t newlen = (size_t)(len/n); 

  if ( i == (n-1) ) 
    // Last column is remainder. 
    newlen = len - i*newlen; 
		
  return newlen;
		
	/* NOT PRODUCTION READY.
	// Is staggered?  Assume yes for dim that has odd integer size. 
	if (len % 2) { 
		// Define size that will allow boundary duplication on contact edges. 
		if (i != (n-1)) { // If this isn't the last partition...
			// Duplicate ending edge of parts except for last partition. 
			newlen += 1; 
		} else { 
			// Inner partition that shares edges with left and right neighbors. 
			newlen = len - i*newlen; 
		} 
	} else { 
		if ( i == (n-1) ) 
			// Last column is remainder. 
			newlen = len - i*newlen; 
	}
	
	return newlen;
	*/
}
//----------------------------------------------------------------------------
/* scatterdims must be preallocated.  Sets each array element to
	 NOSCATTER, SCATTERX, or SCATTERY. These tags denote if the dimension
	 should be scattered.
*/
void getscatterdims(int nc, int ndims, int nvars, int *scatterdims, mnsopts *opt)
{
	int dimid, varid, status;
	char name[NC_MAX_NAME];
	char att[256];
	nc_type type;
	int foundunits = 0;

	/* Condition for scatter dimension:
		 Coordvar is x/y based on metadata, and shares dim name. */
	for(dimid=0; dimid < ndims; ++dimid) {
		scatterdims[dimid] = NOSCATTER;
		
		status = nc_inq_dimname(nc, dimid, name);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query number dim name from input file.\n");
			exit(-1);
		}
		
		// First check user's x/y dim lists from command-line.
		if (instringlist(opt->xdims, name, opt->xdims_len)) {
			scatterdims[dimid] = SCATTERX;
			goto DIMDONE;
		} else if (instringlist(opt->ydims, name, opt->ydims_len)) {
			scatterdims[dimid] = SCATTERY;
			goto DIMDONE;
		}
		
		// Try to auto-detect if not in list from command-line.
		// Check to see if coordvar is X or Y. 
		status = nc_inq_varid(nc, name, &varid);
		if (status != NC_NOERR) {
			/* No variable that shares dimension name. */
			scatterdims[dimid] = NOSCATTER;
			goto DIMDONE;
		}
		
		// If not found in user's explicit list of x/y dims, infer via metadata. 
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
				scatterdims[dimid] = SCATTERX;
				foundunits = 1;
			}	else if	(
								(!strncasecmp(att, "degrees_north",13)) || 
								(!strncasecmp(att, "degree_north",12)) || 
								(!strncasecmp(att, "degrees_n",9)) || 
								(!strncasecmp(att, "degree_n",8)) ||
								(!strncasecmp(att, "degreesn",8)) || 
								(!strncasecmp(att, "degreen",7)) ) {
				scatterdims[dimid] = SCATTERY;
				foundunits = 1;				
			} else {
				scatterdims[dimid] = NOSCATTER;
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
				scatterdims[dimid] = SCATTERX;
			else if ( !strncasecmp(att, "y",1) ) 
				scatterdims[dimid] = SCATTERY;
			else
				scatterdims[dimid] = NOSCATTER;
		}
		
	DIMDONE:
		if (opt->verbose)
			fprintf(stdout, "Dimension \"%s\" will be scattered? %s\n", name,
							scatterdims[dimid] != NOSCATTER ? "Yes" : "No");
	}
}
/* Define dimensions with optional scattering in new output files. */
void defdim_even(int nc, int *ncids, int ndims, int *scatterdims, mnsopts *opt)
{
	int dimid, xi, yi, i, status, dummy, recid;
	size_t len, newlen;
	char name[NC_MAX_NAME];

	status = nc_inq_unlimdim(nc, &recid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to query input file for record id.\n");
		exit(-1);
	}

	for(dimid=0; dimid < ndims; ++dimid) {
		status = nc_inq_dim(nc, dimid, name, &len);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to get dim name and length from input file.\n");
			exit(-1);
		}
		
		if (dimid == recid)
			len = NC_UNLIMITED;

		for(yi=0; yi < opt->ny; ++yi) {
			for(xi=0; xi < opt->nx; ++xi) {
				// Index into complete enumerated list of partitions. Used in new filename.
				i = yi*opt->nx + xi;
				
				switch(scatterdims[dimid]) {
				case NOSCATTER:
					if (opt->verbose) 
						fprintf(stderr, "INFO : %d : DEFDIM ncid = %d, dim = %s, len = %d\n", __LINE__, ncids[i], name, (int)len);
						
					status = nc_def_dim(ncids[i], name, len, &dummy);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to define dimension %s in output file %04d.\n", name, i);
						exit(-1);
					}
					break;
					
				case SCATTERX:
					newlen = dimlen_even(xi, len, opt->nx);
					status = nc_def_dim(ncids[i], name, newlen, &dummy);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to define scattered x dim \"%s\" with length %d in output file %04d.\n", name, (int)newlen, i);
						exit(-1);
					}
					
					if (opt->verbose) 
						fprintf(stderr, "INFO : %d : DEFDIM ncid = %d, dim = %s, len = %d\n", __LINE__, ncids[i], name, (int)newlen);
						
					break;
					
				case SCATTERY:
					newlen = dimlen_even(yi, len, opt->ny);
					status = nc_def_dim(ncids[i], name, newlen, &dummy);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to define scattered y dim with length %d in output file %04d.\n", (int)newlen, i);
						exit(-1);
					}
					
					if (opt->verbose) 
						fprintf(stderr, "INFO : %d : DEFDIM ncid = %d, dim = %s, len = %d\n", __LINE__, ncids[i], name, (int)newlen);
					
					break;
				default:
					break;
				}
			}
		}
	}
}
void defvar_even(int nc, int *ncids, int nvars, int ndims, int *scatterdims, mnsopts *opt)
{
	int varid, xi, yi, i, j, status, varidnew, natt, ndimvar;
	size_t len, newlen;
	char varname[NC_MAX_NAME], attname[NC_MAX_NAME], dimname[NC_MAX_NAME];
	nc_type type;
	int dimids[NC_MAX_DIMS];
	// Scatter attribute.
	int attdata[4];

	if (opt->verbose)	
		fprintf(stdout, "DEBUG : %d\n", __LINE__);
							
	/* The first element in the scatter attribute is always 1. */
	attdata[0] = 1;

	for(varid=0; varid < nvars; ++varid) {
		status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query var from input file.\n");
			exit(-1);
		}

		for(yi=0; yi < opt->ny; ++yi) {
			for(xi=0; xi < opt->nx; ++xi) {
				i = yi*opt->nx + xi;

				status = nc_def_var(ncids[i], varname, type, ndimvar, dimids, &varidnew);
				if (status != NC_NOERR) {
					fprintf(stderr, "Error. Failed to define var %s in output file %d.\n", varname, i);
					exit(-1);
				}

				if (opt->verbose)
					fprintf(stdout, "DEBUG : %d : Defining variable %s\n", __LINE__, varname);

				/* Copy atts. */
				for(j=0; j < natt; ++j) {
					status = nc_inq_attname(nc, varid, j, attname);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to get name for att id %d.\n", j);
						exit(-1);
					}
					
					status = nc_copy_att(nc, varid, attname, ncids[i], varidnew);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to copy att %s to output file %d.\n", varname, i);
						exit(-1);
					}
				}
				
				// Need to add new attribute "domain_decomposition" for
				// dimension variables that are 1D.
				if (ndimvar == 1) {
					if (opt->verbose)
						fprintf(stdout, "DEBUG : %d :\tscatterdims[dimids[0]] = %d\n", __LINE__, scatterdims[dimids[0]]);		
									
					status = nc_inq_dimname(nc, dimids[0], dimname);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to get dim name from input file for dim id %d.\n", dimids[0]);
						exit(-1);
					}
					
					status = nc_inq_dimlen(nc, dimids[0], &len);
					if (status != NC_NOERR) {
						fprintf(stderr, "Error. Failed to get dim length from input file.\n");
						exit(-1);
					}

					// Vars only dimension is scattered and var name is also the dim name.
					if ( (scatterdims[dimids[0]] == SCATTERX) && 
							 (strcmp(varname, dimname)==0) ) {
						attdata[1] = (int)len;
						
						if (xi == 0) {
							attdata[2] = 1;
						} else {
							newlen = dimlen_even(xi-1, len, opt->nx);
							attdata[2] = (int)(xi * newlen) + 1; // 1-based syntax. 
						}
						
						newlen = dimlen_even(xi, len, opt->nx);
						attdata[3] = (int)newlen;

						if (opt->verbose) {
							fprintf(stdout, "Adding attribute:\t%s:domain_decomposition = %d %d %d %d\n", varname, attdata[0], attdata[1], attdata[2], attdata[3]);
						}

						status = nc_put_att_int(ncids[i], varid, "domain_decomposition", NC_INT, 4, attdata);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to set domain_decomposition attribute data for output file %d.\n", i);
							exit(-1);
						}								
					} else if ( (scatterdims[dimids[0]] == SCATTERY) && 
											(strcmp(varname, dimname)==0) ) {
						attdata[1] = (int)len;
						
						if (yi == 0) {
							attdata[2] = 1;
						} else {
							newlen = dimlen_even(yi-1, len, opt->ny);
							attdata[2] = (int)(yi * newlen) + 1; // 1-based syntax.
						}
						
						newlen = dimlen_even(yi, len, opt->ny);
						attdata[3] = (int)newlen;
						
            if (opt->verbose) {
							fprintf(stdout, "Adding attribute:\t%s:domain_decomposition = %d %d %d %d\n", varname, attdata[0], attdata[1], attdata[2], attdata[3]);
						}

						status = nc_put_att_int(ncids[i], varid, "domain_decomposition", NC_INT, 4, attdata);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to set domain_decomposition attribute data for variable %s output file %d.\n", varname, i);
							exit(-1);
						}
					} 
				} 
			}
		}
	}
}
void putvar_even(int nc, int *ncids, int ndims, int nvars, int *scatterdims, mnsopts *opt)
{
	int varid, xi, yi, i, j, status, varidnew, natt, ndimvar, recid, reci, dimi;
	size_t dimlen[NC_MAX_DIMS], size;
	size_t instart[4];
	size_t outstart[4] = {0,0,0,0};
	size_t count[4], newlen;
	size_t nrec = 0;
	size_t maxsize[5] = {0,0,0,0,0}; /* Array of maxsizes for each datatype. */
	enum maxsizeindex {CHAR=0,SHORT,INT,FLOAT,DOUBLE} ;
	char varname[NC_MAX_NAME], attname[NC_MAX_NAME], dimname[NC_MAX_NAME];
	nc_type type;
	int dimids[NC_MAX_DIMS];
	char *tp, *otp;
	short *sp, *osp;
	int *ip, *oip;
	float *fp, *ofp;
	double *dp, *odp;

	status = nc_inq_unlimdim(nc, &recid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to query input file for record id.\n");
		exit(-1);
	}

	if (recid != -1) {
		status = nc_inq_dimlen(nc, recid, &nrec);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to number of records from input file.\n");
			exit(-1);
		}
	}

	for(i=0; i < ndims; ++i) {
		status = nc_inq_dimlen(nc, i, &dimlen[i]);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query dim length for dim id %d in input file.\n", i);
			exit(-1);
		}
	}

	/* Find largest var to limit alloc to single call because of overhead. */
	for(varid=0; varid < nvars; ++varid) {
		status = nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to query var id %d in input file.\n", varid);
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
			fprintf(stderr, "Error. Failed to query var id %d in input file.\n", varid);
			exit(-1);
		}
		
		/* Skip record var. */
		if ( (recid != -1) &&
				 (ndimvar > 0) && 
				 (dimids[0] == recid) )
			continue;
		
		if (opt->verbose) {
			fprintf(stdout, "Reading data for static variable %s.\n", varname);
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

		for(yi=0; yi < opt->ny; ++yi) {
			for(xi=0; xi < opt->nx; ++xi) {
				i = yi*opt->nx + xi;

				if (ndimvar > 0) {
					for(dimi=0; dimi < ndimvar; ++dimi) {
						switch(scatterdims[dimids[dimi]]) {
						case NOSCATTER:
							instart[dimi] = 0;
							count[dimi] = dimlen[dimids[dimi]];
							break;
						case SCATTERX:
							if (xi == 0) {
								newlen = dimlen_even(xi, dimlen[dimids[dimi]], opt->nx);
								instart[dimi] = 0;
							} else {
								newlen = dimlen_even(xi-1, dimlen[dimids[dimi]], opt->nx);
								instart[dimi] = newlen*xi;
								
							}
							count[dimi] = newlen;
							break;
							
						case SCATTERY:
							if (yi == 0) {
								newlen = dimlen_even(yi, dimlen[dimids[dimi]], opt->ny);
								instart[dimi] = 0;
							} else {
								newlen = dimlen_even(yi-1, dimlen[dimids[dimi]], opt->ny);
								instart[dimi] = newlen*yi;
								
							}
							count[dimi] = newlen;
							break;
						}
					}
					
					if (opt->verbose) {
            fprintf(stdout, "\tvar = %s\n", varname);
            fprintf(stdout, "\tstart = ");
            printsizetarray(instart, ndimvar);
            fprintf(stdout, "\tcount = ");
            printsizetarray(count, ndimvar);
            fprintf(stdout, "\toutstart = ");
            printsizetarray(outstart, ndimvar);				
            fprintf(stdout, "Performing hyperslab copy into tile %d.\n", i);
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
				} else {
					/* Scalar variable. */
					switch(type) {
					case NC_BYTE: case NC_CHAR:
						status = nc_put_var_text(ncids[i], varid, tp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put scalar var %s text data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_SHORT:
						status = nc_put_var_short(ncids[i], varid, sp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put scalar var %s short data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_INT:
						status = nc_put_var_int(ncids[i], varid, ip);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put scalar var %s int data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_FLOAT:
						status = nc_put_var_float(ncids[i], varid, fp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put scalar var %s float data into output file %d.\n", varname, i);
							exit(-1);
						}
						break;
					case NC_DOUBLE:
						status = nc_put_var_double(ncids[i], varid, dp);
						if (status != NC_NOERR) {
							fprintf(stderr, "Error. Failed to put scalar var %s double data into output file %d.\n", varname, i);
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

			if (opt->verbose)
				fprintf(stdout, "Reading variable %s, record %d.\n", varname, reci);
			
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
			
			for(yi=0; yi < opt->ny; ++yi) {
				for(xi=0; xi < opt->nx; ++xi) {
					i = yi*opt->nx + xi;
					
					for(dimi=1; dimi < ndimvar; ++dimi) {
						switch(scatterdims[dimids[dimi]]) {
						case NOSCATTER:
							instart[dimi] = 0;
							count[dimi] = dimlen[dimids[dimi]];
							break;
						case SCATTERX:
							if (xi == 0) {
								newlen = dimlen_even(xi, dimlen[dimids[dimi]], opt->nx);
								instart[dimi] = 0;
							} else {
								newlen = dimlen_even(xi-1, dimlen[dimids[dimi]], opt->nx);
								instart[dimi] = newlen*xi;
								
							}
							count[dimi] = newlen;
							break;
							
						case SCATTERY:
							if (yi == 0) {
								newlen = dimlen_even(yi, dimlen[dimids[dimi]], opt->ny);
								instart[dimi] = 0;
							} else {
								newlen = dimlen_even(yi-1, dimlen[dimids[dimi]], opt->ny);
								instart[dimi] = newlen*yi;
								
							}
							count[dimi] = newlen;
							break;

						}
					}
					
					if (opt->verbose)
						fprintf(stdout, "Dicing variable %s, record %d\n", varname, reci);

					hyperslabcopy(type, dimlen, dimids, instart, count, ndimvar == 1?0:ndimvar, tp, sp, ip, fp, dp, otp, osp, oip, ofp, odp); 						

					if (opt->verbose)
						fprintf(stdout, "Writing variable %s, record %d\n", varname, reci);
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
int mppncscatter(mnsopts* popts) 
{
	int status = 0;
	int nfiles = popts->nx * popts->ny; 
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
	int scatterdims[NC_MAX_DIMS];
	int *scatterlenx = NULL;
	int *scatterleny = NULL; /* Stores subdomain dim sizes for each output file. */
	char name[NC_MAX_NAME];

	ncids = (int*)malloc(sizeof(int)*nfiles);

	/* Strip path in file name for output file names. */
	prefix = popts->filein;
	pchar = strstr(prefix, "/");
	if (pchar != NULL) {
		do {
			prefix = pchar + 1;
		} while(pchar = strstr(prefix, "/"));
	}

	status = nc_open(popts->filein, NC_NOWRITE, &nc);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to open input file.\n");
		return -1;
	}

	status = nc_inq(nc, &ndims, &nvars, &ngatts, &unlimdimid);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to query general info for input file.\n");
		return -1;
	}

	status = nc_inq_format(nc, &format);
	if (status != NC_NOERR) {
		fprintf(stderr, "Error. Failed to query input file format.\n");
		return -1;
	}

	for(i=0; i < nfiles; ++i) {
		if (sprintf(output, "%s.%04d", prefix, i+popts->start) < 1) {
			fprintf(stderr, "Error. Failed to create output file name.\n");
			return -1;
		}
		status = nc_create(output, NC_CLOBBER | format, &ncids[i]);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to create output file %s.\n", output);
			return -1;
		}
		status = nc_set_fill(ncids[i], NC_NOFILL, &dummy);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to disable prefill for output file %s.\n", output);
			return -1;
		}

		status = nc_put_att_int(ncids[i], NC_GLOBAL, "NumFilesInSet", NC_INT, 1, &nfiles);
		if (status != NC_NOERR) {
			fprintf(stderr, "Error. Failed to add global attribute \"NumFilesInSet\" to output file %s.\n", output);
			return -1;
		}
	}	
 
	/* Copy all global attributes. */
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

	getscatterdims(nc, ndims, nvars, scatterdims, popts);

	if (popts->verbose)
		fprintf(stdout, "Defining output dimensions.\n");

	defdim_even(nc, ncids, ndims, scatterdims, popts);

	if (popts->verbose)
		fprintf(stdout, "Defining output variables.\n");

	defvar_even(nc, ncids, nvars, ndims, scatterdims, popts);

	if (popts->verbose)
		fprintf(stdout, "Ending define mode.\n");

	for(i=0; i < nfiles; ++i)
		nc_enddef(ncids[i]);

	putvar_even(nc, ncids, ndims, nvars, scatterdims, popts);

	nc_close(nc);
	for(i=0; i < nfiles; ++i)
		nc_close(ncids[i]);

	if (ncids != NULL)
		free(ncids);

	return 0;
}
int
main(int argc, char** argv)
{
	int status;
	mnsopts opts;

	status = 0;

	initmnsopts(&opts);
	
	/* parse command-line args.  */
	status = getmnsopts(argc, argv, &opts);
	
	status = status ? status : mppncscatter(&opts);               
	
	freemnsopts(&opts);
	
	return status;
} 

