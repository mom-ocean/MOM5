#ifndef MPPNCSCATTER_OPT_H
#define MPPNCSCATTER_OPT_H

#include "getopt.h"
#include "common.h"

/*-------------------------------------------------------------------*/
#define USAGE "\
mppncscatter -- Decomposes single NetCDF file into many files (converse to mppnccombine).  The output files are created in the current working directory by default and prefixed with the input file name, i.e. in.nc.0000 ...\n\
\n\
Usage: mppncscatter [OPTION...] in.nc\n\
\n\
  -h, --help           Give this usage message.\n\
      --usage          \n\
  -i, --io-layout-x N  Set io-layout for X dimension.\n\
  -j, --io-layout-y N  Set io-layout for Y dimension.\n\
  -n, --dry-run        Run without writing files.\n\
  -p, --prefix PATH    Prefix path for output files.\n\
  -s, --start N        Start file name suffix numbers from N (default is 0).\n\
  -V, --verbose        Progressively output messages to stdout.\n\
  -v, --version        Print program version.\n\
  -w, --width N        Width for output filename suffix digits.\n\
                       -w 4 (default) creates in.nc.0000 ...\n\
  -x, --npx N          Try to split domain evenly into N columns.\n\
  -X, --xdims d1,...   List of X dimension names to scatter\n\
                       (for those not detectable through metadata).\n\
  -y, --npy N          Try to split domain evenly into N rows.\n\
  -Y, --ydims d1,...   List of Y dimension names to scatter\n\
                       (for those not detectable through metadata).\n\
\n\
Report bugs to Remik . Ziemlinski @ noaa . gov.\n\
"
/*-------------------------------------------------------------------*/
typedef struct MNSOPTS {
  char    dryrun;   /* If should do dry run. */
	char*   filein;		/* Input filename (allocated). */
	int     help;     /* give usage insructions */
	int 		nx;				/* Number of columns to split file into. (Required) */
	int 		ny;				/* Number of rows to split file into. (Required) */
  int     nxio;     /* io-layout for x dims. */
  int     nyio;     /* io-layout for y dims. */
  char*   prefix;   /* Output prefix. */
	int     start;    /* Start filename number suffix from this. */
	int     version;  /* give program version */
	int     verbose;  /* Verbose echos to stdout. */
  int     width;    /* Width of output filename digit suffix. */
	char**	xdims;		/* List of xdim names to scatter. */
	int			xdims_len;/* Number of names in above list. */
	char**	ydims;		/* List of xdim names to scatter. */
	int			ydims_len;/* Number of names in above list. */
} mnsopts;

int getmnsopts(int argc, char** argv, mnsopts* popts);
void printusage(void);
void printversion(void);
void initmnsopts(mnsopts* popts);
void freemnsopts(mnsopts* popts);

/*-------------------------------------------------------------------*/
static struct option const long_options[] =
{ 
	{"help",        no_argument,            0, 'h'},
  {"io-layout-x", required_argument,      0, 'i'},
  {"io-layout-y", required_argument,      0, 'j'},
	{"usage",       no_argument,            0, 'h'},
	{"dry-run",     no_argument,            0, 'n'},
	{"start",       required_argument,      0, 's'},
	{"version",     no_argument,            0, 'v'},
	{"verbose",     no_argument,            0, 'V'},
  {"width",       required_argument,      0, 'w'},
	{"npx",         required_argument,      0, 'x'},
	{"xdims",       required_argument,      0, 'X'},
	{"npy",         required_argument,      0, 'y'},
	{"ydims",       required_argument,      0, 'Y'},
  {"prefix",      required_argument,      0, 'p'},
	{0, 0, 0, 0}
};

#endif /* MPPNCSCATTER_OPT_H */

