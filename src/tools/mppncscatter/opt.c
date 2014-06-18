#include "opt.h"

/*-------------------------------------------------------------------*/
void initmnsopts(mnsopts* popts)
{
	if (popts == NULL) 
		return;

	popts->dryrun = 0;
	popts->filein = NULL;
	popts->help = 0;
	popts->nx = 0;
	popts->ny = 0;
	popts->nxio = 0;
	popts->nyio = 0;
  popts->prefix = 0;
	popts->start = 0;
	popts->verbose = 0;
	popts->version = 0;
  popts->width = 4;
	popts->xdims = NULL;
	popts->xdims_len = 0;
	popts->ydims = NULL;
	popts->ydims_len = 0;
}
/*-------------------------------------------------------------------*/
void freemnsopts(mnsopts* popts)
{
	if (popts == NULL) return;

	if (popts->filein != NULL) {
		free(popts->filein);
		popts->filein = NULL;
	}

	if (popts->prefix != NULL) {
		free(popts->prefix);
		popts->prefix = NULL;
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
/*-------------------------------------------------------------------*/
void printusage()
{
	fprintf(stderr, USAGE);
	printversion();
}
/*-------------------------------------------------------------------*/
void printversion()
{
	fprintf(stderr, "mppncscatter %s\n", MPPNCSCATTER_VERSION);
	fprintf(stderr, "Built with NetCDF %s\n", nc_inq_libvers());
	fprintf(stderr, "Copyright (C) 2007,2009-2010,2013 Remik Ziemlinski\n\
\n\
This program comes with NO WARRANTY, to the extent permitted by law.\n\
You may redistribute copies of this program\n\
under the terms of the GNU General Public License.\n"); 
}
/*-------------------------------------------------------------------*/
int getmnsopts(int argc, char** argv, mnsopts* popts)
{
	int c;       
	char *token, *cp;
  const char delimiters[] = ",";
  int len = 0;
  
	if (popts == NULL) return -1;
  
	if (newstringlist(&popts->xdims, &c, NC_MAX_DIMS)) {
		printf("ERROR: Failed to allocate memory for X dimension list.\n");
		exit(1);
	}

	if (newstringlist(&popts->ydims, &c, NC_MAX_DIMS)) {
		printf("ERROR: Failed to allocate memory for Y dimension list.\n");
		exit(1);
	}
	
	while ( (c = getopt_long(argc, argv, "hi:j:np:s:vVw:x:y:X:Y:", long_options, 0))
					!= -1
					)
		switch (c)
			{
			case 'i':
				popts->nxio = atoi(optarg);
				break;

			case 'j':
				popts->nyio = atoi(optarg);
				break;
        
		  case 'n':
		    popts->dryrun = 1;
		    break;
		    
      case 'p':
        len = strlen(optarg);
        popts->prefix = (char*)malloc(len+1);
        strcpy(popts->prefix, optarg);
        
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

			case 'w':
				popts->width = atoi(optarg);
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