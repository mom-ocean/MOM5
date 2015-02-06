#include "opt.h"
#include "mppncscatter.h"

/*-------------------------------------------------------------------*/
int main(int argc, char** argv)
{
	int status;
	mnsopts opts;

	status = 0;

	initmnsopts(&opts);
	
	/* parse command-line args.  */
	status = getmnsopts(argc, argv, &opts);
	
	status = (status ? status : mppncscatter(&opts));               
	
	freemnsopts(&opts);
	
 	if (opts.verbose) {
		fprintf(stdout, "Info: Done.\n");
		fflush(stdout);
	}

	return status;
} 