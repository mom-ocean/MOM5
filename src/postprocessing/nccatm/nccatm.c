/*
  nccatm

  This is an extensive modification of the original "nccat" program which
  allows more than two netCDF files to be concatenated at once.

  V1.2: Added an option to limit the number of dimensions per variable which
        are processed at once.  The output file is now created from the first
        input file if it doesn't exist.
  V1.1: Changed loop order for increased I/O efficiency; records are now the
        innermost loop then the variables loop.  Added progress information
        options.
  V1.0: Original release.

  Modified by Hans Vahlenkamp
  Geophysical Fluid Dynamics Laboratory/NOAA
  Princeton University Forrestal Campus
  Last Updated: 3/22/00
*/

/*
 *	nccat:	Concatenate two netCDF files.
 *
 *	The NetCDF files must be identical in variable
 *	names and dimensions.  Each variable must have
 *	a leftmost NC_UNLIMITED record dimension.
 *
 *	Copyright (C) 1991 Charles R. Denham, Zydeco.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <netcdf.h>
#include <udunits.h>

/* Auxiliary function prototypes */
void usage();
int clean_up(int, int);
int close_file(int);
char *ltoac(unsigned long);


int main(int argc, char *argv[])
  {
   unsigned char netcdfffiospec=0;  /* Use the shell's NETCDF_FFIOSPEC? */
   unsigned char verbose=0;  /* Verbose mode; print extensive progress info */
   unsigned char quiet=0;  /* Quiet mode; no progress output */
   unsigned char checkrecmem=0;  /* Show memory usage for a record and quit */
   int maxdims=(-1);  /* Max. # of variable dimensions to process at once */
   int outfid=(-1), infid;  /* Argument #s of output and first input files */
   struct stat outstatbuf;  /* Information about output file */
   char syscmd[2048];  /* Shell copy command */
   struct utsname minfo;  /* Machine information */
   long a, d, v, v2, l, f, r, rn=0;  /* Loop variables */
   int outcdfid=(-1), incdfid=(-1);  /* IDs of the I/O netCDF files */
   int outndims, inndims;  /* Number of dimensions */
   int outnvars, innvars;  /* Number of variables */
   int natts;  /* Number of attributes (ignored) */
   int outrecdim, inrecdim;  /* ID of the record dimensions */
   char dimname[MAX_NC_NAME];                             /* Names of record */
   char inrecname[MAX_NC_NAME], outrecname[MAX_NC_NAME];  /* dimensions      */
   long outdimsize[MAX_NC_DIMS];  /* Sizes of the dimensions */
   long indimsize[MAX_NC_DIMS];   /*            "            */
   char outvarname[MAX_NC_VARS][MAX_NC_NAME];  /* Names of the variables */
   char invarname[MAX_NC_VARS][MAX_NC_NAME];   /*           "            */
   nc_type outdatatype[MAX_NC_VARS];  /* Data types of the variables */
   nc_type indatatype[MAX_NC_VARS];   /*              "              */
   nc_type outrecdimtype=(-1);  /* Output record dimension's data type */
   int outvarndims[MAX_NC_VARS];  /* Number of dimensions for each variable */
   int invarndims[MAX_NC_VARS];   /*                   "                    */
   int outvardim[MAX_NC_VARS][MAX_NC_DIMS];  /* Dimensions for each variable */
   int invardim[MAX_NC_VARS][MAX_NC_DIMS];   /*              "               */
   int outrecvar=(-1), inrecvar=(-1);  /* IDs of the record dimensions */
   unsigned long recmemsize, maxmemsize=0;  /* Memory size of var. record */
   char outunitsstr[1024], inunitsstr[1024];  /* "units" attribute value */
   utUnit outunit, inunit;                    /*            "            */
   unsigned char useoutunits=1, useinunits;  /* Are the "units" attributes  */
                                             /* of the record dimensions    */
                                             /* valid according to udunits? */
   double inslope, inint;  /* Conversion values for input to output units */
   long totalnrecords;  /* Number of records processed */
   long outrecstart;  /* Starting record value for output */
   double *outrecvals=NULL;  /* Output record values */
   long nrecords;  /* Number of records in an input file */
   long incoord[MAX_NC_DIMS], outcoord[MAX_NC_DIMS];  /* Ranges for reading  */
   long incount[MAX_NC_DIMS], outcount[MAX_NC_DIMS];  /* and writing records */
   int status;  /* Flag for indicating record dimension errors */
   long datasize;  /* Number of data array elements */
   int outvarid;  /* ID of input variable in output file */
   int ndimloops;  /* Number of dimensions to loop over per variable record */
   long outerloopcount;  /* Number of iterations over non-record dimensions */
   unsigned char update;  /* Flag for updating the hyperslab coordinates */
   unsigned char *values;  /* I/O transfer buffer */
   double *valueptr;  /* Pointer for adjustment to record variable values */

   /* Check the command-line arguments */
   if (argc < 3)
     {
      usage(); return(-1);
     }
   for (a=1; a < argc; a++)
     {
      if (!strcmp(argv[a],"-n")) netcdfffiospec=1;
      else if (!strcmp(argv[a],"-v")) verbose=1;
      else if (!strcmp(argv[a],"-q")) quiet=1;
      else if (!strcmp(argv[a],"-c")) checkrecmem=1;
      else if (!strcmp(argv[a],"-d"))
        {
         if (++a > argc-1)
           {
            usage(); return(-1);
           }
         maxdims=atoi(argv[a]);
         if (maxdims < 1)
           {
            usage(); return(-1);
           }
        }
      else 
        {
         outfid=a; break;
        }
     }
   infid=outfid+1;
   if (outfid==(-1) || (infid > argc-1))
     {
      usage(); return(-1);
     }
   if (!strcmp(argv[outfid],argv[infid]))
     {
      fprintf(stderr,"Cannot concatenate a file to itself!\n");
      return(-1);
     }

   /* Initialize the udunits library with the definitions file */
   if (utInit("/usr/local/etc/udunits.dat") != 0) return(-1);

   /* Copy the first input file to the output if the output doesn't exist */
   if (stat(argv[outfid],&outstatbuf)!=0)
     {
      if (!quiet) printf("Creating output file from: %s\n",argv[infid]);
      sprintf(syscmd,"cp %s %s",argv[infid],argv[outfid]);
      if (system(syscmd)) return(-1);
      if (++infid > argc-1) return(0);
     }

   /* If nccatm is running on a Cray system then set the "NETCDF_FFIOSPEC" */
   /* environment variable for greater read/write efficiency               */
   uname(&minfo);
   if (strstr(minfo.sysname,"CRAY")!=NULL ||
       strstr(minfo.machine,"CRAY")!=NULL)
     {
      if (netcdfffiospec && getenv("NETCDF_FFIOSPEC")!=NULL)
        {
         if (!quiet)
           printf("Using the current value of the \"NETCDF_FFIOSPEC\" environment variable.\n");
        }
      else
        {
         if (!quiet)
           printf("Setting the \"NETCDF_FFIOSPEC\" environment variable to \"cachea:896:2\".\n");
         if (putenv("NETCDF_FFIOSPEC=cachea:896:2"))
           fprintf(stderr,"Error setting the \"NETCDF_FFIOSPEC\" environment variable!\n");
        }
     }

   /* Disable fatal returns from netCDF library functions */
   ncopts=0;

   /* Open the output file */
   if ((outcdfid=ncopen(argv[outfid],NC_WRITE))==(-1))
     {
      fprintf(stderr,"\"ncopen\" failure on the output file!\n"); utTerm();
      return(-1);
     }

   /* Get some information about the output file */
   if (ncinquire(outcdfid,&outndims,&outnvars,&natts,&outrecdim)==(-1))
     {
      fprintf(stderr,"\"ncinquire\" failure on the output file!\n");
      clean_up(outcdfid,-1); return(-1);
     }
   if (outrecdim==(-1))
     {
      fprintf(stderr,"No record dimension exists in the output file!\n");
      clean_up(outcdfid,-1); return(-1);
     }
   if (outnvars > MAX_NC_VARS)
     {
      fprintf(stderr,"Too many variables in the output file, limit = %d!\n",
              MAX_NC_VARS);
      clean_up(outcdfid,-1); return(-1);
     }
   if (outndims > MAX_NC_DIMS)
     {
      fprintf(stderr,"Too many dimensions in the output file, limit = %d!\n",
              MAX_NC_DIMS);
      clean_up(outcdfid,-1); return(-1);
     }
   for (d=0; d < outndims; d++)
     {
      if ((ncdiminq(outcdfid,d,dimname,&outdimsize[d]))==(-1))
        {
         fprintf(stderr,"\"ncdiminq\" failure on the output file!\n");
         clean_up(outcdfid,-1); return(-1);
        }
      if (d==outrecdim) strcpy(outrecname,dimname);
     }
   for (v=0; v < outnvars; v++)
     {
      if ((ncvarinq(outcdfid,v,outvarname[v],&outdatatype[v],&outvarndims[v],
                    outvardim[v],&natts))==(-1))
        {
         fprintf(stderr,"\"ncvarinq\" failure on the output file!\n");
         clean_up(outcdfid,-1); return(-1);
        }
      if (!strcmp(outvarname[v],outrecname))
        {
         outrecdimtype=outdatatype[v]; outrecvar=v;
        }
     }
   if (outrecvar==(-1))
     {
      fprintf(stderr,"Record dimension is not defined as a variable!\n");
      clean_up(outcdfid,-1); return(-1);
     }

   /* Show the memory required to process a variable record */
   if (checkrecmem)
     {
      for (v=0; v < outnvars; v++)
       {
        recmemsize=1;
        for (d=0; d < outvarndims[v]; d++)
          {
           if (outvardim[v][d]!=outrecdim)
             recmemsize*=outdimsize[outvardim[v][d]];
          }
        recmemsize*=nctypelen(outdatatype[v]);
        if (recmemsize > maxmemsize) maxmemsize=recmemsize;
       }
      printf("At least %s bytes will be required to process the records.\n",
             ltoac(maxmemsize));
      return(clean_up(outcdfid,-1));
     }

   /* Record dimension must have a datatype of double */
   if (outrecdimtype!=NC_DOUBLE)
     {
      fprintf(stderr,"Output file's record dimension must be \"double\" (64-bit floating point) type!\n");
      clean_up(outcdfid,-1); return(-1);
     }

   /* Look for the record dimension's "units" attribute */
   if (ncattget(outcdfid,outrecvar,"units",(void *)outunitsstr)==(-1))
     useoutunits=0;
   else
     {
      if (utScan(outunitsstr,&outunit)!=0) useoutunits=0;
     }
   if (!useoutunits)
     printf("Output file has no valid \"units\" attribute for the record dimension!\n");

   /* Get the size and values of the record coordinates */
   totalnrecords=outdimsize[outrecdim]; outrecstart=0;
   if (useoutunits)
     {
      if ((outrecvals=(double *)malloc(sizeof(double)*totalnrecords))==NULL)
        {
         fprintf(stderr,"Cannot allocate %s bytes of memory!\n",
                 ltoac(sizeof(double)*totalnrecords));
         clean_up(outcdfid,-1); return(-1);
        }
      if ((ncvarget(outcdfid,outrecvar,&outrecstart,&totalnrecords,
                    (void *)outrecvals))==(-1))
        {
         fprintf(stderr,"\"ncvarget\" failure on the output file!\n");
         clean_up(outcdfid,-1); return(-1);
        }
     }

   /* Should new variable space be prefilled (probably not for speed)? */
#ifdef NO_NCPREFILL
   if (!quiet) printf("Setting NC_NOFILL mode...\n");
   ncsetfill(outcdfid,NC_NOFILL);
#endif

   /* Loop over all the input files */
   for (f=infid; f < argc; f++)
     {
      if (!quiet) printf("Appending file: %s\n",argv[f]);

      /* Open an input file */
      if ((incdfid=ncopen(argv[f],NC_NOWRITE))==(-1))
        {
         fprintf(stderr,"\"ncopen\" failure on the input file!\n");
         clean_up(outcdfid,-1);
         if (outrecvals!=NULL) free(outrecvals);
         return(-1);
        }

      /* Get some information about the input file */
      if (ncinquire(incdfid,&inndims,&innvars,&natts,&inrecdim)==(-1))
        {
         fprintf(stderr,"\"ncinquire\" failure on the input file!\n");
         clean_up(outcdfid,incdfid);
         if (outrecvals!=NULL) free(outrecvals);
         return(-1);
        }	 
      if (inrecdim==(-1))
        {
         fprintf(stderr,"No record dimension exists in the input file!\n");
         clean_up(outcdfid,incdfid);
         if (outrecvals!=NULL) free(outrecvals);
         return(-1);
        }
      for (d=0; d < inndims; d++)
        {
         if ((ncdiminq(incdfid,d,dimname,&indimsize[d]))==(-1))
           {
            fprintf(stderr,"\"ncdiminq\" failure on the input file!\n");
            clean_up(outcdfid,incdfid); return(-1);
           }
         if (d==inrecdim) strcpy(inrecname,dimname);
        }
      for (v=0; v < innvars; v++)
        {
         if ((ncvarinq(incdfid,v,invarname[v],&indatatype[v],&invarndims[v],
                       invardim[v],&natts))==(-1))
           {
            fprintf(stderr,"\"ncvarinq\" failure on the input file!\n");
            clean_up(outcdfid,incdfid); return(-1);
           }
         if (!strcmp(invarname[v],inrecname)) inrecvar=v;
        }
      if (inrecvar==(-1))
        {
         fprintf(stderr,"Record dimension is not defined as a variable!\n");
         clean_up(outcdfid,-1); return(-1);
        }

      /* Look for the record dimension's "units" attribute */
      if (useoutunits)
        {
         useinunits=1;
         if (ncattget(incdfid,inrecvar,"units",(void *)inunitsstr)==(-1))
           useinunits=0;
         else
           {
            if (utScan(inunitsstr,&inunit)!=0) useinunits=0;
           }
         if (useinunits)
           {
            if (utConvert(&inunit,&outunit,&inslope,&inint)!=0) useinunits=0;
           }
        }
      else useinunits=0;
      if (!useinunits && useoutunits)
        printf("Input file has no valid \"units\" attribute for the record dimension!\n");

      /* Update the number of records */
      nrecords=indimsize[inrecdim]; totalnrecords+=nrecords;
      if (!quiet)
        {
         printf("Records to transfer:  %ld\n",nrecords);
         printf("Records in result:    %ld\n",totalnrecords);
         printf("Variables to process: %d\n",innvars);
        }

      /* Process all the records in each file */
      if (!quiet && !verbose) printf("Record =");
      for (r=0; r < nrecords; r++)
        {
         if (!quiet)
           {
            if (verbose) printf("  Record = %ld\n",r+1);
            else
              {
               if (rn < 10)
                 {
                  printf(" %ld",r+1); rn++;
                 }
               else
                 {
                  rn=1; printf("\n        %ld",r+1);
                 }
               if (r < nrecords-1) printf(",");
              }
            fflush(stdout);
           }

         /* Process all the variables at each record */
         for (v=0; v < innvars; v++)
           {
            /* Skip if the variable is a non-record dimension */
            if (ncdimid(incdfid,invarname[v])!=(-1))
              {
               if (v!=inrecvar)
                 {
                  if (!quiet && verbose)
                    printf("    Variable = %s; is a non-record dimension, skipping\n",
                           invarname[v]);
                  continue;
                 }
              }

            /* Does the input variable exist in the output file? */
            status=0;
            for (v2=0; v2 < outnvars; v2++)
              {
               if (!strcmp(outvarname[v2],invarname[v]))
                 {
                  outvarid=v2; status=1; break;
                 }
              }
            if (!status)
              {
               if (!quiet && verbose)
                 printf("    Variable = %s; not in output file\n",
                        invarname[v]);
               continue;
              }

            /* Is it a record variable? */
            if (invardim[v][0]!=inrecdim)
              {
               if (!quiet && verbose)
                 printf("    Variable = %s; no record dimension\n",
                        invarname[v]);
               continue;
              }
            if (outvardim[outvarid][0]!=outrecdim)
              {
               if (!quiet && verbose)
                 printf("    Variable = %s; no record dimension in output file\n",
                        invarname[v]);
               continue;
              }

            /* Check some input variable info against the output */
            if (outdatatype[outvarid]!=indatatype[v])
              {
               fprintf(stderr,"\nVariable = %s; different input and output data types\n",
                       invarname[v]);
               fprintf(stderr,"Aborting...output file may have corrupted data!\n");
               clean_up(outcdfid,incdfid);
               if (outrecvals!=NULL) free(outrecvals);
               return(-1);
              }
            if (outvarndims[outvarid]!=invarndims[v])
              {
               fprintf(stderr,"\nVariable = %s; different number of input and output dimensions\n",
                       invarname[v]);
               fprintf(stderr,"Aborting...output file may have corrupted data!\n");
               clean_up(outcdfid,incdfid);
               if (outrecvals!=NULL) free(outrecvals);
               return(-1);
              }   
            for (d=0; d < invarndims[v]; d++)
              {
               if (invardim[v][d]!=inrecdim &&
                   indimsize[invardim[v][d]]!=outdimsize[outvardim[outvarid][d]])
                 {
                  fprintf(stderr,"\nVariable = %s; different input and output dimension sizes\n",
                          invarname[v]);
                  fprintf(stderr,"Aborting...output file may have corrupted data!\n");
                  clean_up(outcdfid,incdfid);
                  if (outrecvals!=NULL) free(outrecvals);
                  return(-1);
                 }
              }

            if (!quiet && verbose)
              {
               printf("    Variable = %s; %ld to go ...\n",invarname[v],
                      (innvars-v-1));
               fflush(stdout);
              }

            /* Determine the size of a read/write data block */
            if (maxdims==(-1) || maxdims >= invarndims[v]) ndimloops=1;
            else ndimloops=invarndims[v]-maxdims;
            datasize=1;
            for (d=0; d < invarndims[v]; d++)
              {
               incoord[d]=0; outcoord[d]=0;
               if (d < ndimloops)
                 {
                  incount[d]=1; outcount[d]=1;
                 }
               else
                 {
                  incount[d]=indimsize[invardim[v][d]];
                  outcount[d]=outdimsize[outvardim[outvarid][d]];
                 }
               datasize*=incount[d];
              }

            /* Allocate a transfer buffer */
            if ((values=(unsigned char *)malloc(nctypelen(indatatype[v])*datasize))==NULL)
              {
               fprintf(stderr,"\nCannot allocate %s bytes of memory!\n",
                       ltoac(nctypelen(indatatype[v])*datasize));
               fprintf(stderr,"Aborting...output file may have corrupted data!  Try reducing the # of\n");
               fprintf(stderr,"dimensions being processed per variable record with the \"-d\" option.\n");
               clean_up(outcdfid,incdfid);
               if (outrecvals!=NULL) free(outrecvals);
               return(-1);
              }

            /* Perform the read/write iterations */
            incoord[0]=r; outcoord[0]=totalnrecords-nrecords+r;
            outerloopcount=1; update=1;
            for (l=1; l <= ndimloops-1; l++)
              outerloopcount*=indimsize[invardim[v][l]];
            for (l=0; l < outerloopcount; l++)
              {
               if ((ncvarget(incdfid,v,incoord,incount,(void *)values))==(-1))
                 {
                  fprintf(stderr,"\n\"ncvarget\" failure on the input file!\n");
                  clean_up(outcdfid,incdfid); return(-1);
                 }
               if (useinunits && v==inrecvar)
                 {
                  valueptr=(double *)values;
                  *valueptr=(*valueptr)*inslope+inint;
                 }
               if ((ncvarput(outcdfid,outvarid,outcoord,outcount,(void *)values))==(-1))
                 {
                  fprintf(stderr,"\n\"ncvarput\" failure on the input file!\n");
                  clean_up(outcdfid,incdfid); return(-1);
                 }
               /* Update the hyperslab coordinates */
               for (d=ndimloops-1; d >=1; d--)
                 {
                  if (update)
                    {
                     incoord[d]++; outcoord[d]++;
                    }
                  if (incoord[d] > indimsize[invardim[v][d]]-1)
                    {
                     incoord[d]=0; outcoord[d]=0; update=1;
                    }
                  else update=0;
                 }
               update=1;
              }

            /* Deallocate the transfer buffer */
            free(values);
           }
        }
      if (!quiet && !verbose) printf("\n");
      close_file(incdfid);
     }

   /* Clean up */
   if (outrecvals!=NULL) free(outrecvals);
   return(clean_up(outcdfid,-1));
  }


/* Usage message for the program */
void usage()
  {
   printf("Usage:\n");
   printf("nccatm [-n] [-v | -q] [-c] [-d #] infile1 infile2 ...\n");
   printf("   -n: Use the shell's setting for \"NETCDF_FFIOSPEC\", otherwise it will use\n");
   printf("       the built-in default of \"cachea:896:2\" (Cray only).\n");
   printf("   -v: Verbose mode; print extensive progress information.\n");
   printf("   -q: Quiet mode; no progress output.\n");
   printf("   -c: Show the memory required to process a variable record then quit.\n");
   printf("   -d #: Specify the maximum number of dimensions to process in each read/write\n");
   printf("         operation (>=1); otherwise all the dimensions of a variable's single\n");
   printf("         record are assumed.  This is useful if a single record of a variable\n");
   printf("         is too large for memory allocation.  Example: If a variable has 3\n");
   printf("         dimensions (X,Y,Z) in addition to the record dimension then specifying\n");
   printf("         \"-d 2\" will process the X,Y chunk and loop over Z.\n");
   printf("   infile1: NetCDF filename for input and output.\n");
   printf("   infile2 ...: NetCDF filenames for input.\n");
   printf("Purpose:\n");
   printf("   Concatenate NetCDF record variables.  The NetCDF files must be identical in\n");
   printf("   variable names, dimensions, and data types.  Each variable must have a\n");
   printf("   leftmost NC_UNLIMITED record dimension.  The record dimension must also\n");
   printf("   be type double (64-bit floating point).  The udunits library is used in\n");
   printf("   order to try and resolve time coordinates that have different units between\n");
   printf("   files.  If the output file doesn't exist then it will be a copy of the first\n");
   printf("   input file.\n");
   printf("Example:\n");
   printf("   nccatm foo.nc bar.nc\n");
  }


/* Close the udunits library and any open files */
int clean_up(int cdfid1, int cdfid2)
  {
   int status1=0, status2=0;

   utTerm();
   if (cdfid1 != -1) status1=close_file(cdfid1);
   if (cdfid2 != -1) status2=close_file(cdfid2);
   if (status1 == -1 || status2 == -1) return(-1);
   else return(0);
  }


/* Close an open netCDF file */
int close_file(int cdfid)
  {
   int status;

   if (cdfid != -1)
     {
      if ((status=ncclose(cdfid)) == -1)
        fprintf(stderr,"\"ncclose\" failure!\n");
     }
   return(status);
  }


/* Convert an integer into a string with commas */
char *ltoac(unsigned long value)
  {
   char invaluestr[30], *outvaluestr;
   int invaluelen, ncommas, i, j=0, c=0;

   sprintf(invaluestr,"%ld",value); invaluelen=strlen(invaluestr);
   outvaluestr=(char *)calloc(30,1);
   ncommas=(invaluelen-1)/3; j=invaluelen+ncommas-1;
   for (i=invaluelen-1; i>=0; i--)
     {
      outvaluestr[j--]=invaluestr[i]; c++;
      if (c==3 && i!=0)
        {
         outvaluestr[j--]=','; c=0;
        }
     }
   return(outvaluestr);
  }
