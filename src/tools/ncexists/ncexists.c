/*
  Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
  This program is distributed under the terms of the GNU General Public
  License. See the file COPYING contained in this directory
*/
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

void handle_error(int status);

int main (int argc, char **argv)
{
  
  int status;
  char *filename = NULL;
  char *attr = NULL;
  char *gattr = NULL;
  int ncid;
  char *var_name = NULL;
  int vr_len;
  int var_id;
  int t_len;
  int global = 0;
  int c;
  char *usage = "Usage: ncexists -f filename [ -g global_attribute || -v variable || -v variable -a attribute ]\n       Returns 1 if variable or attribute is found, 0 if not found.\n";
  nc_type vr_type, t_type;
  
  
  
 if (argc == 1)
  {
  	printf (usage);
	//printf ("No arguments: exiting\n");
	return 0;
  }
  else
  {
	while ((c = getopt (argc, argv, "f:v:g:a:")) != -1)
	{
		switch (c)
            	{
           		case 'f':
             		filename = optarg;
             		break;
           		
			case 'g':
         	        gattr = optarg;
			global = 1;
       		        break;
       		        
			case 'v':
			var_name = optarg;
			break;		
			
			case 'a':
        	        attr= optarg;
        	        break;
			
			
			   
			case '?':
             			printf (usage);
             			return 1;
                        default:
				printf (usage);
             			abort ();
           	}	
	 }
	 
	 if (global == 1)
	 {

		status = nc_open(filename, 0, &ncid);	
	 	if (status != NC_NOERR) handle_error(status);
	 
	 	status = nc_inq_att (ncid, NC_GLOBAL, gattr, &t_type, &t_len);
     	 	if (status == NC_NOERR)
		{	
			printf ("1\n");
		}
		else
		{
			printf ("0\n");
		}
	 }
	 else
	 {
	 	status = nc_open(filename, 0, &ncid);	
	 	if (status != NC_NOERR) handle_error(status);
	 
		status = nc_inq_varid (ncid, var_name, &var_id);
     	 	if (status == NC_NOERR) 
		{
                    if ( attr == NULL )
                    {
                         printf ("1\n");
                    } 
                    else 
                    {
	 	 	status = nc_inq_att (ncid, var_id, attr, &vr_type, &vr_len);
         		if (status == NC_NOERR)
			{
				 printf ("1\n");
			}
			else
			{
				printf ("0\n");
			}
                    }
		}
                else
                {
                    printf ("0\n");
                }
	 }
	 
	 
	 status = nc_close(ncid);       /* close netCDF dataset */
     	 if (status != NC_NOERR) handle_error(status);
	
} 

  
  exit(0);
}

void handle_error(int status) 
{
     if (status != NC_NOERR) 
     {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(-1);
     }
}
