/*
   Copyright (C) 2004, Remik Ziemlinski <first d0t surname att n0aa d0t g0v>

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
*/

#include "strlist.h"

int newstringlist(char*** list, int* n, int size)
{
	int i;
	if (*list != NULL) return EXIT_FAILURE;
  
	*list = XMALLOC(char*, size);
	for(i=0; i<size; ++i) {
		(*list)[i] = NULL;
	}

	*n = size;
	
	return EXIT_SUCCESS;
}
int clearstringlist(char** list, int n)
{        
	if (    (list == NULL) ||
					(n <= 0)
					)
		return 1;
	
	for(; --n >= 0;) {
		if (list[n] == NULL) continue;

		XFREE(list[n]);	
		list[n] = NULL;
	}
  
	return EXIT_SUCCESS;
}
void getstringlist(char *optarg, char*** list, int* nitems)
{
	char *p;
	char **plist;
	int i=0; 
	p = optarg;
	*nitems = 1;
	
	/* find # comma-delimited vars */
	while (*p++) { if (*p == ',') ++(*nitems); }

	/*list = XMALLOC(char*, *nitems);*/
	
	if (*list == NULL) 
	{
					return;
	}

	/*plist = *list;*/

	/* copy var names into list */
	for (p = (char*)strtok(optarg, ",");
				p != NULL;
				p = (char*)strtok((char*) NULL, ",")) 
	{
		(*list)[i] = XMALLOC(char, strlen(p) + 1);
		
		if ((*list)[i] == NULL) { /* Abnormal failure. */
			fprintf(stderr, "ERROR: Failed to allocate memory for string.\n");
			exit(1);
		}

		strcpy((*list)[i], p);
		++i;
	}
}
void freestringlist(char*** list, int nitems)
{
        if (*list == NULL) return;
        
        for(; nitems > 0; --nitems)
        {
					XFREE((*list)[nitems-1]);
        }
        
        XFREE(*list);
				*list = NULL;
}
void printstrlist(char** list, int n, FILE* f)
{
    int i;
    for(i = 0; i < n; ++i) 
        fprintf(f, "%s ", list[i]);
    
    fprintf(f, "\n");
}
int instringlist(char** list, char* str, int nitems)
{       
	int i;
	if ((list == NULL) || (str == NULL)) return 0;

	/* printf("instringlist entered\n"); */

	for(i=0; i < nitems; ++i)
	{
		/* printf("instringlist: testing if list[%d] is NULL.\n", i);
		*/
		if (list[i] == NULL)
			continue;

		/*
		printf("instringlist: list[i] not NULL\n");
		printf("instringlist: comparing \"%s\" with \"%s\".\n", list[i], str);
		*/
		if(strcmp(list[i],str) == 0)
			return 1;        
	}
	return 0;
}
int addstringtolist(char** list, char* string, int nitems)
{
	int i;
	if (string == NULL) return EXIT_FAILED;
	if (list == NULL) return EXIT_FAILED;
	
	/* printf("Adding string to list");*/
	for(i=0; i<nitems; ++i) {
		if (list[i] == NULL) {
			/* Empty slot, so fill it. */
			list[i] = XMALLOC(char, strlen(string) + 1);
			strcpy(list[i], string);
			/* printf(" to slot %d.\n", i); */
			break;
		}
	}

	if (i >= nitems)
		return EXIT_FAILED;
	else
		return EXIT_SUCCESS;
}
/* Returns new pointer with adjusted size. */
int appendstringtolist(char*** list, char* string, int *nitems)
{
	int i, newsize;
	char** newlist = NULL;
	if (string == NULL) return EXIT_FAILED;
	if (*list == NULL) return EXIT_FAILED;

	newstringlist(&newlist, &newsize, *nitems+1);
	copystrlist(*list, newlist, *nitems, newsize);

	newlist[newsize-1] = (char*)malloc(sizeof(char)*(strlen(string)+1));
	strcpy(newlist[newsize-1], string);

	freestringlist(list, *nitems);
	*list = newlist;
	*nitems = newsize;

	return EXIT_SUCCESS;
}
int strlistu(char** list1, char** list2, char** listunion, int n1, int n2, int nu)
{
	int i;

	/* printf("strlistu entered.\n"); */

	if(listunion == NULL) 
	{
					return EXIT_FAILED;
	}
	
	for(i = 0; i < n1; ++i)
	{
		if (list1[i] == NULL) continue;
		/* printf("strlistu: testing if string \"%s\" already in union.\n", list1[i]); */
		if (! instringlist(listunion, list1[i], nu) ) {
			/* printf("strlistu: string \"%s\" not found, trying to add.\n", list1[i]); */
			addstringtolist(listunion, list1[i], nu);        
		}
	}                
	for(i = 0; i < n2; ++i)
	{
		if (list2[i] == NULL) continue;
		if (! instringlist(listunion, list2[i], nu) ) {
			addstringtolist(listunion, list2[i], nu);        
		}
	}

	return EXIT_SUCCESS;
}
int strlistsd(char** list1, char** list2, char** listdiff, int n1, int n2, int nsd)
{
        int i;
       
        if(listdiff == NULL) 
        {
                return EXIT_FAILED;
        }
        
        if(n1 == 0)
        {
                return EXIT_SUCCESS;
        } else if (n2 == 0)
        {
                /* exclude list is empty, return complete list1 */
                for(i = 0; i < n1; ++i)
                        addstringtolist(listdiff, list1[i], nsd);
        } else
        {
                for(i = 0; i < n1; ++i)
                {
										if(! instringlist(list2, list1[i], n2) )
											{
												addstringtolist(listdiff, list1[i], nsd);
											}
                }
        }
        return EXIT_SUCCESS;
}
int copystrlist(char** listsrc, char** listdst, int nsrc, int ndst)
{
	int i;

	if (nsrc > ndst)
		return EXIT_FAILED;

	for(i=0; i < ndst; ++i) {
		XFREE(listdst[i]);	
		listdst[i] = NULL;
	}

	for(i=0; i < nsrc; ++i) {
		if (listsrc[i] == NULL) continue;

		listdst[i] = XMALLOC(char, strlen(listsrc[i]) + 1);
		strcpy(listdst[i], listsrc[i]);
	}

	return EXIT_SUCCESS;
}
int getnumstrlist(char** list, int nlist)
{
	int n = 0;
	int i;

	for(i=0; i < nlist; ++i) {
		if ((list[i] != NULL) && (strlen(list[i])>0))
			++n;
	}
	
	return n;
}
