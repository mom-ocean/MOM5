/*
   Copyright (C) 2004-2007, Remik Ziemlinski @ noaa gov

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

#ifndef STRLIST_H
#define STRLIST_H 1

#include "common.h"

/* 0 success, 1 fail; allocates memory for 'size' elements and sets n=0 */
int newstringlist(char*** list, int* n, int size);
/* 0 success, 1 fail; frees memory for each string, but not the list of strings pointer. */
int clearstringlist(char** list, int n);
/* traverses comma-delimited string, returns list and number of list items */
void getstringlist(char *optarg, char*** list, int* nitems);
/* free each string in list. Internally sets *list to NULL. */
void freestringlist(char*** list, int nitems);
/* check membership; returns 0 - false, 1 - true */
int instringlist(char** list, char* str, int nitems);
/* Does deep copy of string. Fails if list is full (no NULL entries). */
int addstringtolist(char** list, char* string, int nitems);
/* Grows the list and sets new number of string pointers. */
int appendstringtolist(char*** list, char* string, int *nitems);
void printstrlist(char** list, int n, FILE* f);

/* find string list union; "union" must already be allocated to 
   accomodate size(set1) + size(set2)
   does not assume listunion is empty, therefore can use as a copy function
*/
int strlistu(char** list1, char** list2, char** listunion, int n1, int n2, int nu);
/* find string list simple difference; 
   list1 items are returned unless exist in list2
   "diff" must already be allocated to 
   accomodate size(set1) + size(set2) worst case.
*/
int strlistsd(char** list1, char** list2, char** listdiff, int n1, int n2, int nsd);

/* Does deep copy. Requires ndst >= nsrc. Strings in listdst after nsrc position are left intact.*/
int copystrlist(char** listsrc, char** listdst, int nsrc, int ndst);
/* Get number of strings in list (and not number of available slots), which don't have to occupy every array position contiguously. */
int getnumstrlist(char** list, int nlist);

#endif
