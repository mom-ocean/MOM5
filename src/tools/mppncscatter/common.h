/*
   Copyright (C) 2010, Remik Ziemlinski <first d0t surname att n0aa d0t g0v>

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

#ifndef COMMON_H
#define COMMON_H 1

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <memory.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <errno.h>
#include <netcdf.h>
#include <math.h>

#ifndef errno
extern int errno;
#endif

#include "xmalloc.h"

#include <limits.h>
#include <inttypes.h>

#define EXIT_SUCCESS        0     
#define EXIT_FAILURE        1
#define EXIT_FATAL          2
#define EXIT_FAILED         3

#define DEBUG 0

#define MPPNCSCATTER_VERSION "0.3.2"

#endif /* !COMMON_H */
