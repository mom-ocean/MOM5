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

#ifndef XMALLOC_H
#define XMALLOC_H 1

#include "common.h"

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS   }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif 

#define XCALLOC(type, num) \
((type*) xcalloc ((size_t)(num), (size_t)sizeof(type)))

#define XMALLOC(type, num) \
((type*) xmalloc((size_t)((num) * sizeof(type))))

#define XREALLOC(type, p, num) \
((type*) xrealloc ((p), (size_t)((num) * sizeof(type))))

#define XFREE(stale) \
do { \
if (stale) { free(stale); stale=0; } \
} while  (0)

BEGIN_C_DECLS

extern void* xcalloc(size_t num, size_t size);
extern void* xmalloc(size_t num);
extern void* xrealloc(void* p, size_t num);

END_C_DECLS

#endif
