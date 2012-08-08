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

#include "common.h"

void* 
xmalloc(size_t num)
{
  void* new = (void*)malloc(num);
  if(!new)
    {
      fprintf(stderr, "Memory exhausted\n");
      exit(EXIT_FATAL);
    }

  return new;
}

void*
xrealloc(void* p, size_t num)
{
  void* new;

  if(!p)
    return xmalloc(num);

  new = (void*) realloc(p, num);
  if(!new)
    {
      fprintf(stderr, "Memory exhausted\n");
      exit(EXIT_FATAL);
    }

  return new;
}

void*
xcalloc(size_t num, size_t size)
{
  void* new = (void*) xmalloc(num*size);
  bzero(new, num*size);
  return new;
}

