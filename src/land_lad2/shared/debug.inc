! -*-f90-*-
! $Id: debug.inc,v 20.0 2013/12/13 23:29:47 fms Exp $

use land_debug_mod, only : dpri

#if defined(__GFORTRAN__)

#define __DEBUG1__(x) call dpri("x",x);write(*,*)
#define __DEBUG2__(x1,x2) call dpri("x1",x1);call dpri("x2",x2);write(*,*) 
#define __DEBUG3__(x1,x2,x3) call dpri("x1",x1);call dpri("x2",x2);call dpri("x3",x3);write(*,*) 
#define __DEBUG4__(x1,x2,x3,x4) call dpri("x1",x1);call dpri("x2",x2);call dpri("x3",x3);call dpri("x4",x4);write(*,*) 
#define __DEBUG5__(x1,x2,x3,x4,x5) call dpri("x1",x1);call dpri("x2",x2);call dpri("x3",x3);call dpri("x4",x4);call dpri("x5",x5);write(*,*) 

#else

#define __DEBUG1__(x) call dpri(#x,x);write(*,*)
#define __DEBUG2__(x1,x2) call dpri(#x1,x1);call dpri(#x2,x2);write(*,*) 
#define __DEBUG3__(x1,x2,x3) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);write(*,*) 
#define __DEBUG4__(x1,x2,x3,x4) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);call dpri(#x4,x4);write(*,*) 
#define __DEBUG5__(x1,x2,x3,x4,x5) call dpri(#x1,x1);call dpri(#x2,x2);call dpri(#x3,x3);call dpri(#x4,x4);call dpri(#x5,x5);write(*,*) 

#endif
