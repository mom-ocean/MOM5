! -*-f90-*-
! $Id: debug.inc,v 17.0 2009/07/21 03:02:31 fms Exp $

#define __DEBUG1__(x) write(*,'(a12,99g)')#x,x
#define __DEBUG2__(x1,x2) write(*,'(99(a12,g))')#x1,x1,#x2,x2 
#define __DEBUG3__(x1,x2,x3) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3 
#define __DEBUG4__(x1,x2,x3,x4) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3,#x4,x4 
#define __DEBUG5__(x1,x2,x3,x4,x5) write(*,'(99(a12,g))')#x1,x1,#x2,x2,#x3,x3,#x4,x4,#x5,x5 
