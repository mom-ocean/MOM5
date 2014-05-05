
#if defined(__GFORTRAN__)

! gfortran does not support ## 
#define CONCAT1(op) op
#define CONCAT2(op,T) CONCAT1(op)T
#define CONCAT3(op,T,D) CONCAT2(op,T)D

#else

#define CONCAT3(op,T,D) op##T##D

#endif
