/***********************************************************************
                       creater_vgrid.h
    This header file contains interface to create vertical grid on supergrid. 
    refinement =2 is assumed in this routine.
    contact: Zhi.Liang
************************************************************************/

#ifndef CREATE_VGRID_H_
#define CREATE_VGRID_H_
void create_vgrid(int nbnds, double *bnds, int *nz, double *dbnds, double stretch, 
                  int use_legacy, double *zeta, int *np, const char *center);
#endif
