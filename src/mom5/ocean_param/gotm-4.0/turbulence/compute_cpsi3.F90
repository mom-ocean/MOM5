#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate c3 from steady-state Richardson number\label{sec:c3}
!
! !INTERFACE:
   REALTYPE function compute_cpsi3(c1,c2,Ri)
!
! !DESCRIPTION:
! Numerically computes $c_{\psi 3}$ for two-equation models from  given
! steady-state Richardson-number $Ri_{st}$ and parameters
! $c_{\psi 1}$ and $c_{\psi 2}$ according to \eq{Ri_st}.
! A Newton-iteration is used to solve the resulting
! implicit non-linear equation.
!
! !USES:
   use turbulence, only:           an,as,cmue1,cmue2
   use turbulence, only:           cm0,cm0_fix,Prandtl0_fix
   use turbulence, only:           turb_method,stab_method
   use turbulence, only:           Constant
   use turbulence, only:           MunkAnderson
   use turbulence, only:           SchumGerz
   use turbulence, only:           EiflerSchrimpf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)            :: c1,c2,Ri
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf
!
! $Log: compute_cpsi3.F90,v $
! Revision 20.0  2013/12/14 00:13:35  fms
! Merged revision 1.1.2.1 onto trunk
!
! Revision 1.1.2.1  2012/05/15 16:00:53  smg
! initial cvs ci for these modules to mom5.
! AUTHOR:Griffies
! REVIEWERS:
! TEST STATUS:
! CHANGES PUBLIC INTERFACES?
! CHANGES ANSWERS?
!
! Revision 1.1.2.1.390.1  2012/04/23 20:30:29  smg
! updated to the gotm-2012.03.09 CVS tag.
! AUTHOR:Martin Schmidt
! REVIEWERS:Griffies
! TEST STATUS:
! CHANGES PUBLIC INTERFACES?
! CHANGES ANSWERS?
!
! Revision 1.3  2010-09-17 12:53:52  jorn
! extensive code clean-up to ensure proper initialization and clean-up of all variables
!
! Revision 1.2  2007-01-06 11:49:15  kbk
! namelist file extension changed .inp --> .nml
!
! Revision 1.1  2005/06/27 10:54:33  kbk
! new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
     integer                       :: i
     integer,parameter             :: imax=100
     REALTYPE                      :: fc,fp,step,ann
     REALTYPE,parameter            :: e=1.e-8
!
!-----------------------------------------------------------------------
!BOC
   ann=5.
   do i=0,imax
      an(1)=ann
      as(1)=an(1)/Ri
      if (turb_method.eq.2) then
         select case(stab_method)
            case(Constant)
               cmue1=cm0_fix
               cmue2=cm0_fix/Prandtl0_fix
            case(MunkAnderson)
               call cmue_ma(2)
            case(SchumGerz)
               call cmue_sg(2)
            case(EiflerSchrimpf)
               call cmue_rf(2)
         end select
      else
         call cmue_d(2)
      end if
      fc=cmue1(1)*an(1)/Ri-cmue2(1)*an(1)-cm0**(-3)
      an(1)=ann+e
      as(1)=an(1)/Ri
      if (turb_method.eq.2) then
         select case(stab_method)
            case(Constant)
               cmue1=cm0_fix
               cmue2=cm0_fix/Prandtl0_fix
            case(MunkAnderson)
               call cmue_ma(2)
            case(SchumGerz)
               call cmue_sg(2)
            case(EiflerSchrimpf)
               call cmue_rf(2)
         end select
      else
         call cmue_d(2)
      end if
      fp=cmue1(1)*an(1)/Ri-cmue2(1)*an(1)-cm0**(-3)
      step=-fc/((fp-fc)/e)
      ann=ann+0.5*step
      if (abs(step).gt.100.) then
         STDERR 'Method for calculating c3 does not converge.'
         STDERR 'Probably, the prescribed steady-state Richardson number'
         STDERR 'is outside the range of the chosen stability function.'
         STDERR 'Please change gotmturb.nml accordingly.'
         STDERR 'If the problem persists, then use another'
         STDERR 'stability function or Algebraic Stress Model.'
         STDERR 'Program aborts now in turbulence.F90.'
         stop
      endif
      if (abs(step).lt.1.e-10) goto 111
   end do
111   an(1)=ann
   as(1)=an(1)/Ri
   if (turb_method.eq.2) then
      select case(stab_method)
         case(Constant)
            cmue1=cm0_fix
            cmue2=cm0_fix/Prandtl0_fix
         case(MunkAnderson)
            call cmue_ma(2)
         case(SchumGerz)
            call cmue_sg(2)
         case(EiflerSchrimpf)
            call cmue_rf(2)
      end select
   else
      call cmue_d(2)
   end if

   compute_cpsi3 = c2+(c1-c2)/Ri*cmue1(1)/cmue2(1)

   return
   end function compute_cpsi3

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
