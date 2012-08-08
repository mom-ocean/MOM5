module topo_drag_mod

!==========================================================================
! TOPOGRAPHIC DRAG CLOSURE FOR GENERAL CIRCULATION MODELS -- Garner (2005)
!==========================================================================

!--------------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!--------------------------------------------------------------------------

  use       Fms_Mod, only: ERROR_MESG, FATAL, &
                           mpp_pe, mpp_root_pe, &
                           write_version_number, stdlog

  implicit none

  private

  character(len=128) :: version = '$Id: topo_drag.F90,v 17.0 2009/07/21 02:58:21 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201207 $'
  logical            :: module_is_initialized = .false.

  public topo_drag, topo_drag_init, topo_drag_end
  public topo_drag_restart

contains

!#############################################################################      
subroutine topo_drag                                                   &
                                     ( is, js, delt, uwnd, vwnd, atmp, &
                                           pfull, phalf, zfull, zhalf, &
                          dtaux, dtauy, dtemp, taux, tauy, taus, kbot )

integer, intent(in) :: is, js
real,    intent(in) :: delt
integer, intent(in), optional, dimension(:,:) :: kbot

! INPUT
! -----

! UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
! VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
! ATMP     Temperature at full levels (IDIM x JDIM x KDIM)
! PFULL    Pressure at full levels (IDIM x JDIM x KDIM)
! PHALF    Pressure at half levels (IDIM x JDIM x KDIM+1)
! ZFULL    Height at full levels (IDIM x JDIM x KDIM)
! ZHALF    Height at half levels (IDIM x JDIM x KDIM+1)

real, intent(in), dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:,:) :: pfull, phalf, zfull, zhalf


! OUTPUT
! ------


! DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
! DTEMP        Tendency of the temperature in K/s (IDIM x JDIM x KDIM)
! TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
! TAUS         clipped saturation momentum flux (IDIM x JDIM x KDIM) for diagnostics

real, intent(out), dimension(:,:)   :: taux, tauy
real, intent(out), dimension(:,:,:) :: dtaux, dtauy, dtemp, taus



!---------------------------------------------------------------------

      call error_mesg('topo_drag', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag

  !=====================================================================

  subroutine topo_drag_init(lonb,latb)

    real,    intent(in), dimension(:,:) :: lonb,latb


!------- write version number and namelist ---------

    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
    endif

    module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('topo_drag_init', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag_init

  !=====================================================================

  subroutine topo_drag_end

      module_is_initialized = .false.
 
!---------------------------------------------------------------------

      call error_mesg('topo_drag_end', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine topo_drag_end

!#############################################################################      

!#######################################################################
! <SUBROUTINE NAME="topo_drag_restart">
!
! <DESCRIPTION>
!   dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine topo_drag_restart(timestamp)
   character(len=*), intent(in), optional :: timestamp

end subroutine topo_drag_restart
! </SUBROUTINE>

end module topo_drag_mod
