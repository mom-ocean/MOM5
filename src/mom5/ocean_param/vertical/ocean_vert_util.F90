module ocean_vert_util_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies 
! </CONTACT>
!
!<OVERVIEW>
! This module contains routines for use in vertical mixing. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Routines for vertical mixing schemes. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_tracer_util_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!
!  <DATA NAME="smooth_n2" TYPE="logical">
!  For vertical smoothing the N2 calculation for Richardson number.
!  Default smooth_n2 = .true.
!  </DATA> 
!  <DATA NAME="num_n2_smooth" TYPE="integer">
!  For vertical smoothing N2 for Richardson number.
!  Default num_n2_smooth = 1.
!  </DATA> 
!
!  <DATA NAME="smooth_ri_number" TYPE="logical">
!  For vertical smoothing richardson number.
!  Default smooth_ri_number = .true.
!  </DATA> 
!  <DATA NAME="num_ri_smooth" TYPE="integer">
!  For vertical smoothing richardson number.
!  Default num_ri_smooth = 1.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,        only: epsln
use diag_manager_mod,     only: register_diag_field
use fms_mod,              only: open_namelist_file, check_nml_error, close_file
use mpp_mod,              only: input_nml_file, stdout, stdlog, FATAL
use mpp_mod,              only: mpp_error, mpp_pe, mpp_min, mpp_max
  
use ocean_density_mod,    only: density_delta_z
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value, rho0r, rho0, grav, onehalf, onefourth
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,      only: ocean_density_type, ocean_velocity_type
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk5
use ocean_workspace_mod,  only: wrk1_v, wrk2_v
use ocean_util_mod,       only: diagnose_3d, diagnose_3d_u

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6

#include <ocean_memory.h>
logical :: module_is_initialized  = .false.

! diagnostics 
logical :: used 
integer :: id_ri_num_dudz  =-1
integer :: id_ri_num_n2    =-1
integer :: id_ri_num_cgrid =-1
integer :: id_rit          =-1
integer :: id_riu          =-1

! nml parameters 

! for cgrid Richardson number calculation 
logical :: smooth_n2        = .true.
logical :: smooth_ri_number = .true.
integer :: num_n2_smooth    = 1 
integer :: num_ri_smooth    = 1 


! for debugging 
logical :: debug_this_module = .false. 


type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_vert_util_init
public ri_for_cgrid
public ri_for_bgrid 

namelist /ocean_vert_util_nml/ debug_this_module, &
         smooth_n2, smooth_ri_number, num_n2_smooth, num_ri_smooth


contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_util_init">
!
! <DESCRIPTION>
! Initialize vertical mixing utilities.
! </DESCRIPTION>
!
subroutine ocean_vert_util_init (Grid, Domain, Time)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type),   intent(in)         :: Time
  
  integer :: ioun, io_status, ierr
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


  if (module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_vert_util_mod (ocean_vert_util_init): module already initialized')
  endif 

  module_is_initialized = .true.
  stdlogunit=stdlog()
  write( stdlogunit,'(/a/)') trim(version)

  Dom => Domain
  Grd => Grid

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_vert_util_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_util_nml')
#else
  ioun = open_namelist_file()
  read(ioun, ocean_vert_util_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_vert_util_nml')
  call close_file(ioun)
#endif
  write (stdlogunit, ocean_vert_util_nml)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_vert_util_nml)


  id_ri_num_dudz = register_diag_field ('ocean_model', 'ri_num_dudz',                & 
     Grid%tracer_axes_wt(1:3), Time%model_time,                                      &
     '(du/dz)^2 for Richardson number computed on T-cell bottom using Cgrid routine',&
     '1/s^2', missing_value=missing_value, range=(/-1e1,1e15/))

  id_ri_num_n2 = register_diag_field ('ocean_model', 'ri_num_n2',              & 
     Grid%tracer_axes_wt(1:3), Time%model_time,                                &
     'N^2 for Richardson number computed on T-cell bottom using Cgrid routine',&
     '1/s^2', missing_value=missing_value, range=(/-1e15,1e15/))

  id_rit = register_diag_field ('ocean_model', 'rit',& 
     Grid%tracer_axes_wt(1:3), Time%model_time,      &
     'Richardson number computed on T-cell bottom',  &
     'dimensionless', missing_value=missing_value, range=(/-1e15,1e15/))

  id_riu = register_diag_field ('ocean_model', 'riu',& 
     Grid%tracer_axes_wt(1:3), Time%model_time,      &
     'Richardson number computed on U-cell bottom',  &
     'dimensionless', missing_value=missing_value, range=(/-1e15,1e15/))


end subroutine ocean_vert_util_init
! </SUBROUTINE> NAME="ocean_vert_util_init">


!#######################################################################
! <SUBROUTINE NAME="ri_for_bgrid">
!
! <DESCRIPTION>
!
! Compute Richardson number assuming horizontal B-grid layout.   
!
! Richardson number rit is centered at T-cell.  
! Richardson number riu is centered at U-cell.  
!
! This calculation differs from that in ocean_vert_kpp_mom4p1 
! since here we compute N^2 using locally referenced potential density, 
! as done for tide mixing scheme and as done in ri_for_cgrid.  Other
! features of the calculation, such as the horizontal averaging, 
! agree with ocean_vert_kpp_mom4p1. 
!
! </DESCRIPTION>
!
subroutine ri_for_bgrid (Time, dzwt, drhodT, drhodS, theta, salinity, velu, velv, rit, riu)

  type(ocean_time_type),        intent(in)    :: Time
  real, dimension(isd:,jsd:,:), intent(in)    :: dzwt
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity 
  real, dimension(isd:,jsd:,:), intent(in)    :: velu
  real, dimension(isd:,jsd:,:), intent(in)    :: velv
  real, dimension(isd:,jsd:,:), intent(inout) :: rit
  real, dimension(isd:,jsd:,:), intent(inout) :: riu
  
  integer :: i, j, k, kbot, kp1, m
  real    :: prev, tmp1, tmp3, active_cells

  rit   = 0.0
  riu   = 0.0
  wrk1  = 0.0 ! dT/dz 
  wrk2  = 0.0 ! dS/dz
  wrk3  = 0.0 ! N2
  wrk4  = 0.0 ! (du/dz)**2
  wrk5  = 0.0 ! 1/dzwt

  ! vertical derivative of theta and salinity 
  do k=1,nk-1
     kp1 = k+1
     do j=jsd,jed
        do i=isd,ied
           wrk5(i,j,k) = 1.0/dzwt(i,j,k)
           wrk1(i,j,k) = (theta(i,j,k)    - theta(i,j,kp1)   )*wrk5(i,j,k)
           wrk2(i,j,k) = (salinity(i,j,k) - salinity(i,j,kp1))*wrk5(i,j,k)
        enddo
     enddo
  enddo

  ! N^2 < 0 is a signature of vertically unstable regions.  
  ! N^2 = 0 at bottom of ocean bottom.
  do k=1,nk-1
     kp1 = k+1
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = -Grd%tmask(i,j,kp1)*grav*rho0r*(drhodT(i,j,k)*wrk1(i,j,k) + drhodS(i,j,k)*wrk2(i,j,k))
        enddo
     enddo
  enddo

  ! smooth N2 to reduce vertical noise
  if(smooth_n2) then
      do m=1,num_n2_smooth
         do j=jsd,jed
            do i=isd,ied
               prev = onefourth*wrk3(i,j,1) 
               kbot = Grd%kmt(i,j)
               if(kbot > 3) then
                   do k=2,kbot-2
                      tmp1        = wrk3(i,j,k)
                      wrk3(i,j,k) = prev + onehalf*wrk3(i,j,k) + onefourth*wrk3(i,j,k+1)
                      prev        = onefourth*tmp1
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute richardson numbers on bottom of U cells.
  ! compute halos to allow later averaging to rit. 
  do k=1,nk-1
     do j=jsd,jec
        do i=isd,iec
           wrk4(i,j,k) = (velu(i,j,k) - velu(i,j,k+1))**2 + (velv(i,j,k) - velv(i,j,k+1))**2 + epsln 
           tmp3 = wrk3(i,j+1,k) + wrk3(i+1,j+1,k) + wrk3(i,j,k) + wrk3(i+1,j,k)
           riu(i,j,k) = Grd%umask(i,j,k+1)*dzwt(i,j,k)*tmp3/wrk4(i,j,k)
        enddo
     enddo
  enddo

  ! smooth Richardson number in the vertical using a 1-2-1 filter
  if (smooth_ri_number) then 
     do m=1,num_ri_smooth
        do j=jsd,jec
           do i=isd,iec
              prev = onefourth*riu(i,j,1)
              kbot = Grd%kmu(i,j)
              if (kbot  > 3) then
                  do k=2,kbot-2
                     tmp1       = riu(i,j,k)
                     riu(i,j,k) = prev + onehalf*riu(i,j,k) + onefourth*riu(i,j,k+1)
                     prev       = onefourth*tmp1
                  enddo
              endif
           enddo
        enddo
     enddo
  endif

  ! compute richardson numbers on bottom of T cells as average
  ! of four nearest richardson numbers on bottom of U cells.
  ! do not consider land cells in the average.
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           active_cells = Grd%umask(i,j,k+1)   + Grd%umask(i-1,j,k+1) &
                        + Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln
           rit(i,j,k)   = (riu(i,j,k)+riu(i-1,j,k)+riu(i,j-1,k)+riu(i-1,j-1,k))/active_cells

           ! make sure no static instability exists (one that is not seen
           ! by the Richardson number).  This may happen due to
           ! horizontal averaging used in calculating the Richardson
           ! number.

           if (rit(i,j,k) > 0.0 .and. wrk3(i,j,k) > 0.0) then
              rit(i,j,k) = -10.0
           endif
        enddo
     enddo
  enddo


  ! diagnostics 
  call diagnose_3d_u(Time, Grd, id_ri_num_dudz, wrk4(:,:,:))
  call diagnose_3d(Time, Grd, id_ri_num_n2, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_rit, rit(:,:,:))
  call diagnose_3d(Time, Grd, id_riu, riu(:,:,:))

end subroutine ri_for_bgrid
! </SUBROUTINE> NAME="ri_for_bgrid"



!#######################################################################
! <SUBROUTINE NAME="ri_for_cgrid">
!
! <DESCRIPTION>
!
! Compute Richardson number assuming horizontal C-grid layout.   
!
! Richardson number rit is centered at T-cell.  
! Richardson number riu is set equal to rit, as there is no 
! separate "U-cell" when working with a Cgrid. 
!
! This calculation differs from that in ocean_vert_kpp_mom4p1 
! since here we compute N^2 using locally referenced potential density, 
! as done for tide mixing scheme and as done in ri_for_cgrid.  
!
! </DESCRIPTION>
!
subroutine ri_for_cgrid (Time, dzwt, drhodT, drhodS, theta, salinity, ut, vt, rit, riu)

  type(ocean_time_type),        intent(in)    :: Time
  real, dimension(isd:,jsd:,:), intent(in)    :: dzwt
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity 
  real, dimension(isd:,jsd:,:), intent(in)    :: ut       ! zonal velocity at t-point
  real, dimension(isd:,jsd:,:), intent(in)    :: vt       ! merid velocity at t-point
  real, dimension(isd:,jsd:,:), intent(inout) :: rit
  real, dimension(isd:,jsd:,:), intent(inout) :: riu
  
  integer :: i, j, k, kbot, kp1, m
  real    :: prev, tmp, tmp1, tmp2

  rit   = 0.0
  riu   = 0.0
  wrk1  = 0.0 ! dT/dz 
  wrk2  = 0.0 ! dS/dz
  wrk3  = 0.0 ! N2
  wrk4  = 0.0 ! (du/dz)**2
  wrk5  = 0.0 ! 1/dzwt

  ! vertical derivative of theta and salinity 
  do k=1,nk-1
     kp1 = k+1
     do j=jsd,jed
        do i=isd,ied
           wrk5(i,j,k) = 1.0/dzwt(i,j,k)
           wrk1(i,j,k) = (theta(i,j,k)-theta(i,j,kp1))*wrk5(i,j,k)
           wrk2(i,j,k) = (salinity(i,j,k)-salinity(i,j,kp1))*wrk5(i,j,k)
        enddo
     enddo
  enddo

  ! N^2 < 0 is a signature of vertically unstable regions.  
  ! N^2 = 0 at bottom of ocean bottom.
  do k=1,nk-1
     kp1 = k+1
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = -Grd%tmask(i,j,kp1)*grav*rho0r*(drhodT(i,j,k)*wrk1(i,j,k)+drhodS(i,j,k)*wrk2(i,j,k))
        enddo
     enddo
  enddo

  ! smooth N2 to reduce vertical noise
  if(smooth_n2) then
      do m=1,num_n2_smooth
         do j=jsd,jed
            do i=isd,ied
               prev = onefourth*wrk3(i,j,1) 
               kbot = Grd%kmt(i,j)
               if(kbot > 3) then
                   do k=2,kbot-2
                      tmp         = wrk3(i,j,k)
                      wrk3(i,j,k) = prev + onehalf*wrk3(i,j,k) + onefourth*wrk3(i,j,k+1)
                      prev        = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! magnitude of vertical shear on bottom of T-cell 
  wrk4(:,:,:) = 0.0
  do k=1,nk-1
     kp1 = k+1
     do j=jsc,jec
        do i=isc,iec
           tmp1        = (ut(i,j,k)-ut(i,j,kp1))*wrk5(i,j,k)
           tmp2        = (vt(i,j,k)-vt(i,j,kp1))*wrk5(i,j,k)
           wrk4(i,j,k) = tmp1**2 + tmp2**2 + epsln 
        enddo
     enddo
  enddo

  ! compute richardson number on bottom of T-cells
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           rit(i,j,k) = wrk3(i,j,k)/wrk4(i,j,k)
        enddo
     enddo
  enddo

  ! smooth Richardson number in the vertical using a 1-2-1 filter:
  if (smooth_ri_number) then 
      do m = 1,num_ri_smooth
         do j=jsc,jec
            do i=isc,iec
               prev =  onefourth*rit(i,j,1)
               kbot = Grd%kmt(i,j)
               if (kbot > 3) then
                   do k=2,kbot-2
                      tmp        =  rit(i,j,k)
                      rit(i,j,k) =  prev + onehalf*rit(i,j,k) + onefourth*rit(i,j,k+1)
                      prev       =  onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! equate riu to rit
  riu(:,:,:) = rit(:,:,:)


  ! diagnostics 
  call diagnose_3d(Time, Grd, id_ri_num_dudz, wrk4(:,:,:))
  call diagnose_3d(Time, Grd, id_ri_num_n2, wrk3(:,:,:))
  call diagnose_3d(Time, Grd, id_rit, rit(:,:,:))
  call diagnose_3d_u(Time, Grd, id_riu, riu(:,:,:))

end subroutine ri_for_cgrid
! </SUBROUTINE> NAME="ri_for_cgrid"


end module ocean_vert_util_mod

