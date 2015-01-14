module auscom_ice_mod
!
!==========================================================================

! AusCOM modules:
use auscom_ice_parameters_mod

!MOM4 modules: 
use fms_mod,         only: file_exist
use fms_mod,         only: open_namelist_file, close_file, check_nml_error
use mpp_domains_mod, only: mpp_update_domains, domain2d
use mpp_mod,         only: mpp_broadcast, mpp_pe, mpp_npes
use mpp_mod,         only: stdlog, stdout
use mpp_domains_mod, only: mpp_get_compute_domain, &
                           mpp_get_data_domain
use ocean_types_mod, only: ice_ocean_boundary_type, &
                           ocean_public_type
use ocean_types_mod, only: ocean_prog_tracer_type, &
                           ocean_diag_tracer_type, &
                           ocean_Thickness_type, &
                           ocean_grid_type, &
                           ocean_time_type, &
                           ocean_time_steps_type

use constants_mod, only: rho_cp &  !(J/m^3/deg C) rho_cp == rho0*cp_ocean
                        ,rho0r  &  !(m^3/kg)  rho0r == 1.0/rho0
                        ,rho0   &  ! 1.035e3 kg/m^3 (rho_sw)
                        ,hlf    &  ! 3.34e5  J/kg
                        ,cp_ocean  ! 3989.24495292815 J/kg/deg

#if defined(UNIT_TESTING)
use dump_field, only: dump_field_2d, dump_field_close
#endif

implicit none

real, dimension(:,:), allocatable :: &
          AQICE,     & ! accumulated ice form/melt flux (C*m)
          QICE         ! total column cooling from ice formation (in C*m) 
                       ! (ICEFLUX can be either time accumulated or averaged!)  

integer :: ATIME       ! accumulated time for ice formation calculation (s)  

integer :: iisc, iiec, jjsc, jjec
!integer :: iisd, iied, jjsd, jjed, iisc, iiec, jjsc, jjec

real :: salref = 34.7           !psu, NOT as used in POP (i.e., msu)
real :: salice = 4.0            !psu, re-set as in CICE

real :: dt_ocean

contains

!============================================================================
subroutine auscom_ice_init(domain, Time_steps)

implicit none

type(domain2d) :: domain
type(ocean_time_steps_type) :: Time_steps
integer :: ioun, io_status, ierr

! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,auscom_ice_nml,IOSTAT=io_status)
  write (stdlog(),auscom_ice_nml)
  write (stdout(),*) 'auscom_ice_nml='
  write (stdout(),'(/)')
  write (stdout(),auscom_ice_nml)
  ierr = check_nml_error(io_status,'auscom_ice_nml')
  call close_file (ioun)

call mpp_get_compute_domain(domain,iisc,iiec,jjsc,jjec)

allocate(QICE (iisc:iiec,jjsc:jjec));     QICE(:,:)  = 0
allocate(AQICE(iisc:iiec,jjsc:jjec));     AQICE(:,:) = 0

ATIME = 0


dt_ocean = Time_steps%dtts

end subroutine auscom_ice_init

!============================================================================
subroutine auscom_ice_formation_new(Time,T_prog,Thickness, Frazil)
!
!calculate ice formation (including top layer and 'sub-surface' layers)
!adapted from POP (ice.F, v 1.11, dated 2003/01/28)
!

! Calculate Frazil as a diagnstic tracer. This should allow us to balance heat budgets

! Using pointers to array sections. Can use array syntax. Much clearer. POTICE now an ARRAY

implicit none

type(ocean_time_type), intent(in)           :: Time
type(ocean_prog_tracer_type), intent(inout),target :: T_prog(:)
type(ocean_Thickness_type), intent(in),target      :: Thickness
type(ocean_diag_tracer_type), intent(inout),target :: Frazil
!
!
real :: cp_over_lhfusion = rho_cp/hlf/1000.0
       !cp_over_lhfusion = rho_sw*cp_sw/(latent_heat_fusion*rho_fw)
       !   (/deg C)        (J/m^3/deg C)      (J/kg)     (1000kg/m^3)
real :: epsilon = 1.0e-20       !as in POP

integer :: i, k, index_temp, index_salt, taup1
integer :: num_prog_tracers

real, pointer :: PTR_TEMP(:,:),PTR_SALT(:,:),PTR_THICK(:,:), PTR_FRAZIL(:,:)

real, dimension(:,:),allocatable ::  POTICE, TEMP_BEFORE

taup1=Time%taup1

num_prog_tracers = size(T_prog)
do i=1, num_prog_tracers
   if (T_prog(i)%name == 'temp') index_temp = i
   if (T_prog(i)%name == 'salt') index_salt = i
enddo

allocate(POTICE(iisc:iiec,jjsc:jjec),TEMP_BEFORE(iisc:iiec,jjsc:jjec))

!  initialize flux to zero
QICE   = 0.0

!-----------------------------------------------------------------------
! compute frazil ice formation for sub-surface layers. if ice
! forms in lower layers but layers above are warm - the heat is
! used to melt the ice. the ice formation occurs at salinity, Si.
! this volume is replaced with an equal volume at the salinity of
! the layer above. the total ice heat flux is accumulated.
!
! (???) WARNING: unless a monotone advection scheme is in place,
! advective errors could lead to temps that are far below freezing
! in some locations and this scheme will form lots of ice.
! ice formation should be limited to the top layer (kmxice=1)
! if the advection scheme is not monotone.
!-----------------------------------------------------------------------

POTICE = 0.0

do k = kmxice, 1, -1
   PTR_TEMP => T_prog(index_temp)%field(iisc:iiec,jjsc:jjec,k,taup1)
   PTR_SALT => T_prog(index_salt)%field(iisc:iiec,jjsc:jjec,k,taup1)
   PTR_THICK => Thickness%dzt(iisc:iiec,jjsc:jjec,k)
   PTR_FRAZIL => Frazil%field(iisc:iiec,jjsc:jjec,k)

#if defined(UNIT_TESTING)
    call dump_field_2d('ice_formation.input.temp', mpp_pe(), PTR_TEMP)
    call dump_field_2d('ice_formation.input.salt', mpp_pe(), PTR_SALT)
    call dump_field_2d('ice_formation.input.thickness', mpp_pe(), PTR_THICK)
    call dump_field_2d('ice_formation.input.frazil', mpp_pe(), PTR_FRAZIL)
#endif

   ! Save original temperature

   TEMP_BEFORE = PTR_TEMP


   !***
   !*** potice is the potential amount of ice formation
   !*** (potice>0) or melting (potice<0) in layer k
   !***

   POTICE = (-0.054* PTR_SALT - PTR_TEMP) * PTR_THICK

   !***
   !*** if POTICE < 0, use the heat to melt any ice from lower layers
   !*** if POTICE > 0, keep on freezing (QICE < 0)
   !***

   POTICE = max(POTICE, QICE)

   !***
   !*** adjust tracer values based on freeze/melt
   !***

   PTR_TEMP = PTR_TEMP + POTICE / PTR_THICK

   if (iceform_adj_salt) then 
      PTR_SALT = PTR_SALT + (salref - salice) * POTICE * cp_over_lhfusion / PTR_THICK
   endif

   QICE = QICE - POTICE  ! accumulate (vertically) freezing potential


 !-----------------------------------------------------------------------
 !
 ! let any residual heat in the upper layer melt previously formed ice
 !
 !-----------------------------------------------------------------------

  if ( k == 1 ) then
     AQICE = AQICE + QICE       !in degC*m

 !-----------------------------------------------------------------------
 !
 ! recalculate freezing potential based on adjusted T, S.  only interested 
 ! in melt potential now (POTICE < 0) -------- use this melt to offset any 
 ! accumulated freezing (AQICE < 0) and adjust T,S to reflect this melting 
 !
 !-----------------------------------------------------------------------

     POTICE = (-0.054* PTR_SALT - PTR_TEMP) * PTR_THICK

 !surface layer ice form[>0]/melt[<0] potential

     POTICE = max(POTICE, AQICE)
 !tricky......
     AQICE = AQICE - POTICE

 !
 ! adjust T, S again:
 !

    PTR_TEMP = PTR_TEMP + POTICE / PTR_THICK

    if (iceform_adj_salt) then 
       PTR_SALT = PTR_SALT + (salref - salice) * POTICE * cp_over_lhfusion / PTR_THICK
    endif

   

 !-----------------------------------------------------------------------
 !
 ! compute the heat flux for input to the sea ice model. 
 ! note that either:
 !   AQICE<0 and SST=Tfreeze => qflux = -AQICE >0 (net ice made), or  
 !   AQICE=0 and SST>Tfreeze => qflux ~ Tfreeze-SST <0 (melt potential)
 !
 ! REMARK: qflux ( QICE below) is needed only if it's time to communicate
 !         with oasis. however, in order to isolate this subroutine from 
 !         any coupling details, qflux is computed at every ocean time 
 !         step (but only the last calculated one is sent to oasis/ice). 
 !         Note qflux will be converted to W/m^2 later in get_o2i_fields
 !         when coupling actually happens.
 !-----------------------------------------------------------------------

    POTICE = (-0.054* PTR_SALT - PTR_TEMP) * PTR_THICK
    QICE = POTICE - AQICE 

  endif

 !-----------------------------------------------------------------------
 !  Calculate total change to Temperature and therefore heat due to Frazil
 !-----------------------------------------------------------------------

  PTR_FRAZIL = rho_cp * frazil_factor * ( PTR_TEMP - TEMP_BEFORE) * PTR_THICK 

#if defined(UNIT_TESTING)
    call dump_field_2d('ice_formation.output.frazil', mpp_pe(), PTR_FRAZIL)
    call dump_field_2d('ice_formation.output.qice', mpp_pe(), QICE)
#endif

enddo  !loop k = kmxice, 1, -1

ATIME = ATIME + dt_ocean

#if defined(UNIT_TESTING)
    call dump_field_2d('ice_formation.output.atime', mpp_pe(), &
                       reshape((/ float(ATIME) /), (/ 1, 1 /)))
#endif

deallocate(TEMP_BEFORE,POTICE)

end subroutine auscom_ice_formation_new

!=========================================================================
subroutine auscom_ice_heatflux_new(Ocean_sfc)
!
! convert the iceflux m*degC into W/m^2 as required by ice model and set 
! the accumulated AQICE and ATIME back to zero. 
! called once in a coupling interval before calling get_o2i_fields.

implicit none
type (ocean_public_type) :: Ocean_sfc

integer :: i, j

do j = jjsc, jjec
do i = iisc, iiec

  Ocean_sfc%frazil(i,j) = QICE(i,j) * frazil_factor * rho_cp / float(ATIME)
  ! here ATIME can be accumulated time (seconds) in a coupling interval 
  ! or just the ocean time step, depending how often the ice_formation_new
  ! is called  (21/07/2008, E. Hunke suggested POP calls ice_foramtion
  !                         once per coupling interval) 

enddo
enddo

AQICE = 0.0
ATIME = 0

#if defined(UNIT_TESTING)
    call dump_field_2d('ice_heatflux.output.frazil', mpp_pe(), Ocean_sfc%frazil)
#endif

end subroutine auscom_ice_heatflux_new

!===========================================================================
end module auscom_ice_mod
