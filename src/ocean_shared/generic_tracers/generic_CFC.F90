!----------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: CFC module
! This module contains the generic version of CFC Tracers and their chemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
! The chemistry calculations in this module are ported from MOM ocmip2_cfc.F90
! released in omsk_2008_03 
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 CFC
!       simulations as outlined in the CFC-HOWTO documentation,
!       revision 1.6, 1999/04/29.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/CFC/HOWTO-CFC.html
! </REFERENCE>
! <DEVELOPER_NOTES>
! nnz: 
! A. Reproducing GOLD results
! The tracers in this module reproduce the corresponding ones
! in non-generic module GOLD_OMIP2_CFC.F90 with branch tag perth_gas_fluxes_nnz.
! The reproducing of non-generic tracers in this case is sufficient evidence for the consistency of vertical diffusion
! (tracer_vertdiff) routine used in the two cases. 
! 
! B. Reproducing MOM results
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_CFC

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, stdout
  use time_manager_mod, only : time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  


  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_CFC'
  character(len=fm_string_len), parameter :: package_name   = 'generic_cfc'

  public do_generic_CFC
  public generic_CFC_register
  public generic_CFC_init
  public generic_CFC_update_from_coupler
  public generic_CFC_update_from_source
  public generic_CFC_set_boundary_values
  public generic_CFC_end

  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_CFC = .false.

  real, parameter :: epsln=1.0e-30

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_CFC_params
     real :: a1_11, a2_11, a3_11, a4_11   ! Coefficients in the calculation of the
     real :: a1_12, a2_12, a3_12, a4_12   ! CFC11 and CFC12 Schmidt numbers, in
     ! units of ND, degC-1, degC-2, degC-3.
     real :: d1_11, d2_11, d3_11, d4_11   ! Coefficients in the calculation of the
     real :: d1_12, d2_12, d3_12, d4_12   ! CFC11 and CFC12 solubilities, in units
     ! of ND, K-1, log(K)^-1, K-2.
     real :: e1_11, e2_11, e3_11          ! More coefficients in the calculation of
     real :: e1_12, e2_12, e3_12          ! the CFC11 and CFC12 solubilities, in
     ! units of PSU-1, PSU-1 K-1, PSU-1 K-2.
     real :: Rho_0
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_CFC_params


  type(generic_CFC_params) :: param

contains

  subroutine generic_CFC_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
    
  end subroutine generic_CFC_register

  ! <SUBROUTINE NAME="generic_CFC_init">
  !  <OVERVIEW>
  !   Initialize the generic CFC module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the CFC Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_CFC_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    !    call user_allocate_arrays !None for CFC module currently

  end subroutine generic_CFC_init

  subroutine user_allocate_arrays
    !Allocate all the private arrays.
    !None for CFC module currently
  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_CFC_params type
    !==============================================================    

    !=============
    !Block Starts: g_tracer_add_param
    !=============
    !Add the known experimental parameters used for calculations
    !in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and g_tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.

    call g_tracer_start_param_list(package_name)

    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !      Use coefficients given by Zheng et al (1998), JGR vol 103, C1
    !         for CFC11 and CFC12
    !-----------------------------------------------------------------------
    !    g_tracer_add_param(name   , variable   ,  default_value)
    call g_tracer_add_param('a1_11', param%a1_11,  3501.8)
    call g_tracer_add_param('a2_11', param%a2_11, -210.31)
    call g_tracer_add_param('a3_11', param%a3_11,  6.1851)
    call g_tracer_add_param('a4_11', param%a4_11, -0.07513)
    call g_tracer_add_param('a1_12', param%a1_12,  3845.4)
    call g_tracer_add_param('a2_12', param%a2_12, -228.95)
    call g_tracer_add_param('a3_12', param%a3_12,  6.1908)
    call g_tracer_add_param('a4_12', param%a4_12, -0.067430)
    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in mol/l/atm
    !      (1) for CFC11, (2) for CFC12
    !     after Warner and Weiss (1985) DSR, vol 32 for CFC11 and CFC12
    !-----------------------------------------------------------------------
    call g_tracer_add_param('d1_11', param%d1_11, -229.9261)
    call g_tracer_add_param('d2_11', param%d2_11,  319.6552)
    call g_tracer_add_param('d3_11', param%d3_11,  119.4471)
    call g_tracer_add_param('d4_11', param%d4_11, -1.39165)
    call g_tracer_add_param('e1_11', param%e1_11, -0.142382)
    call g_tracer_add_param('e2_11', param%e2_11,  0.091459)
    call g_tracer_add_param('e3_11', param%e3_11, -0.0157274)
    call g_tracer_add_param('d1_12', param%d1_12, -218.0971)
    call g_tracer_add_param('d2_12', param%d2_12,  298.9702)
    call g_tracer_add_param('d3_12', param%d3_12,  113.8049)
    call g_tracer_add_param('d4_12', param%d4_12, -1.39165)
    call g_tracer_add_param('e1_12', param%e1_12, -0.143566)
    call g_tracer_add_param('e2_12', param%e2_12,  0.091015)
    call g_tracer_add_param('e3_12', param%e3_12, -0.0153924)

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', param%Rho_0, 1035.0)

    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========



  end subroutine user_add_params

  !
  !   This is an internal sub, not a public interface.
  !   Add all the tracers to be used in this module. 
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list


    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'


    call g_tracer_start_param_list(package_name)!nnz: Does this append?
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip2_cfc.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip2_cfc.res.nc' )
    call g_tracer_add_param('IC_file'       , param%IC_file       , '')
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file=param%ice_restart_file, ocean_restart_file=param%ocean_restart_file )

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !prog_tracers: cfc_11,cfc_12 
    !diag_tracers: none
    !
    !cfc_12
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cfc_12',               &
         longname   = 'cfc_12 Concentration', &
         units      = 'mol/kg',            &
         prog       = .true.,              &
         flux_gas       = .true.,                      &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /), &
         flux_gas_restart_file  = 'ocmip2_cfc_airsea_flux.res.nc' )

    !g_cfc_11
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cfc_11',               &
         longname   = 'cfc_11 Concentration', &
         units      = 'mol/kg',            &
         prog       = .true.,              &
         flux_gas       = .true.,                      &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /), &
         flux_gas_restart_file  = 'ocmip2_cfc_airsea_flux.res.nc' )


  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_CFC_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_CFC_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_update_from_copler'
    !
    !Nothing specific to be done for CFC's
    !
    return
  end subroutine generic_CFC_update_from_coupler

  ! <SUBROUTINE NAME="generic_CFC_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_CFC_update_from_source(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    !
    !No source update for CFC's currently exit in code.
    !
    return
  end subroutine generic_CFC_update_from_source

  ! <SUBROUTINE NAME="generic_CFC_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="SST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature   
  !  </IN>
  !  <IN NAME="SSS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  ! </SUBROUTINE>

  !User must provide the calculations for these boundary values.
  subroutine generic_CFC_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha_11,alpha_12,sc_11,sc_12
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_cfc_11_field,g_cfc_12_field
    real, dimension(:,:), ALLOCATABLE :: g_cfc_11_alpha,g_cfc_11_csurf,g_cfc_12_alpha,g_cfc_12_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no_11,sc_no_12

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_set_boundary_values'


    !nnz: Can we treat these as source and move block to user_update_from_source?
    !
    !=============
    !Block Starts: Calculate the boundary values
    !=============
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'cfc_11','field',g_cfc_11_field)
    call g_tracer_get_pointer(tracer_list,'cfc_12','field',g_cfc_12_field)

    allocate(g_cfc_11_alpha(isd:ied, jsd:jed)); g_cfc_11_alpha=0.0
    allocate(g_cfc_11_csurf(isd:ied, jsd:jed)); g_cfc_11_csurf=0.0
    allocate(g_cfc_12_alpha(isd:ied, jsd:jed)); g_cfc_12_alpha=0.0
    allocate(g_cfc_12_csurf(isd:ied, jsd:jed)); g_cfc_12_csurf=0.0
    allocate(sc_no_11(isd:ied, jsd:jed))
    allocate(sc_no_12(isd:ied, jsd:jed))

    !The atmospheric code needs soluabilities in units of mol/m3/atm
    !
    !MOM
    !       The factor 1.0e+03 is for the conversion  
    !       from mol/(l * atm) to mol/(m3 * atm) 
    conv_fac = 1.0e+03
    !
    !GOLD
    !       The factor 1.e-09 converts 
    !       from mol/(l * atm) to mol/(m3 * pptv).
    !conv_fac = 1.0e-09

    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       ta = (ST(i,j) + 273.15) * 0.01 ! Why is this in dekaKelvin?
       sal = SSS(i,j) ; SST = ST(i,j)

       !---------------------------------------------------------------------
       !     Calculate solubilities
       !       Use Warner and Weiss (1985) DSR, vol 32, final result
       !       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
       !
       !       Use Bullister and Wisegavger for CCl4
       !---------------------------------------------------------------------

       !nnz: MOM hmask=grid_tmask(i,j,1), GOLD hmask=G%hmask 
       alpha_11 = conv_fac * grid_tmask(i,j,1) * &
            exp(param%d1_11 + param%d2_11/ta + param%d3_11*log(ta) + param%d4_11*ta*ta +&
            sal * ((param%e3_11 * ta + param%e2_11) * ta + param%e1_11)&
            )

       alpha_12 = conv_fac * grid_tmask(i,j,1) * &
            exp(param%d1_12 + param%d2_12/ta + param%d3_12*log(ta) + param%d4_12*ta*ta +&
            sal * ((param%e3_12 * ta + param%e2_12) * ta + param%e1_12)&
            )

       !---------------------------------------------------------------------
       !     Calculate Schmidt numbers
       !      use coefficients given by Zheng et al (1998), JGR vol 103, C1
       !---------------------------------------------------------------------
       sc_no_11(i,j) = param%a1_11 + SST * (param%a2_11 + SST * (param%a3_11 + SST * param%a4_11)) * &
            grid_tmask(i,j,1)
       sc_no_12(i,j) = param%a1_12 + SST * (param%a2_12 + SST * (param%a3_12 + SST * param%a4_12)) * &
            grid_tmask(i,j,1)

       !sc_no_term = sqrt(660.0 / (sc_11 + epsln))
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       g_cfc_11_alpha(i,j) = alpha_11              

       g_cfc_11_csurf(i,j) = g_cfc_11_field(i,j,1,taum1) *  param%Rho_0

       g_cfc_12_alpha(i,j) = alpha_12              

       g_cfc_12_csurf(i,j) = g_cfc_12_field(i,j,1,taum1) *  param%Rho_0
 
    enddo; enddo
    !=============
    !Block Ends: Calculate the boundary values
    !=============

    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'cfc_11','alpha',g_cfc_11_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc_11','csurf',g_cfc_11_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc_11','sc_no',sc_no_11,isd,jsd)

    call g_tracer_set_values(tracer_list,'cfc_12','alpha',g_cfc_12_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc_12','csurf',g_cfc_12_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc_12','sc_no',sc_no_12,isd,jsd)

    deallocate(g_cfc_11_alpha,g_cfc_11_csurf,g_cfc_12_alpha,g_cfc_12_csurf,sc_no_11,sc_no_12)

  end subroutine generic_CFC_set_boundary_values

  ! <SUBROUTINE NAME="generic_CFC_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_CFC_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_end'

  end subroutine generic_CFC_end


end module generic_CFC
