!> Contains routines for handling FMS coupler_bc_type tracer flux structures
!! and calculating fluxes when using generic tracers

module gtracer_flux_mod

  use time_manager_mod,          only: time_type
  use coupler_types_mod,         only: coupler_1d_bc_type, coupler_2d_bc_type
  use coupler_types_mod,         only: coupler_type_spawn, coupler_type_set_diags, coupler_type_set_data
  use coupler_types_mod,         only: coupler_type_register_restarts, coupler_type_restore_state
  use coupler_types_mod,         only: ind_flux, ind_deltap, ind_kw
  use coupler_types_mod,         only: ind_pcair, ind_u10, ind_psurf
  use coupler_types_mod,         only: ind_alpha, ind_csurf, ind_sc_no
  use coupler_types_mod,         only: ind_runoff, ind_deposition
  use coupler_types_mod,         only: atmos_ocean_type_fluxes_init => coupler_types_init
  use ocean_model_mod,           only: ocean_public_type, ocean_state_type
  use ocean_model_mod,           only: ocean_model_init_sfc, ocean_model_flux_init
  use ocean_types_mod,           only: ice_ocean_boundary_type
  use atmos_ocean_fluxes_mod,    only: atmos_ocean_fluxes_init
  use mpp_domains_mod,           only: mpp_get_compute_domain
  use field_manager_mod,         only: fm_field_name_len, fm_type_name_len, fm_loop_over_list, fm_change_list
  use fm_util_mod,               only: fm_util_get_real_array
  use fms_io_mod,                only: restart_file_type, save_restart
  use constants_mod,             only: wtmair, rdgas, vonkarm
  use mpp_mod,                   only: mpp_error, FATAL

  implicit none; private

  ! Public member functions
  public :: flux_exchange_init
  public :: gas_fields_restore
  public :: gas_fields_restart
  public :: set_coupler_type_data
  public :: atmos_ocean_fluxes_calc

  character(len=*), parameter      :: mod_name = 'gtracer_flux_mod'
  real, parameter                  :: epsln=1.0e-30

  contains

  !> \brief Initialise the boundary values, including initialising and setting boundary
  !! values in FMS coupler types. Based loosely on the flux_exchange_init subroutine in
  !! flux_exchange_mod (coupler/flux_exchange.F90)
  subroutine flux_exchange_init (Time, Ocean, Ocean_state, Ice_ocean_boundary, atm_fields)
    type(time_type), intent(in)                  :: Time
    type(ocean_public_type), intent(inout)       :: Ocean
    type(ocean_state_type), pointer              :: Ocean_state
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
    type(coupler_2d_bc_type), intent(inout)      :: atm_fields

    ! local variables
    integer                   :: isc,iec,jsc,jec
    type(coupler_1d_bc_type)  :: gas_fields_atm ! tracer fields in atm
    type(coupler_1d_bc_type)  :: gas_fields_ocn ! tracer fields atop the ocean
    type(coupler_1d_bc_type)  :: gas_fluxes ! tracer fluxes between the atm and ocean

    ! Initialise FMS coupler types
    call atmos_ocean_type_fluxes_init( )
    call ocean_model_flux_init(Ocean_state)
    call atmos_ocean_fluxes_init(gas_fluxes, gas_fields_atm, gas_fields_ocn)

    call mpp_get_compute_domain(Ocean%domain, isc, iec, jsc, jec)

    ! Allocate Ice_ocean_boundary fields
    allocate ( &
      Ice_ocean_boundary%u_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%v_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%t_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%q_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%salt_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%lw_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%sw_flux_vis_dir(isc:iec,jsc:jec), &
      Ice_ocean_boundary%sw_flux_vis_dif(isc:iec,jsc:jec), &
      Ice_ocean_boundary%sw_flux_nir_dir(isc:iec,jsc:jec), &
      Ice_ocean_boundary%sw_flux_nir_dif(isc:iec,jsc:jec), &
      Ice_ocean_boundary%lprec(isc:iec,jsc:jec), &
      Ice_ocean_boundary%fprec(isc:iec,jsc:jec), &
      Ice_ocean_boundary%runoff(isc:iec,jsc:jec), &
      Ice_ocean_boundary%calving(isc:iec,jsc:jec), &
      Ice_ocean_boundary%p(isc:iec,jsc:jec), &
      Ice_ocean_boundary%aice(isc:iec,jsc:jec), &
      Ice_ocean_boundary%mh_flux(isc:iec,jsc:jec), &
      Ice_ocean_boundary%wfimelt(isc:iec,jsc:jec), &
      Ice_ocean_boundary%wfiform(isc:iec,jsc:jec), &
      Ice_ocean_boundary%licefw(isc:iec,jsc:jec), &
      Ice_ocean_boundary%liceht(isc:iec,jsc:jec), &
      Ice_ocean_boundary%wnd(isc:iec,jsc:jec))

    Ice_ocean_boundary%u_flux = 0.0
    Ice_ocean_boundary%v_flux = 0.0
    Ice_ocean_boundary%t_flux = 0.0
    Ice_ocean_boundary%q_flux = 0.0
    Ice_ocean_boundary%salt_flux = 0.0
    Ice_ocean_boundary%lw_flux = 0.0
    Ice_ocean_boundary%sw_flux_vis_dir = 0.0
    Ice_ocean_boundary%sw_flux_vis_dif = 0.0
    Ice_ocean_boundary%sw_flux_nir_dir = 0.0
    Ice_ocean_boundary%sw_flux_nir_dif = 0.0
    Ice_ocean_boundary%lprec = 0.0
    Ice_ocean_boundary%fprec = 0.0
    Ice_ocean_boundary%runoff = 0.0
    Ice_ocean_boundary%calving = 0.0
    Ice_ocean_boundary%p = 0.0
    Ice_ocean_boundary%aice = 0.0
    Ice_ocean_boundary%mh_flux = 0.0
    Ice_ocean_boundary% wfimelt = 0.0
    Ice_ocean_boundary% wfiform = 0.0
    Ice_ocean_boundary%licefw = 0.0
    Ice_ocean_boundary%liceht = 0.0
    Ice_ocean_boundary%wnd = 0.0

    ! Spawn 2D Ocean%fields FMS coupler type from 1D gas_fields_ocn
    call coupler_type_spawn(gas_fields_ocn, Ocean%fields, (/isc,isc,iec,iec/), &
      (/jsc,jsc,jec,jec/), suffix = '_ocn')
    call coupler_type_set_diags(Ocean%fields, "ocean_sfc", Ocean%axes(1:2), Time)

    ! This sets the boundary values in Ocean%fields
    call ocean_model_init_sfc(Ocean_state, Ocean)

    ! Spawn 2D Ice_ocean_boundary%fluxes FMS coupler type from 1D gas_fluxes
    ! Annoyingly, spawning doesn't copy param array, so add manually
    call coupler_type_spawn(gas_fluxes, Ice_ocean_boundary%fluxes, (/isc,isc,iec,iec/), &
      (/jsc,jsc,jec,jec/), suffix='_ice_ocn')
    call add_gas_fluxes_param(Ice_ocean_boundary%fluxes)
    call coupler_type_set_diags(Ice_ocean_boundary%fluxes, "ocean_flux", Ocean%axes(1:2), Time)

    ! Spawn 2D atm_fields FMS coupler type from 1D gas_fields_atm
    call coupler_type_spawn(gas_fields_atm, atm_fields, (/isc,isc,iec,iec/), &
          (/jsc,jsc,jec,jec/), suffix='_atm')
    call coupler_type_set_diags(atm_fields, "atmos_sfc", Ocean%axes(1:2), Time)

  end subroutine flux_exchange_init

  !> \brief Restore FMS coupler_bc_type state from ocean restart file. Based loosely on
  !! code from the coupler_init subroutine in coupler_main (coupler/coupler_main.F90 and
  !! FMScoupler)
  subroutine gas_fields_restore(Ocean)
    type(ocean_public_type), intent(inout) :: Ocean

    ! local variables
    type(restart_file_type), dimension(:), pointer :: ocn_bc_restart => NULL()
    integer                                        :: num_ocn_bc_restart

    call coupler_type_register_restarts(Ocean%fields, ocn_bc_restart, num_ocn_bc_restart, &
      Ocean%domain, ocean_restart=.true.)

    call coupler_type_restore_state(Ocean%fields, test_by_field=.true.)

  end subroutine gas_fields_restore

  !> \brief Write ocean restart file for FMS coupler_bc_type state. Based loosely on
  !! code from the coupler_restart subroutine in coupler_main (coupler/coupler_main.F90)
  subroutine gas_fields_restart(Ocean, timestamp)
    type(ocean_public_type), intent(inout) :: Ocean
    character(len=64),  optional           :: timestamp

    ! local variables
    type(restart_file_type), dimension(:), pointer :: ocn_bc_restart => NULL()
    integer                                        :: num_ocn_bc_restart
    integer                                        :: n

    call coupler_type_register_restarts(Ocean%fields, ocn_bc_restart, num_ocn_bc_restart, &
      Ocean%domain, ocean_restart=.true.)

    do n = 1, num_ocn_bc_restart
      call save_restart(ocn_bc_restart(n), timestamp)
    enddo

  end subroutine gas_fields_restart

  !> Retrieve param array from field_manager and add to FMS coupler_bc_type. This is
  !! needed because the coupler_type_spawn routine does not copy the param array into
  !! the spawned type. This routine is based on the atmos_ocean_fluxes_init subroutine
  !! in atmos_ocean_fluxes_mod (shared/coupler/atmos_ocean_fluxes.F90)
  !! https://github.com/NOAA-GFDL/FMS/blob/7f585284f1487c0629f8075be350385e6e75dd2f/coupler/atmos_ocean_fluxes.F90#L448
  subroutine add_gas_fluxes_param(gas_fluxes)
    type(coupler_2d_bc_type), intent(inout) :: gas_fluxes !< FMS coupler_bc_type to add param to

    ! local variables
    integer                                 :: n
    character(len=fm_field_name_len)        :: name
    character(len=fm_type_name_len)         :: typ
    integer                                 :: ind

    character(len=*), parameter             :: sub_name = 'add_gas_fluxes_param'
    character(len=*), parameter             :: error_header =&
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    n = 0
    do while (fm_loop_over_list('/coupler_mod/fluxes', name, typ, ind))
      if (typ .ne. 'list') then
        call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
      endif

      n = n + 1

      if (.not. fm_change_list('/coupler_mod/fluxes/' // trim(name))) then
        call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
      endif

      if (gas_fluxes%bc(n)%name .eq. name) then
        gas_fluxes%bc(n)%param => fm_util_get_real_array('param')
      else
        call mpp_error(FATAL, trim(error_header) // ' Problem setting param array pointer')
      endif
    enddo
  end subroutine add_gas_fluxes_param

  !> Set single field by name in a 2D FMS coupler type from a two-dimensional array.
  !! This is a convenience wrapper on coupler_type_set_data routine in
  !! atmos_ocean_fluxes_mod (shared/coupler/atmos_ocean_fluxes.F90)
  subroutine set_coupler_type_data(array_in, name, field_index, var, scale_factor, halo_size, idim, jdim, override)
    real, dimension(1:,1:), intent(in)          :: array_in     !< The source array for the field; its size must match the
                                                                !! size of the data being copied unless idim and jdim are
                                                                !! supplied.
    character(*), intent(in)                    :: name         !< The name of the boundary condition to set
    integer, intent(in)                         :: field_index  !< The index of the field in the boundary condition that
                                                                !! is being set.
    type(coupler_2d_bc_type), intent(inout)     :: var          !< BC_type structure with the data to set
    real, optional, intent(in)                  :: scale_factor !< A scaling factor for the data that is being added
    integer, optional, intent(in)               :: halo_size    !< The extent of the halo to copy; 0 by default
    integer, dimension(4), optional, intent(in) :: idim         !< The data and computational domain extents of the first
                                                                !! dimension of the output array in a non-decreasing list
    integer, dimension(4), optional, intent(in) :: jdim         !< The data and computational domain extents of the second
                                                                !! dimension of the output array in a non-decreasing list
    logical, optional, intent(in)               :: override     !< If true, set the override flag for the field being set

    integer                                     :: n
    logical                                     :: override_flag = .false.
    logical                                     :: fmatch = .false.

    character(len=*), parameter                 :: sub_name = 'set_coupler_type_data'
    character(len=*), parameter                 :: error_header =&
      '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    if (present(override)) override_flag = override

    do n = 1, var%num_bcs
      if ( var%bc(n)%name .eq. trim(name) ) then
        call coupler_type_set_data(array_in, n, field_index, var, &
          scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)
        if ( override_flag ) then
          var%bc(n)%field(field_index)%override = override_flag
        endif
        fmatch = .true.
        exit
      endif
    enddo
    if (.not. fmatch) then
      call mpp_error(FATAL, trim(error_header) // ' No coupler type field called '//trim(name))
    endif
  end subroutine set_coupler_type_data

  !> \brief Calculate the FMS coupler_bc_type ocean tracer fluxes. Units should be mol/m^2/s.
  !! Upward flux is positive.
  !!
  !! This routine is based on the atmos_ocean_fluxes_calc subroutine in atmos_ocean_fluxes_mod
  !! (shared/coupler/atmos_ocean_fluxes.F90 and FMScoupler) modified in the following ways:
  !! - Operate on 2D inputs, rather than 1D
  !! - Add calculation for 'air_sea_deposition' from atmos_ocean_dep_fluxes_calc subroutine
  !! - Multiply fluxes by ice_fraction input, rather than masking based on seawater input
  !! - Make tsurf input optional, as it is only used by a few implementations
  !! - Use ind_runoff rather than ind_deposition in runoff flux calculation (note, their
  !!   values are equal)
  !! - Rename gas_fields_ice to gas_fields_ocn
  subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ocn, gas_fluxes,&
    ice_fraction, isc, iec, jsc, jec, tsurf, ustar, cd_m)
    type(coupler_2d_bc_type), intent(in)     :: gas_fields_atm ! fields in atm
        !< Structure containing atmospheric surface variables that are used in the calculation
        !! of the atmosphere-ocean tracer fluxes.
    type(coupler_2d_bc_type), intent(in)     :: gas_fields_ocn ! fields atop the ocean
        !< Structure containing ocean surface variables that are used in the calculation of the
        !! atmosphere-ocean tracer fluxes.
    type(coupler_2d_bc_type), intent(inout)  :: gas_fluxes ! fluxes between the atm and ocean
        !< Structure containing the gas fluxes between the atmosphere and the ocean and
        !! parameters related to the calculation of these fluxes.
    real, intent(in)                         :: ice_fraction(isc:iec,jsc:jec) !< sea ice fraction
    integer, intent(in)                      :: isc !< The start i-index of cell centers within
                                                    !! the computational domain
    integer, intent(in)                      :: iec !< The end i-index of cell centers within the
                                                    !! computational domain
    integer, intent(in)                      :: jsc !< The start j-index of cell centers within
                                                    !! the computational domain
    integer, intent(in)                      :: jec !< The end j-index of cell centers within the
                                                    !! computational domain
    real, intent(in), optional               :: tsurf(isc:iec,jsc:jec) !< surface temperature
    real, intent(in), optional               :: ustar(isc:iec,jsc:jec) !< friction velocity, not
                                                                        !! used
    real, intent(in), optional               :: cd_m (isc:iec,jsc:jec) !< drag coefficient, not
                                                                        !! used

    ! local variables
    character(len=*), parameter   :: sub_name = 'atmos_ocean_fluxes_calc'
    character(len=*), parameter   :: error_header =&
        & '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    real, parameter                         :: permeg=1.0e-6

    integer                                 :: n
    integer                                 :: i
    integer                                 :: j
    real, dimension(:,:), allocatable       :: kw
    real, dimension(:,:), allocatable       :: cair
    character(len=128)                      :: error_string

    ! Return if no fluxes to be calculated
    if (gas_fluxes%num_bcs .le. 0) return

    if (.not. associated(gas_fluxes%bc)) then
      if (gas_fluxes%num_bcs .ne. 0) then
        call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
      else
        return
      endif
    endif

    do n = 1, gas_fluxes%num_bcs
      ! only do calculations if the flux has not been overridden
      if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then
        if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then
          if (.not. allocated(kw)) then
            allocate( kw(isc:iec,jsc:jec) )
            allocate ( cair(isc:iec,jsc:jec) )
          elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
            call mpp_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do j = jsc,jec
              do i = isc,iec
                gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                    & (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) * &
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
                cair(i,j) = &
                    gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                    gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                    gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                    & sqrt(660. / (gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j) + epsln)) *&
                    & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
                gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                    & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j)) / &
                  (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
              enddo
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
            if (.not. present(tsurf)) then
              call mpp_error(FATAL, trim(error_header) // ' Implementation ' //&
                  trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                  ' requires input tsurf')
            endif
            do j = jsc,jec
              do i = isc,iec
                gas_fluxes%bc(n)%field(ind_kw)%values(i,j) = &
                    & (1 - ice_fraction(i,j)) * gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) /&
                    & (770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) *&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i,j) *&
                    & max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln))
                !alpha: mol/m3/atm
                cair(i,j) = &
                    gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                    gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                    gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
                cair(i,j) = max(cair(i,j),0.)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                    & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
                gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                    & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                    & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
              enddo
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then
            if (.not. present(tsurf)) then
              call mpp_error(FATAL, trim(error_header) // ' Implementation ' //&
                  trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                  ' requires input tsurf')
            endif
            !f1p: not sure how to pass salinity. For now, just force at 35.
            do j = jsc,jec
              do i = isc,iec
                !calc_kw(tk,p,u10,h,vb,mw,sc_w,ustar,cd_m)
                gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                    & (1 - ice_fraction(i,j)) * calc_kw(tsurf(i,j),&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j),&
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j),&
                    & 101325./(rdgas*wtmair*1e-3*tsurf(i,j)*max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln)),&
                    & gas_fluxes%bc(n)%param(2),&
                    & gas_fluxes%bc(n)%param(1),&
                    & gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j))
                cair(i,j) =&
                    & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
                cair(i,j) = max(cair(i,j),0.)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                    & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                    & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
                gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                    & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                    & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
              enddo
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then
          if (.not. allocated(kw)) then
            allocate( kw(isc:iec,jsc:jec) )
            allocate ( cair(isc:iec,jsc:jec) )
          elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
            call mpp_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
            do j = jsc,jec
              do i = isc,iec
                kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)
                cair(i,j) =&
                    & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                    & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
              enddo
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
            do j = jsc,jec
              do i = isc,iec
                kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                    & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
                cair(i,j) =&
                    & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                    & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
              enddo
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
            do j = jsc,jec
              do i = isc,iec
                kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                    & max(0.0, gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) - gas_fluxes%bc(n)%param(2))
                cair(i,j) =&
                    & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                    & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(3)
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                    & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
              enddo
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call mpp_error(FATAL, 'Bad parameter (' // trim(error_string) //&
                & ') for air_sea_deposition for ' // trim(gas_fluxes%bc(n)%name))
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
            do j = jsc,jec
              do i = isc,iec
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                    gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
              enddo
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
            do j = jsc,jec
              do i = isc,iec
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                    gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
              enddo
            enddo
          else
            call mpp_error(FATAL, 'Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call mpp_error(FATAL, ' Bad parameter (' // trim(error_string) //&
                & ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
          endif

          if (gas_fluxes%bc(n)%implementation .eq. 'river') then
            do j = jsc,jec
              do i = isc,iec
                gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                    & gas_fields_atm%bc(n)%field(ind_runoff)%values(i,j) /&
                    & gas_fluxes%bc(n)%param(1)
              enddo
            enddo
          else
            call mpp_error(FATAL, ' Unknown implementation (' //&
                & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        else
          call mpp_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) //&
              & ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      endif
    enddo

    if (allocated(kw)) then
      deallocate(kw)
      deallocate(cair)
    endif
  end subroutine  atmos_ocean_fluxes_calc

  !> Calculate \f$k_w\f$
  !!
  !! Taken from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010)
  !!
  !! Uses equations defined in Liss[1974],
  !! \f[
  !!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
  !! \f]
  !! where \f$c_g\f$ and \f$C_l\f$ are the bulk gas and liquid concentrations, \f$H\f$
  !! is the Henry's law constant (\f$H = c_{sg}/C_{sl}\f$, where \f$c_{sg}\f$ is the
  !! equilibrium concentration in gas phase (\f$g/cm^3\f$ of air) and \f$C_{sl}\f$ is the
  !! equilibrium concentration of unionised dissolved gas in liquid phase (\f$g/cm^3\f$
  !! of water)),
  !! \f[
  !!    1/K_g = 1/k_g + H/k_l
  !! \f]
  !! and
  !! \f[
  !!    1/K_l = 1/k_l + 1/{Hk_g}
  !! \f]
  !! where \f$k_g\f$ and \f$k_l\f$ are the exchange constants for the gas and liquid
  !! phases, respectively.
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function calc_kw(tk, p, u10, h, vb, mw, sc_w, ustar, cd_m)
    real, intent(in) :: tk !< temperature at surface in kelvin
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
    real, intent(in) :: h !< Henry's law constant (\f$H=c_sg/C_sl\f$) (unitless)
    real, intent(in) :: vb !< Molar volume
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: sc_w
    real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                        !! ustar = \f$u_{10} \sqrt{C_D}\f$.
    real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                        !! ustar is provided.
                                        !! If ustar is not provided,
                                        !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

    real :: ra,rl,tc

    tc = tk-273.15
    ra = 1./max(h*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
    rl = 1./max(calc_kl(tc,u10,sc_w),epsln)
    calc_kw = 1./max(ra+rl,epsln)
  end function calc_kw

  !> Calculate \f$k_a\f$
  !!
  !! See calc_kw
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function calc_ka(t, p, mw, vb, u10, ustar, cd_m)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< molar volume
    real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
    real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                        !! ustar = \f$u_{10} \sqrt{C_D}\f$.
    real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                        !! ustar is provided.
                                        !! If ustar is not provided,
                                        !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

    real             :: sc
    real             :: ustar_t, cd_m_t

    if (.not. present(ustar)) then
      !drag coefficient
      cd_m_t = 6.1e-4 +0.63e-4*u10
      !friction velocity
      ustar_t = u10*sqrt(cd_m_t)
    else
      cd_m_t = cd_m
      ustar_t = ustar
    end if
    sc = schmidt_g(t,p,mw,vb)
    calc_ka = 1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))
  end function calc_ka

  !> Calculate \f$k_l\f$
  !!
  !! See calc_kw, and Nightingale, Global Biogeochemical Cycles, 2000
  !! (https://doi.org/10.1029/1999GB900091)
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function calc_kl(t, v, sc)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: v !< wind speed at surface in m/s
    real, intent(in) :: sc

    calc_kl = (((0.222*v**2)+0.333*v)*(max(sc,epsln)/600.)**(-0.5))/(100.*3600.)
  end function calc_kl

  !> Schmidt number of the gas in air
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function schmidt_g(t, p, mw, vb)
    real, intent(in) :: t !< temperature at surface in C
    real, intent(in) :: p !< pressure at surface in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< molar volume

    real :: d,v

    d = d_air(t,p,mw,vb)
    v = v_air(t)
    schmidt_g = v / d
  end function schmidt_g

  !> From Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007)
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function d_air(t, p, mw, vb)
    real, intent(in) :: t  !< temperature in c
    real, intent(in) :: p  !< pressure in pa
    real, intent(in) :: mw !< molecular weight (g/mol)
    real, intent(in) :: vb !< diffusion coefficient (\f$cm3/mol\f$)

    real, parameter :: ma = 28.97d0 !< molecular weight air in g/mol
    real, parameter :: va = 20.1d0  !< diffusion volume for air (\f$cm^3/mol\f$)

    real            :: pa

    ! convert p to atm
    pa = 9.8692d-6*p
    d_air = 1d-3 *&
        & (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
    ! d_air is in cm2/s convert to m2/s
    d_air = d_air * 1d-4
  end function d_air

  !> kinematic viscosity in air
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function p_air(t)
    real, intent(in) :: t

    real, parameter :: sd_0 = 1.293393662d0,&
        & sd_1 = -5.538444326d-3,&
        & sd_2 = 3.860201577d-5,&
        & sd_3 = -5.2536065d-7
    p_air = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
  end function p_air

  !> Kinematic viscosity in air (\f$m^2/s\f$
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function v_air(t)
    real, intent(in) :: t !< temperature in C
    v_air = n_air(t)/p_air(t)
  end function v_air

  !> dynamic viscosity in air
  !!
  !! This routine was copied from FMScoupler at
  !! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
  real function n_air(t)
    real, intent(in) :: t !< temperature in C

    real, parameter :: sv_0 = 1.715747771d-5,&
        & sv_1 = 4.722402075d-8,&
        & sv_2 = -3.663027156d-10,&
        & sv_3 = 1.873236686d-12,&
        & sv_4 = -8.050218737d-14
    ! in n.s/m^2 (pa.s)
    n_air = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
  end function n_air

end module gtracer_flux_mod
