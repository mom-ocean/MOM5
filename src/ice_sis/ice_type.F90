!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_mod - maintains the sea ice data, reads/writes restarts, reads the  !
!                namelist and initializes diagnostics. - Mike Winton           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_type_mod

  use mpp_mod,          only: mpp_pe, mpp_root_pe, mpp_sum, mpp_clock_id, CLOCK_COMPONENT, &
                              CLOCK_LOOP, CLOCK_ROUTINE, stdout,input_nml_file
  use mpp_domains_mod,  only: domain2D, mpp_update_domains, CORNER, BGRID_NE
  use fms_mod,          only: file_exist, open_namelist_file, check_nml_error, write_version_number,&
                              read_data, close_file, field_exist, &
                              stderr, stdlog, error_mesg, FATAL, WARNING, NOTE, clock_flag_default
  use fms_io_mod,       only: save_restart, restore_state, query_initialized, &
                              register_restart_field, restart_file_type, set_domain, nullify_domain, &
                              parse_mask_table
  use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                              register_static_field, send_data
  use time_manager_mod, only: time_type, get_time
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: Tfreeze, radius, pi
  use ice_grid_mod,     only: set_ice_grid, t_to_uv, ice_grid_end
  use ice_grid_mod,     only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  use ice_grid_mod,     only: geo_lon, geo_lat, cell_area, sin_rot, cos_rot, wett, xb1d, yb1d
  use ice_grid_mod,     only: geo_lonv_ib, geo_latv_ib
  use ice_grid_mod,     only: x_cyclic, tripolar_grid, dtn, dte, wetv
  use ice_grid_mod,     only: reproduce_siena_201303
  use ice_thm_mod,      only: ice_thm_param, DI, DS, e_to_melt
  use ice_dyn_mod,      only: ice_dyn_param
  use constants_mod,    only: LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)
  use ice_bergs,        only: icebergs_init, icebergs_end, icebergs, icebergs_stock_pe
  use ice_bergs,        only: icebergs_save_restart
  use astronomy_mod,    only: astronomy_init, astronomy_end

  implicit none
  private

public :: ice_data_type, ice_model_init, ice_model_end, ice_stock_pe, kmelt,  &
          mom_rough_ice, heat_rough_ice, atmos_winds, hlim, slab_ice,         &
          spec_ice, verbose, ice_bulk_salin, do_ice_restore, do_ice_limit,    &
          max_ice_limit, ice_restore_timescale, do_init, h2o, heat, slp2ocean,&
          cm2_bugs, conservation_check, do_icebergs, ice_model_restart,       &
          add_diurnal_sw, channel_viscosity, smag_ocn, ssh_gravity,           &
          chan_cfl_limit, ice_data_type_chksum
public :: do_sun_angle_for_alb

public  :: id_cn, id_hi, id_hs, id_t1, id_t2, id_ts
public  :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, id_runoff,    &
           id_calving, id_runoff_hflx, id_calving_hflx,                        &
           id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m,            &
           id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna,    &
           id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, id_fay, id_fix,  &
           id_fiy, id_fcx, id_fcy, id_fwx, id_fwy, id_swdn, id_lwdn, id_sn2ic, &
           id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi,&
           id_qfres, id_qflim, id_ix_trans, id_iy_trans,                       &
           id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif,      &
           id_sw_nir_dir, id_sw_nir_dif, id_mib, id_ustar, id_vstar, id_vocean,&
           id_uocean, id_vchan, id_uchan, id_wnd

public  :: id_alb_vis_dir, id_alb_vis_dif,id_alb_nir_dir, id_alb_nir_dif

public  :: iceClock,iceClock1,iceClock2,iceClock3,iceClock4,iceClock5,iceClock6,iceClock7,iceClock8,iceClock9
public  :: iceClocka,iceClockb,iceClockc

  character(len=128) :: version = '$Id: ice_type.F90,v 20.0 2013/12/13 23:28:32 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'

  !---- id for diagnositics -------------------
  integer :: id_xb, id_xt, id_yb, id_yt, id_ct, id_xv, id_yv
  integer :: id_cn, id_hi, id_hs, id_t1, id_t2, id_ts
  integer :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain
  integer :: id_runoff, id_calving, id_runoff_hflx, id_calving_hflx
  integer :: id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m
  integer :: id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna
  integer :: id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, id_fay, id_fix
  integer :: id_fiy, id_fcx, id_fcy, id_fwx, id_fwy, id_swdn, id_lwdn, id_sn2ic
  integer :: id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi
  integer :: id_qfres, id_qflim, id_ix_trans, id_iy_trans
  integer :: id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif
  integer :: id_sw_nir_dir, id_sw_nir_dif, id_mib, id_ustar, id_vstar
  integer :: id_vocean, id_uocean, id_vchan, id_uchan
  integer :: id_alb_vis_dir, id_alb_vis_dif, id_alb_nir_dir, id_alb_nir_dif 
  integer :: id_wnd 

  !--- namelist interface --------------
  real    :: mom_rough_ice  = 1.0e-4     ! momentum same, cd10=(von_k/ln(10/z0))^2
  real    :: heat_rough_ice = 1.0e-4     ! heat roughness length
  real    :: kmelt          = 6e-5*4e6   ! ocean/ice heat flux constant
  real    :: ks             = 0.31       ! snow conductivity (W/mK)
  real    :: alb_sno        = 0.85       ! snow albedo (less if melting)
  real    :: alb_ice        = 0.5826     ! ice albedo (less if melting)
  real    :: pen_ice        = 0.3        ! part unreflected solar penetrates ice
  real    :: opt_dep_ice    = 0.67       ! ice optical depth
  real    :: t_range_melt   = 1.0        ! melt albedos scaled in over T range
  real    :: ice_bulk_salin = 0.0        ! ice bulk salinity (for ocean salt flux)
  real    :: p0             = 2.75e4     ! ice strength parameter
  real    :: c0             = 20.0       ! another ice strength parameter
  real    :: cdw            = 3.24e-3    ! water/ice drag coefficient
  real    :: wd_turn        = 25.0       ! water/ice drag turning angle
  real    :: h_lo_lim       = 0.0        ! min ice thickness for temp. calc.
  integer :: nsteps_dyn     = 432        ! dynamics steps per slow timestep
  integer :: nsteps_adv     = 8          ! advection steps per slow timestep
  integer :: num_part       = 6          ! number of ice grid partitions
                                         ! partition 1 is open water
                                         ! partitions 2 to num_part-1 are
                                         !   thickness limited ice categories
                                         ! partition num_part is unlimited ice
  logical :: atmos_winds = .true.        ! wind stress from atmosphere model over t points and has wrong sign
  logical :: slab_ice    = .false.       ! do old-style GFDL slab ice?
  logical :: spec_ice    = .false.       ! old-style GFDL slab ice with SST, ice thickness and conc. from data
  logical :: do_ice_restore  = .false.   ! restore sea-ice toward climatology
  logical :: do_ice_limit    = .false.   ! limit sea ice to max_ice_limit
  real    :: max_ice_limit   = 4.0       ! maximum sea ice height(m),
                                         ! if do_ice_limit is true
                                         ! TK: default chosen based on observed
                                         !     ice thickness data used by climate
                                         !     group, which range up to 7.31 m
  real    :: ice_restore_timescale = 5.0 ! time scale for restoring ice (days)
  logical :: conservation_check = .true. ! check for heat and h2o conservation
  logical :: slp2ocean          = .false.! apply sea level pressure to ocean surface
  logical :: cm2_bugs           = .false.! keep cm2 bugs for reproducibility        
  logical :: verbose            = .false.! control printing message, will slow model down when turn true
  logical :: do_icebergs        = .false.! call iceberg code to modify calving field
  logical :: add_diurnal_sw     = .false.! apply an additional diurnal cycle to shortwave radiation
  real    :: channel_viscosity  = 0.     ! viscosity used in one-cell wide channels to parameterize transport (m^2/s)
  real    :: smag_ocn           = 0.15   ! Smagorinksy coefficient for viscosity (dimensionless)
  real    :: ssh_gravity        = 9.81   ! Gravity parameter used in channel viscosity parameterization (m/s^2)
  real    :: chan_cfl_limit     = 0.25   ! CFL limit for channel viscosity parameterization (dimensionless)
  logical :: do_sun_angle_for_alb = .false.! find the sun angle for ocean albed in the frame of the ice model
  integer :: layout(2)          = (/0, 0/)
  integer :: io_layout(2)       = (/0, 0/)
  ! mask_table contains information for masking domain ( n_mask, layout and mask_list).
  !   A text file to specify n_mask, layout and mask_list to reduce number of processor
  !   usage by masking out some domain regions which contain all land points. 
  !   The default file name of mask_table is "INPUT/ice_mask_table". Please note that 
  !   the file name must begin with "INPUT/". The first 
  !   line of mask_table will be number of region to be masked out. The second line 
  !   of the mask_table will be the layout of the model. User need to set ice_model_nml
  !   variable layout to be the same as the second line of the mask table.
  !   The following n_mask line will be the position of the processor to be masked out.
  !   The mask_table could be created by tools check_mask. 
  !   For example the mask_table will be as following if n_mask=2, layout=4,6 and 
  !   the processor (1,2) and (3,6) will be masked out. 
  !     2
  !     4,6
  !     1,2
  !     3,6

  character(len=128) :: mask_table = "INPUT/ice_mask_table"


  namelist /ice_model_nml/ mom_rough_ice, heat_rough_ice, p0, c0, cdw, wd_turn,  &
                           kmelt, alb_sno, alb_ice, pen_ice, opt_dep_ice,        &
                           nsteps_dyn, nsteps_adv, num_part, atmos_winds,        &
                           slab_ice, spec_ice, ice_bulk_salin, layout,           &
                           do_ice_restore, do_ice_limit, max_ice_limit,          &
                           ice_restore_timescale, slp2ocean, conservation_check, &
                           t_range_melt, cm2_bugs, ks, h_lo_lim, verbose,        &
                           do_icebergs, add_diurnal_sw, io_layout, channel_viscosity,&
                           smag_ocn, ssh_gravity, chan_cfl_limit, do_sun_angle_for_alb, &
                           mask_table, reproduce_siena_201303

  logical :: do_init = .false.
  real    :: hlim(8) = (/ 0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! thickness limits 1...num_part-1
  real    :: h2o(4), heat(4) ! for conservation analysis
                             ! 1 - initial ice h2o/heat content
                             ! 2 - h2o/heat flux down at top of ice
                             ! 3 - h2o/heat flux down at bottom of ice
                             ! 4 - final ice h2o/heat content

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! This structure contains the ice model data (some used by calling routines);  !
  ! the third index is partition (1 is open water; 2 is ice cover)               !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  type ice_data_type
     type(domain2D)                     :: Domain
     type (time_type)                   :: Time_Init, Time
     type (time_type)                   :: Time_step_fast, Time_step_slow
     integer                            :: avg_kount
     logical                            :: pe
     integer, pointer, dimension(:)     :: pelist              =>NULL() ! Used for flux-exchange.
     logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
     logical, pointer, dimension(:,:,:) :: ice_mask            =>NULL() ! where ice actually is
     real,    pointer, dimension(:,:,:) :: part_size           =>NULL()
     real,    pointer, dimension(:,:,:) :: part_size_uv        =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo              =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_mom           =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_heat          =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_moist         =>NULL()
     real,    pointer, dimension(:,:,:) :: t_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: u_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: v_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: sea_lev             =>NULL()
     real,    pointer, dimension(:,:)   :: s_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ocn               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ocn               =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_u_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_v_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_u_top_bgrid    =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_v_top_bgrid    =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_t_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_q_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lw_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dir_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dif_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_nir_dir_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_nir_dif_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lh_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: lprec_top           =>NULL()
     real,    pointer, dimension(:,:,:) :: fprec_top           =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_u              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_v              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_t              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_q              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lw             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dir     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dif     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_nir_dir     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_nir_dif     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lh             =>NULL()
     real,    pointer, dimension(:,:  ) :: lprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: fprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: p_surf              =>NULL()
     real,    pointer, dimension(:,:  ) :: runoff              =>NULL()
     real,    pointer, dimension(:,:  ) :: calving             =>NULL()
     real,    pointer, dimension(:,:  ) :: runoff_hflx         =>NULL()
     real,    pointer, dimension(:,:  ) :: calving_hflx        =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_salt           =>NULL()
     real,    pointer, dimension(:,:)   :: lwdn                =>NULL()
     real,    pointer, dimension(:,:  ) :: swdn                =>NULL() ! downward long/shortwave
     real,    pointer, dimension(:,:,:) :: pen                 =>NULL()
     real,    pointer, dimension(:,:,:) :: trn                 =>NULL() ! ice optical parameters
     real,    pointer, dimension(:,:,:) :: tmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: bmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: h_snow              =>NULL()
     real,    pointer, dimension(:,:,:) :: h_ice               =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice1              =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice2              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: sig11               =>NULL()
     real,    pointer, dimension(:,:)   :: sig22               =>NULL()
     real,    pointer, dimension(:,:)   :: sig12               =>NULL()
     real,    pointer, dimension(:,:)   :: frazil              =>NULL()
     real,    pointer, dimension(:,:)   :: bheat               =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_lim_ice        =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_res_ice        =>NULL()
     real,    pointer, dimension(:,:)   :: area                =>NULL()
     real,    pointer, dimension(:,:)   :: vmask               =>NULL() ! where ice vels can be non-zero
     real,    pointer, dimension(:,:)   :: mi                  =>NULL() ! This is needed for the wave model. It is introduced here,
                                                                        ! because flux_ice_to_ocean cannot handle 3D fields. This may be
									! removed, if the information on ice thickness can be derived from 
									! eventually from h_ice outside the ice module.
     real,    pointer, dimension(:,:,:) :: wnd                 =>NULL() ! need 10m wind speed to deduce Stokes drift etc in QL_KF_17
     logical, pointer, dimension(:,:)   :: maskmap             =>NULL() ! A pointer to an array indicating which
                                                                        ! logical processors are actually used for
                                                                        ! the ocean code. The other logical
                                                                        ! processors would be all land points and
                                                                        ! are not assigned to actual processors.
                                                                        ! This need not be assigned if all logical
                                                                        ! processors are used
     integer, dimension(3)              :: axes
     type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
     type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
     type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging
     type(icebergs), pointer            :: icebergs

  end type ice_data_type

  integer :: iceClock, iceClock1, iceCLock2, iceCLock3, iceClock4, iceClock5, &
             iceClock6, iceClock7, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc
  type(restart_file_type), save :: Ice_restart

  contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_stock_pe - returns stocks of heat, water, etc. for conservation checks   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_stock_pe(Ice, index, value)

    use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

    use ice_grid_mod, only : all_avg

    type(ice_data_type)    :: Ice
    integer, intent(in) :: index
    real, intent(out)   :: value

    integer :: i, j, k
    real :: icebergs_value

    value = 0.0
    if(.not.Ice%pe) return

    select case (index)

    case (ISTOCK_WATER)

       value = sum(cell_area*all_avg(DI*Ice%h_ice(isc:iec,jsc:jec,:)+   &
            DS*Ice%h_snow(isc:iec,jsc:jec,:),  &
            Ice%part_size(isc:iec,jsc:jec,:))) &
            *4*pi*radius*radius 
  
    case (ISTOCK_HEAT)

       value = 0.0
       do k=2,km
          do j=jsc, jec
             do i=isc, iec
                if ((Ice%part_size(i,j,k)>0.0.and.Ice%h_ice(i,j,k)>0.0)) then
                   if (slab_ice) then
                      value = value - cell_area(i,j) * Ice%part_size(i,j,k)*Ice%h_ice(i,j,2)*DI*LI
                   else
                      value = value - cell_area(i,j) * Ice%part_size(i,j,k)*e_to_melt(Ice%h_snow(i,j,k), &
                                Ice%h_ice(i,j,k)/2, Ice%t_ice1(i,j,k), Ice%h_ice(i,j,k)/2, Ice%t_ice2(i,j,k))
                   end if
                end if
             end do
          end do
       end do
       value = value*4*pi*radius*radius

    case (ISTOCK_SALT)
       !No salt in the h_snow component.
       value =  sum(cell_area*all_avg(DI*Ice%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))) &
              *ice_bulk_salin*4*pi*radius*radius
    case default

       value = 0.0

    end select

    if (do_icebergs) then
      call icebergs_stock_pe(Ice%icebergs, index, icebergs_value)
      value = value + icebergs_value
    endif

  end subroutine ice_stock_pe

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_model_init - initializes ice model data, parameters and diagnostics      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_model_init (Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

    type (ice_data_type), intent(inout) :: Ice
    type (time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
    type (time_type)    , intent(in)    :: Time           ! current time
    type (time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
    type (time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

    integer           :: io, ierr, nlon, nlat, npart, unit, log_unit, k
    integer           :: sc, dy, i, j
    integer           :: id_restart, id_restart_albedo, id_restart_flux_sw
    real              :: dt_slow
    character(len=64) :: restart_file
    integer           :: stdlogunit, stdoutunit

    stdlogunit=stdlog()
    stdoutunit = stdout()
    !
    ! read namelist and write to logfile
    !
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ice_model_nml, iostat=io)
#else
    unit = open_namelist_file()
    read  (unit, ice_model_nml,iostat=io)
    call close_file (unit)
#endif
    ierr = check_nml_error(io,'ice_model_nml')
    write (stdoutunit,'(/)')
    write (stdoutunit, ice_model_nml)
    write (stdlogunit, ice_model_nml)

    call write_version_number(version, tagname)

    if (spec_ice) then
       slab_ice = .true.
       nsteps_dyn = 0
       nsteps_adv = 0
    end if
    if (slab_ice) num_part = 2 ! open water and ice ... but never in same place
    if (num_part>size(hlim(:))+1) &
         call error_mesg ('ice_model_init', 'not enough thickness limits', FATAL)

    call get_time(Time_step_slow, sc, dy); dt_slow=864e2*dy+sc

    if(file_exist(mask_table)) then
       write(stdoutunit, *) '==> NOTE from ice_model_init:  reading maskmap information from '//trim(mask_table)
       if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('ice_model_init', &
          'ice_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

       allocate(Ice%maskmap(layout(1), layout(2)))
       call parse_mask_table(mask_table, Ice%maskmap, "Ice model")
    endif

    if( ASSOCIATED(Ice%maskmap) ) then
       call set_ice_grid(Ice%domain, dt_slow, nsteps_dyn, nsteps_adv, num_part, layout, io_layout, Ice%maskmap  )
    else
       call set_ice_grid(Ice%domain, dt_slow, nsteps_dyn, nsteps_adv, num_part, layout, io_layout )
    end if
    call set_domain(domain)

    call ice_dyn_param(p0, c0, cdw, wd_turn, slab_ice)
    call ice_thm_param(alb_sno, alb_ice, pen_ice, opt_dep_ice, slab_ice, &
                       t_range_melt, cm2_bugs, ks, h_lo_lim)

    allocate ( Ice % mask     (isc:iec, jsc:jec)       , &
         Ice % ice_mask       (isc:iec, jsc:jec, km)   , &
         Ice % t_surf         (isc:iec, jsc:jec, km)   , &
         Ice % s_surf         (isc:iec, jsc:jec)       , &
         Ice % vmask          (isd:ied, jsd:jed)       , &
         Ice % sea_lev        (isd:ied, jsd:jed)       , &
         Ice % part_size      (isd:ied, jsd:jed, km)   , &
         Ice % part_size_uv   (isc:iec, jsc:jec, km)   , &
         Ice % wnd            (isc:iec, jsc:jec, km)   , &
         Ice % u_surf         (isc:iec, jsc:jec, km)   , &
         Ice % v_surf         (isc:iec, jsc:jec, km)   , &
         Ice % u_ocn          (isd:ied, jsd:jed)       , &
         Ice % v_ocn          (isd:ied, jsd:jed)       , &
         Ice % rough_mom      (isc:iec, jsc:jec, km)   , &
         Ice % rough_heat     (isc:iec, jsc:jec, km)   , &
         Ice % rough_moist    (isc:iec, jsc:jec, km)   , &
         Ice % albedo         (isc:iec, jsc:jec, km)   , &                
         Ice % albedo_vis_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_vis_dif (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dif (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u_top   (isd:ied, jsd:jed, km) ,       &
         Ice % flux_v_top         (isd:ied, jsd:jed, km) ,       &
         Ice % flux_t_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_q_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dir_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dif_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_nir_dir_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_nir_dif_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_lw_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_lh_top        (isc:iec, jsc:jec, km) ,       &
         Ice % lprec_top          (isc:iec, jsc:jec, km) ,       &
         Ice % fprec_top          (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u_top_bgrid   (isc:iec, jsc:jec, km) , &
         Ice % flux_v_top_bgrid         (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u    (isc:iec, jsc:jec ) ,       &
         Ice % flux_v          (isc:iec, jsc:jec ) ,       &
         Ice % flux_t          (isc:iec, jsc:jec ) ,       &
         Ice % flux_q          (isc:iec, jsc:jec ) ,       &
         Ice % flux_sw_vis_dir (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_vis_dif (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_nir_dir (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_nir_dif (isc:iec, jsc:jec) ,        &
         Ice % flux_lw         (isc:iec, jsc:jec ) ,       &
         Ice % flux_lh         (isc:iec, jsc:jec ) ,       &
         Ice % lprec           (isc:iec, jsc:jec ) ,       &
         Ice % fprec           (isc:iec, jsc:jec ) ,       &
         Ice % p_surf          (isc:iec, jsc:jec ) ,       &
         Ice % runoff          (isc:iec, jsc:jec ) ,       &
         Ice % calving         (isc:iec, jsc:jec ) ,       &
         Ice % runoff_hflx     (isc:iec, jsc:jec ) ,       &
         Ice % calving_hflx    (isc:iec, jsc:jec ) ,       &
         Ice % flux_salt       (isc:iec, jsc:jec ) ,       &
         Ice % lwdn            (isc:iec, jsc:jec ) ,       &
         Ice % swdn            (isc:iec, jsc:jec )         )
    allocate ( Ice % frazil (isc:iec, jsc:jec), Ice % bheat  (isc:iec, jsc:jec), &
               Ice % u_ice  (isd:ied, jsd:jed), Ice % v_ice  (isd:ied, jsd:jed), &
               Ice % sig11  (isd:ied, jsd:jed), Ice % sig22  (isd:ied, jsd:jed), &
               Ice % sig12  (isd:ied, jsd:jed)                               )
    allocate ( Ice % tmelt  (isc:iec, jsc:jec, 2:km), Ice % bmelt  (isc:iec, jsc:jec, 2:km) , &
               Ice % pen    (isc:iec, jsc:jec, 2:km), Ice % trn    (isc:iec, jsc:jec, 2:km) , &
               Ice % h_snow (isd:ied, jsd:jed, 2:km), Ice % h_ice  (isd:ied, jsd:jed, 2:km) , &
               Ice % t_ice1 (isd:ied, jsd:jed, 2:km), Ice % t_ice2 (isd:ied, jsd:jed, 2:km)   )
    allocate ( Ice % qflx_lim_ice  (isc:iec, jsc:jec) , Ice % qflx_res_ice  (isc:iec, jsc:jec)   )
    allocate ( Ice % area          (isc:iec, jsc:jec) )
    allocate ( Ice % mi            (isc:iec, jsc:jec) )

    Ice % flux_sw_vis_dir =0.
    Ice % flux_sw_vis_dif =0.
    Ice % flux_sw_nir_dir =0.
    Ice % flux_sw_nir_dif =0.
    Ice % flux_lh         =0. 
    Ice % lwdn            =0.
    Ice % swdn            =0.
    Ice % flux_u_top      =0. 
    Ice % flux_v_top      =0.
    Ice % flux_u_top_bgrid=0. 
    Ice % flux_v_top_bgrid=0.
    Ice % sea_lev         =0.
    Ice % part_size       =0.
    Ice % wnd             =0.
    Ice % u_ocn           =0.
    Ice % v_ocn           =0.
    Ice % u_ice           =0.
    Ice % v_ice           =0.
    Ice % sig11           =0.
    Ice % sig12           =0.
    Ice % sig22           =0.
    Ice % h_snow          =0.
    Ice % h_ice           =0.
    Ice % t_ice1          =0.
    Ice % t_ice2          =0.
    Ice % area            = cell_area * 4*PI*RADIUS*RADIUS
    Ice % mi              =0.
    Ice % u_surf          =0.
    Ice % v_surf          =0.
    Ice % s_surf          =0.
    Ice % flux_t_top      =0.
    Ice % flux_q_top      =0.
    Ice % flux_lw_top     =0.
    Ice % flux_sw_vis_dir_top =0.
    Ice % flux_sw_vis_dif_top =0.
    Ice % flux_sw_nir_dir_top =0.
    Ice % flux_sw_nir_dif_top =0.
    Ice % flux_lh_top     =0.
    Ice % lprec_top       =0.
    Ice % fprec_top       =0.
    Ice % flux_salt       =0.
    Ice % pen             =0.
    Ice % trn             =0.
    Ice % bheat           =0.


    do j = jsc, jec
       do i = isc, iec
          if( wett(i,j) > 0.5 ) then
             Ice % mask(i,j) = .true.
          else
             Ice % mask(i,j) = .false.
          end if
          if( wetv(i,j) > 0.5 ) then
             Ice % vmask(i,j) = 1.
          else
             Ice % vmask(i,j) = 0.
          end if
       enddo
    enddo
    if(reproduce_siena_201303) then
       call mpp_update_domains(Ice%vmask, domain=domain )
    else
       call mpp_update_domains(Ice%vmask, domain=domain, position=CORNER )
    endif

    Ice % Time           = Time
    Ice % Time_Init      = Time_Init
    Ice % Time_step_fast = Time_step_fast
    Ice % Time_step_slow = Time_step_slow

    Ice % avg_kount      = 0
    !
    ! read restart
    !
    !! Need to create new restart version including the needed albedos
    !  that have been added ??

    !--- 
    restart_file = 'ice_model.res.nc'
    id_restart = register_restart_field(Ice_restart, restart_file, 'part_size', Ice%part_size, domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'albedo',    Ice%albedo,    domain=domain)
    id_restart_albedo = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dir', Ice%albedo_vis_dir, &
                                               domain=domain, mandatory=.false.)
    id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dir', Ice%albedo_nir_dir, &
                                               domain=domain, mandatory=.false.)
    id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dif', Ice%albedo_vis_dif, &
                                               domain=domain, mandatory=.false.)
    id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dif', Ice%albedo_nir_dif, &
                                               domain=domain, mandatory=.false.)
    id_restart = register_restart_field(Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,        domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,       domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'rough_moist', Ice%rough_moist,      domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 't_surf',      Ice%t_surf,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'h_snow',      Ice%h_snow(:,:,2:km), domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'h_ice',       Ice%h_ice(:,:,2:km),  domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 't_ice1',      Ice%t_ice1(:,:,2:km), domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 't_ice2',      Ice%t_ice2(:,:,2:km), domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'u_ice',       Ice%u_ice,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'v_ice',       Ice%v_ice,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'sig11',       Ice%sig11,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'sig22',       Ice%sig22,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'sig12',       Ice%sig12,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_u',      Ice%flux_u,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_v',      Ice%flux_v,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_t',      Ice%flux_t,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_q',      Ice%flux_q,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_salt',   Ice%flux_salt,        domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'flux_lw',     Ice%flux_lw,          domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'lprec',       Ice%lprec,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'fprec',       Ice%fprec,            domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'runoff',      Ice%runoff,           domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'calving',     Ice%calving,          domain=domain)
    id_restart = register_restart_field(Ice_restart, restart_file, 'runoff_hflx', Ice%runoff_hflx,      domain=domain, mandatory=.false.)
    id_restart = register_restart_field(Ice_restart, restart_file, 'calving_hflx',Ice%calving_hflx,     domain=domain, mandatory=.false.)
    id_restart = register_restart_field(Ice_restart, restart_file, 'p_surf',      Ice%p_surf,           domain=domain)
    id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dir', Ice%flux_sw_vis_dir, &
                                                domain=domain, mandatory=.false.)    
    id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dif', Ice%flux_sw_vis_dif, &
                                                domain=domain, mandatory=.false.)
    id_restart_flux_sw = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dir', Ice%flux_sw_nir_dir, &
                                                domain=domain, mandatory=.false.)
    id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dif', Ice%flux_sw_nir_dif, &
                                                domain=domain, mandatory=.false.)
!
!        Total SW flux is broken into 4 components in Nalanda,
!        it was a single component preNalanda.
!
    restart_file = 'INPUT/ice_model.res.nc'
    if (file_exist(restart_file)) then
       call restore_state(Ice_restart)
       if( .NOT. query_initialized(Ice_restart, id_restart_albedo)) then
          Ice % albedo_vis_dir    = 0.0         
          Ice % albedo_nir_dir    = 0.0       
          Ice % albedo_vis_dif    = 0.0       
          Ice % albedo_nir_dif    = 0.0       
       endif

       if(.NOT. query_initialized(Ice_restart, id_restart_flux_sw) ) then 
          call error_mesg ('ice_model_init', &
            'Restart file does not contain flux_sw_* subcomponents!', WARNING)
          ! for compatibility with preN restarts, we should check for total SW flux 
          if(field_exist(restart_file,'flux_sw')) then
             call read_data( restart_file, 'flux_sw',Ice%flux_sw_vis_dir , domain ) 
             ! simplest way to brake the total flux to 4 components
             Ice % flux_sw_vis_dir = Ice%flux_sw_vis_dir / 4
             Ice % flux_sw_vis_dif = Ice%flux_sw_vis_dir
             Ice % flux_sw_nir_dir = Ice%flux_sw_vis_dir
             Ice % flux_sw_nir_dif = Ice%flux_sw_vis_dir
          else
             call error_mesg ('ice_model_init', &
                  'Restart file does not contain flux_sw total or its components!', FATAL)
          endif
       endif

       !--- update to data domain
       call mpp_update_domains(Ice%part_size, Domain)
       call mpp_update_domains(Ice%h_snow(:,:,2:km), Domain )
       call mpp_update_domains(Ice%h_ice (:,:,2:km), Domain )
       call mpp_update_domains(Ice%t_ice1(:,:,2:km), Domain )
       call mpp_update_domains(Ice%t_ice2(:,:,2:km), Domain )
       if(reproduce_siena_201303) then
          call mpp_update_domains(Ice%u_ice, Domain )
          call mpp_update_domains(Ice%v_ice, Domain )
       else
          call mpp_update_domains(Ice%u_ice, Ice%v_ice, Domain, gridtype=BGRID_NE )
       endif
       call mpp_update_domains(Ice%sig11, Domain )
       call mpp_update_domains(Ice%sig22, Domain )
       call mpp_update_domains(Ice%sig12, Domain )
    else ! no restart => no ice
       Ice % part_size    = 0.0
       !   where (Ice%mask) Ice % part_size (:,:,1) = 1.0  - flux_exchange won't allow
       Ice % part_size (:,:,1) = 1.0
       Ice % albedo            = 0.0
       Ice % albedo_vis_dir    = 0.0         
       Ice % albedo_nir_dir    = 0.0       
       Ice % albedo_vis_dif    = 0.0       
       Ice % albedo_nir_dif    = 0.0       
       Ice % rough_mom         = mom_rough_ice
       Ice % rough_heat        = heat_rough_ice
       Ice % rough_moist       = heat_rough_ice
       Ice % t_surf            = Tfreeze-5.0
       Ice % h_snow            = 0.0
       Ice % h_ice             = 0.0
       Ice % t_ice1            = -5.0
       Ice % t_ice2            = -5.0
       Ice % u_ice             = 0.0
       Ice % v_ice             = 0.0
       Ice % sig11             = 0.0
       Ice % sig22             = 0.0
       Ice % sig12             = 0.0
       Ice % flux_u            = 0.0 
       Ice % flux_v            = 0.0
       Ice % flux_t            = 0.0 
       Ice % flux_q            = 0.0
       Ice % flux_lw           = 0.0
       Ice % flux_salt         = 0.0 
       Ice % lprec             = 0.0
       Ice % fprec             = 0.0
       Ice % p_surf            = 0.0
       Ice % runoff            = 0.0
       Ice % calving           = 0.0
       Ice % runoff_hflx       = 0.0
       Ice % calving_hflx      = 0.0
       Ice % frazil            = 0.0
       Ice % flux_sw_vis_dir   = 0.0 
       Ice % flux_sw_vis_dif   = 0.0 
       Ice % flux_sw_nir_dir   = 0.0 
       Ice % flux_sw_nir_dif   = 0.0 
       do_init = .true. ! done in ice_model
    end if

    Ice % tmelt       = 0.0
    Ice % bmelt       = 0.0

    Ice % qflx_lim_ice = 0.0
    Ice % qflx_res_ice = 0.0

    Ice%part_size_uv(:,:,1) = 1.0
    do k=2,km
       call t_to_uv(Ice%part_size(:,:,k),Ice%part_size_uv(:,:,k))
       Ice%part_size_uv (:,:,1) = Ice%part_size_uv(:,:,1)-Ice%part_size_uv (:,:,k)
    end do

    call ice_diagnostics_init(Ice)
    !Balaji
    iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    iceClock1 = mpp_clock_id( 'Ice: bot to top', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    iceClock2 = mpp_clock_id( 'Ice: update slow (dn)', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    iceClock7 = mpp_clock_id( '  Ice: slow: conservation check', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock4 = mpp_clock_id( '  Ice: slow: dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClocka = mpp_clock_id( '       slow: ice_dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClockb = mpp_clock_id( '       slow: comm/cut check ', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClockc = mpp_clock_id( '       slow: diags', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock5 = mpp_clock_id( '  Ice: slow: thermodynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock6 = mpp_clock_id( '  Ice: slow: restore/limit', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock8 = mpp_clock_id( '  Ice: slow: salt to ocean', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock9 = mpp_clock_id( '  Ice: slow: thermodyn diags', flags=clock_flag_default, grain=CLOCK_LOOP )
    iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )
 
    ! Initialize icebergs
    if (do_icebergs) call icebergs_init(Ice%icebergs, &
             im, jm, layout, io_layout, Ice%axes(1:2), Ice%maskmap, x_cyclic, tripolar_grid, &
             dt_slow, Time, geo_lonv_ib, geo_latv_ib, wett, dtn, dte, cell_area, cos_rot, sin_rot )

   if (add_diurnal_sw .or. do_sun_angle_for_alb) call astronomy_init
   call nullify_domain()

  end subroutine ice_model_init

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_model_end - writes the restart file                                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_model_end (Ice)
    type (ice_data_type), intent(inout) :: Ice
    integer           :: k

    integer           :: unit
    character(len=22) :: restart='RESTART/ice_model.res'

    if (conservation_check) then
       do k=1,4
          call mpp_sum(h2o(k))
          call mpp_sum(heat(k))
       end do
       if (mpp_pe()==mpp_root_pe()) then
          print *
          print '(a10,5a13)',   'ICE MODEL ','   AT START  ', &
               ' TOP FLUX DN.', &
               ' BOT FLUX DN.', &
               '   AT END    ', &
               '   ERROR     '
          print '(a10,5es13.5)','WATER     ', h2o , h2o (4)-(h2o (1)+h2o (2)-h2o (3))
          print '(a10,5es13.5)','HEAT      ', heat, heat(4)-(heat(1)+heat(2)-heat(3))
          print *
       end if
    end if

    call ice_model_restart()

    !--- release memory ------------------------------------------------
    call ice_grid_end()

    deallocate(Ice % mask, Ice % ice_mask, Ice % t_surf, Ice % s_surf, Ice % sea_lev )
    deallocate(Ice % vmask)
    deallocate(Ice % part_size, Ice % part_size_uv, Ice % u_surf, Ice % v_surf )
    deallocate(Ice % u_ocn, Ice % v_ocn ,  Ice % rough_mom, Ice % rough_heat )
    deallocate(Ice % rough_moist, Ice % albedo, Ice % flux_u_top, Ice % flux_v_top )
    deallocate(Ice % flux_u_top_bgrid, Ice % flux_v_top_bgrid )
    deallocate(Ice % flux_t_top, Ice % flux_q_top, Ice % flux_lw_top )
    deallocate(Ice % flux_lh_top, Ice % lprec_top, Ice % fprec_top, Ice % flux_u )
    deallocate(Ice % flux_v, Ice % flux_t, Ice % flux_q, Ice % flux_lw )
    deallocate(Ice % flux_lh, Ice % lprec, Ice % fprec, Ice % p_surf, Ice % runoff ) 
    deallocate(Ice % calving, Ice % runoff_hflx, Ice % calving_hflx )
    deallocate(Ice % flux_salt)
    deallocate( Ice % lwdn)
    deallocate( Ice % swdn)
    deallocate( Ice % frazil )
    deallocate(Ice % bheat, Ice % u_ice, Ice % v_ice, Ice % sig11, Ice % sig22 )
    deallocate(Ice % sig12, Ice % tmelt, Ice % bmelt, Ice % pen, Ice % trn )
    deallocate(Ice % h_snow, Ice % h_ice, Ice % t_ice1, Ice % t_ice2  )
    deallocate(Ice % qflx_lim_ice, Ice % qflx_res_ice )
    deallocate(Ice % flux_sw_vis_dir, Ice % flux_sw_vis_dif )
    deallocate(Ice % flux_sw_nir_dir, Ice % flux_sw_nir_dif )
    deallocate(Ice % wnd )

    ! End icebergs
    if (do_icebergs) call icebergs_end(Ice%icebergs)

    if (add_diurnal_sw .or. do_sun_angle_for_alb) call astronomy_end

  end subroutine ice_model_end

  !#######################################################################
  ! <SUBROUTINE NAME="ice_model_restart">
  ! <DESCRIPTION>
  !  Write out restart files registered through register_restart_file
  ! </DESCRIPTION>
  subroutine ice_model_restart(Ice, time_stamp)
    type (ice_data_type),     intent(inout), optional :: Ice
    character(len=*),         intent(in), optional :: time_stamp

    call save_restart(Ice_restart, time_stamp)
   !call icebergs_save_restart(Ice%icebergs)
   ! This should go here but since "Ice" is not available we have to
   ! rely on the restart written via ice_model_end() -AJA

  end subroutine ice_model_restart
  ! </SUBROUTINE>
  !#######################################################################

  subroutine ice_diagnostics_init(Ice)
    type (ice_data_type), intent(inout) :: Ice

    real, parameter       :: missing = -1e34
    integer, dimension(2) :: axt, axv, axtv, axvt
    integer, dimension(3) :: axt2
    integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
    logical               :: sent

    !
    ! diagnostics MUST use a domain without halos otherwise same as the
    ! regular domain:  Domain (see ice_grid.f90)
    !
    id_xv = diag_axis_init('xv', xb1d(2:im+1), 'degrees_E', 'X','longitude', set_name='ice', Domain2=Domain )
    id_yv = diag_axis_init('yv', yb1d(2:jm+1), 'degrees_N', 'Y','latitude',  set_name='ice', Domain2=Domain )
    id_xb = diag_axis_init('xb', xb1d, 'degrees_E', 'X', 'longitude', set_name='ice', Domain2=Domain )
    id_yb = diag_axis_init('yb', yb1d, 'degrees_N', 'Y', 'latitude', set_name='ice', Domain2=Domain )
    id_xt = diag_axis_init('xt', (xb1d(1:im)+xb1d(2:im+1))/2, 'degrees_E', 'X', &
            'longitude',set_name='ice',edges=id_xb,Domain2=Domain)
    id_yt = diag_axis_init('yt', (yb1d(1:jm)+yb1d(2:jm+1))/2, 'degrees_N', 'Y', &
            'latitude',set_name='ice', edges=id_yb,Domain2=Domain)
    id_ct = diag_axis_init('ct', hlim(1:num_part-1), 'meters','Z', 'thickness')
    axv  = (/ id_xv, id_yv       /)
    axt  = (/ id_xt, id_yt       /)
    axt2 = (/ id_xt, id_yt, id_ct/)
    Ice%axes(:) = axt2(:)
    axtv = (/ id_xt, id_yv /); ! for north faces of t-cells
    axvt = (/ id_xv, id_yt /); ! for east  faces of t-cells

    id_sin_rot   = register_static_field('ice_model', 'SINROT', axt,              &
                   '-SINROT,COSROT points north', 'none')
    id_cos_rot   = register_static_field('ice_model', 'COSROT', axt,              &
                   'COSROT,SINROT points east','none')
    id_geo_lon   = register_static_field('ice_model', 'GEOLON', axt, 'longitude', &
                   'degrees')
    id_geo_lat   = register_static_field('ice_model', 'GEOLAT', axt, 'latitude',  &
                   'degrees')
    id_cell_area = register_static_field('ice_model', 'CELL_AREA', axt,           &
                   'cell area', 'sphere')
    id_ext       = register_diag_field('ice_model', 'MOI', axt, Ice%Time,   &
                   'ice modeled', '0 or 1', missing_value=missing)
    if (id_ext > 0 ) then
       call error_mesg ('ice_model_init', &
            'Diagnostic MOI has been renamed EXT.  Change your diag_table.', WARNING)
    else
       id_ext = register_diag_field('ice_model', 'EXT', axt, Ice%Time, &
                'ice modeled', '0 or 1', missing_value=missing)
    end if
    id_wnd      = register_diag_field('ice_model', 'wnd ', axt, Ice%Time,   &
                   'wind speed', 'm/s', missing_value=missing)
    id_mi       = register_diag_field('ice_model', 'MI', axt, Ice%Time,                  &
                 'ice mass', 'kg/m^2', missing_value=missing)
    id_mib      = register_diag_field('ice_model', 'MIB', axt, Ice%Time,                  &
                 'ice + bergs mass', 'kg/m^2', missing_value=missing)
    id_cn       = register_diag_field('ice_model', 'CN', axt2, Ice%Time,                 &
                 'ice concentration', '0-1', missing_value=missing)
    id_hs       = register_diag_field('ice_model', 'HS', axt, Ice%Time,                  &
                 'snow thickness', 'm-snow', missing_value=missing)
    id_hi       = register_diag_field('ice_model', 'HI', axt, Ice%Time,                  &
                 'ice thickness', 'm-ice', missing_value=missing)
    id_t1       = register_diag_field('ice_model', 'T1', axt, Ice%Time,                  &
                 'upper ice layer temperature', 'C',  missing_value=missing)
    id_t2       = register_diag_field('ice_model', 'T2', axt, Ice%Time,                  &
                 'lower ice layer temperature', 'C',  missing_value=missing)
    id_ts       = register_diag_field('ice_model', 'TS', axt, Ice%Time,                  &
                 'surface temperature', 'C', missing_value=missing)
    id_sh       = register_diag_field('ice_model','SH' ,axt, Ice%Time,                   &
                 'sensible heat flux', 'W/m^2',  missing_value=missing)
    id_lh       = register_diag_field('ice_model','LH' ,axt, Ice%Time,                   &
                 'latent heat flux', 'W/m^2', missing_value=missing)
    id_sw       = register_diag_field('ice_model','SW' ,axt, Ice%Time,                   &
                 'short wave heat flux', 'W/m^2', missing_value=missing)
    id_lw       = register_diag_field('ice_model','LW' ,axt, Ice%Time,                   &
                 'long wave heat flux over ice', 'W/m^2', missing_value=missing)
    id_snofl    = register_diag_field('ice_model','SNOWFL' ,axt, Ice%Time,               &
                 'rate of snow fall', 'kg/(m^2*s)', missing_value=missing)
    id_rain     = register_diag_field('ice_model','RAIN' ,axt, Ice%Time,                 &
                 'rate of rain fall', 'kg/(m^2*s)', missing_value=missing)
    id_runoff   = register_diag_field('ice_model','RUNOFF' ,axt, Ice%Time,               &
                 'liquid runoff', 'kg/(m^2*s)', missing_value=missing)
    id_calving  = register_diag_field('ice_model','CALVING',axt, Ice%Time,               &
                 'frozen runoff', 'kg/(m^2*s)', missing_value=missing)
    id_runoff_hflx   = register_diag_field('ice_model','RUNOFF_HFLX' ,axt, Ice%Time,               &
                 'liquid runoff sensible heat flux', 'W/m^2', missing_value=missing)
    id_calving_hflx  = register_diag_field('ice_model','CALVING_HFLX',axt, Ice%Time,               &
                 'frozen runoff sensible heat flux', 'W/m^2', missing_value=missing)
    id_evap     = register_diag_field('ice_model','EVAP',axt, Ice%Time,                  &
                 'evaporation', 'kg/(m^2*s)', missing_value=missing)
    id_saltf    = register_diag_field('ice_model','SALTF' ,axt, Ice%Time,                &
                 'ice to ocean salt flux', 'kg/(m^2*s)', missing_value=missing)
    id_sn2ic    = register_diag_field('ice_model','SN2IC'  ,axt,Ice%Time,                &
                 'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)
    id_tmelt    = register_diag_field('ice_model','TMELT'  ,axt, Ice%Time,               &
                 'upper surface melting energy flux', 'W/m^2', missing_value=missing)
    id_bmelt    = register_diag_field('ice_model','BMELT'  ,axt, Ice%Time,               &
                 'bottom surface melting energy flux', 'W/m^2', missing_value=missing)
    id_bheat    = register_diag_field('ice_model','BHEAT'  ,axt, Ice%Time,               &
                 'ocean to ice heat flux', 'W/m^2', missing_value=missing)
    id_e2m      = register_diag_field('ice_model','E2MELT' ,axt, Ice%Time,               &
                 'heat needed to melt ice', 'J/m^2', missing_value=missing)
    id_frazil   = register_diag_field('ice_model','FRAZIL' ,axt, Ice%Time,               &
                 'energy flux of frazil formation', 'W/m^2', missing_value=missing)
    id_alb      = register_diag_field('ice_model','ALB',axt, Ice%Time,                   &
                 'surface albedo','0-1', missing_value=missing )
    id_alb_vis_dir = register_diag_field('ice_model','alb_vis_dir',axt, Ice%Time,                &
                 'ice surface albedo vis_dir','0-1', missing_value=missing )
    id_alb_vis_dif = register_diag_field('ice_model','alb_vis_dif',axt, Ice%Time,                &
                 'ice surface albedo vis_dif','0-1', missing_value=missing )
    id_alb_nir_dir = register_diag_field('ice_model','alb_nir_dir',axt, Ice%Time,                &
                 'ice surface albedo nir_dir','0-1', missing_value=missing )
    id_alb_nir_dif = register_diag_field('ice_model','alb_nir_dif',axt, Ice%Time,                &
                 'ice surface albedo nir_dif','0-1', missing_value=missing )
    id_xprt     = register_diag_field('ice_model','XPRT',axt, Ice%Time,                  &
                 'frozen water transport convergence', 'kg/(m^2*yr)', missing_value=missing)
    id_lsrc     = register_diag_field('ice_model','LSRC', axt, Ice%Time,                 &
                 'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
    id_lsnk     = register_diag_field('ice_model','LSNK',axt, Ice%Time,                  &
                 'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
    id_bsnk     = register_diag_field('ice_model','BSNK',axt, Ice%Time,                  &
                 'frozen water local bottom sink', 'kg/(m^2*yr)', missing_value=missing)
    id_qfres    = register_diag_field('ice_model', 'QFLX_RESTORE_ICE', axt, Ice%Time,    &
                 'Ice Restoring heat flux', 'W/m^2', missing_value=missing)
    id_qflim    = register_diag_field('ice_model', 'QFLX_LIMIT_ICE', axt, Ice%Time,      &
                 'Ice Limit heat flux', 'W/m^2', missing_value=missing)
    id_strna    = register_diag_field('ice_model','STRAIN_ANGLE', axt,Ice%Time,          &
                 'strain angle', 'none', missing_value=missing)
    id_sigi     = register_diag_field('ice_model','SIGI' ,axt, Ice%Time,                 &
                 'first stress invariant', 'none', missing_value=missing)
    id_sigii    = register_diag_field('ice_model','SIGII' ,axt, Ice%Time,                &
                 'second stress invariant', 'none', missing_value=missing)
    id_stren    = register_diag_field('ice_model','STRENGTH' ,axt, Ice%Time,             &
                 'ice strength', 'Pa*m', missing_value=missing)
    id_ui       = register_diag_field('ice_model', 'UI', axv, Ice%Time,                  &
                 'ice velocity - x component', 'm/s', missing_value=missing)
    id_vi       = register_diag_field('ice_model', 'VI', axv, Ice%Time,                  &
                 'ice velocity - y component', 'm/s', missing_value=missing)
    id_ix_trans = register_diag_field('ice_model', 'IX_TRANS', axvt, Ice%Time,           &
                 'x-direction ice transport', 'kg/s', missing_value=missing)
    id_iy_trans = register_diag_field('ice_model', 'IY_TRANS', axtv, Ice%Time,           &
                 'y-direction ice transport', 'kg/s', missing_value=missing)
    id_fax      = register_diag_field('ice_model', 'FA_X', axv, Ice%Time,                &
                 'air stress on ice - x component', 'Pa', missing_value=missing)
    id_fay      = register_diag_field('ice_model', 'FA_Y', axv, Ice%Time,                &
                 'air stress on ice - y component', 'Pa', missing_value=missing)
    id_fix      = register_diag_field('ice_model', 'FI_X', axv, Ice%Time,                &
                 'ice internal stress - x component', 'Pa', missing_value=missing)
    id_fiy      = register_diag_field('ice_model', 'FI_Y', axv, Ice%Time,                &
                 'ice internal stress - y component', 'Pa', missing_value=missing)
    id_fcx      = register_diag_field('ice_model', 'FC_X', axv, Ice%Time,                &
                 'coriolis force - x component', 'Pa', missing_value=missing)
    id_fcy      = register_diag_field('ice_model', 'FC_Y', axv, Ice%Time,                &
                 'coriolis force - y component', 'Pa', missing_value=missing)
    id_fwx      = register_diag_field('ice_model', 'FW_X', axv, Ice%Time,                &
                 'water stress on ice - x component', 'Pa', missing_value=missing)
    id_fwy      = register_diag_field('ice_model', 'FW_Y', axv, Ice%Time,                &
                 'water stress on ice - y component', 'Pa', missing_value=missing)
    id_uo       = register_diag_field('ice_model', 'UO', axv, Ice%Time,                  &
                 'surface current - x component', 'm/s', missing_value=missing)
    id_vo       = register_diag_field('ice_model', 'VO', axv, Ice%Time,                  &
                 'surface current - y component', 'm/s', missing_value=missing)
    id_sw_vis   = register_diag_field('ice_model','SW_VIS' ,axt, Ice%Time,               &
                 'visible short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_dir   = register_diag_field('ice_model','SW_DIR' ,axt, Ice%Time,               &
                 'direct short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_dif   = register_diag_field('ice_model','SW_DIF' ,axt, Ice%Time,               &
                 'diffuse short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_vis_dir = register_diag_field('ice_model','SW_VIS_DIR' ,axt, Ice%Time,         &
                 'visible direct short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_vis_dif = register_diag_field('ice_model','SW_VIS_DIF' ,axt, Ice%Time,         &
                 'visible diffuse short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_nir_dir = register_diag_field('ice_model','SW_NIR_DIR' ,axt, Ice%Time,         &
                 'near IR direct short wave heat flux', 'W/m^2', missing_value=missing)
    id_sw_nir_dif = register_diag_field('ice_model','SW_NIR_DIF' ,axt, Ice%Time,         &
                 'near IR diffuse short wave heat flux', 'W/m^2', missing_value=missing)
    id_ustar    = register_diag_field('ice_model', 'U_STAR', axvt, Ice%Time,              &
                 'channel transport velocity - x component', 'm/s', missing_value=missing)
    id_vstar    = register_diag_field('ice_model', 'V_STAR', axtv, Ice%Time,              &
                 'channel transport velocity - y component', 'm/s', missing_value=missing)
    id_uocean   = register_diag_field('ice_model', 'U_CHAN_OCN', axvt, Ice%Time,          &
                 'ocean component of channel transport - x', 'm/s', missing_value=missing)
    id_vocean   = register_diag_field('ice_model', 'V_CHAN_OCN', axtv, Ice%Time,          &
                 'ocean component of channel transport - y', 'm/s', missing_value=missing)
    id_uchan    = register_diag_field('ice_model', 'U_CHAN_VISC', axvt, Ice%Time,         &
                 'viscous component of channel transport - x', 'm/s', missing_value=missing)
    id_vchan    = register_diag_field('ice_model', 'V_CHAN_VISC', axtv, Ice%Time,         &
                 'viscous component of channel transport - y', 'm/s', missing_value=missing)

    !
    ! diagnostics for quantities produced outside the ice model
    !
    id_swdn  = register_diag_field('ice_model','SWDN' ,axt, Ice%Time,       &
               'downward shortwave flux', 'W/m^2', missing_value=missing)
    id_lwdn  = register_diag_field('ice_model','LWDN' ,axt, Ice%Time,       &
               'downward longwave flux', 'W/m^2', missing_value=missing)
    id_ta    = register_diag_field('ice_model', 'TA', axt, Ice%Time,        &
               'surface air temperature', 'C', missing_value=missing)
    id_slp   = register_diag_field('ice_model', 'SLP', axt, Ice%Time,       &
               'sea level pressure', 'Pa', missing_value=missing)
    id_sst   = register_diag_field('ice_model', 'SST', axt, Ice%Time,       &
               'sea surface temperature', 'deg-C', missing_value=missing)
    id_sss   = register_diag_field('ice_model', 'SSS', axt, Ice%Time,       &
               'sea surface salinity', 'psu', missing_value=missing)
    id_ssh   = register_diag_field('ice_model', 'SSH', axt, Ice%Time,       &
               'sea surface height', 'm', missing_value=missing)
    id_obi   = register_diag_field('ice_model', 'OBI', axt, Ice%Time,       &
         'ice observed', '0 or 1', missing_value=missing)

    if (id_sin_rot>0)   sent=send_data(id_sin_rot, sin_rot(isc:iec,jsc:jec), Ice%Time);
    if (id_cos_rot>0)   sent=send_data(id_cos_rot, cos_rot(isc:iec,jsc:jec), Ice%Time);
    if (id_geo_lon>0)   sent=send_data(id_geo_lon, geo_lon, Ice%Time);
    if (id_geo_lat>0)   sent=send_data(id_geo_lat, geo_lat, Ice%Time);
    if (id_cell_area>0) sent=send_data(id_cell_area, cell_area, Ice%Time);

  end subroutine ice_diagnostics_init

subroutine ice_data_type_chksum(id, timestep, data_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_data_type), intent(in) :: data_type
 integer ::   n, m, outunit

    outunit = stdout()
100 FORMAT("   CHECKSUM::",A32," = ",Z20)
    write(outunit,*) "BEGIN CHECKSUM(ice_data_type):: ", id, timestep
    write(outunit,100) 'ice_data_type%part_size          ',mpp_chksum(data_type%part_size(isc:iec,jsc:jec,:)          )
    write(outunit,100) 'ice_data_type%part_size_uv       ',mpp_chksum(data_type%part_size_uv(isc:iec,jsc:jec,:)       )
    write(outunit,100) 'ice_data_type%albedo             ',mpp_chksum(data_type%albedo(isc:iec,jsc:jec,:)             )
    write(outunit,100) 'ice_data_type%albedo_vis_dir     ',mpp_chksum(data_type%albedo_vis_dir(isc:iec,jsc:jec,:)     )
    write(outunit,100) 'ice_data_type%albedo_nir_dir     ',mpp_chksum(data_type%albedo_nir_dir(isc:iec,jsc:jec,:)     )
    write(outunit,100) 'ice_data_type%albedo_vis_dif     ',mpp_chksum(data_type%albedo_vis_dif(isc:iec,jsc:jec,:)     )
    write(outunit,100) 'ice_data_type%albedo_nir_dif     ',mpp_chksum(data_type%albedo_nir_dif(isc:iec,jsc:jec,:)     )
    write(outunit,100) 'ice_data_type%rough_mom          ',mpp_chksum(data_type%rough_mom(isc:iec,jsc:jec,:)          )
    write(outunit,100) 'ice_data_type%rough_heat         ',mpp_chksum(data_type%rough_heat(isc:iec,jsc:jec,:)         )
    write(outunit,100) 'ice_data_type%rough_moist        ',mpp_chksum(data_type%rough_moist(isc:iec,jsc:jec,:)        )
    write(outunit,100) 'ice_data_type%flux_u             ',mpp_chksum(data_type%flux_u(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%flux_v             ',mpp_chksum(data_type%flux_v(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%flux_t             ',mpp_chksum(data_type%flux_t(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%flux_q             ',mpp_chksum(data_type%flux_q(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%flux_lw            ',mpp_chksum(data_type%flux_lw(isc:iec,jsc:jec)              )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dir    ',mpp_chksum(data_type%flux_sw_vis_dir(isc:iec,jsc:jec)      )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dif    ',mpp_chksum(data_type%flux_sw_vis_dif(isc:iec,jsc:jec)      )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dir    ',mpp_chksum(data_type%flux_sw_nir_dir(isc:iec,jsc:jec)      )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dif    ',mpp_chksum(data_type%flux_sw_nir_dif(isc:iec,jsc:jec)      )
    write(outunit,100) 'ice_data_type%lprec              ',mpp_chksum(data_type%lprec(isc:iec,jsc:jec)                )
    write(outunit,100) 'ice_data_type%fprec              ',mpp_chksum(data_type%fprec(isc:iec,jsc:jec)                )
    write(outunit,100) 'ice_data_type%p_surf             ',mpp_chksum(data_type%p_surf(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%runoff             ',mpp_chksum(data_type%runoff(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%calving            ',mpp_chksum(data_type%calving(isc:iec,jsc:jec)              )
    write(outunit,100) 'ice_data_type%flux_salt          ',mpp_chksum(data_type%flux_salt(isc:iec,jsc:jec)            )
    write(outunit,100) 'ice_data_type%h_snow             ',mpp_chksum(data_type%h_snow(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%h_ice              ',mpp_chksum(data_type%h_ice(isc:iec,jsc:jec,:) )
    write(outunit,100) 'ice_data_type%t_ice1             ',mpp_chksum(data_type%t_ice1(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%t_ice2             ',mpp_chksum(data_type%t_ice2(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%u_ice              ',mpp_chksum(data_type%u_ice(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%v_ice              ',mpp_chksum(data_type%v_ice(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%sig11              ',mpp_chksum(data_type%sig11(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%sig22              ',mpp_chksum(data_type%sig22(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%sig12              ',mpp_chksum(data_type%sig12(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%frazil             ',mpp_chksum(data_type%frazil(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%qflx_lim_ice       ',mpp_chksum(data_type%qflx_lim_ice(isc:iec,jsc:jec))
    write(outunit,100) 'ice_data_type%qflx_res_ice       ',mpp_chksum(data_type%qflx_res_ice(isc:iec,jsc:jec))
    write(outunit,*) '   ======The following are not restart variables======'
    write(outunit,100) 'ice_data_type%u_surf             ',mpp_chksum(data_type%u_surf(isc:iec,jsc:jec,:)             )
    write(outunit,100) 'ice_data_type%v_surf             ',mpp_chksum(data_type%v_surf(isc:iec,jsc:jec,:)             )
    write(outunit,100) 'ice_data_type%sea_lev            ',mpp_chksum(data_type%sea_lev(isc:iec,jsc:jec)              )
    write(outunit,100) 'ice_data_type%s_surf             ',mpp_chksum(data_type%s_surf(isc:iec,jsc:jec)               )
    write(outunit,100) 'ice_data_type%u_ocn              ',mpp_chksum(data_type%u_ocn(isc:iec,jsc:jec)                )
    write(outunit,100) 'ice_data_type%v_ocn              ',mpp_chksum(data_type%v_ocn(isc:iec,jsc:jec)                )
    write(outunit,100) 'ice_data_type%flux_u_top         ',mpp_chksum(data_type%flux_u_top(isc:iec,jsc:jec,:)         )
    write(outunit,100) 'ice_data_type%flux_v_top         ',mpp_chksum(data_type%flux_v_top(isc:iec,jsc:jec,:)         )
    write(outunit,100) 'ice_data_type%flux_t_top         ',mpp_chksum(data_type%flux_t_top(isc:iec,jsc:jec,:)         )
    write(outunit,100) 'ice_data_type%flux_q_top         ',mpp_chksum(data_type%flux_q_top(isc:iec,jsc:jec,:)         )
    write(outunit,100) 'ice_data_type%flux_lw_top        ',mpp_chksum(data_type%flux_lw_top(isc:iec,jsc:jec,:)        )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dir_top',mpp_chksum(data_type%flux_sw_vis_dir_top(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%flux_sw_vis_dif_top',mpp_chksum(data_type%flux_sw_vis_dif_top(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%flux_sw_nir_dir_top',mpp_chksum(data_type%flux_sw_nir_dir_top(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%flux_sw_nir_dif_top',mpp_chksum(data_type%flux_sw_nir_dif_top(isc:iec,jsc:jec,:))
    write(outunit,100) 'ice_data_type%flux_lh_top        ',mpp_chksum(data_type%flux_lh_top(isc:iec,jsc:jec,:)        )
    write(outunit,100) 'ice_data_type%lprec_top          ',mpp_chksum(data_type%lprec_top(isc:iec,jsc:jec,:)          )
    write(outunit,100) 'ice_data_type%fprec_top          ',mpp_chksum(data_type%fprec_top(isc:iec,jsc:jec,:)          )
    write(outunit,100) 'ice_data_type%flux_lh            ',mpp_chksum(data_type%flux_lh(isc:iec,jsc:jec)              )
    write(outunit,100) 'ice_data_type%lwdn               ',mpp_chksum(data_type%lwdn(isc:iec,jsc:jec)                 )
    write(outunit,100) 'ice_data_type%swdn               ',mpp_chksum(data_type%swdn(isc:iec,jsc:jec)                 )
    write(outunit,100) 'ice_data_type%pen                ',mpp_chksum(data_type%pen(isc:iec,jsc:jec,:)                )
    write(outunit,100) 'ice_data_type%trn                ',mpp_chksum(data_type%trn(isc:iec,jsc:jec,:)                )
    write(outunit,100) 'ice_data_type%tmelt              ',mpp_chksum(data_type%tmelt(isc:iec,jsc:jec,:)              )
    write(outunit,100) 'ice_data_type%bmelt              ',mpp_chksum(data_type%bmelt(isc:iec,jsc:jec,:)              )
    write(outunit,100) 'ice_data_type%bheat              ',mpp_chksum(data_type%bheat(isc:iec,jsc:jec))

    do n = 1, data_type%ocean_fields%num_bcs  !{
       do m = 1, data_type%ocean_fields%bc(n)%num_fields  !{
          !write(outunit,101) 'ice%', m, n, mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
          write(outunit,101) 'ice%',trim(data_type%ocean_fields%bc(n)%name), &
               trim(data_type%ocean_fields%bc(n)%field(m)%name), &
               mpp_chksum(data_type%ocean_fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("   CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum

end module ice_type_mod
