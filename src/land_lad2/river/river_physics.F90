module river_physics_mod 

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! For the full text of the GNU General Public License,               
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Kirsten Findell </CONTACT> 
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang </CONTACT> 

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file
#endif

  use mpp_mod,         only : mpp_sync_self, mpp_send, mpp_recv, EVENT_RECV, EVENT_SEND
  use mpp_mod,         only : mpp_npes, mpp_error, FATAL, mpp_get_current_pelist
  use mpp_mod,         only : mpp_root_pe, mpp_pe, mpp_max
  use mpp_mod,         only : COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : ZERO, NINETY, MINUS_NINETY, mpp_update_domains 
  use mpp_domains_mod, only : mpp_get_compute_domains
  use mpp_domains_mod, only : mpp_get_num_overlap, mpp_get_overlap
  use mpp_domains_mod, only : mpp_get_update_size, mpp_get_update_pelist
  use fms_mod,         only : stdlog, write_version_number
  use fms_mod,         only : close_file, check_nml_error, file_exist
  use diag_manager_mod,only : register_diag_field, send_data
  use river_type_mod,  only : river_type, Leo_Mad_trios, NO_RIVER_FLAG
  use lake_mod,        only : large_dyn_small_stat
  use lake_tile_mod,   only : num_l
  use constants_mod,   only : tfreeze, hlf, DENS_H2O
  use land_debug_mod,  only : set_current_point, is_watch_cell

  implicit none
  private

  real    :: missing = -1.e8

!--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river_physics.F90,v 20.0 2013/12/13 23:29:43 fms Exp $'
  character(len=128) :: tagname = '$Name: tikal $'


! ---- public interfaces -----------------------------------------------------

  public :: river_physics_init, river_physics_step, river_impedes_lake, river_impedes_large_lake

!----------------------------------------------------------------------
  real               :: clw = 4218.
  real               :: csw = 2106.
  real,    parameter :: sec_in_day = 86400.
  integer :: num_lake_lev

! ---- namelist interface
  character*6 :: algor = 'linear'
  real :: lake_outflow_frac_ceiling = 1.e20
  real :: lake_sfc_w_min = -1.e20
  real :: storage_threshold_for_melt = 1.
  real :: storage_threshold_for_diag = 1.e6
  logical :: ice_frac_from_sfc = .false.
  logical :: use_lake_area_bug = .false.
  logical :: zero_frac_bug     = .false. ! it TRUE, reverts to quebec (buggy)
      ! behavior, where the discharge points with zero land fraction were
      ! missed, resulting in water non-conservation
  real :: ice_frac_factor = 0.
  logical :: prohibit_cold_ice_outflow = .TRUE. ! default retrieves old behavior,
      ! to activate bugfix, set it to FALSE
  logical :: lockstep = .false. ! set to true to recognize that lake level falls
      ! in lockstep with river when integrating river storage
  logical :: river_impedes_lake = .false.
  logical :: river_impedes_large_lake = .true.

  namelist /river_physics_nml/ algor, lake_outflow_frac_ceiling, &
                               lake_sfc_w_min, storage_threshold_for_melt, &
                               storage_threshold_for_diag, &
                               ice_frac_from_sfc, ice_frac_factor, &
                               use_lake_area_bug, zero_frac_bug, &
                               prohibit_cold_ice_outflow, lockstep, &
                               river_impedes_lake, river_impedes_large_lake

  integer, parameter, dimension(8) :: di=(/1,1,0,-1,-1,-1,0,1/)
  integer, parameter, dimension(8) :: dj=(/0,-1,-1,-1,0,1,1,1/)
  integer                          :: isc, iec, jsc, jec  ! compute domain
  integer                          :: isd, ied, jsd, jed  ! data domain
  integer                          :: maxtravel
  integer                          :: npes
  integer                          :: num_species
 
  type comm_type
     integer          :: count
     integer          :: pe
     integer, pointer :: i(:) => NULL()
     integer, pointer :: j(:) => NULL()
     integer, pointer :: k(:) => NULL()
  end type comm_type

  type halo_update_type
     type(comm_type), pointer :: send(:) => NULL();
     type(comm_type), pointer :: recv(:) => NULL();
  end type halo_update_type

  type(halo_update_type),  allocatable :: halo_update(:)
  integer                              :: nsend_update, nrecv_update
  real, dimension(:),      allocatable :: send_buffer, recv_buffer
  logical, dimension(:,:), allocatable :: in_domain
  integer, dimension(:,:), allocatable :: nlev

  ! ---- diag field IDs
  integer :: id_temp, id_ice
  integer :: id_temp_old, id_ice_old ! for compatibility with older diagTables

contains

!#######################################################################

  subroutine river_physics_init(River, domain, id_lon, id_lat )
    type(river_type), intent(inout) :: River
    type(domain2d),   intent(inout) :: domain
    integer, intent(in) :: id_lon, id_lat ! diag field IDs

    integer                         :: unit, io_status, ierr
    integer                         :: i, j


!--- read namelist -------------------------------------------------
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=river_physics_nml, iostat=io_status)
    ierr = check_nml_error(io_status, 'river_physics_nml')
#else
    if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr = 1;
      do while ( ierr/=0 )
         read  (unit, river_physics_nml, iostat=io_status, end=10)
         ierr = check_nml_error(io_status,'river_physics_nml')
      enddo
10    continue
      call close_file (unit)
    endif
#endif

!--- write version and namelist info to logfile --------------------
    call write_version_number(version, tagname)
    unit=stdlog()
    write (unit, river_physics_nml)  

    npes     = mpp_npes()

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    num_lake_lev = num_l
    num_species = size(River%outflow_c,3)
    maxtravel = maxval(River%travel)
    call mpp_max(maxtravel)

!--- set up the halo update 
    call setup_halo_update(River, domain)

    River%i_tocell = NO_RIVER_FLAG; River%j_tocell = NO_RIVER_FLAG
    do j = jsc, jec
       do i = isc, iec
          if(River%tocell(i,j) > 0) then
             River%i_tocell(i,j) = i + di(River%tocell(i,j))
             River%j_tocell(i,j) = j + dj(River%tocell(i,j))
          end if
       end do
    end do
 
    ! ---- register diagnostic fields
    id_ice  = register_diag_field ( 'river', 'rv_ice', (/id_lon, id_lat/), &
         River%Time, 'river ice mass fraction', '-', missing_value=missing, &
         mask_variant=.TRUE. )
    id_temp = register_diag_field ( 'river', 'rv_T', (/id_lon, id_lat/), &
         River%Time, 'river temperature', 'K', missing_value=missing, &
         mask_variant=.TRUE. )

    id_ice_old = register_diag_field ( 'river', 'ice', (/id_lon, id_lat/), &
         River%Time, 'obsolete, pls use rv_ice', '-', missing_value=missing, &
         mask_variant=.TRUE. )
    id_temp_old = register_diag_field ( 'river', 'temp', (/id_lon, id_lat/), &
         River%Time, 'obsolete, pls use rv_T', 'K', missing_value=missing, &
         mask_variant=.TRUE. )
  end subroutine river_physics_init

!#####################################################################

  subroutine river_physics_step(River, cur_travel, &
         lake_sfc_A, lake_sfc_bot, lake_depth_sill, lake_width_sill, &
         lake_whole_area, lake_T, lake_wl, lake_ws )

    type(river_type),     intent(inout) :: River
    integer,                 intent(in) :: cur_travel

    real, dimension(isd:ied,jsd:jed), intent(in) :: &
                             lake_sfc_A, lake_sfc_bot
    real, dimension(isd:ied,jsd:jed,num_lake_lev), intent(inout) :: &
                             lake_wl, lake_ws
    real, dimension(isc:iec,jsc:jec), intent(in) :: &
                lake_depth_sill, lake_width_sill, lake_whole_area
    real, dimension(isc:iec,jsc:jec,num_lake_lev), intent(inout) :: &
                             lake_T
! ---- local vars ----------------------------------------------------------
    integer   :: i, j, to_i, to_j, i_species, lev
    real      :: Q0, dQ_dV, dh_dQ, avail, out_frac, qmelt
    real      :: liq_to_flow, ice_to_flow, liq_this_lev, ice_this_lev
    real      :: lake_area, h, ql, qs, qh, qt, h0, t_scale
    real      :: influx
    real      :: influx_c(River%num_species)
    real      :: v_r_d(River%num_species-River%num_c+1:River%num_species)
    real      :: conc(1:River%num_species)
    logical, dimension(isc:iec,jsc:jec) :: &
         diag_mask ! mask of valid ice and temperature values fo diagnostics
    real, dimension(isc:iec,jsc:jec) :: &
         ice, temperature ! variables for diag output (were in River_type) 
    logical :: used ! flag returned by the send_data

    ! invalidate diag_mask everywhere
    diag_mask = .FALSE.

    ! do for all cells at current number of steps from river mouth
    do j = jsc, jec 
      do i = isc, iec
        call set_current_point(i,j,1) ! for debug output
        if (River%travel(i,j)==cur_travel.and.&
            ((.not.zero_frac_bug).or.(River%landfrac(i,j).gt.0))) then
            ! if zero_frac_bug is FALSE, the second line of condition is
            ! always TRUE, so we revert to bugfix
            ! if zero_frac_bug is TRUE, the second line is simply
            ! River%landfrac(i,j).gt.0, so we get quebec (buggy) condition

            ! FIRST COMPUTE LAKE MASS BALANCE (FROM INFLOC AND INFLOW TO LAKE_OUTFLOW)
          
            lake_area = lake_sfc_A(i,j)
            influx   =(River%inflow  (i,j)  +River%infloc  (i,j))  *DENS_H2O*River%dt_slow
            influx_c =(River%inflow_c(i,j,:)+River%infloc_c(i,j,:))*DENS_H2O*River%dt_slow
            if (River%tocell(i,j).eq.0 .and. River%landfrac(i,j).ge.1.) then
                ! terminal, all-land cell (must have lake)
                h = (clw*lake_wl(i,j,1)+csw*lake_ws(i,j,1))*(lake_T(i,j,1)-tfreeze)
                lake_wl(i,j,1) = lake_wl(i,j,1) + (influx-influx_c(1))/lake_area
                lake_ws(i,j,1) = lake_ws(i,j,1) +         influx_c(1) /lake_area
                lake_T (i,j,1) = tfreeze + &
                   (h+influx_c(2)/lake_area)/(clw*lake_wl(i,j,1)+csw*lake_ws(i,j,1))
                ! LAKE_SFC_C(I,J,:) = LAKE_SFC_C(I,J,:) + INFLUX_C / LAKE_AREA
              else
                ! non-terminal all-land cell (possible lake), or terminal coastal cell (possible lake)
                if (lake_area.gt.0.) then
                     if (is_watch_cell()) then
                          write(*,*) 'lake_wl(1):', lake_wl(i,j,1)
                          write(*,*) 'lake_ws(1):', lake_ws(i,j,1)
                          write(*,*) 'lake_T (1):', lake_T (i,j,1)
                     endif
                     h = (clw*lake_wl(i,j,1)+csw*lake_ws(i,j,1))*(lake_T(i,j,1)-tfreeze)
                     lake_wl(i,j,1) = lake_wl(i,j,1) + (influx-influx_c(1))/lake_area
                     lake_ws(i,j,1) = lake_ws(i,j,1) +         influx_c(1) /lake_area
                     lake_T (i,j,1) = tfreeze + &
                        (h+influx_c(2)/lake_area)/(clw*lake_wl(i,j,1)+csw*lake_ws(i,j,1))
                     if (is_watch_cell()) then
                          write(*,*) 'lake_wl(1):', lake_wl(i,j,1)
                          write(*,*) 'lake_ws(1):', lake_ws(i,j,1)
                          write(*,*) 'lake_T (1):', lake_T (i,j,1)
                     endif
                     ! LAKE_SFC_C(I,J,:) = LAKE_SFC_C(I,J,:) + INFLUX_C / LAKE_AREA
                     h0 = lake_sfc_bot(i,j) + (lake_wl(i,j,1)+lake_ws(i,j,1))/DENS_H2O &
                                           -lake_depth_sill(i,j)
                     qt = lake_area * h0 * DENS_H2O
                     ! qt is mass of water stored transiently above sill
                     ! now reduce it to amount that discharges this time step
                     if (qt.gt.0.) then
                         IF (large_dyn_small_stat) THEN
                             if (is_watch_cell()) write(*,*) 'qt[1]/A', qt/lake_area
                             if (lake_width_sill(i,j) .gt. 0.) then
                                 t_scale = lake_whole_area(i,j)/(0.9*lake_width_sill(i,j)*sqrt(h0))
                                 qt = qt * (1. - (1.+River%dt_slow/t_scale)**(-2) )
                                 if (.not.use_lake_area_bug) qt = qt * lake_whole_area(i,j)/lake_area
                               endif
                             if (is_watch_cell()) write(*,*) 'qt[2]/A', qt/lake_area
                             qt = min(qt, lake_outflow_frac_ceiling * lake_area &
                                          * max(0.,(lake_wl(i,j,1)+lake_ws(i,j,1))))
                             if (is_watch_cell()) write(*,*) 'qt[3]/A', qt/lake_area
                             qt = min(qt, (lake_wl(i,j,1)+lake_ws(i,j,1)-lake_sfc_w_min)*lake_area )
                             if (is_watch_cell()) write(*,*) 'qt[4]/A', qt/lake_area
                           ELSE
                             t_scale = lake_whole_area(i,j)/(0.9*lake_width_sill(i,j)*sqrt(h0))
                             qt = qt * (1. - (1.+River%dt_slow/t_scale)**(-2) )
                             if (.not.use_lake_area_bug) qt = qt * lake_whole_area(i,j)/lake_area
                             qt = min(qt, lake_outflow_frac_ceiling * lake_area &
                                          * max(0.,(lake_wl(i,j,1)+lake_ws(i,j,1))))
                             qt = min(qt, (lake_wl(i,j,1)+lake_ws(i,j,1)-lake_sfc_w_min)*lake_area )
                           ENDIF
                         if (ice_frac_from_sfc) then
                             out_frac = lake_wl(i,j,1)/(lake_wl(i,j,1)+lake_ws(i,j,1))
                           else
                             out_frac = max (sum(lake_wl(i,j,:))/sum(lake_wl(i,j,:)+lake_ws(i,j,:)), &
                                           lake_wl(i,j,1)/(lake_wl(i,j,1)+lake_ws(i,j,1)))
                           endif
                         out_frac = min(1., max(0., out_frac))
                         if (ice_frac_factor.lt.1.) then
                             ql = (1.-ice_frac_factor*(1.-out_frac)) * qt
                             ql = min (ql, lake_area*sum(lake_wl(i,j,:)))
                           else
                             ql = out_frac * qt
                           endif
                         qs = qt - ql
                         liq_to_flow = ql
                         ice_to_flow = qs
                         qh = 0.
                         if (is_watch_cell()) &
                              write(*,*) 'ql/A,qs/A,A',ql/lake_area,qs/lake_area,lake_area
                         do lev = 1, num_lake_lev
                           if (is_watch_cell() .and. lev.le.10) &
                                write(*,'(a,i3,99(x,a,g23.16))')'l=',lev,&
                                    'wl(1)=',lake_wl(i,j,1),'ws(1)=',lake_ws(i,j,1), &
                                    'wl(l)=',lake_wl(i,j,lev),'ws(l)=',lake_ws(i,j,lev)
                           liq_this_lev = max(0.,min(liq_to_flow, lake_area*lake_wl(i,j,lev)))
                           ice_this_lev = max(0.,min(ice_to_flow, lake_area*lake_ws(i,j,lev)))
                           lake_wl(i,j,lev) = lake_wl(i,j,lev) - liq_this_lev/lake_area
                           lake_ws(i,j,lev) = lake_ws(i,j,lev) - ice_this_lev/lake_area
                           liq_to_flow = liq_to_flow - liq_this_lev
                           ice_to_flow = ice_to_flow - ice_this_lev
                           qh = qh + (clw*liq_this_lev+csw*ice_this_lev)*(lake_T(i,j,lev)-tfreeze)
                           if (lev.gt.1) then
                             ! replensih (liquid) water lost from depth, using ice from surface,
                             ! so as to preserve thickness of deeper layer
                             h = (clw*lake_wl(i,j,lev)+csw*lake_ws(i,j,lev)) &
                                                            *(lake_T(i,j,lev)-tfreeze)
                             lake_ws(i,j,lev) = lake_ws(i,j,lev) + liq_this_lev/lake_area
                             lake_ws(i,j,1)   = lake_ws(i,j,1)   - liq_this_lev/lake_area
                             lake_T (i,j,lev) = tfreeze + &
                                (h +(liq_this_lev/lake_area)*csw*(lake_T(i,j,1)-tfreeze))  &
                                            /(clw*lake_wl(i,j,lev)+csw*lake_ws(i,j,lev))
                           endif
                           if (is_watch_cell() .and. lev.le.10) &
                                write(*,'(a,i3,99(x,a,g23.16))')'l=',lev,&
                                    'wl(1)=',lake_wl(i,j,1),'ws(1)=',lake_ws(i,j,1), &
                                    'wl(l)=',lake_wl(i,j,lev),'ws(l)=',lake_ws(i,j,lev)
                           if (liq_to_flow.eq.0..and.ice_to_flow.eq.0.) exit
                           enddo
                         River%lake_outflow  (i,j)   = qt
                         River%lake_outflow_c(i,j,1) = qs
                         River%lake_outflow_c(i,j,2) = qh
                       endif
                     if (is_watch_cell()) then
                          write(*,*) 'lake_wl(1):', lake_wl(i,j,1)
                          write(*,*) 'lake_ws(1):', lake_ws(i,j,1)
                          write(*,*) 'lake_T (1):', lake_T (i,j,1)
                     endif
                   else
                     River%lake_outflow  (i,j  ) = influx
                     River%lake_outflow_c(i,j,1) = influx_c(1)
                     River%lake_outflow_c(i,j,2) = influx_c(2)
                   endif
              endif

            ! NEXT COMPUTE RIVER-REACH MASS BALANCE (FROM LAKE_OUTFLOW TO OUTFLOW)
          
            if (River%tocell(i,j).gt.0 .or. River%landfrac(i,j).lt.1.) then
                ! avail is volume to be split between outflow and new storage
                avail = River%storage(i,j) + River%lake_outflow(i,j) / DENS_H2O
                ! determine total water storage at end of step
                if (River%reach_length(i,j) .gt. 0.) then
                    if (algor.eq.'linear') then   ! assume outflow = Q0+dQ_dV*dS
                        if (River%storage(i,j) .le. 0.) then
                            Q0 = 0.; dQ_dV = 0.
                          else
                            Q0=River%o_coef(i,j)*River%storage(i,j)**River%o_exp
                            dQ_dV=River%o_exp*Q0/River%storage(i,j)
                          endif
                        if (.not.river_impedes_lake.or..not.lockstep) then
                            River%storage(i,j) = River%storage(i,j) + River%dt_slow *   &
                             (River%lake_outflow(i,j)/(DENS_H2O*River%dt_slow)-Q0) &
                             /(1.+River%dt_slow*dQ_dV)
                        else
                            if (River%storage(i,j) .le. 0.) then
                                dh_dQ = 0.
                            else
                                dh_dQ = River%d_coef(i,j)*River%d_exp*Q0**(River%d_exp-1)
                            endif
                            River%storage(i,j) = River%storage(i,j) + River%dt_slow *   &
                             (River%lake_outflow(i,j)/(DENS_H2O*River%dt_slow)-Q0) &
                             /(1.+dQ_dV*(River%dt_slow+lake_whole_area(i,j)*dh_dQ))
                        endif
                      else if (algor.eq.'nonlin') then   ! assume all inflow at start of step 
                        if (avail .gt. 0.) then
                            River%storage(i,j) = (avail**(1.-River%o_exp) &
                                 + River%o_coef(i,j)*(River%o_exp-1.)*River%dt_slow) &
                                 **(1./(1.-River%o_exp))
                          else
                            River%storage(i,j) = avail
                          endif
                      endif
                  endif
                ! determine total water outflow during step
                River%outflow(i,j) = (avail - River%storage(i,j)) / River%dt_slow
                ! given outflow, determine flow width, depth, velocity
                if (River%outflow(i,j) .le. 0.) then
                    River%depth(i,j) = 0.
                    River%width(i,j) = 0.
                    River%vel(i,j)   = 0.
                  else
                    River%depth(i,j) = River%d_coef(i,j) &
                         * River%outflow(i,j)**River%d_exp
                    River%width(i,j) = River%w_coef(i,j) &
                         * River%outflow(i,j)**River%w_exp
                    River%vel(i,j) = River%outflow(i,j) /                   &
                                        (River%width(i,j) * River%depth(i,j))
                  endif
                ! given water outflow and storage, split other tracked stuff same way
                out_frac = 0.
                if (avail .gt. 0.) out_frac = River%outflow(i,j)/avail
                River%outflow_c(i,j,:) = out_frac * (River%storage_c(i,j,:) &
                                         +River%lake_outflow_c(i,j,:)/DENS_H2O)
                ! 2011/05/13 PCM: fix ice outflow temperature bug
                if (prohibit_cold_ice_outflow) then
                  River%outflow_c(i,j,:) = max(River%outflow_c(i,j,:), 0.)
                else
                  River%outflow_c(i,j,1) = max(River%outflow_c(i,j,1), 0.)
                  if(River%num_phys+1 <= River%num_species) then
                     River%outflow_c(i,j,River%num_phys+1:River%num_species) = &
                       max(River%outflow_c(i,j,River%num_phys+1:River%num_species), 0.)
                  endif
                endif
                River%outflow_c(i,j,1) = min(River%outflow_c(i,j,1), River%outflow(i,j))
                River%storage_c(i,j,:) = River%storage_c(i,j,:)       &
                      + River%lake_outflow_c(i,j,:)/DENS_H2O       &
                      - River%outflow_c(i,j,:)*River%dt_slow
                ! define intensive variables for diagnostics and for use in transformations.
                ! along the way, melt swept snow as necessary. freeze will be a separate
                ! process, added later; it will be different in that frozen river water will
                ! be stationary, thus a different species

                if (River%storage(i,j) .gt. storage_threshold_for_melt) then
                    conc(1) = River%storage_c(i,j,1)/River%storage(i,j)
                    conc(2) = tfreeze + River%storage_c(i,j,2) /  &
                       ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
                    if (River%storage_c(i,j,1).gt.0. .and. conc(2).gt.tfreeze) then
!                    if (River%storage_c(i,j,1).gt.0. .and. River%storage_c(i,j,2).gt.0.) then
                        qmelt = min(hlf*River%storage_c(i,j,1), River%storage_c(i,j,2))
                        River%melt(i,j) = qmelt
                        River%storage_c(i,j,1) = River%storage_c(i,j,1) - qmelt/hlf
                        River%storage_c(i,j,2) = River%storage_c(i,j,2) - qmelt
!                        conc(2) = tfreeze + River%storage_c(i,j,2) /  &
!                           ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
                      endif
                  endif

                if (River%storage(i,j) .gt. storage_threshold_for_diag) then
                    conc(1) = River%storage_c(i,j,1)/River%storage(i,j)
                    conc(2) = tfreeze + River%storage_c(i,j,2) /  &
                       ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
                    diag_mask(i,j) = .TRUE. 
                  else
                    conc(1) = missing
                    conc(2) = missing
                  endif

                ice(i,j)=conc(1)
                temperature(i,j)=conc(2)

                if (River%do_age) then
                    River%removal_c(i,j,River%num_phys+1) = -River%storage(i,j)/sec_in_day
                    River%storage_c(i,j,River%num_phys+1) = River%storage_c(i,j,River%num_phys+1) &
                       - River%removal_c(i,j,River%num_phys+1)*River%dt_slow
                  endif

                if (River%storage(i,j) .gt. 0.) then
                    conc(River%num_phys+1:River%num_species) = &
                       River%storage_c(i,j,River%num_phys+1:River%num_species)/River%storage(i,j)
                  else
                    conc(River%num_phys+1:River%num_species) = 0.
                  endif

                if(River%num_c.gt.0) then
                    if (River%depth(i,j).gt.0. .and. conc(2).gt.100.) then
                        v_r_d = River%vf_ref * River%Q10**((conc(2)-River%t_ref)/10.)&
                           / ((1+River%kinv*conc(River%num_species-River%num_c+1:River%num_species)) &
                           *River%depth(i,j))
                        ! next should not be necessary if storage_c is positive, but maybe it's not.
                        v_r_d = River%vf_ref * River%Q10**((conc(2)-River%t_ref)/10.)&
                           / ((1+River%kinv*max(0.,conc(River%num_species-River%num_c+1:River%num_species)))*River%depth(i,j))
                      else
                        v_r_d = 0.
                      endif
                    River%removal_c(i,j,River%num_species-River%num_c+1:River%num_species) = &
                       River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) &
                       * (1-exp( -v_r_d * River%dt_slow)) &
                       / River%dt_slow
                    River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) = &
                       River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) &
                       - River%removal_c(i,j,River%num_species-River%num_c+1:River%num_species)* River%dt_slow
                  endif
              endif

            ! FINALLY, REDEFINE OUTFLOW AS DISCHARGE IF WE HAVE OCEAN HERE
          
            if (River%landfrac(i,j).lt.1.) then
                River%disw2o(i,j) = River%outflow(i,j)
                River%outflow(i,j) = 0.
                do i_species = 1, num_species
                  River%disc2o(i,j,i_species) = River%outflow_c(i,j,i_species)
                  River%outflow_c(i,j,i_species) = 0.
                  enddo
              endif

          endif
        enddo
      enddo

    

    if (cur_travel .gt. 0) call do_halo_update(River, halo_update(cur_travel))
    
    ! ---- diagnostic section
    if (id_ice > 0) used = send_data (id_ice, &
         ice(isc:iec,jsc:jec), River%Time, mask=diag_mask)
    if (id_temp > 0) used = send_data (id_temp, &
         temperature(isc:iec,jsc:jec), River%Time, mask=diag_mask)
    ! for compatibility with old diag table
    if (id_ice_old > 0) used = send_data (id_ice_old, &
         ice(isc:iec,jsc:jec), River%Time, mask=diag_mask)
    if (id_temp_old > 0) used = send_data (id_temp_old, &
         temperature(isc:iec,jsc:jec), River%Time, mask=diag_mask)

  end subroutine river_physics_step

!#####################################################################
  subroutine setup_halo_update(River, domain)
    type(river_type),      intent(inout) :: River
    type(domain2d),        intent(inout) :: domain

    integer, parameter                   :: MAXCOMM      = 8  ! should be no larger than 8.
    integer                              :: travelnow, nsend2, p, toc
    integer                              :: spos, rpos, n, m, l
    integer                              :: buffer_pos, pos, msgsize
    integer                              :: i, j, i1, j1, i2, j2, i3, j3, i4, j4, k, kk
    integer                              :: send_size, recv_size, siz, i_dest, j_dest
    logical                              :: is_my_recv, is_my_send
    integer                              :: my_recv_index, my_send_index, roff, soff, pe, total_size
    integer                              :: nsend, nrecv, total_send, total_recv, max_send, max_recv
    integer, allocatable, dimension(:,:) :: tocell
    integer, allocatable, dimension(:,:) :: is_recv, ie_recv, js_recv, je_recv
    integer, allocatable, dimension(:,:) :: is1_send, ie1_send, js1_send, je1_send
    integer, allocatable, dimension(:,:) :: is2_send, ie2_send, js2_send, je2_send
    integer, allocatable, dimension(:,:) :: rot_send, rot_recv, dir_send, dir_recv
    integer, allocatable, dimension(:)   :: send_pelist, recv_pelist, pelist_r, pelist_s
    integer, allocatable, dimension(:)   :: send_count, recv_count, recv_size2
    integer, allocatable, dimension(:)   :: isl, iel, jsl, jel  
    integer, allocatable, dimension(:)   :: sbuf, rbuf
    type(comm_type), pointer             :: send => NULL()
    integer, allocatable, dimension(:,:,:) :: i_send, j_send, t_send, p_send, n_send

    call mpp_get_data_domain   (domain, isd, ied, jsd, jed)
    
    !--- first get the travel and tocell information onto data domain
    allocate(tocell(isd:ied,jsd:jed))
    allocate(isl(0:npes-1), iel(0:npes-1), jsl(0:npes-1), jel(0:npes-1) )
    tocell(isc:iec,jsc:jec) = River%tocell(isc:iec,jsc:jec)
    call mpp_update_domains(tocell, domain)
    call mpp_get_compute_domains(domain, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)

    !--- first get the halo update information for send and recv.
    call mpp_get_update_size(domain, nsend, nrecv)

    if(nsend>0) then
       allocate(send_count(nsend))
       do p = 1, nsend
          send_count(p) = mpp_get_num_overlap(domain, EVENT_SEND, p)
       enddo
       if(ANY(send_count .LE. 0)) call mpp_error(FATAL, &
                  "river_mod: send_count should be positive for any entry")
       total_send = sum(send_count)
       allocate(rbuf(4*total_send))
    endif

    if(nrecv>0) then
       allocate(recv_count(nrecv))
       do p = 1, nrecv
          recv_count(p) = mpp_get_num_overlap(domain, EVENT_RECV, p)
       enddo
       if(ANY(recv_count .LE. 0)) call mpp_error(FATAL, &
            "river_mod: recv_count should be positive for any entry")
       total_recv = sum(recv_count)
       allocate(sbuf(4*total_recv))
    endif
 
    !--- pre-post recv
    rpos = 0
    if(nsend > 0) then
       allocate(pelist_s(nsend))
       max_send = maxval(send_count)
       allocate(is1_send(nsend,max_send), ie1_send(nsend,max_send) )
       allocate(js1_send(nsend,max_send), je1_send(nsend,max_send) )
       allocate(is2_send(nsend,max_send), ie2_send(nsend,max_send) )
       allocate(js2_send(nsend,max_send), je2_send(nsend,max_send) )
       allocate(dir_send(nsend,max_send), rot_send(nsend,max_send) )
       call mpp_get_update_pelist(domain, EVENT_SEND, pelist_s)
       do p = 1, nsend    
          call mpp_get_overlap(domain, EVENT_SEND, p, is1_send(p,1:send_count(p)), ie1_send(p,1:send_count(p)), &
               js1_send(p,1:send_count(p)), je1_send(p,1:send_count(p)), dir_send(p,1:send_count(p)), &
               rot_send(p,1:send_count(p)) )
          call mpp_recv(rbuf(rpos+1), glen=4*send_count(p), from_pe=pelist_s(p), block=.FALSE., tag=COMM_TAG_1)
          rpos = rpos + 4*send_count(p)
       enddo
    endif

    spos = 0
    if(nrecv>0) then
       allocate(pelist_r(nrecv))
       max_recv = maxval(recv_count)
       allocate(is_recv (nrecv,max_recv), ie_recv (nrecv,max_recv) )
       allocate(js_recv (nrecv,max_recv), je_recv (nrecv,max_recv) )
       allocate(rot_recv(nrecv,max_recv), dir_recv(nrecv,max_recv) )
       call mpp_get_update_pelist(domain, EVENT_RECV, pelist_r)
       do p = 1, nrecv
          call mpp_get_overlap(domain, EVENT_RECV, p, is_recv(p,1:recv_count(p)), ie_recv(p,1:recv_count(p)), &
               js_recv(p,1:recv_count(p)), je_recv(p,1:recv_count(p)), dir_recv(p,1:recv_count(p)), &
               rot_recv(p,1:recv_count(p)))
          !--- send the information to the process that send data.
          do n = 1, recv_count(p)
             sbuf(spos+(n-1)*4+1) = is_recv(p,n)
             sbuf(spos+(n-1)*4+2) = ie_recv(p,n)
             sbuf(spos+(n-1)*4+3) = js_recv(p,n)
             sbuf(spos+(n-1)*4+4) = je_recv(p,n)
          end do
          call mpp_send(sbuf(spos+1), plen = 4*recv_count(p), to_pe = pelist_r(p), tag=COMM_TAG_1)
          spos = spos + 4*recv_count(p)
       end do
    endif

    call mpp_sync_self(check=EVENT_RECV)
    !--- unpack
    do p = nsend, 1, -1
       rpos = rpos - 4*send_count(p)
       do n = 1, send_count(p)
          is2_send(p,n) = rbuf(rpos+(n-1)*4+1)
          ie2_send(p,n) = rbuf(rpos+(n-1)*4+2)
          js2_send(p,n) = rbuf(rpos+(n-1)*4+3)
          je2_send(p,n) = rbuf(rpos+(n-1)*4+4)
       end do
    end do
        
    call mpp_sync_self()

    is_my_recv = .false.
    do p = 1, nsend
       if(pelist_s(p) == mpp_pe()) then
          is_my_recv = .true.
          my_recv_index = p
       endif
    enddo
    is_my_send = .false.
    do p = 1, nrecv
       if(pelist_r(p) == mpp_pe()) then
          is_my_send = .true.
          my_send_index = p
       endif
    enddo
    roff = 0
    soff = 0
    if( is_my_recv) then
       nrecv_update = nsend
       if(nsend > 0) then
          allocate(recv_pelist(nrecv_update))
          recv_pelist = pelist_s
       endif
    else
       nrecv_update = nsend + 1
       allocate(recv_pelist(nrecv_update))
       my_recv_index = 1
       roff = 1
       recv_pelist(1) = mpp_pe()
       do p = 1, nsend
          recv_pelist(p+1) = pelist_s(p)
       enddo
    endif
    if( is_my_send ) then
       nsend_update = nrecv
       if(nrecv>0) then
          allocate(send_pelist(nsend_update))
          send_pelist = pelist_r
       endif
    else
       nsend_update = nrecv + 1
       allocate(send_pelist(nsend_update))
       my_send_index = 1
       soff = 1
       send_pelist(1) = mpp_pe()
       do p = 1, nrecv
          send_pelist(p+1) = pelist_r(p)
       enddo
    endif

    allocate(halo_update(maxtravel) )
    if(nsend_update>0) then
       do travelnow = 1, maxtravel
          allocate(halo_update(travelnow)%send(nsend_update))
          halo_update(travelnow)%send(:)%count = 0
          do p = 1, nsend_update
             halo_update(travelnow)%send(p)%count = 0
             halo_update(travelnow)%send(p)%pe    = send_pelist(p)
          end do
       end do
    endif

    if(nrecv_update>0) then
       do travelnow = 1, maxtravel
          allocate(halo_update(travelnow)%recv(nrecv_update))
          halo_update(travelnow)%recv(:)%count = 0
          do p = 1, nrecv_update
             halo_update(travelnow)%recv(p)%count = 0
             halo_update(travelnow)%recv(p)%pe    = recv_pelist(p)
          end do
       end do
    endif

    do p = 1, nsend
       pe = pelist_s(p) - mpp_root_pe()
       !--- configure points need to receive from other pe.
       !--- (i,j) --- halo index on the pe sent to, one neighbor pe data domain
       !--- (i1,j1) --- neighbor index of (i,j) on the pe sent to, on neighbor pe compute domain
       !--- (i2,j2) --- my index corresponding to (i,j), on my compute domain
       !--- (i3,j3) --- neighbor index of (i2,j2), on my data domain
       !--- (i4,j4) --- index of (i1,j1) tocell. 
       do n = 1, send_count(p)
          select case ( dir_send(p,n) )
          case(1)  ! east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             i1 = i - 1
             do j = js2_send(p,n), je2_send(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsl(pe) .OR. j1 > jel(pe) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! w->e
                      i2 = is1_send(p,n)
                      i3 = i2 -1
                      j2 = js1_send(p,n) + j  - js2_send(p,n)
                      j3 = js1_send(p,n) + j1 - js2_send(p,n)
                   case (NINETY) ! s->e
                      i2 = is1_send(p,n) + (je2_send(p,n) - j )
                      i3 = is1_send(p,n) + (je2_send(p,n) - j1)
                      j2 = js1_send(p,n)
                      j3 = j2 - 1
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                      end if
                   end if
                end do
             end do
          case(2)  ! south east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i - 1
             j1 = j + 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! nw->se 
                i3 = i2 - 1
                j3 = j2 + 1
             case (NINETY)
                i3 = i2 - 1
                j3 = j2 - 1
             case (MINUS_NINETY)
                i3 = i2 + 1
                j3 = j2 + 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                end if
             end if
          case(3)  ! south
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             j1 = j + 1
             do i = is2_send(p,n), ie2_send(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isl(pe) .OR. i1 > iel(pe) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! n->s
                      i2 = is1_send(p,n) + i  - is2_send(p,n)
                      i3 = is1_send(p,n) + i1 - is2_send(p,n)
                      j2 = js1_send(p,n)
                      j3 = j2 + 1
                   case (MINUS_NINETY) ! e->s
                      i2 = is1_send(p,n)
                      i3 = i2 + 1
                      j2 = js1_send(p,n) + (ie2_send(p,n) - i )
                      j3 = js1_send(p,n) + (ie2_send(p,n) - i1)
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                      end if
                   end if
                end do
             end do
          case(4)  ! south west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i + 1
             j1 = j + 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! ne->sw 
                i3 = i2 + 1
                j3 = j2 + 1
             case (NINETY) !  
                i3 = i2 - 1
                j3 = j2 + 1
             case (MINUS_NINETY)
                i3 = i2 + 1
                j3 = j2 - 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                end if
             end if
          case(5)  ! west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             i1 = i + 1
             do j = js2_send(p,n), je2_send(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsl(pe) .OR. j1 > jel(pe) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! e->w
                      i2 = is1_send(p,n)
                      i3 = i2 + 1
                      j2 = js1_send(p,n) + j  - js2_send(p,n)
                      j3 = js1_send(p,n) + j1 - js2_send(p,n)
                   case (NINETY) ! n->w
                      i2 = is1_send(p,n) + (je2_send(p,n) - j )
                      i3 = is1_send(p,n) + (je2_send(p,n) - j1)
                      j2 = js1_send(p,n)
                      j3 = j2 + 1
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                      end if
                   end if
                end do
             end do
          case(6)  ! north west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i + 1
             j1 = j - 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! se->nw 
                i3 = i2 + 1
                j3 = j2 - 1
             case (NINETY) !  
                i3 = i2 + 1
                j3 = j2 + 1
             case (MINUS_NINETY)
                i3 = i2 - 1
                j3 = j2 - 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                end if
             end if
          case(7)  ! north
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             j1 = j - 1
             do i = is2_send(p,n), ie2_send(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isl(pe) .OR. i1 > iel(pe)) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! s->n
                      i2 = is1_send(p,n) + i  - is2_send(p,n)
                      i3 = is1_send(p,n) + i1 - is2_send(p,n)
                      j2 = js1_send(p,n)
                      j3 = j2 - 1
                   case (MINUS_NINETY) ! w->n
                      i2 = is1_send(p,n)
                      i3 = i2 - 1
                      j2 = js1_send(p,n) + (ie2_send(p,n) - i )
                      j3 = js1_send(p,n) + (ie2_send(p,n) - i1)
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                      end if
                   end if
                end do
             end do
          case(8)  ! north east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i - 1
             j1 = j - 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! sw->ne 
                i3 = i2 - 1
                j3 = j2 - 1
             case (NINETY) !  
                i3 = i2 + 1
                j3 = j2 - 1
             case (MINUS_NINETY)
                i3 = i2 - 1
                j3 = j2 + 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p+roff), i2, j2)
                end if
             end if
          end select
       end do       
    enddo

    do p = 1, nrecv
       !--- configure points need to send to other pe.
       !--- (i,j) ---  index on my data domain
       !--- (i1,j1) --- index on my compute domain corresponding to (i,j)
       !--- (i2,j2) --- index of (i1,j1) tocell
       do n = 1, recv_count(p)
          select case ( dir_recv(p,n) )
          case(1)  ! east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             i1 = i - 1
             do j = js_recv(p,n), je_recv(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsc .OR. j1 > jec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                   end if
                end do
             end do
          case(2)  ! south east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i - 1
             j1 = j + 1
             if(River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                end if
             end if
          case(3)  ! south
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             j1 = j + 1
             do i = is_recv(p,n), ie_recv(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isc .OR. i1 > iec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                   end if
                end do
             end do
          case(4)  ! south west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i + 1
             j1 = j + 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                end if
             end if
          case(5)  ! west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             i1 = i + 1
             do j = js_recv(p,n), je_recv(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsc .OR. j1 > jec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                   end if
                end do
             end do
          case(6)  ! north west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i + 1
             j1 = j - 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                end if
             end if
          case(7)  ! north
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             j1 = j - 1
             do i = is_recv(p,n), ie_recv(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isc .OR. i1 > iec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                   end if
                end do
             end do
          case(8)  ! north east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i - 1
             j1 = j - 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p+soff), i1, j1)
                end if
             end if
          end select
       end do       
    end do

    allocate(in_domain(isc:iec,jsc:jec))
    in_domain = .true.
    do j = jsc, jec
       do i = isc, iec
          if(River%tocell(i,j) > 0) then
             i_dest = i + di(River%tocell(i,j))
             j_dest = j + dj(River%tocell(i,j))
             if(i_dest < isc .OR. i_dest > iec .OR. j_dest < jsc .OR. j_dest > jec) then 
                LOOP_TRAVEL: do travelnow = 1, maxtravel
                   do p = 1, nrecv
                      send => halo_update(travelnow)%send(p+soff)
                      do n = 1, send%count
                         if(send%i(n) == i .AND. send%j(n) == j) then
                            in_domain(i,j) = .false.
                            exit LOOP_TRAVEL 
                         end if
                      end do
                   end do
                end do LOOP_TRAVEL
                if(in_domain(i,j)) then
                   i_dest = i
                   j_dest = j
                end if
             end if
             River%i_tocell(i,j) = i_dest
             River%j_tocell(i,j) = j_dest
          end if
       end do
    end do

    !--- add points that sent to self.
    do j = jsc, jec
       do i = isc, iec
          m = River%travel(i,j)
          if(m >0 .and. in_domain(i,j) ) then
             call add_single_overlap(halo_update(m)%send(my_send_index), i, j)
             call add_single_overlap(halo_update(m)%recv(my_recv_index), River%i_tocell(i, j), River%j_tocell(i, j))
          end if
       end do
    end do

    !--- the following is for the purpose of bitwise reproduce between processor count
    if(nrecv_update>0) allocate(recv_size2(nrecv_update))
    do p=1, nrecv_update
       call mpp_recv(recv_size2(p), glen = 1, from_pe = recv_pelist(p), block=.FALSE., tag=COMM_TAG_2 )
    enddo

    do p= 1, nsend_update
       msgsize = 0
       do m = 1, maxtravel
          msgsize = msgsize + 2*halo_update(m)%send(p)%count
       enddo
       call mpp_send(msgsize, plen = 1, to_pe = send_pelist(p), tag=COMM_TAG_2)
    enddo

    call mpp_sync_self(check=EVENT_RECV)
    do p=1, nrecv_update 
       recv_size = 0   
       do m = 1, maxtravel
          recv_size = recv_size + halo_update(m)%recv(p)%count
       end do
       recv_size = recv_size*2
       if(recv_size2(p) .NE. recv_size) then
          print*, "At pe = ", mpp_pe()," p = ", p, " from_pe = ", recv_pelist(p), ", send_size = ", recv_size2(p), "recv_size = ", recv_size
          call mpp_error(FATAL, "river_physics_mod: mismatch at send size and recv size")
       endif
    enddo
    call mpp_sync_self()

    total_size = 0
    do p = 1, nrecv_update
       total_size = total_size + recv_size2(p)
    enddo

    if(total_size >0) allocate(recv_buffer(total_size))
    pos = 0
    do p=1, nrecv_update
       if(recv_size2(p) >0) then
          call mpp_recv(recv_buffer(pos+1), glen = recv_size2(p), from_pe = recv_pelist(p), block=.FALSE., tag=COMM_TAG_3 )
          pos = pos + recv_size2(p)
       endif
    enddo

    send_size = 0
    do p = 1, nsend_update
       do m = 1, maxtravel
          send_size = send_size + halo_update(m)%send(p)%count
       end do
    end do
    send_size = send_size*2
    if(send_size>0) allocate(send_buffer(send_size))
    pos = 0
    do p= 1, nsend_update
       buffer_pos = pos
       do m = 1, maxtravel
          do n = 1, halo_update(m)%send(p)%count
             send_buffer(pos+1) = halo_update(m)%send(p)%i(n)
             send_buffer(pos+2) = halo_update(m)%send(p)%j(n)
             pos = pos + 2
          end do
       end do
       msgsize = pos - buffer_pos
       if(msgsize >0) then
          call mpp_send(send_buffer(buffer_pos+1), plen = msgsize, to_pe = send_pelist(p), tag=COMM_TAG_3)
       end if      
    end do

    call mpp_sync_self(check=EVENT_RECV)

    !--- unpack buffer
    allocate(i_send(isc:iec,jsc:jec,8), j_send(isc:iec,jsc:jec,8) )
    allocate(p_send(isc:iec,jsc:jec,8), t_send(isc:iec,jsc:jec,8) )
    allocate(n_send(isc:iec,jsc:jec,8), nlev(isc:iec,jsc:jec))
    nlev = 0
    pos = 0
    do p=1, nrecv_update
       if(recv_size2(p) >0) then
          do m = 1, maxtravel
             do n = 1, halo_update(m)%recv(p)%count
                i = halo_update(m)%recv(p)%i(n)
                j = halo_update(m)%recv(p)%j(n)
                i1 = recv_buffer(pos+1)
                j1 = recv_buffer(pos+2)
                pos = pos + 2
                do k = 1, nlev(i,j)
                   if( j1 < j_send(i,j,k) .OR. (j1 == j_send(i,j,k) .AND. i1 < i_send(i,j,k) ) ) then
                      do kk = nlev(i,j)+1, k+1, -1
                         i_send(i,j,kk) = i_send(i,j,kk-1)
                         j_send(i,j,kk) = j_send(i,j,kk-1)
                         p_send(i,j,kk) = p_send(i,j,kk-1)
                         t_send(i,j,kk) = t_send(i,j,kk-1)
                         n_send(i,j,kk) = n_send(i,j,kk-1)
                         halo_update(t_send(i,j,kk))%recv(p_send(i,j,kk))%k(n_send(i,j,kk)) = &
                             halo_update(t_send(i,j,kk))%recv(p_send(i,j,kk))%k(n_send(i,j,kk)) + 1
                      end do
                      exit
                   end if
                end do
                nlev(i,j) = nlev(i,j) + 1
                i_send(i,j,k) = i1
                j_send(i,j,k) = j1
                p_send(i,j,k) = p
                t_send(i,j,k) = m
                n_send(i,j,k) = n         
                halo_update(m)%recv(p)%k(n) = k               
             end do
          end do
       end if
    end do

    call mpp_sync_self()    
    if(allocated(send_buffer)) deallocate(send_buffer)
    if(allocated(recv_buffer)) deallocate(recv_buffer)
    if(allocated(recv_size2 )) deallocate(recv_size2 )

    !--- set up buffer for send and recv.
    send_size = 0
    do m = 1, maxtravel
       siz = 0
       do p = 1, nsend_update
          siz = siz + halo_update(m)%send(p)%count
       end do
       send_size = max(send_size, siz)
    end do
    send_size = send_size*(num_species+1)
    if(send_size > 0) allocate(send_buffer(send_size))

    recv_size = 0
    do m = 1, maxtravel
       siz = 0
       do p = 1, nrecv_update
          siz = siz + halo_update(m)%recv(p)%count
       end do
       recv_size = max(recv_size, siz)
    end do
    recv_size = recv_size*(num_species+1)
    if(recv_size > 0) allocate(recv_buffer(recv_size))

    deallocate(tocell)
    deallocate(isl, iel, jsl, jel )
    deallocate(sbuf, rbuf)
    deallocate(is_recv, ie_recv, js_recv, je_recv)
    deallocate(is1_send, ie1_send, js1_send, je1_send)
    deallocate(is2_send, ie2_send, js2_send, je2_send)
    deallocate(rot_send, rot_recv, send_count, recv_count)
    if(ALLOCATED(pelist_r)) deallocate(pelist_r)
    if(ALLOCATED(pelist_s)) deallocate(pelist_s)
    if(ALLOCATED(send_pelist)) deallocate(send_pelist)
    if(ALLOCATED(recv_pelist)) deallocate(recv_pelist)
    return

  end subroutine setup_halo_update

!###############################################################################
  subroutine do_halo_update(River, update)
     type(river_type),    intent(inout) :: River
     type(halo_update_type), intent(in) :: update
     type(comm_type), pointer           :: send=>NULL()
     type(comm_type), pointer           :: recv=>NULL()
     integer                            :: buffer_pos, pos, recv_buffer_pos
     integer                            :: p, n, i, j, count, l, k
     real                               :: wrk_c(isc:iec,jsc:jec,num_species, 8)
     real                               :: wrk  (isc:iec,jsc:jec, 8)

     !--- pre-post recv data
     pos = 0
     do p = 1, nrecv_update
        recv => update%recv(p)
        count = recv%count
        if(count == 0) cycle
        call mpp_recv(recv_buffer(pos+1), glen=count*(num_species+1), from_pe=recv%pe, block=.FALSE., tag=COMM_TAG_4 )            
        pos = pos + count*(num_species+1)
     enddo
     recv_buffer_pos = pos

     !--- send the data
     pos = 0
     do p = 1, nsend_update
        send => update%send(p)
        count = send%count
        if(count == 0) cycle
        buffer_pos = pos
        do n = 1, count
           i = send%i(n)
           j = send%j(n)
           pos = pos + 1
           send_buffer(pos)   = River%outflow(i,j)
           do l = 1, num_species
              pos = pos + 1
              send_buffer(pos) = River%outflow_c(i,j,l)
           end do
        end do
        call mpp_send(send_buffer(buffer_pos+1), plen=count*(num_species+1), to_pe = send%pe, tag=COMM_TAG_4 ) 
     end do

     call mpp_sync_self(check=EVENT_RECV)

     !--- update the buffer in reverse order
     nlev = 0
     pos = recv_buffer_pos
     do p = nrecv_update, 1, -1
        recv => update%recv(p)
        count = recv%count
        if(count == 0) cycle
        pos = recv_buffer_pos - count*(num_species+1)
        recv_buffer_pos = pos
        do n = 1, count
           i = recv%i(n)
           j = recv%j(n)
           k = recv%k(n)
           pos = pos + 1
           wrk(i,j,k) = recv_buffer(pos)
           nlev(i,j) = nlev(i,j)+1
           do l = 1, num_species
              pos = pos + 1
              wrk_c(i,j,l,k) = recv_buffer(pos)
           end do
        end do
     enddo

     do j = jsc, jec
        do i = isc, iec
           do k = 1, nlev(i,j)
              River%inflow(i,j)   = River%inflow(i,j) + wrk(i,j,k)
              River%inflow_c(i,j,:) = River%inflow_c(i,j,:) + wrk_c(i,j,:,k)
           end do
        end do
     end do

     call mpp_sync_self()


     return

  end subroutine do_halo_update


!###############################################################################
!  This routine will add one point for send/recv into the data type, allocate 
!  memory or expand memory if needed 

  subroutine add_single_overlap(comm, i, j)
    type(comm_type), intent(inout) :: comm
    integer,         intent(in)    :: i, j
    integer, parameter             :: DEFAULT_SIZE = 40 ! too big or too small
    integer                        :: count, maxcount
    integer, allocatable           :: tmp(:)

    count = comm%count 
    count = count + 1
    if(count == 1) then ! assign default space to hold index
       allocate(comm%i(DEFAULT_SIZE))
       allocate(comm%j(DEFAULT_SIZE))
       allocate(comm%k(DEFAULT_SIZE))
    end if
    maxcount = size(comm%i)

    if(count > maxcount) then ! need to expend the size to hold the index
       allocate(tmp(maxcount))
       tmp = comm%i
       deallocate(comm%i)
       allocate(comm%i(2*maxcount))
       comm%i(1:maxcount) = tmp
       tmp = comm%j
       deallocate(comm%j)
       allocate(comm%j(2*maxcount))
       comm%j(1:maxcount) = tmp
       deallocate(tmp)
       deallocate(comm%k)
       allocate(comm%k(2*maxcount))
    end if
    comm%i(count) = i
    comm%j(count) = j    
    comm%k(count) = 0
    comm%count    = count

    return    

  end subroutine add_single_overlap

!#####################################################################

end module river_physics_mod
