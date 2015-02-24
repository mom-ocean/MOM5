module ocean_blob_util_mod
#define COMP isc:iec,jsc:jec
!
!<CONTACT EMAIL="m.bates@student.unsw.edu.au"> Michael L. Bates
!</CONTACT>
!
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module contains subroutines that are common (or are likely to 
! need be common to future implementations) to some or all of the 
! various modules that run the Lagrangian blob scheme.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains subroutines that are common (or may be common
! in a future implementation) to some of the modules that make up the
! blobs framework.
!
! Some of the subroutines contained herein perform tasks such as 
! performing checksums, checking a linked lists for very small blobs 
! (and deleting them), inserting blobs into a list, a bunch of routines 
! for writing blob restart and history files, grid cell search algorithms, 
! buffer manipulations, computations etc.
!
! The module has no namelist.  All potential namelist variables are
! controlled through the ocean_blob_mod namelist.
!</DESCRIPTION>
!
!<INFO>
! <REFERENCE>
! Cormen, T. H, Leiserson, C. E. Rivest, R. L., Stein, C. (2001) Introduction
! to Algorithms.  The MIT Press.
! </REFERENCE>
!
! <REFERENCE>
! Shepard, D. (1968) A two-dimensional interpolation function for 
! irregularly-spaced data.  In: Proceedings of the 1968 23rd ACM national
! conference. ACM '68. ACM, New York, NY, USA, pp. 517-524.
! </REFERENCE>
!
! <REFERENCE>
! Murray, R. J. (1996) Explicit generation of orthogonal grids for
! ocean models.  Journal of Computational Physics 126, 251-273.
! </REFERENCE>
!
!</INFO>


use constants_mod,   only: rad_to_deg, epsln
use fms_mod,         only: error_mesg, FATAL, WARNING, stdout, stderr, mpp_error
use fms_mod,         only: read_data
use mpp_domains_mod, only: mpp_global_sum, mpp_get_neighbor_pe, mpp_update_domains
use mpp_domains_mod, only: mpp_get_current_ntile, mpp_get_tile_id
use mpp_domains_mod, only: NORTH, SOUTH, EAST, WEST
use mpp_mod,         only: mpp_sum, NULL_PE
use grid_mod,        only: get_grid_cell_vertices, get_grid_size

use ocean_parameters_mod, only: onehalf, rho0r, grav, omega_earth
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, DEPTH_BASED
use ocean_types_mod,      only: ocean_thickness_type, ocean_lagrangian_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_time_type
use ocean_types_mod,      only: ocean_blob_type, ocean_grid_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_domain_type, ocean_density_type, blob_grid_type
use ocean_workspace_mod,  only: wrk1, wrk2
use ocean_util_mod,       only: write_chksum_3d, write_chksum_3d_int

implicit none

include 'netcdf.inc'

private

integer :: vert_coordinate_class
integer :: vert_coordinate

! set module types
type(ocean_grid_type),   pointer :: Grd   => NULL()
type(ocean_domain_type), pointer :: Dom   => NULL()
type(ocean_domain_type), pointer :: Bdom  => NULL()
type(blob_grid_type),    pointer :: Info  => NULL()
real, dimension(:,:,:),  pointer :: ht    => NULL()
real, dimension(:,:,:),  pointer :: hu    => NULL()

real, dimension(:,:,:,:), allocatable :: vert_t
real, dimension(:,:,:), allocatable :: ij_im1j, im1jm1_im1j, ij_ijm1, im1jm1_ijm1
real, dimension(:,:,:), allocatable :: t_im1jm1, t_ijm1, t_ij, t_im1j
real, dimension(:,:), allocatable :: xt, xu
real, dimension(:,:), allocatable :: yt, yu
real, dimension(:,:), allocatable :: datdtime_r ! 1./(Grd%dat*dtime)
integer, dimension(9) :: im, jm
integer, dimension(5) :: it, jt
integer, dimension(4) :: iu, ju
integer :: nig,njg,ni,nip1
real :: grav_rho0r
real :: two_omega
real :: dtimer
real :: grav_dtimer

public blob_util_init
public blob_chksum
public lagrangian_system_chksum
public E_and_L_totals
public blob_delete
public insert_blob
public count_blob
public put_att
public inq_var
public get_double
public get_int
public put_double
public put_int
public def_var
public blob_util_end
public write_blobs
public check_ijcell
public check_kcell
public kill_blob
public free_blob_memory
public hashfun
public unlink_blob
public allocate_interaction_memory
public reallocate_interaction_memory
public interp_tcoeff
public interp_ucoeff
public check_cyclic

! variables that are read in during init
integer :: index_temp
integer :: index_salt
integer :: global_sum_flag
integer :: num_prog_tracers
logical :: bitwise_reproduction
logical :: debug_this_module
logical :: really_debug
integer :: isd,ied,jsd,jed
integer :: isc,iec,jsc,jec
integer :: isg,ieg,jsg,jeg
integer :: nk
integer :: isbd, iebd, jsbd, jebd
real    :: blob_small_mass

contains

!######################################################################
! <SUBROUTINE NAME="blob_util_init">
!
! <DESCRIPTION>
! Initialises this module.
! </DESCRIPTION>
!
subroutine blob_util_init(Grid, Domain, PE_info, Blob_domain, &
                          sum_flag, num_tracers, itemp, isalt,   &
                          small_mass, dtime, bitwise, ver_coordinate_class,&
                          ver_coordinate, debug, debug_lots)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_domain_type), intent(in), target :: Blob_domain
  type(blob_grid_type),    intent(in), target :: PE_info
  integer,                 intent(in)         :: sum_flag
  integer,                 intent(in)         :: num_tracers
  integer,                 intent(in)         :: itemp
  integer,                 intent(in)         :: isalt
  real,                    intent(in)         :: small_mass
  real,                    intent(in)         :: dtime
  logical,                 intent(in)         :: bitwise
  integer,                 intent(in)         :: ver_coordinate_class
  integer,                 intent(in)         :: ver_coordinate
  logical,                 intent(in)         :: debug
  logical,                 intent(in)         :: debug_lots

  real, dimension(:,:), allocatable :: dxte, dxue, dytn, dyun, verticies
  real, dimension(:,:), allocatable :: lon_vert, lat_vert
  integer, dimension(:), allocatable :: tile_ids
  integer :: tile, nlon, nlat
  integer :: stdoutunit
  integer :: nfstatus, gsfile, x_vert_t_id, y_vert_t_id, ni, nj
  integer :: i,j,m, iscii,jscjj
  integer :: ii(4), jj(4)
  logical :: do_wst_bnd, do_est_bnd, do_nth_bnd, do_sth_bnd
  logical :: other_vars
  character(len=128) :: filename

  stdoutunit = stdout()

  write (stdoutunit, '(/,a)') 'Note from ocean_blob_util_mod: Initialising ocean_blob_util_mod'

  Grd  => Grid
  Dom  => Domain
  Bdom => Blob_domain
  Info => PE_info

  isd=Dom%isd; ied=Dom%ied; jsd=Dom%jsd; jed=Dom%jed
  isc=Dom%isc; iec=Dom%iec; jsc=Dom%jsc; jec=Dom%jec
  isg=Dom%isg; ieg=Dom%ieg; jsg=Dom%jsg; jeg=Dom%jeg
  ni = Grd%ni; nk = Grd%nk
  nip1 = ni+1

  nig = ieg-(isg-1); njg = jeg-(jsg-1)
  isbd=Bdom%isd; iebd=Bdom%ied; jsbd=Bdom%jsd; jebd=Bdom%jed

  index_temp            = itemp
  index_salt            = isalt
  num_prog_tracers      = num_tracers
  global_sum_flag       = sum_flag
  bitwise_reproduction  = bitwise
  blob_small_mass       = small_mass
  debug_this_module     = debug
  really_debug          = debug_lots
  vert_coordinate_class = ver_coordinate_class
  vert_coordinate       = ver_coordinate

  ! handy constants
  grav_rho0r  = grav*rho0r
  two_omega   = 2 * omega_earth
  dtimer      = 1./dtime
  grav_dtimer = grav*dtimer

  !1: (i,j)
  !2: (i+1,j), 3: (i+1,j+1), 4: (i,j+1), 5: (i-1,j+1)
  !6: (i-1,j), 7: (i-1,j-1), 8: (i,j-1), 9: (i+1,j-1)
  im(1)= 0; jm(1)= 0
  im(2)= 1; jm(2)= 0
  im(3)= 1; jm(3)= 1
  im(4)= 0; jm(4)= 1
  im(5)=-1; jm(5)= 1
  im(6)=-1; jm(6)= 0
  im(7)=-1; jm(7)=-1
  im(8)= 0; jm(8)=-1
  im(9)= 1; jm(9)=-1

  allocate(xt(isbd:iebd, jsbd:jebd))
  allocate(yt(isbd:iebd, jsbd:jebd))
  allocate(xu(isbd:iebd, jsbd:jebd))
  allocate(yu(isbd:iebd, jsbd:jebd))

  xt(isc:iec, jsc:jec) = Grd%xt(isc:iec, jsc:jec)
  yt(isc:iec, jsc:jec) = Grd%yt(isc:iec, jsc:jec)
  xu(isc:iec, jsc:jec) = Grd%xu(isc:iec, jsc:jec)
  yu(isc:iec, jsc:jec) = Grd%yu(isc:iec, jsc:jec)
  call mpp_update_domains(xt(:,:), Bdom%domain2d, complete=.false.)
  call mpp_update_domains(yt(:,:), Bdom%domain2d, complete=.false.)
  call mpp_update_domains(xu(:,:), Bdom%domain2d, complete=.false.)
  call mpp_update_domains(yu(:,:), Bdom%domain2d, complete=.true.)

  ! For cyclic boundary conditions, we need to adjust the halo values
  ! for lon and lat so that the lon and lat are always monotonic
  ! on the one processor to make the calculations for the point location
  ! scheme a bit easier.
  if (Grd%cyclic_x) then 
     do_wst_bnd = .false.
     do_est_bnd = .false.
     do j=jsc,jec
        if (xt(iec+1,j) < xt(iec,  j)) do_est_bnd=.true.
        if (xt(isc,  j) < xt(isc-1,j)) do_wst_bnd=.true.
     enddo

     allocate(dxte(isbd:iebd,jsbd:jebd))
     allocate(dxue(isbd:iebd,jsbd:jebd))
     dxte(isc:iec,jsc:jec) = Grd%dxte(isc:iec,jsc:jec)
     dxue(isc:iec,jsc:jec) = Grd%dxue(isc:iec,jsc:jec)
     call mpp_update_domains(dxte(:,:), Bdom%domain2d)
     call mpp_update_domains(dxue(:,:), Bdom%domain2d)
        
     if (do_wst_bnd) then
        xt(isc-1,jsc:jec) = &
          xt(isc  ,jsc:jec) - rad_to_deg*2*dxte(isc-1,jsc:jec)/(Info%ht(isc  ,jsc:jec,1)+Info%ht(isc-1,jsc:jec,1))
        xt(isc-2,jsc:jec) = &
          xt(isc-1,jsc:jec) - rad_to_deg*2*dxte(isc-2,jsc:jec)/(Info%ht(isc-1,jsc:jec,1)+Info%ht(isc-2,jsc:jec,1))
        xu(isc-1,jsc:jec) = &
          xu(isc  ,jsc:jec) - rad_to_deg*2*dxue(isc-1,jsc:jec)/(Info%hu(isc  ,jsc:jec,2)+Info%hu(isc-1,jsc:jec,2))
        xu(isc-2,jsc:jec) = & 
          xu(isc-1,jsc:jec) - rad_to_deg*2*dxue(isc-2,jsc:jec)/(Info%hu(isc-1,jsc:jec,2)+Info%hu(isc-2,jsc:jec,2))
        
        ! Take care of the corners
        if(Info%pe_NW /= NULL_PE) then
           xt(isc-1,jed:jebd) = &
             xt(isc  ,jed:jebd) - rad_to_deg*2*dxte(isc-1,jed:jebd)/(Info%ht(isc  ,jed:jebd,1)+Info%ht(isc-1,jed:jebd,1))
           xt(isc-2,jed:jebd) = &
             xt(isc-1,jed:jebd) - rad_to_deg*2*dxte(isc-2,jed:jebd)/(Info%ht(isc-1,jed:jebd,1)+Info%ht(isc-2,jed:jebd,1))
           xu(isc-1,jed:jebd) = &
             xu(isc  ,jed:jebd) - rad_to_deg*2*dxue(isc-1,jed:jebd)/(Info%hu(isc  ,jed:jebd,2)+Info%hu(isc-1,jed:jebd,2))
           xu(isc-2,jed:jebd) = &
             xu(isc-1,jed:jebd) - rad_to_deg*2*dxue(isc-2,jed:jebd)/(Info%hu(isc-1,jed:jebd,2)+Info%hu(isc-2,jed:jebd,2))
        endif
        if(Info%pe_SW /= NULL_PE) then
           xt(isc-1,jsbd:jsd) = &
             xt(isc  ,jsbd:jsd) - rad_to_deg*2*dxte(isc-1,jsbd:jsd)/(Info%ht(isc  ,jsbd:jsd,1)+Info%ht(isc-1,jsbd:jsd,1))
           xt(isc-2,jsbd:jsd) = &
             xt(isc-1,jsbd:jsd) - rad_to_deg*2*dxte(isc-2,jsbd:jsd)/(Info%ht(isc-1,jsbd:jsd,1)+Info%ht(isc-2,jsbd:jsd,1))
           xu(isc-1,jsbd:jsd) = &
             xu(isc  ,jsbd:jsd) - rad_to_deg*2*dxue(isc-1,jsbd:jsd)/(Info%hu(isc  ,jsbd:jsd,2)+Info%hu(isc-1,jsbd:jsd,2))
           xu(isc-2,jsbd:jsd) = &
             xu(isc-1,jsbd:jsd) - rad_to_deg*2*dxue(isc-2,jsbd:jsd)/(Info%hu(isc-1,jsbd:jsd,2)+Info%hu(isc-2,jsbd:jsd,2))
        endif
     endif

     if (do_est_bnd) then
        xt(iec+1,jsc:jec) = &
          xt(iec  ,jsc:jec) + rad_to_deg*2*dxte(iec  ,jsc:jec)/(Info%ht(iec  ,jsc:jec,1)+Info%ht(iec+1,jsc:jec,1))
        xt(iec+2,jsc:jec) = &
          xt(iec+1,jsc:jec) + rad_to_deg*2*dxte(iec+1,jsc:jec)/(Info%ht(iec+1,jsc:jec,1)+Info%ht(iec+2,jsc:jec,1))
        xu(iec+1,jsc:jec) = &
          xu(iec  ,jsc:jec) + rad_to_deg*2*dxue(iec  ,jsc:jec)/(Info%hu(iec  ,jsc:jec,2)+Info%hu(iec+1,jsc:jec,2))
        xu(iec+2,jsc:jec) = &
          xu(iec+1,jsc:jec) + rad_to_deg*2*dxue(iec+1,jsc:jec)/(Info%hu(iec+1,jsc:jec,2)+Info%hu(iec+2,jsc:jec,2))

        ! Take care of the corners
        if(Info%pe_NE /= NULL_PE) then
           xt(iec+1,jed:jebd) = &
             xt(iec  ,jed:jebd) + rad_to_deg*2*dxte(iec  ,jed:jebd)/(Info%ht(iec  ,jed:jebd,1)+Info%ht(iec+1,jed:jebd,1))
           xt(iec+2,jed:jebd) = &
             xt(iec+1,jed:jebd) + rad_to_deg*2*dxte(iec+1,jed:jebd)/(Info%ht(iec+1,jed:jebd,1)+Info%ht(iec+2,jed:jebd,1))
           xu(iec+1,jed:jebd) = &
             xu(iec  ,jed:jebd) + rad_to_deg*2*dxue(iec  ,jed:jebd)/(Info%hu(iec  ,jed:jebd,2)+Info%hu(iec+1,jed:jebd,2))
           xu(iec+2,jed:jebd) = &
             xu(iec+1,jed:jebd) + rad_to_deg*2*dxue(iec+1,jed:jebd)/(Info%hu(iec+1,jed:jebd,2)+Info%hu(iec+2,jed:jebd,2))
        endif
        if(Info%pe_SE /= NULL_PE) then
           xt(iec+1,jsbd:jsd) = &
             xt(iec  ,jsbd:jsd) + rad_to_deg*2*dxte(iec  ,jsbd:jsd)/(Info%ht(iec  ,jsbd:jsd,1)+Info%ht(iec+1,jsbd:jsd,1))
           xt(iec+2,jsbd:jsd) = &
             xt(iec+1,jsbd:jsd) + rad_to_deg*2*dxte(iec+1,jsbd:jsd)/(Info%ht(iec+1,jsbd:jsd,1)+Info%ht(iec+2,jsbd:jsd,1))
           xu(iec+1,jsbd:jsd) = &
             xu(iec  ,jsbd:jsd) + rad_to_deg*2*dxue(iec  ,jsbd:jsd)/(Info%hu(iec  ,jsbd:jsd,2)+Info%hu(iec+1,jsbd:jsd,2))
           xu(iec+2,jsbd:jsd) = &
             xu(iec+1,jsbd:jsd) + rad_to_deg*2*dxue(iec+1,jsbd:jsd)/(Info%hu(iec+1,jsbd:jsd,2)+Info%hu(iec+2,jsbd:jsd,2))
        endif
     endif

     deallocate(dxte)
     deallocate(dxue)
  endif

  if (Grd%cyclic_y .or. Grd%tripolar) then
     do_nth_bnd = .false.
     do_sth_bnd = .false.
     do i=isc,iec
        if (yt(i,jec+1) < yt(i,jec)  ) do_nth_bnd=.true.
        if (yt(i,jsc)   < yt(i,jsc-1)) do_sth_bnd=.true.
     enddo
     
     allocate(dytn(isbd:iebd,jsbd:jebd))
     allocate(dyun(isbd:iebd,jsbd:jebd))
     dytn(isc:iec,jsc:jec) = Grd%dytn(isc:iec,jsc:jec)
     dyun(isc:iec,jsc:jec) = Grd%dyun(isc:iec,jsc:jec)
     call mpp_update_domains(dytn(:,:), Bdom%domain2d)
     call mpp_update_domains(dyun(:,:), Bdom%domain2d)

     if (do_sth_bnd .and. .not. Grd%tripolar) then
        ! The tripolar grid has a solid southern boundary, so, we only want to alter the halos at the northern boundary
        yt(isc:iec,jsc-1) = &
          yt(isc:iec,jsc  ) - rad_to_deg*2*dytn(isc:iec,jsc-1)/(Info%ht(isc:iec,jsc  ,1)+Info%ht(isc:iec,jsc-1,1))
        yt(isc:iec,jsc-2) = &
          yt(isc:iec,jsc-1) - rad_to_deg*2*dytn(isc:iec,jsc-2)/(Info%ht(isc:iec,jsc-1,1)+Info%ht(isc:iec,jsc-2,1))
        yu(isc:iec,jsc-1) = &
          yu(isc:iec,jsc  ) - rad_to_deg*2*dyun(isc:iec,jsc-1)/(Info%hu(isc:iec,jsc  ,2)+Info%hu(isc:iec,jsc-1,2))
        yu(isc:iec,jsc-2) = &
          yu(isc:iec,jsc-1) - rad_to_deg*2*dyun(isc:iec,jsc-2)/(Info%hu(isc:iec,jsc-1,2)+Info%hu(isc:iec,jsc-2,2))

        ! Do the corners
        if (Info%pe_SW /= NULL_PE) then
           yt(isbd:isd,jsc-1) = &
             yt(isbd:isd,jsc  ) - rad_to_deg*2*dytn(isbd:isd,jsc-1)/(Info%ht(isbd:isd,jsc  ,1)+Info%ht(isbd:isd,jsc-1,1))
           yt(isbd:isd,jsc-2) = &
             yt(isbd:isd,jsc-1) - rad_to_deg*2*dytn(isbd:isd,jsc-2)/(Info%ht(isbd:isd,jsc-1,1)+Info%ht(isbd:isd,jsc-2,1))
           yu(isbd:isd,jsc-1) = &
             yu(isbd:isd,jsc  ) - rad_to_deg*2*dyun(isbd:isd,jsc-1)/(Info%hu(isbd:isd,jsc  ,2)+Info%hu(isbd:isd,jsc-1,2))
           yu(isbd:isd,jsc-2) = &
             yu(isbd:isd,jsc-1) - rad_to_deg*2*dyun(isbd:isd,jsc-2)/(Info%hu(isbd:isd,jsc-1,2)+Info%hu(isbd:isd,jsc-2,2))
        endif
        if (Info%pe_SE /= NULL_PE) then
           yt(ied:iebd,jsc-1) = &
             yt(ied:iebd,jsc  ) - rad_to_deg*2*dytn(ied:iebd,jsc-1)/(Info%ht(ied:iebd,jsc  ,1)+Info%ht(ied:iebd,jsc-1,1))
           yt(ied:iebd,jsc-2) = &
             yt(ied:iebd,jsc-1) - rad_to_deg*2*dytn(ied:iebd,jsc-2)/(Info%ht(ied:iebd,jsc-1,1)+Info%ht(ied:iebd,jsc-2,1))
           yu(ied:iebd,jsc-1) = &
             yu(ied:iebd,jsc  ) - rad_to_deg*2*dyun(ied:iebd,jsc-1)/(Info%hu(ied:iebd,jsc  ,2)+Info%hu(ied:iebd,jsc-1,2))
           yu(ied:iebd,jsc-2) = &
             yu(ied:iebd,jsc-1) - rad_to_deg*2*dyun(ied:iebd,jsc-2)/(Info%hu(ied:iebd,jsc-1,2)+Info%hu(ied:iebd,jsc-2,2))
        endif

     endif
     if (do_nth_bnd) then
        yt(isc:iec,jec+1) = &
          yt(isc:iec,jec  ) + rad_to_deg*2*dytn(isc:iec,jsc  )/(Info%ht(isc:iec,jec  ,1)+Info%ht(isc:iec,jec+1,1))
        yt(isc:iec,jec+2) = &
          yt(isc:iec,jec+1) + rad_to_deg*2*dytn(isc:iec,jsc+1)/(Info%ht(isc:iec,jec+1,1)+Info%ht(isc:iec,jec+2,1))
        yu(isc:iec,jec+1) = &
          yu(isc:iec,jec  ) + rad_to_deg*2*dyun(isc:iec,jsc  )/(Info%hu(isc:iec,jec  ,2)+Info%hu(isc:iec,jec+1,2))
        yu(isc:iec,jec+2) = &
          yu(isc:iec,jec+1) + rad_to_deg*2*dyun(isc:iec,jsc+1)/(Info%hu(isc:iec,jec+1,2)+Info%hu(isc:iec,jec+2,2))

        ! Do the corners
        if (Info%pe_NW /= NULL_PE) then
           yt(isbd:isd,jec+1) = &
             yt(isbd:isd,jec  ) + rad_to_deg*2*dytn(isbd:isd,jsc  )/(Info%ht(isbd:isd,jec  ,1)+Info%ht(isbd:isd,jec+1,1))
           yt(isbd:isd,jec+2) = &
             yt(isbd:isd,jec+1) + rad_to_deg*2*dytn(isbd:isd,jsc+1)/(Info%ht(isbd:isd,jec+1,1)+Info%ht(isbd:isd,jec+2,1))
           yu(isbd:isd,jec+1) = &
             yu(isbd:isd,jec  ) + rad_to_deg*2*dyun(isbd:isd,jsc  )/(Info%hu(isbd:isd,jec  ,2)+Info%hu(isbd:isd,jec+1,2))
           yu(isbd:isd,jec+2) = &
             yu(isbd:isd,jec+1) + rad_to_deg*2*dyun(isbd:isd,jsc+1)/(Info%hu(isbd:isd,jec+1,2)+Info%hu(isbd:isd,jec+2,2))
        endif
        if (Info%pe_NE /= NULL_PE) then
           yt(ied:iebd,jec+1) = &
             yt(ied:iebd,jec  ) + rad_to_deg*2*dytn(ied:iebd,jsc  )/(Info%ht(ied:iebd,jec  ,1)+Info%ht(ied:iebd,jec+1,1))
           yt(ied:iebd,jec+2) = &
             yt(ied:iebd,jec+1) + rad_to_deg*2*dytn(ied:iebd,jsc+1)/(Info%ht(ied:iebd,jec+1,1)+Info%ht(ied:iebd,jec+2,1))
           yu(ied:iebd,jec+1) = &
             yu(ied:iebd,jec  ) + rad_to_deg*2*dyun(ied:iebd,jsc  )/(Info%hu(ied:iebd,jec  ,2)+Info%hu(ied:iebd,jec+1,2))
           yu(ied:iebd,jec+2) = &
             yu(ied:iebd,jec+1) + rad_to_deg*2*dyun(ied:iebd,jsc+1)/(Info%hu(ied:iebd,jec+1,2)+Info%hu(ied:iebd,jec+2,2))
        endif
     endif
     deallocate(dytn)
     deallocate(dyun)
  endif
     
  allocate( datdtime_r(isd:ied,jsd:jed) )
  datdtime_r(:,:) = Grd%datr(:,:)/dtime

  ! Calculate the vectors of the sides of grid cells (in 
  ! the horizontal) for the point location scheme.
  allocate( vert_t(2,4,isd:ied,jsd:jed) )

  allocate(tile_ids(mpp_get_current_ntile(Dom%domain2d)))
  tile_ids = mpp_get_tile_id(Dom%domain2d)
  tile = tile_ids(1)    ! Assume one tile per PE
  deallocate(tile_ids)

  call get_grid_size('OCN', tile, nlon, nlat)
  allocate(lon_vert(nlon+1, nlat+1))
  allocate(lat_vert(nlon+1, nlat+1))

  call get_grid_cell_vertices('OCN', tile, lon_vert, lat_vert)

  ! In this grid configuration, we derive the verticies from a 2d configuration
  ! with shape (isc:iec+1,jsc:jec+1).
  ! The verticies with the bottom left hand corner corresponding to i,j
  !
  !     4     3
  !     +-----+
  !     | i,j |
  !     +-----+
  !     1     2
  !
  ! 1==(i,j); 2==(i+1,j); 3==(i+1,j+1); 4==(i,j+1)

  call mpp_update_domains(lon_vert(:,:), Dom%domain2d)
  vert_t(1, 1, isc:iec, jsc:jec) = lon_vert(isc:iec, jsc:jec)
  vert_t(1, 2, isc:iec, jsc:jec) = lon_vert((1+isc):(1+iec), jsc:jec)
  vert_t(1, 3, isc:iec, jsc:jec) = lon_vert((1+isc):(1+iec), (1+jsc):(1+jec))
  vert_t(1, 4, isc:iec, jsc:jec) = lon_vert(isc:iec, (1+jsc):(1+jec))

  call mpp_update_domains(lat_vert(:,:), Dom%domain2d)
  vert_t(2, 1, isc:iec, jsc:jec) = lat_vert(isc:iec, jsc:jec)
  vert_t(2, 2, isc:iec, jsc:jec) = lat_vert((1+isc):(1+iec), jsc:jec)
  vert_t(2, 3, isc:iec, jsc:jec) = lat_vert((1+isc):(1+iec), (1+jsc):(1+jec))
  vert_t(2, 4, isc:iec, jsc:jec) = lat_vert(isc:iec, (1+jsc):(1+jec))

  deallocate(lon_vert)
  deallocate(lat_vert)

  ! Now fill/get the values in the halos
  if (Info%pe_E==NULL_PE)  vert_t(:,:,ied,jsc:jec) = 0.0
  if (Info%pe_N==NULL_PE)  vert_t(:,:,isc:iec,jed) = 0.0
  if (Info%pe_W==NULL_PE)  vert_t(:,:,isd,jsc:jec) = 0.0
  if (Info%pe_S==NULL_PE)  vert_t(:,:,isc:iec,jsd) = 0.0
  if (Info%pe_NE==NULL_PE) vert_t(:,:,ied,jed)     = 0.0
  if (Info%pe_NW==NULL_PE) vert_t(:,:,isd,jed)     = 0.0
  if (Info%pe_SW==NULL_PE) vert_t(:,:,isd,jsd)     = 0.0
  if (Info%pe_SE==NULL_PE) vert_t(:,:,ied,jsd)     = 0.0

  ! The verticies are indexed from southwest, anticlockwise
  !                  4-----3
  !                  |     |
  !                  |  T  |
  !                  |     |
  !                  1-----2

  !Pre calculate some of the vectors for the point location scheme
  !See notes for details
  allocate( ij_im1j(    2,isd:ied,jsd:jed) )
  allocate( im1jm1_ijm1(2,isd:ied,jsd:jed) )
  allocate( ij_ijm1(    2,isd:ied,jsd:jed) )
  allocate( im1jm1_im1j(2,isd:ied,jsd:jed) )

  allocate( t_im1j(     2,isd:ied,jsd:jed) )
  allocate( t_ij(       2,isd:ied,jsd:jed) )
  allocate( t_ijm1(     2,isd:ied,jsd:jed) )
  allocate( t_im1jm1(   2,isd:ied,jsd:jed) )

  !North vector: U(i  ,j  )=>U(i-1,j  ), aka 3=>4
  ij_im1j(    :,isd:ied,jsd:jed) = vert_t(:,4,isd:ied,jsd:jed) - vert_t(:,3,isd:ied,jsd:jed)

  !South vector: U(i-1,j-1)=>U(i-1,j  ), aka 1=>2
  im1jm1_ijm1(:,isd:ied,jsd:jed) = vert_t(:,2,isd:ied,jsd:jed) - vert_t(:,1,isd:ied,jsd:jed)

  !East vector:  U(i  ,j  )=>U(i  ,j-1), aka 3=>2
  ij_ijm1(    :,isd:ied,jsd:jed) = vert_t(:,2,isd:ied,jsd:jed) - vert_t(:,3,isd:ied,jsd:jed)

  !West vector:  U(i-1,j-1)=>U(i-1,j  ), aka 1=>4
  im1jm1_im1j(:,isd:ied,jsd:jed) = vert_t(:,4,isd:ied,jsd:jed) - vert_t(:,1,isd:ied,jsd:jed)

  !T(i,j)=>U(i-1,j-1)
  t_im1jm1(1,isd:ied,jsd:jed) = vert_t(1,1,isd:ied,jsd:jed) - Grd%xt(isd:ied,jsd:jed)
  t_im1jm1(2,isd:ied,jsd:jed) = vert_t(2,1,isd:ied,jsd:jed) - Grd%yt(isd:ied,jsd:jed)

  !T(i,j)=>U(i  ,j-1)
  t_ijm1(1,isd:ied,jsd:jed) = vert_t(1,2,isd:ied,jsd:jed) - Grd%xt(isd:ied,jsd:jed)
  t_ijm1(2,isd:ied,jsd:jed) = vert_t(2,2,isd:ied,jsd:jed) - Grd%yt(isd:ied,jsd:jed)

  !T(i,j)=>U(i  ,j  )
  t_ij(1,isd:ied,jsd:jed) = vert_t(1,3,isd:ied,jsd:jed) - Grd%xt(isd:ied,jsd:jed)
  t_ij(2,isd:ied,jsd:jed) = vert_t(2,3,isd:ied,jsd:jed) - Grd%yt(isd:ied,jsd:jed)

  !T(i,j)=>U(i-1,j  )
  t_im1j(1,isd:ied,jsd:jed) = vert_t(1,4,isd:ied,jsd:jed) - Grd%xt(isd:ied,jsd:jed)
  t_im1j(2,isd:ied,jsd:jed) = vert_t(2,4,isd:ied,jsd:jed) - Grd%yt(isd:ied,jsd:jed)

  allocate(Info%maxlon(jsc:jec))
  allocate(Info%minlon(jsc:jec))
  allocate(Info%maxlat(isc:iec))
  allocate(Info%minlat(isc:iec))

  do j=jsc,jec
     Info%minlon(j) = minval(vert_t(1,:,isc:iec,j))
     Info%maxlon(j) = maxval(vert_t(1,:,isc:iec,j))
  enddo
  do i=isc,iec
     Info%minlat(i) = minval(vert_t(2,:,i,jsc:jec))
     Info%maxlat(i) = maxval(vert_t(2,:,i,jsc:jec))
  enddo

end subroutine blob_util_init
! </SUBROUTINE>  NAME="blob_util_init"


!######################################################################
! <SUBROUTINE NAME="blob_chksum">
!
! <DESCRIPTION>
! Performs global sums and checksums for all blob types (for diagnostic
! purposes).  
! </DESCRIPTION>
!
subroutine blob_chksum(T_prog, head_static_free, head_static_bott, &
                       head_dynamic_free, head_dynamic_bott, blob_counter)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  type(ocean_blob_type),        pointer    :: head_static_free
  type(ocean_blob_type),        pointer    :: head_static_bott
  type(ocean_blob_type),        pointer    :: head_dynamic_free
  type(ocean_blob_type),        pointer    :: head_dynamic_bott
  integer, dimension(isc:iec,jsc:jec,nk), intent(in) :: blob_counter

  type(ocean_blob_type), pointer                       :: this=>NULL()
  real, dimension(isd:ied,jsd:jed,nk,num_prog_tracers) :: grdtracer, grdfield
  real, dimension(isd:ied,jsd:jed,nk)    :: grdlat, grdlon, grddepth, grdgeodepth
  real, dimension(isd:ied,jsd:jed,nk)    :: grdu, grdv, grdw, grdent
  real, dimension(isd:ied,jsd:jed,nk)    :: grdmass, grddens, grdvolume
  real, dimension(isd:ied,jsd:jed,nk)    :: grdh1, grdh2, grdstep, grdst
  integer, dimension(isd:ied,jsd:jed,nk) :: grdhash, grdnum 
  integer, dimension(isd:ied,jsd:jed,nk) :: grdnsteps, grdmsteps
  integer, dimension(isd:ied,jsd:jed,nk) :: grdi, grdj, grdk
  integer           :: tmpi
  integer           :: i, j, k, p, n, nblobs
  character(len=28) :: tname
  character(len=28) :: fname
  character(len=*), parameter :: fmti="(a,x,i25)"
  character(len=*), parameter :: fmte="(a,x,es25.18)"
  integer :: stdoutunit

  real :: tmpr

  stdoutunit=stdout()

  grdlat(:,:,:)      = 0.0
  grdlon(:,:,:)      = 0.0
  grddepth(:,:,:)    = 0.0
  grdgeodepth(:,:,:) = 0.0
  grdst(:,:,:)       = 0.0
  grdmass(:,:,:)     = 0.0
  grddens(:,:,:)     = 0.0
  grdvolume(:,:,:)   = 0.0
  grdu(:,:,:)        = 0.0
  grdv(:,:,:)        = 0.0
  grdw(:,:,:)        = 0.0
  grdent(:,:,:)      = 0.0
  grdh1(:,:,:)       = 0.0
  grdh2(:,:,:)       = 0.0
  grdstep(:,:,:)     = 0.0
  grdtracer(:,:,:,:) = 0.0
  grdfield(:,:,:,:)  = 0.0
  grdhash(:,:,:)     = 0
  grdnum(:,:,:)      = 0
  grdnsteps(:,:,:)   = 0
  grdmsteps(:,:,:)   = 0
  grdi(:,:,:)        = 0
  grdj(:,:,:)        = 0
  grdk(:,:,:)        = 0
  nblobs             = 0 !number of blobs

  do p=1,4
     nullify(this)
     if(p==1 .and. associated(head_static_free))  this=>head_static_free
     if(p==2 .and. associated(head_static_bott))  this=>head_static_bott
     if(p==3 .and. associated(head_dynamic_free)) this=>head_dynamic_free
     if(p==4 .and. associated(head_dynamic_bott)) this=>head_dynamic_bott

     if (associated(this)) then
        fullcycle: do
           i = this%i - Dom%ioff
           j = this%j - Dom%joff
           k = this%k

           if (isc<=i .and. i<=iec .and. jsc<=j .and. j<=jec) then
              grdlat(i,j,k)      = grdlat(i,j,k)      + this%lat
              grdlon(i,j,k)      = grdlon(i,j,k)      + this%lon
              grddepth(i,j,k)    = grddepth(i,j,k)    + this%depth
              grdgeodepth(i,j,k) = grdgeodepth(i,j,k) + this%geodepth
              grdst(i,j,k)       = grdst(i,j,k)       + this%st
              grdmass(i,j,k)     = grdmass(i,j,k)     + this%mass
              grddens(i,j,k)     = grddens(i,j,k)     + this%density
              grdvolume(i,j,k)   = grdvolume(i,j,k)   + this%volume
              grdu(i,j,k)        = grdu(i,j,k)        + this%v(1)
              grdv(i,j,k)        = grdv(i,j,k)        + this%v(2)
              grdw(i,j,k)        = grdw(i,j,k)        + this%v(3)
              grdent(i,j,k)      = grdent(i,j,k)      + this%ent
              grdh1(i,j,k)       = grdh1(i,j,k)       + this%h1
              grdh2(i,j,k)       = grdh2(i,j,k)       + this%h2
              grdhash(i,j,k)     = grdhash(i,j,k)     + this%hash
              grdnum(i,j,k)      = grdnum(i,j,k)      + this%number
              grdstep(i,j,k)     = grdstep(i,j,k)     + this%step
              grdnsteps(i,j,k)   = grdnsteps(i,j,k)   + this%nsteps
              grdmsteps(i,j,k)   = grdmsteps(i,j,k)   + this%model_steps
              grdi(i,j,k)        = grdi(i,j,k)        + this%i
              grdj(i,j,k)        = grdj(i,j,k)        + this%j
              grdk(i,j,k)        = grdk(i,j,k)        + this%k
              nblobs             = nblobs + 1
              do n=1,num_prog_tracers
                 grdtracer(i,j,k,n) = grdtracer(i,j,k,n) + this%tracer(n)
                 if (p==3.or.p==4) grdfield(i,j,k,n)  = grdfield(i,j,k,n) + this%field(n)
              enddo
           endif

           this=>this%next
           if(.not. associated(this)) exit fullcycle
        enddo fullcycle
     endif
  enddo

  
  call mpp_sum(nblobs)

  write(stdoutunit, fmti)    'total number of blobs               =', nblobs

  tmpi = sum(blob_counter(:,:,:))
  call mpp_sum(tmpi)
  write(stdoutunit, fmti) 'global number of blobs for all time =', tmpi

  tmpr = mpp_global_sum(Dom%domain2d, grdmass(isc:iec,jsc:jec,1:nk),global_sum_flag)
  write(stdoutunit, fmte) 'global mass of blobs                =', tmpr
  do n=1,num_prog_tracers
     tmpr = mpp_global_sum(Dom%domain2d, grdtracer(isc:iec,jsc:jec,1:nk,n),global_sum_flag)
     if(n==index_temp) then
             write(stdoutunit, fmte) 'global heat content of blobs        =', tmpr
     else
        tname = trim(T_prog(n)%name)//' content of blobs'
        fname = trim(T_prog(n)%name)//' concentration of blobs'
        write(stdoutunit, fmte) 'global '//tname(:)//' =', tmpr
     endif
  enddo

  call write_chksum_3d_int('blob counter', blob_counter(COMP,1:nk))
  call write_chksum_3d('blob latitude', grdlat(COMP,1:nk))  
  call write_chksum_3d('blob longitude', grdlon(COMP,1:nk))
  call write_chksum_3d('blob depth', grddepth(COMP,1:nk))
  call write_chksum_3d('blob geodepth', grdgeodepth(COMP,1:nk))
  call write_chksum_3d('blob vertical coord. (st)', grdst(COMP,1:nk))
  call write_chksum_3d('blob mass', grdmass(COMP,1:nk))
  call write_chksum_3d('blob density', grddens(COMP,1:nk))
  call write_chksum_3d('blob volume', grdvolume(COMP,1:nk))
  call write_chksum_3d('blob zonal velocity', grdu(COMP,1:nk))
  call write_chksum_3d('blob meridional velocity', grdv(COMP,1:nk))
  call write_chksum_3d('blob vertical velocity', grdw(COMP,1:nk))
  call write_chksum_3d('blob entrainment velocity', grdent(COMP,1:nk))
  call write_chksum_3d('stretching function 1', grdh1(COMP,1:nk))
  call write_chksum_3d('stretching function 2', grdh2(COMP,1:nk))
  call write_chksum_3d_int('blob number', grdnum(COMP,1:nk))
  call write_chksum_3d('blob step size', grdstep(COMP,1:nk))
  call write_chksum_3d_int('number of blob steps', grdnsteps(COMP,1:nk))
  call write_chksum_3d_int('number of E system steps', grdmsteps(COMP,1:nk))
  call write_chksum_3d_int('zonal cell index (i)', grdi(COMP,1:nk))
  call write_chksum_3d_int('meridional cell index (j)', grdj(COMP,1:nk))
  call write_chksum_3d_int('depth cell index (k)', grdk(COMP,1:nk))
  do n=1,num_prog_tracers    
     if (n==index_temp) then
        call write_chksum_3d('blob heat content', grdtracer(COMP,1:nk,n))
        call writE_chksum_3d('blob temp concentration', grdfield(COMP,1:nk,n))
     else
        tname = trim(T_prog(n)%name)//' content'
        call write_chksum_3d('blob '//tname(1:19), grdtracer(COMP,1:nk,n))
        call write_chksum_3d(fname(1:19), grdfield(COMP,1:nk,n))
     endif
  enddo
  write(stdoutunit, '(a)') 'end blob chksums'

end subroutine blob_chksum
! </SUBROUTINE>  NAME="blob_chksum"

!######################################################################
! <SUBROUTINE NAME="lagrangian_system_chksum">
!
! <DESCRIPTION>
!
! Performs checksums for the Lagrangian_system derived type.  This is
! the derived type that stores all of the "gridded" blob variables, 
! and is essential for the accounting required to interact with the 
! Eulerian model in a conservative manner.  The checksums are for
! diagnostic purposes.
!
! </DESCRIPTION>
!
subroutine lagrangian_system_chksum(L_system)

  type(ocean_lagrangian_type), intent(in) :: L_system

  integer :: stdoutunit
  character(len=*), parameter :: fmt="(a,x,i20)"

  stdoutunit=stdout()

  call write_chksum_3d('T-grid upper cell rho_dzt (taup1)', L_system%rho_dztup(COMP,1:nk))
  call write_chksum_3d('T-grid lower cell rho_dzt (taup1)', L_system%rho_dztlo(COMP,1:nk))
  call write_chksum_3d('T-grid blob convergence (taup1)', L_system%conv_blob(COMP,1:nk))
  write(stdoutunit, fmt) 'end Lagrangian system chksums'

end subroutine lagrangian_system_chksum
! </SUBROUTINE>  NAME="lagrangian_system_chksum"

!######################################################################
! <SUBROUTINE NAME="E_and_L_totals">
!
! <DESCRIPTION>
! Gives a brief summary of the total mass, volume and tracer content
! of the E, L and total systems.  Usually used for debuggin purposes.
! </DESCRIPTION>
!
subroutine E_and_L_totals(L_system, Thickness, T_prog, idx)
  type(ocean_lagrangian_type),  intent(in) :: L_system
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer,                      intent(in) :: idx

  real, dimension(isd:ied,jsd:jed) :: tmpE, tmpL1, tmpL2, tmpT
  real :: tmpE_total, tmpL1_total, tmpL2_total, tmpT_total
  character(len=128) :: tname
  integer :: n, k
  integer :: stdoutunit

  stdoutunit=stdout()

  tmpE  = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%rho_dzt( :,:,:,idx),3)
  tmpL1 = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%rho_dztL(:,:,:,idx),3)
  do k=1,nk
     tmpL2 = Grd%tmask(:,:,k)*Grd%dat(:,:)*(L_system%rho_dztup(:,:,k)+L_system%rho_dztlo(:,:,k))
  enddo
  tmpT  = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%rho_dztT(:,:,:,idx),3)
  call maketotal(.true.)
  write(stdoutunit,'(a,x,es25.18,x,a)') 'E mass                  =',tmpE_total, 'kg'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'L mass (from thickness) =',tmpL1_total,'kg'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'L mass (from L_system)  =',tmpL2_total,'kg'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'E + L mass              =',tmpE_total+tmpL1_total,'kg'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'T mass (from rho_dztT)  =',tmpT_total,'kg'

  tmpE  = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%dzt( :,:,:),3)
  tmpL1 = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%dztL(:,:,:),3)
  tmpT  = Grd%dat(:,:)*sum(Grd%tmask(:,:,:)*Thickness%dztT(:,:,:,idx),3)
  call maketotal(.true.)
  write(stdoutunit,'(a,x,es25.18,x,a)') 'E volume                  =',tmpE_total, 'm^3'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'L volume (from thickness) =',tmpL1_total,'m^3'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'E + L volume              =',tmpE_total+tmpL1_total,'m^3'
  write(stdoutunit,'(a,x,es25.18,x,a)') 'T volume (from dztT)      =',tmpT_total,'m^3'

  do n=1,num_prog_tracers
     tmpE  = T_prog(n)%conversion*Grd%dat(:,:)                   &
             *sum( Grd%tmask(:,:,:)*Thickness%rho_dzt(:,:,:,idx)*T_prog(n)%field(:,:,:,idx), 3)
     tmpL1 = T_prog(n)%conversion*sum(Grd%tmask(:,:,:)*T_prog(n)%sum_blob(:,:,:,idx), 3)
     tmpT  = tmpE + tmpL1
     call maketotal(.false.)
     if(n==index_temp) then
        write(stdoutunit,'(a,x,es25.18,x,a)') 'E heat       =',tmpE_total, 'J'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'L heat       =',tmpL1_total,'J'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'T heat       =',tmpT_total, 'J'
     elseif(T_prog(n)%name(1:3)=='age') then
        tname = trim(T_prog(n)%name)
        write(stdoutunit,'(a,x,es25.18,x,a)') 'E '//tname(1:10)//' =',tmpE_total, 'yr'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'L '//tname(1:10)//' =',tmpL1_total,'yr'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'T '//tname(1:10)//' =',tmpT_total, 'yr'
     else
        tname = trim(T_prog(n)%name)
        write(stdoutunit,'(a,x,es25.18,x,a)') 'E '//tname(1:10)//' =',tmpE_total, 'kg'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'L '//tname(1:10)//' =',tmpL1_total,'kg'
        write(stdoutunit,'(a,x,es25.18,x,a)') 'T '//tname(1:10)//' =',tmpT_total, 'kg'
     endif
  enddo

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! This is a nested subroutine that does the global sum of an array for !
  ! each of the E system, L system and the combined system.              !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine maketotal(do_L2)
    logical, intent(in) :: do_L2
    tmpE_total  = mpp_global_sum(Dom%domain2d, tmpE,  global_sum_flag)
    tmpL1_total = mpp_global_sum(Dom%domain2d, tmpL1, global_sum_flag)
    tmpT_total  = mpp_global_sum(Dom%domain2d, tmpT,  global_sum_flag)
    if (do_L2) then 
       tmpL2_total = mpp_global_sum(Dom%domain2d, tmpL2, global_sum_flag)
    endif
  end subroutine maketotal

end subroutine E_and_L_totals
! </SUBROUTINE>  NAME="E_and_L_totals"

!######################################################################
! <SUBROUTINE NAME="write_blobs">
!
! <DESCRIPTION>
! Dumps most of the information carried around by blobs, for all blobs
! in a particular list.  Useful for debugging.
! </DESCRIPTION>
!
subroutine write_blobs(head, head_name, time)
  type(ocean_blob_type), pointer    :: head
  character(len=*),      intent(in) :: head_name
  character(len=*),      intent(in) :: time

  type(ocean_blob_type), pointer :: this=>NULL()
  integer :: i,j,n
  integer :: stdoutunit

  stdoutunit = stdout()

  write(stdoutunit, '(a)') ' '
  write(stdoutunit, '(a)') 'Summary of all blobs in '//trim(head_name)//' list ('//trim(time)//')'
  write(stdoutunit,'(4(a3,x),2(a6,x),4(a6,x),2(a10,x),1(a8,x),2(a10,x),2(a6,x),3(a9,x),(a10))') &
       'pe','i','j','k', &
       'hash','number', &
       'lat','lon','depth', 'gdepth',&
       'mass','volume','density', &
       'heat','S (kg)', &
       'temp','sal', &
       'u','v','w', &
       'dzt'
  n=0
  if(associated(head)) then
     this=>head
     blobcycle: do
        i=this%i; j=this%j
        print('(4(i3,x),2(i6,x),4(f6.1,x),2(es10.3,x),(f8.2,x),2(es10.3,x),2(f6.2,x),3(es9.2,x),(es10.3))'), &
             Info%pe_this, i, j, this%k,&
             this%hash, this%number, &
             this%lat, this%lon, this%depth, this%geodepth,&
             this%mass, this%volume, &
             this%density,&
             this%tracer(index_temp), this%tracer(index_salt), &
             this%field(index_temp), this%field(index_salt),&
             this%v(1), this%v(2), this%v(3), &
             this%volume*Grd%datr(i,j)
        n=n+1
        this=>this%next
        if(.not.associated(this)) exit blobcycle
     enddo blobcycle
  endif
  print('(a,i3,a,i10)'), 'total '//trim(head_name)//' blobs on pe ',Info%pe_this,' = ', n

end subroutine write_blobs
! </SUBROUTINE>  NAME="write_blobs"

!######################################################################
! <SUBROUTINE NAME="blob_delete">
!
! <DESCRIPTION>
! Deletes all (nearly) zero mass blob objects from the linked list.
! The size of the blobs that are deleted is controlled by the variable
! blob_small_mass in the ocean_blob_nml.
! </DESCRIPTION>
!
subroutine blob_delete(Time, Thickness, T_prog, head)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_blob_type),        pointer       :: head

  ! local variables
  type(ocean_blob_type), pointer :: prev => NULL()
  type(ocean_blob_type), pointer :: this => NULL()
  type(ocean_blob_type), pointer :: next => NULL()
  integer :: i,j,k,tau,taup1

  taup1 = Time%taup1
  tau   = Time%tau

  ! point this to the head blob
  this => head
  
  !====================================================================!
  ! if the head is not associated, the list is empty and there is      !
  ! nothing to delete if it is associated, then cycle through the list,!
  ! deleting all zero-mass blobs until the end (tail) of the list is   !
  ! reached                                                            !
  !====================================================================!

  ! check if the head is associated
  if(associated(this)) then
     next => this%next

     blobdel: do
        i = this%i
        j = this%j
        k = this%k

        if(this%mass < blob_small_mass) then
           
           ! return any blob properties to the E system
           call kill_blob(Thickness, T_prog(:), this, i, j, k)
           
           ! Deallocate the blob from memory
           call free_blob_memory(this)
           
           ! Unlink the blob from the list. Point the previous blob 
           ! to the next blob and vice versa.
           if(associated(prev)) then
              prev%next=>next
           else
              ! if the prev is not associated, we are at the head
              ! and we need to point the head towards the next
              head => next
              if(associated(next)) next%prev=>NULL()
           endif
           
           ! this is the vice versa bit
           if(associated(next)) next%prev=>prev
           
        endif ! blob mass small?

        ! check to see whether we are at the end of the list.  If so, exit
        ! the loop.  If not, then move our current, prev and next blobs
        ! further down the list and go back to the beginning of the do loop
        if(associated(next)) then
           this => next
           prev => this%prev
           next => this%next
        else
           exit blobdel
        endif

     enddo blobdel

     nullify(this)
     nullify(prev)
     nullify(next)

  endif !head associated?

end subroutine blob_delete
! </SUBROUTINE>  NAME="blob_delete"

!######################################################################
! <SUBROUTINE NAME="unlink_blob">
!
! <DESCRIPTION>
! Unlinks a blob from a doubly linked list.  It returns pointers to 
! the blob, the head of the list, the (formerly) previous blob in the 
! list and the (formerly) next blob in the list.
! </DESCRIPTION>
!
subroutine unlink_blob(blob, head, prev, next)
  type(ocean_blob_type), pointer :: blob
  type(ocean_blob_type), pointer :: head
  type(ocean_blob_type), pointer :: prev
  type(ocean_blob_type), pointer :: next
  prev=>blob%prev
  next=>blob%next
  ! Unlink the list from the blob
  if(associated(prev)) then
     prev%next=>next
  else
     head=>next
     if(associated(next)) next%prev=>NULL()
  endif
  if(associated(next)) next%prev=>prev
  ! Unlink the blob from the list
  blob%next=>NULL()
  blob%prev=>NULL()
end subroutine unlink_blob
! </SUBROUTINE>  NAME="unlink_blob"

!######################################################################
! <SUBROUTINE NAME="insert_blob">
!
! <DESCRIPTION>
! Inserts a blob to the linked list.  The relative order of blobs in 
! a linked list determines whether bitwise reproduction is possible.
!
! Regardless of bitwise reproducability or not, we must ensure that 
! blobs always appear in the same relative order when we are using 
! dynamic blobs because if we have a situation where dztL>dztT, we 
! start destroying blobs to enforce dztL<dztT.  In order that we do 
! not significantly change answers, we must always destroy the same 
! blob, regardless of domain decomposition, restarts, etc. So, we must 
! always sort blobs so they appear in the linked list in the same 
! relative order.
! </DESCRIPTION>
!
subroutine insert_blob(blob, head)

  type(ocean_blob_type), pointer :: blob
  type(ocean_blob_type), pointer :: head

  ! local variables
  type(ocean_blob_type), pointer :: prev
  type(ocean_blob_type), pointer :: next
  logical :: order, eol
  integer :: stdoutunit 

  stdoutunit = stdout()

  if (associated(head)) then
     prev => NULL()
     next => head

     ! cycle through hash, looking for a spot
     checkhash: do
        call inorder(blob%hash, next%hash, order, eol)
        if (eol .or. order) exit checkhash
     enddo checkhash
     if (eol) return
     
     ! cycle through the blob number, ensuring we remain
     ! in the correct hash
     checknumber: do
        if (blob%hash > next%hash) then
           order = .true.
        elseif(blob%hash == next%hash) then
           call inorder(blob%number, next%number, order, eol)
        else
           write(stdoutunit,'(a)') 'ocean_blob_util_mod, insert_blob: list out of order'
           call debugoutput('blob information (hash and number)', prev, blob, next)
           call mpp_error(FATAL, 'ocean_blob_util_mod, insert_blob: list out of order')
        endif
        if (eol .or. order) exit checknumber
     enddo checknumber
     if (eol) return
     
     ! if we have made it this far, then we have not reached the end of the 
     ! list and the blob should be in its correct place in the list, so 
     ! insert it into the list between prev and next, or if we are at the 
     ! top of the list, insert it at the head.
     if (associated(prev)) then
        ! insert the blob after the previous blob
        prev%next => blob
        blob%prev => prev
     else
        ! insert the blobs at the head of the list
        head      => blob
        blob%prev => NULL()
     endif
     next%prev => blob
     blob%next => next

  else  !head associated

     ! the list is empty, so make the blob the only item in the list
     head      => blob
     blob%next => NULL()
     blob%prev => NULL()

  endif !head associated

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that checks if the blob is in the right spot
! in the linked list
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine inorder(var1, var2, foundslot, eolist)

    integer, intent(in)  :: var1, var2
    logical, intent(out) :: foundslot, eolist

    foundslot = .false.
    eolist    = .false.

    if (var1 >= var2) then
       foundslot = .true.
    else
       call checkeol(eolist)
       if (.not. eolist) then
          ! move down the list
          prev => next
          next => next%next
       endif
    endif

  end subroutine inorder

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that checks if we are at the end of the 
! linked list
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine checkeol(endoflist)

    logical :: endoflist
    
    endoflist = .false.

    if(.not. associated(next%next)) then
       ! we are at the end of the list
       next%next => blob
       blob%prev => next
       blob%next => NULL()
       endoflist = .true.
    endif

  end subroutine checkeol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a nested subroutine that outputs useful debugging information
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine debugoutput(message, pr, th, ne)

    type(ocean_blob_type), pointer :: pr, th, ne
    character(len=*), intent(in) :: message

    print ('(/,a)'), trim(message)
    if(associated(pr)) then
       print('(a5,2i6)'), 'prev:', pr%hash, pr%number
    else
       print('(a19)'), 'prev not associated'
    endif
    if(associated(th)) then
       print('(a5,2i6)'), 'this:', th%hash, th%number
    else
       print('(a19)'), 'this not associated'
    endif
    if(associated(ne)) then
       print('(a5,2i6)'), 'next:', ne%hash, ne%number
    else
       print('(a19)'), 'next not associated'
    endif

  end subroutine debugoutput

end subroutine insert_blob
! </SUBROUTINE>  NAME="insert_blob"

!######################################################################
! <SUBROUTINE NAME="count_blob">
!
! <DESCRIPTION>
! Allocates a blob its hash and a number.  These two numbers can
! uniquely identify any blob.  The hash and number is based on the grid
! cell of origin.  Each grid cell has a unique hash.  We have an array 
! which keeps track of the number of blobs formed in a grid cell.  These
! two numbers give the unique identifier.  So, we also need to increment
! the counter array.
! </DESCRIPTION>
!
subroutine count_blob(blob, blob_counter)

  type(ocean_blob_type) :: blob
  integer, dimension(isc:iec,jsc:jec,nk), intent(inout) :: blob_counter

  integer :: ii, jj, kk

  blob%hash = hashfun(blob%i, blob%j, blob%k)

  ii = blob%i
  jj = blob%j
  kk = blob%k

  blob%number            = blob_counter(ii, jj, kk) + 1
  blob_counter(ii,jj,kk) = blob%number
  
end subroutine count_blob
! </SUBROUTINE>  NAME="count_blob"


!######################################################################
! <SUBROUTINE NAME="put_att">
!
! <DESCRIPTION>
! Writes an attribute to a netcdf file.
! </DESCRIPTION>
!
subroutine put_att(ncid, id, att, attval)
  integer,           intent(in) :: ncid, id
  character (len=*), intent(in) :: att, attval
  integer :: vallen, mret
  integer :: stderrunit

  stderrunit=stderr()

  ! Define the attribute for the netcdf file
  vallen=len_trim(attval)
  mret = nf_put_att_text(ncid, id, att, vallen, attval)
  ! If something goes wrong give error and bring the model down
  if (mret .ne. NF_NOERR) then
     write(stderrunit,'(a)') 'ocean_blob_util_mod, putt_att: '&
          //'nf_put_att_text failed adding '//trim(att)//' = '//trim(attval)
     write(stderrunit,'(3x,a,i3)') 'error code = ',mret
     call give_error_code(mret)
     call error_mesg('ocean_blob_util_mod, put_att', &
          'netcdf function returned a failure!', FATAL)
 endif
end subroutine put_att
! </SUBROUTINE>  NAME="put_att"

!######################################################################
! <FUNCTION NAME="inq_var">
!
! <DESCRIPTION>
! Gets the variable identifier from a netcdf file.
! </DESCRIPTION>
!
integer function inq_var(ncid, var)
  integer,           intent(in) :: ncid
  character(len=*),  intent(in) :: var
  integer :: iret
  integer :: stderrunit 

  stderrunit = stderr()
  iret = nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
     write(stderrunit,*) 'ocean_blob_util_mod, inq_var: nf_inq_varid ',&
          var,' failed'
     call error_mesg('ocean_blob_util_mod, inq_var', &
          'netcdf function returned a failure!', FATAL)
   endif
end function inq_var
! </FUNCTION>  NAME="inq_var"

!######################################################################
! <FUNCTION NAME="get_double">
!
! <DESCRIPTION>
! Gets the value of a "double" variable from a netcdf file
! </DESCRIPTION>
!
real function get_double(ncid, id, m)
  integer, intent(in) :: ncid, id, m
  integer :: iret
  integer :: stderrunit 

  stderrunit = stderr()
  iret=nf_get_var1_double(ncid, id, m, get_double)
  if (iret .ne. NF_NOERR) then
     write(stderrunit,'(a)') 'ocean_blob_util_mod, get_double: ' &
          //'nf_get_var1_double failed reading'
     call error_mesg('ocean_blob_util_mod, get_double', &
          'netcdf function returned a failure!', FATAL)
  endif
end function get_double
! </FUNCTION>  NAME="get_double"

!######################################################################
! <FUNCTION NAME="get_int">
!
! <DESCRIPTION>
! Gets the value of an integer variable from a netcdf file
! </DESCRIPTION>
!
integer function get_int(ncid, id, m)
  integer, intent(in) :: ncid, id, m
  integer :: iret
  integer :: stderrunit 

  stderrunit = stderr()
  iret=nf_get_var1_int(ncid, id, m, get_int)
  if (iret .ne. NF_NOERR) then
     write(stderrunit,'(a)') 'ocean_blob_util_mod, get_int: '&
          //'nf_get_var1_double failed reading'
     call error_mesg('ocean_blob_util_mod, get_int', &
          'netcdf function returned a failure!', FATAL)
  endif
end function get_int
! </FUNCTION>  NAME="get_int"

!######################################################################
! <SUBROUTINE NAME="put_double">
!
! <DESCRIPTION>
! Writes the value of a "double" variable to a netcdf file
! </DESCRIPTION>
!
subroutine put_double(ncid, varid, start, val)
  integer, intent(in) :: ncid, varid, start
  real,    intent(in) :: val
  integer :: mret, mret1
  character(len=31) :: varname
  integer :: stderrunit 

  stderrunit = stderr()
  mret = nf_put_vara_double(ncid, varid, start, 1, val)
  if (mret .ne. NF_NOERR) then
     mret1 = nf_inq_varname(ncid, varid, varname)
     write(stderrunit,*) 'ocean_blob_util_mod, put_double: '&
          //'nf_put_vara_double failed writing'//trim(varname)
     write(stderrunit, '(a,i10)') 'Failed with error code: ', mret
     print *, 'ocean_blob_util_mod, put_double: '&
          //'nf_put_vara_double failed writing '//trim(varname)
     print *, 'Failed with error code: ', mret
     call give_error_code(mret)
     call error_mesg('ocean_blob_util_mod, put_double', &
          'netcdf function returned a failure!', FATAL)
  endif
end subroutine put_double
! </SUBROUTINE>  NAME="put_double"


!######################################################################
! <SUBROUTINE NAME="put_int">
!
! <DESCRIPTION>
! Writes the value of an integer variable to a netcdf file
! </DESCRIPTION>
!
subroutine put_int(ncid, varid, start, val)
  integer, intent(in) :: ncid, varid, start, val
  integer :: mret, mret1
  character(len=31) :: varname
  integer :: stderrunit 

  stderrunit = stderr()
  mret = nf_put_vara_int(ncid, varid, start, 1, val)
  if (mret .ne. NF_NOERR) then
     mret1 = nf_inq_varname(ncid, varid, varname)
     write(stderrunit,*) 'ocean_blob_util_mod, put_int: '&
          //'nf_put_vara_int failed writing '//trim(varname)
     write(stderrunit, '(a,i10)') 'Failed with error code: ', mret
     print *, 'ocean_blob_util_mod, put_int: '&
          //'nf_put_vara_int failed writing '//trim(varname)
     print *, 'Failed with error code: ', mret
     call give_error_code(mret)
     call error_mesg('ocean_blob_util_mod, put_int', &
          'netcdf function returned a failure!', FATAL)
  endif
end subroutine put_int
! </SUBROUTINE>  NAME="put_int"

!######################################################################
! <FUNCTION NAME="def_var">
!
! <DESCRIPTION>
! Defines a netcdf variable
! </DESCRIPTION>
!
integer function def_var(ncid, var, ntype, idim)
  integer,           intent(in) :: ncid
  integer,           intent(in) :: ntype 
  integer,           intent(in) ::idim
  character(len=*),  intent(in) :: var
  integer :: mret
  character(len=30) :: problem
  integer :: stderrunit 

  stderrunit = stderr()
  mret = nf_def_var(ncid, var, ntype, 1, idim, def_var)

  if (mret .ne. NF_NOERR) then
     write(stderrunit,'(a,i4)') 'ocean_blob_util_mod, def_var: '//&
          'nf_def_var failed for '//trim(var)//'with error code ',mret
     write(stderrunit,'(3x,a,i3)') 'error code = ',mret
     call give_error_code(mret)
     call error_mesg('ocean_blob_util_mod, '//trim(problem), &
          'netcdf function returned a failure!', FATAL)
  endif
end function def_var
! </FUNCTION>  NAME="def_var"

!######################################################################
! <SUBROUTINE NAME="give_error_code">
!
! <DESCRIPTION>
! Gives error descriptions for netcdf calls.
! </DESCRIPTION>
!
subroutine give_error_code(retval)
  integer :: retval
  print *, ' '
  if(retval == NF_EBADDIM)    print *, 'Bad dimension'
  if(retval == NF_ENAMEINUSE) print *, 'Name already in use'
  if(retval == NF_ENOTVAR)    print *, 'Not a variable'
end subroutine give_error_code
! </SUBROUTINE>  NAME="give_error_code"

!######################################################################
! <FUNCTION NAME="hashfun">
!
! <DESCRIPTION>
! Calculates the hash
! </DESCRIPTION>
!
integer function hashfun(i,j,k)
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  hashfun = k + nk*j + Info%nk_nj*i
end function hashfun
! </FUNCTION>  NAME="hashfun"

!######################################################################
! <SUBROUTINE NAME="blob_util_end">
!
! <DESCRIPTION>
! Does what is necessary to shut down the module.
! </DESCRIPTION>
!
subroutine blob_util_end()
  
  nullify(Info)
  nullify(Grd)
  nullify(Dom)
end subroutine blob_util_end
! </SUBROUTINE>  NAME="blob_util_end"

!######################################################################
! <SUBROUTINE NAME="check_ijcell">
!
! <DESCRIPTION>
! Checks whether a blob (horizontally) resides in a grid cell or not.  
! If it does not it figures out which direction the blob is in and 
! checks the neighbouring grid cell, until it finds which grid cell 
! the blob resides in.
!
! It uses a cross product technique from computational geometry
! (Cormen et al., 2001).
! </DESCRIPTION>
!
subroutine check_ijcell(dx, dy, i, j, h, a, lon, lat, off)
  real,                  intent(in)    :: dx
  real,                  intent(in)    :: dy
  integer,               intent(inout) :: i
  integer,               intent(inout) :: j
  real,    dimension(2), intent(inout) :: h
  real,    dimension(2), intent(in)    :: a
  real,                  intent(out)   :: lon
  real,                  intent(out)   :: lat
  logical, dimension(2), intent(out)   :: off

  real, dimension(2) :: b, a_b, b_ij, b_im1jm1, t_b
  logical :: check_north, check_south, check_east, check_west, doublecheck
  integer :: new_i, new_j
  integer :: stdoutunit

  stdoutunit = stdout()

  ! set default values
  check_north = .true. 
  check_south = .true. 
  check_east  = .true. 
  check_west  = .true.
  doublecheck = .true.

  ! Update the blob longitude and latitude and calculate the required vectors
  a_b(1) = dx/h(1)
  a_b(2) = dy/h(2)
  a_b(:) = rad_to_deg*a_b(:)
  b(:)   = a(:) + a_b(:)

  new_i = i
  new_j = j

  ! The verticies are indexed from southwest, anticlockwise
  ! 4-----3    U(i-1,j  ) == im1j   == 4
  ! |     |    U(i  ,j  ) ==   ij   == 3
  ! |  T  |     tracer cell, T(i,j)
  ! |   a-+-b  U(i   j-1) ==   ijm1 == 2
  ! 1-----2    U(i-1,j-1) == im1jm1 == 1
  ! a is the start position of the blob and b is the end position

  checkcell: do while (doublecheck)
          
     doublecheck = .false.

     b_im1jm1(:) = vert_t(:,1,i,j) - b(:)
     b_ij(:)     = vert_t(:,3,i,j) - b(:)
     t_b(1)      = b(1) - xt(i,j)
     t_b(2)      = b(2) - yt(i,j)

     if (check_north) then
        if( cross( t_im1j(:,i,j), t_b(:)         ) < 0.0 .and. & !T=>i-1j x T =>b
            cross( t_ij(:,i,j),   t_b(:)         ) > 0.0 .and. & !T=>ij   x T =>b
            cross( b_ij(:),       ij_im1j(:,i,j) ) < 0.0 ) then  !b=>ij   x ij=>i-1j
           new_j=new_j+1
           check_south = .false.
           doublecheck = .true.
        endif
     endif

     if (check_south) then
        if( cross( t_ijm1(:,i,j),   t_b(:)             ) < 0.0 .and. & !T=>ij-1   x T     =>b
            cross( t_im1jm1(:,i,j), t_b(:)             ) > 0.0 .and. & !b=>i-1j-1 x T     =>b
            cross( b_im1jm1(:),     im1jm1_ijm1(:,i,j) ) < 0.0 ) then  !b=>i-1j-1 x i-1j-1=>ij-1
           new_j=new_j-1
           check_north = .false.
           doublecheck = .true.
        endif
     endif

     if (check_east) then
        if( cross( t_ij(:,i,j),   t_b(:)         ) < 0.0 .and. & !T=>ij   x T =>b
            cross( t_ijm1(:,i,j), t_b(:)         ) > 0.0 .and. & !T=>ij-1 x T =>b
            cross( b_ij(:),       ij_ijm1(:,i,j) ) > 0.0 ) then  !b=>ij   x ij=>ij-1
           new_i=new_i+1
           check_west  = .false.
           doublecheck = .true.
        endif
     endif

     if (check_west) then
        if( cross( t_im1jm1(:,i,j), t_b(:)             ) < 0.0 .and. & !T=>i-1j-1 x T     =>b
            cross( t_im1j(:,i,j),   t_b(:)             ) > 0.0 .and. & !T=>i-1j   x T     =>b
            cross( b_im1jm1(:),     im1jm1_im1j(:,i,j) ) > 0.0) then   !b=>i-1j-1 x i-1j-1=>i-1j
           new_i=new_i-1
           check_east  = .false.
           doublecheck = .true.
        endif
     endif

     if (new_i<isd .or. ied<new_i .or. new_j<jsd .or. jed<new_j) then
        if (bitwise_reproduction) then
           write (stdoutunit, '(a,i4,a,a,2(i4,a),a,4(i4,a))')                       &
                'Error on pe ',Info%pe_this, ' a blob has gone outside the halo. ',&
                'blobs (i,j) location = (',new_i,',',new_j,') ',                   &
                '(isd:ied,jsd:jed) = (',isd,':',ied,',',jsd,':',jed,')'
           call mpp_error(FATAL, &
                '==>Error in ocean_blob_util_mod (check_ijcell): ' &
                //'blob has gone outside of halo.')
        else
           b(:) = a(:)
           exit checkcell
        endif
        
     else
        i = new_i
        j = new_j
     endif

  enddo checkcell

  off(1:2) = .false.
  if (i<isc .or. iec<i) off(1) = .true.
  if (j<jsc .or. jec<j) off(2) = .true.
  lon = b(1) 
  lat = b(2)

contains 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! A nested function that does the cross product of two vectors of length!
! two.                                                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pure function cross(vec1,vec2)
    real             :: cross
    real, intent(in) :: vec1(2)
    real, intent(in) :: vec2(2)
    cross = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end function cross

end subroutine check_ijcell
! </SUBROUTINE>  NAME="check_ijcell"

!######################################################################
! <SUBROUTINE NAME="check_kcell">
!
! <DESCRIPTION>
! Searches for which (vertical) grid cell a blob resides in).
! </DESCRIPTION>
!
subroutine check_kcell(Time,Ext_mode,Thickness,geodepth,w,i,j,k,off)
  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_thickness_type),     intent(in)    :: Thickness
  real,                           intent(in)    :: geodepth
  real,                           intent(in)    :: w
  integer,                        intent(in)    :: i
  integer,                        intent(in)    :: j
  integer,                        intent(inout) :: k
  logical,                        intent(out)   :: off

  integer :: tau, increment
  logical :: foundk

  tau    = Time%tau
  off    = .false.
  foundk = .false.

  ! Check if the blob has gone out of bounds.  This can happen even if
  ! a blob does not change from one k level to another.
  if (geodepth < -Ext_mode%eta_t(i,j,tau)) then
     k=1
     off=.true.
     return
  elseif (geodepth > Grd%ht(i,j)) then
     k=Grd%kmt(i,j)
     off=.true.
     return
  endif

  ! Next, check if the blob remains in the same vertical k-level.
  if (k>1) then
     if (Thickness%geodepth_zwt(i,j,k-1) < geodepth .and. geodepth <= Thickness%geodepth_zwt(i,j,k)) foundk=.true.
  else
     if (-Ext_mode%eta_t(i,j,tau) < geodepth .and. geodepth <= Thickness%geodepth_zwt(i,j,k)) foundk=.true.
  endif

  if (foundk) then
     ! Because of the way that Thickness%geodepth_zwt is calculated, there 
     !can be a very small discrepency with Grd%ht.  So, we double check 
     ! that the cell that we are in is actually a water cell.  If it is not, 
     ! and we have made it here, we know that 
     ! Thickness%geodepth_zwt(k)<geodepth<=Grd%ht for some k<kmt. However, 
     ! there can be circumstances where, due to the discrepency between 
     ! geodepth_zwt and ht, we can have geodepth<Grd%ht and 
     ! Thickness%geodepth_zwt(kmt)<geodepth<=Thickness%geodepth_zwt(kmt+1). 
     ! In these cases we need to ensure that the blob is in the kmt cell and 
     ! not kmt+1.
     if (k>Grd%kmt(i,j)) k=Grd%kmt(i,j)
     return
  endif

  ! Guess which direction to search based on the vertical velocity
  if(w>=0.0) then !up
     increment = -1
  else !down
     increment = +1
  endif

  kcycle: do
     k = k+increment
     if (k<2 .or. Grd%kmt(i,j)<k) exit kcycle !note: we need to treat the surface grid cell differently
     if (Thickness%geodepth_zwt(i,j,k-1) < geodepth .and. geodepth <= Thickness%geodepth_zwt(i,j,k)) return
  enddo kcycle

  ! Treating the surface cell
  k=1
  if ( -Ext_mode%eta_t(i,j,tau) < geodepth .and. geodepth <= Thickness%geodepth_zwt(i,j,1)) return

  ! If we have made it this far then the search has not been successful, so, we take a brute
  ! force approach and do a full sweep of the water column from top to bottom
  column: do
     k=k+1
     if (k>Grd%kmt(i,j)) exit column !something is wrong
     if (Thickness%geodepth_zwt(i,j,k-1) < geodepth .and. geodepth <= Thickness%geodepth_zwt(i,j,k)) return
  enddo column

  ! Somtimes, there can be really small differences between Thickness%geodepth_zwt(kmt) and Grd%ht. 
  ! So, we check against Grd%ht before raising an error
  if (Thickness%geodepth_zwt(i,j,Grd%kmt(i,j))<geodepth .and. geodepth<=Grd%ht(i,j)) then
     k = Grd%kmt(i,j)
     return
  endif
     
  ! By now, we have covered all bases, so if we have made it this far, there is a problem with the logic
  ! of the search algorithm, so, we raise a fatal error
  call mpp_error(FATAL, 'ocean_blob_util_mod, check_kcell: Cannot find vertical cell for blob!')

end subroutine check_kcell
! </SUBROUTINE>  NAME="check_kcell"

!######################################################################
! <SUBROUTINE NAME="kill_blob">
!
! <DESCRIPTION>
! Kills a blob by returning all of its remaining properties to the E
! system.
! </DESCRIPTION>
!
subroutine kill_blob(Thickness, T_prog, this, i, j, k)
  type(ocean_thickness_type),   intent(inout) :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_blob_type),        pointer       :: this
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k

  integer :: n

  Thickness%blob_source(i,j) = Thickness%blob_source(i,j) &
                               + datdtime_r(i,j)*this%mass
  this%mass = 0.0

  do n=1,num_prog_tracers
     T_prog(n)%tend_blob(i,j,k) = T_prog(n)%tend_blob(i,j,k) &
                                  + this%tracer(n)*datdtime_r(i,j)
     this%tracer(n) = 0.0
  enddo

end subroutine kill_blob
! </SUBROUTINE>  NAME="kill_blob"


!######################################################################
! <SUBROUTINE NAME="free_blob_memory">
!
! <DESCRIPTION>
! Frees the heap memory taken up by a blob.
! </DESCRIPTION>
!
subroutine free_blob_memory(blob)
  type(ocean_blob_type), pointer :: blob

  deallocate(blob)
  nullify(blob)

end subroutine free_blob_memory
! </SUBROUTINE>  NAME="free_blob_memory"


!######################################################################
! <SUBROUTINE NAME="allocate_interaction_memory">
!
! <DESCRIPTION>
! Allocates the history arrays for a blob (only used when
! bitwise_reproduction=.true. in the ocean_blob_nml).
! </DESCRIPTION>
!
subroutine allocate_interaction_memory(blob, total_ns)
  type(ocean_blob_type), pointer   :: blob
  integer,              intent(in) :: total_ns

  allocate(blob%dtracer(num_prog_tracers))
  if (bitwise_reproduction) then
     allocate(blob%di(0:total_ns))
     allocate(blob%dj(0:total_ns))
     allocate(blob%dk(0:total_ns))
     allocate(blob%entrainment(1:total_ns,0:num_prog_tracers)) 
     allocate(blob%detrainment(1:total_ns,0:num_prog_tracers)) 
     allocate(blob%mass_in(1:total_ns))
     allocate(blob%mass_out(0:total_ns))
  endif
end subroutine allocate_interaction_memory
! </SUBROUTINE>  NAME="allocate_interaction_memory"

!######################################################################
! <SUBROUTINE NAME="reallocate_interaction_memory">
!
! <DESCRIPTION>
! Different blobs can have different history memory requirements.  
! When they change from one type of blob to another, we need to change
! the memory allocated to a blob to reflect the new requirements.  This
! is only necessary if bitwise_reproduction=.true. in the ocean_blob_nml.
! </DESCRIPTION>
!
subroutine reallocate_interaction_memory(blob, head, total_ns)
  type(ocean_blob_type), pointer :: blob
  type(ocean_blob_type), pointer :: head
  integer, intent(in) :: total_ns

  type(ocean_blob_type), pointer :: new_blob

  ! Allocate memory to the new blob
  allocate(new_blob)
  allocate(new_blob%tracer(num_prog_tracers))
  allocate(new_blob%field( num_prog_tracers))
  call allocate_interaction_memory(new_blob, total_ns)

  ! Copy the data from the old blob to the new blob
  new_blob%i = blob%i
  new_blob%j = blob%j
  new_blob%k = blob%k

  new_blob%m   = blob%m
  new_blob%kdw = blob%kdw
  new_blob%kup = blob%kup

  new_blob%hash   = blob%hash
  new_blob%number = blob%number

  new_blob%model_steps = blob%model_steps
  new_blob%nsteps      = blob%nsteps

  new_blob%sink = blob%sink
  new_blob%new  = blob%new

  new_blob%h1 = blob%h1
  new_blob%h2 = blob%h2

  new_blob%lat = blob%lat
  new_blob%lon = blob%lon

  new_blob%depth    = blob%depth
  new_blob%geodepth = blob%geodepth

  new_blob%st        = blob%st
  new_blob%mass      = blob%mass
  new_blob%density   = blob%density
  new_blob%densityr  = blob%densityr
  new_blob%volume    = blob%volume
  new_blob%tracer(:) = blob%tracer(:)
  new_blob%field(:)  = blob%field(:)

  new_blob%step        = blob%step
  new_blob%nfrac_steps = blob%nfrac_steps

  new_blob%v(:) = blob%v(:)

  ! Unlink the blob from the list, and link in the new blob
  new_blob%next => blob%next
  if(associated(new_blob%next)) new_blob%next%prev=>new_blob
  new_blob%prev => blob%prev
  if(associated(new_blob%prev)) then
     new_blob%prev%next=>new_blob
  else
     head=>new_blob
  endif

  ! Deallocate memory from the old blob
  call free_blob_memory(blob)

  ! Now point to the new bit of memory
  blob => new_blob

end subroutine reallocate_interaction_memory
! </SUBROUTINE>  NAME="reallocate_interaction_memory"


!######################################################################
! <SUBROUTINE NAME="interp_tcoeff">
!
! <DESCRIPTION>
! Used for the horizontal interpolation of T grid variables.  The
! routine returns coefficients required for inverse distance
! weighting (Shephard, 1968).
! </DESCRIPTION>
!
subroutine interp_tcoeff(i, j, h, lon, lat, dsq_r)
  integer,               intent(in)  :: i
  integer,               intent(in)  :: j
  real,    dimension(2), intent(in)  :: h
  real,                  intent(in)  :: lon
  real,                  intent(in)  :: lat
  real,    dimension(9), intent(out) :: dsq_r

  real,    dimension(9) :: distance
  integer :: m, iit, jjt

  distance(:) = 0.0

     ! We have special treatment for solid boundaries
     tcoeff: do m=1,Info%tidx(0,i,j)
        iit = i+Info%it(Info%tidx(m,i,j))
        jjt = j+Info%jt(Info%tidx(m,i,j))
        ! Calculate the distance from each point to the blob
        ! We do not need to include deg_to_rad because it cancels anyway
        distance(m) = onehalf                                              &
             *sqrt( abs(lon-xt(iit,jjt))**2 * (h(1)+Info%ht(iit,jjt,1))**2 &
             +abs(lat-yt(iit,jjt))**2 * (h(2)+info%ht(iit,jjt,2))**2 )
        
        if (distance(m)<epsln) then
           dsq_r(:) = 0.0
           dsq_r(m) = 1.0
           exit tcoeff
        endif
        ! Accumulate that distance 
        dsq_r(m) = 1.0/distance(m)**2
     enddo tcoeff

end subroutine interp_tcoeff
! </SUBROUTINE>  NAME="interp_tcoeff"


!######################################################################
! <SUBROUTINE NAME="interp_ucoeff">
!
! <DESCRIPTION>
! Used for the horizontal interpolation of U grid variables.  The
! routine returns coefficients required for inverse distance
! weighting (Shephard, 1968).
! </DESCRIPTION>
!
subroutine interp_ucoeff(i, j, h, lon, lat, dsq_r)
  integer,               intent(in)  :: i
  integer,               intent(in)  :: j
  real,    dimension(2), intent(in)  :: h
  real,                  intent(in)  :: lon
  real,                  intent(in)  :: lat
  real,    dimension(4), intent(out) :: dsq_r
  
  real, dimension(4) :: distance
  integer :: m, iiu, jju
  
  !Initialise some variables
  distance(:) = 0.0

     ucoeff: do m=1,Info%uidx(0,i,j)
        iiu = i+Info%iu(Info%uidx(m,i,j))
        jju = j+Info%ju(Info%uidx(m,i,j))
        ! Calculate the distance from each point to the blob
        ! We do not need to include deg_to_rad because it cancels anyway
        distance(m) = onehalf                                  &
             *sqrt( abs(lon-xu(iiu,jju))**2 * (h(1)+Info%hu(iiu,jju,1))**2  &
             +abs(lat-yu(iiu,jju))**2 * (h(2)+Info%hu(iiu,jju,2))**2 )

        if (distance(m)<epsln) then
           dsq_r(:) = 0.0
           dsq_r(m) = 1.0
           exit ucoeff
        endif
        
        ! Accumulate that distance 
        dsq_r(m) = 1.0/distance(m)**2
     enddo ucoeff

end subroutine interp_ucoeff
! </SUBROUTINE>  NAME="interp_ucoeff"


!######################################################################
! <SUBROUTINE NAME="check_cyclic">
!
! <DESCRIPTION>
! Checks and adjusts blob position and grid cell index
! for cylclic/periodic domains, as well as
! the Murray (1996) tripolar grid.
! </DESCRIPTION>
!
subroutine check_cyclic(blob, i, j, adjust_latlon)
  type(ocean_blob_type), pointer    :: blob
  integer,               intent(inout) :: i
  integer,               intent(inout) :: j
  logical,               intent(in) :: adjust_latlon

  logical :: change_ij

  ! If we have a cyclic grid and a blob goes from one PE to another
  ! we need to reset its (i,j) coordinates and its (lon,lat).
  ! Same with the tripolar grid.  Things get a bit more 
  ! complicated if we cross the arctic bipolar fold.
  change_ij = .false.
  
  if (Grd%cyclic_x) then
     if (i<isg) then
        i = nig - i
        change_ij=.true.
     elseif (i>ieg) then
        i = i - nig
        change_ij=.true.
     endif
  endif
  
  if (Grd%cyclic_y) then
     if (j<jsg) then
        j = njg - j
        change_ij=.true.
     elseif (j>jeg) then
        j = j - njg
        change_ij=.true.
     endif
  elseif(Grd%tripolar) then
     if (j>jeg) then
        ! We have crossed the bipolar Arctic fold.
        ! We reset the i vlue to correspond to the
        ! opposing i value, and reset the j value 
        ! to be jeg.
        ! See figure 4.6 of Griffies et al (2004)
        i = nip1 - i
        j = jeg
        
        ! If we cross the geographic north pole.
        if (blob%lat>90.) then
           blob%lat  = 180-blob%lat
           blob%v(2) = -blob%v(2)
        endif
        
        if (adjust_latlon) then
           ! Handle the x-cyclic bit, if required
           if    (blob%lon<Info%minlon(j)) then
              blob%lon = Info%maxlon(j) - modulo(Info%minlon(j),blob%lon)
           elseif(blob%lon>Info%maxlon(j)) then
              blob%lon = Info%minlon(j) + modulo(blob%lon,Info%maxlon(j))
           endif
        endif
     endif
  endif

  if (adjust_latlon) then
     
     if (change_ij) then
        if    (blob%lon<Info%minlon(j)) then
           blob%lon = Info%maxlon(j) - modulo(Info%minlon(j),blob%lon)
        elseif(blob%lon>Info%maxlon(j)) then
           blob%lon = Info%minlon(j) + modulo(blob%lon,Info%maxlon(j))
        endif
        
        if    (blob%lat<Info%minlat(i)) then
           blob%lat = Info%maxlat(i) - modulo(Info%minlat(i),blob%lat)
        elseif(blob%lat>Info%maxlat(i)) then
           blob%lat = Info%minlat(i) + modulo(blob%lat,Info%maxlat(i))
        endif
     
     endif
  endif


end subroutine check_cyclic
! </SUBROUTINE>  NAME="check_cyclic"

end module ocean_blob_util_mod
