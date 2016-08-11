
module bgrid_prog_var_mod

!-----------------------------------------------------------------------
!
!     allocates memory for B-grid core prognostics variables
!    and provides several routines for handling these variables
!
!     the module has routines for:
!       1) setting up the prognostic variables in a data structure
!       2) performing simple arithmetic with this data structure
!       3) applying forward explicit time differencing 
!       4) reading and writing prognostic variable restart files
!
!-----------------------------------------------------------------------

use      bgrid_horiz_mod, only: horiz_grid_type
use       bgrid_vert_mod, only: vert_grid_type
use      bgrid_masks_mod, only: grid_mask_type
use       bgrid_halo_mod, only: update_halo, TEMP, UWND, VWND
use bgrid_cold_start_mod, only: cold_start_resol, cold_start
use              fms_mod, only: file_exist, open_restart_file, mpp_error, &
                                FATAL, close_file, mpp_pe, mpp_root_pe,   &
                                set_domain, read_data, write_data,        &
                                write_version_number, nullify_domain,     &
                                field_size, NOTE, mpp_chksum,             &
                                error_mesg, FATAL, stdout, stdlog,        &
                                uppercase
use    field_manager_mod, only: MODEL_ATMOS
use   tracer_manager_mod, only: get_tracer_names, set_tracer_profile

use         platform_mod, only: I8_KIND

implicit none
private

public :: prog_var_type, prog_var_init, var_init,  &
          prog_var_time_diff, prog_var_times_scalar, &
          prog_var_equals_scalar
public :: open_prog_var_file, read_prog_var, write_prog_var

!-----------------------------------------------------------------------
! type(prog_var_type)
!     data structure that contains all prognostic fields and tracers
!
!     nlon = size of the global x-axis compute grid
!             (i.e., number of longitude points)
!     nlat = size of the global y-axis compute grid
!             (i.e., number of latitude points)
!     nlev = number of vertical levels
!     ntrace = number of tracers (in the current data structure)
!
!     ilb, iub = lower and upper index bounds of x-axis (data domain)
!     jlb, jub = lower and upper index bounds of y-axis (data domain)
!     klb, kub = lower and upper index bounds of k-axis
!
! arrays with dimension(ilb:iub,jlb:jub)
!     ps   = surface pressure
!     pssl = surface pressure adjusted to eta=1. (for eta coordinate)
!
! arrays with dimension(ilb:iub,jlb:jub,klb:kub)
!     u    = zonal wind component
!     v    = meridional wind component
!     t    = temperature
!
! arrays with dimension(ilb:iub,jlb:jub,klb:kub,1:ntrace)
!     r    = arbitrary number of tracers
!
!-----------------------------------------------------------------------

type prog_var_type
     integer       :: nlon, nlat, nlev, ntrace
     integer       :: ilb, iub, jlb, jub, klb, kub
     real, pointer :: ps(:,:) =>NULL(), &
                      pssl(:,:) =>NULL()
     real, pointer :: u(:,:,:) =>NULL(), &
                      v(:,:,:) =>NULL(), &
                      t(:,:,:) =>NULL(), &
                      r(:,:,:,:) =>NULL()
end type prog_var_type

! overloaded interface for initializing real model arrays
interface var_init
    module procedure var_init_type_4d, var_init_bound_4d, &
                     var_init_type_3d, var_init_bound_3d, &
                     var_init_type_2d, var_init_bound_2d
end interface

!-----------------------------------------------------------------------
! private data

logical :: do_log = .true.
character(len=128) :: version='$Id: bgrid_prog_var.F90,v 13.0 2006/03/28 21:05:18 fms Exp $'
character(len=128) :: tagname='$Name: tikal $'

integer :: unit_in
character(len=128) :: directory_in
logical :: read_pssl
logical :: old_restart_format, do_cold_start

character(len=64) :: res_file_name = 'bgrid_prog_var.res'

character(len=80) :: restart_format = &
              'bgrid grid atmospheric dynamical core: restart format 05'

contains

!#######################################################################
! creates a prog_var_type variable

 subroutine prog_var_init (Hgrid, nlev, ntrs, Vars)

  type(horiz_grid_type), intent(in)    :: Hgrid
  integer,               intent(in)    :: nlev, ntrs
  type(prog_var_type)  , intent(inout) :: Vars
!-----------------------------------------------------------------------
! write version info to logfile
  if (do_log) then
    call write_version_number(version, tagname)
    do_log = .false.
  endif

! all arrays have the same horizontal dimensions regardless of
! whether the field is on the temperature or velocity grid

    Vars % ilb = Hgrid % ilb
    Vars % iub = Hgrid % iub
    Vars % jlb = Hgrid % jlb
    Vars % jub = Hgrid % jub
    Vars % klb = 1
    Vars % kub = nlev

    Vars % nlon = Hgrid % nlon
    Vars % nlat = Hgrid % nlat
    Vars % nlev = nlev
    Vars % ntrace = ntrs

    Vars % ps   => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % pssl => var_init_bound_2d (Vars % ilb, Vars % iub, &
                                      Vars % jlb, Vars % jub)

    Vars % u => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % v => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % t => var_init_bound_3d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev)
    Vars % r => var_init_bound_4d (Vars % ilb, Vars % iub, &
                                   Vars % jlb, Vars % jub, &
                                   nlev, ntrs)

 end subroutine prog_var_init

!#######################################################################
!##### overloaded functions that allocate a single real variable #######
!#######################################################################
!
!      variables must be declard as pointers
!      real, pointer :: field(:,:,:)
!      field => var_init (Hgrid,nlev)
!
!#######################################################################

 function var_init_bound_2d (ilb, iub, jlb, jub) result (var)

  integer, intent(in)           :: ilb, iub, jlb, jub
  real, dimension(:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub) )
    var = 0.0

 end function var_init_bound_2d

!#######################################################################

 function var_init_type_2d (Hgrid) result (var)

  type(horiz_grid_type), intent(in) :: Hgrid
  real, dimension(:,:), pointer     :: var

    var => var_init_bound_2d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub)

 end function var_init_type_2d

!#######################################################################

 function var_init_bound_3d (ilb, iub, jlb, jub, kdim) result (var)

  integer, intent(in)             :: ilb, iub, jlb, jub, kdim
  real, dimension(:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim) )
    var = 0.0

 end function var_init_bound_3d

!#######################################################################

 function var_init_type_3d (Hgrid, kdim) result (var)

  type(horiz_grid_type), intent(in) :: Hgrid
  integer, intent(in)               :: kdim
  real, dimension(:,:,:), pointer   :: var

    var => var_init_bound_3d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim)

 end function var_init_type_3d

!#######################################################################

 function var_init_bound_4d (ilb, iub, jlb, jub, kdim, ntrace) result (var)

  integer, intent(in)               :: ilb, iub, jlb, jub, kdim, ntrace
  real, dimension(:,:,:,:), pointer :: var

    allocate ( var (ilb:iub, jlb:jub, 1:kdim, 1:ntrace) )
    var = 0.0

 end function var_init_bound_4d

!#######################################################################

 function var_init_type_4d (Hgrid, kdim, ntrace) result (var)

  type(horiz_grid_type), intent(in)   :: Hgrid
  integer, intent(in)                 :: kdim, ntrace
  real, dimension(:,:,:,:), pointer   :: var

    var => var_init_bound_4d (Hgrid % ilb, Hgrid % iub, &
                              Hgrid % jlb, Hgrid % jub, kdim, ntrace)

 end function var_init_type_4d

!#######################################################################
!#######################################################################
! sets all prognostic variables to a scalar

 subroutine prog_var_equals_scalar (Var, scalar)

  type(prog_var_type), intent(inout) :: Var
  real               , intent(in)    :: scalar

     Var % u    = scalar
     Var % v    = scalar
     Var % t    = scalar
     Var % r    = scalar
     Var % ps   = scalar
     Var % pssl = scalar

 end subroutine prog_var_equals_scalar

!#######################################################################
! multiplies all prognostic variables by a scalar

 subroutine prog_var_times_scalar (Var, scalar)

  type(prog_var_type), intent(inout) :: Var
  real               , intent(in)    :: scalar

     Var % u    = Var % u    * scalar
     Var % v    = Var % v    * scalar
     Var % t    = Var % t    * scalar
     Var % r    = Var % r    * scalar
     Var % ps   = Var % ps   * scalar
     Var % pssl = Var % pssl * scalar

 end subroutine prog_var_times_scalar

!#######################################################################
! performs time differencing on all prognostic variables
!             Var = Var + dt * Var_dt
! all tracers are used unless argument nt is supplied

 subroutine prog_var_time_diff (dt, Masks, Var_dt, Var, nt)

  real,                 intent(in)    :: dt
  type(grid_mask_type), intent(in)    :: Masks
  type(prog_var_type),  intent(inout) :: Var_dt, Var
  integer, optional,    intent(in)    :: nt

  integer :: ntp, i, j, k

!----- explicit differencing with two time levels -----

   ntp = Var_dt % ntrace
   if (present(nt)) ntp = min(Var_dt%ntrace, nt)

   Var % ps   = Var % ps   + dt * Var_dt % ps
   Var % pssl = Var % pssl + dt * Var_dt % pssl

   Var % u = Var % u + dt * Var_dt % u
   Var % v = Var % v + dt * Var_dt % v
   Var % t = Var % t + dt * Var_dt % t

   Var % r(:,:,:,1:ntp) = Var % r(:,:,:,1:ntp) + &
                                dt * Var_dt % r(:,:,:,1:ntp)

  ! apply mask for step-mountain
   if (.not.Masks%sigma) then
       do k = Var%klb, Var%kub
       do j = Var%jlb, Var%jub
       do i = Var%ilb, Var%iub
          if (Masks%Vel%mask(i,j,k) < .01) then
             Var % u (i,j,k) = 0.
             Var % v (i,j,k) = 0.
          endif
        ! temp & tracers
          if (Masks%Tmp%mask(i,j,k) < .01) then
             Var % t (i,j,k) = 0.
             Var % r (i,j,k,1:ntp) = 0.
          endif
       enddo
       enddo
       enddo
   endif

!----- zero out tendencies -----

   Var_dt % ps   = 0.0
   Var_dt % pssl = 0.0

   Var_dt % u = 0.0
   Var_dt % v = 0.0
   Var_dt % t = 0.0

   Var_dt % r(:,:,:,1:ntp) = 0.0


 end subroutine prog_var_time_diff

!#######################################################################
!########## routines for reading and writing restart files #############
!#######################################################################

subroutine open_prog_var_file ( ix, jx, kx, dir )

!-----------------------------------------------------------------------
! This routine can open for reading either a NATIVE or NETCDF restart
! file.  Only the model resolution is returned. A subsequent call to 
! read_prog_var is needed to read the data.
! The restart file is called "INPUT/bgrid_prog_var.res".
! For NetCDF files a ".nc" suffix is added.
!
!   ix, jx, kx = global resolution read from the restart file
!   dir        = directory where input restart files reside
!                  for the current directory use: "" or "."
!                  default: "INPUT"
!-----------------------------------------------------------------------
 integer,          intent(out)          :: ix, jx, kx
 character(len=*), intent(in), optional :: dir

 integer :: ic, vers, day, sec, ntsd, nt, ntp, siz(4)
 character(len=128) :: filename_in
 character(len=80)  :: control
 character(len=2)   :: avers

! write version info to logfile
  if (do_log) then
    call write_version_number(version, tagname)
    do_log = .false.
  endif

! set-up restart directory and file name
  directory_in = 'INPUT/'
  if (present(dir)) then
      if (len_trim(dir) > 0) then
          directory_in = trim(dir)//"/"
      else
          directory_in = '' ! null string
      endif
  endif

  filename_in = trim(directory_in)//trim(res_file_name)

! when restart file does not exist
! set up simple initial conditions

  if (file_exist(trim(filename_in)//'.nc') ) then
       call field_size (trim(filename_in)//'.nc', 't', siz )
       ix = siz(1); jx = siz(2); kx = siz(3)
       if (min(ix,jx,kx) <= 0) call mpp_error ('bgrid_prog_var_mod', &
                        'problem reading field size; siz <= 0', FATAL)
       old_restart_format = .false.
       do_cold_start      = .false.
       if (mpp_pe() == mpp_root_pe()) call mpp_error ('bgrid_prog_var_mod', &
                              'Reading NetCDF formatted restart file.', NOTE)

! old native format restart files
  else if (file_exist(trim(filename_in)) ) then

      ! open restart file and get restart version number
      ! if control record cannot be read then file uses older format (1 or 2)
        unit_in = open_restart_file ( trim(filename_in), 'read' )
        old_restart_format = .true.
        do_cold_start      = .false.
        if (mpp_pe() == mpp_root_pe()) call mpp_error ('bgrid_prog_var_mod', &
                               'Reading native formatted restart file.', NOTE)
        read  (unit_in,err=2)  control 

      ! extract version number
        ic = index(control,'restart format ')
        if (ic == 0) call mpp_error ('bgrid_prog_var_mod', &
                     'problem extracting restart version number', FATAL)
        avers = control(ic+15:ic+16)
        read (avers,'(i2.2)') vers
        go to 3 

      ! read version number from old format (first rewind file)
      2 rewind (unit_in)
        read   (unit_in) vers
        write  (avers,'(i2.2)') vers

      ! read first (non-control) record of restart file
      ! note: ntsd,day,sec are no longer used and
      !       number of time levels (nvar) is not read or used
      3 continue
        select case (vers)
          case (1:2)
               read  (unit_in) ntsd, day, sec, ix, jx, kx, nt, ntp
               read_pssl = .false.
          case (3)
               read  (unit_in) ntsd, day, sec, ix, jx, kx, nt, ntp
               read_pssl = .true.
          case (4)
               read  (unit_in) ix, jx, kx, nt, ntp
               read_pssl = .true.
          case (5)
               read  (unit_in) ix, jx, kx
               read_pssl = .true.
          case default
             call mpp_error ('bgrid_prog_var_mod', &
                             'cannot not read old restart version '//avers, FATAL)
        end select

! no restart --> self start
  else
       call cold_start_resol ( ix, jx, kx )
       do_cold_start = .true.
       read_pssl = .false.
       return
  endif

end subroutine open_prog_var_file

!#######################################################################

subroutine read_prog_var (Hgrid, Var, eta, peta, fis, res)

!-----------------------------------------------------------------------
! This routine can read either a NATIVE or NETCDF restart file.
!
! Hgrid  = horizontal grid constants
! Var    = prognostic variables
! eta    = sigma/eta/bk values at model layer interfaces (half levels)
! peta   = reference pressures (pk) at model layer interfaces
! fis    = geopotential height of the surface
! res    = reciprocal of eta at the surface
!-----------------------------------------------------------------------
 type(horiz_grid_type), intent(inout) :: Hgrid
 type  (prog_var_type), intent(inout) :: Var
 real, intent(out), dimension(:)      :: eta, peta
 real, intent(out), dimension(Hgrid%ilb:Hgrid%iub, &
                              Hgrid%jlb:Hgrid%jub) :: fis, res
 integer :: n, unit, siz(4)
 integer :: isd,ied,vsd,ved
 character(len=64)  :: tr_name
 character(len=128) :: filename_in
 real :: tr_surf, tr_mult
 logical :: found

! set-up file name
  filename_in = trim(directory_in)//trim(res_file_name)
  if (.not.old_restart_format) filename_in = trim(filename_in)//'.nc'

  if (do_cold_start) then

     ! set up simple initial conditions
     ! when restart file does not exist

       call cold_start ( Hgrid, eta, peta, fis, res, Var%ps, Var%pssl, &
                                           Var%u, Var%v, Var%t )
  else

     ! must pass fields to read_data on data domain
     ! mass fields are on data domain
     ! set up indexing for velocity fields on data domain

       isd = Hgrid%Vel%isd;  ied = Hgrid%Vel%ied
       vsd = Hgrid%Vel%jsd;  ved = Hgrid%Vel%jed

     ! read non-distributed data from root pe (vertical coordinate info)
       if (old_restart_format) then
           read (unit_in) eta, peta
       else
           call nullify_domain ()
           call read_data ( filename_in,  'eta',  eta )
           call read_data ( filename_in, 'peta', peta )
       endif

     !  --- read variables ---

     ! initialize domain for temperature grid 
     ! read surf pres, topog, and more

       if (old_restart_format) then
           call set_domain ( Hgrid%Tmp%Domain )
           call read_data ( unit_in, Var%ps  )
           if (read_pssl) &
           call read_data ( unit_in, Var%pssl)
           call read_data ( unit_in,     res )
           call read_data ( unit_in,     fis )
       else
           call read_data ( filename_in, 'ps'  , Var%ps  , Hgrid%Tmp%Domain )
           if (read_pssl) &
           call read_data ( filename_in, 'pssl', Var%pssl, Hgrid%Tmp%Domain )
           call read_data ( filename_in, 'res' ,     res , Hgrid%Tmp%Domain )
           call read_data ( filename_in, 'fis' ,     fis , Hgrid%Tmp%Domain )
       endif
     
     ! initialize domain for velocity grid 
     ! read u and v wind components

     ! pass velocity fields on data domain
       if (old_restart_format) then
           call set_domain ( Hgrid%Vel%Domain )
           call read_data ( unit_in, Var%u(isd:ied,vsd:ved,:) )
           call read_data ( unit_in, Var%v(isd:ied,vsd:ved,:) )
       else
           call read_data ( filename_in, 'u', Var%u(isd:ied,vsd:ved,:), Hgrid%Vel%Domain )
           call read_data ( filename_in, 'v', Var%v(isd:ied,vsd:ved,:), Hgrid%Vel%Domain )
       endif

     ! re-initialize domain for temperature grid 
     ! read temperature and tracers
       if (old_restart_format) then
           call set_domain ( Hgrid%Tmp%Domain )
           call read_data ( unit_in, Var%t )
           call close_file (unit_in) ! done reading B-grid dynamics restart
       else
           call read_data ( filename_in, 't', Var%t, Hgrid%Tmp%Domain )
       endif

  endif

! read single tracer restart file
  
  filename_in = trim(directory_in)//'atmos_tracers.res.nc'

  if (file_exist(trim(filename_in))) then
      do n = 1, Var%ntrace
         call get_tracer_names ( MODEL_ATMOS, n, tr_name )
         call field_size (trim(filename_in), tr_name, siz, field_found=found)
         if (found) then
             if (siz(1) == Hgrid%nlon .and. siz(2) == Hgrid%nlat .and. siz(3) == size(eta(:))-1) then
                call read_data (trim(filename_in), tr_name, Var%r(:,:,:,n), Hgrid%Tmp%Domain)
             else
                call error_mesg ('bgrid_prog_var_mod','can not read tracers with wrong size', FATAL)
             endif
         else
          ! initialize new tracers (apply surface value only)
            call set_tracer_profile ( MODEL_ATMOS, n, Var%r(:,:,:,n) )
         endif
      enddo
  else

    ! old format
    ! read separate tracer restart files
      do n = 1, Var%ntrace
         call get_tracer_names ( MODEL_ATMOS, n, tr_name )
         filename_in = trim(directory_in)//'tracer_'//trim(tr_name)//'.res'

         if (file_exist(trim(filename_in))) then
             unit = open_restart_file( trim(filename_in), 'read' )
             call set_domain ( Hgrid%Tmp%Domain )
             call read_data ( unit, Var%r(:,:,:,n) )
             call close_file (unit)
         else
           ! initialize new tracers (apply surface value only)
             call set_tracer_profile ( MODEL_ATMOS, n, Var%r(:,:,:,n) )
         endif
      enddo

  endif

! update all boundaries for restart variables

  call update_halo (Hgrid, TEMP, res)
  call update_halo (Hgrid, TEMP, fis)
  call update_halo (Hgrid, TEMP, Var%ps)
  if (read_pssl) &
  call update_halo (Hgrid, TEMP, Var%pssl)
  call update_halo (Hgrid, TEMP, Var%t)
  call update_halo (Hgrid, TEMP, Var%r)
  call update_halo (Hgrid, UWND, Var%u)
  call update_halo (Hgrid, VWND, Var%v)

! for old restart formats, initialize pssl
  if (.not.read_pssl) Var%pssl = Var%ps * res

! check sum input data
  call print_check_sum ('Check sums for B-grid input data:', &
                        Hgrid, Var, eta, peta, fis, res)

end subroutine read_prog_var

!#######################################################################
! writes a B-grid core restart file named "RESTART/bgrid_prog_var.res"
! For NetCDF files a ".nc" suffix is added.

 subroutine write_prog_var (Var, Hgrid, Vgrid, fis, res, dir, format)

!-----------------------------------------------------------------------
! Var    = prognostic variables
! Hgrid  = horizontal grid constants
! Vgrid  = vertical grid constants
! fis    = geopotential height of the surface
! res    = reciprocal of eta at the surface
! dir    = directory where output restart files will be written
!            for the current directory use: "" or "."
!            default: "RESTART"
! format = file format, either: NATIVE or NETCDF
!-----------------------------------------------------------------------
 type  (prog_var_type), intent(in) :: Var
 type(horiz_grid_type), intent(in) :: Hgrid
 type (vert_grid_type), intent(in) :: Vgrid
 real, intent(in), dimension(Hgrid%ilb:Hgrid%iub, &
                             Hgrid%jlb:Hgrid%jub) :: fis, res
 character(len=*), intent(in), optional :: dir, format

 integer :: n, unit
 integer :: isd,ied, vsd,ved
 character(len=8)  :: oform
 character(len=64) :: tr_name
 character(len=128) :: directory, restart_name
 logical :: use_old_format

!-----------------------------------------------------------------------
! output format
  oform = 'NETCDF'; if (present(format)) oform = format
  if (trim(uppercase(oform)) == 'NATIVE') then
      use_old_format = .true.
  else if (trim(uppercase(oform)) == 'NETCDF') then
      use_old_format = .false.
  else
      call error_mesg ('write_prog_var', &
                       'invalid value for optional argument format', FATAL)
  endif

! check sum output data
  call print_check_sum ('Check sums for B-grid output data:', &
                        Hgrid, Var, Vgrid%eta, Vgrid%peta, fis, res)

! set-up restart file name
  directory = 'RESTART/'
  if (present(dir)) then
      if (len_trim(dir) > 0) then
          directory = trim(dir)//'/'  ! append /
      else
          directory = ''              ! cwd
      endif
  endif

  restart_name = trim(directory)//trim(res_file_name)

! open output restart file
  if (use_old_format) then
      unit = open_restart_file ( trim(restart_name), 'write' )
    ! write non-distributed data from root pe
      if ( mpp_pe() == mpp_root_pe() ) then
           write (unit)  restart_format
           write (unit)  Hgrid%nlon, Hgrid%nlat, Vgrid%nlev
           write (unit)  Vgrid%eta, Vgrid%peta ! vertical coordinate info
      endif
  else
! netcdf format
      restart_name = trim(restart_name)//'.nc' ! add netcdf file suffix
      call nullify_domain ()
      call write_data ( restart_name,  'eta', Vgrid%eta  )
      call write_data ( restart_name, 'peta', Vgrid%peta )
  endif

! must pass fields to write_data on data domain
! mass fields are on data domain
! set up indexing for velocity fields on data domain

  isd = Hgrid%Vel%isd;  ied = Hgrid%Vel%ied
  vsd = Hgrid%Vel%jsd;  ved = Hgrid%Vel%jed

! output 2d-fields on mass grid (surf pres, topog)
  if (use_old_format) then
      call set_domain ( Hgrid%Tmp%Domain )
      call write_data ( unit, Var%ps  )
      call write_data ( unit, Var%pssl)
      call write_data ( unit,     res )   
      call write_data ( unit,     fis )   
  else
      call write_data ( restart_name, 'ps',   Var%ps  , Hgrid%Tmp%Domain )
      call write_data ( restart_name, 'pssl', Var%pssl, Hgrid%Tmp%Domain )
      call write_data ( restart_name, 'res',  res     , Hgrid%Tmp%Domain )
      call write_data ( restart_name, 'fis',  fis     , Hgrid%Tmp%Domain )
  endif

! output 3d-fields on velocity grid (u,v)
! pass velocity fields on data domain
  if (use_old_format) then
      call set_domain ( Hgrid%Vel%Domain )
      call write_data ( unit, Var%u(isd:ied,vsd:ved,:) )
      call write_data ( unit, Var%v(isd:ied,vsd:ved,:) )
  else
      call write_data ( restart_name, 'u', Var%u(isd:ied,vsd:ved,:), Hgrid%Vel%Domain )
      call write_data ( restart_name, 'v', Var%v(isd:ied,vsd:ved,:), Hgrid%Vel%Domain )
  endif

! output 3d temperature field
  if (use_old_format) then
    ! re-initialize domain for temperature grid (save temp and tracers)
      call set_domain ( Hgrid%Tmp%Domain )
      call write_data ( unit, Var%t   )
    ! done writing B-grid dynamics restart
      call close_file ( unit )
  else
      call write_data ( restart_name, 't', Var%t, Hgrid%Tmp%Domain )
  endif

! write tracer restart file(s)
  do n = 1, Var%ntrace
     call get_tracer_names ( MODEL_ATMOS, n, tr_name )
     if (use_old_format) then
       ! write separate files when using old format
         restart_name = trim(directory)//'tracer_'//trim(tr_name)//'.res'
         unit = open_restart_file( trim(restart_name), 'write' )
         call write_data ( unit, Var%r(:,:,:,n) )
         call close_file ( unit )
     else
       ! write single netcdf tracer file
         restart_name = trim(directory)//'atmos_tracers.res.nc'
         call write_data ( trim(restart_name), trim(tr_name), &
                           Var%r(:,:,:,n), Hgrid%Tmp%Domain )
     endif
  enddo

 end subroutine write_prog_var

!#######################################################################

subroutine print_check_sum (label, Hgrid, Var, eta, peta, fis, res, unit)

 character(len=*),      intent(in) :: label
 type(horiz_grid_type), intent(in) :: Hgrid
 type  (prog_var_type), intent(in) :: Var
 real, intent(in), dimension(:)    :: eta, peta
 real, intent(in), dimension(Hgrid%ilb:Hgrid%iub, &
                             Hgrid%jlb:Hgrid%jub) :: fis, res
 integer, intent(in), optional :: unit
 integer :: n, chksum_unit
 integer :: is, ie, hs, he, vs, ve
 integer(I8_KIND) :: zsum

  chksum_unit = stdout()
  if (present(unit)) chksum_unit = unit
 
  is = Hgrid%Tmp%is;  ie = Hgrid%Tmp%ie
  hs = Hgrid%Tmp%js;  he = Hgrid%Tmp%je
  vs = Hgrid%Vel%js;  ve = Hgrid%Vel%je

  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,'(/a)') trim(label)
! sum all restart variables (except eta,peta)
  zsum = mpp_chksum(fis (is:ie,hs:he))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(fis)',zsum
  zsum = mpp_chksum(res (is:ie,hs:he))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(res)',zsum
  zsum = mpp_chksum(Var%ps  (is:ie,hs:he))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(ps)',zsum
  zsum = mpp_chksum(Var%pssl(is:ie,hs:he))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(pssl)',zsum
  zsum = mpp_chksum(Var%t(is:ie,hs:he,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(t)',zsum
  zsum = mpp_chksum(Var%u(is:ie,vs:ve,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(u)',zsum
  zsum = mpp_chksum(Var%v(is:ie,vs:ve,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(v)',zsum
  zsum = 0
! sum of all tracers
  do n = 1, Var%ntrace
    zsum = zsum + mpp_chksum(Var%r(is:ie,hs:he,:,n))
  enddo
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(r)',zsum
10 format ('chksum',a6,' = ',z16)

end subroutine print_check_sum

!#######################################################################

end module bgrid_prog_var_mod

