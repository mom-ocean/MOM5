module ocean_operators_mod
!   
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S.M. Griffies
!</CONTACT>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> A. Rosati
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Zhi Liang
!</REVIEWER>
!
!<REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Alexander Pletzer
!</REVIEWER>
!
!<OVERVIEW>
! Operators for MOM
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes discrete operators used by MOM.  
!</DESCRIPTION>
!
! <INFO>
!
!<NOTE>
! All operators will be replaced by generic forms when Fortran
! can properly support functions of allocatable arrays.  The 
! problems presently with this replacement are are as follows:
! Allocatable arrays cannot be inside of derived types.
! Only pointers to allocatable arrays can be inside derived types.
! Supposedly the former will be allowed in Fortran 95
! Also, functions cannot be typed as a derived type without conflicts
! which preclude using function as general operators operating on
! derived types.
! </NOTE>
!
! <NOTE>
! Mnemonics for simple operators
!
! 1st letter (direction of operation)
! <BR/>
! F => Forward direction with respect to the index.
! <BR/>
! B => Backward direction with respect to the index.
!
! 2nd letter (operation)
! <BR/>
! D => Derivative
! <BR/> 
! A => Average
! <BR/>
! M => Minimum
!
! 3rd letter (axis)
!
! X => along the X axis
! <BR/>
! Y => along the Y axis
! <BR/>
! Z => along the Z axis
!
! 4th letter (placement of quantity being operated on)
!
! E => East face
! <BR/>
! N => North face
! <BR/>
! B => Bottom face
! <BR/>
! P => Point (grid point within cell)
!
! 5th letter (type of grid cell)
!
! U => U-cell
! <BR/>
! T => T-cell
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_operators_nml">
!  <DATA NAME="use_legacy_DIV_UD" TYPE="logical">
!  Set use_legacy_DIV_UD=.true. to reproduce Riga results for 
!  DIV_UD on Bgrid. For the case that the model grid is tripolar grid, 
!  when barotropic_halo > 1 in ocean_barotropic.F90, then 
!  we must set use_legacy_DIV_UD=.false., since will not reproduce between
!  different number of processors if set use_legacy_DIV_UD=.true.
!
!  Tests indicate that with wider barotropic halos, there are
!  some performance enhancements for use_legacy_DIV_UD=.false.  
!  Hence, the default is use_legacy_DIV_UD=.false.
!
!  For the case that the model grid is regular lat-lon grid,
!  use_legacy_DIV_UD could be set to .true. or .false. for 
!  any positive value of barotropic_halo.
!
!  Note that the only difference between the new and old DIV_UD
!  is order of operations induced by parentheses, which occurs in the  
!  tripolar fold region in the Arctic: 
!  old: DIV_UD(i,j) = (uh_bay  - uhim_bay  + vh_bax  - vhjm_bax)*datr_bt(i,j)
!  new: DIV_UD(i,j) = ((uh_bay - uhim_bay) + (vh_bax - vhjm_bax))*datr_bt(i,j) 
!  </DATA> 
!</NAMELIST>

  use mpp_domains_mod,   only: mpp_update_domains, CGRID_NE, BGRID_NE
  use mpp_domains_mod,   only: mpp_get_data_domain, mpp_get_compute_domain
  use mpp_domains_mod,   only: EAST, NORTH, CORNER
  use mpp_mod,           only: mpp_error, FATAL, stdlog, stdout, input_nml_file
  use fms_mod,           only: open_namelist_file, check_nml_error, close_file, file_exist

  use ocean_domains_mod,    only: get_local_indices, get_halo_sizes    
  use ocean_domains_mod,    only: set_ocean_domain
  use ocean_parameters_mod, only: onefourth, oneeigth, onesixteenth
  use ocean_parameters_mod, only: MOM_BGRID, MOM_CGRID 
  use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, ocean_thickness_type 

  implicit none

  private

  public FAX, FAY ! forward averaging operators
  public BAX, BAY ! backward averaging operators
  public FMX, FMY ! forward minimum operators

  public FDZ_T    ! forward derivative operator in vertical  

  ! derivative operators along constant depth surfaces 
  ! account taken of general sloping of the k-surfaces 
  public FDX_ZT, FDY_ZT  ! forward 

  ! derivative operators along constant k-levels 
  ! equaly, derivatives along constant vertical coordinate surfaces 
  public FDX_T,  FDY_T   ! forward for fields at T-points
  public FDX_U,  FDY_U   ! forward for fields at U-points
  public FDX_NT, FDY_ET  ! forward for fields at T-cell faces 
  public BDX_ET, BDY_NT  ! backward for fields at T-cell faces
  public BDX_EU, BDY_NU  ! backward for fields at U-cell faces 

  ! Not-so-simple operators (no mnemonics) 
  public REMAP_NT_TO_NU    ! remap thickness weighted advective velocity on north face of T-cells to U-cells 
  public REMAP_ET_TO_EU    ! remap thickness weighted advective velocity on east face of T-cells to U-cells 
  public REMAP_BT_TO_BU    ! remap T-cell thickness or vertical velocity on base of T-cells to U-cells
  public DIV_UD            ! divergence of thickness weighted barotropic velocity
  public GRAD_BAROTROPIC_P ! Gradient of surface pressure or bottom Montgomery potential
  public S2D               ! Smoothing operator: 2D version of [1/4, 1/2, 1/4] filter
  public LAP_T             ! Lateral Laplacian of T-cell fields weighted by a diffusivity 
  public set_barotropic_domain
  public get_use_legacy_DIV_UD

  !namelist
  logical :: use_legacy_DIV_UD = .false. 
  namelist /ocean_operators_nml/ use_legacy_DIV_UD

  integer :: barotropic_halo, isd_bt, ied_bt, jsd_bt, jed_bt
  real, dimension(:,:),   allocatable :: dxtn_bt,  dyte_bt
  real, dimension(:,:),   allocatable :: dxte_bt,  dytn_bt
  real, dimension(:,:),   allocatable :: dxter_bt, dytnr_bt
  real, dimension(:,:),   allocatable :: dxu_bt,   dyu_bt
  real, dimension(:,:),   allocatable :: dxur_bt,  dyur_bt
  real, dimension(:,:),   allocatable :: dat_bt,   datr_bt
  real, dimension(:,:),   allocatable :: umask_bt
  real, dimension(:,:,:), allocatable :: tmasken_bt
  type(ocean_domain_type), pointer    :: Dom_bt => NULL()

  character (len=128) :: version = &
       '$Id: ocean_operators.F90,v 20.0 2013/12/14 00:10:53 fms Exp $'

  character (len=128) :: tagname = &
       '$Name: tikal $'

  type(ocean_grid_type), pointer      :: Grd =>NULL()
  type(ocean_domain_type), pointer    :: Dom =>NULL()
  type(ocean_thickness_type), pointer :: Thk =>NULL()
  type(ocean_domain_type), save       :: Dom_flux

  !for Bgrid or Cgrid
  integer :: horz_grid

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS
  real, dimension(isd:ied,jsd:jed,nk) :: fmx_tmask
  real, dimension(isd:ied,jsd:jed,nk) :: fmy_tmask
  integer, parameter                  :: halo = 1
#else
  real, dimension(:,:,:), allocatable :: fmx_tmask
  real, dimension(:,:,:), allocatable :: fmy_tmask
  integer                             :: halo
#endif

  logical :: module_is_initialized = .FALSE.

  public ocean_operators_init

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_operators_init">
!
! <DESCRIPTION>
! Initialize the operator module
! </DESCRIPTION>
!
subroutine ocean_operators_init(Grid, Domain, Thickness, hor_grid)

    type(ocean_grid_type), intent(in), target      :: Grid
    type(ocean_domain_type), intent(in), target    :: Domain
    type(ocean_thickness_type), intent(in), target :: Thickness
    integer, intent(in)                            :: hor_grid 

    integer :: xhalo
    integer :: yhalo
    integer :: k
    integer :: stdoutunit, stdlogunit
    integer :: ioun, io_status, ierr

    stdoutunit=stdout()
    stdlogunit = stdlog()

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ocean_operators_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'ocean_operators_nml')
#else
    ioun = open_namelist_file()
    read  (ioun, ocean_operators_nml,iostat=io_status)
    ierr = check_nml_error(io_status, 'ocean_operators_nml')
    call close_file(ioun)
#endif
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_operators_nml)
    write (stdlogunit, ocean_operators_nml)

    Grd => Grid
    Dom => Domain
    Thk => Thickness
    horz_grid = hor_grid 

    call set_ocean_domain(Dom_flux, Grid, name='horz diff flux', maskmap=Dom%maskmap)

#ifndef MOM_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    call get_halo_sizes(Domain,xhalo, yhalo)
    if (xhalo /= yhalo) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_operators_mod (ocean_operators_init): with static memory, xhalo must equal yhalo')
    endif 
    nk = Grd%nk
    allocate(fmx_tmask(isd:ied,jsd:jed,nk))
    allocate(fmy_tmask(isd:ied,jsd:jed,nk))
    halo = xhalo
#endif
    do k=1,nk
       fmx_tmask(:,:,k) = FMX(Grd%tmask(:,:,k))
       fmy_tmask(:,:,k) = FMY(Grd%tmask(:,:,k))
    enddo

! set up data for the operator used in barotropic domain.
   if( .NOT. Associated(Dom_bt) ) then
      call mpp_error(FATAL, 'ocean_operators_mod(ocean_operators_init): &
                &set_barotropic_domain must be called before calling ocean_operators_init')
   endif

   call get_halo_sizes(Dom_bt,xhalo, yhalo)
   if (xhalo /= yhalo) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_operators_mod (ocean_operators_init): xhalo must equal yhalo for barotropic domain')
   endif 
   barotropic_halo = xhalo
   ! when barotropic_halo > 1, use_legacy_DIV_UD must be false.
   if( barotropic_halo > 1 .AND. use_legacy_DIV_UD .AND. Grd%tripolar ) then
       call mpp_error(FATAL, &
      '==>Error in ocean_operators_mod (ocean_operators_init):  use_legacy_DIV_UD must be false '// &
        'when barotropic_halo > 1 and the grid is tripolar')
   endif 

   isd_bt = isc - barotropic_halo
   ied_bt = iec + barotropic_halo
   jsd_bt = jsc - barotropic_halo
   jed_bt = jec + barotropic_halo

   allocate(dxtn_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dyte_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dxte_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dytn_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dxter_bt(isd_bt:ied_bt, jsd_bt:jed_bt))
   allocate(dytnr_bt(isd_bt:ied_bt, jsd_bt:jed_bt))
   allocate(dxu_bt(isd_bt:ied_bt,   jsd_bt:jed_bt))
   allocate(dyu_bt(isd_bt:ied_bt,   jsd_bt:jed_bt))
   allocate(dat_bt(isd_bt:ied_bt,   jsd_bt:jed_bt))
   allocate(datr_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dxur_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(dyur_bt(isd_bt:ied_bt,  jsd_bt:jed_bt))
   allocate(umask_bt(isd_bt:ied_bt, jsd_bt:jed_bt))
   allocate(tmasken_bt(isd_bt:ied_bt, jsd_bt:jed_bt, 2))
   dxtn_bt    = 0
   dyte_bt    = 0
   dxte_bt    = 0
   dytn_bt    = 0
   dxter_bt   = 0
   dytnr_bt   = 0
   dxu_bt     = 0
   dyu_bt     = 0
   dat_bt     = 0
   datr_bt    = 0
   dxur_bt    = 0
   dyur_bt    = 0
   umask_bt   = 0
   tmasken_bt = 0

   dxtn_bt(isd:ied, jsd:jed)       = Grid%dxtn
   dyte_bt(isd:ied, jsd:jed)       = Grid%dyte
   dxte_bt(isd:ied, jsd:jed)       = Grid%dxte
   dytn_bt(isd:ied, jsd:jed)       = Grid%dytn
   dxter_bt(isd:ied, jsd:jed)      = Grid%dxter
   dytnr_bt(isd:ied, jsd:jed)      = Grid%dytnr
   dxu_bt(isd:ied, jsd:jed)        = Grid%dxu
   dyu_bt(isd:ied, jsd:jed)        = Grid%dyu
   dat_bt(isd:ied, jsd:jed)        = Grid%dat
   datr_bt(isd:ied, jsd:jed)       = Grid%datr
   dxur_bt(isd:ied, jsd:jed)       = Grid%dxur
   dyur_bt(isd:ied, jsd:jed)       = Grid%dyur
   umask_bt(isd:ied, jsd:jed)      = Grid%umask(:,:,1)
   tmasken_bt(isd:ied, jsd:jed, 1) = Grid%tmasken(:,:,1,1)
   tmasken_bt(isd:ied, jsd:jed, 2) = Grid%tmasken(:,:,1,2)
   if(barotropic_halo > 1 ) then
      call mpp_update_domains(dxtn_bt,  Dom_bt%domain2d, position=NORTH  )
      call mpp_update_domains(dyte_bt,  Dom_bt%domain2d, position=EAST   )
      call mpp_update_domains(dxte_bt,  Dom_bt%domain2d, position=EAST   )
      call mpp_update_domains(dytn_bt,  Dom_bt%domain2d, position=NORTH  )
      call mpp_update_domains(dxter_bt, Dom_bt%domain2d, position=EAST   )
      call mpp_update_domains(dytnr_bt, Dom_bt%domain2d, position=NORTH  )
      call mpp_update_domains(dxu_bt,   Dom_bt%domain2d, position=CORNER )
      call mpp_update_domains(dyu_bt,   Dom_bt%domain2d, position=CORNER )
      call mpp_update_domains(dxur_bt,  Dom_bt%domain2d, position=CORNER )
      call mpp_update_domains(dyur_bt,  Dom_bt%domain2d, position=CORNER )
      call mpp_update_domains(dat_bt,   Dom_bt%domain2d )
      call mpp_update_domains(datr_bt,  Dom_bt%domain2d )
      call mpp_update_domains(umask_bt, Dom_bt%domain2d, position=CORNER )
      call mpp_update_domains(tmasken_bt(:,:,1), Dom_bt%domain2d, position=EAST )
      call mpp_update_domains(tmasken_bt(:,:,2), Dom_bt%domain2d, position=NORTH)
   endif


    write( stdlogunit,'(/a/)') trim(version)

    return

  end subroutine ocean_operators_init
! </SUBROUTINE> NAME="ocean_operators_init"

!#######################################################################
! <SUBROUTINE NAME="set_barotropic_domain">
!
! <DESCRIPTION>
! Set the barotropic domain used in barotropic time step.
! </DESCRIPTION>
!
! <IN NAME="Domain_in" TYPE="type(ocean_domain_type)">
! Store the barotropic domain.
! </IN>
!
subroutine set_barotropic_domain( Domain_in )
   type(ocean_domain_type), intent(inout), target :: Domain_in

   Dom_bt => Domain_in

end subroutine set_barotropic_domain
! </SUBROUTINE> NAME="set_barotropic_domain"


!#######################################################################
! <FUNCTION NAME="get_use_legacy_DIV_UD">
!
! <DESCRIPTION>
! Return the value of ocean_operators_nml variable use_legacy_DIV_UD
! </DESCRIPTION>
!
function get_use_legacy_DIV_UD()
   logical :: get_use_legacy_DIV_UD

   get_use_legacy_DIV_UD = use_legacy_DIV_UD

   return

end function get_use_legacy_DIV_UD
! </FUNCTION>


!#######################################################################
! <FUNCTION NAME="REMAP_NT_TO_NU">
!
! <DESCRIPTION>
! REMAP_NT_TO_NU remaps a normal flux at the north 
! face of T-cells to the north face of U-cells
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be remapped 
! </IN>
!
function REMAP_NT_TO_NU(a) 
    real, intent(in), dimension(isd:,jsd:) :: a
    real, dimension(isd:ied,jsd:jed)       :: REMAP_NT_TO_NU
    integer :: i, j

    do j=jsc-halo,jec+halo-1
       do i=isc-halo,iec+halo-1
          REMAP_NT_TO_NU(i,j) = ((a(i,j)*Grd%duw(i,j) + a(i+1,j)*Grd%due(i,j))*Grd%dus(i,j+1) +&
               ((a(i,j+1)*Grd%duw(i,j+1) + a(i+1,j+1)*Grd%due(i,j+1)))*Grd%dun(i,j))*Grd%dater(i,j+1)
       enddo
    enddo
    REMAP_NT_TO_NU(iec+halo,:) = 0.0
    REMAP_NT_TO_NU(:,jec+halo) = 0.0

end function REMAP_NT_TO_NU
! </FUNCTION> NAME="REMAP_NT_TO_NU"


!#######################################################################
! <FUNCTION NAME="REMAP_ET_TO_EU">
!
! <DESCRIPTION>
! REMAP_ET_TO_EU remaps a normal flux at the east 
! face of T-cells to the east face of U-cells
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be remapped 
! </IN>
!
function REMAP_ET_TO_EU(a) 

    real, intent(in), dimension(isd:,jsd:) :: a
    real, dimension(isd:ied,jsd:jed)       :: REMAP_ET_TO_EU
    integer :: i, j

    do j=jsc-halo,jec+halo-1
       do i=isc-halo,iec+halo-1
          REMAP_ET_TO_EU(i,j) = ((a(i,j)*Grd%dus(i,j) + a(i,j+1)*Grd%dun(i,j))*Grd%duw(i+1,j) +&
               ((a(i+1,j)*Grd%dus(i+1,j) + a(i+1,j+1)*Grd%dun(i+1,j)))*Grd%due(i,j))*Grd%datnr(i+1,j)
       enddo
    enddo
    REMAP_ET_TO_EU(iec+halo,:) = 0.0
    REMAP_ET_TO_EU(:,jec+halo) = 0.0

end function REMAP_ET_TO_EU
! </FUNCTION> NAME="REMAP_ET_TO_EU"


!#######################################################################
! <FUNCTION NAME="REMAP_BT_TO_BU">
!
! <DESCRIPTION>
! REMAP_BT_TO_BU remaps a T-cell thickness or 
! vertical velocity on the base of T-cells to U-cells
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be remapped 
! </IN>
function REMAP_BT_TO_BU(a) 

    real, intent(in), dimension(isd:,jsd:) :: a
    real, dimension(isd:ied,jsd:jed)       :: REMAP_BT_TO_BU
    integer :: i, j

    do j=jsc-halo,jec+halo-1
       do i=isc-halo,iec+halo-1
          REMAP_BT_TO_BU(i,j) = (a(i,j)*Grd%dte(i,j)*Grd%dus(i,j) + a(i+1,j)*Grd%dtw(i+1,j)*Grd%dus(i,j)&
               + a(i,j+1)*Grd%dte(i,j+1)*Grd%dun(i,j) + a(i+1,j+1)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
       enddo
    enddo
    REMAP_BT_TO_BU(iec+halo,:) = 0.0
    REMAP_BT_TO_BU(:,jec+halo) = 0.0

end function REMAP_BT_TO_BU
! </FUNCTION> NAME="REMAP_BT_TO_BU"


!#######################################################################
! <FUNCTION NAME="DIV_UD">
!
! <DESCRIPTION>
! Compute divergence of vertically integrated velocity.  
! For MOM with generalized vertical coordinates, ud has an extra
! density factor built in.  
!
! The Bgrid uh and vh fields are located on the corner points, so require
! some backward averaging.  
!
! The Cgrid uh and vh fields are located on the T-cell faces, so need no
! spatial averaging. 
!
! Bgrid code is a speedier version of
! <BR/>
! uhy(:,:) = BAY(ud(:,:,1)*dyu(:,:))
! <BR/>
! vhx(:,:) = BAX(ud(:,:,2)*dxu(:,:))
! <BR/>
! DIV_UD(ud) = BDX_ET(uhy(:,:)/dyte(:,:)) + BDY_NT(vhx(:,:)/dxtn(:,:))
!
! Cgrid divergence is straightforward divergence operator.
!
! </DESCRIPTION>
!
function DIV_UD (ud, halo_in, halo_out )

    integer, intent(in)                                                  :: halo_in, halo_out
    real, intent(in), dimension(isc-halo_in:,jsc-halo_in:,:)             :: ud
    real, dimension(isc-halo_out:iec+halo_out,jsc-halo_out:jec+halo_out) :: DIV_UD
    integer :: i, j, jstart, jend, istart, iend
    real    :: uh, uhim, uhjm, uhimjm, vh, vhim, vhjm, vhimjm
    real    :: uh_bay, vh_bax, uhim_bay, vhjm_bax

    if( halo_in < halo_out) call mpp_error(FATAL, &
         'ocean_operators(DIV_UD): halo_in must be no less than halo_out')

    istart = isc - halo_out
    iend   = iec + halo_out
    jstart = jsc - halo_out
    jend   = jec + halo_out

    DIV_UD = 0.d0

    if(horz_grid == MOM_CGRID) then

          do j=jstart+1,jend
             do i=istart+1,iend
                uhim = ud(i-1,j,1)*dat_bt(i-1,j)*dxter_bt(i-1,j)
                vhjm = ud(i,j-1,2)*dat_bt(i,j-1)*dytnr_bt(i,j-1)
                uh   = ud(i,j,1)*dat_bt(i,j)*dxter_bt(i,j)
                vh   = ud(i,j,2)*dat_bt(i,j)*dytnr_bt(i,j)
                DIV_UD(i,j) = ((uh - uhim) + (vh - vhjm))*datr_bt(i,j)
             enddo
          enddo

    else ! Bgrid 

       if( use_legacy_DIV_UD) then
          do j=jstart+1,jend
             i = istart+1
             uhim = ud(i-1,j,1)*dyu_bt(i-1,j)
             vhim = ud(i-1,j,2)*dxu_bt(i-1,j)
             uhimjm = ud(i-1,j-1,1)*dyu_bt(i-1,j-1)
             vhimjm = ud(i-1,j-1,2)*dxu_bt(i-1,j-1)
             uhim_bay = 0.5*(uhim+uhimjm)
             do i=istart+1,iend
                uh  = ud(i,j,1)*dyu_bt(i,j)
                vh  = ud(i,j,2)*dxu_bt(i,j)

                uhjm = ud(i,j-1,1)*dyu_bt(i,j-1)
                vhjm = ud(i,j-1,2)*dxu_bt(i,j-1)

                uh_bay = 0.5*(uh+uhjm)
                vh_bax = 0.5*(vh+vhim)

                vhjm_bax = 0.5*(vhjm+vhimjm)

                DIV_UD(i,j) = (uh_bay - uhim_bay + vh_bax - vhjm_bax)*datr_bt(i,j)

                uhim = uh
                vhim = vh
                uhimjm = uhjm
                vhimjm = vhjm
                uhim_bay = uh_bay
             enddo
          enddo

       else

          do j=jstart+1,jend
             i = istart+1
             uhim = ud(i-1,j,1)*dyu_bt(i-1,j)
             vhim = ud(i-1,j,2)*dxu_bt(i-1,j) 
             uhimjm = ud(i-1,j-1,1)*dyu_bt(i-1,j-1)
             vhimjm = ud(i-1,j-1,2)*dxu_bt(i-1,j-1)
             uhim_bay = 0.5*(uhim+uhimjm)
             do i=istart+1,iend
                uh  = ud(i,j,1)*dyu_bt(i,j)
                vh  = ud(i,j,2)*dxu_bt(i,j)

                uhjm = ud(i,j-1,1)*dyu_bt(i,j-1)
                vhjm = ud(i,j-1,2)*dxu_bt(i,j-1)

                uh_bay = 0.5*(uh+uhjm)
                vh_bax = 0.5*(vh+vhim)

                vhjm_bax = 0.5*(vhjm+vhimjm)

                DIV_UD(i,j) = ((uh_bay - uhim_bay) + (vh_bax - vhjm_bax))*datr_bt(i,j)

                uhim = uh
                vhim = vh
                uhimjm = uhjm
                vhimjm = vhjm
                uhim_bay = uh_bay
             enddo
          enddo

       endif  ! endif for legacy operator 

    endif ! endif for horz_grid 

    DIV_UD(:,:jstart) = 0.0
    DIV_UD(:istart,:) = 0.0

end function DIV_UD
! </FUNCTION> NAME="DIV_UD"


!#######################################################################
! <FUNCTION NAME="GRAD_BAROTROPIC_P">
!
! <DESCRIPTION>
! Compute horizontal gradient of the pressure field associated with 
! either the free surface height or the bottom pressure. 
!
! Account taken here for either Bgrid or Cgrid operators. 
! For the Bgrid, the gradient is centered onto the U-cell 
! for use in updating barotropic velocity. For the Cgrid, the 
! two components are centred on the T-cell faces.  
!
! The Bgrid algorithm is a speedier version of
! <BR/>
! grad_barotropic_p(:,:,1) = FDX_NT(FAY(press(:,:)))
! <BR/>
! grad_barotropic_p(:,:,2) = FDY_ET(FAX(press(:,:)))
! </DESCRIPTION>
!
function GRAD_BAROTROPIC_P (press, halo_in, halo_out )

    integer, intent(in)                                    :: halo_in, halo_out
    real, intent(in), dimension(isc-halo_in:,jsc-halo_in:) :: press
    real, dimension(isc-halo_out:iec+halo_out,jsc-halo_out:jec+halo_out,2) :: GRAD_BAROTROPIC_P
    integer :: i, j, jstart, jend, istart, iend
    real    :: p_fay, p_fax, pip_fay, pjp_fax

    if( halo_in < halo_out) then
         call mpp_error(FATAL,'ocean_operators(GRAD_BAROTROPIC_P): halo_in must be no less than halo_out')
    endif 

    istart = isc - halo_out
    iend   = iec + halo_out
    jstart = jsc - halo_out
    jend   = jec + halo_out

    if(horz_grid == MOM_BGRID) then 

       do j=jstart,jend-1
          do i=istart,iend-1
             p_fay = 0.5*(press(i,j) + press(i,j+1))
             p_fax = 0.5*(press(i,j) + press(i+1,j))
             pip_fay = 0.5*(press(i+1,j) + press(i+1,j+1))
             pjp_fax = 0.5*(press(i,j+1) + press(i+1,j+1))

             GRAD_BAROTROPIC_P(i,j,1) = (pip_fay - p_fay)*dxur_bt(i,j)*umask_bt(i,j)
             GRAD_BAROTROPIC_P(i,j,2) = (pjp_fax - p_fax)*dyur_bt(i,j)*umask_bt(i,j)
          enddo
       enddo

    else ! cgrid 

       do j=jstart,jend-1
          do i=istart,iend-1
             GRAD_BAROTROPIC_P(i,j,1) = (press(i+1,j) - press(i,j))*dxter_bt(i,j)*tmasken_bt(i,j,1)
             GRAD_BAROTROPIC_P(i,j,2) = (press(i,j+1) - press(i,j))*dytnr_bt(i,j)*tmasken_bt(i,j,2)
          enddo
       enddo

    endif 

    GRAD_BAROTROPIC_P(iend:,:,:) = 0.0
    GRAD_BAROTROPIC_P(:,jend:,:) = 0.0

end function GRAD_BAROTROPIC_P
! </FUNCTION> NAME="GRAD_BAROTROPIC_P"

!#######################################################################
! <FUNCTION NAME="S2D">
!
! <DESCRIPTION>
! Smooth a 2D field with a 2D version of a 1D filter with weights (1/4, 1/2, 1/4)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be smoothed 
! </IN>
!
function S2D(a) 

    real, intent(in), dimension(isd:,jsd:) :: a
    real, dimension(isd:ied,jsd:jed)       :: S2D
    real, dimension(-1:1,-1:1)             :: fir
    integer :: i, j, iq, jq

    fir(-1,-1) = onesixteenth
    fir(0,-1)  = oneeigth
    fir(1,-1)  = onesixteenth
    fir(-1,0)  = oneeigth
    fir(0,0)   = onefourth
    fir(1,0)   = oneeigth
    fir(-1,1)  = onesixteenth
    fir(0,1)   = oneeigth
    fir(1,1)   = onesixteenth         
    do j=jsc-halo+1,jec+halo-1
       do i=isc-halo+1,iec+halo-1
          S2D(i,j) = 0.0
          do jq=-1,1
             do iq=-1,1
                S2D(i,j) = S2D(i,j) + fir(iq,jq)*a(i+iq,j+jq)
             enddo
          enddo
       enddo
    enddo
    S2D(iec+halo,:) = 0.0
    S2D(isc-halo,:) = 0.0
    S2D(:,jec+halo) = 0.0
    S2D(:,jsc-halo) = 0.0

end function S2D
! </FUNCTION> NAME="S2D"


!#######################################################################
! <FUNCTION NAME="LAP_T">
!
! <DESCRIPTION>
!  Compute horizontal 5-point Laplacian operator on eta_t.
!  Result lives at T-cell center. 
!
!  Redundancy update for tripolar is needed to conserve  
!  total volume and tracer.  It is likely unimportant 
!  when call LAP_T from within the barotropic loop. Yet it 
!  is essential when call LAP_T from ocean_surface_smooth.
!
!  Mixing coefficient is assumed to be centred on the T-cell.
!  It is averaged to compute its value on the i-face and j-face
!  for computing fluxes. 
!
! </DESCRIPTION>
!
function LAP_T (a, mix, update)

    real, intent(in),    dimension(isd:ied,jsd:jed) :: a
    real, intent(in),    dimension(isd:ied,jsd:jed) :: mix
    logical, intent(in),                   optional :: update

    real, dimension(isd:ied,jsd:jed) :: fx, fy, LAP_T
    real                             :: chg

    integer :: i, j, k
    logical :: update_domains

    k = 1

    if (PRESENT(update)) then 
       update_domains = update
    else 
       update_domains = .false.
    endif  

    do j = jsc-halo, jec+halo
       do i = isc-halo,iec+halo-1
          fx(i,j) = 0.5*(mix(i,j)+mix(i+1,j))*((a(i+1,j)-a(i,j))*Grd%dxter(i,j))*fmx_tmask(i,j,k)
       enddo
    enddo
    fx(iec+halo,:) = 0.0

    do j = jsc-halo, jec+halo-1
       do i = isc-halo,iec+halo
          fy(i,j) = 0.5*(mix(i,j)+mix(i,j+1))*((a(i,j+1)-a(i,j))*Grd%dytnr(i,j))*fmy_tmask(i,j,k)
       enddo
    enddo
    fy(:,jec+halo) = 0.0

    if(Grd%tripolar .and. update_domains) then
       call mpp_update_domains(fx, fy, Dom_flux%domain2d, gridtype=CGRID_NE)
    endif 

    LAP_T(:,:) = 0.0

    do j=jsc,jec
       do i=isc,iec
          chg        = (Grd%dyte(i,j)*fx(i,j) - Grd%dyte(i-1,j)*fx(i-1,j))*Grd%datr(i,j) &
                     + (Grd%dxtn(i,j)*fy(i,j) - Grd%dxtn(i,j-1)*fy(i,j-1))*Grd%datr(i,j)
          LAP_T(i,j) = Grd%tmask(i,j,k)*chg
       enddo
    enddo
 
end function LAP_T
! </FUNCTION> NAME="LAP_T">


!#######################################################################
! <FUNCTION NAME="FAX">
!
! <DESCRIPTION>
! Forwards average in the i-direction on the X-axis.
! If input is a(i,j) then output is defined at (i+1/2,j)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be averaged  
! </IN>
function FAX(a) 

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FAX
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1
        FAX(i,j) = 0.5*(a(i,j) + a(i+1,j))
     enddo
     FAX(iec+halo,j) = 0.0
  enddo
end function FAX
! </FUNCTION> NAME="FAX"


!#######################################################################
! <FUNCTION NAME="BAX">
!
! <DESCRIPTION>
! Backwards average in the i-direction along the X-axis.
! If input is a(i,j) then output is defined at (i-1/2,j)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be averaged  
! </IN>
!
function BAX(a) 

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: BAX
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo+1,iec+halo
        BAX(i,j) = 0.5*(a(i,j) + a(i-1,j))
     enddo
     BAX(isc-halo,j) = 0.0
  enddo
end function BAX
! </FUNCTION> NAME="BAX"


!#######################################################################
! <FUNCTION NAME="FAY">
!
! <DESCRIPTION>
! Forwards average in the j-direction on the Y-axis
! If input is a(i,j) then output is defined at (i,j+1/2)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be averaged  
! </IN>
function FAY(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FAY
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FAY(i,j) = 0.5*(a(i,j) + a(i,j+1))
     enddo
  enddo
  FAY(:,jec+halo) = 0.0

end function FAY
! </FUNCTION> NAME="FAY"


!#######################################################################
! <FUNCTION NAME="BAY">
!
! <DESCRIPTION>
! Backwards average in the j-direction along the Y-axis
! If input is a(i,j) then output is defined at (i,j-1/2)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be averaged  
! </IN>
function BAY(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed) :: BAY
  integer :: i, j

  do j=jsc-halo+1,jec+halo
     do i=isc-halo,iec+halo
        BAY(i,j) = 0.5*(a(i,j) + a(i,j-1))
     enddo
  enddo
  BAY(:,jsc-halo) = 0.0

end function BAY
! </FUNCTION> NAME="BAY"


!#######################################################################
! <FUNCTION NAME="BDX_EU">
!
! <DESCRIPTION>
! Backwards Derivative in X of a quantity defined on the East face of a U-cell
! If input is a(i,j) then output is defined at (i-1/2,j)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function BDX_EU(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDX_EU
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo+1,iec+halo
        BDX_EU(i,j) = (Grd%dyue(i,j)*a(i,j) - Grd%dyue(i-1,j)*a(i-1,j))*Grd%daur(i,j)
     enddo
     BDX_EU(isc-halo,j) = 0.0
  enddo

end function BDX_EU
! </FUNCTION> NAME="BDX_EU"


!#######################################################################
! <FUNCTION NAME="BDX_ET">
!
! <DESCRIPTION>
! Backwards derivative in X of a quantity defined on the East face of a T-cell.
! If input is a(i,j) then output is defined at (i-1/2,j)
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function BDX_ET(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDX_ET
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo+1,iec+halo
        BDX_ET(i,j) = (Grd%dyte(i,j)*a(i,j) - Grd%dyte(i-1,j)*a(i-1,j))*Grd%datr(i,j) 
     enddo
     BDX_ET(isc-halo,j) = 0.0
  enddo

end function BDX_ET
! </FUNCTION> NAME="BDX_ET"


!#######################################################################
! <FUNCTION NAME="FDX_U">
!
! <DESCRIPTION>
! Forward Derivative in X of a quantity defined on the grid point of a U-cell.
! If input is a(i,j) then output is defined at (i+1/2,j).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDX_U(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDX_U
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1
        FDX_U(i,j) = (a(i+1,j) - a(i,j))*Grd%dxuer(i,j)
     enddo
     FDX_U(iec+halo,j) = 0.0
  enddo

end function FDX_U
! </FUNCTION> NAME="FDX_U"


!#######################################################################
! <FUNCTION NAME="FDX_ZT">
!
! <DESCRIPTION>
!
! Forward Derivative in X of a quantity defined on a tracer grid point
! where it is necessary to take derivative with depth held constant.  
! When grid points live at different depths, then have an extra 
! contribution to the derivative.  
!
! Input a(i,j,1) is at the grid point of a T-cell at level k-1.
! Input a(i,j,2) is at the grid point of a T-cell at level k.
!
! Output is defined at (i+1/2,j) which is at the east face in a T-cell at level k.
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
! <IN NAME="k" TYPE="integer">
! Depth level 
! </IN>
function FDX_ZT(a,k)

  real, intent(in), dimension(isd:,jsd:,:) :: a
  real, dimension(isd:ied,jsd:jed) :: FDX_ZT
  integer, intent(in)  :: k
  integer :: i, j

  do i=isc-halo,iec+halo-1
     FDX_ZT(i,:) = (a(i+1,:,2) - a(i,:,2))*Grd%dxter(i,:)
  enddo
  FDX_ZT(iec+halo,:) = 0.0

!   expanded form of 
!   FDX_ZT(:,:) = FDX_ZT(:,:) - FAX(FDZ_T(a(:,:,1:2),k-1))*FDX_T(Thk%depth_zt(:,:,k))

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1   
        FDX_ZT(i,j) = FDX_ZT(i,j) - 0.5*( -(a(i,  j,2) - a(i,  j,1))/Thk%dzwt(i,  j,k-1)   &
                                          -(a(i+1,j,2) - a(i+1,j,1))/Thk%dzwt(i+1,j,k-1) ) &
                                   *((Thk%depth_zt(i+1,j,k) - Thk%depth_zt(i,j,k))*Grd%dxter(i,j))
     enddo
  enddo

end function FDX_ZT
! </FUNCTION> NAME="FDX_ZT"


!#######################################################################
! <FUNCTION NAME="FDX_T">
!
! <DESCRIPTION>
!
! Forward Derivative in X of a quantity defined on the grid point
! of a T-cell.
!
! For lateral derivatives where vertical coordinate is held constant.  
!
! Input a(i,j) is at the grid point of a T-cell
! Output is defined at (i+1/2,j) which is at the east face in a T-cell.
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDX_T(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDX_T
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1
        FDX_T(i,j) = (a(i+1,j) - a(i,j))*Grd%dxter(i,j) 
     enddo
     FDX_T(iec+halo,j) = 0.0
  enddo

end function FDX_T
! </FUNCTION> NAME="FDX_T"


!#######################################################################
! <FUNCTION NAME="FDY_ZT">
!
! <DESCRIPTION>
! Forward Derivative in Y of a quantity defined on a tracer grid point
! where it is necessary to take derivative with depth held constant.  
! When grid points live at different depths, then have an extra 
! contribution to the derivative.  
!
! Input a(i,j,1) is at the grid point of a T-cell at level k-1.
! Input a(i,j,2) is at the grid point of a T-cell at level k.
!
! Output is defined at (i,j+1/2) which is at the north face in a T-cell at level k.
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
! <IN NAME="k" TYPE="integer">
! Depth level 
! </IN>
function FDY_ZT(a,k)

  real, intent(in), dimension(isd:,jsd:,:) :: a
  real, dimension(isd:ied,jsd:jed) ::  FDY_ZT
  integer, intent(in)  :: k
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FDY_ZT(i,j) = (a(i,j+1,2) - a(i,j,2))*Grd%dytnr(i,j)
     enddo
  enddo
  FDY_ZT(:,jec+halo) = 0.0

! expanded form of 
! FDY_ZT(:,:) = FDY_ZT(:,:) - FAY(FDZ_T(a(:,:,1:2),k-1))*FDY_T(Grd%depth_zt(:,:,k))

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FDY_ZT(i,j) = FDY_ZT(i,j) - 0.5*( -(a(i,  j,2) - a(i,  j,1))/Thk%dzwt(i,  j,k-1)   &
                                          -(a(i,j+1,2) - a(i,j+1,1))/Thk%dzwt(i,j+1,k-1) ) &
                                    *((Thk%depth_zt(i,j+1,k) - Thk%depth_zt(i,j,k))*Grd%dytnr(i,j))
     enddo
  enddo


end function FDY_ZT
! </FUNCTION> NAME="FDY_ZT"


!#######################################################################
! <FUNCTION NAME="FDY_T">
!
! <DESCRIPTION>
! Forward Derivative in Y of a quantity defined on the grid Point
! of a T-cell.
!
! For lateral derivatives where vertical coordinate is held constant.  
!
! Input a(i,j) is at the grid point of a T-cell.
! Output is defined at (i,j+1/2) which is at the north face in a T-cell.
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDY_T(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDY_T
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FDY_T(i,j) = (a(i,j+1) - a(i,j))*Grd%dytnr(i,j)
     enddo
  enddo
  do i=isc-halo,iec+halo
    FDY_T(i,jec+halo) = 0.0
  enddo

end function FDY_T
! </FUNCTION> NAME="FDY_T"


!#######################################################################
! <FUNCTION NAME="FDX_NT">
!
! <DESCRIPTION>
! Forward Derivative in X of a quantity defined on the North face of a T-cell.
! If input is a(i,j) then output is defined at (i+1/2,j).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDX_NT(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDX_NT
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1
        FDX_NT(i,j) = (a(i+1,j) - a(i,j))*Grd%dxur(i,j)
     enddo
     FDX_NT(iec+halo,j) = 0.0
  enddo

end function FDX_NT
! </FUNCTION> NAME="FDX_NT"


!#######################################################################
! <FUNCTION NAME="BDY_NU">
!
! <DESCRIPTION>
! Backward Derivative in Y of a quantity defined on the North face of a U-cell.
! If input is a(i,j) then output is defined at (i,j-1/2).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function BDY_NU(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDY_NU
  integer :: i, j

  do j=jsc-halo+1,jec+halo
     do i=isc-halo,iec+halo
        BDY_NU(i,j) = (Grd%dxun(i,j)*a(i,j) - Grd%dxun(i,j-1)*a(i,j-1))*Grd%daur(i,j)
     enddo
  enddo
  BDY_NU(:,jsc-halo) = 0.0

end function BDY_NU
! </FUNCTION> NAME="BDY_NU"


!#######################################################################
! <FUNCTION NAME="BDY_NT">
!
! <DESCRIPTION>
! Backward Derivative in Y of a quantity defined on the North face of a T-cell.
! If input is a(i,j) then output is defined at (i,j-1/2).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function BDY_NT(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDY_NT
  integer :: i, j

  do j=jsc-halo+1,jec+halo
     do i=isc-halo,iec+halo
        BDY_NT(i,j) = (Grd%dxtn(i,j)*a(i,j) - Grd%dxtn(i,j-1)*a(i,j-1))*Grd%datr(i,j)  
     enddo
  enddo
  BDY_NT(:,jsc-halo) = 0.0

end function BDY_NT
! </FUNCTION> NAME="BDY_NT"


!#######################################################################
! <FUNCTION NAME="FDY_U">
!
! <DESCRIPTION>
! Forward Derivative in Y of a quantity defined on the grid Point of a U-cell.
! If input is a(i,j) then output is defined at (i,j+1/2).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDY_U(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDY_U
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FDY_U(i,j) = (a(i,j+1) - a(i,j))*Grd%dyunr(i,j)
     enddo
  enddo
  FDY_U(:,jec+halo) = 0.0

end function FDY_U
! </FUNCTION> NAME="FDY_U"


!#######################################################################
! <FUNCTION NAME="FDY_ET">
!
! <DESCRIPTION>
! Forward Derivative in Y of a quantity defined on the East face of a T-cell.
! If input is a(i,j) then output is defined at (i,j+1/2).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
function FDY_ET(a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FDY_ET
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FDY_ET(i,j) = (a(i,j+1) - a(i,j))*Grd%dyur(i,j)
     enddo
  enddo
  FDY_ET(:,jec+halo) = 0.0

end function FDY_ET
! </FUNCTION> NAME="FDY_ET"


!#######################################################################
! <FUNCTION NAME="FDZ_T">
!
! <DESCRIPTION>
! Forward Derivative in Z of a field on the T-cell point at level k.
!
! input a(i,j,1) is at the grid point of a T-cell at level k.
! input a(i,j,2) is at the grid point of a T-cell at level k+1.
!
! output is at (i,j,k+3/2) which is bottom face of T-cell at level k.
!
! minus sign due to convention that z-increases upwards, whereas
! k increases downward.  
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to be finite differenced 
! </IN>
! <IN NAME="k" TYPE="integer">
! Depth level index 
! </IN>
!
function FDZ_T(a, k)

  real, intent(in), dimension(isd:,jsd:,:) :: a
  real, dimension(isd:ied,jsd:jed)         :: FDZ_T
  integer, intent(in) :: k

  FDZ_T(:,:) =  -(a(:,:,2) - a(:,:,1))/Thk%dzwt(:,:,k) 

end function FDZ_T
! </FUNCTION> NAME="FDZ_T"


!#######################################################################
! <FUNCTION NAME="FMX">
!
! <DESCRIPTION>
! Forwards Minimum in the X direction.
! If input is a(i,j) then output is defined at (i+1/2,j).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to find minimum
! </IN>
!
function FMX (a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FMX
  integer :: i, j

  do j=jsc-halo,jec+halo
     do i=isc-halo,iec+halo-1
        FMX(i,j) = min(a(i+1,j), a(i,j))
     enddo
     FMX(iec+halo,j) = 0.0
  enddo

end function FMX
! </FUNCTION> NAME="FMX"


!#######################################################################
! <FUNCTION NAME="FMY">
!
! <DESCRIPTION>
! Forwards Minimum in the Y direction.
! If input is a(i,j) then output is defined at (i,j+1/2).
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! Field to find minimum
! </IN>
function FMY (a)

  real, intent(in), dimension(isd:,jsd:) :: a
  real, dimension(isd:ied,jsd:jed)       :: FMY
  integer :: i, j

  do j=jsc-halo,jec+halo-1
     do i=isc-halo,iec+halo
        FMY(i,j) = min(a(i,j+1), a(i,j))
     enddo
  enddo
  FMY(:,jec+halo) = 0.0

end function FMY
! </FUNCTION> NAME="FMY"



end module ocean_operators_mod

