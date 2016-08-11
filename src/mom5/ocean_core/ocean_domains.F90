module ocean_domains_mod
!  
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matthew Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Set the ocean domain parameters. 
!</OVERVIEW>
!
!<DESCRIPTION>
! The module computes the horizontal domain parameters needed to 
! run MOM in a parallel computational environment. 
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison,  R.C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_domains_nml">
!  <DATA NAME="halo" UNITS="dimensionless" TYPE="integer">
!  For specifying the halo size by hand.
!  </DATA> 
!  <DATA NAME="max_tracers" UNITS="dimensionless" TYPE="integer">
!  temporary - need to call domains_init before tracer_init 
!  Used for computing mpp_stack_size.
!  </DATA> 
!  <DATA NAME="x_cyclic_offset" UNITS="dimensionless" TYPE="integer">
!  offset to be applied on x-direction boundary condition. Its value could
!  be positive or negative and the default value is 0. When the y-direction 
!  boundary condition is folded-north(tripolar grid), x_cyclic_offset must
!  be 0. For torus (cyclic in x and y-direction), at least one of 
!  x_cyclic_offset and y_cyclic_offset must be 0.   
!  </DATA> 
!  <DATA NAME="y_cyclic_offset" UNITS="dimensionless" TYPE="integer">
!  offset to be applied on y-direction boundary condition. Its value could
!  be positive or negative and the default value is 0. For torus (cyclic 
!  in x and y-direction), at least one of x_cyclic_offset and 
!  y_cyclic_offset must be 0.   
!  </DATA> 
!</NAMELIST>

use fms_mod,         only: open_namelist_file, write_version_number, check_nml_error, close_file
use mpp_domains_mod, only: mpp_define_layout, mpp_define_domains, mpp_domains_set_stack_size
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_mod, only: FOLD_NORTH_EDGE, CYCLIC_GLOBAL_DOMAIN, mpp_define_io_domain
use mpp_mod,         only: input_nml_file, mpp_max, mpp_min, mpp_npes, mpp_error, FATAL, stdout, stdlog

use ocean_types_mod, only: ocean_domain_type, ocean_grid_type
use ocean_types_mod, only: ocean_prog_tracer_type, ocean_velocity_type

implicit none

private

integer, dimension(2) :: domain_layout=(/1,1/)
integer, dimension(2) :: io_domain_layout=(/0,0/)

character(len=128) :: version='$Id: ocean_domains.F90,v 20.0 2013/12/14 00:10:42 fms Exp $'
character(len=128) :: tagname='$Name: tikal $'

character(len=32) :: name_default = 'mom_domain'

public ocean_domain_init
public set_ocean_domain
public get_local_indices
public get_halo_sizes
public get_global_indices
public get_active_indices
public reduce_active_domain
public get_domain_offsets

integer     :: halo=1            ! size (cells) of halo region
integer     :: max_tracers = 10  ! temporary - need to call domains_init before tracer_init 
integer     :: x_cyclic_offset=0 ! offset applied to x-direction cyclic boundary condition
integer     :: y_cyclic_offset=0 ! offset applied to y-direction cyclic boundary condition
namelist /ocean_domains_nml/ halo, max_tracers, x_cyclic_offset, y_cyclic_offset

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_domain_init">
!
! <DESCRIPTION>
! Initialise the domain module.
! </DESCRIPTION>
!
subroutine ocean_domain_init()

  integer :: ierr,ioun,io_status,stdoutunit,stdlogunit
  stdoutunit=stdout();stdlogunit=stdlog() 

  call write_version_number(version, tagname)

  ! provide for namelist over-ride
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=ocean_domains_nml, iostat=io_status)
ierr = check_nml_error(io_status,'ocean_domains_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_domains_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_domains_nml')
  call close_file(ioun)
#endif
  write (stdlogunit,ocean_domains_nml)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_domains_nml)

end subroutine ocean_domain_init
! </SUBROUTINE> NAME="set_ocean_domain"

!#######################################################################
! <SUBROUTINE NAME="set_ocean_domain">
!
! <DESCRIPTION>
! For setting the ocean domain layout and associated parameters. 
! </DESCRIPTION>
!
subroutine set_ocean_domain (Domain, Grid, xhalo, yhalo, name, layout, io_layout, maskmap)

  type(ocean_grid_type),   intent(in)           :: Grid
  type(ocean_domain_type), intent(inout)        :: Domain
  integer,                 intent(in), optional :: xhalo, yhalo, layout(2), io_layout(2)
  logical,                 intent(in), optional :: maskmap(:,:)
  character(len=*),        intent(in), optional :: name

  real              :: ph, pc
  integer           :: n, nhp, ncp, ncp_max, ncp_min, ncpx, ncpy, lay_out(2), xsiz, ysiz
  integer           :: mpp_stack_size=-1 
  character(len=32) :: name_
  character(len=4)  :: char_lay1, char_lay2, char_npes, char_xsiz, char_ysiz


  if (PRESENT(layout)) domain_layout = layout
  if (PRESENT(io_layout)) io_domain_layout = io_layout
  if (domain_layout(1) == 1 .and. domain_layout(2) == 1) then
     call mpp_define_layout ((/1,Grid%ni,1,Grid%nj/), mpp_npes(), lay_out)
  else
     lay_out = domain_layout
  endif

  Domain%layout = lay_out
  Domain%io_layout = io_domain_layout

  allocate(Domain%maskmap(lay_out(1), lay_out(2)) )

  if(PRESENT(maskmap)) then
     Domain%maskmap = maskmap
  else
     Domain%maskmap = .TRUE.
  end if

  if (PRESENT(name)) then
     name_ = trim(name)
  else
     name_ = trim(name_default)
  endif

  if (PRESENT(xhalo)) then
     Domain%xhalo=xhalo
  else
     Domain%xhalo=halo
  endif

  if (PRESENT(yhalo)) then
     Domain%yhalo=yhalo
  else
     Domain%yhalo=halo
  endif


! setup 2D domains for local x-y regions (including halos) on each processor.

  if (Grid%cyclic_x) then

      ! spherical with tripolar Arctic
      if (Grid%tripolar) then
          if(x_cyclic_offset .NE. 0) call  mpp_error(FATAL, &
             '==>Error in ocean_domains_mod: x_cyclic_offset must be 0 for tripolar grid')
          call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), lay_out, Domain%domain2d, maskmap = Domain%maskmap&
               , xflags = cyclic_global_domain, yflags = FOLD_NORTH_EDGE, xhalo=Domain%xhalo, yhalo=Domain%yhalo,name=name_)
          Domain%xflags = cyclic_global_domain
          Domain%yflags = FOLD_NORTH_EDGE

      ! torus 
      elseif(Grid%cyclic_y) then 
          if(x_cyclic_offset*y_cyclic_offset .NE. 0) call  mpp_error(FATAL, &
             '==>Error in ocean_domains_mod: at least one of x_cyclic_offset and y_cyclic_offset must be 0 for torus')
          call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), lay_out, Domain%domain2d, maskmap = Domain%maskmap&
               , xflags = cyclic_global_domain, yflags = cyclic_global_domain, xhalo=Domain%xhalo, yhalo=Domain%yhalo, name=name_&
               , x_cyclic_offset=x_cyclic_offset, y_cyclic_offset=y_cyclic_offset)
          Domain%xflags = cyclic_global_domain
          Domain%yflags = cyclic_global_domain
 
      ! cyclic in x and solid wall in y
      else 
          call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), lay_out, Domain%domain2d, maskmap = Domain%maskmap&
               , xflags = cyclic_global_domain, xhalo=Domain%xhalo, yhalo=Domain%yhalo,name=name_, x_cyclic_offset=x_cyclic_offset)
          Domain%xflags = cyclic_global_domain
          Domain%yflags = 0
      endif

  else

      ! cyclic in y and solid wall in x (only physically relevant on beta or f plane)
      if(Grid%cyclic_y) then 
          call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), lay_out, Domain%domain2d, maskmap = Domain%maskmap&
               , yflags = cyclic_global_domain, xhalo=Domain%xhalo, yhalo=Domain%yhalo,name=name_, y_cyclic_offset=y_cyclic_offset)
          Domain%xflags = 0
          Domain%yflags = cyclic_global_domain
      
      ! four solid walls
      else 
          call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), lay_out, Domain%domain2d, maskmap = Domain%maskmap&
               , xhalo=Domain%xhalo, yhalo=Domain%yhalo,name=name_)
          Domain%xflags = 0
          Domain%yflags = 0

      endif

  endif
  call mpp_define_io_domain(Domain%domain2d, io_domain_layout)

  Domain%x_cyclic_offset = x_cyclic_offset
  Domain%y_cyclic_offset = y_cyclic_offset

  ! define data and compute indices for local domain on this processor

  call mpp_get_compute_domain (Domain%domain2d, Domain%isc, Domain%iec, &
       Domain%jsc, Domain%jec)

  call mpp_get_data_domain (Domain%domain2d, Domain%isd, Domain%ied, &
       Domain%jsd, Domain%jed)

  call mpp_get_global_domain (Domain%domain2d, Domain%isg, Domain%ieg, &
       Domain%jsg, Domain%jeg)

  Domain%isa = Domain%isc
  Domain%iea = Domain%iec
  Domain%jsa = Domain%jsc
  Domain%jea = Domain%jec

  Domain%ioff = 0
  Domain%joff = 0
  
#ifdef MOM_STATIC_ARRAYS

  if (Domain%xhalo == 1 .AND. Domain%yhalo == 1) then
     xsiz = Domain%iec - Domain%isc + 1
     ysiz = Domain%jec - Domain%jsc + 1
     Domain%ioff = Domain%isc - 1     
     Domain%isc = 1
     Domain%isd = 0
     Domain%iec = xsiz
     Domain%ied = xsiz+1
     Domain%joff = Domain%jsc - 1          
     Domain%jsc = 1
     Domain%jsd = 0
     Domain%jec = ysiz
     Domain%jed = ysiz+1
     
     if (xsiz /= NI_LOCAL_ .OR. ysiz /= NJ_LOCAL_) then
        write( char_xsiz,'(i4)' ) xsiz
        write( char_ysiz,'(i4)' ) ysiz
        call mpp_error(FATAL,'==>Error in ocean_domains_mod: domain size (xsiz,ysiz) = ('&
        //trim(char_xsiz)//','//trim(char_ysiz)// ') does not match size set by preprocessor macro.&
          & Disable MOM_STATIC_ARRAYS cpp-preprocessor option to use dynamic allocation. ')
     endif
 endif

#endif


  ncp_max = (Domain%iec-Domain%isc+1)*(Domain%jec-Domain%jsc+1)
  ncp_min = ncp_max
  call mpp_max(ncp_max)
  call mpp_min(ncp_min)

  ncp = (Domain%iec-Domain%isc+1)*(Domain%jec-Domain%jsc+1)
  nhp = 4*halo**2 + 2*(Domain%iec-Domain%isc+1 + Domain%jec-Domain%jsc+1)*halo
  ph = 100.0*nhp/(nhp+ncp)
  pc = 100.0 - ph

  ! estimate stack size
  if (mpp_stack_size <= 0) then 
    ncpx = Domain%iec-Domain%isc+1
    ncpy = Domain%jec-Domain%jsc+1
    call mpp_max (ncpx)
    call mpp_max (ncpy)
    mpp_stack_size = 2*(ncpx + 2*max(Domain%xhalo,Domain%yhalo))*&
                    (ncpy + 2*max(Domain%xhalo,Domain%yhalo))*Grid%nk*(max_tracers) 
  endif

  call mpp_domains_set_stack_size(mpp_stack_size)

end subroutine set_ocean_domain
! </SUBROUTINE> NAME="set_ocean_domain"


!#######################################################################
! <SUBROUTINE NAME="get_local_indices">
!
! <DESCRIPTION>
! For getting local indices from domain derived type. 
! </DESCRIPTION>
!
subroutine get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)

  type(ocean_domain_type), intent(in) :: Domain
  integer, intent(out)                :: isd, ied, jsd, jed, isc, iec, jsc, jec
  
  isd=Domain%isd;ied=Domain%ied;jsd=Domain%jsd;jed=Domain%jed
  isc=Domain%isc;iec=Domain%iec;jsc=Domain%jsc;jec=Domain%jec

  
end subroutine get_local_indices
! </SUBROUTINE> NAME="get_local_indices"


!#######################################################################
! <SUBROUTINE NAME="get_domain_offsets">
!
! <DESCRIPTION>
! For getting domain offsets from domain derived type. 
! </DESCRIPTION>
!
subroutine get_domain_offsets(Domain, ioff, joff)

  type(ocean_domain_type), intent(in)  :: Domain
  integer,                 intent(out) :: ioff, joff
  
  ioff = Domain%ioff
  joff = Domain%joff
  
end subroutine get_domain_offsets
! </SUBROUTINE> NAME="get_domain_offsets"


!#######################################################################
! <SUBROUTINE NAME="get_active_indices">
!
! <DESCRIPTION>
! For getting active domain indices from domain derived type. 
! </DESCRIPTION>
!
subroutine get_active_indices(Domain, isa, iea, jsa, jea)

  type(ocean_domain_type), intent(in) :: Domain
  integer, intent(out)                :: isa, iea, jsa, jea

  isa=Domain%isa;iea=Domain%iea;jsa=Domain%jsa;jea=Domain%jea

end subroutine get_active_indices
! </SUBROUTINE> NAME="get_active_indices"


!#######################################################################
! <SUBROUTINE NAME="get_global_indices">
!
! <DESCRIPTION>
! For getting global indices from domain derived type. 
! </DESCRIPTION>
!
subroutine get_global_indices(Domain, isg, ieg, jsg, jeg)

  type(ocean_domain_type), intent(in)  :: Domain
  integer,                 intent(out) :: isg, ieg, jsg, jeg

  isg=Domain%isg;ieg=Domain%ieg;jsg=Domain%jsg;jeg=Domain%jeg
  
end subroutine get_global_indices
! </SUBROUTINE> NAME="get_global_indices"



!#######################################################################
! <SUBROUTINE NAME="reduce_active_domain">
!
! <DESCRIPTION>
! For getting reducing the active domain 
! </DESCRIPTION>
!
subroutine reduce_active_domain(Domain,inc_x,inc_y)

  type(ocean_domain_type), intent(inout)        :: Domain
  integer,                 intent(in), optional :: inc_x, inc_y
  integer                                       :: inc

  if (PRESENT(inc_x)) then
     if (inc_x < 1) call mpp_error(FATAL,'==>Error: active domain decrement must be positive')
     inc = inc_x
  else
     inc = 1
  endif

  Domain%isa = Domain%isa + inc
  Domain%iea = Domain%iea - inc

  if (PRESENT(inc_y)) then
     if (inc_y < 1) call mpp_error(FATAL,'==>Error: active domain decrement must be positive')
     inc = inc_y
  else
     inc = 1
  endif

  Domain%jsa = Domain%jsa + inc
  Domain%jea = Domain%jea - inc

end subroutine reduce_active_domain
! </SUBROUTINE> NAME="reduce_active_domain"



!#######################################################################
! <SUBROUTINE NAME="get_halo_sizes">
!
! <DESCRIPTION>
! For getting halo sizes from domain derived type. 
! </DESCRIPTION>
!
subroutine get_halo_sizes(Domain, xhalo, yhalo)

  type(ocean_domain_type), intent(in)            :: Domain
  integer,                 intent(out), optional :: xhalo, yhalo

  if (PRESENT(xhalo)) xhalo = Domain%xhalo
  if (PRESENT(yhalo)) yhalo = Domain%yhalo

end subroutine get_halo_sizes
! </SUBROUTINE> NAME="get_halo_sizes"

  
end module ocean_domains_mod
