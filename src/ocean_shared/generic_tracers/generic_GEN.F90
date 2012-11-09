!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL=""> 
! </REVIEWER>
!
!<OVERVIEW>
! This module contains the generic version of "Generic" Tracers.
! It is a very light-weight module, with minimal FMS dependence.
! It is designed to run with MOM in applications, such as GEOS-5,
! where MOM is being hosted by a non-fms system. It does not depend
! on any ESMF or GEOS-5 concepts. It assumes that all 
! internal chemistry or biology calculations are done by the host 
! system and that the host system has already added all other boundary
! sources. It also assumes that the host system is in charge of initialization
! of the tracers. MOM is only responsible for transport calculations.
! Note that, because of MOM's allocation strategy, it is not currently possible
! to fully share pointers; therefore, using this module may require 
! having multiple copies of the tracers.
!
! The module's only public method is generic_GEN_register, which must
! be called by the generic tracer handler after it has read the namelist.
! The only MOM dependencies are g_tracer_type and the g_tracer_add method.
!
!</OVERVIEW>
!
!<DESCRIPTION>
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! </REFERENCE>
!
! <DEVELOPER_NOTES>
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_GEN

  use mpp_mod,           only : stdout
  use g_tracer_utils,    only : g_tracer_type,  g_tracer_add

  implicit none
  private

  public do_generic_GEN, Num_generic_GEN_Tracers
  public generic_GEN_register

  logical, save :: do_generic_GEN  = .false.
  integer, save :: Num_generic_GEN_Tracers = 0

contains

  subroutine generic_GEN_register(tracer_list)

    type(g_tracer_type), pointer :: tracer_list

    integer :: i

    character(len=13) :: name="GENtracer_xxx"

    if(Num_generic_GEN_Tracers>0) then

       write (stdout(),'(/,A,i4,A)')  "Will register ",Num_generic_GEN_Tracers, &
            " genericGEN tracers "
       do i=1,Num_generic_GEN_Tracers
          write(name(11:13),'(I3.3)') i
          call g_tracer_add(tracer_list,'generic_gen',&
               name       = name,                &
               longname   = name//"_longname",   &
               units      = name//"_units",      &
               prog       = .true.)

          write (stdout(),'(/,A)') &
               "  Registered genericGEN tracer -> "//trim(name)
       end do
    else
          write (stdout(),'(/,A,I4)') &
               "  Something strange--Trying to register genericGEN tracers, but Num_generic_GEN_Tracers is ", Num_generic_GEN_Tracers               
    end if
  end subroutine generic_GEN_register

end module generic_GEN
