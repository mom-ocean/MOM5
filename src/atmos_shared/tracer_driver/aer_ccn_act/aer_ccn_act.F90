        module aer_ccn_act_mod

use fms_mod,             only: error_mesg, FATAL, open_namelist_file, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               file_exist, write_version_number, &
                               check_nml_error, close_file
use mpp_mod,             only: input_nml_file, get_unit
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act2_k, &
                               aer_ccn_act_wpdf_k, aer_ccn_act_k_init, &
                               aer_ccn_act_k_end, aer_ccn_act_wpdf_m_k

implicit none
private
    private Loading
      
    public aer_ccn_act, aer_ccn_act2, aer_ccn_act_wpdf, &
           aer_ccn_act_wpdf_m, aer_ccn_act_init, aer_ccn_act_end

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: aer_ccn_act.F90,v 19.0 2012/01/06 20:31:34 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!---------------- private data -------------------


!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

logical  :: nooc = .false.   ! include organic aerosols as ccns ?
real     :: sul_concen = 0.1
real     :: low_concen = 0.1
real     :: high_concen = 1.
 !Parameters for look-up tables
 
real ::  lowup=0.3 !m/s
real ::  highup=10.

! earlier values: lowup2 = 0.001, highmass2 = 1000., highmass3 = 1000.
!real ::  lowup2=0.0001 !m/s
real ::  lowup2=0.01   !m/s
real ::  highup2=0.3
real ::  lowmass2=0.01 !ug m-3
!real ::  highmass2=1000.
real ::  highmass2=100.
real ::  lowmass3=0.01 !ug m-3
!real ::  highmass3=1000.
real ::  highmass3=100.
real ::  lowmass4=0.01 !ug m-3
real ::  highmass4=100.
real ::  lowmass5=0.01 !ug m-3
real ::  highmass5=100.
real :: lowT2=243.15 !K
real :: highT2=308.15

namelist /aer_ccn_act_nml/ nooc, sul_concen, low_concen, high_concen, &
                           lowup, highup, lowup2, highup2, lowmass2, &
                           highmass2, lowmass3, highmass3,  &
                           lowmass4, highmass4, lowmass5, highmass5, &
                           lowT2, highT2


logical :: module_is_initialized  = .false.
 
contains

subroutine aer_ccn_act (T1, P1, Updraft1, TotalMass, Drop)
real, dimension(:), intent(inout) :: TotalMass
real, intent(in) :: T1, P1, Updraft1
real, intent(inout) :: Drop
    
  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()

  tym = size (totalmass,1)

  call aer_ccn_act_k (T1, P1, Updraft1, TotalMass, tym, Drop, ier,  &
                      ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act', ermesg, FATAL)

  
end subroutine aer_ccn_act

subroutine aer_ccn_act2 (T1, P1, Updraft1, TotalMass, mu,airdens,Nc,qc,qt,qe,tc,te,Drop)

!T1 temperature (K)
!P1 pressure (Pa)
!Updraft1 updraft velocity (m/s)
!TotalMass aerosol mass ()
!mu entrainment coef. (/s)
!airdens air density (kg/m3 air)
!Nc droplet mixing ratio (#/kg air)
!qc in-cloud vapor mixing ratio (kg water/kg air)
!qt in-cloud total water mixing ratio qc + ql (kg water/kg air)
!qe environment vapor mixing ratio (kg water/kg air)
!tc in-cloud temperature (K)
!te environment temperature (K)
!Drop droplet number concentration (#/cc)

real, dimension(:), intent(in) :: TotalMass
real, intent(in) :: T1, P1, Updraft1, mu,airdens, Nc, qc, qt, qe, tc, te
real, intent(inout) :: Drop

  integer :: tym, ier
  character(len=256) :: ermesg
        
  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

  call aer_ccn_act2_k (T1, P1, Updraft1, TotalMass, tym, mu,  &
                       airdens,Nc,qc,qt,qe,tc,te,Drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act2', ermesg, FATAL)

end subroutine aer_ccn_act2

!-->cjg: addition
!
! Additional subroutines to compute CCN activation by integrating
! over an assumed subgrid-scale PDF of w

subroutine aer_ccn_act_wpdf(T, p, wm, wp2, totalmass, drop)

! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

real, intent(in)    :: T, p, wm, wp2
real, intent(inout) :: totalmass(4)
real, intent(out)   :: drop

  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

   call aer_ccn_act_wpdf_k (T, p, wm, wp2, totalmass, tym,           &
                            drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act_wpdf', ermesg, FATAL)

end subroutine aer_ccn_act_wpdf

!----------------------------------------------------------------------

subroutine aer_ccn_act_wpdf_m(T, p, wm, wp2, offs, totalmass, drop)


! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

real, intent(in)    :: T, p, wm, wp2
integer, intent(in) :: offs
real, intent(inout) :: totalmass(4)
real, intent(out)   :: drop

  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

  call aer_ccn_act_wpdf_m_k (T, p, wm, wp2, offs, totalmass, tym,       &
                             drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act_wpdf_m', ermesg, FATAL)

end subroutine aer_ccn_act_wpdf_m

!------------------------------------------------------------------------

subroutine aer_ccn_act_init ()

!--------------------------------------------------------------------  
!  local variables:
      
      integer   ::   unit, ierr, io, logunit
      integer, parameter :: res = 20 !
      real, dimension(res,res,res,res,res) :: droplets

      integer, parameter :: res2 = 20 !
      real, dimension(res2,res2,res2,res2,res2) :: droplets2

      if (module_is_initialized) return

!--------------------------------------------------------------------- 
!    read namelist.
!--------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=aer_ccn_act_nml, iostat=io)
        ierr = check_nml_error(io,'aer_ccn_act_nmliostat=io')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aer_ccn_act_nml, iostat=io,  &
               end=10)
        ierr = check_nml_error(io,'aer_ccn_act_nml')
        end do
10      call close_file (unit)   
#endif
      endif                      
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!--------------------------------------------------------------------
       call write_version_number (version, tagname)
       logunit=stdlog()
       if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=aer_ccn_act_nml)

       call Loading( droplets, droplets2)

       call aer_ccn_act_k_init (droplets,   &
                      droplets2, res, res2, nooc,  &
                       sul_concen, low_concen, high_concen, &
                       lowup, highup, lowup2, highup2, lowmass2, &
                       highmass2, lowmass3, highmass3,  &
                       lowmass4, highmass4, lowmass5, highmass5, &
                      lowT2, highT2  )
       module_is_initialized  = .true.

end subroutine aer_ccn_act_init


subroutine Loading(droplets, droplets2)

real, dimension(:,:,:,:,:), intent(out) :: droplets, droplets2
real xx
integer i, j, k, l, m, unit
integer res, res2

  res = size(droplets,1)
  res2 = size(droplets2,1)
  unit = get_unit()
  open(unit, FILE='INPUT/droplets.dat')
  do k=1,res
    do i=1,res
      do j=1, res
        do l=1, res
        do m=1, res
          read(unit,*) xx
          droplets(m,l,j,i,k)=xx
        end do
        end do
      end do
    end do
  end do
  close(unit)

  unit = get_unit()
  open(unit, FILE='INPUT/droplets2.dat')
  do k=1,res2
    do i=1,res2
      do j=1, res2
        do l=1, res2
        do m=1, res2
          read(unit,*) xx
          droplets2(m,l,j,i,k)=xx
        end do
        end do
      end do
    end do
  end do
  close(unit)

end subroutine Loading


subroutine aer_ccn_act_end()

  call aer_ccn_act_k_end 
  module_is_initialized  = .false.

end subroutine aer_ccn_act_end

end module aer_ccn_act_mod
