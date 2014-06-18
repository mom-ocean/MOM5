module detr_ice_num_mod

use mg_const_mod, only : mg_const_init, rhoi 
use mpp_mod,      only : input_nml_file
use fms_mod,      only : mpp_pe, mpp_root_pe, file_exist, stdlog, &
                         open_namelist_file, check_nml_error, close_file, &
                         write_version_number, error_mesg, FATAL

implicit none
private

!-----------------------------------------------------------------------
!----public interfaces-------------------------------------------------

public  detr_ice_num,  detr_ice_num_init,  detr_ice_num_end


!-----------------------------------------------------------------------
!----version number-----------------------------------------------------
character(len=128) :: version = '$Id: detr_ice_num.F90,v 19.0 2012/01/06 20:10:40 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

!------------------------------------------------------------------------
!--namelist--------------------------------------------------------------

real :: ice_detr_tm = 213.15 !  temperature threshold for constraining 
                             !  max. detrained ice num
real :: ice_detr_rv = 25.e-6 !  mean volume radius for detrained ice 
                             !  used for ice_detr_opt = 2
character(len=64) :: detr_ice_radius = 'kristjansson'  
                             !  ice size assumption to be used for the 
                             !  detrained ice number calculation
                             !  valid values: 'kristjansson', 'constant'

namelist / detr_ice_num_nml /  ice_detr_tm, ice_detr_rv, detr_ice_radius

!-----------------------------------------------------------------------

logical            :: use_kristjansson_radius = .false.
logical            :: use_constant_radius = .false.
logical            :: module_is_initialized = .false.



contains


!##########################################################################

SUBROUTINE detr_ice_num_init

      integer :: unit, io, ierr, logunit

      IF (module_is_initialized) return

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=detr_ice_num_nml, iostat=io)
      ierr = check_nml_error(io,'detr_ice_num_nml')
#else
      if ( file_exist('input.nml')) then
 
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=detr_ice_num_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'detr_ice_num_nml')
        enddo
10      call close_file (unit)
      endif
#endif

!--------- write version and namelist to standard log ------------

        call write_version_number ( version, tagname )
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
                      write ( logunit, nml=detr_ice_num_nml )

!--------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!--------------------------------------------------------------------------
      call mg_const_init

!--------------------------------------------------------------------------
!    check for valid values.
!--------------------------------------------------------------------------
      if (trim(detr_ice_radius)  == 'kristjansson') then
        use_kristjansson_radius = .true.
      else if (trim(detr_ice_radius) == 'constant') then   
        use_constant_radius = .true.
      else
        call error_mesg ('detr_ice_num_mod', &
             'value of nml variable detr_ice_radius not acceptable', FATAL)
      END IF

!------------------------------------------------------------------------
!    mark module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------


END SUBROUTINE  detr_ice_num_init



!#########################################################################

SUBROUTINE detr_ice_num (tin, dqi, dni)

!-------------------------------------------------------------------------

real, intent (in),     dimension(:,:,:) :: tin, dqi
real, intent (out),    dimension(:,:,:) :: dni

!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      real    :: dl, asp, R_equiv
      integer :: i, j, k

!-------------------------------------------------------------------------
!    check for initialization.
!-------------------------------------------------------------------------
      if (.not. module_is_initialized) call  detr_ice_num_init

!-------------------------------------------------------------------------
!    calculate number of ice particles detrained given mass of ice detrained
!    and mean particle size. two options are given. rhoi is density of ice 
!    in morrison-gettleman microphysics.
!-------------------------------------------------------------------------
      IF (use_kristjansson_radius) THEN

!-------------------------------------------------------------------------
!    option 1:
!    compute L(T) from Kristjansson et al., JGR, 2000, Eq. 3 based on 
!    Ryan, Mc Farquhar and Heymsfield . assume hexagonal ice: 
!    Volume = 3*SQRT(3) / 8 * D^2 * L  (as Fu + Liou 1992)
!    with aspect ratio a = D/L (a from Fu, 1996, JClim)
!              => R_eqiv = 3/(4 !Pi)  3*SQRT(3) / 8 ( a^2 * L^3 ) 
!
!    this could be replaced by lookup table
!-------------------------------------------------------------------------
        do k=1,size(tin,3) 
          do j=1,size(tin,2)
            do i=1,size(tin,1)
              if (dqi(i,j,k) /= 0.0) then
                dl = 1030.7*  &
                         exp(0.05522*(MAX(tin(i,j,k),ice_detr_tm) - 279.5))
                IF ( dl .LE. 30. )  then
                  asp = 1.
                else IF ( dl .GT. 30. .AND. dl .LE. 80. )  then
                  asp =0.8
                else IF ( dl .GT. 80. .AND. dl .LE. 200. )  then
                  asp=0.5
                else IF ( dl .GT. 200. .AND. dl .LE. 500. ) then
                  asp= 0.34      
                else IF ( dl .GT. 500. ) then
                  asp=0.22
                endif
                R_equiv = 1.E-6*  &
                        (3./(4.*3.14)*3.*SQRT(3.)/8.*asp**2*dl**3 )**(1./3.)
                dni(i,j,k) = dqi(i,j,k)/rhoi*3./(4.*3.14*R_equiv**3 )
              else
                dni(i,j,k) = 0.
              endif
            end do
          end do
        end do
      ELSE IF (use_constant_radius) THEN

!-------------------------------------------------------------------------
!    option 2:
!    assume ice_detr_rv micron mean volume radius for detrained ice
!-------------------------------------------------------------------------
        do k=1,size(tin, 3) 
          do j=1,size(tin, 2)
            do i=1,size(tin, 1)
              dni(i,j,k) = dqi(i,j,k)/rhoi*3./(4.*3.14*ice_detr_rv**3)
            end do
          end do
        end do
      END IF

!-------------------------------------------------------------------------


END SUBROUTINE detr_ice_num



!#######################################################################

 SUBROUTINE  detr_ice_num_end

!------------------------------------------------------------------------
!    mark the module as uninitialized.
!------------------------------------------------------------------------
      module_is_initialized = .FALSE.

!-----------------------------------------------------------------------


 END SUBROUTINE  detr_ice_num_end




!########################################################################



end module detr_ice_num_mod
