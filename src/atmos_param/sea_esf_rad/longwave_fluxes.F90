                  module longwave_fluxes_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Fei Liu
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Dan Schwartzkopf
! </REVIEWER>
! <OVERVIEW>
!  This code is a helper module that provides various operations on 
!  longwave flux variables.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!
!  shared modules:

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, fms_init, &
                              mpp_pe, mpp_root_pe, stdlog, &
                              file_exist, write_version_number, &
                              check_nml_error, error_mesg, &
                              FATAL, close_file

!  shared radiation package modules:

use rad_utilities_mod, only:  Rad_control, &
                              rad_utilities_init, lw_diagnostics_type

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    longwave_fluxes calculates the longwave fluxes between each model
!    level and every other modle level for each of the longwave
!    spectral parameterization bands.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: longwave_fluxes.F90,v 19.0 2012/01/06 20:15:17 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!-------  interfaces --------

public    &
       longwave_fluxes_init, &
       longwave_fluxes_ks, longwave_fluxes_k_down, &
       longwave_fluxes_KE_KEp1, longwave_fluxes_diag, &
       longwave_fluxes_sum, longwave_fluxes_end



!---------------------------------------------------------------------
!-------- namelist  ---------

real      ::  dummy = 1.0                    

namelist / longwave_fluxes_nml /        &
                                    dummy

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized = .false. ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains

! <SUBROUTINE NAME="longwave_fluxes_init">
!  <OVERVIEW>
!   Subroutine to initialize longwave fluxes namelist
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize longwave fluxes namelist
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_fluxes_init 

!---------------------------------------------------------------------
!    longwave_fluxes_init is the constructor for longwave_fluxes_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
     integer    ::  unit, ierr, io, logunit

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=longwave_fluxes_nml, iostat=io)
      ierr = check_nml_error(io,"longwave_fluxes_nml")
#else
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_fluxes_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_fluxes_nml')
        end do
10      call close_file (unit)
      endif
#endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=longwave_fluxes_nml)

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

end subroutine longwave_fluxes_init



!#####################################################################
! <SUBROUTINE NAME="longwave_fluxes_ks">
!  <OVERVIEW>
!   Subroutine to calculate longwave diagnostic fluxes
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate longwave diagnostic fluxes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_ks ( source, trans, source2, trans2,  &
!                             cld_trans, cld_ind, Lw_diagnostics)
!  </TEMPLATE>
!  <IN NAME="source" TYPE="real">
!   source is longwave source function.
!  </IN>
!  <IN NAME="trans" TYPE="real">
!   trans is longwve transmittance function
!  </IN>
!  <IN NAME="source2" TYPE="real">
!   source2 is longwave source function
!  </IN>
!  <IN NAME="trans2" TYPE="real">
!   trans2 is longwve transmittance function
!  </IN>
!  <IN NAME="cld_trans" TYPE="real">
!   cld_trans is longwave cloud transmittance function
!  </IN>
!  <IN NAME="cld_ind" TYPE="real">
!   cld_ind is a lookup table to translate longwave band index to cloud index
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   Lw_diagnostics contains the longwave diagnostics flux values
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_fluxes_ks (source, trans, source2, trans2,  &
                               cld_trans, cld_ind, Lw_diagnostics)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
integer, dimension(:),      intent(in)    :: cld_ind
real, dimension (:,:,:,:),  intent(in)    :: source
real, dimension (:,:,:,:),  intent(in)    :: source2
real, dimension (:,:,:,:),  intent(in)    :: trans2, trans
real, dimension (:,:,:,:),  intent(in)    :: cld_trans
type(lw_diagnostics_type),  intent(inout) :: Lw_diagnostics
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     cld_ind
!     source
!     source2
!     trans
!     trans2
!     cld_trans
!
!  intent(inout) variables:
!
!     Lw_diagnostics
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension (size(source2,1), &
                       size(source2,2), &
                       size(source2,3) ) ::    flux_tmp, flux_tmp2

      integer   ::   k, ks, ke, nbands, m

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      flux_tmp2
!      k
!      ks
!      ke
!      nbands
!      m
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ks =1
      ke = size(source2,3)-1
      nbands = size(source,4)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m = 1, nbands
        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,KS,m)*trans(:,:,k,m    )
          flux_tmp2(:,:,k) = source2(:,:,k,m)*trans2(:,:,k ,m    )
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if ((m  == 1) .or. (m  >= 7)) then
          Lw_diagnostics%fluxn(:,:,KS,m) =    &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     source(:,:,KS,m)*trans(:,:,KS,m)
        else
          Lw_diagnostics%fluxn(:,:,KS,m) =    &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     source(:,:,KS,m)
        endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) =   &
                                     Lw_diagnostics%fluxn(:,:,k,m) +  &
                                     flux_tmp(:,:,k)* & 
                                     cld_trans(:,:,k,cld_ind(m))
          Lw_diagnostics%fluxn(:,:,KS,m) =  &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     flux_tmp2(:,:,k)*   &
                                     cld_trans(:,:,k,cld_ind(m))
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          if ((m  == 1) .or. (m  >= 7)) then
            Lw_diagnostics%fluxncf(:,:,KS,m) = source(:,:,KS,m)*   &
                                               trans(:,:,KS,m)
          else
            Lw_diagnostics%fluxncf(:,:,KS,m) =  source(:,:,KS,m)
          endif
          do k=KS+1,KE+1
            Lw_diagnostics%fluxncf(:,:,k,m) =   &
                                    Lw_diagnostics%fluxncf(:,:,k,m) + &
                                    flux_tmp(:,:,k)
            Lw_diagnostics%fluxncf(:,:,KS,m) =  &
                                    Lw_diagnostics%fluxncf(:,:,KS,m) + &
                                    flux_tmp2(:,:,k)
          end do
        endif
     end do   ! (m loop)

!---------------------------------------------------------------------


end subroutine longwave_fluxes_ks



!####################################################################
! <SUBROUTINE NAME="longwave_fluxes_k_down">
!  <OVERVIEW>
!   Subroutine to calculate longwave diagnostic fluxes
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate longwave diagnostic fluxes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_k_down (klevel, source, trans, trans2,   &
!                                cld_trans, cld_ind,   Lw_diagnostics)
!  </TEMPLATE>
!  <IN NAME="klevel" TYPE="integer">
!   klevel is the starting vertical level to calculate longwave fluxes
!  </IN>
!  <IN NAME="source" TYPE="real">
!   source is longwave flux source function
!  </IN>
!  <IN NAME="trans" TYPE="real">
!   trans is longwave flux transmittance function
!  </IN>
!  <IN NAME="trans2" TYPE="real">
!   trans2 is longwave flux transmittance function
!  </IN>
!  <IN NAME="cld_trans" TYPE="real">
!   cld_trans is longwave cloud transmittance function
!  </IN>
!  <IN NAME="cld_ind" TYPE="real">
!   cld_ind is a lookup table to translate longwave band index to cloud index
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   Lw_diagnostics contains the longwave diagnostics flux values
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_fluxes_k_down (klevel, source, trans, trans2,   &
                                   cld_trans, cld_ind, Lw_diagnostics)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

integer,                      intent(in)     ::  klevel
real,    dimension (:,:,:,:), intent(in)     ::  source, trans, &
                                                 trans2, cld_trans
type(lw_diagnostics_type),    intent(inout)  ::  Lw_diagnostics
integer, dimension(:),        intent(in)     ::  cld_ind

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     klevel
!     source
!     trans
!     trans2
!     cld_trans
!     cld_ind
!
!  intent(inout) variables:
!
!     Lw_diagnostics
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension (size(source,1), size(source,2)) :: flux4, flux4a

      real    ::  flux_tmp, flux_tmp2
      integer ::  kp, i, j, israd, ierad, jsrad, jerad
      integer :: ke
      integer :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux4
!      flux4a
!      flux3a
!      flux_tmp
!      flux_tmp2
!      kp
!      i,j,k
!      nn
!      ntot
!      israd,ierad
!      jsrad,jerad
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ierad  = size(source,1)
      jerad  = size(source,2)
      israd  = 1
      jsrad  = 1
      ke     = size(source,3)-1
      nbands = size(trans,4)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m=1, nbands
        do kp=klevel+1,KE+1
          do j=jsrad,jerad
            do i=israd,ierad   
              flux_tmp = source(i,j,klevel,m)*trans(i,j,kp,m)
              Lw_diagnostics%fluxn(i,j,kp,m) =    &
                          Lw_diagnostics%fluxn(i,j,kp,m) + flux_tmp*   &
                          cld_trans(i,j,kp, cld_ind(m))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
              if (Rad_control%do_totcld_forcing) then
                Lw_diagnostics%fluxncf(i,j,kp,m) =   &
                            Lw_diagnostics%fluxncf(i,j,kp,m) + flux_tmp 

              endif
            end do
          end do
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        flux4(:,:)  = 0.0
        flux4a(:,:) = 0.0
        do kp=klevel+1,KE+1
          do j=jsrad,jerad
            do i=israd,ierad   
              flux_tmp2 = source(i,j,kp,m)*trans2(i,j,kp,m)
              flux4(i,j) = flux4(i,j) + flux_tmp2*  &
                           cld_trans(i,j,kp, cld_ind(m))
              if (Rad_control%do_totcld_forcing) then
                flux4a (i,j) = flux4a (i,j) + flux_tmp2        
              endif
            end do
          end do
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do j=jsrad,jerad
          do i=israd,ierad   
            Lw_diagnostics%fluxn  (i,j,klevel,m) =   &
                               Lw_diagnostics%fluxn  (i,j,klevel,m) +  &
                               flux4    (i,j       )
            if (Rad_control%do_totcld_forcing) then
              Lw_diagnostics%fluxncf(i,j,klevel,m) =  &
                               Lw_diagnostics%fluxncf(i,j,klevel,m) +  &
                               flux4a   (i,j       )
            endif
          end do
        end do
      end do  ! (nbands loop)

!---------------------------------------------------------------------


end subroutine longwave_fluxes_k_down



!####################################################################
! <SUBROUTINE NAME="longwave_fluxes_KE_KEp1">
!  <OVERVIEW>
!   Subroutine to calculate longwave diagnostic fluxes
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate longwave diagnostic fluxes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_KE_KEp1 (source, trans, trans2,   &
!                                cld_trans, cld_ind,   Lw_diagnostics)
!  </TEMPLATE>
!  <IN NAME="source" TYPE="real">
!   source is longwave flux source function
!  </IN>
!  <IN NAME="trans" TYPE="real">
!   trans is longwave flux transmittance function
!  </IN>
!  <IN NAME="trans2" TYPE="real">
!   trans2 is longwave flux transmittance function
!  </IN>
!  <IN NAME="cld_trans" TYPE="real">
!   cld_trans is longwave cloud transmittance function
!  </IN>
!  <IN NAME="cld_ind" TYPE="real">
!   cld_ind is a lookup table to translate longwave band index to cloud index
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   Lw_diagnostics contains the longwave diagnostics flux values
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_fluxes_KE_KEp1 (source, trans, trans2, cld_trans,&
                                       cld_ind, Lw_diagnostics)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real,    dimension (:,:,:,:),   intent(in)    :: source, cld_trans
real,    dimension (:,:,:),     intent(in)    :: trans, trans2
integer, dimension(:),          intent(in)    :: cld_ind
type(lw_diagnostics_type),      intent(inout) :: Lw_diagnostics

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     source
!     cld_trans
!     trans
!     trans2
!     cld_ind
!
!  intent(inout) variables:
!
!     Lw_diagnostics
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(trans,1), size(trans,2)) ::  &
                                                   flux_tmp, flux_tmp2

      integer :: ke
      integer :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      flux_tmp2
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ke     = size(source,3) - 1
      nbands = size(trans,3) 

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m=1,nbands
        flux_tmp(:,:) = source(:,:,KE+1,m)*trans(:,:,m)
        flux_tmp2(:,:) = source(:,:,KE,m)*trans2(:,:,m)
        Lw_diagnostics%fluxn(:,:,KE,m) =    &
                       Lw_diagnostics%fluxn(:,:,KE,m) + &
                       flux_tmp(:,:)*cld_trans(:,:,KE+1,cld_ind(m))
        Lw_diagnostics%fluxn(:,:,KE+1,m) =  &
                       Lw_diagnostics%fluxn(:,:,KE+1,m) +   &
                       flux_tmp2(:,:)*cld_trans(:,:,KE+1,cld_ind(m))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          Lw_diagnostics%fluxncf(:,:,KE,m) =   &
                       Lw_diagnostics%fluxncf(:,:,KE,m) + flux_tmp(:,:)
          Lw_diagnostics%fluxncf(:,:,KE+1,m) =  &
                       Lw_diagnostics%fluxncf(:,:,KE+1,m) +  &
                       flux_tmp2(:,:)
        endif
      end do  ! (nbands loop)

!---------------------------------------------------------------------

end subroutine longwave_fluxes_KE_KEp1



!####################################################################
! <SUBROUTINE NAME="longwave_fluxes_diag">
!  <OVERVIEW>
!   Subroutine to calculate longwave diagnostic fluxes
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to calculate longwave diagnostic fluxes
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_diag (source, trans,  &
!                               cld_trans, cld_ind,   Lw_diagnostics)
!  </TEMPLATE>
!  <IN NAME="source" TYPE="real">
!   source is longwave flux source function
!  </IN>
!  <IN NAME="trans" TYPE="real">
!   trans is longwave flux transmittance function
!  </IN>
!  <IN NAME="cld_trans" TYPE="real">
!   cld_trans is longwave cloud transmittance function
!  </IN>
!  <IN NAME="cld_ind" TYPE="real">
!   cld_ind is a lookup table to translate longwave band index to cloud index
!  </IN>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   Lw_diagnostics contains the longwave diagnostics flux values
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_fluxes_diag (source, trans, cld_trans, cld_ind, &
                                 Lw_diagnostics)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
real, dimension (:,:,:,:), intent(in)    :: source, trans, cld_trans
integer, dimension(:),     intent(in)    :: cld_ind
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
!-------------------------------------------------------------------
!  intent(in) variables:
!
!     source
!     trans
!     cld_trans
!     cld_ind
!
!  intent(inout) variables:
!
!     Lw_diagnostics
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      real, dimension (size(trans,1), &
                       size(trans,2), &
                       size(trans,3)) ::   flux_tmp

      integer   :: k, ks, ke
      integer   :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      k
!      ks
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ks     = 1
      ke     = size(trans,3) - 1
      nbands = size(trans,4)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m=1,nbands

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,k,m)*trans(:,:,k,m)
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) =   &
                             Lw_diagnostics%fluxn(:,:,k,m) +  &
                             flux_tmp(:,:,k)*cld_trans(:,:,k,cld_ind(m))
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          do k=KS+1,KE+1
            Lw_diagnostics%fluxncf(:,:,k,m) =    &
                               Lw_diagnostics%fluxncf(:,:,k,m) +   &
                               flux_tmp(:,:,k)
          end do
        endif
      end do ! (m loop)

!---------------------------------------------------------------------



end subroutine longwave_fluxes_diag




!###################################################################
! <SUBROUTINE NAME="longwave_fluxes_sum">
!  <OVERVIEW>
!   Subroutine to compute summation of diagnostic longwave fluxes over
!   all bands
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute summation of diagnostic longwave fluxes over
!   all bands
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_fluxes_sum (is, ie, js, je, flux, NBTRGE,         &
!                             Lw_diagnostics, fluxcf)
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   Obsolete
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   Obsolete
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   Obsolete
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   Obsolete
!  </IN>
!  <OUT NAME="flux" TYPE="real">
!   all sky total longwave flux
!  </OUT>
!  <IN NAME="NBTRGE" TYPE="integer">
!   number of longwave flux bands 
!  </IN>
!  <IN NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   longwave flux diagnostics
!  </IN>
!  <OUT NAME="fluxcf" TYPE="real">
!   clear sky total longwave flux
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_fluxes_sum (is, ie, js, je, flux, nbtrge,         &
                                Lw_diagnostics, fluxcf)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

integer,                          intent(in)    :: is, ie, js, &
                                                   je, nbtrge
real, dimension(:,:,:),           intent(out)   :: flux
real, dimension(:,:,:), optional, intent(out)   :: fluxcf
type(lw_diagnostics_type),        intent(in)    :: Lw_diagnostics

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     is,ie,js,je
!     nbtrge
!     Lw_diagnostics
!
!  intent(out) variables:
!
!     flux
!     fluxcf
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
    integer       ::   m

!---------------------------------------------------------------------
!  local variables:
!
!      j,m        
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      flux = 0.
      do m= 1, 6+NBTRGE               
        flux(:,:,:) = flux(:,:,:) + Lw_diagnostics%fluxn(:,:,:,m)
      end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then 
        fluxcf = 0.
        do m= 1, 6+NBTRGE               
          fluxcf(:,:,:) = fluxcf(:,:,:) +    &
                          Lw_diagnostics%fluxncf(:,:,:,m)
        end do
      endif

!--------------------------------------------------------------------


end subroutine longwave_fluxes_sum


!#####################################################################

subroutine longwave_fluxes_end

!--------------------------------------------------------------------
!    longwave_fluxes_end is the destructor for the longwave_fluxes_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine longwave_fluxes_end


!#####################################################################


                end module longwave_fluxes_mod
