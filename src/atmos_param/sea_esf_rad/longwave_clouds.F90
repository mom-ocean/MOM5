                     module longwave_clouds_mod
 
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  This code calculates longwave cloud radiative parameters, i.e.
!  cloud optical depth, flux, and heating rate.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!

!   shared modules:

use mpp_mod,           only: input_nml_file
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, close_file
use constants_mod,     only: constants_init, radcon

!  shared radiation package modules:

use rad_utilities_mod, only: rad_utilities_init, lw_output_type, &
                             Cldrad_control, lw_clouds_type, &
                             cld_specification_type,  & 
                             cldrad_properties_type

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    longwave_clouds_mod determines lw cloud transmission for each
!    longwave cloud band.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: longwave_clouds.F90,v 19.0 2012/01/06 20:17:31 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public       &
            longwave_clouds_init, &
            cldtau,  cloud,  thickcld, &
            lw_clouds_dealloc,  &
            longwave_clouds_end


!---------------------------------------------------------------------
!---- namelist   -----

logical    :: dummy         = .false.


namelist / longwave_clouds_nml /   &
                                     dummy

!----------------------------------------------------------------------
!--- public data ---------


!----------------------------------------------------------------------
!---   private ---------

integer  :: NLWCLDB
logical  :: module_is_initialized = .false.    ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
!
! <SUBROUTINE NAME="longwave_clouds_init">
!  <OVERVIEW>
!   The constructor method of longwave_clouds module.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This method does the initialization of longwave cloud module. It
!   reads the longwave clouds namelist from input namelist file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_clouds_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_clouds_init

!--------------------------------------------------------------------
!    longwave_clouds_init is the constructor for longwave_clouds_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer               :: unit, ierr, io, logunit

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
      call constants_init
      call rad_utilities_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=longwave_clouds_nml, iostat=io)
      ierr = check_nml_error(io,'longwave_clouds_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_clouds_nml')
        end do
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=longwave_clouds_nml)

!--------------------------------------------------------------------
!    mark the module as initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------


end subroutine longwave_clouds_init



!####################################################################
! <SUBROUTINE NAME="cldtau">
!  <OVERVIEW>
!   Subroutine to calculate cloud optical depth
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine calculates cloud transmission function from cloud
!   emissivity.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cldtau(Cldrad_props, Cld_spec, Lw_clouds)
!  </TEMPLATE>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification input data to cloud optical depth calculation
!  </IN>
!  <INOUT NAME="Lw_clouds" TYPE="Lw_clouds">
!   cloud longwave parameters
!  </INOUT>
! </SUBROUTINE>
!
subroutine cldtau (nprofile, Cldrad_props, Cld_spec, Lw_clouds)
 
!--------------------------------------------------------------------
!    cldtau claculates the cloud transmission function for max overlap
!    and weighted random overlap clouds in each of the lw cloud
!    parameterization bands.  
!--------------------------------------------------------------------

integer,                      intent(in)    :: nprofile
type(cldrad_properties_type), intent(in)    :: Cldrad_props
type(cld_specification_type), intent(in)    :: Cld_spec       
type(lw_clouds_type),         intent(inout) :: Lw_clouds

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Cldrad_props
!     Cld_spec
!
!   intent(inout) variables:
!
!     Lw_clouds
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer  :: is, ie, js, je, ks, ke
      integer   :: n, k, i, j
      integer   :: emiss_index, profile_index

!---------------------------------------------------------------------
!  local variables:
!
!      is,ie,js,je,ks,ke
!      i,j,k,n
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      is = 1
      ie = size(Cld_spec%cmxolw, 1)
      js = 1
      je = size(Cld_spec%cmxolw, 2)
      ks = 1
      ke = size(Cld_spec%cmxolw, 3)
      NLWCLDB = Cldrad_control%NLWCLDB

!--------------------------------------------------------------------
!    allocate the components of the lw_clouds derived type variable.
!--------------------------------------------------------------------
      allocate (Lw_clouds%taucld_rndlw (IS:IE, JS:JE,KS:KE, NLWCLDB) )
      allocate (Lw_clouds%taunbl_mxolw (IS:IE, JS:JE,KS:KE, NLWCLDB) )
      allocate (Lw_clouds%taucld_mxolw (IS:IE, JS:JE,KS:KE, NLWCLDB) )

!----------------------------------------------------------------------
!    define max overlap layer transmission function over layers KS,KE
!----------------------------------------------------------------------
      do n = 1,NLWCLDB
        do k = KS,KE
          Lw_clouds%taucld_mxolw(:,:,k,n) = 1.0E+00 -    &
                                Cldrad_props%emmxolw(:,:,k,n, nprofile)
        end do
      end do
 
!----------------------------------------------------------------------
!    define "weighted random cloud" layer transmission function
!    over layers KS,KE
!----------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
        do n = 1,NLWCLDB
          if (Cldrad_control%do_ica_calcs) then
            profile_index = nprofile
            emiss_index = nprofile
          else
            profile_index = n
            emiss_index = 1
          endif
          do k = KS,KE
            do j = JS,JE
              do i = IS,IE
                if (Cld_spec%crndlw_band(i,j,k,profile_index) >   &
                                                          0.0E+00) then
                  Lw_clouds%taucld_rndlw(i,j,k,n) =   &
                        (Cld_spec%crndlw_band(i,j,k,profile_index)/ &
                        (1.0E+00 - Cld_spec%cmxolw(i,j,k)))*  &
                        (1.0E+00 -   &
                          Cldrad_props%emrndlw(i,j,k,n,emiss_index)) + &
                        1.0E+00 -    &
                          Cld_spec%crndlw_band(i,j,k,profile_index)/   &
                           (1.0E+00 - Cld_spec%cmxolw(i,j,k))
                else
                  Lw_clouds%taucld_rndlw(i,j,k,n) = 1.0E+00
                endif
              end do
            end do
          end do
        end do
      else
        do n = 1,NLWCLDB
          do k = KS,KE
            do j = JS,JE
              do i = IS,IE
                if (Cld_spec%crndlw(i,j,k) > 0.0E+00) then
                  Lw_clouds%taucld_rndlw(i,j,k,n) =   &
                         (Cld_spec%crndlw(i,j,k)/(1.0E+00 - &
                          Cld_spec%cmxolw(i,j,k)))*  &
                    (1.0E+00 - Cldrad_props%emrndlw(i,j,k,n,1)) + &
                         1.0E+00 - Cld_spec%crndlw(i,j,k)/   &
                         (1.0E+00 - Cld_spec%cmxolw(i,j,k))
                else
                  Lw_clouds%taucld_rndlw(i,j,k,n) = 1.0E+00
                endif
              end do
            end do
          end do
        end do
      endif
 
!--------------------------------------------------------------------
!    define "nearby layer" cloud transmission function for max
!    overlapped clouds (if emissivity not equal to one)
!--------------------------------------------------------------------
      do n = 1,NLWCLDB
        do k=KS,KE
          Lw_clouds%taunbl_mxolw(:,:,k,n) = 0.0E+00
        end do
      end do
 
!------------------------------------------------------------------

  
end subroutine cldtau



!######################################################################
! <SUBROUTINE NAME="cloud">
!  <OVERVIEW>
!   Subroutine to calculate cloud transmission function
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine calculates cloud transmission functions above certain
!   level.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud (kl, Cldrad_props, Cld_spec, Lw_clouds, cldtf)
!  </TEMPLATE>
!  <IN NAME="kl" TYPE="integer">
!   the vertical level above which cloud transmission functions are desired.
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification input data to cloud optical depth calculation
!  </IN>
!  <IN NAME="Lw_clouds" TYPE="lw_clouds_type">
!   cloud longwave radiative properties
!  </IN>
!  <OUT NAME="cldtf" TYPE="real">
!   cloud transmission functions
!  </OUT>
! </SUBROUTINE>
!
subroutine cloud (kl, Cldrad_props, Cld_spec,  Lw_clouds, cldtf)

!---------------------------------------------------------------------
!    
!---------------------------------------------------------------------

integer,                      intent(in)  :: kl
type(cldrad_properties_type), intent(in)  :: Cldrad_props
type(cld_specification_type), intent(in)  :: Cld_spec
type(lw_clouds_type),         intent(in)  :: Lw_clouds
real, dimension(:,:,:,:),     intent(out) :: cldtf        

!--------------------------------------------------------------------
!  intent(in) variables;
!
!     kl
!     Cldrad_props
!     Cld_spec 
!     Lw_clouds
!
!  intent(out) variables:
!
!     cldtf
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cld_spec%cmxolw,1), &
                       size(Cld_spec%cmxolw,2), & 
                       size(Cld_spec%cmxolw,3) + 1) ::   &
                                                       cldtfmo, cldtfrd

      integer   :: is, ie, js, je, ks, ke
      integer   :: i, j, kp, n

!---------------------------------------------------------------------
!   local variables:
!
!       cldtfmo
!       cldtfrd
!       is,ie,js,je,ks,ke
!       i,j,kp,n
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
! 
!---------------------------------------------------------------------
      is = 1
      ie = size (Cld_spec%cmxolw, 1)
      js = 1
      je = size(Cld_spec%cmxolw,2)
      ks = 1
      ke = size(Cld_spec%cmxolw, 3)

!---------------------------------------------------------------------
!    the definition of "within a max overlapped cloud" is:
!    at pressure level k (separating layers k and (k-1)), the max
!    overlap cloud amounts for layers k and (k-1) must be 1) nonzero
!    (else no such cloud) and 2) equal (else one such cloud ends at
!    level k and another begins). Another way to define this is: if
!    considering the transmission across layer kp (between levels
!    kp and (kp+1)) the max overlap cloud amounts for layers kp and
!    (kp-1) must be nonzero and equal.
!---------------------------------------------------------------------
      do n = 1,NLWCLDB

!--------------------------------------------------------------------
!   cloud "nearby layer" transmission functions
!--------------------------------------------------------------------
        cldtfmo(:,:,kl) = 0.0
        cldtfrd(:,:,kl) = 1.0

!--------------------------------------------------------------------
!   if level kl is within a maximum overlapped cloud, the cloud
!   "nearby layer" transmission function may be non-unity. Exception:
!   at levels KS,KE+1  the function must be unity.
!--------------------------------------------------------------------
        if (kl > KS .AND. kl < KE+1) then
          do j=JS,JE
            do i=IS,IE
              if (Cld_spec%cmxolw(i,j,kl-1) /= 0.0 .and.     &
                  Cld_spec%cmxolw(i,j,kl) ==   &
                  Cld_spec%cmxolw(i,j,kl-1)) then
                cldtfmo(i,j,kl) = Cld_spec%cmxolw(i,j,kl)*  &
                                  Lw_clouds%taunbl_mxolw(i,j,kl,n)
                cldtfrd(i,j,kl) = 1.0 - Cld_spec%cmxolw(i,j,kl)
              endif
            enddo
          enddo
        endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
        cldtf(:,:,kl,n) = cldtfmo(:,:,kl) + cldtfrd(:,:,kl)

!--------------------------------------------------------------------
!     cloud transmission functions between level kl and higher
!     levels ( when kl le KE+1)
!--------------------------------------------------------------------
        if (kl .LT. KE+1) then
          cldtfmo(:,:,kl) = 0.0
          cldtfrd(:,:,kl) = 1.0

!--------------------------------------------------------------------
!    for first layer below  level kl, assume flux at level kl
!   is unity and is apportioned between (cmxolw) max. overlap cld,
!   (crndlw) rnd overlap cld, and remainder as clear sky.
!--------------------------------------------------------------------
          cldtfmo(:,:,kl+1) = Cld_spec%cmxolw(:,:,kl)*   &
                              Lw_clouds%taucld_mxolw(:,:,kl,n)
          cldtfrd(:,:,kl+1) = (1.0 - Cld_spec%cmxolw(:,:,kl))*  &
                               Lw_clouds%taucld_rndlw(:,:,kl,n)
          cldtf(:,:,kl+1,n) = cldtfmo(:,:,kl+1) + cldtfrd(:,:,kl+1)

!--------------------------------------------------------------------
!    if layers above and below level (kp-1) have no max overlap cloud,
!    or their amounts differ (ie, either top of new max overlap cld or
!    no max overlap cld at all), then apportion total "flux" (or,
!    cloud tf (cldtf)) between any max overlap cloud in layer(kp-1),
!    any rnd overlap cloud and clear sky.
!--------------------------------------------------------------------
          do kp = kl+2, KE+1
            do j=JS,JE
              do i=IS,IE
                if (Cld_spec%cmxolw(i,j,kp-2) .eq. 0. .or.    &
                    Cld_spec%cmxolw(i,j,kp-2) .ne.    &
                    Cld_spec%cmxolw(i,j,kp-1)) then
                  cldtfmo(i,j,kp) = cldtf(i,j,kp-1,n)*   &
                                    Cld_spec%cmxolw(i,j,kp-1)* &
                                    Lw_clouds%taucld_mxolw(i,j,kp-1,n)
                  cldtfrd(i,j,kp) = cldtf(i,j,kp-1,n)*   &
                                    (1.0 - Cld_spec%cmxolw(i,j,kp-1))*&
                                    Lw_clouds%taucld_rndlw(i,j,kp-1,n)
                  cldtf(i,j,kp,n) = cldtfmo(i,j,kp) + cldtfrd(i,j,kp)

!--------------------------------------------------------------------
!    if layer above level (kp-1) has a max overlap cloud, and layer
!    layer below level (kp-1) also does (ie, within max overlap cld)
!    obtain separate cloud tfs for max overlap cloud and for 
!    remainder (which may contain a random overlap cloud).
!--------------------------------------------------------------------
                else 
                  cldtfmo(i,j,kp) = cldtfmo(i,j,kp-1)*   &
                                    Lw_clouds%taucld_mxolw(i,j,kp-1,n)
                  cldtfrd(i,j,kp) = cldtfrd(i,j,kp-1)*   &
                                    Lw_clouds%taucld_rndlw(i,j,kp-1,n)
                  cldtf(i,j,kp,n) = cldtfmo(i,j,kp) + cldtfrd(i,j,kp)
                endif
              end do
            end do
          end do
        endif
      end do

!---------------------------------------------------------------------


end subroutine cloud



!####################################################################

! <SUBROUTINE NAME="thickcld">
!  <OVERVIEW>
!   Subroutine to calculate longwave cloud flux
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine calculates longwave cloud flux at model pressure
!   levels and heating rate.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call thickcld (pflux_in, Cldrad_props, Cld_spec, Lw_output)
!  </TEMPLATE>
!  <IN NAME="pflux_in" TYPE="real">
!   pressure at flux levels of model
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification input data to cloud optical depth calculation
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!   cloud longwave radiative flux
!  </INOUT>
! </SUBROUTINE>
!
subroutine thickcld (pflux_in, Cldrad_props, Cld_spec, Lw_output)

!------------------------------------------------------------------
!    thickcld recomputes cloud fluxes in "thick" clouds assuming
!    that df/dp is constant. the effect is to reduce top-of-cloud
!    cooling rates, thus performing a "pseudo-convective adjustment"
!    by heating (in a relative sense) the cloud top.
!    NOTE: this subroutine cannot handle a frequency-dependent 
!    emissivity. therefore, it assumes that emissivity quantities 
!    (emmxolw) are from frequency band 1 (normally unity).
!---------------------------------------------------------------------

real,   dimension (:,:,:),    intent(in)    ::  pflux_in
type(cldrad_properties_type), intent(in)    :: Cldrad_props
type(cld_specification_type), intent(in)    :: Cld_spec      
type(lw_output_type),         intent(inout) :: Lw_output      

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     pflux   =  pressure at flux levels of model.
!     Cldrad_props
!     Cld_spec
!
!  intent(inout) variables:
!
!     Lw_output
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
 
      real, dimension (size(pflux_in,1),   &
                       size(pflux_in,2))  :: &
                                        delptc, fbtm, ftop, pbtm, ptop

      integer, dimension (size(pflux_in,1),   &
                          size(pflux_in,2))  :: &
                                               itopmxo, ibtmmxo

      integer, dimension (size(pflux_in,1),    &
                          size(pflux_in,2),  &
                          size(pflux_in,3)-1)  :: &
                                                  ktopmxo, kbtmmxo

      real, dimension (size(pflux_in,1),   &
                       size(pflux_in,2),  &
                       size(pflux_in,3)-1)  :: &
                                                  tmp1, pdfinv

      real, dimension (size(pflux_in,1),  &
                       size(pflux_in,2),  &
                       size(pflux_in,3)  )  :: pflux


      integer   :: is, ie, js, je, ks, ke
      integer   ::  kmxolw
      integer   :: i,j, k, kc, kc1, kc2

!---------------------------------------------------------------------
!   local variables:
!
!     delptc
!     fbtm
!     ftop
!     pbtm
!     ptop
!     itopmxo
!     ibtmmxo
!     ktopmxo
!     kbtmmxo
!     tmp1
!     pdfinv 
!     pflux
!     is,ie,js,je,ks,ke
!     kmxolw
!     i,j,k
!     kc
!     kc1,kc2
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      kmxolw = MAXVAL(Cld_spec%nmxolw)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
    is = 1
    ie = size (pflux_in, 1)
    js = 1
    je = size(pflux_in,2)
    ks = 1
    ke = size(pflux_in, 3) - 1

!--------------------------------------------------------------------
!    convert pflux to cgs.
!--------------------------------------------------------------------
      pflux = 10.0*pflux_in
      pdfinv(:,:,ks:ke) = 1.0/(pflux(:,:,ks+1:ke+1) - pflux(:,:,ks:ke))

!--------------------------------------------------------------------
!    determine levels at which max overlap clouds start and stop
!--------------------------------------------------------------------
      itopmxo(:,:) = 0
      ibtmmxo(:,:) = 0
      ktopmxo(:,:,:) = 0
      kbtmmxo(:,:,:) = 0

!--------------------------------------------------------------------
!    max overlap cloud in first layer (not likely)
!--------------------------------------------------------------------
      do j = JS,JE
        do i = IS,IE
          if ( Cld_spec%cmxolw(i,j,KS) .GT. 0.0E+00) then
            itopmxo(i,j) = itopmxo(i,j) + 1
            ktopmxo(i,j,itopmxo(i,j)) = KS
          endif
        end do
      end do

!--------------------------------------------------------------------
!    k-level for which top of max overlap cloud is defined
!--------------------------------------------------------------------
      do k = KS+1,KE
        do j = JS,JE
          do i = IS,IE
            if (Cld_spec%cmxolw(i,j,k) .GT. 0.0E+00 .AND.   &
                Cld_spec%cmxolw(i,j,k-1) /= Cld_spec%cmxolw(i,j,k)) then
              itopmxo(i,j) = itopmxo(i,j) + 1
              ktopmxo(i,j,itopmxo(i,j)) = k
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!    k-level for which bottom of max overlap cloud is defined
!--------------------------------------------------------------------
      do k=KS,KE-1
        do j=JS,JE
          do i=IS,IE
            if (CLd_spec%cmxolw(i,j,k) > 0.0E+00 .AND.    &
                Cld_spec%cmxolw(i,j,k+1) /= Cld_spec%cmxolw(i,j,k)) then
              ibtmmxo(i,j) = ibtmmxo(i,j) + 1
              kbtmmxo(i,j,ibtmmxo(i,j)) = k+1
            endif
          enddo
        enddo
      enddo

!--------------------------------------------------------------------
!    bottom of max overlap cloud in KE'th level
!--------------------------------------------------------------------
      do j = JS,JE
        do i = IS,IE
          if (Cld_spec%cmxolw(i,j,KE) .GT. 0.0E+00) then
            ibtmmxo(i,j) = ibtmmxo(i,j) + 1
            kbtmmxo(i,j,ibtmmxo(i,j)) = KE+1
          endif
        enddo
      enddo
 
!---------------------------------------------------------------------- 
!    obtain the pressures and fluxes of the top and bottom of the cloud.
!---------------------------------------------------------------------- 
      if (kmxolw .NE. 0) then
        do kc=1,kmxolw 
          do j=JS,JE
            do i=IS,IE
              if (kbtmmxo(i,j,kc) > ktopmxo(i,j,kc)) then
                kc1 = ktopmxo(i,j,kc)
                kc2 = kbtmmxo(i,j,kc)
                ptop(i,j) = pflux (i,j,kc1) 
                pbtm(i,j) = pflux (i,j,kc2)
                ftop(i,j) = Lw_output%flxnet(i,j,kc1)
                fbtm(i,j) = Lw_output%flxnet(i,j,kc2)

!-----------------------------------------------------------------------
!      compute the "flux derivative" df/dp delptc.
!-----------------------------------------------------------------------
                delptc(i,j) = (ftop(i,j) - fbtm(i,j))/   &
                              (ptop(i,j) - pbtm(i,j))
!-----------------------------------------------------------------------
!      compute the total flux change from the top of the cloud.
!-----------------------------------------------------------------------
                do k=kc1+1,kc2-1
                  tmp1(i,j,k) = ftop(i,j) + (pflux(i,j,k) - ptop(i,j))*&
                                delptc(i,j) 
                  Lw_output%flxnet(i,j,k) =   &
                        Lw_output%flxnet(i,j,k)*(1.0E+00 -    &
                        Cld_spec%cmxolw(i,j,k)*   &
                        Cldrad_props%emmxolw(i,j,k,1,1)) +  &
                        tmp1(i,j,k)*Cld_spec%cmxolw(i,j,k)*   &
                        Cldrad_props%emmxolw(i,j,k,1,1)
                end do
              endif
            end do
          end do
        end do
      endif

!-----------------------------------------------------------------------
!     recompute the heating rates based on the revised fluxes.
!-----------------------------------------------------------------------
      Lw_output%heatra(:,:,KS:KE) = radcon* &
                                    (Lw_output%flxnet(:,:,KS+1:KE+1) - &
                                     Lw_output%flxnet(:,:,KS:KE))* &
                                     pdfinv(:,:,KS:KE)

!----------------------------------------------------------------------


end subroutine thickcld
  
   
!#####################################################################
! <SUBROUTINE NAME="lw_clouds_dealloc">
!  <OVERVIEW>
!    subroutine to deallocate the array components of the
!    lw_clouds_type variable that is input.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This subroutine deallocates the array components of the
!    lw_clouds_type variable that is input.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_clouds_dealloc (Lw_clouds)
!  </TEMPLATE>
!  <INOUT NAME="Lw_clouds" TYPE="lw_clouds_type">
!   lw_clouds_type variable containing cloud trans-
!   mission function information
!  </INOUT>
! </SUBROUTINE>

subroutine lw_clouds_dealloc (Lw_clouds)

!------------------------------------------------------------------
!    lw_clouds_dealloc deallocates the array components of the
!    lw_clouds_type variable that is input.
!-------------------------------------------------------------------

type(lw_clouds_type), intent(inout)  :: Lw_clouds

!---------------------------------------------------------------------
!  intent(inout) variables:
!
!     Lw_clouds      lw_clouds_type variable containing cloud trans-
!                    mission function information
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the array components of Lw_clouds.
!--------------------------------------------------------------------
      deallocate (Lw_clouds%taucld_rndlw)
      deallocate (Lw_clouds%taucld_mxolw)
      deallocate (Lw_clouds%taunbl_mxolw)

!-------------------------------------------------------------------

end subroutine lw_clouds_dealloc 


!#####################################################################
!
! <SUBROUTINE NAME="longwave_clouds_end">
!  <OVERVIEW>
!   The destructor for longwave_clouds module.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine closes the longwave cloud module. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_clouds_end 
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_clouds_end

!---------------------------------------------------------------------
!    longwave_clouds_end is the destructor for longwave_clouds_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine longwave_clouds_end


!######################################################################


                   end module longwave_clouds_mod

