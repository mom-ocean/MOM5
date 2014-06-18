                     module longwave_driver_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Code to set up longwave radiation calculation
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
! 

use mpp_mod,            only: input_nml_file
use fms_mod,            only: open_namelist_file, fms_init, &
                              mpp_pe, mpp_root_pe, stdlog, &
                              file_exist, write_version_number, &
                              check_nml_error, error_mesg, &
                              FATAL, close_file
use time_manager_mod,   only: time_type

! shared radiation package modules:

use rad_utilities_mod,  only: rad_utilities_init, Rad_control, &
                              cldrad_properties_type, &
                              cld_specification_type, lw_output_type, &
                              atmos_input_type, radiative_gases_type, &
                              aerosol_type, aerosol_properties_type,  &
                              aerosol_diagnostics_type, &
                              Lw_control, assignment(=), &
                              lw_table_type, lw_diagnostics_type

!   radiation package module:

use sealw99_mod,        only: sealw99_init,sealw99_time_vary, sealw99, &
                              sealw99_endts, sealw99_end

!------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    longwave_driver_mod is the driver for the longwave radiation
!    component of the sea_esf_rad radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: longwave_driver.F90,v 19.0 2012/01/06 20:18:03 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'

!---------------------------------------------------------------------
!-------  interfaces --------

public      &
   longwave_driver_init, longwave_driver_time_vary, longwave_driver,   &
   longwave_driver_endts, longwave_driver_end

private      &

! called from longwave_driver:
         longwave_driver_alloc


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16) :: lwform= 'sealw99'
 

namelist / longwave_driver_nml /    &
                                 lwform

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized =  .false.   ! module initialized ?
logical :: do_sealw99 = .false.               ! sealw99 parameter-
                                              ! ization active ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                    PUBLIC SUBROUTINES
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="longwave_driver_init">
!  <OVERVIEW>
!   longwave_driver_init is the constructor for longwave_driver_mod
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine initializes longwave radiation package
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_init (latb, lonb, pref, Lw_tables)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   2d array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   2d array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   array containing two reference pressure profiles [pascals]
!  </IN>
!  <INOUT NAME="Lw_tables" TYPE="lw_table_type">
!   lw_tables_type variable containing various longwave
!                 table specifiers needed by radiation_diag_mod.
!  </INOUT>
! </SUBROUTINE>
!
subroutine longwave_driver_init (latb, lonb, pref, Lw_tables)
 
!---------------------------------------------------------------------
!    longwave_driver_init is the constructor for longwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:),   intent(in)    :: latb, lonb
real, dimension(:,:),   intent(in)    :: pref
type(lw_table_type),    intent(inout) :: Lw_tables

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      2d array of model latitudes at cell corners 
!                 [ radians ]
!       lonb      2d array of model longitudes at cell corners 
!                 [ radians ]
!       pref      array containing two reference pressure profiles 
!                 [ Pa ]
!
!  intent(out) variables:
!
!       Lw_tables lw_tables_type variable containing various longwave
!                 table specifiers needed by radiation_diag_mod.
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables

      integer     :: unit, ierr, io, logunit

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

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=longwave_driver_nml, iostat=io)
      ierr = check_nml_error(io,'longwave_driver_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_driver_nml')
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
                          write (logunit, nml=longwave_driver_nml)

!--------------------------------------------------------------------
!    determine if valid specification of lw radiation has been made.
!    if optional packages are provided at some later time, this is where
!    the choice of package will be made.
!---------------------------------------------------------------------
      if (trim(lwform) == 'sealw99') then
        do_sealw99 = .true.
        call sealw99_init ( latb, lonb, pref, Lw_tables)
      else
        call error_mesg ( 'longwave_driver_mod', &
                 'invalid longwave radiation form specified', FATAL)
      endif

!---------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------


end  subroutine longwave_driver_init


!#####################################################################

subroutine longwave_driver_time_vary (Time, Rad_gases_tv)
 
!-------------------------------------------------------------------- 
type(time_type), intent(in) :: Time
type(radiative_gases_type), intent(inout) :: Rad_gases_tv

 

      call sealw99_time_vary (Time, Rad_gases_tv)
 
end subroutine longwave_driver_time_vary      

        
 
!#####################################################################

subroutine longwave_driver_endts (Rad_gases_tv)
          
type(radiative_gases_type), intent(in) :: Rad_gases_tv
 

     call sealw99_endts (Rad_gases_tv)


end subroutine longwave_driver_endts 



!#####################################################################
! <SUBROUTINE NAME="longwave_driver">
!  <OVERVIEW>
!   Subroutine to set up and execute longwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver (is, ie, js, je, Atmos_input, Rad_gases, &
!                         Aerosol, Cldrad_props, Cld_spec, Lw_output, &
!                         Lw_diagnostics)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!   starting subdomain i indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   ending subdomain i indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   starting subdomain j indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   ending subdomain j indice of data in the physics_window being
!       integrated
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!  </IN>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </INOUT>
!  <INOUT NAME="Lw_diagnostics" TYPE="lw_diagnostics_type">
!   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  </INOUT>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol input data to longwave radiation
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud specification input data to longwave radiation
!  </IN>
! </SUBROUTINE>
!
subroutine longwave_driver (is, ie, js, je, Rad_time, Atmos_input,  &
                            Rad_gases, Aerosol, Aerosol_props,   &
                            Cldrad_props, Cld_spec, Aerosol_diags, &
                            Lw_output, Lw_diagnostics)

!--------------------------------------------------------------------
!    longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!--------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
type(time_type),              intent(in)     :: Rad_time
type(atmos_input_type),       intent(in)     :: Atmos_input  
type(radiative_gases_type),   intent(inout)  :: Rad_gases   
type(aerosol_type),           intent(in)     :: Aerosol     
type(aerosol_properties_type),intent(inout)  :: Aerosol_props
type(aerosol_diagnostics_type),intent(inout)  :: Aerosol_diags
type(cldrad_properties_type), intent(in)     :: Cldrad_props
type(cld_specification_type), intent(in)     :: Cld_spec     
type(lw_output_type), dimension(:),  intent(inout)  :: Lw_output
type(lw_diagnostics_type),    intent(inout)  :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Rad_time       time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol 
!                     fields that are seen by the longwave radiation 
!                     package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(inout) variables:
!
!      Aerosol_props  aerosol_properties_type variable containing the 
!                     aerosol radiative properties needed by the rad-
!                     iation package 
!      Lw_output      lw_output_type variable containing longwave 
!                     radiation output data 
!      Lw_diagnostics lw_diagnostics_type variable containing diagnostic
!                     longwave output used by the radiation diagnostics
!                     module
!  
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      type(lw_output_type)  :: Lw_output_std, Lw_output_ad
      logical :: calc_includes_aerosols
      integer  :: ix, jx, kx  ! dimensions of current physics window

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_driver_mod',   &
          'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    call longwave_driver_alloc to allocate component arrays of a
!    lw_output_type variable.
!----------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1
!**************************************
      ! This is a temporary fix! Lw_output needs to be allocated at a higher level!
      ! Constructor and destructor for lw_output_type needs to be provided through
      ! rad_utilities
!**************************************
      call longwave_driver_alloc (ix, jx, kx, Lw_output(1))
      call longwave_driver_alloc (ix, jx, kx, Lw_output_std)
      if (Rad_control%do_lwaerosol_forcing) then
      ! This is a temporary fix! Lw_output needs to be allocated at a higher level!
        call longwave_driver_alloc (ix, jx, kx, Lw_output(Rad_control%indx_lwaf))
        call longwave_driver_alloc (ix, jx, kx, Lw_output_ad)
      endif

!--------------------------------------------------------------------
!    calculate the longwave radiative heating rates and fluxes.
!--------------------------------------------------------------------
      if (do_sealw99) then

             
!--------------------------------------------------------------------
!    call sealw99 to use the simplified-exchange-approximation (sea)
!    parameterization.
!----------------------------------------------------------------------
         if (Rad_control%do_lwaerosol_forcing) then
           if (Lw_control%do_lwaerosol) then
             calc_includes_aerosols = .false.
           else
             calc_includes_aerosols = .true.
           endif

!----------------------------------------------------------------------
!    call sealw99 with aerosols (if model is being run without) and 
!    without aerosols (if model is being run with). save the radiation
!    fluxes to Lw_output_ad (which does not feed back into the model),
!    but which may be used to define the aerosol forcing.
!----------------------------------------------------------------------
           call sealw99 (is, ie, js, je, Rad_time, Atmos_input,  &
                     Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
                     Cld_spec, Aerosol_diags, Lw_output_ad, &
                     Lw_diagnostics, calc_includes_aerosols)
           Lw_output(Rad_control%indx_lwaf) = Lw_output_ad
         endif
 
!----------------------------------------------------------------------
!    standard call, where radiation output feeds back into the model.
!----------------------------------------------------------------------
        call sealw99 (is, ie, js, je, Rad_time, Atmos_input,  &
                      Rad_gases, Aerosol, Aerosol_props, Cldrad_props, &
                      Cld_spec, Aerosol_diags, Lw_output_std,  &
                      Lw_diagnostics, Lw_control%do_lwaerosol)
        Lw_output(1) = Lw_output_std
      else

!--------------------------------------------------------------------
!    at the current time sealw99 is the only longwave parameterization 
!    available.
!----------------------------------------------------------------------
        call error_mesg ('longwave_driver_mod', &
         'invalid longwave radiation parameterization selected', FATAL)
      endif

      call longwave_driver_dealloc (Lw_output_std)
      if (Rad_control%do_lwaerosol_forcing) then
        call longwave_driver_dealloc (Lw_output_ad)
      endif

!---------------------------------------------------------------------

end subroutine longwave_driver


!#####################################################################
! <SUBROUTINE NAME="longwave_driver_end">
!  <OVERVIEW>
!   Subroutine to end longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine end longwave calculation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_driver_end                  

!--------------------------------------------------------------------
!    longwave_driver_end is the destructor for longwave_driver_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_driver_mod',   &
          'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    call sealw99_end to close sealw99_mod.
!-------------------------------------------------------------------
      if (do_sealw99) then
        call sealw99_end
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------


end subroutine longwave_driver_end                  


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="longwave_driver_alloc">
!  <OVERVIEW>
!   Subroutine to allocate output variables from longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_alloc (ix, jx, kx, Lw_output)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   Dimension 1 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   Dimension 2 length of radiation arrays to be allocated
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   Dimension 3 length of radiation arrays to be allocated
!  </IN>
!  <OUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_driver_alloc (ix, jx, kx, Lw_output)

!--------------------------------------------------------------------
!    longwave_driver_alloc allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!--------------------------------------------------------------------

integer,                   intent(in)    :: ix, jx, kx
type(lw_output_type),      intent(inout) :: Lw_output

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      ix,jx,kx     (i,j,k) dimensions of current physics window 
!
!
!   intent(inout) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!  
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    allocate and initialize arrays to hold net longwave fluxes and 
!    the longwave heating rate at each gridpoint. if the
!    cloud-forcing calculation is to be done, also allocate and init-
!    ialize arrays for fluxes and heating rates without clouds.
!-------------------------------------------------------------------
      allocate (Lw_output%flxnet( ix, jx, kx+1) )
      allocate (Lw_output%heatra( ix, jx, kx  ) )
      allocate (Lw_output%netlw_special   &
                                ( ix, jx, Rad_control%mx_spec_levs  ) )
      allocate (Lw_output%bdy_flx         &
                                ( ix, jx, 4) )
      Lw_output%flxnet(:,:,:) = 0.0
      Lw_output%heatra(:,:,:) = 0.0
      Lw_output%netlw_special(:,:,:) = 0.0
      Lw_output%bdy_flx (:,:,:) = 0.0      
      if (Rad_control%do_totcld_forcing)  then
        allocate (Lw_output%flxnetcf( ix, jx, kx+1) )
        allocate (Lw_output%heatracf( ix, jx, kx  ) )
        allocate (Lw_output%netlw_special_clr  &
                                ( ix, jx, Rad_control%mx_spec_levs  ) )
        allocate (Lw_output%bdy_flx_clr         &
                                ( ix, jx, 4) )
        Lw_output%flxnetcf(:,:,:) = 0.0
        Lw_output%heatracf(:,:,:) = 0.0
        Lw_output%netlw_special_clr(:,:,:) = 0.0
        Lw_output%bdy_flx_clr (:,:,:) = 0.0      
      endif
    
!--------------------------------------------------------------------

end subroutine longwave_driver_alloc


!#####################################################################
! <SUBROUTINE NAME="longwave_driver_dealloc">
!  <OVERVIEW>
!   Subroutine to deallocate output variables from longwave calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_driver_alloc (Lw_output)
!  </TEMPLATE>
!  <OUT NAME="Lw_output" TYPE="lw_output_type">
!   lw_output_type variable containing longwave 
!                   radiation output data
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_driver_dealloc (Lw_output)

!--------------------------------------------------------------------
!    longwave_driver_alloc deallocates the components
!    of the lw_output_type variable Lw_output.
!--------------------------------------------------------------------

type(lw_output_type),      intent(inout) :: Lw_output

!--------------------------------------------------------------------
!
!   intent(inout) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!  
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    deallocate arrays to hold net longwave fluxes and 
!    the longwave heating rate at each gridpoint
!-------------------------------------------------------------------
      deallocate (Lw_output%flxnet)
      deallocate (Lw_output%heatra)
      deallocate (Lw_output%netlw_special)
      deallocate (Lw_output%bdy_flx)
      if (Rad_control%do_totcld_forcing)  then
        deallocate (Lw_output%flxnetcf)
        deallocate (Lw_output%heatracf)
        deallocate (Lw_output%netlw_special_clr)
        deallocate (Lw_output%bdy_flx_clr)
      endif
    
!--------------------------------------------------------------------

end subroutine longwave_driver_dealloc

!###################################################################



                end module longwave_driver_mod
