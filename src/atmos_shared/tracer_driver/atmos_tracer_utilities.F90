module atmos_tracer_utilities_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Bruce Wyman
! </REVIEWER>


! <OVERVIEW>
!     This code provides some utility routines for atmospheric tracers in the FMS framework.
! </OVERVIEW>
! <DESCRIPTION>
!    This module gives utility routines which can be used to provide 
!    consistent removal mechanisms for atmospheric tracers. 
!
!    In particular it provides schemes for wet and dry deposiiton that 
!    can be easily utilized.
!
! </DESCRIPTION>


use            fms_mod, only : lowercase, &
                               write_version_number, &
                               stdlog, &
                               mpp_pe, &
                               mpp_root_pe, &
                               error_mesg, &
                               NOTE, FATAL
use   time_manager_mod, only : time_type
use   diag_manager_mod, only : send_data, &
                               register_diag_field
use tracer_manager_mod, only : query_method, &
                               get_tracer_names, &
                               get_number_tracers, &
                               MAX_TRACER_FIELDS
use  field_manager_mod, only : MODEL_ATMOS, parse
use   horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
                               horiz_interp_new, horiz_interp, horiz_interp_del
use  monin_obukhov_mod, only : mo_profile
use      constants_mod, only : GRAV, &     ! acceleration due to gravity [m/s2]
                               RDGAS, &    ! gas constant for dry air [J/kg/deg]
                               vonkarm, &
                               PI, &
                               DENS_H2O, & ! Water density [kg/m3]
                               WTMH2O, &   ! Water molecular weight [g/mole]
                               WTMAIR, &   ! Air molecular weight [g/mole]
                               AVOGNO      ! Avogadro's number
use   interpolator_mod, only : interpolator,  &
                               obtain_interpolator_time_slices, &
                               unset_interpolator_time_flag, &
                               interpolate_type
use      astronomy_mod, only : universal_time

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  wet_deposition,    &
        dry_deposition,    &
        dry_deposition_time_vary,    &
        dry_deposition_endts,        &
        interp_emiss,      &
        atmos_tracer_utilities_end, &
        atmos_tracer_utilities_init, &
        get_wetdep_param, &
        get_rh,   &
        get_w10m, &
        get_cldf, &
        sjl_fillz

!---- version number -----
character(len=128) :: version = '$Id: atmos_tracer_utilities.F90,v 20.0 2013/12/13 23:24:13 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

logical :: module_is_initialized = .FALSE.

character(len=7), parameter :: mod_name = 'tracers'
integer, parameter :: max_tracers = MAX_TRACER_FIELDS
!-----------------------------------------------------------------------
!--- identification numbers for  diagnostic fields and axes ----
integer :: id_tracer_ddep(max_tracers), id_tracer_dvel(max_tracers), &
           id_tracer_wdep_ls(max_tracers),   id_tracer_wdep_cv(max_tracers),  &
           id_tracer_wdep_lsin(max_tracers), id_tracer_wdep_cvin(max_tracers),&
           id_tracer_wdep_lsbc(max_tracers), id_tracer_wdep_cvbc(max_tracers),&
           id_tracer_reevap_ls(max_tracers), id_tracer_reevap_cv(max_tracers),&
           id_tracer_wdep_ls_3d(max_tracers)
integer :: id_tracer_ddep_cmip(max_tracers)
integer :: id_w10m, id_delm
integer :: id_u_star, id_b_star, id_rough_mom, id_z_pbl,  &
           id_mo_length_inv, id_vds
character(len=32),  dimension(max_tracers) :: tracer_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_longnames = ' '
character(len=32),  dimension(max_tracers) :: tracer_wdep_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_wdep_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_wdep_longnames = ' '
character(len=32),  dimension(max_tracers) :: tracer_ddep_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_dvel_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_ddep_units     = ' '
character(len=32),  dimension(max_tracers) :: tracer_dvel_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_ddep_longnames = ' '
character(len=128), dimension(max_tracers) :: tracer_dvel_longnames = ' '
real, allocatable :: blon_out(:,:), blat_out(:,:)
!----------------parameter values for the diagnostic units--------------
real, parameter :: mw_air = WTMAIR/1000.  ! Convert from [g/mole] to [kg/mole]
real, parameter :: mw_h2o = WTMH2O/1000.  ! Convert from [g/mole] to [kg/mole]
real, parameter :: twopi = 2*PI

type wetdep_type
   character (len=500) :: scheme, text_in_scheme, control
   real  :: Henry_constant
   real  :: Henry_variable
   real  :: frac_in_cloud
   real  :: frac_in_cloud_snow
   real  :: alpha_r
   real  :: alpha_s
   logical :: Lwetdep, Lgas, Laerosol, Lice
end type wetdep_type

type(wetdep_type), dimension(:), allocatable :: Wetdep


type drydep_type
   character (len=500) :: scheme, name, control
   real  :: land_dry_dep_vel
   real  :: sea_dry_dep_vel
   real  :: ice_dry_dep_vel
   real  :: snow_dry_dep_vel
   real  :: vegn_dry_dep_vel
   logical :: Ldrydep
end type drydep_type

type(drydep_type), dimension(:), allocatable :: Drydep


contains

!
! ######################################################################
!
!<SUBROUTINE NAME="atmos_tracer_utilities_init">
!<OVERVIEW>
! This is a routine to create and register the dry and wet deposition 
! fields of the tracers.
!</OVERVIEW>
!<DESCRIPTION>
!  This routine creates diagnostic names for dry and wet deposition fields of the tracers.
!  It takes the tracer name and appends "ddep" for the dry deposition field and "wdep" for 
!  the wet deposition field. This names can then be entered in the diag_table for 
!  diagnostic output of the tracer dry and wet deposition. The module name associated with
!  these fields in "tracers". The units of the deposition fields are assumed to be kg/m2/s.
!</DESCRIPTION>
!<TEMPLATE>
! call atmos_tracer_utilities_init(lonb,latb, mass_axes, Time)
!</TEMPLATE>
!   <IN NAME="lonb" TYPE="real" DIM="(:,:)">
!     The longitude corners for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real" DIM="(:,:)">
!     The latitude corners for the local domain.
!   </IN>
!   <IN NAME="mass_axes" TYPE="integer" DIM="(3)">
!     The axes relating to the tracer array.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>

subroutine atmos_tracer_utilities_init(lonb, latb, mass_axes, Time)

! Routine to initialize the tracer identification numbers. 
! This registers the 2D fields for the wet and dry deposition.
real, dimension(:,:),  intent(in) :: lonb, latb
integer, dimension(3), intent(in) :: mass_axes
type(time_type),       intent(in) :: Time

integer :: ntrace
character(len=20) :: units =''
!
integer :: n, logunit
character(len=128) :: name

logical  :: flag

! Make local copies of the local domain dimensions for use 
! in interp_emiss.
      allocate ( blon_out(size(lonb,1),size(lonb,2)))
      allocate ( blat_out(size(latb,1),size(latb,2)))
!      allocate ( data_out(size(lonb(:))-1, size(latb(:))-1))
      blon_out = lonb
      blat_out = latb
      
      do n = 1, max_tracers
         write ( tracer_names(n),     100 ) n
         write ( tracer_longnames(n), 102 ) n
         tracer_units(n) = 'none'
      enddo
  100 format ('tr',i3.3)
  102 format ('tracer ',i3.3)

call get_number_tracers(MODEL_ATMOS, num_tracers= ntrace)

   if (ntrace > 0) then
     allocate (Wetdep(ntrace))
     allocate (Drydep(ntrace))
   endif
   do n = 1, ntrace
!--- set tracer tendency names where tracer names have changed ---

call get_tracer_names(MODEL_ATMOS,n,tracer_names(n),tracer_longnames(n),tracer_units(n))
      write (name,100) n
      if (trim(tracer_names(n)) /= name) then
          tracer_ddep_names(n) = trim(tracer_names(n)) //'_ddep'
          tracer_dvel_names(n) = trim(tracer_names(n)) //'_dvel'
          tracer_wdep_names(n) = trim(tracer_names(n)) //'_wdep'
      endif
      write (name,102) n
      if (trim(tracer_longnames(n)) /= name) then
          tracer_wdep_longnames(n) = &
                  trim(tracer_longnames(n)) // ' wet deposition for tracers'
          tracer_ddep_longnames(n) = &
                  trim(tracer_longnames(n)) // ' dry deposition for tracers'
          tracer_dvel_longnames(n) = &
                  trim(tracer_longnames(n)) // ' dry deposition velocity for tracers'
      endif

      select case (trim(tracer_units(n)))
        case ('mmr')
          units = 'kg/m2/s'
        case ('kg/kg')
          units = 'kg/m2/s'
        case ('vmr')
          units = 'mole/m2/s'
        case ('mol/mol')
          units = 'mole/m2/s'
        case ('mole/mole')
          units = 'mole/m2/s'
        case default
          units = trim(tracer_units(n))//' kg/(m2 s)'
          call error_mesg('atmos_tracer_utilities_init',&
          ' Dry dep units set to '//trim(units)//' in atmos_tracer_utilities for '//trim(tracer_names(n)),&
           NOTE)
      end select

    
      flag = query_method ('wet_deposition',MODEL_ATMOS,n, &
                            Wetdep(n)%text_in_scheme,Wetdep(n)%control)
      call get_wetdep_param(Wetdep(n)%text_in_scheme,  &
                            Wetdep(n)%control,&
                            Wetdep(n)%scheme, &
                            Wetdep(n)%Henry_constant,  &
                            Wetdep(n)%Henry_variable, &
                            Wetdep(n)%frac_in_cloud, &
                            Wetdep(n)%frac_in_cloud_snow, &
                            Wetdep(n)%alpha_r, Wetdep(n)%alpha_s, &
                            Wetdep(n)%Lwetdep, Wetdep(n)%Lgas, &
                            Wetdep(n)%Laerosol, Wetdep(n)%Lice )

      Drydep(n)%Ldrydep = query_method ('dry_deposition', MODEL_ATMOS,&
                                  n,Drydep(n)%name, Drydep(n)%control)


      call get_drydep_param(Drydep(n)%name,Drydep(n)%control,  &
                         Drydep(n)%scheme,Drydep(n)%land_dry_dep_vel,  &
                                             Drydep(n)%sea_dry_dep_vel)
! When formulation of dry deposition is resolved perhaps use the following?
!      call get_drydep_param    &
!           (Drydep(n)%name, Drydep(n)%control, Drydep(n)%scheme,  &
!            Drydep(n)%land_dry_dep_vel, Drydep(n)%sea_dry_dep_vel,  & 
!            Drydep(n)%ice_dry_dep_vel, Drydep(n)%snow_dry_dep_vel, &
!                                             Drydep(n)%vegn_dry_dep_vel )

! Register the dry deposition of the n tracers
     id_tracer_ddep(n) = register_diag_field ( mod_name,                    &
            trim(tracer_ddep_names(n)), mass_axes(1:2), Time,               &
            trim(tracer_ddep_longnames(n)),                                 &
            trim(units), missing_value=-999.     )
     id_tracer_ddep_cmip(n) = register_diag_field ( mod_name,               &
            trim(tracer_ddep_names(n))//'_cmip', mass_axes(1:2), Time,      &
            trim(tracer_ddep_longnames(n)),                                 &
           'kg/m2/s', missing_value=-999.     )
! Register the dry deposition of the n tracers
     id_tracer_dvel(n) = register_diag_field ( mod_name,                    &
            trim(tracer_dvel_names(n)), mass_axes(1:2), Time,               &
            trim(tracer_dvel_longnames(n)),                                 &
            'm/s', missing_value=-999.     )
! Register the wet deposition of the n tracers by large scale clouds
     id_tracer_wdep_ls(n) = register_diag_field ( mod_name,                 &
            trim(tracer_wdep_names(n))//'_ls', mass_axes(1:2), Time,        &
            trim(tracer_wdep_longnames(n))//' in large scale',              &
            trim(units), missing_value=-999.    )

! Register the wet deposition of the n tracers by large scale clouds 3d
     id_tracer_wdep_ls_3d(n) = register_diag_field ( mod_name,         &
            trim(tracer_wdep_names(n))//'_ls_3d', mass_axes(1:3), Time, &
            trim(tracer_wdep_longnames(n))//' in large scale 3D',       &
            trim(units), missing_value=-999.    )

! Register the wet deposition of the n tracers by convective clouds
     id_tracer_wdep_cv(n) = register_diag_field ( mod_name,               &
              trim(tracer_wdep_names(n))//'_cv', mass_axes(1:2), Time, &
              trim(tracer_wdep_longnames(n))//' in convective scheme',                   &
              trim(units), missing_value=-999.    )
! Register in-cloud rainout by large scale clouds
     id_tracer_wdep_lsin(n) = register_diag_field ( mod_name,               &
            trim(tracer_wdep_names(n))//'_lsin', mass_axes(1:2), Time,      &
            trim(tracer_wdep_longnames(n))//' in_cloud by lscale precip',   &
            trim(units), missing_value=-999.    )
! Register below-cloud washout by large scale clouds
     id_tracer_wdep_lsbc(n) = register_diag_field ( mod_name,               &
            trim(tracer_wdep_names(n))//'_lsbc', mass_axes(1:2), Time,      &
            trim(tracer_wdep_longnames(n))//' below_cloud by lscale precip',&
            trim(units), missing_value=-999.  )
! Register in-cloud re-evaporation by large scale clouds
     id_tracer_reevap_ls(n) = register_diag_field ( mod_name,               &
            trim(tracer_names(n))//'_reevap_ls', mass_axes(1:3), Time,      &
            trim(tracer_longnames(n))//' re-evap by lscale clouds',         &
            trim(units), missing_value=-999.    )
! Register in-cloud rainout by convective clouds
! Register in-cloud rainout of the n tracers by convective clouds
     id_tracer_wdep_cvin(n) = register_diag_field ( mod_name,               &
            trim(tracer_wdep_names(n))//'_cvin', mass_axes(1:2), Time,      &
            trim(tracer_wdep_longnames(n))//' in_cloud by conv precip',     &
            trim(units), missing_value=-999.    )
! Register below-cloud washout by convective clouds
     id_tracer_wdep_cvbc(n) = register_diag_field ( mod_name,               &
            trim(tracer_wdep_names(n))//'_cvbc', mass_axes(1:2), Time,      &
            trim(tracer_wdep_longnames(n))//' below_cloud by conv precip',  &
            trim(units), missing_value=-999.  )
! Register re-evaporation by convective clouds
     id_tracer_reevap_cv(n) = register_diag_field ( mod_name,               &
            trim(tracer_names(n))//'_reevap_cv', mass_axes(1:3), Time,      &
            trim(tracer_longnames(n))//' re-evap by conv precip',         &
            trim(units), missing_value=-999.    )
   enddo
! Register scaling factor to calculate wind speed at 10 meters
   id_delm   = register_diag_field ( mod_name,                &
               'delm', mass_axes(1:2),Time,                   &
               'Scaling factor', 'none',                      &
                missing_value=-999.                           )
! Register the wind speed at 10 meters
   id_w10m   = register_diag_field ( mod_name,                &
               'w10m', mass_axes(1:2),Time,                   &
               'Wind speed at 10 meters', 'm/s',              &
                missing_value=-999.                           )

     id_u_star = register_diag_field ( mod_name,                    &
            'u_star_atm', mass_axes(1:2), Time,               &
            'u star',                                 &
            'm/s', missing_value=-999.     )
     id_b_star = register_diag_field ( mod_name,                    &
            'b_star_atm', mass_axes(1:2), Time,               &
            'b star',                                 &
            'm/s2', missing_value=-999.     )
     id_rough_mom = register_diag_field ( mod_name,                    &
            'rough_mom_atm', mass_axes(1:2), Time,               &
            'rough length z0',                                 &
            'm', missing_value=-999.     )
     id_z_pbl = register_diag_field ( mod_name,                    &
            'z_pbl_atm', mass_axes(1:2), Time,               &
            'z pbl',                                 &
            'm', missing_value=-999.     )
     id_mo_length_inv = register_diag_field ( mod_name,                &
            'mo_length_inv_atm', mass_axes(1:2), Time,               &
            'monin Obukhov length',                                 &
            'm', missing_value=-999.     )
     id_vds = register_diag_field ( mod_name,                    &
            'vds_atm', mass_axes(1:2), Time,               &
            'vds',                                 &
            'm/s', missing_value=-999.     )

 
      call write_version_number (version, tagname)

    if ( mpp_pe() == mpp_root_pe() ) then
         logunit=stdlog()
         call write_namelist_values (logunit,ntrace)
    endif

      module_is_initialized = .TRUE.

end subroutine atmos_tracer_utilities_init


!####################################################################

subroutine dry_deposition_time_vary (drydep_data, Time)

type(time_type), intent(in) :: Time
type(interpolate_type), dimension(:), intent(inout) :: drydep_data

      integer :: n


      do n=1,size(drydep_data,1)
        if (Drydep(n)%Ldrydep .and. Drydep(n)%scheme == 'file') then
          call obtain_interpolator_time_slices (drydep_data(n), Time)
        endif
      end do

end subroutine dry_deposition_time_vary 



!####################################################################

subroutine dry_deposition_endts (drydep_data)                          

type(interpolate_type), dimension(:), intent(inout) :: drydep_data

      integer :: n

      do n=1, size(drydep_data,1)     
        call unset_interpolator_time_flag (drydep_data(n))
      end do


end subroutine dry_deposition_endts       



!####################################################################

!</SUBROUTINE>
!
!#######################################################################
!
subroutine write_namelist_values (unit, ntrace)
    integer, intent(in) :: unit, ntrace
    integer :: n

    write (unit,10)
    do n = 1, ntrace
       write (unit,11) trim(tracer_wdep_names(n)),     &
                       trim(tracer_wdep_longnames(n)), &
                       trim(tracer_wdep_units(n))
       write (unit,11) trim(tracer_ddep_names(n)),     &
                       trim(tracer_ddep_longnames(n)), &
                       trim(tracer_ddep_units(n))
       write (unit,11) trim(tracer_dvel_names(n)),     &
                       trim(tracer_dvel_longnames(n)), &
                       'm/s'
    enddo

 10 format (' &TRACER_DIAGNOSTICS_NML', &
          /,'    TRACER:  names  longnames  (units)')
 11 format (a16,2x,a,2x,'(',a,')')

 end subroutine write_namelist_values

!
!#######################################################################
!
!<SUBROUTINE NAME = "dry_deposition">
subroutine dry_deposition( n, is, js, u, v, T, pwt, pfull, dz, &
                           u_star, landmask, dsinku, tracer, Time, &
                           Time_next, lon, half_day, drydep_data)
! When formulation of dry deposition is resolved perhaps use the following?
!                           landfr, seaice_cn, snow_area, & 
!                           vegn_cover, vegn_lai, & 
!                           b_star, z_pbl, rough_mom )
!
!<OVERVIEW>
! Routine to calculate the fraction of tracer to be removed by dry 
! deposition.
!</OVERVIEW>
!<DESCRIPTION>
! There are three types of dry deposition coded.
!
! 1) Wind driven derived dry deposition velocity.
!
! 2) Fixed dry deposition velocity.
! 
! 3) Dry deposition velocities read in from input file

! There are an addition three types of dry deposition coded 
! but presently commented out.
!
! 4) Wind driven derived dry deposition velocity, surface dependent.
!
! 5) Wind driven derived dry deposition velocity, surface and boundary
!    layer stability dependent.
!
! 6) Fixed dry deposition velocity, difference between land, snow-covered
!    land, sea and ice-covered sea.
! 
! The theory behind the wind driven dry deposition velocity calculation
! assumes that the deposition can be modeled as a parallel resistance type 
! problem.
!
!  Total resistance to HNO3-type dry deposition, 
!<PRE>       R = Ra + Rb
!  resisa = aerodynamic resistance
!  resisb = surface resistance (laminar layer + uptake)
!         = 5/u*  [s/cm]        for neutral stability
!      Vd = 1/R
!</PRE>
! For the fixed dry deposition velocity, there is no change in the 
! deposition velocity but the variation of the depth of the surface 
! layer implies that there is variation in the amount deposited.
!
! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "dry_deposition","wind_driven","surfr=XXX"
!     where XXX is the total resistance defined above.
!
! "dry_deposition","fixed","land=XXX, sea=YYY"
!     where XXX is the dry deposition velocity (m/s) over land
!       and YYY is the dry deposition velocity (m/s) over sea.
!
! "dry_deposition","file","FILENAME.NC"
!     where FILENAME.NC is the NetCDF file name.
!</PRE>
!</DESCRIPTION>
!<TEMPLATE>
! call dry_deposition( n, is, js, u, v, T, pwt, pfull, dz,
!                      u_star, landmask, dsinku, tracer, Time, drydep_data)
!</TEMPLATE>
!
!  <IN NAME="n" TYPE="integer">
!    The tracer number.
!  </IN>
!  <IN NAME="is, js" TYPE="integer">
!    Start indices for array (computational indices).
!  </IN>
!  <IN NAME="u" TYPE="real" DIM="(:,:)">
!    U wind field.
!  </IN>
!  <IN NAME="v" TYPE="real" DIM="(:,:)">
!    V wind field.
!  </IN>
!  <IN NAME="T" TYPE="real" DIM="(:,:)">
!    Temperature.
!  </IN>
!  <IN NAME="pwt" TYPE="real" DIM="(:,:)">
!     Pressure differential of half levels.
!  </IN>
!  <IN NAME="pfull" TYPE="real" DIM="(:,:)">
!     Full pressure levels.
!  </IN>
!  <IN NAME="u_star" TYPE="real" DIM="(:,:)">
!     Friction velocity.
!  </IN>
!  <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude.
!  </IN>
!  <IN NAME="landmask" TYPE="logical">
!     Land - sea mask.
!  </IN>
!  <INOUT NAME="drydep_data" TYPE="interpolate_type">
!     Dry deposition data interpolated from input file.
!  </INOUT>
!
!  <OUT NAME="dsinku" TYPE="real" DIM="(:,:)">
!    The amount of tracer in the surface layer which is dry deposited per second.
!  </OUT>
!
integer, intent(in)                 :: n, is, js
real, intent(in), dimension(:,:)    :: u, v, T, pwt, pfull, u_star, tracer, dz
real, intent(in), dimension(:,:)    :: lon, half_day
logical, intent(in), dimension(:,:) :: landmask
! When formulation of dry deposition is resolved perhaps use the following?
!real, intent(in), dimension(:,:)    :: landfr, z_pbl, b_star, rough_mom
!real, intent(in), dimension(:,:)    :: seaice_cn, snow_area, vegn_cover,  &
!                                       vegn_lai
type(time_type), intent(in)         :: Time, Time_next
type(interpolate_type),intent(inout)  :: drydep_data
real, intent(out), dimension(:,:)   :: dsinku

real,dimension(size(u,1),size(u,2))   :: hwindv,frictv,resisa,drydep_vel
!real,dimension(size(u,1),size(u,2))   :: mo_length_inv, vds, rs, k1, k2
integer :: i,j, flagsr, id, jd
real    :: land_dry_dep_vel, sea_dry_dep_vel, ice_dry_dep_vel,  &
           snow_dry_dep_vel, vegn_dry_dep_vel,   &
           surfr, sear, icer,  snowr, vegnr
real    :: diag_scale
real    :: factor_tmp, gmt, dv_on, dv_off, dayfrac, vd_night, vd_day, loc_angle
logical :: used, diurnal
integer :: flag_species, flag_diurnal
character(len=10) ::units,names
character(len=500) :: name,control,scheme, speciesname,dummy

! Default zero
dsinku = 0.0
if (.not. Drydep(n)%Ldrydep) return
name =Drydep(n)%name
control = Drydep(n)%control
scheme = Drydep(n)%scheme
land_dry_dep_vel = Drydep(n)%land_dry_dep_vel
sea_dry_dep_vel = Drydep(n)%sea_dry_dep_vel
!ice_dry_dep_vel = Drydep(n)%ice_dry_dep_vel
!snow_dry_dep_vel = Drydep(n)%snow_dry_dep_vel
!vegn_dry_dep_vel = Drydep(n)%vegn_dry_dep_vel

! delta z = dp/(rho * grav)
! delta z = RT/g*dp/p    pwt = dp/g
!dz(:,:) = pwt(:,:)*rdgas*T(:,:)/pfull(:,:)
id=size(pfull,1); jd=size(pfull,2)


  select case(lowercase(scheme))
  
   case('wind_driven')
! Calculate horizontal wind velocity and aerodynamic resistance:
!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2), 
!   and  u*=sqrt(momentum flux)  is friction velocity.
!
!****  Compute dry sinks (loss frequency, need modification when 
!****    different vdep values are to be used for species)
        flagsr=parse(control,'surfr',surfr)
        if(flagsr == 0) surfr=500.
        hwindv=sqrt(u**2+v**2)
        frictv=u_star
        resisa=hwindv/(u_star*u_star)
        where (frictv .lt. 0.1) frictv=0.1
        dsinku = (1./(surfr/frictv + resisa))/dz
        drydep_vel(:,:) = 0.

!    case('sfc_dependent_wind_driven')
!! Calculate horizontal wind velocity and aerodynamic resistance:
!!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2), 
!!   and  u*=sqrt(momentum flux)  is friction velocity.
!!
!!****  Compute dry sinks (loss frequency, need modification when 
!!****    different vdep values are to be used for species)
!        flagsr=parse(control,'surfr',surfr)
!        if(flagsr == 0) surfr=500.
!
!        flagsr=parse(control,'sear',sear)
!        if(flagsr == 0) sear=surfr
!        
!        flagsr=parse(control,'icer',icer)
!        if(flagsr == 0) icer=surfr
!        
!        flagsr=parse(control,'snowr',snowr)
!        if(flagsr == 0) snowr=surfr
!        
!        flagsr=parse(control,'vegnr',vegnr)
!        if(flagsr == 0) vegnr=surfr
!
!        hwindv=sqrt(u**2+v**2)
!        frictv=u_star
!        resisa=hwindv/(u_star*u_star)
!        where (frictv .lt. 0.1) frictv=0.1
!        drydep_vel(:,:) = (1./(surfr/frictv + resisa))*  &
!                                         (landfr(:,:) - snow_area(:,:)) + &
!                          (1./(snowr/frictv + resisa))*snow_area(:,:)  +&
!                          (1./(sear/frictv + resisa))*   &
!                                 (1. - landfr(:,:) - seaice_cn(:,:))   + &
!                          (1./(icer/frictv + resisa))*seaice_cn(:,:)
!        dsinku(:,:) = drydep_vel(:,:) / dz(:,:)
!
!    case('sfc_BL_dependent_wind_driven')
!! Calculate horizontal wind velocity and aerodynamic resistance:
!!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2), 
!!   and  u*=sqrt(momentum flux)  is friction velocity.
!!
!!****  Compute dry sinks (loss frequency, need modification when 
!!****    different vdep values are to be used for species)
!        flagsr=parse(control,'surfr',surfr)
!        if(flagsr == 0) surfr=500.
!
!        flagsr=parse(control,'sear',sear)
!        if(flagsr == 0) sear=surfr
!        
!        flagsr=parse(control,'icer',icer)
!        if(flagsr == 0) icer=surfr
!        
!        flagsr=parse(control,'snowr',snowr)
!        if(flagsr == 0) snowr=surfr
!        
!        flagsr=parse(control,'vegnr',vegnr)
!        if(flagsr == 0) vegnr=surfr
!
!        hwindv=sqrt(u**2+v**2)
!        frictv=u_star
!        resisa=hwindv/(u_star*u_star)
!        where (frictv .lt. 0.1) frictv=0.1
!        mo_length_inv = - vonkarm * b_star/(frictv*frictv)      
!        where (rough_mom > 0.005)
!           k1 =  0.001222 * log10(rough_mom) + 0.003906
!        elsewhere
!           k1 =  0.001222 * log10(0.005) + 0.003906
!        endwhere        
!        where(mo_length_inv < 0)
!           k2 = 0.0009 * ( - z_pbl * mo_length_inv)**(2.0/3.0)
!        elsewhere
!           k2 = 0.0
!        endwhere        
!        vds = frictv * (k1 + k2)
!        rs  = 1.0 / vds
!        drydep_vel(:,:) = (1./(rs + resisa))*      &
!                                   (landfr(:,:) - snow_area(:,:))  +  &
!                          (1./(snowr/frictv + resisa))*snow_area(:,:)  + &
!                          (1./(sear/frictv + resisa))*   &
!                                (1. - landfr(:,:) - seaice_cn(:,:))  + &
!                          (1./(icer/frictv + resisa))*seaice_cn(:,:)
!        dsinku(:,:) = drydep_vel(:,:) / dz(:,:)

    case('fixed')
! For the moment let's try to calculate the delta-z of the bottom 
! layer and using a simple dry deposition velocity times the 
! timestep, idt, calculate the fraction of the lowest layer which 
! deposits.
       where (landmask(:,:))
! dry dep value over the land surface
         drydep_vel(:,:) = land_dry_dep_vel
      elsewhere
! dry dep value over the sea surface
         drydep_vel(:,:) = sea_dry_dep_vel
      endwhere
      dsinku(:,:) = drydep_vel(:,:) / dz(:,:)

!    case('sfc_dependent_fixed')
!      drydep_vel(:,:) = land_dry_dep_vel*    &
!                                     (landfr(:,:) - snow_area(:,:))   + &
!                        snow_dry_dep_vel*snow_area(:,:)   +  &
!                        sea_dry_dep_vel*  &
!                                  (1. - landfr(:,:) - seaice_cn(:,:))  + &
!                        ice_dry_dep_vel*seaice_cn(:,:)
!      dsinku(:,:) = drydep_vel(:,:) / dz(:,:) 

    case('file')
        flag_species = parse(control,'name',speciesname)
        if(flag_species>0) then
           name = trim(speciesname)
        else
           call get_tracer_names(MODEL_ATMOS,n,name)
        endif
        flag_diurnal = parse(control,'diurnal',dummy)
        diurnal = (flag_diurnal > 0)
        call interpolator( drydep_data, Time, drydep_vel, trim(name), is, js )

        if (diurnal) then
           do j = 1,jd
           do i = 1,id
! half_day is between 0 and pi, so dv_off btwn 0 to pi, dv_on btwn -pi and 0
              dv_off = MIN( 1.2*half_day(i,j), PI )
              dv_on = -dv_off
              dayfrac = dv_off/PI
! apply the mean dep vel during polar day or polar night (or nearby)
              if (dv_off > 0 .and. dv_off < PI  ) then
                 vd_night = MIN(0.001, 0.5*drydep_vel(i,j))
                 vd_day = ( drydep_vel(i,j)-vd_night*(1.-dayfrac) ) / dayfrac
                 gmt = universal_time(Time)
                 loc_angle = gmt + lon(i,j) - PI
                 if (loc_angle >= PI) loc_angle = loc_angle - twopi
                 if (loc_angle < -PI) loc_angle = loc_angle + twopi
                 if( loc_angle >= dv_off .or. loc_angle <= dv_on ) then
                    drydep_vel(i,j) = vd_night
                 else
                    factor_tmp = loc_angle - dv_on
                    factor_tmp = factor_tmp / MAX(2*dv_off,1.e-6)
                    drydep_vel(i,j) = 0.5*PI*sin(factor_tmp*PI)*(vd_day-vd_night) + vd_night
                 end if
              end if
           end do
           end do
        end if !(diurnal)

        dsinku(:,:) = drydep_vel(:,:) / dz(:,:)
    case('default')
        drydep_vel(:,:) = 0.
  end select

dsinku(:,:) = MAX(dsinku(:,:), 0.0E+00)
where(tracer>0)
  dsinku=dsinku*tracer
elsewhere
  dsinku=0.0
endwhere

! Now save the dry deposition to the diagnostic manager
! delta z = dp/(rho * grav)
! delta z *rho  = dp/g
! tracer(kgtracer/kgair) * dz(m)* rho(kgair/m3) = kgtracer/m2
! so rho drops out of the equation
    if (id_tracer_ddep(n) > 0 ) then
      call get_tracer_names(MODEL_ATMOS,n,names,units=units)
      select case (trim(units))
        case ('vmr')
          diag_scale = mw_air
        case ('mol/mol')
          diag_scale = mw_air
        case ('mole/mole')
          diag_scale = mw_air
        case default
          diag_scale = 1.
      end select
      used = send_data ( id_tracer_ddep(n), dsinku*pwt/diag_scale, Time_next, &
          is_in =is,js_in=js)
    endif
    if (id_tracer_ddep_cmip(n) > 0 ) then
      call get_tracer_names(MODEL_ATMOS,n,names,units=units)
      select case (trim(names))
        case ('so2')
          diag_scale = mw_air/0.064
        case ('so4')
          diag_scale = mw_air/0.096
        case ('dms')
          diag_scale = mw_air/0.062
        case ('nh3')
          diag_scale = mw_air/0.017
        case default
          diag_scale = 1.
        end select
       used = send_data ( id_tracer_ddep_cmip(n), dsinku*pwt/diag_scale,Time_next, &
           is_in =is,js_in=js)
    endif
    if (id_tracer_dvel(n) > 0 ) then
      used = send_data ( id_tracer_dvel(n), drydep_vel, Time_next, &
          is_in =is,js_in=js)
    end if
end subroutine dry_deposition
!</SUBROUTINE>
!
!#######################################################################
!
!<SUBROUTINE NAME = "wet_deposition">
!<TEMPLATE>
!CALL wet_deposition( n, T, pfull, phalf, zfull, zhalf, &
!                     rain, snow, qdt, cloud, rain3d, snow3d, &
!                     tracer, tracer_dt, Time, cloud_param, is, js, dt )
!</TEMPLATE>
subroutine wet_deposition( n, T, pfull, phalf, zfull, zhalf, &
                           rain, snow, qdt, cloud, cloud_frac, &
                           f_snow_berg, rain3d, snow3d, &
                           tracer, tracer_dt, Time, cloud_param, &
                           is, js, dt, sum_wdep_out )
!      
!<OVERVIEW>
! Routine to calculate the fraction of tracer removed by wet deposition
!</OVERVIEW>
!
!<IN NAME="n" TYPE="integer">
!   Tracer number
!</IN>
!<IN NAME="is, js" TYPE="integer">
!   start indices for array (computational indices)
!</IN>
!<IN NAME="T" TYPE="real" DIM="(:,:,:)">
!   Temperature
!</IN>
!<IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!   Full level pressure field (Pa)
!</IN>
!<IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!   Half level pressure field (Pa)
!</IN>
!<IN NAME="zfull" TYPE="real" DIM="(:,:,:)">
!   Full level height field (m)
!</IN>
!<IN NAME="zhalf" TYPE="real" DIM="(:,:,:)">
!   Half level height field (m)
!</IN>
!<IN NAME="rain" TYPE="real" DIM="(:,:)">
!   Precipitation in the form of rain
!</IN>
!<IN NAME="snow" TYPE="real" DIM="(:,:)">
!   Precipitation in the form of snow
!</IN>
!<IN NAME="qdt" TYPE="real" DIM="(:,:,:)">
!   The tendency of the specific humidity (+ condenstate) due to the cloud parametrization (kg/kg/s)
!</IN>
!<IN NAME="cloud" TYPE="real" DIM="(:,:,:)">
!   Cloud amount (liquid + ice) (kg/kg)
!</IN>
!<IN NAME="cloud_frac" TYPE="real" DIM="(:,:,:)">
!   Cloud area fraction
!</IN>
!<IN NAME="rain3d" TYPE="real" DIM="(:,:,:)">
!   Precipitation in the form of rain (kg/m2/s)
!</IN>
!<IN NAME="snow3d" TYPE="real" DIM="(:,:,:)">
!   Precipitation in the form of snow (kg/m2/s)
!</IN>
!<IN NAME="tracer" TYPE="real" DIM="(:,:,:)">
!   The tracer field 
!</IN>
!<IN NAME="Time" TYPE="type(time_type)">
!   The time structure for submitting wet deposition as a diagnostic
!</IN>
!<IN NAME="cloud_param" TYPE="character">
!   Is this a convective (convect) or large scale (lscale) cloud parametrization?
!</IN>
!<IN NAME="dt" TYPE="real">
!   The model timestep (in seconds)
!</IN>
!<OUT NAME="tracer_dt" TYPE="real" DIM="(:,:,:)">
!   The tendency of the tracer field due to wet deposition
!</OUT>
!<DESCRIPTION>
! Schemes allowed here are 
!
! 1) Deposition removed in the same fractional amount as the modeled precipitation rate is to 
!    a standardized precipitation rate.
!    Basically this scheme assumes that a fractional area of the gridbox is affected by 
!    precipitation and that this precipitation rate is due to a cloud of standardized cloud 
!    liquid water content. Removal is constant throughout the column where precipitation is occuring.
!
! 2) Removal according to Henry's Law. This law states that the ratio of the concentation in 
!    cloud water and the partial pressure in the interstitial air is a constant. If tracer
!    is in VMR, the units for Henry's constant are mole/L/Pa (normally it is mole/L/atm).
!    Parameters for a large number of species can be found at
!    http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!
! 3) Aerosol removal, using specified in-cloud tracer fraction

! 4) Similar as 3) with some lwh modifications
!
! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "wet_deposition","henry","henry=XXX, dependence=YYY"
!     where XXX is the Henry's constant for the tracer in question
!       and YYY is the temperature dependence of the Henry's Law constant.
!
! "wet_deposition","fraction","lslwc=XXX, convlwc=YYY"
!     where XXX is the liquid water content of a standard large scale cloud
!       and YYY is the liquid water content of a standard convective cloud.
!</PRE>

!</DESCRIPTION>

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
integer,          intent(in)                     :: n, is, js
real,             intent(in),  dimension(:,:,:)  :: T, pfull,phalf, zfull, zhalf, qdt, cloud, tracer
real,             intent(in),  dimension(:,:,:)  :: cloud_frac
real,             intent(in),  dimension(:,:,:)  :: f_snow_berg     
                                     ! snow production by Bergeron process
real,             intent(in),  dimension(:,:)    :: rain, snow
character(len=*), intent(in)                     :: cloud_param
type (time_type), intent(in)                     :: Time
real,             intent(out), dimension(:,:,:)  :: tracer_dt
real,             intent(in)                     :: dt
real,             intent(in),  dimension(:,:,:)  :: rain3d, snow3d
real,             intent(out),  dimension(:,:), optional :: sum_wdep_out

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
real, dimension(size(T,1),size(T,2),size(pfull,3)) :: &
      Htemp, xliq, n_air, rho_air, pwt, zdel, precip3d, scav_fact3d, &
      precip3ds, precip3dr
real, dimension(size(T,1),size(T,2)) :: &
      temp_factor, scav_factor, washout, sum_wdep, &
      w_h2o, K1, K2, beta, f_a,  scav_factor_s,  &
      wdep_in, wdep_bc, fluxr,fluxs, tracer_flux
real, dimension(size(T,1),size(T,2),size(pfull,3)) :: &
      in_temp, bc_temp, dt_temp, reevap_fraction, reevap_diag
integer, dimension(size(T,1),size(T,2)) :: &
      ktopcd, kendcd
real, dimension(size(rain3d,1),size(rain3d,2),size(rain3d,3)) :: rainsnow3d

integer :: i, j, k, kk, id, jd, kd, flaglw

real, dimension(size(T,3)) :: conc

real :: conc_rain, conc_rain_total, conc_sat

real, parameter ::  DENS_SNOW = 500.    ! Snow density [kg/m3]
real, parameter ::  RGAS      = 8.3143  ! ideal gas constant Pa m3/mol/K

real    :: &
      Henry_constant, Henry_variable, &
      clwc, wash, premin, prenow, hwtop, &
      diag_scale

real, parameter :: &
      inv298p15 = 1./298.15, &     ! 1/K
      kboltz = 1.38E-23,         & ! J/K
      rain_diam  = 1.89e-3,     &  ! mean diameter of rain drop (m)
      rain_vterm = 7.48,        &  ! rain drop terminal velocity (m/s)
      vk_air = 6.18e-6,         &  ! kinematic viscosity of air (m^2/s)
      d_g = 1.12e-5,            &  ! diffusive coefficient (m^2/s)
      geo_fac = 6.,             &  ! geometry factor (surface area/volume = geo_fac/diameter)
      cm3_2_m3 = 1.e-6             ! m3/cm3

real :: &
      k_g,                       & ! mass transfer coefficient (m/s)
      stay,                      & ! fraction
      fall_time                    ! fall time through layer (s)

real :: f_a0, scav_factor0, sa_drop0, fgas0
real :: frac_in_cloud, frac_in_cloud_snow, frac_int, ph
real , parameter :: &
      R_r = 0.001, &               ! radius of cloud-droplets for rain
      R_s = 0.001, &               ! radius of cloud-droplets for snow
      frac_int_gas = 1.0,   &
      frac_int_aerosol= 0.5

real :: alpha_r, alpha_s

logical :: &
      used, &
      Lwetdep, Lgas, Laerosol, Lice
character(len=500) :: &
      tracer_name, control, scheme, units, &
      text_in_scheme

!-----------------------------------------------------------------------

    ktopcd = 0
    kendcd = 0

    tracer_dt   = 0.
    wdep_in     = 0.
    wdep_bc     = 0.
    beta        = 0.
    reevap_fraction = 0.
    reevap_diag = 0.
    tracer_flux = 0.

    sum_wdep = 0.

    id = size(T,1)
    jd = size(T,2)
    kd = size(T,3)

    call get_tracer_names(MODEL_ATMOS,n,tracer_name, units = units)
    if ( .not. Wetdep(n)%Lwetdep) return
     text_in_scheme = Wetdep(n)%text_in_scheme
     control = Wetdep(n)%control
     scheme = Wetdep(n)%scheme
     Henry_constant = Wetdep(n)%Henry_constant
     Henry_variable = Wetdep(n)%Henry_variable
     frac_in_cloud = Wetdep(n)%frac_in_cloud
     frac_in_cloud_snow = Wetdep(n)%frac_in_cloud_snow
     alpha_r = Wetdep(n)%alpha_r
     alpha_s = Wetdep(n)%alpha_s
     Lwetdep = Wetdep(n)%Lwetdep
     Lgas    = Wetdep(n)%Lgas
     Laerosol = Wetdep(n)%Laerosol
     Lice    = Wetdep(n)%Lice

    rho_air(:,:,:) = pfull(:,:,:) / ( T(:,:,:)*RDGAS ) ! kg/m3
!   Lice = .not. (scheme=='henry_noice' .or. scheme=='henry_below_noice' .or. &
!                 scheme=='aerosol_noice' .or. scheme=='aerosol_below_noice' )
    if (Lice) then
       rainsnow3d(:,:,:) = rain3d(:,:,:) + snow3d(:,:,:)
    else
       rainsnow3d(:,:,:) = rain3d(:,:,:)
    end if
    do k=1,kd
       precip3d(:,:,k) = rainsnow3d(:,:,k+1)-rainsnow3d(:,:,k)
       precip3dr(:,:,k) = rain3d(:,:,k+1)-rain3d(:,:,k)
       pwt(:,:,k)  = ( phalf(:,:,k+1) - phalf(:,:,k) )/GRAV ! kg/m2
       zdel(:,:,k) = zhalf(:,:,k) - zhalf(:,:,k+1) ! m
    end do
    if (Lice) then
      do k=1,kd
        precip3ds(:,:,k) = snow3d(:,:,k+1)-snow3d(:,:,k)
      end do
    else
      do k=1,kd
        precip3ds(:,:,k) = 0.
      end do
    endif
!
!+++ pag: 
!
if(lowercase(scheme)=='wdep_gas' .or. lowercase(scheme)=='wdep_aerosol') then
      ph = 5.0
!  cloud liquid water content
      xliq(:,:,:) = 0.5E-3                           !default value
      if(trim(cloud_param) .eq. 'convect') then
         xliq(:,:,:) = 1.0e-3
      elseif(trim(cloud_param) .eq. 'lscale') then
         xliq(:,:,:) = 0.5e-3
      endif

!                    
! tracer == gas      
!             
      if(lowercase(scheme)=="wdep_gas") THEN
          frac_int = frac_int_gas
          do k = 1, kd
            do j = 1, jd
              do i = 1, id
                Htemp(i,j,k)=henry_constant &
                          *exp(-Henry_variable*(1./298.-1./T(i,j,k)))
                K1(i,j)=1.2e-2*exp(-2010*(1/298.-1/T(i,j,k)))
                K2(i,j)=6.6e-8*exp(-1510*(1/298.-1/T(i,j,k)))
                HTemp(i,j,k)=Htemp(i,j,k)*(1 + K1(i,j) &
                           /10.**(-ph) + K1(i,j)*K2(i,j)/(10.**(-ph))**2)
                f_a(i,j) = Htemp(i,j,k)/101.325*RGAS &
                           *T(i,j,k)*xliq(i,j,k)*rho_air(i,j,k)/DENS_H2O
                scav_fact3d(i,j,k)=f_a(i,j)/(1.+f_a(i,j))
              enddo
            enddo
          enddo
      elseif(lowercase(scheme)=="wdep_aerosol") THEN
          frac_int = frac_int_aerosol
          scav_fact3d(:,:,:)=frac_in_cloud
      else
      print *,' Aerosol number =',n,' tracer_name=',tracer_name,' scheme=',text_in_scheme
      print *, 'Please check "am2p12.ft'
      call ERROR_MESG('wet_deposition', 'Tracer is neither aerosol NOR gas.', FATAL )
    endif
!                    
!in cloud scavenging
!
        do k=1,kd
          do j = 1, jd
            do i = 1, id
              beta(i,j)  = MAX( 0.0, precip3d(i,j,k)/pwt(i,j,k)/xliq(i,j,k))
              in_temp(i,j,k) = (exp(-beta(i,j)*scav_fact3d(i,j,k)*dt)-1.0)
              if ( tracer(i,j,k) .gt. 0.) then
                wdep_in(i,j)=wdep_in(i,j) &
                             - in_temp(i,j,k)*tracer(i,j,k)*pwt(i,j,k)
                tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                                   - in_temp(i,j,k)*tracer(i,j,k)/dt
              endif
!                    
!--reevaporation     
!    calculation of fracion of aerosols to-be
!    reevaporated to the atmosphere:

              beta(i,j)=precip3d(i,j,k)
              if (beta(i,j) < 0.) then
                beta(i,j) = beta(i,j)/rainsnow3d(i,j,k)
              endif
              if (rainsnow3d(i,j,k+1) == 0. ) then 
!--reevaporation total
                beta(i,j)=MIN(MAX(0.,-beta(i,j)),1.)
              else
                beta(i,j)=MIN(MAX(0.,-beta(i,j))*frac_int,1.)
              endif
! reevporating to atmosphere
              reevap_diag(i,j,k)=beta(i,j)*wdep_in(i,j)
              wdep_in(i,j) = wdep_in(i,j)*(1.-beta(i,j))
              tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                                 - reevap_diag(i,j,k)/pwt(i,j,k)/dt
            enddo
          enddo
        enddo
! Below cloud scavenging
        do k=1,kd
          do j = 1, jd
            do i = 1, id
              fluxs(i,j) = (snow3d(i,j,k+1)+snow3d(i,j,k))/2.0
              fluxr(i,j) = (rain3d(i,j,k+1)+rain3d(i,j,k))/2.0
              bc_temp(i,j,k) = 3./4.*dt* &
                       (fluxr(i,j)*alpha_r/R_r/DENS_H2O &
                        + fluxs(i,j)*alpha_s/R_s/DENS_SNOW)
              if ( tracer(i,j,k) .gt. 0. ) then
                wdep_bc(i,j)=wdep_bc(i,j) &
                             + bc_temp(i,j,k)*tracer(i,j,k)*pwt(i,j,k)
                tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                                   + bc_temp(i,j,k)*tracer(i,j,k)/dt
              endif
            enddo
          enddo
        enddo
!
!  end wdep_gas or wdep_aerosol
!
    else  

! Calculate fraction of precipitation reevaporated in layer
     do k=1,kd
       where( rainsnow3d(:,:,k) > 0. .and. precip3d(:,:,k) < 0. )
          reevap_fraction(:,:,k) = &
             -precip3d(:,:,k) / (rainsnow3d(:,:,k)) ! fraction
       end where
! Assume that the tracer reevaporation fraction is 50% of the precip
! reevaporation fraction, except when fraction = 100%
!      where( reevap_fraction(:,:,k) < 1. )
!         reevap_fraction(:,:,k) = 0.5*reevap_fraction(:,:,k)
!      end where
    end do
!  cloud liquid water content
!   xliq = 0.5E-3                           !default value
!   if(trim(cloud_param) .eq. 'convect') then
!      xliq = 1.0e-3
!   elseif(trim(cloud_param) .eq. 'lscale') then
!      xliq = 0.5e-3
!   endif

! Lgas = lowercase(scheme)=='henry' .or. lowercase(scheme)=='henry_below' .or. &
!        lowercase(scheme)=='henry_noice' .or. lowercase(scheme)=='henry_below_noice'
! Laerosol = lowercase(scheme)=='aerosol' .or. lowercase(scheme)=='aerosol_below' .or. &
!            lowercase(scheme)=='aerosol_noice' .or. lowercase(scheme)=='aerosol_below_noice'
! Assume that the aerosol reevaporation fraction is 50% of the precip
! reevaporation fraction, except when fraction = 100%
if( Lgas ) then
   frac_int = frac_int_gas
elseif( Laerosol ) then
   frac_int = frac_int_aerosol
else
   frac_int = 1.
end if

if( Lgas .or. Laerosol ) then
! units = VMR
!
! Henry_constant (mole/L/Pa) = [X](aq) / Px(g) 
! where [X](aq) is the concentration of tracer X in precipitation (mole/L)
!       Px(g) is the partial pressure of the tracer in the air (Pa)
!
! VMR (total) = VMR (gas) + VMR (aq)
!             = VMR (gas) + [X] * L
!
! where L = cloud liquid amount (kg H2O/mole air)
!
! Using Henry's Law, [X] = H * Px = H * VMR(gas) * Pfull
!
! So, VMR (total) =  VMR(gas) * [ 1 + H * Pfull * L ]
! 
! VMR(gas) = VMR(total) / [1 + H * Pfull * L]
!
! [X] = H * Pfull * VMR(total) / [ 1 + H * Pfull * L]
!
! Following Giorgi and Chameides, JGR, 90(D5), 1985, the first-order loss
! rate constant (s^-1) of X due to wet deposition equals:
!
! k = W_X / n_X
!
! where W_x = the loss rate (molec/cm3/s), and n_X = the number density (molec/cm3)
! 
! W_X = [X] * W_H2O / (55 mole/L)
! n_x = VMR(total) * n_air (molec/cm3) = VMR(total) * P/(kT) * 1E-6 m3/cm3
! 
! where P = atmospheric pressure (Pa)
!       k = Boltzmann's constant = 1.38E-23 J/K
!       T = temperature (K)
!       W_H2O = removal rate of water (molec/cm3/s)
! 
!             [X] * W_H2O / 55         
! So, k = ------------------------------
!         VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O    H * VMR(total) * P / [ 1 + H * P *L ]
!       = ----- * ---------------------------------------
!          55          VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O     H * kT * 1E6
!       = ----- *  -------------    
!          55      1 + H * P * L 
!
!         W_H2O     1     1     H * P * L
!       = ----- * ----- * - * -------------
!          55     n_air   L   1 + H * P * L
!
! where W_H2O = precip3d (kg/m2/s) * (AVOGNO/mw_h2o) (molec/kg) / zdel (m) * 1E-6 m3/cm3
!
   if( (Lgas .and. Henry_constant > 0) .or. Laerosol ) then
      in_temp(:,:,:) = 0.
      bc_temp(:,:,:) = 0.
      do k=1,kd
! Calculate the temperature dependent Henry's Law constant
         scav_factor(:,:) = 0.0
         xliq(:,:,k)  = MAX( cloud(:,:,k) * mw_air, 0. ) ! (kg H2O)/(mole air)
         n_air(:,:,k) = pfull(:,:,k) / (kboltz*T(:,:,k)) * cm3_2_m3 ! molec/cm3
         if (Lgas) then
            temp_factor(:,:) = 1/T(:,:,k)-inv298p15
            Htemp(:,:,k) = Henry_constant * &
                           exp( Henry_variable*temp_factor )
            f_a(:,:) = Htemp(:,:,k) * pfull(:,:,k) * xliq(:,:,k) ! / cloud_frac
            scav_factor(:,:) = f_a(:,:) / ( 1.+f_a(:,:) )
            scav_factor_s(:,:) = scav_factor(:,:)
         else if (Laerosol) then
            scav_factor(:,:) = frac_in_cloud
            if (frac_in_cloud == frac_in_cloud_snow) then
              scav_factor_s(:,:) = scav_factor(:,:)
            else
              scav_factor_s(:,:) = f_snow_berg(:,:,k)*frac_in_cloud_snow +&
                                   (1.-f_snow_berg(:,:,k))*frac_in_cloud
            endif 
         end if
!        where (precip3d(:,:,k) > 0.0)
         where (precip3d(:,:,k) > 0. .and. xliq(:,:,k) > 0.)
            w_h2o(:,:) = precip3d(:,:,k) * (AVOGNO/mw_h2o) / zdel(:,:,k) * cm3_2_m3 ! molec/cm3/s
            beta(:,:) = w_h2o(:,:) * mw_h2o  / (n_air(:,:,k) * xliq(:,:,k))
            where (precip3ds(:,:,k) > 0.0 .and. precip3dr(:,:,k) > 0.0)
              where (scav_factor(:,:) /= scav_factor_s(:,:)) 
                scav_factor(:,:) = ( scav_factor(:,:)*precip3dr(:,:,k) +  &
                   scav_factor_s(:,:) * precip3ds(:,:,k)) / precip3d(:,:,k)
              end where
            elsewhere
              where (precip3ds(:,:,k) > 0.0)
                 scav_factor(:,:) = scav_factor_s(:,:)
              endwhere
            endwhere
            in_temp(:,:,k) = beta(:,:) * scav_factor(:,:) ! 1/s
        endwhere
      enddo 
!-----------------------------------------------------------------
! Below-cloud wet scavenging
!-----------------------------------------------------------------
      if( lowercase(scheme)=='henry_below' .or. lowercase(scheme)=='henry_below_noice') then
         k_g = d_g/rain_diam * &
               ( 2. + 0.6 * sqrt( rain_diam*rain_vterm/vk_air ) * (vk_air/d_g)**(1./3.) )
         do i = 1,id
         do j = 1,jd
            conc(:) = tracer(i,j,:) * n_air(i,j,:) / cm3_2_m3 ! Convert from VMR to molec/m3
            do kk = 1,kd      
               stay = 1.
               if( precip3d(i,j,kk) > 0. ) then
                  conc_rain_total = 0.
                  stay = zfull(i,j,kk) / (rain_vterm * dt)
                  stay = min( stay, 1. )
                  do k = kk,kd
                     f_a0 = Htemp(i,j,k) * pfull(i,j,k) * xliq(i,j,kk) * n_air(i,j,kk)/n_air(i,j,k)
                     scav_factor0 = f_a0 / ( 1.+f_a0 )
                     conc_sat = conc(k) * scav_factor0 ! molec/m3 <== (xeqca1)
                     sa_drop0 = geo_fac / rain_diam * xliq(i,j,kk) * n_air(i,j,kk) / &
                                ( DENS_H2O * AVOGNO * cm3_2_m3 ) ! (m2 H2O) / (m3 air)
                     fgas0 = conc(k) * k_g ! molec/m2/s
                     fall_time = zdel(i,j,k) / rain_vterm ! sec
                     conc_rain = fgas0 * sa_drop0 * fall_time ! molec/m3 <== (xca1)
                     conc_rain_total = conc_rain_total + conc_rain ! molec/m3 <== (all1)
                     if ( conc_rain_total < conc_sat ) then
                        conc(k) = max( conc(k)-conc_rain, 0. )
                     end if
                  end do
                  conc(kk) = conc(kk) / n_air(i,j,kk) * cm3_2_m3 ! Convert to VMR
                  conc(kk) = tracer(i,j,kk) - conc(kk)
                  if ( conc(kk) /= 0. .and. tracer(i,j,kk) /= 0. ) then
                     fall_time = zdel(i,j,kk)/rain_vterm
                     bc_temp(i,j,kk) = bc_temp(i,j,kk) + &
                                       conc(kk) / (tracer(i,j,kk) * fall_time) * stay ! 1/s
                  end if
               end if         
            end do
         end do
         end do

      else if ( lowercase(scheme) == 'aerosol_below' .or. lowercase(scheme) == 'aerosol_below_noice') then

        do k=1,kd
           fluxs = (snow3d(:,:,k+1)+snow3d(:,:,k))/2.0
           fluxr = (rain3d(:,:,k+1)+rain3d(:,:,k))/2.0
           bc_temp(:,:,k) = 3./4. * &
                       (fluxr(:,:)*alpha_r/R_r/DENS_H2O + &
                        fluxs(:,:)*alpha_s/R_s/DENS_SNOW)
         end do

      end if


      do k = 1,kd
         wdep_in(:,:) = wdep_in(:,:) - &
                        in_temp(:,:,k)*tracer(:,:,k)*pwt(:,:,k)*  &
                                                         cloud_frac(:,:,k)
         wdep_bc(:,:) = wdep_bc(:,:) - &
                        bc_temp(:,:,k)*tracer(:,:,k)*pwt(:,:,k)
      enddo
      dt_temp(:,:,:) = 1. - exp( -bc_temp(:,:,:)*dt ) & ! fractional loss/timestep
                          * ( cloud_frac(:,:,:)*exp( -in_temp(:,:,:)*dt ) + (1-cloud_frac(:,:,:)) )
      tracer_dt(:,:,:) = dt_temp(:,:,:) / dt !+ve loss frequency (1/sec)
    endif

else if(lowercase(scheme)=='fraction') then
   tracer_dt = 0.0
!-----------------------------------------------------------------------
!
!     Compute areal fractions experiencing wet deposition:
!
!     Set minimum precipitation rate below which no wet removal
!     occurs to 0.01 cm/day ie 1.16e-6 mm/sec (kg/m2/s)
   premin=1.16e-6
!
!     Large scale cloud liquid water content (kg/m3)
!     and below cloud washout efficiency (cm-1):
   flaglw =parse(control,'lslwc',clwc)
   if (flaglw == 0 ) clwc=0.5e-3
   wash=1.0  
!
!     When convective adjustment occurs, use convective cloud liquid water content:
!
    if(trim(cloud_param) .eq. 'convect') then
      flaglw = parse(control,'convlwc',clwc)
      if (flaglw == 0) clwc=2.0e-3
      wash=0.3 
   end if
!
   do j=1,size(rain,2)
   do i=1,size(rain,1)
      tracer_dt(i,j,:)=0.0
      washout(i,j)=0.0
      prenow = rain(i,j) + snow(i,j)
      if(prenow .gt. premin) then      
!
! Assume that the top of the cloud is where the highest model level 
! specific humidity is reduced. And the the bottom of the cloud is the
! lowest model level where specific humidity is reduced.
!
         ktopcd(i,j) = 0
         do k = kd,1,-1
            if (qdt(i,j,k) < 0.0 ) ktopcd(i,j) = k
         enddo
         kendcd(i,j) = 0
         do k = 1,kd
            if (qdt(i,j,k) < 0.0 ) kendcd(i,j) = k
         enddo
!
!     Thickness of precipitating cloud deck:
!
         if(ktopcd(i,j).gt.1) then
            hwtop = 0.0
            do k=ktopcd(i,j),kendcd(i,j)
               hwtop=hwtop+(phalf(i,j,k+1)-phalf(i,j,k))*rdgas*T(i,j,k)/grav/pfull(i,j,k)
            enddo
            do k=ktopcd(i,j),kendcd(i,j)
!     Areal fraction affected by precip clouds (max = 0.5):
               tracer_dt(i,j,k)=prenow/(clwc*hwtop)
            end do  
          endif

         washout(i,j)=prenow*wash
          endif
   end do
   end do
          endif

! Now multiply by the tracer mixing ratio to get the actual tendency.
tracer_dt(:,:,:) = MIN( MAX(tracer_dt(:,:,:), 0.0E+00), 0.5/dt)
where (tracer > 0.)
   tracer_dt = tracer_dt*tracer
else where
   tracer_dt = 0.
end where

!++lwh
!
! Re-evaporation
!
do k = 1,kd
   where (reevap_fraction(:,:,k) > 0.) 
      reevap_diag(:,:,k) = reevap_fraction(:,:,k) * tracer_flux(:,:)
! tracer reevaporation fraction is reduced from precip reevaporation,
! except when complete reevaporation occurs
      where( reevap_fraction(:,:,k) < 1. )
         reevap_diag(:,:,k) = reevap_diag(:,:,k) * frac_int
      end where
      tracer_dt(:,:,k) = tracer_dt(:,:,k) - reevap_diag(:,:,k) / pwt(:,:,k)
   end where
   tracer_flux(:,:) = tracer_flux(:,:) + tracer_dt(:,:,k)*pwt(:,:,k)
end do
!--lwh
!
endif ! End branching pag/lwh
!
! Output diagnostics in kg/m2/s (if MMR) or mole/m2/s (if VMR)
if(trim(units) .eq. 'mmr') then
   diag_scale = 1.
        elseif(trim(units) .eq. 'vmr') then
   diag_scale = mw_air ! kg/mole
else
   write(*,*) ' Tracer number =',n,' tracer_name=',tracer_name
   write(*,*) ' scheme=',text_in_scheme
   write(*,*) ' control=',control
   write(*,*) ' scheme=',scheme
   write(*,*) 'Please check field table'
   write(*,*) 'tracers units =',trim(units),'it should be either  mmr or vmr!'
!  <ERROR MSG="Unsupported tracer units" STATUS="FATAL">
!     Tracer units must be either VMR or MMR
!  </ERROR>
   call error_mesg('wet_deposition', 'Unsupported tracer units.', FATAL )
          endif


! Column integral of wet deposition
sum_wdep = 0.
do k=1,kd
   sum_wdep = sum_wdep + tracer_dt(:,:,k)*pwt(:,:,k)/diag_scale
end do

if (present (sum_wdep_out))  sum_wdep_out = -sum_wdep

if(trim(cloud_param) == 'lscale') then
   if (id_tracer_reevap_ls(n) > 0 ) then
      used = send_data ( id_tracer_reevap_ls(n), reevap_diag/diag_scale, Time ,is,js,1)
   endif
   if (id_tracer_wdep_ls(n) > 0 ) then
      used = send_data ( id_tracer_wdep_ls(n), sum_wdep, Time, is_in =is, js_in=js )
   endif
   if (id_tracer_wdep_ls_3d(n) > 0 ) then
     used = send_data ( id_tracer_wdep_ls_3d(n), tracer_dt*pwt/diag_scale, Time, is_in =is, js_in=js ,ks_in=1)
   endif
   if (id_tracer_wdep_lsin(n) > 0 ) then
       used = send_data ( id_tracer_wdep_lsin(n), wdep_in/diag_scale, Time, is_in=is, js_in=js )
   endif
   if (id_tracer_wdep_lsbc(n) > 0 ) then
      used = send_data ( id_tracer_wdep_lsbc(n), wdep_bc/diag_scale, Time, is_in=is, js_in=js )
   endif

else if(trim(cloud_param) == 'convect') then
   if (id_tracer_reevap_cv(n) > 0 ) then
      used = send_data ( id_tracer_reevap_cv(n), reevap_diag/diag_scale, Time ,is,js,1)
   endif
   if(id_tracer_wdep_cv(n) > 0) then
      used = send_data( id_tracer_wdep_cv(n), sum_wdep, Time, is_in=is, js_in=js)
   endif
   if (id_tracer_wdep_cvin(n) > 0 ) then
       used = send_data ( id_tracer_wdep_cvin(n), wdep_in/diag_scale, Time, is_in=is, js_in=js)
   endif
   if (id_tracer_wdep_cvbc(n) > 0 ) then
      used = send_data ( id_tracer_wdep_cvbc(n), wdep_bc/diag_scale, Time, is_in=is, js_in=js)
   endif
endif

end subroutine wet_deposition
!</SUBROUTINE>
!
!#######################################################################
!
subroutine get_drydep_param(text_in_scheme,text_in_param,scheme,land_dry_dep_vel,sea_dry_dep_vel)
!subroutine get_drydep_param(text_in_scheme, text_in_param, scheme,  &
!                            land_dry_dep_vel, sea_dry_dep_vel, &
!                            ice_dry_dep_vel, snow_dry_dep_vel,  &
!                            vegn_dry_dep_vel)
!
! Subroutine to initialiize the parameters for the dry deposition scheme.
! If the dry dep scheme is a "fixed" form then the 
! dry_deposition velocity value has to be set.
! If the dry dep scheme is a "wind_driven" form then the dry_deposition
! velocity value will be calculated. So set to a dummy value of 0.0
! INTENT IN
!  text_in_scheme   : The text that has been parsed from tracer table as 
!                     the dry deposition scheme to be used.
!  text_in_param    : The parameters that are associated with the dry 
!                     deposition scheme.
! INTENT OUT
!  scheme           : The scheme that is being used.
!  land_dry_dep_vel : Dry deposition velocity over the land
!  sea_dry_dep_vel  : Dry deposition velocity over the sea
!
character(len=*), intent(in)    :: text_in_scheme, text_in_param
character(len=*), intent(out)   :: scheme
real, intent(out)               :: land_dry_dep_vel, sea_dry_dep_vel!, &
!                                   ice_dry_dep_vel, snow_dry_dep_vel, &
!                                   vegn_dry_dep_vel

integer :: flag

!Default
scheme                  = 'None'
land_dry_dep_vel=0.0
sea_dry_dep_vel=0.0
!ice_dry_dep_vel=0.0
!snow_dry_dep_vel=0.0
!vegn_dry_dep_vel=0.0

if(lowercase(trim(text_in_scheme(1:4))).eq.'wind') then
 scheme                  = 'Wind_driven'
endif

!if(lowercase(trim(text_in_scheme(1:15))).eq.'sfc_dependent_w') then
!scheme                 = 'sfc_dependent_wind_driven'
!endif
!
!if(lowercase(trim(text_in_scheme(1:6))).eq.'sfc_bl') then
!scheme                 = 'sfc_BL_dependent_wind_driven'
!endif

if(lowercase(trim(text_in_scheme(1:5))).eq.'fixed') then
scheme                 = 'fixed'
flag=parse(text_in_param,'land',land_dry_dep_vel)
flag=parse(text_in_param,'sea', sea_dry_dep_vel)
endif

!if(lowercase(trim(text_in_scheme(1:15))).eq.'sfc_dependent_f') then
!scheme                 = 'sfc_dependent_fixed'
!flag=parse(text_in_param,'land',land_dry_dep_vel)
!flag=parse(text_in_param,'sea', sea_dry_dep_vel)
!flag=parse(text_in_param,'ice', ice_dry_dep_vel)
!if (flag == 0) then
!   ice_dry_dep_vel = sea_dry_dep_vel
!end if
!flag=parse(text_in_param,'snow', snow_dry_dep_vel)
!if (flag == 0) then
!   snow_dry_dep_vel = land_dry_dep_vel
!end if
!flag=parse(text_in_param,'vegn', vegn_dry_dep_vel)
!if (flag == 0) then
!   vegn_dry_dep_vel = land_dry_dep_vel
!end if
!endif

if(lowercase(trim(text_in_scheme(1:4))).eq.'file') then
   scheme = 'file'
endif

end subroutine get_drydep_param
!
!#######################################################################
!
!<SUBROUTINE NAME="get_wetdep_param">
!<TEMPLATE>
!CALL get_wetdep_param(text_in_scheme, text_in_param, scheme,&
!                      henry_constant, henry_temp, &
!                      frac_in_cloud, alpha_r, alpha_s)
!</TEMPLATE>
subroutine get_wetdep_param(text_in_scheme,text_in_param,scheme,&
                            henry_constant, henry_temp, &
                            frac_in_cloud, frac_in_cloud_snow,  &
                            alpha_r,alpha_s, &
                            Lwetdep, Lgas, Laerosol, Lice, &
                            frac_in_cloud_uw, frac_in_cloud_donner)
!<OVERVIEW>
! Routine to initialize the parameters for the wet deposition scheme.
!</OVERVIEW>
!
! shm has modified this subroutine to include additional parameters:
!      frac_in_cloud, alpha_r, and alpha_s
!
! INTENT IN
!<IN NAME="text_in_scheme" TYPE="character">
!   Text read from the tracer table which provides information on which
!                   wet deposition scheme to use.
!</IN>
!<IN NAME="text_in_param" TYPE="character">
!   Parameters associated with the wet deposition scheme. These will be
!                   parsed in this routine.
!</IN>
!<OUT NAME="scheme" TYPE="character">
!   Wet deposition scheme to use.
!   Choices are: None, Fraction, Henry, Henry_below, Aerosol, Aerosol_below,
!                wdep_aerosol, wdep_gas
!</OUT>
!<OUT NAME="henry_constant" TYPE="real">
!   Henry's Law constant for the tracer (see wet_deposition for explanation of Henry's Law)
!</OUT>
!<OUT NAME="henry_temp" TYPE="real">
!   The temperature dependence of the Henry's Law constant.
!</OUT>
!<OUT NAME="frac_in_cloud" TYPE="real">
!   In-cloud fraction for aerosols
!</OUT>
!<OUT NAME="alpha_r" TYPE="real">
!   Controls below-cloud aerosol scavenging by rain
!</OUT>
!<OUT NAME="alpha_s" TYPE="real">
!   Controls below-cloud aerosol scavenging by snow
!</OUT>
!<OUT NAME="Lwetdep" TYPE="logical">
!   Does tracer have wet removal?
!</OUT>
!<OUT NAME="Lgas" TYPE="logical">
!   Is tracer a gas?
!</OUT>
!<OUT NAME="Laerosol" TYPE="logical">
!   Is tracer an aerosol?
!</OUT>
!<OUT NAME="Lice" TYPE="logical">
!   Is tracer removed by snow (or just rain)?
!</OUT>


character(len=*), intent(in)    :: text_in_scheme, text_in_param
character(len=*), intent(out)   :: scheme
real, intent(out)               :: henry_constant, henry_temp
real, intent(out)               :: frac_in_cloud, frac_in_cloud_snow
real, intent(out)               :: alpha_r, alpha_s
logical, intent(out)            :: Lwetdep, Lgas, Laerosol, Lice
real, intent(out), optional     :: frac_in_cloud_uw, frac_in_cloud_donner

integer :: flag

!Default
scheme                  = 'None'
henry_constant= 0.
henry_temp    = 0.
frac_in_cloud = 0.
frac_in_cloud_snow = 0.
alpha_r       = 0.
alpha_s       = 0.
Lwetdep = .false.
Lgas = .false.
Laerosol = .false.

if (present(frac_in_cloud_uw))     frac_in_cloud_uw = 0.
if (present(frac_in_cloud_donner)) frac_in_cloud_donner = 0.

if( trim(lowercase(text_in_scheme)) == 'fraction' ) then
   scheme                 = 'Fraction'
else if( trim(lowercase(text_in_scheme)) == 'henry' .or. &
         trim(lowercase(text_in_scheme)) == 'henry_below' .or. &
         trim(lowercase(text_in_scheme)) == 'henry_noice' .or. &
         trim(lowercase(text_in_scheme)) == 'henry_below_noice' ) then
   if( trim(lowercase(text_in_scheme)) == 'henry' ) then
      scheme                 = 'henry'
   else if ( trim(lowercase(text_in_scheme)) == 'henry_below' ) then
      scheme                 = 'henry_below'
   else if ( trim(lowercase(text_in_scheme)) == 'henry_noice' ) then
      scheme                 = 'henry_noice'
   else if ( trim(lowercase(text_in_scheme)) == 'henry_below_noice' ) then
      scheme                 = 'henry_below_noice'
   end  if
   flag=parse(text_in_param,'henry',     henry_constant)
   flag=parse(text_in_param,'dependence',henry_temp    )
   Lgas = .true.
else if( trim(lowercase(text_in_scheme)) == 'aerosol' .or. &
         trim(lowercase(text_in_scheme)) == 'aerosol_below' .or. &
         trim(lowercase(text_in_scheme)) == 'aerosol_noice' .or. &
         trim(lowercase(text_in_scheme)) == 'aerosol_below_noice' ) then
   if( trim(lowercase(text_in_scheme)) == 'aerosol' ) then
      scheme                 = 'aerosol'
   else if ( trim(lowercase(text_in_scheme)) == 'aerosol_below' ) then
      scheme                 = 'aerosol_below'
   else if ( trim(lowercase(text_in_scheme)) == 'aerosol_noice' ) then
      scheme                 = 'aerosol_noice'
   else if ( trim(lowercase(text_in_scheme)) == 'aerosol_below_noice' ) then
      scheme                 = 'aerosol_below_noice'
   end if
   flag=parse(text_in_param,'frac_incloud',frac_in_cloud)

   flag=parse(text_in_param,'frac_incloud_snow',frac_in_cloud_snow)
   if (flag == 0) then
      frac_in_cloud_snow = frac_in_cloud
   end if

   if (present(frac_in_cloud_uw)) then
      flag=parse(text_in_param,'frac_incloud_uw',frac_in_cloud_uw)
      if (flag == 0) then
         frac_in_cloud_uw = frac_in_cloud
      end if
   end if
   if (present(frac_in_cloud_donner)) then
      flag=parse(text_in_param,'frac_incloud_donner',   &
                                                frac_in_cloud_donner)
      if (flag == 0) then
         frac_in_cloud_donner = frac_in_cloud
      end if
   end if
   flag=parse(text_in_param,'alphar',alpha_r)
   flag=parse(text_in_param,'alphas',alpha_s)
   Laerosol = .true.
end if
if( trim(lowercase(text_in_scheme)) == 'wdep_aerosol') scheme= 'wdep_aerosol'
if ( trim(lowercase(text_in_scheme)) == 'wdep_gas' ) scheme= 'wdep_gas'
if (scheme .eq. 'wdep_aerosol' .or. scheme .eq. 'wdep_gas') then
   flag=parse(text_in_param,'frac_incloud',frac_in_cloud)
   flag=parse(text_in_param,'alphar',alpha_r)
   flag=parse(text_in_param,'alphas',alpha_s)
end if

Lice = .not. ( scheme=='henry_noice' .or. scheme=='henry_below_noice' .or. &
               scheme=='aerosol_noice' .or. scheme=='aerosol_below_noice' )
Lwetdep = scheme /= 'None'

end subroutine get_wetdep_param
!</SUBROUTINE>
!
!#######################################################################
!
!<SUBROUTINE NAME="interp_emiss">
subroutine interp_emiss(global_source, start_lon, start_lat, &
                        lon_resol, lat_resol, data_out)
!
!<OVERVIEW>
! A routine to interpolate emission fields of arbitrary resolution onto the 
! resolution of the model.
!</OVERVIEW>
!<DESCRIPTION>
! Routine to interpolate emission fields (or any 2D field) to the model 
! resolution. The local section of the global field is returned to the 
! local processor.
!</DESCRIPTION>
! 
!<TEMPLATE>
! call interp_emiss(global_source, start_lon, start_lat, &
!                        lon_resol, lat_resol, data_out)
!</TEMPLATE>
! INTENT IN
!<IN NAME="global_source" TYPE="real" DIM="(:,:)">
!  Global emission field.
!</IN>
!<IN NAME="start_lon" TYPE="real">
!  Longitude of starting point of emission field 
!  (in radians). This is the westernmost boundary of the 
!  global field.
!</IN>
!<IN NAME="start_lat" TYPE="real">
!  Latitude of starting point of emission field
!  (in radians). This is the southern boundary of the 
!  global field.
!</IN>
!<IN NAME="lon_resol" TYPE="real">
!  Longitudinal resolution of the emission data (in radians).
!</IN>
!<IN NAME="lat_resol" TYPE="real">
!  Latitudinal resolution of the emission data (in radians).
!</IN>
! 
! INTENT OUT
!<OUT NAME="data_out" TYPE="real" DIM="(:,:)">
!  Interpolated emission field on the local PE. 
!</OUT>

real, intent(in)  :: global_source(:,:)
real, intent(in)  :: start_lon,start_lat,lon_resol,lat_resol
real, intent(out) :: data_out(:,:)

integer :: i, j, nlon_in, nlat_in
real :: blon_in(size(global_source,1)+1)
real :: blat_in(size(global_source,2)+1)
type (horiz_interp_type) :: Interp
! Set up the global surface boundary condition longitude-latitude boundary values

   nlon_in = size(global_source,1)
   nlat_in = size(global_source,2)
! For some reason the input longitude needs to be incremented by 180 degrees.
   do i = 1, nlon_in+1
      blon_in(i) = start_lon + float(i-1)*lon_resol + PI
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-twopi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+twopi

   do j = 2, nlat_in
      blat_in(j) = start_lat + float(j-1)*lat_resol
   enddo
      blat_in(1)         = -0.5*PI
      blat_in(nlat_in+1) =  0.5*PI

! Now interpolate the global data to the model resolution
   call horiz_interp_init
   call horiz_interp_new (Interp, blon_in, blat_in, &
                                  blon_out, blat_out)
   call horiz_interp (Interp, global_source, data_out)
   call horiz_interp_del ( Interp )


end subroutine interp_emiss
!</SUBROUTINE>
! ######################################################################
!

      subroutine GET_RH (T,Q,P,RH,mask)
!***********************************************************************
!  SUBROUTINE GET_RH
!  PURPOSE
!     VECTOR COMPUTATION OF RELATIVE HUMIDITY
!  DESCRIPTION OF PARAMETERS
!     T        TEMPERATURE VECTOR (DEG K)
!     P        PRESSURE VECTOR (Pa)
!     Q        SPECIFIC HUMIDITY (kg/kg)
!     RH       RELATIVE HUMIDITY
!     MASK     EARTH SURFACE BELOW GROUND (in pressure coordinates)
!
!***********************************************************************
!

      Real, parameter :: ONE    = 1.
      Real, parameter :: ZP622  = 0.622
      Real, parameter :: Z1P0S1 = 1.00001
      Real, parameter :: Z1P622 = 1.622
      Real, parameter :: Z138P9 = 138.90001
      Real, parameter :: Z198P9 = 198.99999
      Real, parameter :: Z200   = 200.0
      Real, parameter :: Z337P9 = 337.9

      real, intent(in), dimension(:,:,:) :: T !temp at curr time step [ deg k ]
      real, intent(in), dimension(:,:,:) :: P !pressure at full levels [ Pa ]
      real, intent(in), dimension(:,:,:) :: Q !specific humidity at current time step  [ kg / kg ]
      real, intent(in), dimension(:,:,:), optional :: mask !
      real, intent(out), dimension(:,:,:) :: RH !relative humidity [0-1]

! Dynamic Work Space
! ------------------
      real :: A1622
      real, dimension(size(q,1),size(q,2),size(q,3)) :: e1, e2, tq, qs
      integer, dimension(size(q,1),size(q,2),size(q,3)) :: i1, i2
      integer :: i,j,k

!
      real, dimension(67) :: EST1
      data EST1/       0.31195E-02, 0.36135E-02, 0.41800E-02, &
          0.48227E-02, 0.55571E-02, 0.63934E-02, 0.73433E-02, &
          0.84286E-02, 0.96407E-02, 0.11014E-01, 0.12582E-01, &
          0.14353E-01, 0.16341E-01, 0.18574E-01, 0.21095E-01, &
          0.23926E-01, 0.27096E-01, 0.30652E-01, 0.34629E-01, &
          0.39073E-01, 0.44028E-01, 0.49546E-01, 0.55691E-01, &
          0.62508E-01, 0.70077E-01, 0.78700E-01, 0.88128E-01, &
          0.98477E-01, 0.10983E+00, 0.12233E+00, 0.13608E+00, &
          0.15121E+00, 0.16784E+00, 0.18615E+00, 0.20627E+00, &
          0.22837E+00, 0.25263E+00, 0.27923E+00, 0.30838E+00, &
          0.34030E+00, 0.37520E+00, 0.41334E+00, 0.45497E+00, &
          0.50037E+00, 0.54984E+00, 0.60369E+00, 0.66225E+00, &
          0.72589E+00, 0.79497E+00, 0.86991E+00, 0.95113E+00, &
          0.10391E+01, 0.11343E+01, 0.12372E+01, 0.13484E+01, &
          0.14684E+01, 0.15979E+01, 0.17375E+01, 0.18879E+01, &
          0.20499E+01, 0.22241E+01, 0.24113E+01, 0.26126E+01, &
          0.28286E+01, 0.30604E+01, 0.33091E+01, 0.35755E+01/
!
      real, dimension(72)  :: EST2
      data EST2/ &
          0.38608E+01, 0.41663E+01, 0.44930E+01, 0.48423E+01, &
          0.52155E+01, 0.56140E+01, 0.60394E+01, 0.64930E+01, &
          0.69767E+01, 0.74919E+01, 0.80406E+01, 0.86246E+01, &
          0.92457E+01, 0.99061E+01, 0.10608E+02, 0.11353E+02, &
          0.12144E+02, 0.12983E+02, 0.13873E+02, 0.14816E+02, &
          0.15815E+02, 0.16872E+02, 0.17992E+02, 0.19176E+02, &
          0.20428E+02, 0.21750E+02, 0.23148E+02, 0.24623E+02, &
          0.26180E+02, 0.27822E+02, 0.29553E+02, 0.31378E+02, &
          0.33300E+02, 0.35324E+02, 0.37454E+02, 0.39696E+02, &
          0.42053E+02, 0.44531E+02, 0.47134E+02, 0.49869E+02, &
          0.52741E+02, 0.55754E+02, 0.58916E+02, 0.62232E+02, &
          0.65708E+02, 0.69351E+02, 0.73168E+02, 0.77164E+02, &
          0.81348E+02, 0.85725E+02, 0.90305E+02, 0.95094E+02, &
          0.10010E+03, 0.10533E+03, 0.11080E+03, 0.11650E+03, &
          0.12246E+03, 0.12868E+03, 0.13517E+03, 0.14193E+03, &
          0.14899E+03, 0.15634E+03, 0.16400E+03, 0.17199E+03, &
          0.18030E+03, 0.18895E+03, 0.19796E+03, 0.20733E+03, &
          0.21708E+03, 0.22722E+03, 0.23776E+03, 0.24871E+03/
!
      real, dimension(139) :: EST
      EQUIVALENCE (EST(1)  , EST1(1)), (EST(68),EST2(1))
!***********************************************************************
!
      A1622   = ONE  / Z1P622
      TQ = T - Z198P9
      I1(:,:,:) = 1
      I2(:,:,:) = 1
      where ( T < Z200 ) TQ = Z1P0S1
      where ( T > Z337P9 ) TQ = Z138P9
      IF (present(mask)) THEN
        where ( mask > 0. )
          I1 = int(TQ)
          I2 = I1 + 1
        end where
      else
        I1 = int(TQ)
        I2 = I1 + 1
      endif
      do i=1,size(q,1)
        do j=1,size(q,2)
          do k=1,size(q,3)
            E1(i,j,k) =  EST( I1(i,j,k) )
            E2(i,j,k) =  EST( I2(i,j,k) )
          enddo
        enddo
     enddo
      QS(:,:,:) = TQ(:,:,:) - float(I1(:,:,:))
      QS(:,:,:) = E1(:,:,:) + QS(:,:,:) * ( E2(:,:,:)-E1(:,:,:) )
      E1(:,:,:) = (0.01 * P(:,:,:)) * A1622
      where ( E1 < QS ) QS = E1
      if (present(mask)) then
        where ( mask > 0. )  QS = ZP622 * QS / ( P * 0.01)
      else
        QS(:,:,:) = ZP622 * QS(:,:,:) / ( P(:,:,:) * 0.01)
      endif
      RH(:,:,:) = Q(:,:,:)/QS(:,:,:)

end subroutine GET_RH

! ######################################################################
!
subroutine get_w10m(z_full, u, v, rough_mom,u_star, b_star, q_star, &
       w10m_ocean, w10m_land, Time_next, is,js)

real, intent(in),    dimension(:,:) :: z_full, u, v
real, intent(in),    dimension(:,:)   :: rough_mom
real, intent(in),    dimension(:,:)   :: u_star, b_star, q_star
type(time_type), intent(in)           :: Time_next
integer, intent(in)                   :: is,js

logical :: used

real, intent(out),   dimension(:,:)   :: w10m_ocean, w10m_land
real, dimension(size(u,1),size(u,2)) ::  del_m
real, dimension(size(u,1),size(u,2)) ::  del_h
real, dimension(size(u,1),size(u,2)) ::  del_q
! Reference heights for momentum and heat [m]
real, parameter :: zrefm = 10.
real, parameter :: zrefh = 2.
real, parameter :: scaling_factor=1.

     w10m_ocean(:,:)    = 0.0
     w10m_land (:,:)    = 0.0
     del_m(:,:)   = 0.0

      call mo_profile(zrefm, zrefh, z_full, &
           rough_mom, rough_mom, rough_mom, &
           u_star, b_star, q_star, &
           del_m, del_h, del_q )
!-----------------------------------------------------------------
!       ... Wind speed at anemometer level (10 meters above ground)
!-----------------------------------------------------------------
      w10m_ocean(:,:)=sqrt(u(:,:)**2 +v(:,:)**2 )*del_m(:,:)
      w10m_land (:,:)=sqrt(u(:,:)**2 +v(:,:)**2 )*del_m(:,:)*scaling_factor

! Send the scaling factor
      if (id_delm > 0 ) then
        used = send_data ( id_delm, del_m, Time_next, is_in=is,js_in=js )
      endif

! Send the 10m wind speed data to the diag_manager for output.
      if (id_w10m > 0 ) then
        used = send_data ( id_w10m, w10m_land, Time_next, is_in=is,js_in=js )
      endif

end subroutine get_w10m

! ######################################################################
!
subroutine get_cldf(ps, pfull, rh, cldf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This subroutine estimates the cloud fraction "cldf" for
!!! each grid box using an empirical function of the relative
!!! humidity in that grid box, following Sundqvist et al., Mon. Weather Rev.,
!!! v117, 164101657, 1989:
!!!
!!!             cldf = 1 - sqrt[ 1 - (RH - RH0)/(1 - RH0) ]
!!!
!!! where RH is the relative humidity and RH0 is the threshold relative
!!! humidity for condensation specified as a function of pressure based
!!! on Xu and Krueger, Mon. Weather Rev., v119, 342-367, 1991.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real, intent(in),  dimension(:,:)   :: ps
      real, intent(in),  dimension(:,:,:) :: pfull
      real, intent(in),  dimension(:,:,:) :: rh
      real, intent(out), dimension(:,:,:) :: cldf
      real, parameter :: zrt = 0.6
      real, parameter :: zrs = 0.99
      integer           :: i,j,k, id, jd, kd
      real              :: p, r, r0, b0
      id=size(pfull,1); jd=size(pfull,2); kd=size(pfull,3)

      do k = 1, kd
      do j = 1, jd
      do i = 1, id
       P = pfull(i,j,k)
       R = RH(i,j,k)
       R0 = ZRT + (ZRS-ZRT) * exp(1.-(PS(i,j)/P)**2.5)
       B0 = (R-R0) / (1.-R0)
       if (R .lt.R0) B0 = 0.
       if (B0.gt.1.) B0 = 1.
       CLDF(i,j,k) = 1.-sqrt(1.-B0)
      end do
      end do
      end do

end subroutine get_cldf

!######################################################################
!<SUBROUTINE NAME="tracer_utilities_end">
!<OVERVIEW>
!  The destructor routine for the tracer utilities module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>

subroutine atmos_tracer_utilities_end
 

   deallocate(blon_out, blat_out)
   module_is_initialized = .FALSE.

 end subroutine atmos_tracer_utilities_end
!</SUBROUTINE>

! ######################################################################
! !IROUTINE: sjl_fillz --- Fill from neighbors below and above
!
! !INTERFACE:
 subroutine sjl_fillz(im, km, nq, q, dp)

 implicit none

! !INPUT PARAMETERS:
   integer,  intent(in):: im                ! No. of longitudes
   integer,  intent(in):: km                ! No. of levels
   integer,  intent(in):: nq                ! Total number of tracers

   real, intent(in)::  dp(im,km)       ! pressure thickness
! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout) :: q(im,km,nq)   ! tracer mixing ratio

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !BUGS:
!   Currently this routine only performs the east-west fill algorithm.
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   00.04.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i, k, ic
   real qup, qly, dup

   do ic=1,nq
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0.) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( 0.75*qly, qup )        !borrow no more than 75% from top
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1)
             q(i,k  ,ic) = 0.
          endif
          enddo
      enddo

! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic) < 0.) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( qly, qup )
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
             q(i,k,ic) = 0.
          endif
      enddo
   enddo
end subroutine sjl_fillz

end module atmos_tracer_utilities_mod
