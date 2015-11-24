module atmos_ch3i_mod


use mpp_mod, only: input_nml_file 
use            fms_mod, only : file_exist,   &
                               write_version_number, &
                               error_mesg, &
                               FATAL, &
                               NOTE, &
                               mpp_pe,  &
                               mpp_root_pe, &
                               lowercase,   &
                               open_namelist_file, &
                               check_nml_error, &
                               close_file,   &
                               stdlog
use  field_manager_mod, only : MODEL_ATMOS,          &
                               parse
use tracer_manager_mod, only : get_tracer_index,     &
                               get_tracer_names,     &
                               query_method
use   time_manager_mod, only : time_type
use   diag_manager_mod, only : send_data,            &
                               register_diag_field
use   interpolator_mod, only : interpolate_type,     &
                               interpolator_init,    &
                               obtain_interpolator_time_slices, &
                               unset_interpolator_time_flag, &
                               interpolator_end,     &
                               interpolator,         &
                               query_interpolator,   &
                               CONSTANT,             &
                               INTERP_WEIGHTED_P  
use      constants_mod, only : WTMAIR, &
                               AVOGNO, &
                               SECONDS_PER_DAY
implicit none

private
public :: atmos_ch3i_init, atmos_ch3i_time_vary, atmos_ch3i,  &
          atmos_ch3i_endts, atmos_ch3i_end

!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
character(len=128) :: conc_filename = ''

namelist /atmos_ch3i_nml/    &
   conc_filename

logical :: has_emissions                      
character(len=128) :: emis_filename
type(interpolate_type), save :: ch3i_emissions, input_conc

character(len=128), allocatable :: field_names(:)

integer :: id_emissions, id_loss, id_j_ch3i, ind_ch3i
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4        !conversion factor (cm2/m2)
real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO        
real, parameter :: boltz = 1.38044e-16      ! Boltzmann's Constant (erg/K)


character(len=7), parameter :: module_name = 'tracers'
!---- version number -----
character(len=128) :: version = '$Id: atmos_ch3i.F90,v 20.0 2013/12/13 23:23:47 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
logical :: module_is_initialized = .FALSE.

contains

subroutine atmos_ch3i_init( lonb_mod, latb_mod, axes, Time, mask )

!-----------------------------------------------------------------------
!
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   real, intent(in), dimension(:,:) :: lonb_mod
   real, intent(in), dimension(:,:) :: latb_mod
   type(time_type), intent(in) :: Time
   integer        , intent(in) :: axes(4)
   real, intent(in),    dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

   integer ::  unit, nfields, flag_file
   integer :: ierr, io, logunit

   character(len=128) :: tracer_name, tracer_units, name, control

   if (module_is_initialized) then
      return
   end if

   ind_ch3i = get_tracer_index(MODEL_ATMOS, 'ch3i')

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)

!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=atmos_ch3i_nml, iostat=io)
     ierr = check_nml_error(io, 'atmos_ch3i_nml')
#else
     unit = open_namelist_file ()
     ierr=1; do while (ierr /= 0)
     read  (unit, nml=atmos_ch3i_nml, iostat=io, end=10)
     ierr = check_nml_error(io, 'atmos_ch3i_nml')
     enddo
10   call close_file (unit)
#endif
   endif

  
   logunit=stdlog()
   if(mpp_pe() == mpp_root_pe()) then       
      write(logunit, nml=atmos_ch3i_nml)
   endif
     
!-----------------------------------------------------------------------
!     ... set up emissions
!-----------------------------------------------------------------------
   has_emissions = .false.
   if( query_method('emissions',MODEL_ATMOS,ind_ch3i,name,control) ) then
      if( trim(name) == 'file' ) then
         flag_file = parse(control, 'file',emis_filename)
         if (flag_file > 0) then
            has_emissions = .true.
            call interpolator_init( ch3i_emissions, trim(emis_filename), &
                                    lonb_mod, latb_mod,  &
                                    data_out_of_bounds=(/CONSTANT/) )
            call query_interpolator(ch3i_emissions,nfields=nfields)
            allocate( field_names(nfields) )
            call query_interpolator(ch3i_emissions,field_names=field_names)
            id_emissions = register_diag_field( module_name, 'ch3i_emis', axes(1:2), &
                                                Time, 'ch3i_emis', 'molec/cm2/s' )
         else
            call error_mesg('atmos_ch3i_init','CH3I emission file not specified in field table',FATAL)
         end if
      else 
         call error_mesg('atmos_ch3i_init','CH3I emission file not specified in field table',FATAL)
      end if
   else
      call error_mesg('atmos_ch3i_init','No emissions specified for CH3I in field table',NOTE)
   end if

!-----------------------------------------------------------------------
!     ... set up tracer concentrations
!-----------------------------------------------------------------------
   if (conc_filename /= '') then
      call interpolator_init( input_conc,trim(conc_filename), lonb_mod, latb_mod,&
                              data_out_of_bounds=(/CONSTANT/), &
                              vert_interp=(/INTERP_WEIGHTED_P/) )
   end if


   call get_tracer_names( MODEL_ATMOS,ind_ch3i,tracer_name,units=tracer_units)
   id_loss =register_diag_field( module_name, 'ch3i_loss', axes(1:3), Time, 'ch3i_loss', TRIM(tracer_units)//'/s' )
   id_j_ch3i =register_diag_field( module_name, 'j_ch3i', axes(1:3), Time, 'j_ch3i', '/s' )

   module_is_initialized = .true.

end subroutine atmos_ch3i_init


!#####################################################################

subroutine atmos_ch3i_time_vary (Time)


type(time_type), intent(in) :: Time

      if (has_emissions) then
        call obtain_interpolator_time_slices (ch3i_emissions, Time)  
      endif
      if (conc_filename /= '') then
        call obtain_interpolator_time_slices (input_conc, Time)  
      endif

end subroutine atmos_ch3i_time_vary


!######################################################################

subroutine atmos_ch3i_endts              


      if (has_emissions) then
        call unset_interpolator_time_flag (ch3i_emissions)  
      endif
      if (conc_filename /= '') then
        call unset_interpolator_time_flag (input_conc)  
      endif

end subroutine atmos_ch3i_endts


!-----------------------------------------------------------------------

subroutine atmos_ch3i( lon, lat, land, pwt, ch3i, ch3i_dt,       &
                       Time, phalf, pfull, t, is, js, dt,    &
                       z_half, z_full, q, tsurf, albedo, coszen, &
                       Time_next, kbot)

real, intent(in),    dimension(:,:)            :: lon, lat
real, intent(in),    dimension(:,:)            :: land
real, intent(in),    dimension(:,:,:)          :: pwt
real, intent(in),    dimension(:,:,:)          :: ch3i
real, intent(out),   dimension(:,:,:)          :: ch3i_dt
type(time_type), intent(in)                    :: Time, Time_next     
integer, intent(in)                            :: is,js
real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
real, intent(in)                               :: dt !to be passed into chemdr
real, intent(in),    dimension(:,:,:)          :: z_half !height in meters at half levels
real, intent(in),    dimension(:,:,:)          :: z_full !height in meters at full levels
real, intent(in),    dimension(:,:,:)          :: q !specific humidity at current time step in kg/kg
real, intent(in),    dimension(:,:)            :: tsurf !surface temperature
real, intent(in),    dimension(:,:)            :: albedo
real, intent(in),    dimension(:,:)            :: coszen
integer, intent(in),  dimension(:,:), optional :: kbot

!-----------------------------------------------------------------------

logical :: used
integer :: i, j, k, id, jd, kd, kb
real, dimension(size(ch3i,1),size(ch3i,2),size(ch3i,3)) :: &
   conc_o3, conc_oh, j_ch3i, k_ch3i_oh, &
   emis_source, ch3i_loss, &
   air_dens
real, dimension(size(ch3i,1),size(ch3i,2)) :: emis,temp_data

!   <ERROR MSG="tracer_driver_init must be called first." STATUS="FATAL">
!     Tracer_driver_init needs to be called before tracer_driver.
!   </ERROR>
if (.not. module_is_initialized)  &
   call error_mesg ('atmos_ch3i','atmos_ch3i_init must be called first1', FATAL)

id=size(ch3i,1); jd=size(ch3i,2); kd=size(ch3i,3)

!-----------------------------------------------------------------------
!     ... Get emissions
!-----------------------------------------------------------------------
emis_source(:,:,:) = 0.
emis(:,:) = 0.
if (has_emissions) then
   do k = 1,size(field_names)
      call interpolator( ch3i_emissions, Time, temp_data, field_names(k), is, js )
      emis(:,:) = emis(:,:) + temp_data(:,:)
   end do
   if (present(kbot)) then
      do j=1,jd
         do i=1,id
            kb=kbot(i,j)
            emis_source(i,j,kb) = emis(i+is,j+js)/pwt(i,j,kb) * emis_cons
         enddo
      enddo
   else
      emis_source(:,:,kd) = emis(:,:)/pwt(:,:,kd) * emis_cons
   end if
   used = send_data(id_emissions,emis,Time_next,is_in=is,js_in=js)
end if

!-----------------------------------------------------------------------
!     ... Get tracer concentrations
!-----------------------------------------------------------------------
call interpolator(input_conc, Time, phalf, conc_o3, 'ox')
call interpolator(input_conc, Time, phalf, conc_oh, 'oh')

air_dens(:,:,:) = 10. * pfull(:,:,:) / (boltz*t(:,:,:)) ! molec/cm3
conc_oh(:,:,:) = conc_oh(:,:,:) * air_dens(:,:,:) ! convert VMR to molec/cm3

!-----------------------------------------------------------------------
!     ... reaction rates
!-----------------------------------------------------------------------
k_ch3i_oh(:,:,:) = 2.9e-12 * exp( -1100./t(:,:,:) ) ! s^-1
do k = 1,kd
   where( coszen(:,:) > 0. )
      j_ch3i(:,:,k) = 2./(4.*SECONDS_PER_DAY) ! s^-1
   elsewhere
      j_ch3i(:,:,k) = 0.
   endwhere
end do
ch3i_loss(:,:,:) = ( j_ch3i(:,:,:) + k_ch3i_oh(:,:,:)*conc_oh(:,:,:) ) &
                 * ch3i(:,:,:) ! VMR/s
used = send_data( id_loss, ch3i_loss, Time_next, is_in=is, js_in=js )
used = send_data( id_j_ch3i, j_ch3i, Time_next, is_in=is, js_in=js )

ch3i_dt(:,:,:) = emis_source(:,:,:) - ch3i_loss(:,:,:)


end subroutine atmos_ch3i

!-----------------------------------------------------------------------

subroutine atmos_ch3i_end

   deallocate( field_names )
   module_is_initialized = .FALSE.

end subroutine atmos_ch3i_end

end module atmos_ch3i_mod
