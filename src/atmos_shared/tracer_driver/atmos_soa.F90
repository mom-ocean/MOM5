module atmos_soa_mod
! <DESCRIPTION>
!   This module is an implementation of Secondary organic aerosols (SOA)
!   from anthropogenic activities, and is based on Tie et al. (JGR, 2003).
!   The only souce of SOA is due to the oxydation of C4H10 by OH.
!   The concentrations of these 2 gas species are read as input.
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory. 
!  This implies that the "atmos_SOA_init" should be executed at the begining 
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use mpp_mod, only: input_nml_file 
use                    fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL
use           time_manager_mod, only : time_type
use           diag_manager_mod, only : send_data,               &
                                       register_diag_field,     &
                                       register_static_field
use         tracer_manager_mod, only : get_tracer_index,        &
                                       set_tracer_atts
use          field_manager_mod, only : MODEL_ATMOS
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR
use           interpolator_mod, only:  interpolate_type,  &
                                       interpolator_init, &
                                       obtain_interpolator_time_slices,&
                                       unset_interpolator_time_flag, &
                                       interpolator, interpolator_end, &
                                       CONSTANT, INTERP_WEIGHTED_P

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_SOA_init, atmos_SOA_end, atmos_SOA_chem, &
        atmos_SOA_time_vary, atmos_soa_endts

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSOA = 0  ! tracer number for Secondary Organic Aerosol 

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH_conc            = 0
integer ::   id_C4H10_conc         = 0
integer ::   id_SOA_chem           = 0
integer ::   id_SOA_chem_col       = 0

type(interpolate_type),save         ::  gas_conc_interp
character(len=32)  :: gas_conc_filename = 'gas_conc_3D.nc'
character(len=32), dimension(2) :: gas_conc_name
data gas_conc_name/'OH','C4H10'/

namelist /secondary_organics_nml/ gas_conc_filename, gas_conc_name

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_soa.F90,v 20.0 2013/12/13 23:24:02 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_SOA_init">
!<OVERVIEW>
! The constructor routine for the soa module.
!</OVERVIEW>
 subroutine atmos_SOA_init ( lonb, latb, nlev, axes, Time, mask)
!-----------------------------------------------------------------------
real,             intent(in), dimension(:,:)        :: lonb, latb
integer,          intent(in)                        :: nlev
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
!
!-----------------------------------------------------------------------
!
      integer  unit,io,ierr, logunit
      character(len=3) :: SOA_tracer
!
      data SOA_tracer/'SOA'/

!
      if (module_is_initialized) return
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=secondary_organics_nml, iostat=io)
        ierr = check_nml_error(io,'secondary_organics_nml')
#else
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=secondary_organics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'secondary_organics_nml')
        end do
10      call close_file (unit)
#endif
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=secondary_organics_nml)

!----- set initial value of soa ------------

      nSOA = get_tracer_index(MODEL_ATMOS,'SOA')
      if (nSOA > 0) then
         call set_tracer_atts(MODEL_ATMOS,'SOA','SOA','mmr')
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (*,30) SOA_tracer,nsoa
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (logunit,30) SOA_tracer,nsoa
      endif


  30   format (A,' was initialized as tracer number ',i2)

     call interpolator_init (gas_conc_interp, trim(gas_conc_filename),  &
                             lonb, latb,&        
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gas_conc_name, & 
                             vert_interp=(/INTERP_WEIGHTED_P/) )
      if (id_OH_conc .eq. 0 ) &
        id_OH_conc    = register_diag_field ( mod_name,           &
                      'OH_SOA_conc',axes(1:3),Time,                        &
                      'Hydroxyl radical concentration',           &
                      'molec.cm-3')

      id_C4H10_conc    = register_diag_field ( mod_name,           &
                      'C4H10_mmr',axes(1:3),Time,                        &
                      'nButane concentration',           &
                      'mmr')

      id_SOA_chem    = register_diag_field ( mod_name,       &
                      'SOA_chem',axes(1:3),Time,            &
                      'SOA production by C4H10 + OH',        &
                      'kg/m2/s')

      id_SOA_chem_col= register_diag_field ( mod_name,       &
                      'SOA_chem_col',axes(1:2),Time,            &
                      'column SOA production by C4H10 + OH',        &
                      'kg/m2/s')

      call write_version_number (version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
 end subroutine atmos_SOA_init




!#####################################################################

subroutine atmos_SOA_time_vary (Time)

type(time_type), intent(in) :: Time


      call obtain_interpolator_time_slices (gas_conc_interp, Time)

end subroutine atmos_SOA_time_vary


!#####################################################################

subroutine atmos_SOA_endts             


      call unset_interpolator_time_flag (gas_conc_interp)


end subroutine atmos_SOA_endts



!#####################################################################

!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_SOA_end">
!<OVERVIEW>
!  The destructor routine for the soa module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_SOA_end
!</TEMPLATE>
 subroutine atmos_SOA_end

      call interpolator_end (gas_conc_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_SOA_end
!</SUBROUTINE>
!-----------------------------------------------------------------------
      SUBROUTINE atmos_SOA_chem(pwt,temp,pfull, phalf, dt, &
                          jday,hour,minute,second,lat,lon, &
                          SOA, SOA_dt, Time,Time_next,is,ie,js,je,kbot)

! ****************************************************************************
      real, intent(in),    dimension(:,:,:)          :: pwt
      real, intent(in),    dimension(:,:,:)          :: temp,pfull,phalf
      real, intent(in)                               :: dt
      integer, intent(in)                            :: jday, hour,minute,second
      real, intent(in),  dimension(:,:)              :: lat, lon  ! [radian]
      real, intent(in),    dimension(:,:,:)          :: SOA
      real, intent(out),   dimension(:,:,:)          :: SOA_dt
      type(time_type), intent(in)                    :: Time, Time_next
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      real, dimension(size(SOA,1),size(SOA,2),size(SOA,3)) :: &
               SOA_chem, OH_conc, C4H10_conc
      real, dimension(size(SOA,1),size(SOA,2)) :: &
               SOA_prod, &
               xu, dayl, h, hl, hc, hred, fac_OH, fact_OH
      real, parameter                            :: wtm_C = 12.
      real, parameter                            :: wtm_C4H10 = 58.
      real, parameter                            :: yield = 0.1
      real, parameter                            :: small_value=1.e-21
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real                                       :: decl, hd, x
      integer :: i,j,k,id,jd,kd
      integer                                    :: istep, nstep
! Local grid sizes
      id=size(SOA,1); jd=size(SOA,2); kd=size(SOA,3)

      OH_conc(:,:,:)=0.  ! molec/cm3
      call interpolator(gas_conc_interp, Time, phalf, OH_conc, &
                       trim(gas_conc_name(1)), is, js)

      C4H10_conc(:,:,:)=0.0
      call interpolator(gas_conc_interp, Time, phalf, C4H10_conc, &
                       trim(gas_conc_name(2)), is, js)
      C4H10_conc(:,:,:)=C4H10_conc(:,:,:)*WTM_C4H10/WTMAIR

      x = 2. *pi *float(jday-1)/365.
      decl = A0 - A1*cos(  X) + B1*sin(  X) - A2*cos(2.*X) + B2*sin(2.*X) &
           - A3*cos(3.*X) + B3*sin(3.*X)
      xu(:,:) = -tan(lat(:,:))*tan(decl)
      where ( xu > -1 .and. xu < 1 ) dayl=acos(xu)/pi
      where ( xu <= -1 ) dayl = 1.
      where ( xu >= 1 ) dayl = 0.
!   Calculate normalization factors for OH and NO3 such that
!   the diurnal average respect the monthly input values.
      hd=0.
      fact_OH(:,:)  = 0.
      nstep = int(24.*3600./dt)
      do istep=1,nstep
        hd=hd+dt/3600./24.
        hl(:,:) = pi*(1.-dayl(:,:))
        hc(:,:) = pi*(1.+dayl(:,:))
        h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
        where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fact_OH  = fact_OH + amax1(0.,sin(pi*hred)/2.)/nstep
        endwhere
      enddo


      hd=amax1(0.,amin1(1.,(hour+minute/60.+second/3600.)/24.))
      hl(:,:) = pi*(1.-dayl(:,:))
      hc(:,:) = pi*(1.+dayl(:,:))
      h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
      fac_OH(:,:)  = 0.
      where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fac_OH  = amax1(0.,sin(pi*hred)/2.)/fact_OH
      elsewhere
! Nightime
          fac_OH  = 0.
      endwhere

      do i=1,id
        do j=1,jd
          do k=1,kd 
            SOA_dt(i,j,k) = 1.55E-11 * exp( -540./temp(i,j,k) ) *yield &
                * C4H10_conc(i,j,k)*OH_conc(i,j,k)*fac_oh(i,j)
          enddo
        enddo
      enddo

      SOA_chem(:,:,:)=SOA_dt(:,:,:)*pwt(:,:,:)

      if (id_SOA_chem > 0) then
        used = send_data ( id_SOA_chem, &
              SOA_chem, Time_next,is_in=is,js_in=js,ks_in=1)
      endif

! column production of SOA 


      SOA_prod = 0.
      do k=1,kd
        SOA_prod = SOA_prod +  SOA_chem(:,:,k)
      end do

      if (id_SOA_chem_col > 0) then
        used = send_data ( id_SOA_chem_col, &
                           SOA_prod, Time_next,is_in=is,js_in=js)
      endif


end subroutine atmos_SOA_chem


end module atmos_SOA_mod
