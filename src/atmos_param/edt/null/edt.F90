module edt_mod

!=======================================================================
!
!
!
!      EDT (Entrainment and Diagnostic Turbulence) MODULE
!
!
!      February 2002
!      Contact person: Steve Klein
!
!
!      These routines calculate the diffusivity coefficients for
!      momentum and temperature-moisture-scalars using the moist
!      thermodynamcs modules based on:
!
!      H. Grenier and C. Bretherton, 2001: A moist PBL parameterization
!      for large-scale models and its application to subtropical
!      cloud-topped marine boundary layers. Mon. Wea. Rev., 129,
!      357-377.
!
!      The actual routine is not described in this paper but is
!      a simplified extension of the parameterization discussed
!      here.  The original code, given to Steve Klein from 
!      Chris Bretherton in May 2001, was tested in the NCAR 
!      atmospheric model, formerly known as CCM. The code has 
!      been adapted for the FMS system by Steve Klein and Paul
!      Kushner.
!
!
!      To quote the Bretherton and Grenier description:
!
!      Driver routine to compute eddy diffusion coefficients for 
!      momentum, moisture, trace constituents and static energy.  Uses 
!      first order closure for stable turbulent layers. For convective 
!      layers, an entrainment closure is used, coupled to a diagnosis 
!      of layer-average TKE from the instantaneous thermodynamic and 
!      velocity profiles. Convective layers are diagnosed by extending 
!      layers of moist static instability into adjacent weakly stably 
!      stratified interfaces, stopping if the stability is too strong.  
!      This allows a realistic depiction of dry convective boundary 
!      layers with a downgradient approach."
! 
!      Authors:  Herve Grenier, 06/2000, Chris Bretherton 09/2000
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! outside modules used
!

use      constants_mod, only: grav,vonkarm,cp_air,rdgas,rvgas,hlv,hls, &
                              tfreeze,radian

use            fms_mod, only: file_exist, open_namelist_file, error_mesg, FATAL,&
                              NOTE, mpp_pe, mpp_root_pe, close_file, &
                              write_version_number, stdlog

use   diag_manager_mod, only: register_diag_field, send_data
        
use   time_manager_mod, only: time_type, get_date, month_name
 
use  monin_obukhov_mod, only: mo_diff

use sat_vapor_pres_mod, only: lookup_es, lookup_des

implicit none
private

!-----------------------------------------------------------------------
!
!      public interfaces

public edt, edt_init, edt_end, edt_on, qaturb, qcturb,tblyrtau

!-----------------------------------------------------------------------
!
! declare version number 
!

character(len=128) :: Version = '$Id: edt.F90,v 17.0 2009/07/21 02:55:02 fms Exp $'
character(len=128) :: Tagname = '$Name: siena_201207 $'
logical            :: module_is_initialized = .false.
logical            :: edt_on = .false.

!-----------------------------------------------------------------------
!
!      global storage variable
!

real, allocatable, dimension(:,:,:) :: qaturb ! cloud fraction diagnosed
                                              ! from turbulence model
      ! (fraction)
real, allocatable, dimension(:,:,:) :: qcturb ! cloud condensate 
                                              ! diagnosed from turb.
      ! model (kg liq/kg air)
real, allocatable, dimension(:,:,:) :: tblyrtau  ! turbulent layer
                                                 ! time scale


contains

!======================================================================= 
!
!      subroutine edt_init 
!        
!
!      this subroutine reads the namelist file and restart data
!      and initializes some constants.
!        

subroutine edt_init(lonb, latb, axes,time,idim,jdim,kdim)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
! 
!      idim,jdim,kdim    size of the first 3 dimensions 
!      axes, time        variables needed for netcdf diagnostics
!      latb, lonb        latitudes and longitudes at grid box boundaries
!
!
!      --------
!      internal
!      --------
! 
!      unit              unit number for namelist and restart file
!      io                internal variable for reading of namelist file
!      full              indices for full level axes coordinates
!      half              indices for half level axes coordinates
!
!-----------------------------------------------------------------------

integer,         intent(in) :: idim,jdim,kdim,axes(4)
type(time_type), intent(in) :: time
real, dimension(:,:),intent(in) :: lonb, latb


!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif


!-----------------------------------------------------------------------
!
!      initialize edt_on

       module_is_initialized = .true.
end subroutine edt_init

!
!======================================================================= 




!======================================================================= 
!
!      subroutine edt
!        
!
!       this subroutine is the main driver program to the routines
!       provided by Chris Bretherton
!        

subroutine edt(is,ie,js,je,dt,time,tdtlw_in, u_star,b_star,q_star,t,qv,ql,qi,qa, &
               u,v,z_full,p_full,z_half,p_half,stbltop,k_m,k_t,pblh,  &
               kbot,tke)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      is,ie,js,je  i,j indices marking the slab of model working on
!      dt        physics time step (seconds)
!      time      variable needed for netcdf diagnostics
!      u_star    friction velocity (m/s)
!      b_star    buoyancy scale (m/(s**2))
!      q_star    moisture scale (kg vapor/kg air)
!
!      three dimensional fields on model full levels, reals dimensioned
!      (:,:,pressure), third index running from top of atmosphere to 
!      bottom
!          
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      ql        liquid water specific humidity (kg cond/kg air)
!      qi        ice water specific humidity (kg cond/kg air)
!      qa        cloud fraction 
!      u         zonal wind (m/s)
!      v         meridional wind (m/s) 
!      z_full    height of full levels (m)
!      p_full    pressure (Pa)
!
!      the following two fields are on the model half levels, with
!      size(z_half,3) = size(t,3) +1, z_half(:,:,size(z_half,3)) 
!      must be height of surface (if you are not using eta-model)
!
!      z_half    height at half levels (m)
!      p_half    pressure at half levels (Pa)
!        
!      ------
!      output
!      ------
!
!      stbltop   maximum altitude the very stable boundary layer
!                is permitted to operate
!
!      The following variables are defined at half levels and are
!      dimensions 1:nlev+1.
!
!      k_m       diffusivity for momentum (m**2/s)
!      k_t       diffusivity for temperature and scalars (m**2/s)
!
!      k_m and k_t are defined at half-levels so that size(k_m,3) 
!      should be at least as large as size(t,3). Note, however, that 
!      only the returned values at levels 2 to size(t,3) are 
!      meaningful; other values will be returned as zero.
!
!      --------------
!      optional input
!      --------------
!
!      kbot      integer indicating the lowest true layer of atmosphere
!
!      ---------------
!      optional output
!      ---------------
!
!      pblh      depth of planetary boundary layer (m)
!      tke       turbulent kinetic energy (m*m)/(s*s)
!
!-----------------------------------------------------------------------

integer,         intent(in)                            :: is,ie,js,je
real,            intent(in)                            :: dt
type(time_type), intent(in)                            :: time
real,            intent(in),  dimension(:,:,:)         :: tdtlw_in
real,            intent(in),  dimension(:,:)           :: u_star,b_star
real,            intent(in),  dimension(:,:)           :: q_star
real,            intent(in),  dimension(:,:,:)         :: t,qv,ql,qi,qa
real,            intent(in),  dimension(:,:,:)         :: u, v
real,            intent(in),  dimension(:,:,:)         :: z_full, p_full
real,            intent(in),  dimension(:,:,:)         :: z_half, p_half
real,            intent(out), dimension(:,:)           :: stbltop
real,            intent(out), dimension(:,:,:)         :: k_m,k_t
integer,         intent(in),  dimension(:,:), optional :: kbot
!real,            intent(out), dimension(:,:), optional :: pblh
real,            intent(out), dimension(:,:)           :: pblh
real,            intent(out), dimension(:,:,:),optional:: tke

end subroutine edt

!
!======================================================================= 

!======================================================================= 
!
!      subroutine edt_end
!        
!
!      this subroutine writes out the restart field
!        

subroutine edt_end()

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
!      subroutine end
!
       module_is_initialized = .false.

end subroutine edt_end

!
!=======================================================================

end module edt_mod
