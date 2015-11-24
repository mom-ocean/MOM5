!FDOC_TAG_GFDL

module cu_mo_trans_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Isaac Held
! </CONTACT>
! <OVERVIEW>
!    A simple module that computes a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
! </OVERVIEW>
! <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in improving the simulation of tropical
!     precipitation in some models.  The diffusivity computed here is
!     simply 
!<PRE>
! diffusivity = c*W*L 
! W = M/rho  (m/sec) 
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air <p>
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g 
!   (diff_norm is a namelist parameter;
!      the factor of g = acceleration of gravity here is an historical artifact) <p>
! for further discussion see 
!     <LINK SRC="cu_mo_trans.pdf">cu_mo_trans.pdf</LINK>
!</PRE>
! </DESCRIPTION>


!=======================================================================
!
!                 DIFFUSIVE CONVECTIVE MOMENTUM TRANSPORT MODULE
!
!=======================================================================

  use   constants_mod, only:  GRAV, RDGAS, RVGAS, CP_AIR
 

  use         mpp_mod, only: input_nml_file
  use         fms_mod, only: file_exist, check_nml_error,    &
                             open_namelist_file, close_file, &
                             write_version_number,           &
                             mpp_pe, mpp_root_pe, stdlog,    &
                             error_mesg, FATAL, NOTE

  use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use  Time_Manager_Mod, ONLY: time_type

implicit none
private


! public interfaces
!=======================================================================
public :: cu_mo_trans_init, &
          cu_mo_trans,      &
          cu_mo_trans_end

!=======================================================================

! form of interfaces
!=======================================================================

      
logical :: module_is_initialized = .false.


!---------------diagnostics fields------------------------------------- 

integer :: id_diff_cmt, id_utnd_cmt, id_vtnd_cmt, id_ttnd_cmt, &
           id_massflux_cmt, id_detmf_cmt

character(len=11) :: mod_name = 'cu_mo_trans'

real :: missing_value = -999.
logical  ::  do_diffusive_transport = .false.
logical  ::  do_nonlocal_transport = .false.


!--------------------- namelist variables with defaults -------------

real    :: diff_norm =   1.0
logical :: limit_mass_flux = .false.  ! when true, the mass flux 
                                      ! out of a grid box is limited to
                                      ! the mass in that grid box
character(len=64) :: transport_scheme = 'diffusive'
integer :: non_local_iter = 2  ! iteration count for non-local scheme
logical :: conserve_te = .true.  ! conserve total energy ?
real    ::  gki = 0.7  ! Gregory et. al. constant for p-gradient param
real    :: amplitude = 1.0 ! Tuning parameter (1=full strength)

namelist/cu_mo_trans_nml/ diff_norm, &
                          limit_mass_flux, &
                          non_local_iter, conserve_te, gki, &
                          amplitude,  &
                          transport_scheme


!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: cu_mo_trans.F90,v 19.0 2012/01/06 20:05:02 fms Exp $'
character(len=128) :: tagname = '$Name: tikal $'

contains

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_init">
!  <OVERVIEW>
!   initializes module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Reads namelist and registers one diagnostic field
!     (diff_cmt:  the kinematic diffusion coefficient)
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_init( axes, Time )
!
!  </TEMPLATE>
!  <IN NAME=" axes" TYPE="integer">
!    axes identifier needed by diag manager
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    time at initialization needed by diag manager
!  </IN>
! </SUBROUTINE>
!
subroutine cu_mo_trans_init( axes, Time, doing_diffusive )

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 logical,         intent(out) :: doing_diffusive

integer :: unit, ierr, io, logunit
integer, dimension(3)  :: half =  (/1,2,4/)

!------ read namelist ------

   if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=cu_mo_trans_nml, iostat=io)
      ierr = check_nml_error(io,'cu_mo_trans_nml')
#else   
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=cu_mo_trans_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'cu_mo_trans_nml')
      enddo
 10   call close_file (unit)
#endif
   endif

!--------- write version number and namelist ------------------

      call write_version_number(version, tagname)
      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) &
        write ( logunit, nml=cu_mo_trans_nml )

!----------------------------------------------------------------------
!    define logicals indicating momentum transport scheme to use.
!----------------------------------------------------------------------
      if (trim(transport_scheme) == 'diffusive') then
        do_diffusive_transport = .true.
      else if (trim(transport_scheme) == 'nonlocal') then
        do_nonlocal_transport = .true.
      else
        call error_mesg ('cu_mo_trans', &
         'invalid specification of transport_scheme', FATAL)
      endif

      doing_diffusive = do_diffusive_transport

! --- initialize quantities for diagnostics output -------------

   if (do_diffusive_transport) then
     id_diff_cmt = &
      register_diag_field ( mod_name, 'diff_cmt', axes(1:3), Time,    &
                        'cu_mo_trans coeff for momentum',  'm2/s', &
                         missing_value=missing_value               )
     id_massflux_cmt = &
      register_diag_field ( mod_name, 'massflux_cmt', axes(half), Time, &
                        'cu_mo_trans mass flux',  'kg/(m2 s)', &
                         missing_value=missing_value               )
    else if (do_nonlocal_transport) then 
     id_utnd_cmt = &
      register_diag_field ( mod_name, 'utnd_cmt', axes(1:3), Time,    &
                        'cu_mo_trans u tendency',  'm/s2', &
                         missing_value=missing_value               )
     id_vtnd_cmt = &
      register_diag_field ( mod_name, 'vtnd_cmt', axes(1:3), Time,    &
                        'cu_mo_trans v tendency',  'm/s2', &
                         missing_value=missing_value               )
     id_ttnd_cmt = &
      register_diag_field ( mod_name, 'ttnd_cmt', axes(1:3), Time,    &
                        'cu_mo_trans temp tendency',  'deg K/s', &
                         missing_value=missing_value               )
     id_massflux_cmt = &
      register_diag_field ( mod_name, 'massflux_cmt', axes(half), Time, &
                        'cu_mo_trans mass flux',  'kg/(m2 s)', &
                         missing_value=missing_value               )
     id_detmf_cmt = &
      register_diag_field ( mod_name, 'detmf_cmt', axes(1:3), Time,  &
                      'cu_mo_trans detrainment mass flux',  'kg/(m2 s)',&
                         missing_value=missing_value               )
    endif

!--------------------------------------------------------------

  module_is_initialized = .true.


end subroutine cu_mo_trans_init

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans_end">
!  <OVERVIEW>
!   terminates module
!  </OVERVIEW>
!  <DESCRIPTION>
!   This is the destructor for cu_mo_trans
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cu_mo_trans_end()

  module_is_initialized = .false.

end subroutine cu_mo_trans_end

!#######################################################################

! <SUBROUTINE NAME="cu_mo_trans">
!  <OVERVIEW>
!   picks one of the available cumulus momentum transport parameteriz-
!    ations based on namelist-supplied information (currently diffusive 
!    or non-local options are available).  For the diffusive scheme, it
!    returns a diffusivity proportional to the convective mass 
!    flux. For the non-local scheme temperature, specific humidity and
!    momentum tendencies are returned.
!  </OVERVIEW>
!  <DESCRIPTION>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cu_mo_trans (is, js, Time, mass_flux, t,           &
!                p_half, p_full, z_half, z_full, diff)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! ! horizontal domain on which computation to be performed is
!    (is:is+size(t,1)-1,ie+size(t,2)-1) in global coordinates
!   (used by diag_manager only)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! current time, used by diag_manager
!  </IN>
!  <IN NAME="mass_flux" TYPE="real">
! convective mass flux (Kg/(m**2 s)), dimension(:,:,:), 3rd dimension is
!    vertical level (top down) -- defined at interfaces, so that
!    size(mass_flux,3) = size(p_half,3); entire field processed;
!   all remaining fields are 3 dimensional;
!  size of first two dimensions must confrom for all variables
!  </IN>
!  <IN NAME="t" TYPE="real">
! temperature (K) at full levels, size(t,3) = size(p_full,3)
!  </IN>
!  <IN NAME="p_half" TYPE="real">
! pressure at interfaces (Pascals) 
!  size(p_half,3) = size(p_full,3) + 1
!  </IN>
!  <IN NAME="p_full" TYPE="real">
! pressure at full levels (levels at which temperature is defined)
!  </IN>
!  <IN NAME="z_half" TYPE="real">
! height at half levels (meters); size(z_half,3) = size(p_half,3)
!  </IN>
!  <IN NAME="z_full" TYPE="real">
! height at full levels (meters); size(z_full,3) = size(p_full,3)
!  </IN>
!  <OUT NAME="diff" TYPE="real">
! kinematic diffusivity (m*2/s); defined at half levels 
!   size(diff,3) = size(p_half,3)
!  </OUT>
! </SUBROUTINE>
!




subroutine cu_mo_trans (is, js, Time, mass_flux, t,           &
                        p_half, p_full, z_half, z_full, dt, uin, vin,&
                        tracer, pmass, det0, utnd, vtnd, ttnd, &
                        qtrcumo, diff)

type(time_type), intent(in) :: Time
integer,         intent(in) :: is, js

real,   intent(in) :: dt
real, intent(inout)   , dimension(:,:,:,:) :: tracer           
real, intent(inout)   , dimension(:,:,:) :: uin, vin, t        
real, intent(in)   , dimension(:,:,:) :: mass_flux, &
                                         pmass, det0, &
                                         p_half, p_full, z_half, z_full
real, intent(out), dimension(:,:,:) :: utnd, vtnd          
real, intent(out), dimension(:,:,:) :: ttnd          
real, intent(out), dimension(:,:,:,:) :: qtrcumo                   
real, intent(inout), dimension(:,:,:) :: diff                      


      integer :: im, jm, km, nq, nq_skip

!-----------------------------------------------------------------------

 if (.not.module_is_initialized) call error_mesg ('cu_mo_trans',  &
                      'cu_mo_trans_init has not been called.', FATAL)

!-----------------------------------------------------------------------
!   utnd = 0.
!   vtnd = 0.
!   ttnd = 0.
!   qtrcumo = 0.

    if (do_diffusive_transport) then
      call diffusive_cu_mo_trans (is, js, Time, mass_flux, t,      &
                      p_half, p_full, z_half, z_full, diff)
    utnd = 0.
    vtnd = 0.
    ttnd = 0.
    qtrcumo = 0.
      
    else if (do_nonlocal_transport) then
!      call error_mesg ('cu_mo_trans',  &
!        'non-local transport not currently available', FATAL)
       im = size(uin,1)
       jm = size(uin,2)
       km = size(uin,3)
       nq = size(tracer,4)
       nq_skip = nq
       qtrcumo(:,:,:,1:nq_skip) = 0.0
       call non_local_mot (im, jm, km, is, js, Time, dt, t, uin, vin,  &
                           nq, nq_skip,   &
                 tracer, pmass, mass_flux, det0, utnd, vtnd, ttnd,  &
                 qtrcumo)
    endif

end subroutine cu_mo_trans


!#######################################################################

! <SUBROUTINE NAME="diffusive_cu_mo_trans">
!  <OVERVIEW>
!   returns a diffusivity proportional to the 
!    convective mass flux, for use with diffusive 
!    convective momentum transport closure
!  </OVERVIEW>
!  <DESCRIPTION>
!   A diffusive approximation to convective momentum transport is crude but
!    has been found to be useful in inproving the simulation of tropical
!    precipitation in some models.  The diffusivity computed here is
!    simply 
!<PRE>
! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! rho - density of air
! L = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (diff_norm is a namelist parameter;
!      the factor of g here is an historical artifact)
! for further discussion see cu_mo_trans.ps
!</PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diffusive_cu_mo_trans (is, js, Time, mass_flux, t,      &
!                p_half, p_full, z_half, z_full, diff)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! ! horizontal domain on which computation to be performed is
!    (is:is+size(t,1)-1,ie+size(t,2)-1) in global coordinates
!   (used by diag_manager only)
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! current time, used by diag_manager
!  </IN>
!  <IN NAME="mass_flux" TYPE="real">
! convective mass flux (Kg/(m**2 s)), dimension(:,:,:), 3rd dimension is
!    vertical level (top down) -- defined at interfaces, so that
!    size(mass_flux,3) = size(p_half,3); entire field processed;
!   all remaining fields are 3 dimensional;
!  size of first two dimensions must confrom for all variables
!  </IN>
!  <IN NAME="t" TYPE="real">
! temperature (K) at full levels, size(t,3) = size(p_full,3)
!  </IN>
!  <IN NAME="p_half" TYPE="real">
! pressure at interfaces (Pascals) 
!  size(p_half,3) = size(p_full,3) + 1
!  </IN>
!  <IN NAME="p_full" TYPE="real">
! pressure at full levels (levels at which temperature is defined)
!  </IN>
!  <IN NAME="z_half" TYPE="real">
! height at half levels (meters); size(z_half,3) = size(p_half,3)
!  </IN>
!  <IN NAME="z_full" TYPE="real">
! height at full levels (meters); size(z_full,3) = size(p_full,3)
!  </IN>
!  <OUT NAME="diff" TYPE="real">
! kinematic diffusivity (m*2/s); defined at half levels 
!   size(diff,3) = size(p_half,3)
!  </OUT>
! </SUBROUTINE>
!
subroutine diffusive_cu_mo_trans (is, js, Time, mass_flux, t,     &
                        p_half, p_full, z_half, z_full, diff)

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: is, js

real, intent(in)   , dimension(:,:,:) :: mass_flux, t, &
                                         p_half, p_full, z_half, z_full
real, intent(out), dimension(:,:,:) :: diff                      

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho
real, dimension(size(t,1),size(t,2))           :: zbot, ztop

integer :: k, nlev
logical :: used

!-----------------------------------------------------------------------

 if (.not.module_is_initialized) call error_mesg ('cu_mo_trans',  &
                      'cu_mo_trans_init has not been called.', FATAL)

!-----------------------------------------------------------------------

nlev = size(t,3)

zbot = z_half(:,:,nlev+1)
ztop = z_half(:,:,nlev+1)
  
do k = 2, nlev
  where(mass_flux(:,:,k) .ne. 0.0 .and. mass_flux(:,:,k+1) == 0.0) 
    zbot = z_half(:,:,k)
  endwhere
  where(mass_flux(:,:,k-1) == 0.0 .and. mass_flux(:,:,k) .ne. 0.0) 
    ztop = z_half(:,:,k)
  endwhere
end do

rho  = p_full/(RDGAS*t)  ! density 
   ! (including the virtual temperature effect here might give the 
   ! impression that this theory is accurate to 2%!)

! diffusivity = c*W*L
! W = M/rho  (m/sec)
! M = convective mass flux (kg/(m2 sec)) 
! L = ztop - zbot = depth of convecting layer (m)
! c = normalization constant = diff_norm/g
!   (the factor of g here is an historical artifact)

diff(:,:,1) = 0.0
do k = 2, nlev
  diff(:,:,k) = diff_norm*mass_flux(:,:,k)*(ztop-zbot)/(rho(:,:,k)*GRAV)
end do


! --- diagnostics
     if ( id_diff_cmt > 0 ) then
        used = send_data ( id_diff_cmt, diff, Time, is, js, 1 )
     endif
     if ( id_massflux_cmt > 0 ) then
        used = send_data ( id_massflux_cmt, mass_flux, Time, is, js, 1 )
     endif

end subroutine diffusive_cu_mo_trans


!#######################################################################

 subroutine non_local_mot(im, jm, km, is, js, Time, dt, tin, uin, vin, nq, nq_skip, qin, pmass, mc,       &
                          det0, utnd, vtnd, ttnd, qtnd)
!
! This is a non-local cumulus transport algorithm based on the given cloud mass fluxes (mc).
! Detrainment fluxes are computed internally by mass (or momentum) conservation.
! Simple upwind algorithm is used. It is possible replace the large-scale part by
! a higher-order upwind scheme (such as van Leer or PPM). But the cost is perhaps
! not worth it becasue there is so much uncertainty in the accuracy of the cloud
! mass fluxes.
! Contact Shian-Jiann Lin for more information and a tech-note.
!
   implicit none
  integer, intent(in):: im, jm, km                   ! dimensions
  integer, intent(in):: is, js                                   
  type(time_type), intent(in) :: Time      ! timestamp for diagnmostics
  integer, intent(in):: nq                           ! tracer dimension
  integer, intent(in):: nq_skip                      ! # of tracers to skip
                                ! Here the assumption is that sphum is the #1 tracer
!  integer, intent(in):: iter                         ! Number of iterations
  real,    intent(in):: dt                           ! model time step (seconds)
  real,    intent(inout):: tin (im,jm,km)  ! temperature    
  real,    intent(inout):: uin(im,jm,km), vin(im,jm,km) ! input winds (m/s)
  real,    intent(inout):: qin(im,jm,km,nq)
  real,    intent(in):: pmass(im,jm,km)              ! layer mass (kg/m**2)
  real,    intent(in):: mc(im,jm,km+1)               ! cloud mass flux [kg/(s m**2)]
                                                     ! positive -> upward
  real,    intent(in):: det0(im,jm,km)               ! detrained mass fluxes

!  real,    intent(in):: cp                           ! 
!  real,    intent(in):: gki                          ! Gregory et al constant (0.7)
!  logical, intent(in):: conserve                     ! Conserve Total Energy?
! Output
! added on tendencies
! real,    intent(inout)::ttnd(im,jm,km)             ! temperature due to TE conservation
  real,    intent(out)::ttnd(im,jm,km)             ! temperature due to TE conservation
! real,    intent(inout)::qtnd(im,jm,km,nq)
  real,    intent(out)::qtnd(im,jm,km,nq)

  real,    intent(out)::utnd(im,jm,km), vtnd(im,jm,km)   ! m/s**2
!
! Local 
  real dm1(km), u1(km), v1(km), u2(km), v2(km)
  real q1(km,nq), q2(km,nq), qc(km,nq)
  real uc(km), vc(km)
  real mc1(km+1)
  real mc2(km+1)
  real det1(km), det2(km)

  integer i, j, k, it, iq
  integer ktop, kbot, kdet
  real rdt, ent, fac_mo, fac_t
  real,  parameter:: eps = 1.E-15       ! Cutoff value
! real,  parameter:: amplitude = 1.0    ! Tuning parameter (1=full strength)
  real x_frac
  real rtmp
  logical :: used

  utnd = 0.
  vtnd = 0.
  ttnd =0.
  qtnd = 0.

  rdt   = 1./dt
  fac_mo = amplitude * dt
!  fac_t = 0.5/(dt*cp)
  fac_t = 0.5/(dt*CP_AIR)

!  write(*,*) 'Within non_local_mot, num_tracers=', nq

  do j=1,jm
     do i=1,im

! Copy to 1D arrays to better utilize the cache
        do k=1,km
           dm1(k) = pmass(i,j,k)
        enddo

        do k=1,km
          if (limit_mass_flux) then
 ! Limit mass flux by available layer mass
            mc1(k) = min(dm1(k),   mc(i,j,k)*fac_mo) 
           det1(k) = min(dm1(k), det0(i,j,k)*fac_mo)  
          else
            mc1(k) =  mc(i,j,k)*fac_mo
           det1(k) =  det0(i,j,k)*fac_mo  
          endif
        enddo

!---------------------------
! Locate cloud base
!---------------------------
        kbot = 1
        do k = km,2,-1
           if( mc1(k) > eps ) then
                 kbot = k
                 go to 1111
           endif
        end do
1111    continue

!---------------------------
! Locate cloud top
!---------------------------
        ktop = km
        do k = 1, km-1 
           if( mc1(k+1) > eps ) then
                 ktop = k
                 go to 2222
           endif
        end do
2222    continue

        if ( kbot > ktop+1 ) then    ! ensure there are at least 3 layers to work with

           kdet = ktop
           do k=kbot-1,ktop,-1
              if ( det1(k) > eps ) then
                 kdet = k
                 go to 3333
              endif
           enddo
3333       continue
                                     ! cloudbase, interior, and cloudtop layers
           do k=ktop, kbot
              u1(k) = uin(i,j,k)
              u2(k) = u1(k)
              v1(k) = vin(i,j,k)
              v2(k) = v1(k)
           enddo

        if ( nq_skip < nq ) then
           do iq=nq_skip+1,nq
           do k=ktop, kbot
              q1(k,iq) = qin(i,j,k,iq)
              q2(k,iq) = q1(k,iq)
           enddo
           enddo
        endif

!       do it=1, iter
        do it=1, non_local_iter

!          x_frac = real(it) / real(iter)
           x_frac = real(it) / real(non_local_iter)
           do k=ktop, kbot
                 mc2(k) = x_frac * mc1(k)
                det2(k) = x_frac * det1(k)
           enddo

!----------------------------------------------------------
! In-cloud fields: Cloud base
!----------------------------------------------------------
           rtmp = 1. / (dm1(kbot)+mc2(kbot))
           uc(kbot) = (dm1(kbot)*u1(kbot) + mc2(kbot)*u2(kbot-1)) * rtmp
           vc(kbot) = (dm1(kbot)*v1(kbot) + mc2(kbot)*v2(kbot-1)) * rtmp

           if ( nq_skip < nq ) then
              do iq=nq_skip+1,nq
                 qc(kbot,iq) = (dm1(kbot)*q1(kbot,iq) + mc2(kbot)*q2(kbot-1,iq)) * rtmp
              enddo
           endif
!----------------------------------------------------------
! In-cloud fields: interior
!----------------------------------------------------------
!
!
! Below the detrainment level:  (det=0)
           if ( kdet < kbot-1 ) then
!-----------------------------------------------------------
! The in-cloud fields are modified (diluted) by entrainment of
! environment air, and the GKI effect will not be added here.
!-----------------------------------------------------------
              do k=kbot-1, kdet+1,-1
                 ent = mc2(k) - mc2(k+1) + det2(k)
                uc(k) = (mc2(k+1)*uc(k+1)+ent*u2(k)) / mc2(k)  
                vc(k) = (mc2(k+1)*vc(k+1)+ent*v2(k)) / mc2(k) 
              enddo
              if ( nq_skip < nq ) then
                 do iq=nq_skip+1,nq
                    do k=kbot-1, kdet+1,-1
                       ent = mc2(k) - mc2(k+1) + det2(k)
                       qc(k,iq) = (mc2(k+1)*qc(k+1,iq)+ent*q2(k,iq)) / mc2(k) 
                    enddo
                 enddo
              endif
           endif
!
! Pressure-gradient effect of Gregory et al 1997 added when
! the clouds are detraining.
! Entrained mass fluxes are diagnosed by mass conservation law
           if ( kdet > ktop ) then
              do k=kdet, ktop+1, -1
                 ent = mc2(k) - mc2(k+1) + det2(k)
                 uc(k) = (mc2(k+1)*uc(k+1)+ent*u2(k))/(mc2(k+1)+ent)     &
                        + gki*(u2(k)-u2(k+1))
                 vc(k) = (mc2(k+1)*vc(k+1)+ent*v2(k))/(mc2(k+1)+ent)     &
                        + gki*(v2(k)-v2(k+1))
              enddo
              if ( nq_skip < nq ) then
                 do iq=nq_skip+1,nq
                    do k=kdet, ktop+1, -1
                       ent = mc2(k) - mc2(k+1) + det2(k)
                       qc(k,iq) = (mc2(k+1)*qc(k+1,iq)+ent*q2(k,iq))/(mc2(k+1)+ent)
                    enddo
                 enddo
              endif
           endif

!----------------
! Update fields:
!----------------
! Cloud top
          rtmp = 1. / (dm1(ktop)+mc2(ktop+1))
          u2(ktop) = (dm1(ktop)*u1(ktop)+mc2(ktop+1)*uc(ktop+1)) * rtmp
          v2(ktop) = (dm1(ktop)*v1(ktop)+mc2(ktop+1)*vc(ktop+1)) * rtmp
          if ( nq_skip < nq ) then
               do iq=nq_skip+1,nq
                  q2(ktop,iq) = (dm1(ktop)*q1(ktop,iq)+mc2(ktop+1)*qc(ktop+1,iq)) * rtmp
               enddo
          endif

!---------------------------------------------------------------------
! Interior (this loop can't be vectorized due to k to k-1 dependency)
!---------------------------------------------------------------------
           do k=ktop+1,kbot-1
              u2(k) = (dm1(k)*u1(k) + mc2(k)*(u2(k-1)-uc(k)) + mc2(k+1)*uc(k+1)) /  &
                      (dm1(k) + mc2(k+1))
              v2(k) = (dm1(k)*v1(k) + mc2(k)*(v2(k-1)-vc(k)) + mc2(k+1)*vc(k+1)) /  &
                      (dm1(k) + mc2(k+1))
           enddo

          if ( nq_skip < nq ) then
           do iq=nq_skip+1,nq
           do k=ktop+1,kbot-1
              q2(k,iq) = (dm1(k)*q1(k,iq)+mc2(k)*(q2(k-1,iq)-qc(k,iq))+mc2(k+1)*qc(k+1,iq))  &
                       / (dm1(k) + mc2(k+1))
           enddo
           enddo
          endif

!---------------------------------
! Update fields in the Cloud base
!---------------------------------
               rtmp = mc2(kbot)/dm1(kbot)
           u2(kbot) = u1(kbot) + rtmp*(u2(kbot-1)-uc(kbot))
           v2(kbot) = v1(kbot) + rtmp*(v2(kbot-1)-vc(kbot))
          if ( nq_skip < nq ) then
            do iq=nq_skip+1,nq
               q2(kbot,iq) = q1(kbot,iq) + rtmp*(q2(kbot-1,iq)-qc(kbot,iq))
            enddo
          endif
        enddo         ! end iteration

!--------------------
! Compute tendencies:
!--------------------
           do k=ktop,kbot
              utnd(i,j,k) = (u2(k) - u1(k)) * rdt
              vtnd(i,j,k) = (v2(k) - v1(k)) * rdt
! Update winds:
              uin(i,j,k) = u2(k)
              vin(i,j,k) = v2(k)
           enddo
           if ( nq_skip < nq ) then
              do iq=nq_skip+1,nq
                 do k=ktop,kbot
! qtnd is the total tendency
!                   qtnd(i,j,k,iq) = qtnd(i,j,k,iq) + (q2(k,iq) - q1(k,iq)) * rdt
                    qtnd(i,j,k,iq) =                (q2(k,iq) - q1(k,iq)) * rdt
                    qin(i,j,k,iq) = q2(k,iq)
                 enddo
              enddo
           endif
  
!          if ( conserve ) then
           if ( conserve_te ) then
           do k=ktop,kbot
! ttnd is the total tendency containing contribution from RAS
              ttnd(i,j,k) = ((u1(k)+u2(k))*(u1(k)-u2(k)) + (v1(k)+v2(k))*(v1(k)-v2(k))) &
!                          * fac_t + ttnd(i,j,k)
                           * fac_t 
              tin(i,j,k) = tin(i,j,k) + ttnd(i,j,k)*dt
           enddo
            else
             do k=ktop,kbot
               ttnd(i,j,k) = 0.0
             enddo
           endif
        endif
     enddo
  enddo

! --- diagnostics
     if ( id_utnd_cmt > 0 ) then
        used = send_data ( id_utnd_cmt, utnd, Time, is, js, 1 )
     endif
     if ( id_vtnd_cmt > 0 ) then
        used = send_data ( id_vtnd_cmt, vtnd, Time, is, js, 1 )
     endif
     if (conserve_te) then
     if ( id_ttnd_cmt > 0 ) then
        used = send_data ( id_ttnd_cmt, ttnd, Time, is, js, 1 )
     endif
     endif
     if ( id_massflux_cmt > 0 ) then
        used = send_data ( id_massflux_cmt, mc, Time, is, js, 1 )
     endif
     if ( id_detmf_cmt > 0 ) then
        used = send_data ( id_detmf_cmt, det0, Time, is, js, 1 )
     endif

 end subroutine non_local_mot

end module cu_mo_trans_mod

! <INFO>

! </INFO>


