module microphysics_mod

use fms_mod,                      only :  FATAL, error_mesg,  &
                                          write_version_number
use rotstayn_klein_mp_mod,        only :  rotstayn_klein_microp, &
                                          rotstayn_klein_microp_init,  &
                                          rotstayn_klein_microp_end
use morrison_gettelman_microp_mod, only : morrison_gettelman_microp,  &
                                          morrison_gettelman_microp_init, &
                                          morrison_gettelman_microp_end
use strat_cloud_utilities_mod,    only :  strat_cloud_utilities_init, &
                                          diag_id_type, diag_pt_type, &
                                          strat_nml_type, atmos_state_type,&
                                          cloud_state_type, particles_type,&
                                          strat_constants_type, &
                                          cloud_processes_type, &
                                          precip_state_type

implicit none
private

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

public  microphysics, microphysics_init, microphysics_end    



!------------------------------------------------------------------------
!---version number-------------------------------------------------------

character(len=128) :: version = '$Id: microphysics.F90,v 19.0 2012/01/06 20:26:13 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'


logical :: module_is_initialized = .false.




CONTAINS




!##########################################################################

subroutine microphysics_init (Nml)

type(strat_nml_type), intent(in) :: Nml

      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    write version number to output file.
!-------------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-------------------------------------------------------------------------
      call strat_cloud_utilities_init
      call rotstayn_klein_microp_init
      call morrison_gettelman_microp_init (Nml%do_pdf_clouds)

      module_is_initialized = .true.

end subroutine microphysics_init


!##########################################################################

subroutine microphysics &         
                    (idim, jdim, kdim, Nml, Constants, N3D, Atmos_state, &
                     Cloud_state, Cloud_processes, Particles, n_diag_4d, &
                     diag_4d, diag_id, diag_pt, n_diag_4d_kp1, diag_4d_kp1,&
                     ST_out, SQ_out, Precip_state, otun, ncall, &
                     qa_upd_0, SA_0, nrefuse, isamp, jsamp, ksamp, &
                     debugo, debugo0, debugo1)    

!------------------------------------------------------------------------
type(strat_constants_type),        intent(inout) :: Constants     
type(atmos_state_type),            intent(inout) :: Atmos_state
type(cloud_state_type),            intent(inout) :: Cloud_state
type(precip_state_type),           intent(inout) :: Precip_state
type(cloud_processes_type),        intent(inout) :: Cloud_processes
type(particles_type),              intent(inout) :: Particles
integer,                           intent(inout) :: nrefuse         
integer,                           intent(in)    :: idim, jdim, kdim, &
                                                    n_diag_4d, otun, ncall,&
                                                    n_diag_4d_kp1, isamp, &
                                                    jsamp, ksamp
type(strat_nml_type),              intent(in)    :: Nml
logical,                           intent(in)    :: debugo, debugo0, debugo1
real, dimension (idim,jdim,kdim),  intent(in)    :: N3D 
 real, dimension(idim,jdim,kdim),  intent(inout) :: ST_out,SQ_out,    &
                                                    qa_upd_0, SA_0   
real, dimension(idim,jdim,kdim,0:n_diag_4d),                &
                                   intent(inout) :: diag_4d
real, dimension(idim,jdim,kdim+1,0:n_diag_4d),               &
                                   intent(inout) :: diag_4d_kp1
type(diag_id_type),                intent(in)    :: diag_id
type(diag_pt_type),                intent(inout) :: diag_pt

!------------------------------------------------------------------------
!---local variables------------------------------------------------------

      integer :: i,j,k


!------------------------------------------------------------------------
!    call selected microphysics scheme.
!------------------------------------------------------------------------
      if (Constants%do_rk_microphys ) then

!------------------------------------------------------------------------
!    Rotstayn-Klein microphysics
!------------------------------------------------------------------------
        call rotstayn_klein_microp ( &
                         idim, jdim, kdim,  Nml, N3D,     &
                         Constants%overlap, Constants%dtcloud,  &
                         Constants%inv_dtcloud, Atmos_state%pfull,&
                         Atmos_state%deltpg, Atmos_state%airdens,     &
                         Constants%mask_present, Constants%mask, &
                         Atmos_state%esat0, Cloud_state%ql_in,  &
                         Cloud_state%qi_in, Cloud_state%qa_in,   &
                         Cloud_state%ql_mean, Cloud_state%qa_mean, &
                         Cloud_state%qn_mean, Atmos_state%omega,  &
                         Atmos_state%T_in, Atmos_state%U_ca, &
                         Atmos_state%qv_in, Atmos_state%qs,  &
                         Cloud_processes%D_eros, Cloud_processes%dcond_ls, &
                         Cloud_processes% dcond_ls_ice,         &
                         Cloud_processes%qvg, Atmos_state%gamma,   &
                         Cloud_processes%tmp5, Particles%drop1,    &
                         Particles%concen_dust_sub, Cloud_state%ql_upd,   &
                         Cloud_state%qi_upd, Cloud_state%qn_upd,       & 
                         Cloud_state%qi_mean, Cloud_state%qa_upd,   &
                         Atmos_state%ahuco, n_diag_4d, diag_4d, diag_id, &
                         diag_pt, n_diag_4d_kp1, diag_4d_kp1,  &
                         Constants%limit_conv_cloud_frac, &
                         Cloud_state%SA_out, Cloud_state%SN_out,        & 
                         ST_out, SQ_out, Cloud_state%SL_out,  &
                         Cloud_state%SI_out, Precip_state%rain3d,  &
                         Precip_state%snow3d, Precip_state%snowclr3d,    &
                         Precip_state%surfrain, Precip_state%surfsnow,  &
                         Cloud_processes%f_snow_berg, otun)                 
      else if (Constants%do_mg_microphys) then

!--------------------------------------------------------------------------
!     for morrison-gettelman, some additional fields are needed. for droplet
!     activation Yi's drop1 is used. 
!--------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              Particles%drop2(i,j,k) = Particles%drop1(i,j,k)*  &
                                          1.e6/Atmos_state%airdens(i,j,k)
              Cloud_processes%dcond_ls_tot(i,j,k) =   &
                               Cloud_processes%dcond_ls(i,j,k) +   &
                                      Cloud_processes%dcond_ls_ice(i,j,k) 
              Atmos_state%tn(i,j,k) = Atmos_state%T_in(i,j,k) +   &
                                                             ST_out(i,j,k)
              Atmos_state%qvn(i,j,k) = Atmos_state%qv_in(i,j,k) + &
                                                              SQ_out(i,j,k)
            end do
          end do   
        end do   

        do j=1,jdim

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency prior to
!    microphysics.
!------------------------------------------------------------------------
          if (debugo) then
            if ( j .eq. jsamp) then
              write(otun, *) " ST samp bef mg ", ST_out(isamp,jsamp,ksamp)
            end if
          endif

!-------------------------------------------------------------------------
!    call morrison-gettelman microphysics package.
!-------------------------------------------------------------------------
          call morrison_gettelman_microp( &
                  ncall, j ,idim, jdim, kdim, Nml, &
                  Constants%dtcloud, Atmos_state%pfull(:,j,:),  &
                  Atmos_state%delp(:,j,:), Atmos_state%zhalf(:,j,:),&
                  Atmos_state%tn(:,j,:),  Atmos_state%T_in(:,j,:),    &
                  Atmos_state%qvn(:,j,:), Atmos_state%qv_in(:,j,:),  &
                  Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:), &
                  Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                  Cloud_state%qa_upd(:,j,:), Atmos_state%ahuco(:,j,:),  &
                  Constants%limit_conv_cloud_frac,   &
                  Cloud_processes%dcond_ls_tot(:,j,:), &
                  Particles%drop2(:,j,:), Particles%crystal1(:,j,:), &
                  Particles%rbar_dust(:,j,:), Particles%ndust(:,j,:),    &
                  Cloud_processes%tmp5(:,j,:), Cloud_state%qa_upd(:,j,:), &
                  qa_upd_0(:,j,:), SA_0(:,j,:),   &
                  Cloud_processes%D_eros(:,j,:), Atmos_state%gamma(:,j,:), &
                  Constants%inv_dtcloud, Cloud_state%ql_in(:,j,:), &
                  Cloud_state%qi_in(:,j,:), Cloud_state%qa_in(:,j,:),     &
                  Cloud_state%qn_in(:,j,:), Cloud_state%qni_in(:,j,:),&
                  Atmos_state%rh_crit(:,j,:), ST_out(:,j,:), SQ_out(:,j,:),&
                  Cloud_state%SL_out(:,j,:), Cloud_state%SI_out(:,j,:),  &
                  Cloud_state%SN_out(:,j,:), Cloud_state%SNI_out(:,j,:),&
                  Cloud_state%SA_out(:,j,:), Precip_state%rain3d,   &
                  Precip_state%snow3d, Precip_state%surfrain,   &
                  Precip_state%surfsnow, Precip_state%qrout3d_mg(:,j,:),   &
                  Precip_state%qsout3d_mg(:,j,:), Precip_state%lsc_snow, &
                  Precip_state%lsc_rain, Precip_state%lsc_snow_size, &
                  Precip_state%lsc_rain_size,   &
                  Cloud_processes%f_snow_berg(:,j,:), &
                  n_diag_4d, diag_4d, diag_id, &
                  diag_pt, nrefuse, debugo0, debugo1, otun)    

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency after microphysics.
!------------------------------------------------------------------------
          IF (debugo) THEN
            if ( j .eq. jsamp) then
              write(otun, *) " ST samp aft mg ", ST_out(isamp,jsamp,ksamp)
            end if
          END IF
        end do

!-------------------------------------------------------------------------
!    exit with error if no valid microphysics scheme was specified.
!-------------------------------------------------------------------------
      else
        call error_mesg ('strat_cloud/microphysics', &
           'invalid strat_cloud_nml microphys_scheme option', FATAL)
      endif   

!------------------------------------------------------------------------


end subroutine microphysics


!#######################################################################

subroutine microphysics_end (Nml)

type(strat_nml_type), intent(in) :: Nml

      if (.not. module_is_initialized) return

      call rotstayn_klein_microp_end 
      call morrison_gettelman_microp_end (Nml%do_pdf_clouds)

      module_is_initialized = .false.


end subroutine microphysics_end

!#######################################################################



end module microphysics_mod
