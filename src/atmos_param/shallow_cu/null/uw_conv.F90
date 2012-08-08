
MODULE UW_CONV_MOD

  use   Time_Manager_Mod, ONLY: time_type
  use           fms_mod, only : write_version_number, ERROR_MESG, FATAL

  use  rad_utilities_mod, only : aerosol_type
  

!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: uw_conv.F90,v 18.0 2010/03/02 23:33:14 fms Exp $'
  character(len=128) :: tagname = '$Name: siena_201207 $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end, calculate_uw_closure

  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

contains

!#####################################################################
!#####################################################################

  SUBROUTINE UW_CONV_INIT( do_strat, axes, Time, kd, tracers_in_uw )
    logical,         intent(in) :: do_strat
    integer,         intent(in) :: axes(4), kd
    type(time_type), intent(in) :: Time
    logical,         intent(in) :: tracers_in_uw(:)
    
    call write_version_number (version, tagname)
    module_is_initialized = .true.
      call error_mesg('UW_CONV_INIT', &
      'This module is not supported as part of the public release', FATAL)
  end SUBROUTINE UW_CONV_INIT

!#####################################################################
!#####################################################################

  subroutine uw_conv_end
    module_is_initialized = .FALSE.
     call error_mesg('uw_conv_end', &
      'This module is not supported as part of the public release', FATAL)
  end subroutine uw_conv_end

!#####################################################################
!#####################################################################

  SUBROUTINE uw_conv(is, js, Time, tb, qv, ub, vb, pmid, pint,zmid,  & !input
       zint, q, omega, delt, pblht, ustar, bstar, qstar, land, coldT,& !input
       asol,                                                         & !input
       cush, do_strat,  skip_calculation, max_available_cf,          & !input
       tten, qvten, qlten, qiten, qaten, qnten,                      & !output
       uten, vten, rain, snow,                                       & !output
       cmf, hlflx, qtflx, pflx, liq_pflx, ice_pflx, cldql, cldqi, cldqa,cldqn, cbmfo,  & !output
        tracers, trtend, uw_wetdep)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    real, intent(in), dimension(:,:,:)   :: max_available_cf !  largest
                                     ! realizable value for uw cld frac
                                   ! after accounting for deep cld frac
    logical,intent(in), dimension(:,:)   :: skip_calculation ! do not
                                                 ! calculate where .true.
    logical,intent(in)                   :: do_strat !logical flag
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar !pbl height...
    real, intent(inout), dimension(:,:)  :: cush  ! convective scale height (m) 

    type(aerosol_type),  intent (in)     :: asol
   
    real, intent(out), dimension(:,:,:)  :: tten,qvten              ! T,qv tendencies
    real, intent(out), dimension(:,:,:)  :: qlten,qiten,qaten,qnten ! q tendencies
    real, intent(out), dimension(:,:,:)  :: uten,vten               ! u,v tendencies
   
    real, intent(out), dimension(:,:,:)  :: cldql,cldqi,cldqa, cldqn!in-updraft q
    real, intent(out), dimension(:,:,:)  :: cmf    ! mass flux at level above layer (kg/m2/s)
    real, intent(out), dimension(:,:,:)  :: pflx   ! precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: liq_pflx   ! liq precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: ice_pflx   ! solid precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: hlflx ! theta_l flux
    real, intent(out), dimension(:,:,:)  :: qtflx  ! qt  flux
    real, intent(out), dimension(:,:)    :: rain, snow
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies
    real, intent(out), dimension(:,:,:)  :: uw_wetdep       ! calculated wet depostion for tracers

      call error_mesg('UW_CONV', &
      'This module is not supported as part of the public release', FATAL)

  END SUBROUTINE UW_CONV

!#####################################################################
!#####################################################################


  SUBROUTINE calculate_uw_closure(is, js, Time, tb, qv, ub, vb, pmid, & !input
       pint, zmid, zint, q, omega, delt, pblht, ustar, bstar, qstar, & !input
       land, coldT, asol, cush,                                      & !input
                                                             cbmfo,  & !output
        tracers        )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar !pbl height...
    type(aerosol_type),  intent (in)     :: asol
    real, intent(in   ), dimension(:,:)  :: cush  ! convective scale height (m) 

   
   
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers

      call error_mesg('calculate_uw_closure', &
      'This module is not supported as part of the public release', FATAL)
  END SUBROUTINE calculate_uw_closure


end MODULE UW_CONV_MOD
