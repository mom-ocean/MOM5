                    module optical_path_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  ds
! </REVIEWER>
! <OVERVIEW>
!  Module that set up optical depth calculaiton
! </OVERVIEW>
! <DESCRIPTION>
!  radiative fluxes
! </DESCRIPTION>

!   shared modules:

use mpp_mod,             only: input_nml_file
use fms_mod,             only: open_namelist_file, fms_init, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               file_exist, write_version_number, &
                               check_nml_error, error_mesg, &
                               FATAL, close_file
use constants_mod,       only: RDGAS, RVGAS, GRAV, wtmair, &
                               avogno, pstd, diffac, tfreeze, &
                               constants_init

!   shared radiation package modules:

use rad_utilities_mod,   only: looktab, longwave_tables3_type, &
                               rad_utilities_init,  &
                               radiative_gases_type, &
                               aerosol_type,  &
                               aerosol_diagnostics_type,&
                               aerosol_properties_type, &
                               atmos_input_type, &
                               Lw_parameters,  Lw_control, &
                               Rad_control, &
                               optical_path_type, &
                               gas_tf_type, &
                               table_alloc
use longwave_params_mod, only: longwave_params_init, NBCO215,&
                               NBLY_RSB

!   radiation package modules:

use lw_gases_stdtf_mod,  only: lw_gases_stdtf_init, cfc_exact,&
                               cfc_overod, cfc_overod_part,   &
                               cfc_exact_part

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    optical_path_mod computes the optical depths and associated
!    transmission functions for various atmospheric components 
!    including radiative gases and aerosols.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

   character(len=128)  :: &
   version =  '$Id: optical_path.F90,v 19.0 2012/01/06 20:20:47 fms Exp $'
   character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public     &
         optical_path_init, optical_path_setup,     &
         optical_trans_funct_from_KS,    &
         optical_trans_funct_k_down, &
         optical_trans_funct_KE,  &
         optical_trans_funct_diag, &
         get_totch2o, get_totch2obd, &
         get_totvo2, optical_dealloc, &
         optical_path_end

private    &

!   called from optical_path_init:
         optical_ckd_init,     &

!   called from optical_path_setup:
         optical_path_ckd, optical_o3, optical_rbts,   &
         optical_h2o, cfc_optical_depth, optical_depth_aerosol


!---------------------------------------------------------------------
!---- namelist   -----

logical   :: tmp_dpndnt_h2o_lines = .true.  ! the 1200-1400 cm(-1)
                                            ! band h2o line intensities
                                            ! are temperature dependent?



namelist / optical_path_nml /    &
                                 tmp_dpndnt_h2o_lines

!-------------------------------------------------------------------
!-----  public data ----


!---------------------------------------------------------------------
!----- private  data  ----

!--------------------------------------------------------------------
!    data from former block data bs296 for self-broadened continuum
!    at 296K, band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!-------------------------------------------------------------------
real       :: v1sh2o_296, v2sh2o_296, dvsh2o_296,   &
              ssh2o_296(2000)
integer    :: nptsh2o_296

!--------------------------------------------------------------------
!  data from former block data bfh2o for foreign-broadened continuum
!    band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!--------------------------------------------------------------------
real        ::  v1fh2o, v2fh2o, dvfh2o, sfh2o(2000)
integer     ::  nptfh2o

!--------------------------------------------------------------------
!    array sfac is the frequency-dependent multiplicative factor used
!    to change the original self-broadened continuum coefficients
!    to those used in ckd2.1 or ckd2.4 (including intermediate changes).
!
!    array fscal is the frequency-dependent multiplicative factor used
!    to change the original foreign-broadened continuum coefficients
!    to those used in ckd2.1 or ckd2.4 (including intermediate changes).
!
!    array tmpfctrs is the logarithmic temperature dependence (per K)
!    of the self-broadened continuum coefficient, as a function of
!    frequency, used in all ckd AFGL continuum models.
!    the frequency ranges and intervals are as in sh2o_296.
!----------------------------------------------------------------------
real         :: sfac(2000), fscal(2000), tmpfctrs(2000)

!----------------------------------------------------------------------
!         the radfunc function (1 - exp(-h*nu/kt))/(1 + exp(-h*nu/kt))
!    is tabulated from 5 to 2995 cm-1 at intervals of 10 cm-1,
!    and from 100K to 490K at 10K intervals. note that the
!    radfn function used in ckd models equals the radfunc function
!    defined above, multiplied by nu (in cm-1).
!        the temperature derivative (at 105K to 485K, with the final
!    array value set to zero) is obtained from radfunc, and stored in
!    radfuncderiv.
!        tktab and vjtab are the respective temperature and frequency
!    points at which tabulations occurred.
!----------------------------------------------------------------------
type (longwave_tables3_type),save  :: radfunc
integer                            :: ioffh2o, nptch2o 
real                               :: vvj(2000)

!---------------------------------------------------------------------
!        fvj = foreign-broadened ckd 2.1 (ckd2.4) coefficient (including
!              all corrections), averaged over 7 specified wide
!              frequency bands in the 560-1200 cm-1 range. The average
!              is weighted by the frequency of the individual 10 cm-1
!              bands used in the averaging process.
!     fvjinw = band-averaged foreign coefficient (as in fvj) over
!              the 900-990,1070-1200 cm-1 range.
!      fvjwd = band-averaged foreign coefficient (as in fvj) over
!              the 560-800 cm-1 range.
!        svj = self-broadened ckd 2.1 (ckd2.4) coefficient (including
!              all corrections), averaged over 7 specified wide
!              frequency bands in the 560-1200 cm-1 range. The average
!              is weighted by the frequency of the individual 10 cm-1
!              bands used in the averaging process.
!     fvjinw = band-averaged self coefficient (as in svj) over
!              the 900-990,1070-1200 cm-1 range.
!      svjwd = band-averaged self coefficient (as in svj) over
!              the 560-800 cm-1 range.
!    radfnbd = the radiation function (radfn) averaged over each of
!              the 7 frequency bands: assumed to be altitude-independent
! radfnbdinw = same as radfnbd, but for the 560-800 cm-1 range.
!  radfnbdwd = same as radfnbd, but for the 900-990,1070-1200 cm-1 range
!----------------------------------------------------------------------
real      ::    fvj(7), fvjinw, fvjwd, svj(7), svjinw, svjwd,    &
                radfnbd(7), radfnbdinw, radfnbdwd

real      ::    ao3rnd(3), bo3rnd(3)

real      ::    ab15wd 

!---------------------------------------------------------------------
!  define continuum coefficients over special bands, the choices
!  depend on model architecture. the program gasbnd is used.
!
!    1) 560-800 as 1 band
!----------------------------------------------------------------------
real      ::    betawd 

!----------------------------------------------------------------------
!    3) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb on.
!    4) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
!--------------------------------------------------------------------
real, dimension (NBLY_RSB)           :: betacm

!---------------------------------------------------------------------

real, allocatable, dimension (:,:)     ::             csfah2o

!---------------------------------------------------------------------
!   the values of the molecular weights of f11 and f12 are derived
!   from elemental atomic weights adopted by the International Union of 
!   Pure and Applied Chemistry in 1961. These values are also used in 
!   the US Standard Atmosphere, 1976.
!   some previous radiative calculations at gfdl have used the
!   values  137.5, 121.0 for the molecular weights of f11 and f12.
!---------------------------------------------------------------------
real       ::  wtmf11  = 137.36855
real       ::  wtmf12  = 120.91395
real       ::  wtmf113 = 187.3765
real       ::  wtmf22  =  86.46892

!---------------------------------------------------------------------
real, dimension(2,10)  :: cpf10h2o, csf10h2o
real, dimension(2, 4)  :: cpf4h2o, csf4h2o
real, dimension(2, 2)  :: cpf2h2o, csf2h2o
real, dimension(2   )  :: cpf1h2o, csf1h2o


real      :: d622 = RDGAS/RVGAS
integer   :: NBTRG, NBTRGE
!!$integer   :: n

integer   :: ks, ke

logical   :: module_is_initialized      = .false. ! module has been
                                                  ! initialized ?


!----------------------------------------------------------------------




                              contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
!#####################################################################
! <SUBROUTINE NAME="optical_path_init">
!  <OVERVIEW>
!   Subroutine to initialize optical depth calculation and read
!   optical path namelist from input file.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize optical depth calculation and read
!   optical path namelist from input file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_path_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine optical_path_init(pref)

!--------------------------------------------------------------------
!    optical_path_init is the constructor for optical_path_mod.
!--------------------------------------------------------------------

      real, dimension(:,:), intent(in) :: pref
!--------------------------------------------------------------------
!  local variables:

      real                    :: dum
      real, dimension (NBLY_RSB)  :: dummy_n
      real, dimension (Lw_parameters%NBTRGE) :: dummy_ch4n2o
      real                    :: awide_c, bwide_c, awide_n, bwide_n, &
                                 awide, bwide
      integer, dimension(5)   :: no_h2o12001400bands = &
                                  (/ 1, 2, 4, 10, 20 /)
      real, dimension(20)     :: arndm_12001400, brndm_12001400,    &
                                 ap_12001400, bp_12001400,          &
                                 atp_12001400, btp_12001400,        &
                                 fbdlo_12001400, fbdhi_12001400
      integer                 ::  unit, ierr, io, logunit
      integer                 :: inrad, k, m
      integer                 :: subb

!---------------------------------------------------------------------
!  local variables:
!
!       dum
!       dummy
!       dummy_n
!       dummy_ch4n2o
!       ap
!       bp
!       atp
!       btp
!    define random band parameters for special bands. the choices 
!    depend on model architecture. the program gasbnd is used.
!    1)  560-800 as 1 band
!       awide_c
!       bwide_c
!       awide_n
!       bwide_n
!    end comment for above
!       awide
!       bwide
!       no_h2o12001400bands
!       arndm_12001400
!       brndm_12001400
!       ap_12001400
!       bp_12001400
!       atp_12001400
!       btp_12001400
!       fbdlo_12001400
!       fbdhi_12001400
!       unit
!       ierr
!       io
!       inrad
!       k,m
!       subb
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
      call longwave_params_init
      call lw_gases_stdtf_init(pref)

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=optical_path_nml, iostat=io)
      ierr = check_nml_error(io,"optical_path_nml")
#else
!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=optical_path_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'optical_path_nml')
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
                        write (logunit, nml=optical_path_nml)

!--------------------------------------------------------------------
!    verify that Lw_parameters%NBTRG and NBTRGE have been initialized.
!--------------------------------------------------------------------
      if (Lw_parameters%NBTRG_iz) then
        NBTRG  = Lw_parameters%NBTRG
      else
        call error_mesg ('optical_path_mod',  &
           ' Lw_parameters%NBTRG not yet initialized', FATAL)
      endif
      if (Lw_parameters%NBTRGE_iz) then
        NBTRGE = Lw_parameters%NBTRGE
      else
        call error_mesg ('optical_path_mod',  &
           ' Lw_parameters%NBTRGE not yet initialized', FATAL)
      endif

      if (nbtrge == 0 .and. tmp_dpndnt_h2o_lines) then
        call error_mesg ('optical_path_mod', &
        'cannot have temperature-dependent h2o line intensities &
             &without having separate 1200-1400 cm(-1) band(s)', FATAL) 
      endif
        
!---------------------------------------------------------------------
!    read needed data from raduiation input files.
!---------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_1992' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          inrad = open_namelist_file('INPUT/h2ocoeff_ckd_speccombwidebds_hi92')
          read (inrad,9000) awide_c   ! ckd rndm coeff for 560-800 band
          read (inrad,9000) bwide_c   ! ckd rndm coeff for 560-800 band
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          inrad = open_namelist_file('INPUT/h2ocoeff_rsb_speccombwidebds_hi92')
          read (inrad,9000) awide_n   ! rsb rndm coeff for 560-800 band
          read (inrad,9000) bwide_n   ! rsb rndm coeff for 560-800 band
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) betawd    ! rsb cont coeff for 560-800 band
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
!  rsb cont coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (betacm(k),k=1,NBLY_RSB)
        endif
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          inrad = open_namelist_file('INPUT/h2ocoeff_ckd_speccombwidebds_hi00')
          read (inrad,9000) awide_c   ! ckd rndm coeff for 560-800 band
          read (inrad,9000) bwide_c   ! ckd rndm coeff for 560-800 band
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          inrad = open_namelist_file('INPUT/h2ocoeff_rsb_speccombwidebds_hi00')
          read (inrad,9000) awide_n   ! rsb rndm coeff for 560-800 band
          read (inrad,9000) bwide_n   ! rsb rndm coeff for 560-800 band
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) betawd    ! rsb cont coeff for 560-800 band
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
!  rsb cont coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (betacm(k),k=1,NBLY_RSB)
        endif
      endif
9000  format(5e14.6)
      call close_file (inrad)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        awide = awide_c
        bwide = bwide_c
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        awide = awide_n
        bwide = bwide_n
      endif

!---------------------------------------------------------------------
!    compute a*b for computational frequency bands for the 15 um
!    region, as 1 band (ab15wd)
!---------------------------------------------------------------------
      ab15wd = awide*bwide

      if (trim(Lw_control%linecatalog_form) == 'hitran_1992') then
        inrad = open_namelist_file('INPUT/o39001200_hi92_data')
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000') then
        inrad = open_namelist_file('INPUT/o39001200_hi00_data')
      endif
      read (inrad,2001) (ao3rnd(k),k=1,3)
      read (inrad,2001) (bo3rnd(k),k=1,3)

!---------------------------------------------------------------------
!    verify that Lw_control%do_ch4 has been initialized.
!--------------------------------------------------------------------
      if (Lw_control%do_ch4_iz) then
      else
        call error_mesg ( 'optical_path_mod',  &
                      ' do_ch4 not yet initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    verify that Lw_control%do_n2o has been initialized.
!--------------------------------------------------------------------
      if (Lw_control%do_n2o_iz) then
      else
        call error_mesg ( 'optical_path_mod',  &
                      ' do_n2o not yet initialized', FATAL)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( csfah2o(2, NBTRGE) )
        if (trim(Lw_control%linecatalog_form) == 'hitran_1992') then
          inrad = open_namelist_file('INPUT/h2o12001400_hi92_data')
        else if (trim(Lw_control%linecatalog_form) ==    &
                                                    'hitran_2000') then
          inrad = open_namelist_file('INPUT/h2o12001400_hi00_data')
        endif

!----------------------------------------------------------------------
!     read in random coefficients for 1200-1400 freq region, spacing
!     through the data until  those appropriate for NBTRGE h2o bands
!     are reached. note: unless a continuum is inserted beyond 1200
!     cm-1, the band coefficients are independent of continuum type.
!---------------------------------------------------------------------
        do subb = 1,5    ! 5 = no. band divisions in h2o 1200-1400 data
          if (NBTRGE == no_h2o12001400bands(subb)) then
!  read and process data for sub-band number from data matching NBTRGE
!  then exit subb loop
            read (inrad,2001) (arndm_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (brndm_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (ap_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (bp_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (atp_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (btp_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (fbdlo_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (fbdhi_12001400(k),k=1,NBTRGE)
            do m=1,NBTRGE
              csfah2o(1,m) =          atp_12001400(m)
              csfah2o(2,m) =          btp_12001400(m)
            end do
            exit
          else if (subb < 5) then 
!  read data for sub-band number from  data not matching NBTRGE
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
            read (inrad,2001)    &
                       (dummy_ch4n2o(k),k=1,no_h2o12001400bands(subb))
          else
!    failure of any sub-band number to match NBTRGE
            call error_mesg ('optical_path_mod',  &
              'NBTRGE is inconsistent with available data', FATAL)
          endif
        end do
2001  format(5e14.6)
        call close_file(inrad)
      endif

!------------------------------------------------------------------
!
!------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        call optical_ckd_init
      endif

!------------------------------------------------------------------
!    mark the module as initialized.
!------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------



end subroutine optical_path_init




!###################################################################
! <SUBROUTINE NAME="optical_path_setup">
!  <OVERVIEW>
!   Subroutine to prepare optical path calculation, such as memory
!   allocation.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to prepare optical path calculation, such as memory
!   allocation.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_path_setup (is, ie, js, je,  Atmos_input, &
!                            Rad_gases, Aerosol, Aerosol_props, Optical)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   Latitude and longitude bound of model physics window.
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases input data
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol radiative properties input data
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties output (extinction coefficient,
!   single scattering albedo and asymmetry parameter in different
!   bands)
!  </INOUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   optical path output
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_path_setup (is, ie, js, je, Atmos_input, &
                               Rad_gases, Aerosol, Aerosol_props,  &
                               Aerosol_diags, Optical, &
                               including_aerosols)  

!------------------------------------------------------------------
!
!------------------------------------------------------------------

integer, intent(in)                          :: is, ie, js, je
type(atmos_input_type),        intent(in)    :: Atmos_input
type(radiative_gases_type),    intent(in)    :: Rad_gases
type(aerosol_type),            intent(in)    :: Aerosol      
type(aerosol_properties_type), intent(inout) :: Aerosol_props      
type(aerosol_diagnostics_type), intent(inout) :: Aerosol_diags      
type(optical_path_type),       intent(inout) :: Optical     
logical,                   intent(in)            :: including_aerosols  

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je
!      Atmos_input
!      Rad_gases
!      Aerosol
!
!  intent(inout) variables:
!
!      Aerosol_props
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:
 
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3) )  :: press, pflux, &
                                                       temp, tflux, &
                                                       atmden, vv

      real, dimension (size(Atmos_input%press,1),   &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3) - 1 )  ::   &
                                                       rh2o, deltaz
      real, dimension (size(Atmos_input%press,3) ) :: bsum

      integer      :: n_aerosol_bands
      integer      :: k, i, j, n
      integer      :: ix, jx, kx
      integer      :: israd, ierad, jsrad, jerad

!--------------------------------------------------------------------
!  local variables:
!
!       press
!       pflux
!       temp
!       tflux
!       atmden
!       vv             layer-mean pressure in atmospheres. due to quad-
!                      rature considerations, this does not equal the 
!                      pressure at the data level (press).
!       rh2o
!       deltaz
!       n_aerosol_bands
!       i,k
!       ix,jx,kx
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
      ix = ie -is + 1
      jx = je -js +1
      israd = 1
      ierad = ix
      jsrad = 1
      jerad = jx
      ks    = 1
      kx = size(Atmos_input%press,3) - 1
      ke    = kx

!  convert press and pflux to cgs.
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
      pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)
      deltaz = Atmos_input%deltaz
      temp = Atmos_input%temp
      rh2o = Atmos_input%rh2o
      tflux = Atmos_input%tflux

!--------------------------------------------------------------------
!    atmden   =  atmospheric density, in gm/cm**2, for each of the
!                KMAX layers.
!-------------------------------------------------------------------
      allocate (Optical%wk       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
      allocate (Optical%rh2os    (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
      allocate (Optical%rfrgn    (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
      allocate (Optical%tfac     (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
      allocate (Optical%avephi   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

      if (NBTRGE > 0) then
        allocate (Optical%avephif(ISRAD:IERAD, JSRAD:JERAD,    &
                                  KS:KE+1, NBTRGE) )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      Optical%wk      = 0.                                      
      Optical%rh2os   = 0.                                      
      Optical%rfrgn   = 0.                                       
      Optical%tfac    = 0.                                       
      Optical%avephi   = 0.

      if (NBTRGE > 0) then
        Optical%avephif   = 0.
      endif
 
!----------------------------------------------------------------------
!    define the layer-mean pressure in atmospheres (vv) and the layer 
!    density (atmden). 
!----------------------------------------------------------------------
      do k=KS,KE
        atmden(:,:,k) = (pflux(:,:,k+1) - pflux(:,:,k))/(1.0E+02*GRAV)
        vv(:,:,k)     = 0.5E+00*(pflux(:,:,k+1) + pflux(:,:,k)  )/pstd
      end do

!----------------------------------------------------------------------
!     compute optical paths.
!----------------------------------------------------------------------
      call optical_h2o (pflux, atmden, vv, press, temp, rh2o,   &
                        tflux, Optical)

!---------------------------------------------------------------------
!    call optical_ckd2p1 to determine self- and foreign-broadened h2o
!    continuum paths, for the given temperature, pressure and mixing
!    ratio, over the predetermined frequency range for the ckd2.1 
!    continuum. call optical_roberts for self-broadened continuum
!    paths for the rsb (Roberts) continuum.
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        call optical_path_ckd  (atmden, press, temp, rh2o, Optical)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        call optical_rbts  (temp, rh2o, Optical)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      call optical_o3 (atmden, Rad_gases%qo3, vv, Optical)

!--------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_optical_depth (atmden, Rad_gases, Optical) 
      endif

!---------------------------------------------------------------------
!    compute aerosol layer transmission functions for all layers.
!    option predlwaer is planned, but not yet available. when  it 
!    becomes available,  aeralb, aerext and aerconc will be additional 
!    arguments going to aertau.
!---------------------------------------------------------------------
      if (including_aerosols .or. Rad_control%volcanic_lw_aerosols) then
        n_aerosol_bands = Lw_parameters%n_lwaerosol_bands

!---------------------------------------------------------------------
!    allocate space for and then retrieve the aerosol mixing ratios and
!    aerosol optical properties from the aerosol module.
!--------------------------------------------------------------------
        allocate (Optical%totaerooptdep (ix,jx,kx+1, N_AEROSOL_BANDS))
        allocate (Optical%aerooptdep_KE_15 (ix, jx ) )
        Optical%totaerooptdep = 0.                              
        Optical%aerooptdep_KE_15 = 0.           
      else
        n_aerosol_bands = 0
      endif
!-------------------------------------------------------------------
!    for each aerosol frequency band, retrieve aerosol optical proper-
!    ties for each aerosol category. then call optical_depth_aerosol 
!    to compute for aerosol optical depth. 
!-------------------------------------------------------------------
      do n=1,n_aerosol_bands  !  loop on aerosol frequency bands
        if (including_aerosols) then
          call optical_depth_aerosol (js, Atmos_input, n, Aerosol,   &
                                      Aerosol_props, Aerosol_diags, &
                                      Optical)
        endif   ! (including_aerosols)

        if (Rad_control%volcanic_lw_aerosols) then
          if (size(Aerosol_props%lw_ext,4) /= 0) then
            do j=1,jx
              do i=1,ix
                bsum(1) = 0.0
                do k=2,kx+1
                  if (n == 5) then
                    Aerosol_diags%lw_extopdep_vlcno(i,j,k,1) =  &
                                   Aerosol_props%lw_ext(i,j,k-1,n)*&
                                   Atmos_input%deltaz(i,j,k-1)
                    Aerosol_diags%lw_absopdep_vlcno(i,j,k,1) =  &
                          Aerosol_diags%lw_extopdep_vlcno(i,j,k,1) 
!! NOT CURRENTLY AVAILABLE IN SEA LW CODE -- lw_ssa not processed
!                     Aerosol_diags%lw_absopdep_vlcno(i,j,k,2) =  &
!                              (1.0-Aerosol_props%lw_ssa(i,j,k-1,n))*  &
!                                  Aerosol_props%lw_ext(i,j,k-1,n)*&
!                                  Atmos_input%deltaz(i,j,k-1)
                  endif
                  if (n == 6) then
                    Aerosol_diags%lw_extopdep_vlcno(i,j,k,2) =  &
                                   Aerosol_props%lw_ext(i,j,k-1,n)*&
                                   Atmos_input%deltaz(i,j,k-1)
                    Aerosol_diags%lw_absopdep_vlcno(i,j,k,2) =  &
                           Aerosol_diags%lw_extopdep_vlcno(i,j,k,2) 
!! NOT CURRENTLY AVAILABLE IN SEA LW CODE -- lw_ssa not processed
!                     Aerosol_diags%lw_absopdep_vlcno(i,j,k,1) =  &
!                           (1.0-Aerosol_props%lw_ssa(i,j,k-1,n))*  &
!                                  Aerosol_props%lw_ext(i,j,k-1,n)*&
!                                  Atmos_input%deltaz(i,j,k-1)
                  endif
                  bsum(k) = bsum(k-1) +    &
                            Aerosol_props%lw_ext(i,j,k-1,n)*&
                            Atmos_input%deltaz(i,j,k-1)
                end do
                Optical%totaerooptdep(i,j,2:kx+1,n) =    &
                        Optical%totaerooptdep(i,j,2:kx+1,n) +   &
                        bsum(2:kx+1)
                if (n == n_aerosol_bands) then
                  Optical%aerooptdep_KE_15(i,j) = &
                                Optical%aerooptdep_KE_15(i,j) +  &
                                Aerosol_props%lw_ext(i,j,kx,n)* &
                                Atmos_input%deltaz(i,j,kx)
                endif
              end do   
            end do
          endif ! (size)
        endif  ! (volcanic_lw_aerosols)
      end do  ! (n_aerosol_bnads)

!---------------------------------------------------------------------
       


end subroutine  optical_path_setup



!####################################################################
! <SUBROUTINE NAME="optical_trans_funct_from_KS">
!  <OVERVIEW>
!   Subroutine to compute transmission function from level KS to another
!   level
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute transmission function from level KS to another
!   level
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_trans_funct_from_KS (Gas_tf, to3cnt, overod, Optical, &
!                                        cnttaub1, cnttaub2, cnttaub3)
!  </TEMPLATE>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission functions
!  </INOUT>
!  <OUT NAME="to3cnt" TYPE="real">
!   Ozone continuum transmission function
!  </OUT>
!  <OUT NAME="overod" TYPE="real">
!   Transmission function due to h2o continuum and aerosol
!  </OUT> 
!  <INOUT NAME="Optical" TYPE="real">
!   Optical depth function
!  </INOUT>
!  <OUT NAME="cnttaub1, cnttaub2, cnttaub3" TYPE="real">
!   Transmission functions of gas continuum
!  </OUT>
! </SUBROUTINE>
!
subroutine optical_trans_funct_from_KS (Gas_tf, to3cnt, overod,   &
                                        Optical, cnttaub1, cnttaub2, &
                                        cnttaub3, including_aerosols)  

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension (:,:,:), intent(out)   ::  to3cnt, overod, &
                                           cnttaub1, cnttaub2, cnttaub3
type(optical_path_type), intent(inout) ::  Optical
type(gas_tf_type),       intent(inout) ::  Gas_tf 
logical,                   intent(in)            :: including_aerosols  
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(inout) variables:
!
!     Optical
!     Gas_tf
!
!   intent(out) variables:
!
!     to3cnt
!     overod
!     cnttaub1
!     cnttaub2
!     cnttaub3
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)) ::   &
                                               tmp1, tmp2, tmp3,   &
                                               totch2o_tmp,  &
                                               totaer_tmp, tn2o17

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)-1) ::    cfc_tf

      integer    :: m

!---------------------------------------------------------------------
!  local variables:
!
!     tmp1
!     tmp2
!     tmp3
!     totch2o_tmp
!     totaer_tmp
!     tn2o17
!     cfc_tf
!     m
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
              'module has not been initialized', FATAL )
      endif
!-----------------------------------------------------------------------
!    compute transmission functions in 990-1070 cm-1 range, including
!    ozone and h2o continuum, from level KS to all other levels. 
!------------------------------------------------------------------
      if (Lw_control%do_o3) then
        tmp1  (:,:,KS:KE) = bo3rnd(2)*Optical%tphio3(:,:,KS+1:KE+1)/  &
                            Optical%toto3(:,:,KS+1:KE+1)
        tmp2(:,:,KS:KE) = 0.5*(tmp1(:,:,KS:KE)*(SQRT(1.0E+00 +   &
                              (4.0E+00*ao3rnd(2)*  &
                               Optical%toto3(:,:,KS+1:KE+1))/  &
                               tmp1(:,:,KS:KE))  - 1.0E+00))
      else
        tmp2(:,:,KS:KE)  = 0.0E+00
      endif

      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        call get_totch2obd(6, Optical, totch2o_tmp)
        tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + diffac*   &
                          totch2o_tmp(:,:,KS+1:KE+1)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + betacm(14)*   &
                          Optical%totvo2(:,:,KS+1:KE+1)
      endif

      if (including_aerosols) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,6)
        tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +    &
                          totaer_tmp   (:,:,KS+1:KE+1)
      endif
      to3cnt(:,:,KS) = 1.0
      to3cnt(:,:,KS+1:KE+1) = EXP(-1.0E+00*tmp2(:,:,KS:KE))
 
!--------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (6, Optical, cfc_tf)
        to3cnt(:,:,KS+1:KE+1) = to3cnt(:,:,KS+1:KE+1)* cfc_tf(:,:,KS:KE)
      endif

!---------------------------------------------------------------------
!    compute transmission function in the 560-800 cm-1 range
!    evaluate  optical depth contributions 
!    add contributions from h2o(lines) and h2o(continuum).
!    h2o(continuum) contributions are either Roberts or CKD2.1
!---------------------------------------------------------------------
      if (Lw_control%do_h2o) then
        tmp1(:,:,KS:KE) = SQRT(ab15wd*Optical%totphi(:,:,KS+1:KE+1)) 
      else
        tmp1(:,:,:) = 0.0
      endif
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*   &
                          Optical%totch2obdwd(:,:,KS+1:KE+1)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + betawd*   &
                          Optical%totvo2     (:,:,KS+1:KE+1)
      endif

!--------------------------------------------------------------------
!    add contribution from longwave aerosols (if desired).
!--------------------------------------------------------------------
      if (including_aerosols) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:, 9)
        tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) +    &
                          totaer_tmp(:,:,KS+1:KE+1)
      endif
 
!----------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
      overod(:,:,KS) = 1.0
      overod(:,:,KS+1:KE+1) = EXP(-1.0E+00*tmp1     (:,:,KS:KE))

!---------------------------------------------------------------------
!    add contribution from the 17 um n2o band (if desired).
!    the expression with tn2o17 retains the 560-630 cm-1 equi-
!    valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_n2o) then
        tn2o17(:,:,ks+1:ke+1) = Gas_tf%tn2o17(:,:,ks+1:ke+1)
        if (NBCO215 .EQ. 2) then
          overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1) *    &
                                  (130./240. + 110./240.*   &
                                  tn2o17(:,:,KS+1:KE+1))
        elseif (NBCO215 .EQ. 3) then
          overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1)*(170./240. +  &
                                  70./240.*tn2o17(:,:,KS+1:KE+1))
        endif
      endif

!--------------------------------------------------------------------- 
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!--------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_overod (Optical, cfc_tf)
        overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1)*cfc_tf(:,:,KS:KE)
      endif 

!----------------------------------------------------------------------
!    compute continuum band transmission functions from level KS to
!    other levels (cnttau). the continuum transmission function from
!    level k to kp (contod) equals cnttau for k=KS, so is not
!    evaluated here. for all other levels k, contod is obtained by
!    division of relevant values of cnttau.
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        call get_totch2obd(4, Optical, totch2o_tmp)
        tmp1(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
        call get_totch2obd(5, Optical, totch2o_tmp)
        tmp2(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
        call get_totch2obd(7, Optical, totch2o_tmp)
        tmp3(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp1(:,:,KS:KE) = betacm(12)*Optical%totvo2(:,:,KS+1:KE+1)
        tmp2(:,:,KS:KE) = betacm(13)*Optical%totvo2(:,:,KS+1:KE+1)
        tmp3(:,:,KS:KE) = betacm(15)*Optical%totvo2(:,:,KS+1:KE+1)
      endif

      if (including_aerosols) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,4)
        tmp1(:,:,KS:KE) =  tmp1(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1  )
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,5)
        tmp2(:,:,KS:KE) =  tmp2(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,7)
        tmp3(:,:,KS:KE) =  tmp3(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
      endif

      cnttaub1(:,:,KS) = 1.0                       
      cnttaub2(:,:,KS) = 1.0                       
      cnttaub3(:,:,KS) = 1.0                       
      cnttaub1(:,:,KS+1:KE+1) = EXP(-1.0*tmp1(:,:,KS:KE))
      cnttaub2(:,:,KS+1:KE+1) = EXP(-1.0*tmp2(:,:,KS:KE))
      cnttaub3(:,:,KS+1:KE+1) = EXP(-1.0*tmp3(:,:,KS:KE))

!---------------------------------------------------------------------
!    if cfcs are included, add transmission functions for f11, f12,
!    f113, and f22.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (4, Optical, cfc_tf)
        cnttaub1(:,:,KS+1:KE+1) = cnttaub1(:,:,KS+1:KE+1)*  &
                                  cfc_tf(:,:,KS:KE)
        call cfc_exact (5, Optical, cfc_tf)
        cnttaub2(:,:,KS+1:KE+1) = cnttaub2(:,:,KS+1:KE+1)*  &
                                  cfc_tf(:,:,KS:KE)
        call cfc_exact (7, Optical, cfc_tf)
        cnttaub3(:,:,KS+1:KE+1) = cnttaub3(:,:,KS+1:KE+1)*   &
                                  cfc_tf(:,:,KS:KE)
      endif 
 
!----------------------------------------------------------------------
!    evaluate h2o (mbar*phibar) between level KS and other levels.
!----------------------------------------------------------------------
      Optical%avephi(:,:,KS:KE) = Optical%totphi(:,:,KS+1:KE+1)
 
!----------------------------------------------------------------------
!    the evaluation of emiss over the layer between data level (KS)
!    and flux level (KE+1) is done by averaging E2 functions referring
!    to the top and bottom of the layer. a special value of (mbar*
!    phibar) is required; it is stored in the (otherwise vacant)
!    KE+1'th position of avephi.
!----------------------------------------------------------------------
      Optical%avephi(:,:,KE+1) = Optical%avephi(:,:,KE-1) +   &
                                 Optical%emx1(:,:)

!----------------------------------------------------------------------
!    if h2o lines in the 1200-1400 range are assumed to have a temp-
!    erature dependent intensity, similar evaluation for (mbar*phibar)
!    is performed, with a special value for the lowest layer
!----------------------------------------------------------------------
      if (NBTRGE > 0) then
        if (tmp_dpndnt_h2o_lines) then
          do m=1,NBTRGE
            Optical%avephif(:,:,KS:KE,m) =     &
                                       Optical%tphfh2o(:,:,KS+1:KE+1,m)
          end do
          do m=1,NBTRGE
            Optical%avephif(:,:,KE+1,m) =   &
                                        Optical%avephif(:,:,KE-1,m) +  &
                                        Optical%emx1f(:,:,m)
          end do
        else 
          do m=1,NBTRGE
            Optical%avephif(:,:,KS:KE,m) =     &
                                       Optical%avephi(:,:,KS:KE)
          end do
          do m=1,NBTRGE
            Optical%avephif(:,:,KE+1,m) = Optical%avephi(:,:,KE+1)
          end do
        endif
      endif

!----------------------------------------------------------------------


end subroutine optical_trans_funct_from_KS




!####################################################################
! <SUBROUTINE NAME="optical_trans_funct_k_down">
!  <OVERVIEW>
!   Subroutine to compute transmission function downward from level k
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute transmission function downward from level k
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_trans_funct_k_down (Gas_tf, k,                     &
!                                    to3cnt, overod, Optical)
!  </TEMPLATE>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission functions
!  </INOUT>
!  <IN NAME="k" TYPE="integer">
!   The data level from which downward transmission functions are computed
!  </IN>
!  <OUT NAME="to3cnt" TYPE="real">
!   Ozone continuum transmission function
!  </OUT>
!  <OUT NAME="overod" TYPE="real">
!   Transmission function due to h2o continuum and aerosol
!  </OUT> 
!  <INOUT NAME="Optical" TYPE="real">
!   Optical depth function
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_trans_funct_k_down (Gas_tf, k, to3cnt, overod,   &
                                       Optical,including_aerosols)  

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

integer,                 intent (in)    :: k
real, dimension (:,:,:), intent(out)    :: to3cnt, overod
type(optical_path_type), intent(inout)  :: Optical
type(gas_tf_type),       intent(inout)  :: Gas_tf 
logical,                   intent(in)            :: including_aerosols  

!---------------------------------------------------------------------
!   intent(in) variable:
!        
!       k
!
!   intent(inout) variables:
!
!       Optical
!       Gas_tf
!
!   intent(out) variables:
!
!       to3cnt
!       overod
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)) ::    &
                                                avmo3, avpho3, tmp1, &
                                                tmp2, avvo2,  &
                                                avckdwd, avckdo3, &
                                                avaero3, totch2o_tmp,  &
                                                totaer_tmp, tn2o17

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)-1) ::     cfc_tf

      integer       :: kp, m

!---------------------------------------------------------------------
!   local variables:
!
!       avmo3
!       avpho3
!       tmp1
!       tmp2
!       avvo2
!       avchdwd
!       avckdo3
!       avaero3  
!       totch2o_tmp
!       totaer_tmp
!       tn2o17
!       cfc_tf
!       kp
!       m
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
              'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        call get_totch2obd(6, Optical, totch2o_tmp)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (including_aerosols) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,6)
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do kp=1,KE+1-k
        avmo3 (:,:,kp+k-1) = Optical%toto3 (:,:,kp+k) -    &
                             Optical%toto3 (:,:,k)
        avmo3 (:,:,kp+k-1) = max(avmo3 (:,:,kp+k-1),1.0e-10)
        avpho3(:,:,kp+k-1) = Optical%tphio3(:,:,kp+k) -    &
                             Optical%tphio3(:,:,k) 
        avpho3 (:,:,kp+k-1) = max(avpho3 (:,:,kp+k-1),1.0e-12)
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          avckdwd(:,:,kp+k-1) = Optical%totch2obdwd(:,:,kp+k) -   &
                                Optical%totch2obdwd(:,:,k)
          avckdo3(:,:,kp+k-1) = totch2o_tmp(:,:,kp+k) -  &
                                totch2o_tmp(:,:,k)
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          avvo2 (:,:,kp+k-1) = Optical%totvo2(:,:,kp+k) -   &
                               Optical%totvo2(:,:,k)
        endif 
        if (including_aerosols) then
          avaero3(:,:,kp+k-1) =  &
                       totaer_tmp   (:,:,kp+k) - totaer_tmp   (:,:,k)
         endif
       end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
       do kp=1,KE+1-k
         Optical%avephi  (:,:,kp+k-1) = Optical%totphi(:,:,kp+k) -  &
                                        Optical%totphi(:,:,k)
       end do
       Optical%avephi (:,:,KE+1) = Optical%avephi(:,:,KE-1) +   &
                                   Optical%emx1(:, :)

!---------------------------------------------------------------------
!    if h2o lines in the 1200-1400 range are assumed to have a temp-
!    erature dependent intensity, similar evaluation for (mbar*phibar)
!    is performed, with a special value for the lowest layer
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        if (tmp_dpndnt_h2o_lines) then
          do m=1,NBTRGE
            do kp=1,KE+1-k
              Optical%avephif(:,:,kp+k-1,m) =   &
                                     Optical%tphfh2o(:,:,kp+k,m) -  &
                                     Optical%tphfh2o(:,:,k,   m)
            end do
            Optical%avephif(:,:,KE+1,m) =   &
                                         Optical%avephif(:,:,KE-1,m) + &
                                         Optical%emx1f(:,:,m)
          end do
        else
          do m=1,NBTRGE
            do kp=1,KE+1-k
              Optical%avephif(:,:,kp+k-1,m) = Optical%avephi(:,:,kp+k-1)
            end do
            Optical%avephif(:,:,KE+1,m) = Optical%avephi(:,:,KE+1) 
          end do
        endif
      endif

!----------------------------------------------------------------------
!    compute transmission function in the 560-800 cm-1 range
!    evaluate  optical depth contributions 
!
!    add contributions from h2o(lines) and h2o(continuum).
!    h2o(continuum) contributions are either Roberts or CKD2.1
!----------------------------------------------------------------------
      if (Lw_control%do_h2o) then 
        tmp1(:,:,k:KE) = SQRT(ab15wd*Optical%avephi(:,:,k:KE)) 
      else
        tmp1(:,:,k:KE) = 0.0
      endif

      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        tmp1(:,:,k:KE) = tmp1(:,:,k:KE) + diffac*   &
                         avckdwd    (:,:,k:KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp1(:,:,k:KE) = tmp1(:,:,k:KE) + betawd*   &
                         avvo2      (:,:,k:KE)
      endif

!-------------------------------------------------------------------
!    add contribution from longwave aerosols (if desired).
!-------------------------------------------------------------------
      if (including_aerosols) then
        totaer_tmp      (:,:,:) = Optical%totaerooptdep   (:,:,:,9)
        do kp=k,KE
          tmp1(:,:,kp) = tmp1(:,:,kp) +    &
                         (totaer_tmp(:,:,kp+1) - totaer_tmp(:,:,k) )
        end do
      endif

!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!--------------------------------------------------------------------
      overod(:,:,k+1:KE+1) = EXP(-1.0E+00*tmp1(:,:,k:KE))

!----------------------------------------------------------------------
!    add contribution from the 17 um n2o band (if desired).
!    the expression with tn2o17 retains the 560-630 cm-1 equi-
!    valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_n2o) then
        tn2o17(:,:,k+1:ke+1) = Gas_tf%tn2o17(:,:,k+1:ke+1)
        if (NBCO215 .EQ. 2) then
          overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1) *(130./240. +  &
                                 110./240.*tn2o17(:,:,k+1:KE+1))
        elseif (NBCO215 .EQ. 3) then
          overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1)*(170./240. + &
                                 70./240.*tn2o17(:,:,k+1:KE+1))
        endif
      endif

!----------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_overod_part ( Optical, cfc_tf, k)
        overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1)*cfc_tf(:,:,k:KE)
      endif

!--------------------------------------------------------------------
!    compute transmission functions in 990-1070 cm-1 range, including
!    ozone and h2o continuum, from level k to all other levels. 
!---------------------------------------------------------------------
      if (Lw_control%do_o3) then
        tmp1  (:,:,k:KE) = bo3rnd(2)*avpho3(:,:,k:KE)/avmo3(:,:,k:KE)
        tmp2(:,:,k:KE) = 0.5*(tmp1(:,:,k:KE)*(SQRT(1.0E+00 + (4.0E+00* &
                           ao3rnd(2)*avmo3(:,:,k:KE))/tmp1(:,:,k:KE))  &
                           - 1.0E+00))
      else
        tmp2(:,:,k:KE) = 0.0
      endif

      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        tmp2(:,:,k:KE) = tmp2(:,:,k:KE) + diffac*   &
                         avckdo3  (:,:,k:KE) 
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp2(:,:,k:KE) = tmp2(:,:,k:KE) + betacm(14)*   &
                         avvo2 (:,:,k:KE)
      endif
      if (including_aerosols) then
        tmp2(:,:,k:KE) = tmp2(:,:,k:KE) +   &
                         avaero3      (:,:,k:KE)
      endif
      to3cnt(:,:,k+1:KE+1) = EXP(-1.0E+00*tmp2(:,:,k:KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact_part (6, Optical, cfc_tf, k)
        to3cnt(:,:,k+1:KE+1) = to3cnt(:,:,k+1:KE+1)*cfc_tf(:,:,k:KE)
      endif 
!---------------------------------------------------------------------


end subroutine optical_trans_funct_k_down



!#################################################################
! <SUBROUTINE NAME="optical_trans_funct_KE">
!  <OVERVIEW>
!   Subroutine to compute transmission function from level KE
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute transmission function from level KE
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_trans_funct_KE (Gas_tf, to3cnt, Optical, overod)
!  </TEMPLATE>
!  <INOUT NAME="Gas_tf" TYPE="gas_tf_type">
!   Gas transmission functions
!  </INOUT>
!  <OUT NAME="to3cnt" TYPE="real">
!   Ozone continuum transmission function
!  </OUT>
!  <OUT NAME="overod" TYPE="real">
!   Transmission function due to h2o continuum and aerosol
!  </OUT> 
!  <INOUT NAME="Optical" TYPE="real">
!   Optical depth function
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_trans_funct_KE (Gas_tf, to3cnt, Optical, overod, &
                                   including_aerosols)  

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, dimension (:,:,:), intent(out)   :: to3cnt, overod
type(optical_path_type), intent(inout) :: Optical
type(gas_tf_type),       intent(inout) :: Gas_tf 
logical,                   intent(in)            :: including_aerosols  

!---------------------------------------------------------------------
!   intent(inout) variables:
!
!     Optical
!     Gas_tf
!
!   intent(out) variables:
!
!     to3cnt
!     overod
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)) ::    &
                                             tmp1, tmp2, tn2o17

      real, dimension (size(to3cnt,1), size(to3cnt,2), &
                       size(to3cnt,3)-1) ::    &
                                             cfc_tf

      real, dimension (size(to3cnt,1), size(to3cnt,2)) :: &
                                             aerooptdep_KE_15

!---------------------------------------------------------------------
!   local variables:
!
!      tmp1
!      tmp2
!      tn2o17
!      cfc_tf
!      aer_tmp
!      aerooptdep_KE_15
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
             'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------------
!    compute transmission function in the 560-800 cm-1 range. evaluate 
!    optical depth contributions. add contributions from h2o(lines) and
!    h2o(continuum). h2o(continuum) contributions are either Roberts 
!    or CKD2.1 or CKD2.4.
!----------------------------------------------------------------------
      if (Lw_control%do_h2o) then
        tmp1     (:,:,KE) = SQRT(ab15wd*Optical%var2  (:,:,KE)) 
      else
        tmp1     (:,:,KE) = 0.0
      endif
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        tmp1(:,:,KE) = tmp1(:,:,KE) + diffac*   &
                       Optical%xch2obdwd   (:,:,KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp1(:,:,KE) = tmp1(:,:,KE) + betawd*  &
                       Optical%cntval     (:,:,KE)
      endif

!---------------------------------------------------------------------
!    add contribution from longwave aerosols (if desired).
!---------------------------------------------------------------------
      if (including_aerosols) then
        aerooptdep_KE_15(:,:) = Optical%aerooptdep_KE_15(:,:)
        tmp1(:,:,KE) = tmp1(:,:,KE) + aerooptdep_KE_15(:,:)  
      endif
 
!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
      overod(:,:,KE+1) = EXP(-1.0E+00*tmp1     (:,:,KE))
 
!---------------------------------------------------------------------
!    add contribution from the 17 um n2o band (if desired).
!    the expression with tn2o17 retains the 560-630 cm-1 equi-
!    valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_n2o) then
        tn2o17(:,:,ke+1    ) = Gas_tf%tn2o17(:,:,ke+1)
        if (NBCO215 .EQ. 2) then
          overod(:,:,KE+1) = overod(:,:,KE+1) *  &
                             (130./240. + 110./240.*tn2o17(:,:,KE+1))
        else if (NBCO215 .EQ. 3) then
          overod(:,:,KE+1) = overod(:,:,KE+1) *   &
                             (170./240. + 70./240.*tn2o17(:,:,KE+1))
        endif
      endif

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_overod_part (Optical, cfc_tf, KE)
        overod(:,:,KE+1) = overod(:,:,KE+1)*cfc_tf(:,:,KE)
      endif 

!-----------------------------------------------------------------------
!    compute transmission functions in 990-1070 cm-1 range, including
!    ozone and h2o continuum, from level KS to all other levels. 
!---------------------------------------------------------------------
      if (Lw_control%do_o3) then
        tmp1  (:,:,KE) = bo3rnd(2)*Optical%var4(:,:,KE)/  &
                         Optical%var3(:,:,KE)
        tmp2(:,:,KE) = 0.5*(tmp1(:,:,KE)*(SQRT(1.0E+00 + (4.0E+00*  &
                       ao3rnd(2)*Optical%var3 (:,:,KE))/  &
                       tmp1(:,:,KE)) - 1.0E+00))
      else
        tmp2(:,:,KE) = 0.0
      endif

      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        tmp2(:,:,KE) = tmp2(:,:,KE) + diffac*Optical%xch2obd  (:,:,KE,6)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        tmp2(:,:,KE) = tmp2(:,:,KE) + betacm(14)*Optical%cntval (:,:,KE)
      endif

      to3cnt(:,:,KE+1) = EXP(-1.0E+00*tmp2(:,:,KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod and to3cnt.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact_part (6, Optical, cfc_tf, KE)
        to3cnt(:,:,KE+1) = to3cnt(:,:,KE+1)*cfc_tf(:,:,KE)
      endif

!-------------------------------------------------------------------


end subroutine optical_trans_funct_KE




!####################################################################
! <SUBROUTINE NAME="optical_trans_funct_diag">
!  <OVERVIEW>
!   Subroutine to compute diagnostic transmission function
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute diagnostic transmission function
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_trans_funct_diag (Atmos_input, contdg, to3dg, &
!                                  Optical)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data
!  </IN>
!  <OUT NAME="to3dg" TYPE="real">
!   Ozone continuum diagnostic transmission function
!  </OUT>
!  <OUT NAME="contdg" TYPE="real">
!   Diagnostic continuum transmission functions
!  </OUT> 
!  <INOUT NAME="Optical" TYPE="real">
!   Optical depth function
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_trans_funct_diag (Atmos_input, contdg, to3dg, &
                                     Optical)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

real, dimension (:,:,:),   intent(out)   :: to3dg                
real, dimension (:,:,:,:), intent(out)   :: contdg               
type(optical_path_type),   intent(inout) :: Optical
type(atmos_input_type),    intent(in)    :: Atmos_input

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Atmos_input
!
!   intent(inout) variables:
!
!     Optical
!
!   intent(out) variables:
!
!     to3dg
!     contdg
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (size(Atmos_input%pflux,1),             &
                       size(Atmos_input%pflux,2),             &
                       size(Atmos_input%pflux,3)-1) ::        &
                                                       pdfinv

      real, dimension (size(Atmos_input%pflux,1),          &
                       size(Atmos_input%pflux,2), &
                       size(Atmos_input%pflux,3)) ::  &
                                    press, pflux, ca, cb, csuba,  &
                                    csubb, ctmp2, ctmp3, delpr1, delpr2

!---------------------------------------------------------------------
!   local variables:
!
!      pdfinv
!      press 
!      pflux 
!      ca        
!      cb      
!      csuba 
!      csubb 
!      ctmp2 
!      ctmp3 
!      delpr1
!      delpr2
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
               'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
!    convert press and pflux to cgs.
!---------------------------------------------------------------------
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
      pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      pdfinv(:,:,ks:ke) = 1.0/(pflux(:,:,ks+1:ke+1) - pflux(:,:,ks:ke))
      delpr1(:,:,KS+1:KE)   = pdfinv (:,:,KS+1:KE)*  &
                              (press(:,:,KS+1:KE) - pflux(:,:,KS+1:KE)) 
      delpr2(:,:,KS+1:KE+1) = pdfinv(:,:,KS:KE)*   &
                              (pflux(:,:,KS+1:KE+1) - press(:,:,KS:KE)) 

!-----------------------------------------------------------------------
!    compute nearby-layer transmissivities for the o3 band and for the
!    one-band continuum band.  the sf function is used.
!    the method is the same as described for co2 in reference(4).
!-----------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'rsb' ) then
        ctmp2(:,:,KS+1:KE)  = Optical%cntval(:,:,KS+1:KE)*  &
                              delpr1(:,:,KS+1:KE) 
        ctmp3(:,:,KS+1:KE)  = Optical%cntval(:,:,KS:KE-1)*   &
                              delpr2(:,:,KS+1:KE) 
      endif
    
!-----------------------------------------------------------------------
!    compute sf2.
!    continuum band 1
!-----------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,4)*  &
                              delpr1(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,4)*  &
                              delpr2(:,:,KS+1:KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        csuba(:,:,KS+1:KE)  = betacm(12)*ctmp2(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = betacm(12)*ctmp3(:,:,KS+1:KE)
      endif
      ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +    &
                            csuba(:,:,KS+1:KE)*(0.166666E+00 -  &
                            csuba(:,:,KS+1:KE)*0.416666E-01))   
      cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +  &
                            csubb(:,:,KS+1:KE)*(0.166666E+00 - &
                            csubb(:,:,KS+1:KE)*0.416666E-01)) 
      contdg(:,:,KE+1,1)    = 1.0E+00 + cb (:,:,KE)
      contdg(:,:,KS+1:KE,1) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +  &
                              cb (:,:,KS+1:KE))

!--------------------------------------------------------------------
!    continuum band 2
!---------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,5)*   &
                              delpr1(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,5)*  &
                              delpr2(:,:,KS+1:KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        csuba(:,:,KS+1:KE)  = betacm(13)*ctmp2(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = betacm(13)*ctmp3(:,:,KS+1:KE)
      endif
      ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +  &
                            csuba(:,:,KS+1:KE)*(0.166666E+00 -   &
                            csuba(:,:,KS+1:KE)*0.416666E-01)) 
      cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +   &
                            csubb(:,:,KS+1:KE)*(0.166666E+00 -   &
                            csubb(:,:,KS+1:KE)*0.416666E-01)) 
      contdg(:,:,KE+1,2)    = 1.0E+00 + cb (:,:,KE)
      contdg(:,:,KS+1:KE,2) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +  &
                              cb (:,:,KS+1:KE))

!--------------------------------------------------------------------
!    continuum band 3
!--------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,7)*   &
                              delpr1(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,7)*  &
                              delpr2(:,:,KS+1:KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        csuba(:,:,KS+1:KE)  = betacm(15)*ctmp2(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = betacm(15)*ctmp3(:,:,KS+1:KE)
      endif
      ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +    &
                            csuba(:,:,KS+1:KE)*(0.166666E+00 -  &
                            csuba(:,:,KS+1:KE)*0.416666E-01)) 
      cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +   &
                            csubb(:,:,KS+1:KE)*(0.166666E+00 -  &
                            csubb(:,:,KS+1:KE)*0.416666E-01)) 
      contdg(:,:,KE+1,3)    = 1.0E+00 + cb (:,:,KE)
      contdg(:,:,KS+1:KE,3) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +   &
                              cb (:,:,KS+1:KE))

!--------------------------------------------------------------------
!    ozone band
!--------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,6)*   &
                              delpr1(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,6)*  &
                              delpr2(:,:,KS+1:KE)
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        csuba(:,:,KS+1:KE)  = betacm(14)*ctmp2(:,:,KS+1:KE)
        csubb(:,:,KS+1:KE)  = betacm(14)*ctmp3(:,:,KS+1:KE)
      endif
      ca   (:,:,KS+1:KE)  = csuba(:,:,KS+1:KE)*(-0.5E+00 +   &
                            csuba(:,:,KS+1:KE)*   &
                            (0.166666E+00 - csuba(:,:,KS+1:KE)*  &
                            0.416666E-01)) 
      cb   (:,:,KS+1:KE)  = csubb(:,:,KS+1:KE)*(-0.5E+00 +  &
                            csubb(:,:,KS+1:KE)*   &
                            (0.166666E+00 - csubb(:,:,KS+1:KE)*   &
                            0.416666E-01)) 
      to3dg (:,:,KE+1)    = 1.0E+00 + cb(:,:,KE)
      to3dg (:,:,KS+1:KE) = 1.0E+00 + 0.5E+00*(ca(:,:,KS+1:KE) +   &
                            cb(:,:,KS+1:KE))

!-------------------------------------------------------------------



end subroutine optical_trans_funct_diag


!###################################################################
! <SUBROUTINE NAME="get_totch2o">
!  <OVERVIEW>
!   Subroutine to compute self broadened temperature dependent
!   water vapor continuum
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute self broadened temperature dependent
!   water vapor continuum
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_totch2o (n, Optical, totch2o, dte1, ixoe1)
!  </TEMPLATE>
!  <IN NAME="n" TYPE="integer">
!   frequency band index
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth output
!  </INOUT>
!  <OUT NAME="totch2o" TYPE="real">
!   self broadened and temperature dependent continuum
!  </OUT>
!  <IN NAME="dte1" TYPE="real">
!   temperature step delta
!  </IN>
!  <IN NAME="ixoe1" TYPE="integer">
!   temperature index array
!  </IN>
! </SUBROUTINE>
!
subroutine get_totch2o (n, Optical, totch2o, dte1, ixoe1)

!------------------------------------------------------------------
!
!------------------------------------------------------------------

real, dimension(:,:,:),    intent(in)      :: dte1    
type(optical_path_type),   intent(inout)   :: Optical
integer, dimension(:,:,:), intent(in)      :: ixoe1   
real, dimension(:,:,:),    intent(out)     :: totch2o
integer,                   intent(in)      :: n

!-----------------------------------------------------------------
!   intent(in) variables:
!
!        dte1
!        ixoe1
!        n
!      
!   intent(inout) variables:
!        Optical
!
!   intent(out) variables:
!
!        totch2o
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!   local variables:

      real, dimension (size(Optical%tfac,1), size(Optical%tfac,2), &
                       size(Optical%tfac,3)) ::     &
                                                 radf, sh2o , tmpexp

      real               ::  fh2o0, sh2o0
      integer            ::  k, nu

!------------------------------------------------------------------
!   local variables:
!
!       radf
!       sh2o
!       tmpexp
!       fh2o0
!       sh2o0
!       k
!       nu
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
              'module has not been initialized', FATAL )
      endif
!--------------------------------------------------------------------
!    compute self-broadened temperature-dependent continuum coefficient
!    using the single coefficient -.013 for all frequencies in
!    the 160-560 cm-1 range. experiments with the mid-latitude
!    summer profile show errors of < .01 W/m**2 (in the net broadband
!    flux, 0-2200 cm-1) using this value. this value is used instead
!    of tmpfctrs at each frequency band.
!--------------------------------------------------------------------
      tmpexp(:,:,KS:KE) = EXP(-.013*Optical%tfac(:,:,KS:KE))

!--------------------------------------------------------------------
!    compute source function for frequency bands (ioffh2o+1 to ioffh2o
!    +nptch2o) at layer temperatures using table lookup.
!    note that ixoe1 can be used for temp index, and dte1 for deltat,
!    as the table extent for radf is the same as for the e1 tables
!    of the model.
!--------------------------------------------------------------------
      nu = n
      call looktab (radfunc, ixoe1, dte1, radf, KS, KE, nu+ioffh2o)
      sh2o0 = ssh2o_296(nu+ioffh2o)*sfac(nu+ioffh2o)

      do k=KS,KE 
        sh2o(:,:,k) = sh2o0*        tmpexp(:,:,k)
      end do
 
!--------------------------------------------------------------------
!    compute h2o self- and foreign- broadened continuum optical path,
!    summed from the top of the atmosphere through layer k.
!--------------------------------------------------------------------
      fh2o0 = sfh2o(nu+ioffh2o)*fscal(nu+ioffh2o)
      totch2o(:,:,1) = 0.0E+00
      do k = KS+1,KE+1
        totch2o(:,:,k) = Optical%wk(:,:,k-1)*1.0e-20*   &
                         (sh2o(:,:,k-1)*Optical%rh2os(:,:,k-1) +    &
                          fh2o0*Optical%rfrgn(:,:,k-1))* &
                          vvj(nu)*radf(:,:,k-1   )    +   &
                          totch2o(:,:,k-1)
      end do

!------------------------------------------------------------------

end subroutine get_totch2o



!#####################################################################
! <SUBROUTINE NAME="get_totch2obd">
!  <OVERVIEW>
!   Subroutine to compute self broadened temperature dependent
!   water vapor continuum
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute self broadened temperature dependent
!   water vapor continuum
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_totch2obd (n, Optical, totch2obd)
!  </TEMPLATE>
!  <IN NAME="n" TYPE="integer">
!   frequency band index
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth output
!  </INOUT>
!  <OUT NAME="totch2obd" TYPE="real">
!   self broadened and temperature dependent h2o continuum
!  </OUT>
! </SUBROUTINE>
!
subroutine get_totch2obd (n, Optical, totch2obd)

!------------------------------------------------------------------
!
!------------------------------------------------------------------

real, dimension(:,:,:), intent(out)     :: totch2obd
integer,                intent(in)      :: n
type(optical_path_type), intent(inout) :: Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      n
!
!   intent(inout) variable:
!
!      Optical
!
!   intent(out) variable:
!
!      totch2obd
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer            ::  k, nu

!--------------------------------------------------------------------
!  local variables:
!
!      k
!      nu
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
              'module has not been initialized', FATAL )
      endif
!---------------------------------------------------------------------
!    compute h2o self- and foreign- broadened continuum optical path 
!    for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!    the top of the atmosphere through layer k (totch2obd,
!    totch2obdinw, totch2obdwd).
!---------------------------------------------------------------------
      nu = n     
      totch2obd(:,:,1) = 0.0E+00
      do k = KS+1,KE+1
        totch2obd(:,:,k) = totch2obd(:,:,k-1) +   &
                           Optical%xch2obd(:,:,k-1,nu)
      end do

!--------------------------------------------------------------------
 
end subroutine get_totch2obd




!#####################################################################
! <SUBROUTINE NAME="get_totvo2">
!  <OVERVIEW>
!   Subroutine to compute continuum coefficients in band n
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute continuum coefficients in band n
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_totvo2 (n, Optical, totvo2_out) 
!  </TEMPLATE>
!  <IN NAME="n" TYPE="integer">
!   frequency band index
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Optical depth output
!  </INOUT>
!  <OUT NAME="totvo2_out" TYPE="real">
!   Continuum coefficients in band n
!  </OUT>
! </SUBROUTINE>
!
subroutine get_totvo2 (n, Optical, totvo2_out) 

!------------------------------------------------------------------
!
!------------------------------------------------------------------

integer,                 intent(in)       :: n
type(optical_path_type), intent(inout)    :: Optical
real, dimension(:,:,:),  intent(out)      :: totvo2_out

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      n
!
!   intent(inout) variable:
!
!      Optical
!
!   intent(out) variable:
!
!      totvo2_out
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
             'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------

      totvo2_out(:,:,:) = betacm(n)*Optical%totvo2(:,:,KS+1:KE+1)

end subroutine get_totvo2 



!####################################################################
! <SUBROUTINE NAME="optical_dealloc">
!  <OVERVIEW>
!   Subroutine to deallocate the array components of the 
!   optical_path_type input variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine deallocates the array components of the 
!   optical_path_type input variable. Dependent on the namelist
!   options chosen, some of the arrays may or may nothave been
!   allocated.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_dealloc (Optical)            
!  </TEMPLATE>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Derived type variable containing information related to
!   the computation of optical depth associated with 
!   different atmospheric constituents. 
!  </INOUT>
! </SUBROUTINE>
!

subroutine optical_dealloc (Optical, including_aerosols)  

!-------------------------------------------------------------------
!    optical_dealloc deallocates the array components of the 
!    optical_path_type input variable.
!--------------------------------------------------------------------

type(optical_path_type), intent(inout) :: Optical
logical,                   intent(in)            :: including_aerosols  

!--------------------------------------------------------------------
! intent(inout) variables:
!
!    Optical       optical_path_type variable containing fields used
!                  in the calculation of optical paths for various
!                  atmospheric constituents
! 
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    deallocate the array elements of Optical.
!--------------------------------------------------------------------
       deallocate (Optical%empl1          )
       deallocate (Optical%empl2          )
       deallocate (Optical%var1           )
       deallocate (Optical%var2           )
       deallocate (Optical%avephi         )
       deallocate (Optical%totphi         )
       deallocate (Optical%emx1           )
       deallocate (Optical%emx2           )

 
       if (NBTRGE > 0) then
         deallocate (Optical%avephif        )
         deallocate (Optical%emx1f          )
         deallocate (Optical%emx2f          )
         deallocate (Optical%empl1f         )
         deallocate (Optical%empl2f         )
         deallocate (Optical%vrpfh2o        )
         deallocate (Optical%tphfh2o         )
       endif
 
       if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
           trim(Lw_control%continuum_form) == 'ckd2.4' ) then
         deallocate (Optical%xch2obd        )
         deallocate (Optical%totch2obdwd    )
         deallocate (Optical%xch2obdwd      )
       else if (trim(Lw_control%continuum_form) == 'rsb' ) then
         deallocate (Optical%cntval         )
         deallocate (Optical%totvo2         )
       endif
 
       deallocate (Optical%toto3          )
       deallocate (Optical%tphio3         )
       deallocate (Optical%var3           )
       deallocate (Optical%var4           )
       deallocate (Optical%wk             )
       deallocate (Optical%rh2os          )
       deallocate (Optical%rfrgn          )
       deallocate (Optical%tfac           )

       if (Lw_control%do_cfc) then
         deallocate (Optical%totf11         )
         deallocate (Optical%totf12         )
         deallocate (Optical%totf113         )
         deallocate (Optical%totf22         )
       endif

       if (including_aerosols) then
         deallocate (Optical%totaerooptdep  )
         deallocate (Optical%aerooptdep_KE_15  )
       endif

!-------------------------------------------------------------------


end subroutine optical_dealloc



!####################################################################
! <SUBROUTINE NAME="optical_path_end">
!  <OVERVIEW>
!   optical_path_end is the destructor for optical_path_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   optical_path_end is the destructor for optical_path_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_depth_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine optical_path_end

!--------------------------------------------------------------------
!    optical_path_end is the destructor for optical_path_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'optical_path_mod',  &
             'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------
!    mark the module as uninitialized.
!-----------------------------------------------------------------
      module_is_initialized = .false.

!------------------------------------------------------------------



end subroutine optical_path_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
                                   
                                  
!###################################################################
! <SUBROUTINE NAME="optical_ckd_init">
!  <OVERVIEW>
!   Subroutine to initialize water vapor self and foreign broadened
!   continuum coefficients. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   Idckdh2o reads ckd2.1 self and foreign-broadened h2o continuum
!   coefficients, corrections, and coefficients for temperature
!   dependence of the self-continuum. these are tabulated at 10
!   cm-1 intervals from 0 to 20000 cm-1
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_ckd_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine optical_ckd_init

!------------------------------------------------------------------
!    optical_ckd_init reads ckd2.1 or ckd2.4 self and foreign-broadened
!    h2o continuum coefficients, corrections, and coefficients for
!    temperature dependence of the self-continuum. these are tabulated
!    at 10 cm-1 intervals from 0 to 20000 cm-1.
!    (the above information is as of 2/12/96).
!
!    references:
!
!    (1) clough, s. a.  et al. "line shape and the water vapor
!        continuum," atmospheric research, 23 (1989) 229-241.
!
!
!    author: m. d. schwarzkopf
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:


!-------------------------------------------------------------------
!   data from former block data bs260 for self-broadened continuum
!    at 260K, band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!---------------------------------------------------------------------
      real    ::  v1sh2o_260, v2sh2o_260, dvsh2o_260,    &
                  ssh2o_260(2000)
      integer ::  nptsh2o_260

!--------------------------------------------------------------------
!        tktab and vjtab are the respective temperature and frequency
!    points at which tabulations occurred.
!---------------------------------------------------------------------
      real   ::   tktab(40),  vjtab(300)

!---------------------------------------------------------------------
      integer  :: inrad, k, j, ihih2o

!--------------------------------------------------------------------
!   local variables:
!
!      v1sh2o_260
!      v2sh2o_260
!      dvsh2o_260
!      ssh2o_260
!      nptsh2o_260
!      tktab
!      vjtab
!      inrad
!      k,j
!      ihih2o
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    call routine to allocate radfunc table
!---------------------------------------------------------------------
      call table_alloc (radfunc, 40, 300)

!--------------------------------------------------------------------
!    read h2o (original) data
!    data are at frequencies 5 - 19995 cm-1, at 10 cm-1 intervals
!-------------------------------------------------------------------
      inrad = open_namelist_file ('INPUT/h2ockd2.1_data')
      read (inrad,9001) v1sh2o_296, v2sh2o_296, dvsh2o_296,  &
                        nptsh2o_296
      read (inrad,9002) (ssh2o_296(k),k=1,2000)
      read (inrad,9001) v1sh2o_260, v2sh2o_260, dvsh2o_260,   &
                        nptsh2o_260
      read (inrad,9002) (ssh2o_260(k),k=1,2000)
      read (inrad,9001) v1fh2o, v2fh2o, dvfh2o, nptfh2o
      read (inrad,9002) (sfh2o(k),k=1,2000)
9001  format (3f12.1,i8)
9002  format (5e14.5)
 
      call close_file (inrad)

!--------------------------------------------------------------------
!    read h2o corrected data
!--------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1') then
        inrad = open_namelist_file ('INPUT/h2ockd2.1_corrdata')
      else if (trim(Lw_control%continuum_form) == 'ckd2.4') then
        inrad = open_namelist_file ('INPUT/h2ockd2.4_corrdata')
      endif
      read (inrad,9007) (sfac(k),k=1,2000)
      read (inrad,9007) (fscal(k),k=1,2000)
      read (inrad,9007) (tmpfctrs(k),k=1,2000)
9007  format (5e13.6)
 
      call close_file (inrad)

!--------------------------------------------------------------------
!    read radfn data
!--------------------------------------------------------------------
      inrad = open_namelist_file ('INPUT/radfn_5-2995_100-490k')
      read (inrad,9000) ((radfunc%vae(k,j),radfunc%td(k,j),k=1,40), &
                                                           j=1,300)
9000  format (8f14.6)
      call close_file (inrad)

!---------------------------------------------------------------------
      do k=1,40
        tktab(k) = 100. + 10.*(k-1)
      end do
      do j=1,300
        vjtab(j) = 5. + 10.*(j-1)
      end do
 
!--------------------------------------------------------------------
!    compute range to use in datasets for actual frequency intervals
!    used in model.
!
!    freqlo = 160.
!    freqhi = 560.
!
!    define initial offset and number of data points to use
!    for the 3 h2o continua over the frequency range of the
!    calculations (freqlo,freqhi). note: we assume no interpolation
!    is needed. if interp. was required, these limits would be
!    expanded. values are put into commons in include file tab.h
!    for transmission into Optical_ckd2.1.F.
!
!    ioff is the offset from the absorption tables (starting at 5)
!    needed for proper freq computations. first index used then
!    is (ioff+1). for calculations with the first band beginning
!    at 160 cm-1, this number is 16, and the index number of the
!    band ending at 560 cm-1 is 56.
!-----------------------------------------------------------------------
      ioffh2o = 16

!---------------------------------------------------------------------
!    the final index number used in the calculation is (ihi)
!--------------------------------------------------------------------
      ihih2o  = 56

!--------------------------------------------------------------------
!    nptc is the number of frequency points used in the calculation.
!    ( = ihi - (ioff+1) + 1)
!---------------------------------------------------------------------
      nptch2o = ihih2o - ioffh2o

!---------------------------------------------------------------------
!    vvj are the frequencies for calculation of h2o coefficients. by
!    assumption, no other frequencies are used.
!----------------------------------------------------------------------
      do j=1,nptch2o
        vvj(j) = v1sh2o_296 + dvsh2o_296*float(j+ioffh2o-1)
      end do

!---------------------------------------------------------------------
!    compute h2o coefficients averaged over the broad bands used
!    in the 560 -1200 cm-1 range. until the frequency bands are read
!    in, I will re-list them here, rather than use rnddta.H variables
!    (where they are stored).
!
!    the required wide bands are:
!        560-630 cm-1
!        630-700   (assuming 3 bands in 15um complex)
!        700-800
!        560-800   (1 band for entire complex)
!        800-900
!        900-990
!        990-1070
!        1070-1200
!        800-900,1070-1200   (until this band is broken into 2)
!    we assume that, for best accuracy:
!    the quantity required is <svj> and <fvj) where angle brackets are
!    averages over frequency, s and f are self- and foreign coeff-
!    icients, including corrections, and vj is frequency (from vjtab).
!    notations for special bands attempt similarity with that
!    previously used in the radiation code.
!    we also assume that one value may be used (at all altitudes)
!    for the radiation correction term radfn, in each frequency band.
!    the values used below result from experimentation.
!---------------------------------------------------------------------
      svj = 0.0
      fvj = 0.0
      svjwd = 0.0
      fvjwd = 0.0
      svjinw = 0.0
      fvjinw = 0.0

!--------------------------------------------------------------------
!    560-630 band:
!--------------------------------------------------------------------
      do j=57,63
        svj(1) = svj(1) + vjtab(j)*ssh2o_296(j)*sfac(j)/7.
        fvj(1) = fvj(1) + vjtab(j)*sfh2o(j)*fscal(j)/7.
      end do
      radfnbd(1) = 0.90

!--------------------------------------------------------------------
!    630-700 band:
!--------------------------------------------------------------------
      do j=64,70
        svj(2) = svj(2) + vjtab(j)*ssh2o_296(j)*sfac(j)/7.
        fvj(2) = fvj(2) + vjtab(j)*sfh2o(j)*fscal(j)/7.
      end do
      radfnbd(2) = 0.92

!--------------------------------------------------------------------
!    700-800 band:
!--------------------------------------------------------------------
      do j=71,80
        svj(3) = svj(3) + vjtab(j)*ssh2o_296(j)*sfac(j)/10.
        fvj(3) = fvj(3) + vjtab(j)*sfh2o(j)*fscal(j)/10.
      end do
      radfnbd(3) = 0.95
!--------------------------------------------------------------------
!    800-900 band:
!--------------------------------------------------------------------
      do j=81,90
        svj(4) = svj(4) + vjtab(j)*ssh2o_296(j)*sfac(j)/10.
        fvj(4) = fvj(4) + vjtab(j)*sfh2o(j)*fscal(j)/10.
      end do
      radfnbd(4) = 0.97

!--------------------------------------------------------------------
!    900-990 band:
!--------------------------------------------------------------------
      do j=91,99
        svj(5) = svj(5) + vjtab(j)*ssh2o_296(j)*sfac(j)/9.
        fvj(5) = fvj(5) + vjtab(j)*sfh2o(j)*fscal(j)/9.
      end do
      radfnbd(5) = 0.98

!--------------------------------------------------------------------
!    990-1070 band:
!--------------------------------------------------------------------
      do j=100,107
        svj(6) = svj(6) + vjtab(j)*ssh2o_296(j)*sfac(j)/8.
        fvj(6) = fvj(6) + vjtab(j)*sfh2o(j)*fscal(j)/8.
      end do
      radfnbd(6) = 0.99

!--------------------------------------------------------------------
!    1070-1200 band:
!--------------------------------------------------------------------
      do j=108,120
        svj(7) = svj(7) + vjtab(j)*ssh2o_296(j)*sfac(j)/13.
        fvj(7) = fvj(7) + vjtab(j)*sfh2o(j)*fscal(j)/13.
      end do
      radfnbd(7) = 0.992

!--------------------------------------------------------------------
!    560-800 combined band:
!-------------------------------------------------------------------
      do j=57,80
        svjwd = svjwd + vjtab(j)*ssh2o_296(j)*sfac(j)/24.
        fvjwd = fvjwd + vjtab(j)*sfh2o(j)*fscal(j)/24.
      end do
      radfnbdwd = 0.92

!--------------------------------------------------------------------
!    800-990,1070-1200 combined band:
!--------------------------------------------------------------------
      do j=81,99
        svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac(j)/22.
        fvjinw = fvjinw + vjtab(j)*sfh2o(j)*fscal(j)/32.
      end do
      do j=108,120
        svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac(j)/32.
        fvjinw = fvjinw + vjtab(j)*sfh2o(j)*fscal(j)/32.
      end do
      radfnbdinw = 0.98

!--------------------------------------------------------------------


end subroutine optical_ckd_init




!###################################################################
! <SUBROUTINE NAME="optical_path_ckd">
!  <OVERVIEW>
!   Subroutine to compute water vapor self and foreign broadened 
!   continuum optical paths
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute water vapor self and foreign broadened 
!   continuum optical paths over the frequency range specified by
!    ioffh2o and nptch2o using the ckd algorithm, modified for 
!    the gcm parameterization.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_path_ckd (atmden, press, temp, rh2o, Optical)
!  </TEMPLATE>
!  <IN NAME="atmden" TYPE="real">
!   Atmospheric density profile
!  </IN>
!  <IN NAME="press" TYPE="real">
!   The pressure coordinate array
!  </IN>
!  <IN NAME="temp" TYPE="real">
!   Temperature
!  </IN> 
!  <IN NAME="rh2o" TYPE="real">
!   mass mixing ratio of h2o at model data levels
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   water vapor continuum optical path otuput
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_path_ckd (atmden, press, temp, rh2o, Optical) 

!------------------------------------------------------------------
!    subroutine optical_ckd computes h2o continuum optical paths
!    (self + foreign) over the frequency range specified by
!    ioffh2o and nptch2o using the ckd algorithm, modified for 
!    the gcm parameterization.
!    (this routine is previously called contnm.F)
!------------------------------------------------------------------

real, dimension (:,:,:), intent(in)       :: atmden, press, temp, rh2o
type(optical_path_type), intent(inout)    :: Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      atmden
!      press
!      temp
!      rh2o
!
!   intent(inout) variable:
!
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(press,1), size(press,2), &
                       size(press,3)) ::                totch2obdinw

      real, dimension (size(press,1), size(press,2), &
                       size(press,3)-1) ::       &
                                    xch2obdinw, tmpexp, rvh2o, rhoave

      real                    ::  t0 = 296.0
      integer                 ::  k, nu
      integer      :: israd, ierad, jsrad, jerad

!---------------------------------------------------------------------
!  local variables:
!
!      totch2obdinw
!      xch2obdinw
!      tmpexp
!      rvh2o
!      rhoave
!      t0
!      n,k
!      nu
!
!--------------------------------------------------------------------
      israd = 1
      ierad = size(press,1)
      jsrad = 1
      jerad = size(press,2)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      allocate (Optical%xch2obd    (ISRAD:IERAD, JSRAD:JERAD,    &
                                                          KS:KE  , 7) )
      allocate (Optical%totch2obdwd(ISRAD:IERAD, JSRAD:JERAD,    &
                                                          KS:KE+1   ) )
      allocate (Optical%xch2obdwd  (ISRAD:IERAD, JSRAD:JERAD,    &
                                                          KS:KE     ) )
      Optical%xch2obd  = 0.                                           
      Optical%totch2obdwd = 0.                                        
      Optical%xch2obdwd  = 0.      

!--------------------------------------------------------------------
!    define the volume mixing ratio of h2o
!---------------------------------------------------------------------
      rvh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)/d622

!---------------------------------------------------------------------
!    define input arguments to optical_ckd
!    wk is column density (molec/cm2) of water vapor
!    rfrgn is partial pressure (Amagat) at 296K from N2+O2+Ar
!    rh2os is partial pressure (Amagat) at 296K from water vapor
!-------------------------------------------------------------------
      Optical%wk(:,:,KS:KE) =  rvh2o(:,:,KS:KE)*avogno/wtmair*   &
                               atmden(:,:,KS:KE)/   &
                               (1.0 + rvh2o(:,:,KS:KE))
      rhoave(:,:,KS:KE) = (press(:,:,KS:KE)/pstd)*   &
                          (tfreeze/temp(:,:,KS:KE))
      Optical%rfrgn(:,:,KS:KE) =  rhoave(:,:,KS:KE)*(t0/tfreeze)/  &
                                  (1.0 + rvh2o(:,:,KS:KE))
      Optical%rh2os(:,:,KS:KE) = Optical%rfrgn(:,:,KS:KE)*   &
                                 rvh2o(:,:,KS:KE)
      Optical%tfac(:,:,KS:KE) = temp(:,:,KS:KE) - t0

!--------------------------------------------------------------------
!    compute self-broadened temperature-dependent continuum coefficient
!    using the single coefficient -.020 for all frequencies in
!    the 560-1200 cm-1 range. experiments with the mid-latitude
!    summer profile show errors of < .01 W/m**2 (in the net broadband
!    flux, 0-2200 cm-1) using this value. this value is used instead
!    of tmpfctrs at each frequency band.
!-------------------------------------------------------------------
      tmpexp(:,:,KS:KE) = EXP(-.020*Optical%tfac(:,:,KS:KE))
 
!-------------------------------------------------------------------
!    compute h2o self- and foreign- broadened continuum optical path 
!    for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!    the top of the atmosphere through layer k (totch2obd,
!    totch2obdinw, totch2obdwd).
!--------------------------------------------------------------------
      do nu = 1,7
        do k = KS,KE 
          Optical%xch2obd(:,:,k,nu) = Optical%wk(:,:,k)*1.0e-20*   &
                                      (svj(nu)*Optical%rh2os(:,:,k)*&
                                      tmpexp(:,:,k) + fvj(nu)*   &
                                      Optical%rfrgn(:,:,k))*radfnbd(nu)
        end do
      end do
 
      do k = KS,KE 
        xch2obdinw(:,:,k) = Optical%wk(:,:,k)*1.0e-20*(svjinw*  &
                            Optical%rh2os(:,:,k)* tmpexp(:,:,k) +   &
                            fvjinw*Optical%rfrgn(:,:,k))*radfnbdinw
        Optical%xch2obdwd(:,:,k) = Optical%wk(:,:,k)*1.0e-20*   &
                                   (svjwd*Optical%rh2os(:,:,k)* &
                                   tmpexp(:,:,k) + fvjwd*  &
                                   Optical%rfrgn(:,:,k))*radfnbdwd
      end do
 
      totch2obdinw(:,:,1) = 0.0E+00
      Optical%totch2obdwd(:,:,1) = 0.0E+00
      do k = KS+1,KE+1
        totch2obdinw(:,:,k) = totch2obdinw(:,:,k-1) +    &
                              xch2obdinw(:,:,k-1)
        Optical%totch2obdwd(:,:,k) = Optical%totch2obdwd(:,:,k-1) + &
                                     Optical%xch2obdwd(:,:,k-1)
      end do

!----------------------------------------------------------------------

 
end subroutine optical_path_ckd
 


!################################################################## 
! <SUBROUTINE NAME="optical_o3">
!  <OVERVIEW>
!   Subroutine to compute optical paths for o3.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute optical paths for o3.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_o3 (atmden, qo3, vv, Optical)
!  </TEMPLATE>
!  <IN NAME="atmden" TYPE="real">
!   Atmospheric density profile
!  </IN>
!  <IN NAME="qo3" TYPE="real">
!   mass mixing ratio of o3 at model data levels
!  </IN>
!  <IN NAME="vv" TYPE="real">
!   Ozone volume mixing atio
!  </IN> 
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   ozone optical path otuput
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_o3 (atmden, qo3, vv, Optical)

!------------------------------------------------------------------
!    optical_o3 computes optical paths for o3.
!------------------------------------------------------------------

real, dimension(:,:,:),  intent(in)    ::  atmden, qo3, vv
type(optical_path_type), intent(inout) ::  Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!     atmden
!     qo3     mass mixing ratio of o3 at model data levels.
!     vv
!
!   intent(inout) variable:
!
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer  ::    k    ! do-loop index
      integer      :: israd, ierad, jsrad, jerad

!---------------------------------------------------------------------
      israd = 1
      ierad = size(qo3,1)
      jsrad = 1
      jerad = size(qo3,2)

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      allocate (Optical%toto3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (Optical%tphio3(ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (Optical%var3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))
      allocate (Optical%var4  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))
      Optical%toto3  = 0.
      Optical%tphio3 = 0.
      Optical%var3  = 0.
      Optical%var4  = 0.                                        

!-----------------------------------------------------------------------
!    compute optical paths for o3, using the diffusivity 
!    approximation 1.66 for the angular integration.  obtain 
!    unweighted values var3 and weighted values  var4.
!    the quantities  0.003 (.003) appearing in the
!    var4 expression are the approximate voigt corrections
!    for o3.
!---------------------------------------------------------------------  
      Optical%var3(:,:,KS:KE) = atmden(:,:,KS:KE)*qo3(:,:,KS:KE)*diffac
      Optical%var4(:,:,KS:KE) = Optical%var3(:,:,KS:KE)*    &
                                (vv(:,:,KS:KE) + 3.0E-03)

!----------------------------------------------------------------------
!    compute summed optical paths for o3.
!----------------------------------------------------------------------
      Optical%toto3 (:,:,KS) = 0.0E+00
      Optical%tphio3(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%toto3 (:,:,k) = Optical%toto3 (:,:,k-1) +    &
                                Optical%var3  (:,:,k-1) 
        Optical%tphio3(:,:,k) = Optical%tphio3(:,:,k-1) +    &
                                Optical%var4  (:,:,k-1) 
      end do

!----------------------------------------------------------------------


end subroutine optical_o3




!#####################################################################
! <SUBROUTINE NAME="optical_rbts">
!  <OVERVIEW>
!   Subroutine to compute optical paths for h2o rbts continuum
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute optical paths for h2o rbts continuum
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_rbts (temp, rh2o, Optical) 
!  </TEMPLATE>
!  <IN NAME="temp" TYPE="real">
!   temperature profile used in continuum calculation
!  </IN>
!  <IN NAME="rh2o" TYPE="real">
!   mass mixing ratio of h2o at model data levels
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   water vapor robert continuum optical path
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_rbts (temp, rh2o, Optical) 

!------------------------------------------------------------------
!    optical_rbts computes optical paths for h2o rbts comtinuum.
!------------------------------------------------------------------

real, dimension(:,:,:),  intent(in)    :: temp, rh2o
type(optical_path_type), intent(inout) :: Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      temp
!      rh2o
!
!   intent(inout) variable:
!
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension(size(temp,1), size(temp,2), &
                                     size(temp,3)) :: texpsl
      integer     :: k
      integer      :: israd, ierad, jsrad, jerad

!--------------------------------------------------------------------
!  local variables:
!
!      texpsl
!      i,k
!
!----------------------------------------------------------------------
      israd = 1
      ierad = size(temp,1)
      jsrad = 1
      jerad = size(temp,2)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      allocate (Optical%cntval(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1   ))
      allocate (Optical%totvo2(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1   ))
      Optical%cntval = 0.                                         
      Optical%totvo2 = 0.                                        

!----------------------------------------------------------------------
!    compute argument for constant temperature coefficient (this is 
!    1.800E+03/(1.0E+00/temp - 1.0E+00/2.960E+02)).
!---------------------------------------------------------------------- 
      texpsl(:,:,KS:KE+1) = EXP(1.800E+03/temp(:,:,KS:KE+1) -   &
                                6.081081081E+00) 

!----------------------------------------------------------------------
!    compute optical path for the h2o continuum, using roberts 
!    coefficients betinw, and temperature correction texpsl. 
!    the diffusivity approximation (which cancels out in this
!    expression) is assumed to be 1.66.  the use of the diffusivity
!    factor has been shown to be a significant source of error in the
!    continuum calculations, however, the time penalty of an angular
!    integration is severe.
!---------------------------------------------------------------------  
      Optical%cntval(:,:,KS:KE) = texpsl(:,:,KS:KE)*rh2o(:,:,KS:KE)*   &
                                  Optical%var2(:,:,KS:KE)/   &
                                  (rh2o(:,:,KS:KE) + d622   )

!----------------------------------------------------------------------
!    compute summed optical paths for h2o roberts continuum.
!----------------------------------------------------------------------
      Optical%totvo2(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%totvo2(:,:,k) = Optical%totvo2(:,:,k-1) +   &
                                Optical%cntval(:,:,k-1) 
      end do

!----------------------------------------------------------------------



end subroutine optical_rbts



!####################################################################
! <SUBROUTINE NAME="optical_h2o">
!  <OVERVIEW>
!   Subroutine to compute water vapor optical paths
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute water vapor optical paths
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_h2o (pflux, atmden, vv, press, temp, rh2o, tflux, &
!                     Optical) 
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   pressure at flux levels of model
!  </IN>
!  <IN NAME="atmden" TYPE="real">
!   Atmospheric density profile
!  </IN>
!  <IN NAME="vv" TYPE="real">
!   volume mixing ratio of h2o at model data levels
!  </IN>
!  <IN NAME="press" TYPE="real">
!   The pressure coordinate array
!  </IN>
!  <IN NAME="temp" TYPE="real">
!   Temperature at data levels of model
!  </IN> 
!  <IN NAME="rh2o" TYPE="real">
!   mass mixing ratio of h2o at model data levels
!  </IN>
!  <IN NAME="tflux" TYPE="real">
!   Temperature at flux levels of model
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   water vapor optical path otuput
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_h2o (pflux, atmden, vv, press, temp, rh2o, tflux, &
                        Optical) 

!----------------------------------------------------------------------
!    optical_h2o computes optical paths for h2o.
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(in)    ::  pflux, atmden, vv, press, &
                                           temp, rh2o, tflux
type(optical_path_type), intent(inout) ::  Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!     pflux     pressure at flux levels of model.
!     atmden
!     vv
!     press     pressure at data levels of model.
!     temp      temperature at data levels of model. 
!     rh2o      mass mixing ratio of h2o at model data levels 
!     tflux
!
!   intent(inout) variable:
!
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(pflux,1), size(pflux,2), &
                       size(pflux,3)) ::        &
                                             tpl1, tpl2, &
                                             qh2o, tdif, tdif2
      integer    ::  m, k
      integer      :: israd, ierad, jsrad, jerad

!--------------------------------------------------------------------
!  local variables:
!
!      tpl1
!      tpl2
!      qh2o       h2o mass mixing ratio, multiplied by the diffusivity
!                 factor diffac.
!      tdif
!      tdif2
!      m,k
!
!-----------------------------------------------------------------------

      israd = 1
      ierad = size(pflux,1)
      jsrad = 1
      jerad = size(pflux,2)
!-------------------------------------------------------------------- 
!    compute mean temperature in the "nearby layer" between a flux
!    level and the first data level below the flux level (tpl1) or the
!    first data level above the flux level (tpl2)
!---------------------------------------------------------------------
      tpl1(:,:,KS   )         = temp(:,:,KE   )
      tpl1(:,:,KS   +1:KE   ) = tflux(:,:,KS   +1:KE   )
      tpl1(:,:,KE   +1)       = 0.5E+00*(tflux(:,:,KE   +1) +   &
                                temp(:,:,KE   ))
      tpl2(:,:,KS   +1:KE   ) = tflux(:,:,KS   +1:KE   )
      tpl2(:,:,KE   +1)       = 0.5E+00*(tflux(:,:,KE   ) +    &
                                temp(:,:,KE   ))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      allocate (Optical%empl1  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1   ))
      allocate (Optical%empl2  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1   ))
      allocate (Optical%totphi (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1   ))
      allocate (Optical%var1   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE     ))
      allocate (Optical%var2   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE     ))
      allocate (Optical%emx1   (ISRAD:IERAD, JSRAD:JERAD              ))
      allocate (Optical%emx2   (ISRAD:IERAD, JSRAD:JERAD              ))
      Optical%empl1   = 0.
      Optical%empl2  =0.
      Optical%totphi  = 0.
      Optical%var1   = 0.
      Optical%var2   = 0.
      Optical%emx1   = 0.
      Optical%emx2   = 0.

!----------------------------------------------------------------------
!    compute optical paths for h2o, using the diffusivity 
!    approximation 1.66 for the angular integration.  obtain 
!    unweighted values var1, and weighted values var2.
!    the quantities 0.0003 (.0003) appearing in the
!    var2 expressions are the approximate voigt corrections
!    for h2o.  vv is the layer-mean pressure (in 
!    atmosphere), which is not the same as the level pressure press.
!---------------------------------------------------------------------  
      qh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)*diffac
      Optical%var1(:,:,KS:KE) = atmden(:,:,KS:KE)*qh2o(:,:,KS:KE)
      Optical%var2(:,:,KS:KE) = Optical%var1(:,:,KS:KE)*   &
                                (vv(:,:,KS:KE) + 3.0E-04)

!----------------------------------------------------------------------
!    compute summed optical paths for h2o.
!----------------------------------------------------------------------
      Optical%totphi(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%totphi(:,:,k) = Optical%totphi(:,:,k-1) +   &
                                Optical%var2  (:,:,k-1) 
      end do

!----------------------------------------------------------------------
!    emx1 is the additional pressure-scaled mass from press(KE) to 
!    pflux(KE).  it is used in nearby layer and emiss calculations.
!    emx2 is the additional pressure-scaled mass from press(KE) to 
!    pflux(KE+1).  it is used in calculations between flux levels k
!    and KE+1.
!----------------------------------------------------------------------
      Optical%emx1(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(press(:,:,KE) - &
                          pflux(:,:,KE))/(1.0E+02*GRAV*pstd)
      Optical%emx2(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(pflux(:,:,KE+1) -&
                          press(:,:,KE))/(1.0E+02*GRAV*pstd)

!----------------------------------------------------------------------
!    empl is the pressure scaled mass from pflux(k) to press(k) or to 
!    press(k+1).
!----------------------------------------------------------------------
      Optical%empl1(:,:,KS) = Optical%var2(:,:,KE)
      Optical%empl1(:,:,KS+1:KE+1) = qh2o(:,:,KS:KE)*    &
                                     pflux(:,:,KS+1:KE+1)*   &
                                     (pflux(:,:,KS+1:KE+1) -   &
                                      press(:,:,KS:KE))/   &
                                      (1.0E+02*GRAV*pstd)
      Optical%empl2(:,:,KS+1:KE) =    &
                 qh2o(:,:,KS+1:KE)*pflux(:,:,KS+1:KE)*   &
                 (press(:,:,KS+1:KE) - pflux(:,:,KS+1:KE))/  &
                 (1.0E+02*GRAV*pstd)
      Optical%empl2(:,:,KE+1) = Optical%empl2(:,:,KE) 

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( Optical%empl1f (ISRAD:IERAD , JSRAD:JERAD ,    & 
                                                  KS:KE+1,  NBTRGE ) ) 
        allocate ( Optical%empl2f (ISRAD:IERAD , JSRAD:JERAD ,     &  
                                                  KS:KE+1,  NBTRGE ) ) 
        allocate ( Optical%tphfh2o(ISRAD:IERAD , JSRAD:JERAD ,     & 
                                                  KS:KE+1,  NBTRGE ) ) 
        allocate ( Optical%vrpfh2o(ISRAD:IERAD , JSRAD:JERAD ,    &
                                                  KS:KE+1,  NBTRGE ) )
        allocate ( Optical%emx1f  (ISRAD:IERAD , JSRAD:JERAD ,   &
                                                            NBTRGE ) )
        allocate ( Optical%emx2f  (ISRAD:IERAD , JSRAD:JERAD ,   &
                                                            NBTRGE ) )
        Optical%empl1f  = 0.
        Optical%empl2f  = 0.
        Optical%tphfh2o  = 0.
        Optical%vrpfh2o = 0.
        Optical%emx1f   = 0.
        Optical%emx2f  = 0.                               

        if (tmp_dpndnt_h2o_lines) then
!----------------------------------------------------------------------
!    compute h2o optical paths for use in the 1200-1400 cm-1 range if
!    temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
          tdif(:,:,KS:KE) = temp(:,:,KS:KE)-2.5E+02

          do m=1,NBTRGE
            Optical%vrpfh2o(:,:,KS:KE,m) = Optical%var2(:,:,KS:KE)*   &
                                           EXP(csfah2o(1,m)*   &
                                               (tdif(:,:,KS:KE)) +   &
                                               csfah2o(2,m)*   &
                                               (tdif(:,:,KS:KE))**2 )
          end do
          do m=1,NBTRGE
            Optical%tphfh2o(:,:,KS,m) = 0.0E+00
            do k=KS+1,KE+1
              Optical%tphfh2o(:,:,k,m) = Optical%tphfh2o(:,:,k-1,m) +  &
                                         Optical%vrpfh2o(:,:,k-1,m)
            end do
          end do

          tdif2(:,:,KS+1:KE+1) = tpl2(:,:,KS+1:KE+1)-2.5E+02
          tdif (:,:,KS+1:KE+1) = tpl1(:,:,KS+1:KE+1)-2.5E+02

!---------------------------------------------------------------------
!    compute this additional mass, for use in the 1200-1400 cm-1 range,
!    if temperature dependence of line intensities is accounted for.
!--------------------------------------------------------------------
          do m=1,NBTRGE
            Optical%emx1f(:,:,m) = Optical%emx1(:,:) *    &
                                   EXP(csfah2o(1,m)*(tdif2(:,:,KE+1)) +&
                                     csfah2o(2,m)*(tdif2(:,:,KE+1))**2 )
            Optical%emx2f(:,:,m) = Optical%emx2(:,:) *    &
                                 EXP(csfah2o(1,m)*(tdif (:,:,KE+1)) + &
                                     csfah2o(2,m)*(tdif (:,:,KE+1))**2 )
          end do

!----------------------------------------------------------------------
!    compute this additional mass, for use in the 1200-1400 cm-1 range,
!    if temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
          do m=1,NBTRGE
            Optical%empl1f(:,:,KS+1:KE+1,m) =     &
                                        Optical%empl1(:,:,KS+1:KE+1)*&
                                        EXP(csfah2o(1,m)*   &
                                            (tdif(:,:,KS+1:KE+1)) + &
                                            csfah2o(2,m)*   &
                                            (tdif(:,:,KS+1:KE+1))**2 )
            Optical%empl2f(:,:,KS+1:KE,m) = Optical%empl2(:,:,KS+1:KE)*&
                                          EXP(csfah2o(1,m)*  &
                                              (tdif2(:,:,KS+1:KE)) +   &
                                              csfah2o(2,m)*  &
                                              (tdif2(:,:,KS+1:KE))**2 )
            Optical%empl1f(:,:,KS ,m) = Optical%vrpfh2o(:,:,KE,m)
            Optical%empl2f(:,:,KE+1,m) = Optical%empl2f(:,:,KE,m)
          end do
        else
          do m=1,NBTRGE
            Optical%empl1f(:,:,ks+1:ke+1,m) =   &
                                           Optical%empl1(:,:,ks+1:ke+1)
            Optical%empl2f(:,:,ks+1:ke,m) = Optical%empl2(:,:,ks+1:ke)
            Optical%emx1f(:,:,m)   = Optical%emx1(:,:)
            Optical%emx2f(:,:,m)  = Optical%emx2(:,:)              
            Optical%tphfh2o(:,:,:,m) = Optical%totphi (:,:,:)
            Optical%vrpfh2o(:,:,KE,m) = Optical%var2(:,:,KE)
            Optical%vrpfh2o(:,:,KS:ke,m) = Optical%var2(:,:,KS:Ke)
            Optical%empl1f(:,:,KS,m) = Optical%vrpfh2o(:,:,KE,1)
            Optical%empl2f(:,:,KE+1,m) = Optical%empl2f (:,:,KE,1)
          end do
        endif
      endif
!---------------------------------------------------------------------



end subroutine optical_h2o



!####################################################################
! <SUBROUTINE NAME="cfc_optical_depth">
!  <OVERVIEW>
!   Subroutine to compute CFC optical depths
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute CFC optical depths. The code assumes
!   a constant mixing ratio throughout the atmosphere.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cfc_optical_depth (density, Rad_gases, Optical)
!  </TEMPLATE>
!  <IN NAME="density" TYPE="real">
!   density profile of CFC in the atmosphere
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases optical properties input data
!  </IN>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   CFC Optical depth output
!  </INOUT>
! </SUBROUTINE>
!
subroutine cfc_optical_depth (density, Rad_gases, Optical)

!------------------------------------------------------------------
!    cfc_optical_depth computes optical paths for cfc. The code assumes
!    a constant mixing ratio throughout the atmosphere.
!------------------------------------------------------------------

real, dimension (:,:,:),    intent(in)     :: density 
type(radiative_gases_type), intent(in)     :: Rad_gases
type(optical_path_type),    intent(inout)  :: Optical 

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      density
!      Rad_gases
!
!   intent(inout) variable:
!
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real          ::  rrf11, rrf12, rrf113, rrf22
      real          ::  rf11air, rf12air, rf113air, rf22air
      integer       ::  k
      integer       ::  kx

!--------------------------------------------------------------------
!  local variables:
!
!      rrf11
!      rrf12
!      rrf113
!      rrf22
!      rf11air
!      rf12air
!      rf113air
!      rf22air
!      k
!      kx
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      allocate ( Optical%totf11 (size(density,1), size(density,2),    &
                                 size(density,3) ) )
      allocate ( Optical%totf12 (size(density,1), size(density,2),    &
                                 size(density,3) ) )
      allocate ( Optical%totf113(size(density,1), size(density,2),    &
                                 size(density,3) ) )
      allocate ( Optical%totf22 (size(density,1), size(density,2),    &
                                 size(density,3) ) )
      Optical%totf11  = 0.
      Optical%totf12  = 0.
      Optical%totf113 = 0.
      Optical%totf22 = 0.
 
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      kx = size (density,3)

!--------------------------------------------------------------------
!    define cfc mixing ratio conversion factors.
!--------------------------------------------------------------------
      rf11air  = wtmf11/wtmair
      rf12air  = wtmf12/wtmair
      rf113air = wtmf113/wtmair
      rf22air  = wtmf22/wtmair

      rrf11 = Rad_gases%rrvf11*rf11air
      rrf12 = Rad_gases%rrvf12*rf12air
      rrf113 = Rad_gases%rrvf113*rf113air
      rrf22 = Rad_gases%rrvf22*rf22air

!----------------------------------------------------------------------
!    compute summed optical paths for f11,f12, f113 and f22  with the 
!    diffusivity factor of 2 (appropriate for weak-line absorption 
!    limit).
!----------------------------------------------------------------------
      Optical%totf11(:,:,1) = 0.0E+00
      Optical%totf12(:,:,1) = 0.0E+00
      Optical%totf113(:,:,1) = 0.0E+00
      Optical%totf22 (:,:,1) = 0.0E+00
      do k=2,kx           
        Optical%totf11(:,:,k) = Optical%totf11(:,:,k-1) +    &
                                density(:,:,k-1)*rrf11*2.0E+00
        Optical%totf12(:,:,k) = Optical%totf12(:,:,k-1) +    &
                                density(:,:,k-1)*rrf12*2.0E+00
        Optical%totf113(:,:,k) = Optical%totf113(:,:,k-1) +  &
                                 density(:,:,k-1)*rrf113*2.0E+00
        Optical%totf22(:,:,k) = Optical%totf22(:,:,k-1) +    &
                                density(:,:,k-1)*rrf22*2.0E+00
      end do
       
!--------------------------------------------------------------------


end subroutine cfc_optical_depth



!#####################################################################
! <SUBROUTINE NAME="optical_depth_aerosol">
!  <OVERVIEW>
!   Subroutine to compute aerosol optical depths
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to compute aerosol optical depths. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call optical_depth_aerosol (Atmos_input, n, Aerosol,    &
!                                  Aerosol_props, Optical)
!  </TEMPLATE>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to model grid point for radiative 
!   properties calculation
!  </IN>
!  <IN NAME="n" TYPE="integer">
!   aerosol optical index
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties
!  </INOUT>
!  <INOUT NAME="Optical" TYPE="optical_path_type">
!   Aerosol Optical depth output
!  </INOUT>
! </SUBROUTINE>
!
subroutine optical_depth_aerosol ( js, Atmos_input, n, Aerosol,    &
                                  Aerosol_props, Aerosol_diags, &
                                  Optical)

!------------------------------------------------------------------
!
!------------------------------------------------------------------

integer,                       intent(in)    :: js
type(atmos_input_type),        intent(in)    :: Atmos_input
integer,                       intent(in)    :: n
type(aerosol_type),            intent(in)    :: Aerosol
type(aerosol_properties_type), intent(inout) :: Aerosol_props
type(aerosol_diagnostics_type),intent(inout) :: Aerosol_diags
type(optical_path_type),       intent(inout) :: Optical

!-----------------------------------------------------------------
!   intent(in) variables:
!
!      Atmos_input
!      n
!      Aerosol
!
!   intent(inout) variable:
!
!      Aerosol_props
!      Optical
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(Aerosol%aerosol,1),  &
                       size(Aerosol%aerosol,2),  &
                       size(Aerosol%aerosol,3), &
                       size(Aerosol%aerosol,4))   :: aerooptdepspec, &
                                                     aerooptdepspec_cn

      real, dimension (size(Aerosol%aerosol,1),  &
                       size(Aerosol%aerosol,2),  &
                       size(Aerosol%aerosol,3))  :: aerooptdep
!yim
      integer, dimension (size(Aerosol%aerosol,1),  &
                          size(Aerosol%aerosol,2),  &
                          size(Aerosol%aerosol,3))  :: opt_index_v1, &
                          opt_index_v2, opt_index_v3, opt_index_v4, &
                          opt_index_v5, opt_index_v6, opt_index_v7,opt_index_v8

      real, dimension (size(Aerosol%aerosol,3)+1) :: bsum

      real      :: asum
      integer   :: nfields, irh
      integer   ::  N_AEROSOL_BANDS 
      integer   :: i,j,k
      integer   :: ix, jx, kx
      integer   :: nsc, opt_index
!--------------------------------------------------------------------
!  local variables:
!
!      aerooptdepspec
!      aerooptdep
!      irh
!      opt_index_v
!      bsum
!      asum
!      nfields
!      n_aerosol_bands
!      i,j,k
!      ix,jx,kx
!      na, nw, ni  
!      nsc 
!      opt_index
!      
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ix = size (Aerosol%aerosol,1)
      jx = size (Aerosol%aerosol,2)
      kx = size (Aerosol%aerosol,3)
      nfields = size (Aerosol%aerosol,4)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      aerooptdep(:,:,:) = 0.0
      Optical%totaerooptdep(:,:,:,n) = 0.0
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            irh = MIN(100, MAX(0,     &
                      NINT(100.*Atmos_input%aerosolrelhum(i,j,k))))
            opt_index_v1(i,j,k) =     &
                        Aerosol_props%sulfate_index (irh, &
                                             Aerosol_props%ivol(i,j,k) )
            opt_index_v2(i,j,k) =     &
                               Aerosol_props%omphilic_index( irh )
            opt_index_v3(i,j,k) =     &
                               Aerosol_props%bcphilic_index( irh )
            opt_index_v4(i,j,k) =     &
                               Aerosol_props%seasalt1_index( irh )
            opt_index_v5(i,j,k) =     &
                               Aerosol_props%seasalt2_index( irh )
            opt_index_v6(i,j,k) =     &
                               Aerosol_props%seasalt3_index( irh )
            opt_index_v7(i,j,k) =     &
                               Aerosol_props%seasalt4_index( irh )
            opt_index_v8(i,j,k) =     &
                               Aerosol_props%seasalt5_index( irh )
          end do
        end do
      end do

!---------------------------------------------------------------------
!    using relative humidity criterion (where necessary) determine the
!    aerosol category (as an index) appropriate for the aerosol species
!---------------------------------------------------------------------
  do nsc=1,nfields  ! loop on aerosol species
      if (Aerosol_props%optical_index(nsc) > 0 ) then   

      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
                opt_index = Aerosol_props%optical_index(nsc)
                if (opt_index == 0 ) then
                   call error_mesg ('optical_path_init', &
                  'Cannot find aerosol optical properties for species = ' // &
                   TRIM( Aerosol%aerosol_names(nsc) ),  FATAL )
                endif
                aerooptdepspec(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                     (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))*&
                            Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                   aerooptdepspec_cn(i,j,k,nsc) =    &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                   (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                   Aerosol_props%aerextbandlw_cn(n,opt_index)
                 end if
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) == &   
                          Aerosol_props%sulfate_flag  ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v1(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
!yim
     else if (Aerosol_props%optical_index(nsc) == &   
                          Aerosol_props%bc_flag  ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v1(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%omphilic_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v2(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
!yim
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%bcphilic_flag ) then
      if (Rad_control%using_im_bcsul) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v1(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
      else ! (using_im_bcsul)
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v3(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
      endif  ! (using_im_bcsul)
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%seasalt1_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v4(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%seasalt2_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v5(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%seasalt3_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v6(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%seasalt4_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v7(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
     else if (Aerosol_props%optical_index(nsc) ==  &
                        Aerosol_props%seasalt5_flag ) then
      do k = 1,kx         
        do j = 1,jx         
          do i = 1,ix           
            opt_index = opt_index_v8(i,j,k)
                aerooptdepspec(i,j,k,nsc) =     &
                   diffac*Aerosol%aerosol(i,j,k,nsc)*&
                   (1.0 - Aerosol_props%aerssalbbandlw(n,opt_index))* &
                          Aerosol_props%aerextbandlw(n,opt_index)
                if (n == 1) then
                  aerooptdepspec_cn(i,j,k,nsc) =    &
                     diffac*Aerosol%aerosol(i,j,k,nsc)*   &
                 (1.0 - Aerosol_props%aerssalbbandlw_cn(n,opt_index))*&
                        Aerosol_props%aerextbandlw_cn(n,opt_index)
              endif
            end do
          end do
        end do
      endif
   end do

!---------------------------------------------------------------------
!    save optical path contributions from each layer for band4 and the
!    continuum band. note that if the lw scheme is changed to allow
!    longwave scattering then the %absopdep must be defined approp-
!    riately.
!---------------------------------------------------------------------
      if (n == 1) then
        Aerosol_diags%extopdep(:,:,:,:,3) = aerooptdepspec_cn(:,:,:,:)
        Aerosol_diags%absopdep(:,:,:,:,3) = aerooptdepspec_cn(:,:,:,:)
      endif

!---------------------------------------------------------------------
!    sum optical depths over all species and obtain column optical depth
!---------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            asum = 0.0
            do nsc=1,nfields
              asum = asum + aerooptdepspec(i,j,k,nsc)
            end do
            aerooptdep(i,j,k) = asum                         
          end do
        end do
      end do

      do j=1,jx
        do i=1,ix
          bsum(1) = 0.0
          do k=2,kx+1         
            bsum(k) = bsum(k-1) + aerooptdep(i,j,k-1)
          end do
          do k=2,kx+1         
            Optical%totaerooptdep(i,j,k,n) = bsum(k)
          end do
        end do
      end do

!---------------------------------------------------------------------
!    continuum band is the last indx:
!---------------------------------------------------------------------
      n_aerosol_bands = Lw_parameters%n_lwaerosol_bands
      if ( n == n_aerosol_bands) then
        Optical%aerooptdep_KE_15(:,:) = aerooptdep(:,:,kx)
      endif
    
!---------------------------------------------------------------------

end subroutine optical_depth_aerosol


!#####################################################################

                   end module optical_path_mod


