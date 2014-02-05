                  module longwave_tables_mod

! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!  Fei Liu
! </CONTACT>
! <REVIEWER EMAIL="ds@gfdl.noaa.gv">
!  Dan Schwarzkopf
! </REVIEWER>
! <OVERVIEW>
!  This code defines longwave radiation tables, it also
!  allocate, compute related parameters based on prescribed
!  tables.
! </OVERVIEW>
! <DESCRIPTION>
! </DESCRIPTION>
!

!    shared modules:

use mpp_mod,               only: input_nml_file
use fms_mod,               only: open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, close_file

!  shared radiation package modules:

use rad_utilities_mod,     only: rad_utilities_init,       &  
                                 longwave_tables1_type,  &
                                 longwave_tables2_type,  &
                                 longwave_tables3_type,  &
                                 lw_table_type, Lw_parameters,&
                                 table_alloc, mass_1, temp_1, &
                                 Lw_control
use longwave_params_mod,   only: longwave_params_init, NBLW, NBLX, &
                                 NBLY_RSB, NBLY_CKD

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    longwave_tables_mod constructs various tables used in the longwave
!    radiation parameterization.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: longwave_tables.F90,v 19.0 2012/01/06 20:18:37 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'


!---------------------------------------------------------------------
!------    interfaces   ------

public      &
          longwave_tables_init, &
          longwave_tables_end

private      &

!  called from longwave_tables_init:
          idrbtsh2o, id2h2o, table


!---------------------------------------------------------------------
!------  namelist  -----

real       :: dummy = 1.0

namelist / longwave_tables_nml /  &
                                   dummy

!---------------------------------------------------------------------
!---- public data -------


!---------------------------------------------------------------------
!---- private data -------

!--------------------------------------------------------------------
!    define continuum coefficients over special bands, the choices 
!    depend on model architecture. the program gasbnd is used.
!--------------------------------------------------------------------
real, dimension(:), allocatable :: afach4, afan2o
              
real, dimension(:), allocatable :: fbdlo_12001400, fbdhi_12001400
real, dimension(:), allocatable :: dummy_ch4n2o

real, dimension(:), allocatable :: bdlahcn, bdhahcn

real, dimension(:), allocatable :: bfach4, bfan2o             

real, dimension(:), allocatable :: dch4, dn2o, ech4, en2o
real                            :: d171n2o, e171n2o

real, dimension(:), allocatable :: acomb, bcomb, apcm, bpcm, atpcm,  &
                                   btpcm, bdlocm, bdhicm

integer, parameter              :: NTTABH2O   = 28
integer, parameter              :: NUTABH2O   = 181

real, dimension (NBLW)          :: bandlo, bandhi, arndm, brndm, betad
integer, dimension(40)          :: iband
real, dimension(3)              :: ao3cm, bo3cm
real, dimension(2)              :: ab15cm

integer                         :: NBTRG, NBTRGE, NBLY
real                            :: apwd, bpwd, atpwd, btpwd, bdlowd, &
                                   bdhiwd 

logical :: module_is_initialized = .false.   !  module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
!#####################################################################

! <SUBROUTINE NAME="longwave_tables_init">
!  <OVERVIEW>
!   Constructor of longwave_tables module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Defines continuum coefficients and random band parameters for longwave
!   gas species.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_tables_init (Lw_tables, tabsr,   &
!                        tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)
!  </TEMPLATE>
!  <IN NAME="Lw_tables" TYPE="lw_table_type">
!   Contains the tables used in longwave radiation
!  </IN>
!  <OUT NAME="tabsr" TYPE="longwave_tables3_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1w" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
! </SUBROUTINE>
!
subroutine longwave_tables_init (Lw_tables, tabsr, tab1, tab2, tab3, &
                                 tab1w, tab1a, tab2a, tab3a)

!--------------------------------------------------------------------
!    longwave_tables_init is the constructor for longwave_tables_mod.
!--------------------------------------------------------------------

type(lw_table_type),          intent(inout) :: Lw_tables
type(longwave_tables3_type),  intent(inout) :: tabsr
type (longwave_tables1_type), intent(inout) :: tab1, tab2, tab3, tab1w
type (longwave_tables2_type), intent(inout) :: tab1a, tab2a, tab3a

!---------------------------------------------------------------------
!   intent(inout) variables:
!
!    Lw_tables
!    tabsr
!    tab1
!    tab2
!    tab3
!    tab1w
!    tab1a
!    tab2a
!    tab3a
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    define continuum coefficients over special bands, the choices 
!    depend on model architecture. the program gasbnd is used.
!---------------------------------------------------------------------
      real                          :: apwd_c, bpwd_c, atpwd_c,    &
                                       btpwd_c, bdlowd_c, bdhiwd_c
      real, dimension (NBLY_CKD)    :: acomb_c, bcomb_c, apcm_c,  &
                                       bpcm_c, atpcm_c,   &
                                       btpcm_c, bdlocm_c,  bdhicm_c
      integer                       :: inrad, k
      integer                       :: subb
      integer, dimension(5)         :: no_h2o12001400bands = &
                                        (/ 1, 2, 4, 10, 20 /)
 
!---------------------------------------------------------------------
!    2) 160-560 (as 40 bands). program gasbnd is used with 10 cm-1
!    bandwidth. iband is straightforward mapping.
!---------------------------------------------------------------------
      integer, dimension(40)        :: iband_c
      data iband_c /    &
          1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   &
         11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  &
         21,  22,  23,  24,  25,  26,  27,  28,  29,  30,   &
         31,  32,  33,  34,  35,  36,  37,  38,  39,  40/ 

!----------------------------------------------------------------------
!    define random band parameters for special bands. the choices 
!    depend on model architecture. the program gasbnd is used.
!    2) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb
!    on. combined bands defined according to index iband.
!----------------------------------------------------------------------
      real                        ::   apwd_n, bpwd_n, atpwd_n,   &
                                       btpwd_n, bdlowd_n, bdhiwd_n
      real, dimension (NBLY_RSB)  ::   acomb_n, bcomb_n, apcm_n,  &
                                       bpcm_n, atpcm_n, btpcm_n,  &
                                       bdlocm_n,  bdhicm_n
      real, dimension(NBLY_RSB)   ::   dummy_n
      real                        ::   dum
 
      integer, dimension(40)      ::   iband_n
      data iband_n /   &
          2,   1,   2,   2,   1,   2,   1,   3,   2,   2,   &
          3,   2,   2,   4,   2,   4,   2,   3,   3,   2,  &
          4,   3,   4,   3,   7,   5,   6,   7,   6,   5,  &
          7,   6,   7,   8,   6,   6,   8,   8,   8,   8/
 
!---------------------------------------------------------------------
!    miscellaneous variables:
!    unit            io unit number used for namelist file
!    ierr            error code
!    io              error status returned from io operation
!    k4
!    n4
!---------------------------------------------------------------------
      integer    :: unit, ierr, io, logunit


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
      call longwave_params_init

!-----------------------------------------------------------------------
!    read namelist.
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=longwave_tables_nml, iostat=io)
      ierr = check_nml_error(io,"longwave_tables_nml")
#else
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=longwave_tables_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'longwave_tables_nml')
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
                          write (logunit, nml=longwave_tables_nml)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_1992' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
          inrad = open_namelist_file ('INPUT/h2ocoeff_ckd_speccombwidebds_hi92')
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) apwd_c   ! ckd capphi coeff for 560-800 band
          read (inrad,9000) bpwd_c   ! ckd cappsi coeff for 560-800 band
          read (inrad,9000) atpwd_c  ! ckd capphi coeff for 560-800 band
          read (inrad,9000) btpwd_c  ! ckd cappsi coeff for 560-800 band
          read (inrad,9000) bdlowd_c  ! lo freq of 560-800 band
          read (inrad,9000) bdhiwd_c  ! hi freq of 560-800 band
!  ckd rndm coeff for 40 bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (acomb_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bcomb_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (apcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bpcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (atpcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (btpcm_c(k),k=1,NBLY_CKD)
!  ckd lo/hi freq for 40 bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (bdlocm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bdhicm_c(k),k=1,NBLY_CKD)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          inrad = open_namelist_file ('INPUT/h2ocoeff_rsb_speccombwidebds_hi92')
          read (inrad,9000) dum     
          read (inrad,9000) dum    
          read (inrad,9000) apwd_n   ! rsb capphi coeff for 560-800 band
          read (inrad,9000) bpwd_n   ! rsb cappsi coeff for 560-800 band
          read (inrad,9000) atpwd_n  ! rsb capphi coeff for 560-800 band
          read (inrad,9000) btpwd_n  ! rsb cappsi coeff for 560-800 band
          read (inrad,9000) bdlowd_n ! lo freq of 560-800 band
          read (inrad,9000) bdhiwd_n ! hi freq of 560-800 band
          read (inrad,9000) dum   
!  rsb rndm coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (acomb_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bcomb_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (apcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bpcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (atpcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (btpcm_n(k),k=1,NBLY_RSB)
!  rsb lo/hi freq for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (bdlocm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bdhicm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
        endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          inrad = open_namelist_file ('INPUT/h2ocoeff_ckd_speccombwidebds_hi00')
          read (inrad,9000) dum
          read (inrad,9000) dum
          read (inrad,9000) apwd_c   ! ckd capphi coeff for 560-800 band
          read (inrad,9000) bpwd_c   ! ckd cappsi coeff for 560-800 band
          read (inrad,9000) atpwd_c  ! ckd capphi coeff for 560-800 band
          read (inrad,9000) btpwd_c  ! ckd cappsi coeff for 560-800 band
          read (inrad,9000) bdlowd_c ! lo freq of 560-800 band
          read (inrad,9000) bdhiwd_c ! hi freq of 560-800 band
!  ckd rndm coeff for 40 bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (acomb_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bcomb_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (apcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bpcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (atpcm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (btpcm_c(k),k=1,NBLY_CKD)
!  ckd lo/hi freq for 40 bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (bdlocm_c(k),k=1,NBLY_CKD)
          read (inrad,9000) (bdhicm_c(k),k=1,NBLY_CKD)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          inrad = open_namelist_file ('INPUT/h2ocoeff_rsb_speccombwidebds_hi00')
          read (inrad,9000) dum     
          read (inrad,9000) dum    
          read (inrad,9000) apwd_n   ! rsb capphi coeff for 560-800 band
          read (inrad,9000) bpwd_n   ! rsb cappsi coeff for 560-800 band
          read (inrad,9000) atpwd_n  ! rsb capphi coeff for 560-800 band
          read (inrad,9000) btpwd_n  ! rsb cappsi coeff for 560-800 band
          read (inrad,9000) bdlowd_n ! lo freq of 560-800 band
          read (inrad,9000) bdhiwd_n ! hi freq of 560-800 band
          read (inrad,9000) dum   
!  rsb rndm coeff for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (acomb_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bcomb_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (apcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bpcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (atpcm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (btpcm_n(k),k=1,NBLY_RSB)
!  rsb lo/hi freq for 8 comb bands (160-560) and 8 wide bands (560-1400)
          read (inrad,9000) (bdlocm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (bdhicm_n(k),k=1,NBLY_RSB)
          read (inrad,9000) (dummy_n(k),k=1,NBLY_RSB)
        endif
      endif
      call close_file (inrad)

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (Lw_parameters%NBTRG_iz) then
        NBTRG  = Lw_parameters%NBTRG
      else
        call error_mesg ('longwave_tables_mod', &
                       ' Lw_parameters%NBTRG not yet defined', FATAL) 
      endif
      if (Lw_parameters%NBTRGE_iz) then
        NBTRGE = Lw_parameters%NBTRGE
      else
        call error_mesg ('longwave_tables_mod', &
                       ' Lw_parameters%NBTRGE not yet defined', FATAL)
      endif

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( fbdlo_12001400 (NBTRGE) )
        allocate ( fbdhi_12001400 (NBTRGE) )
        allocate ( dummy_ch4n2o (NBTRGE) )
      endif
      if (NBTRG  > 0) then
        allocate ( afach4 (NBTRG ) )
        allocate ( afan2o (NBTRG ) )
        allocate ( bdlahcn(NBTRG ) )
        allocate ( bdhahcn(NBTRG ) )
        allocate ( bfach4 (NBTRG ) )
        allocate ( bfan2o (NBTRG ) )
        allocate ( dch4   (NBTRG ) )
        allocate ( dn2o   (NBTRG ) )
        allocate ( ech4   (NBTRG ) )
        allocate ( en2o   (NBTRG ) )
      endif

!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (NBTRGE > 0) then
        if (trim(Lw_control%linecatalog_form) == 'hitran_1992') then
          inrad = open_namelist_file ('INPUT/h2o12001400_hi92_data')
        else if(trim(Lw_control%linecatalog_form) == 'hitran_2000') then
          inrad = open_namelist_file ('INPUT/h2o12001400_hi00_data')
        endif

!----------------------------------------------------------------------
!    read in random coefficients for 1200-1400 freq region, spacing
!    through the data until  those appropriate for NBTRGE h2o bands
!    are reached. note: unless a continuum is inserted beyond 1200
!    cm-1, the band coefficients are independent of continuum type.
!---------------------------------------------------------------------
        do subb = 1,5    ! 5 = no. band divisions in h2o 1200-1400 data
          if (NBTRGE == no_h2o12001400bands(subb)) then

!---------------------------------------------------------------------
!    read and process data for sub-band number from data matching NBTRGE
!    then exit subb loop
!---------------------------------------------------------------------
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (dummy_ch4n2o(k),k=1,NBTRGE)
            read (inrad,2001) (fbdlo_12001400(k),k=1,NBTRGE)
            read (inrad,2001) (fbdhi_12001400(k),k=1,NBTRGE)
            exit
          else if (subb < 5) then 

!---------------------------------------------------------------------
!    read data for sub-band number from  data not matching NBTRGE
!---------------------------------------------------------------------
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
            read (inrad,2001) (dummy_ch4n2o(k),k=1,  &
                                          no_h2o12001400bands(subb))
          else

!---------------------------------------------------------------------
!    failure of any sub-band number to match NBTRGE
!---------------------------------------------------------------------
            call error_mesg ('longwave_tables_mod', &
                 'NBTRGE is inconsistent with available data', FATAL)
          endif
        end do
        call close_file(inrad)
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (trim(Lw_control%linecatalog_form) == 'hitran_1992') then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          NBLY = NBLY_CKD
          call id2h2o ('INPUT/id2h2obdckd2p1')
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          NBLY = NBLY_RSB
          call id2h2o ('INPUT/id2h2obdfull')

!----------------------------------------------------------------------
!  read roberts continuum data for self-broadened h2o continuum
!----------------------------------------------------------------------
          call idrbtsh2o
        endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      else if (trim(Lw_control%linecatalog_form) == 'hitran_2000' ) then
        if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
            trim(Lw_control%continuum_form) == 'ckd2.4' ) then
          NBLY = NBLY_CKD
          call id2h2o ('INPUT/h2ocoeff_ckd_0_3000_10cm_hi00')
        else if (trim(Lw_control%continuum_form) == 'rsb' ) then
          NBLY = NBLY_RSB
          call id2h2o ('INPUT/h2ocoeff_rsb_0_3000_10cm_hi00')

!----------------------------------------------------------------------
!  read roberts continuum data for self-broadened h2o continuum
!----------------------------------------------------------------------
          call idrbtsh2o
        endif
      endif

      allocate  ( acomb(NBLY))
      allocate  ( bcomb(NBLY))
      allocate  ( apcm (NBLY))
      allocate  ( bpcm (NBLY))
      allocate  ( atpcm(NBLY))
      allocate  ( btpcm(NBLY))
      allocate  (bdlocm(NBLY))
      allocate  (bdhicm(NBLY))

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (trim(Lw_control%continuum_form) == 'ckd2.1' .or.     &
          trim(Lw_control%continuum_form) == 'ckd2.4' ) then
        apwd = apwd_c
        bpwd = bpwd_c
        atpwd = atpwd_c
        btpwd = btpwd_c
        bdlowd = bdlowd_c
        bdhiwd = bdhiwd_c
        iband = iband_c
        acomb = acomb_c
        bcomb = bcomb_c
        apcm = apcm_c
        bpcm = bpcm_c
        atpcm = atpcm_c
        btpcm = btpcm_c
        bdlocm = bdlocm_c
        bdhicm = bdhicm_c
      else if (trim(Lw_control%continuum_form) == 'rsb' ) then
        apwd = apwd_n
        bpwd = bpwd_n
        atpwd = atpwd_n
        btpwd = btpwd_n
        bdlowd = bdlowd_n
        bdhiwd = bdhiwd_n
        iband = iband_n
        acomb = acomb_n
        bcomb = bcomb_n
        apcm = apcm_n
        bpcm = bpcm_n
        atpcm = atpcm_n
        btpcm = btpcm_n
        bdlocm = bdlocm_n
        bdhicm = bdhicm_n
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      call table_alloc (tab1 , NTTABH2O, NUTABH2O)
      call table_alloc (tab2 , NTTABH2O, NUTABH2O)
      call table_alloc (tab3 , NTTABH2O, NUTABH2O)
      call table_alloc (tab1w, NTTABH2O, NUTABH2O)
      if (NBTRGE > 0) then
        call table_alloc (tab1a, NTTABH2O, NUTABH2O, NBTRGE)
        call table_alloc (tab2a, NTTABH2O, NUTABH2O, NBTRGE)
        call table_alloc (tab3a, NTTABH2O, NUTABH2O, NBTRGE)
      endif
      call table_alloc (tabsr, NTTABH2O, NBLY    )

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      call table (tabsr, tab1, tab2, tab3, tab1w, &
                  tab1a, tab2a, tab3a )

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      allocate (Lw_tables%bdlocm(NBLY))
      allocate (Lw_tables%bdhicm(NBLY))
      allocate (Lw_tables%iband (40))
      allocate (Lw_tables%bandlo (NBLW))
      allocate (Lw_tables%bandhi (NBLW))
      Lw_tables%bdlocm = bdlocm
      Lw_tables%bdhicm = bdhicm
      Lw_tables%iband  = iband 
      Lw_tables%bandlo = bandlo
      Lw_tables%bandhi = bandhi

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------
2001  format (5e14.6)
9000  format (5e14.6)

!----------------------------------------------------------------------



end subroutine longwave_tables_init



!#####################################################################


! <SUBROUTINE NAME="longwave_tables_end">
!  <OVERVIEW>
!   Destructor of longwave_tables module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Closes out longwave tables module.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call longwave_tables_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine longwave_tables_end

!--------------------------------------------------------------------
!    longwave_tables_end is the destructor for longwave_tables_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_tables_mod', &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

 

end subroutine longwave_tables_end 




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="idrbtsh2o">
!  <OVERVIEW>
!   Subroutine to read h2o roberts continuum quantities used in longwave
!   radiation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine reads h2o roberts continuum quantities used in
!   longwave radiation from an INPUT file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call idrbtsh2o
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine idrbtsh2o

!----------------------------------------------------------------------
!    idrbtsh2o reads h2o roberts continuum quantities used in
!    longwave radiation.
!    author: m. d. schwarzkopf
!    revised: 1/1/96
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      integer   :: inrad  !  unit number for i/o
      integer   :: k      !  do-loop index

!-----------------------------------------------------------------------
!    the following roberts continuum coefficients are computed using the
!    program (gasbnd) over the 0-3000 cm-1 range with 10 cm-1 bandwidth.
!-----------------------------------------------------------------------
      inrad = open_namelist_file ('INPUT/id2h2orbts')
      read (inrad, FMT = '(5e14.6)') (betad(k),k=1,NBLW)
      call close_file (inrad)

!---------------------------------------------------------------------
 
end subroutine idrbtsh2o


!#####################################################################

subroutine id2h2o (filename)

!---------------------------------------------------------------------
!    id2h2o reads h2o random model band parameters used for 
!    longwave radiation
!    references:
!     (1) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!    author: m. d. schwarzkopf
!    revised: 1/1/96
!    certified:  radiation version 1.0
!---------------------------------------------------------------------

character(len=*), intent(in)   :: filename

!---------------------------------------------------------------------
!  intent(in) variable:
!
!    filename
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension (NBLW) :: dummy   ! dummy array

      integer   :: inrad  !  unit number for i/o
      integer   :: k      !  do-loop index

!-----------------------------------------------------------------------
!    the following h2o random band parameters are obtained from the
!    afgl 1992 HITRAN tape. parameters are obtained using an auxi-
!    liary program (gasbnd). values depend on assumptions as to
!    line shape, line strength and width. The inputted values span
!    the 0-3000 cm-1 range, with 10 cm-1 bandwidth. other parameter
!    values used in the program are obtained separately.
!-----------------------------------------------------------------------
      inrad = open_namelist_file (filename)
      read (inrad,9000) (arndm(k),k=1,NBLW)
      read (inrad,9000) (brndm(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (bandlo(k),k=1,NBLW)
      read (inrad,9000) (bandhi(k),k=1,NBLW)
      call close_file (inrad)

!--------------------------------------------------------------------
9000  format(5e14.6)

!--------------------------------------------------------------------


end subroutine id2h2o




!#####################################################################
! <SUBROUTINE NAME="table">
!  <OVERVIEW>
!   Subroutine to compute table entries used in longwave radiation
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine computes the table entries used in longwave radiation
!  </DESCRIPTION>
!  <TEMPLATE>
!   call table  (tabsr, tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)
!  </TEMPLATE>
!  <OUT NAME="tabsr" TYPE="longwave_tables3_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1w" TYPE="longwave_tables1_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab1a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tabs2a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
!  <OUT NAME="tab3a" TYPE="longwave_tables2_type">
!   Contains the tables used in longwave radiation
!  </OUT>
! </SUBROUTINE>
!
subroutine table  (tabsr, tab1, tab2, tab3, tab1w, tab1a, tab2a, tab3a)

!---------------------------------------------------------------------
!    table computes table entries used in longwave radiation.  
!    author: m. d. schwarzkopf
!    revised: 1/1/93
!    certified:  radiation version 1.0
!---------------------------------------------------------------------
 
type(longwave_tables3_type), intent(inout)   :: tabsr
type(longwave_tables1_type), intent(inout)   :: tab1, tab2, tab3, tab1w
type(longwave_tables2_type), intent(inout)   :: tab1a, tab2a, tab3a

!----------------------------------------------------------------------
!  intent(inout) variables:
!
!     tabsr
!     tab1
!     tab2
!     tab3
!     tab1w
!     tab1a
!     tab2a
!     tab3a
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (:,:), allocatable   :: r1a, r2a, s2a, t3a,   &
                                              sum4a, sum6a, sum7a, sum8a
      real, dimension(:,:,:),allocatable   :: suma, sumdbea, sum3a 
      real, dimension (NBLW)               :: alfanb, anb, arotnb,   &
                                              betanb, bnb, centnb, delnb
      real, dimension (30)                 :: cnusb, dnusb
      real, dimension (NTTABH2O,NBLW)      :: dbdtnb, src1nb
      real, dimension (NTTABH2O,NBLX)      :: srcwd        
      real, dimension (NTTABH2O, NUTABH2O) :: sumdbe, sum, sum3, sumwde
      real, dimension (NTTABH2O)           :: ddsc, fortcu, r1, r1wd, &
                                              r2, s2, sc, srcs, sum4, &
                                              sum4wd, sum6, sum7, sum8,&
                                              t3, tfour, x, x1, xtemv
      real, dimension (NUTABH2O)           :: expo, fac, x2, zmass, &
                                              zroot
      integer                              :: n, m, ioffset, itab,   &
                                              jtab, nsubds, nsb,  iter
      real                                 :: zmassincr, cent, del,&
                                              bdlo, bdhi, anu, c1,   &
                                              freq_cutoff

!---------------------------------------------------------------------
!  local variables:
!
!    r1a
!    r2a
!    s2a
!    t3a
!    sum4a
!    sum6a
!    sum7a
!    sum8a
!    suma
!    suma
!    sumdbea
!    sum3a
!    ETC. 

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (NBTRGE > 0) then
        allocate ( r1a     (NTTABH2O,NBTRGE) )
        allocate (r2a      (NTTABH2O,NBTRGE) )
        allocate (s2a      (NTTABH2O,NBTRGE) )
        allocate ( t3a     (NTTABH2O,NBTRGE) )
        allocate ( suma    (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sumdbea (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sum3a   (NTTABH2O,NUTABH2O,NBTRGE) )
        allocate ( sum4a   (NTTABH2O,NBTRGE) )
        allocate ( sum6a   (NTTABH2O,NBTRGE) )
        allocate ( sum7a   (NTTABH2O,NBTRGE) )
        allocate (sum8a    (NTTABH2O,NBTRGE) )
      endif

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      if (Lw_parameters%offset_iz) then
        ioffset = Lw_parameters%offset
      else
        call error_mesg ('longwave_tables_mod', &
                 ' Lw_parameters%offset not yet defined', FATAL)
      endif

!--------------------------------------------------------------------- 
!     compute local quantities and ao3, bo3, ab15 for narrow bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        anb   (n) = arndm(n) 
        bnb   (n) = brndm(n) 
        centnb(n) = 0.5E+00*(bandlo(n) + bandhi(n)) 
        delnb (n) = bandhi(n) - bandlo(n)
        betanb(n) = betad(n)
      enddo

!---------------------------------------------------------------------
!    compute a*b and sqrt(a*b) for all 10 cm-1 frequency bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        alfanb(n) = bnb(n)*anb(n) 
        arotnb(n) = SQRT(alfanb(n))
      enddo

!-------------------------------------------------------------------
!   define critical frequency (cutoff for wide band ?? )
!------------------------------------------------------------------
      if (NBTRGE > 0) then
        freq_cutoff = 1400.
      else
        freq_cutoff = 1200.
      endif

!---------------------------------------------------------------------
!    begin table computations here.  compute temperatures and masses
!    for table entries.
!    note: the dimensioning and initialization of xtemv and other
!    arrays with dimension of NTTABH2O=28 imply a restriction of model 
!    temperatures from 100k to 370k.
!    the dimensioning of zmass, zroot and other arrays with 
!    dimension of NUTABH2O=181 imply a restriction of model h2o amounts
!    such that optical paths are between 10**-16 and 10**2, in cgs 
!    units (index 2-181), plus zero (index 1).
!---------------------------------------------------------------------
      zmass(1) = 0.0
      zmass(2) = 10.0E+00**mass_1%min_val
      zmassincr = 10.0E+00**mass_1%tab_inc
 
!---------------------------------------------------------------------
!    the definition of zmassincr as 10**0.1 is slightly different from
!    all previous versions, in which it is 1.258925411E+00. This
!    produces slightly different answers (fluxes differ by 1.0e-6 W/m2).
!---------------------------------------------------------------------
      do jtab=3,NUTABH2O
        zmass(jtab) = zmass(jtab-1)*zmassincr
      enddo
      zroot(1) = 0.0
      do jtab=2,NUTABH2O
        zroot(jtab) = SQRT(zmass(jtab))
      enddo 
      do itab=1,NTTABH2O
        xtemv (itab) = temp_1%min_val + temp_1%tab_inc*(itab-1)
        tfour (itab) = xtemv(itab)**4
        fortcu(itab) = 4.0E+00*xtemv(itab)**3
      enddo
      
!---------------------------------------------------------------------
!    the computation of source, dsrce is needed only for the combined 
!    wide band case.  to obtain them,  the source must be computed 
!    for each of the NBLX wide bands srcwd then combined using iband
!    into source.
!---------------------------------------------------------------------
      do n=1,NBLY
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = 0.0E+00
        enddo
      enddo
      do n=1,NBLX
        do itab=1,NTTABH2O
          srcwd(itab,n) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    begin frequency loop.
!---------------------------------------------------------------------
      do n=1,NBLX 
  
!---------------------------------------------------------------------
!     the 160-560 cm-1 region
!---------------------------------------------------------------------
        if (n .LE. 40) then
          cent = centnb(n+16) 
          del  = delnb (n+16) 
          bdlo = bandlo(n+16) 
          bdhi = bandhi(n+16) 
 
!---------------------------------------------------------------------
!      the 560-1200 cm-1 region, and the 2270-2380 cm-1 region
!---------------------------------------------------------------------
        else
          cent = 0.5E+00*(bdlocm(n-32+ioffset) + bdhicm(n-32+ioffset))
          del  = bdhicm(n-32+ioffset) - bdlocm(n-32+ioffset)
          bdlo = bdlocm(n-32+ioffset)
          bdhi = bdhicm(n-32+ioffset)
        endif

!---------------------------------------------------------------------
!    for purposes of accuracy, all evaluations of planck functions
!    are made on 10 cm-1 intervals, then summed into the NBLX wide 
!    bands.  the last subband may be narrower than 10 cm-1.
!---------------------------------------------------------------------
        nsubds = (del - 1.0E-03)/10 + 1
        do nsb=1,nsubds 
          if(nsb .NE. nsubds) then 
            cnusb(nsb) = 10.0E+00*(nsb - 1) + bdlo + 5.0E+00
            dnusb(nsb) = 10.0E+00
          else
            cnusb(nsb) = 0.5E+00*(10.0E+00*(nsb - 1) + bdlo + bdhi)
            dnusb(nsb) = bdhi -  (10.0E+00*(nsb - 1) + bdlo)
          endif 
          c1 = 3.7412E-05*cnusb(nsb)**3

!---------------------------------------------------------------------
!    begin temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x    (itab)   = 1.4387E+00*cnusb(nsb)/xtemv(itab)
            x1   (itab)   = EXP(x(itab)) 
            srcs (itab)   = c1/(x1(itab) - 1.0E+00)
            srcwd(itab,n) = srcwd(itab,n) + srcs(itab)*dnusb(nsb)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!    the following loops create the combined wide band quantities 
!    source and dsrce.  the first 40 bands map to bands 1 to 8 in
!    source and dsrce for the bands in the case of using the rsb
!    continuum . the first 40 bands map to bands 1 to 40 if the
!    band structure for the ckd continuum is used.
!---------------------------------------------------------------------
      do n=1,40
        do itab=1,NTTABH2O 
          tabsr%vae  (itab,iband(n)) = tabsr%vae(itab,iband(n)) +      &
                                  srcwd(itab,n)
        enddo
      enddo
      do n=9+ioffset,NBLY  
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = srcwd(itab,n+32-ioffset)
        enddo
      enddo
      do n=1,NBLY
        do itab=1,NTTABH2O-1 
          tabsr%td(itab,n) = (tabsr%vae(itab+1,n) -   &
                      tabsr%vae(itab,n))*0.1E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    first compute planck functions src1nb and derivatives dbdtnb 
!    for use in table evaluations.  these are different from source,
!    dsrce because different frequency points are used in evaluation,
!    the frequency ranges are different, and the derivative algorithm
!    is different.
!---------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)
        del  = delnb (n)

!---------------------------------------------------------------------
!    note: at present, the iter loop is only used for iter=2.  the 
!    loop structure is kept so that in the future, we may use a
!    quadrature scheme for the planck function evaluation, rather
!    than use the mid-band frequency.
!---------------------------------------------------------------------
        do iter=2,2
          anu = cent + 0.5E+00*(iter - 2)*del 
          c1  = (3.7412E-05)*anu**3
!---------------------------------------------------------------------
!    temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x  (itab)      = 1.4387E+00*anu/xtemv(itab)
            x1 (itab)      = EXP(x(itab))
            sc (itab)      = c1/((x1(itab) - 1.0E+00) + 1.0E-20) 
            sc (itab)      = c1/(x1(itab) - 1.0E+00)
            ddsc(itab)     = sc(itab)/(x1(itab)-1.0E+00)*x1(itab)*  &
                             x(itab)/xtemv(itab)
            src1nb(itab,n) = del*sc (itab)
            dbdtnb(itab,n) = del*ddsc(itab)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!    next compute r1, r2, s2, and t3 coefficients used for e3 
!    function when the optical path is less than 10**-4.  in this 
!    case, we assume a different dependence on zmass.  also obtain 
!    r1wd, which is r1 summed over the 160-560 cm-1 range.
!---------------------------------------------------------------------
      do itab=1,NTTABH2O
        sum4  (itab) = 0.0E+00
        sum6  (itab) = 0.0E+00
        sum7  (itab) = 0.0E+00
        sum8  (itab) = 0.0E+00
        sum4wd(itab) = 0.0E+00
      enddo

      if (NBTRGE > 0) then
        sum4a (:,:)    = 0.0E+00
        sum6a (:,:)    = 0.0E+00
        sum7a (:,:)    = 0.0E+00
        sum8a (:,:)    = 0.0E+00
      endif

      do n=1,NBLW 
        cent = centnb(n)
!---------------------------------------------------------------------
!#ifndef ch4n2o
!    perform summations for frequency ranges of 0-560, 1200-2200 cm-1 
!#else   ch4n2o
!    perform summations for frequency ranges of 0-560, 1400-2200 cm-1 
!#endif ch4n2o
!---------------------------------------------------------------------
        if (cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.   &
            cent .LE. 2.2E+03) then
          do itab=1,NTTABH2O 
            sum4(itab) = sum4(itab) + src1nb(itab,n)
            sum6(itab) = sum6(itab) + dbdtnb(itab,n)
            sum7(itab) = sum7(itab) + dbdtnb(itab,n)*arotnb(n)
            sum8(itab) = sum8(itab) + dbdtnb(itab,n)*alfanb(n)
          enddo
        endif

        if (NBTRGE > 0) then

!---------------------------------------------------------------------
!    perform summations for frequency range of 1200-1400 cm-1
!    for sum4a, sum6a, sum7a, and sum8a. the computation depends
!    on the value of NBTRGE.
!---------------------------------------------------------------------
          if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
            do m=1,NBTRGE
              if (cent .GT. fbdlo_12001400(m) .AND.   &
                  cent .LE. fbdhi_12001400(m)) then
                sum4a(:,m) = sum4a(:,m) + src1nb(:,n)
                sum6a(:,m) = sum6a(:,m) + dbdtnb(:,n)
                sum7a(:,m) = sum7a(:,m) + dbdtnb(:,n)*arotnb(n)
                sum8a(:,m) = sum8a(:,m) + dbdtnb(:,n)*alfanb(n)
              endif
            enddo
          endif
        endif

!---------------------------------------------------------------------
!    perform summations over 160-560 cm-1 frequency range for e1 
!    calculations sum4wd.
!---------------------------------------------------------------------
        if (cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then
          do itab=1,NTTABH2O 
            sum4wd(itab) = sum4wd(itab) + src1nb(itab,n)
          enddo
        endif 
      enddo

!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      do itab=1,NTTABH2O
        r1(itab)   = sum4(itab)/tfour (itab)
        r2(itab)   = sum6(itab)/fortcu(itab) 
        s2(itab)   = sum7(itab)/fortcu(itab) 
        t3(itab)   = sum8(itab)/fortcu(itab) 
        r1wd(itab) = sum4wd(itab)/tfour(itab)
      enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          sum   (itab,jtab) = 0.0E+00 
          sumdbe(itab,jtab) = 0.0E+00
          sum3  (itab,jtab) = 0.0E+00
          sumwde(itab,jtab) = 0.0E+00
        enddo
      enddo
      if (NBTRGE > 0) then
        do m=1,NBTRGE
          r1a(:,m)   = sum4a(:,m)/tfour(:)
          r2a(:,m)   = sum6a(:,m)/fortcu(:)
          s2a(:,m)   = sum7a(:,m)/fortcu(:)
          t3a(:,m)   = sum8a(:,m)/fortcu(:)
        enddo
        suma   (:,:,:) = 0.0E+00 
        sumdbea(:,:,:) = 0.0E+00
        sum3a  (:,:,:) = 0.0E+00
      endif

!---------------------------------------------------------------------
!    frequency loop begins.
!--------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)

!---------------------------------------------------------------------
!    perform calculations for frequency ranges of 0-560, 
!#ifndef ch4n2o
!    1200-2200 cm-1.
!#else   ch4n2o
!    1400-2200 cm-1.
!#endif  ch4n2o
!---------------------------------------------------------------------
        if (cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.    &
            cent .LE. 2.2E+03) then
          do jtab=1,NUTABH2O 
            x2  (jtab) = arotnb(n)*zroot(jtab) 
            expo(jtab) = EXP( - x2(jtab))
          enddo
          do jtab=122,NUTABH2O
            fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*expo(jtab))/ &
                        (alfanb(n)*zmass(jtab))
          enddo
          do jtab=1,NUTABH2O 
            do itab=1,NTTABH2O
              sum   (itab,jtab) = sum   (itab,jtab) +   &
                                  src1nb(itab,n)*expo(jtab)
              sumdbe(itab,jtab) = sumdbe(itab,jtab) +    &
                                  dbdtnb(itab,n)*expo(jtab)
            enddo
          enddo 
          do jtab=122,NUTABH2O
            do itab=1,NTTABH2O 
              sum3(itab,jtab) = sum3(itab,jtab) +    &
                                dbdtnb(itab,n)*fac(jtab)
            enddo 
          enddo 
        endif

!-------------------------------------------------------------------
!    perform calculations over the frequency range 1200-1400 cm-1. 
!    the calculations depend on the value of NBTRGE.
!-------------------------------------------------------------------
        if (NBTRGE > 0) then 
          if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
            do m=1,NBTRGE
              if (cent .GT. fbdlo_12001400(m) .AND.   &
                  cent .LE. fbdhi_12001400(m)) then
                x2  (:) = arotnb(n)*zroot(:) 
                expo(:) = EXP( - x2(:))
                do jtab=122,NUTABH2O
                  fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*  &
                               expo(jtab))/(alfanb(n)*zmass(jtab))
                enddo
                do jtab=1,NUTABH2O 
                  suma(:,jtab,m)    = suma(:,jtab,m) +  &
                                      src1nb(:,n)*expo(jtab)
                  sumdbea(:,jtab,m) = sumdbea(:,jtab,m) +   &
                                      dbdtnb(:,n)*expo(jtab)
                enddo
                do jtab=122,NUTABH2O 
                  sum3a(:,jtab,m)   = sum3a(:,jtab,m) +   &
                                      dbdtnb(:,n)*fac(jtab)
                enddo
              endif
            enddo
          endif
        endif

!---------------------------------------------------------------------
!    compute sum over 160-560 cm-1 range for use in e1 calculations
!    sumwde.
!---------------------------------------------------------------------
        if (cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then 
          do jtab=1,NUTABH2O
            do itab=1,NTTABH2O 
              sumwde(itab,jtab) = sumwde(itab,jtab) +     &
                                  src1nb(itab,n)*expo(jtab)
            enddo
          enddo 
        endif 
      enddo

!--------------------------------------------------------------------
!    frequency loop ends
!--------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = sum(itab,jtab)/tfour(itab)
          tab2%vae(itab,jtab) = sumdbe(itab,jtab)/fortcu(itab)
        enddo 
      enddo
      do jtab=122,NUTABH2O
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = sum3(itab,jtab)/fortcu(itab)
        enddo
      enddo
      do jtab=1,3
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = r1(itab)
        enddo
      enddo
      do jtab=1,121
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = r2(itab)/2.0E+00 -    &
                                s2(itab)*zroot(jtab)/3.0E+00 +   &
                                t3(itab)*zmass(jtab)/8.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    compute e1 tables for 160-560 cm-1 bands.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = sumwde(itab,jtab)/tfour(itab)
        enddo
      enddo
      do jtab=1,3
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = r1wd(itab)
        enddo 
      enddo

!---------------------------------------------------------------------
!    initialize all derivative table entries.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%td (itab,jtab) = 0.0E+00
          tab1w%td(itab,jtab) = 0.0E+00
          tab2%td  (itab,jtab) = 0.0E+00
          tab3%td (itab,jtab) = 0.0E+00
          tab1%md (itab,jtab) = 0.0E+00
          tab1w%md(itab,jtab) = 0.0E+00
          tab2%md  (itab,jtab) = 0.0E+00
          tab3%md (itab,jtab) = 0.0E+00
          tab1%cd (itab,jtab) = 0.0E+00
          tab1w%cd(itab,jtab) = 0.0E+00
          tab2%cd  (itab,jtab) = 0.0E+00
          tab3%cd (itab,jtab) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for temperature derivatives.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O-1
          tab1%td  (itab,jtab) =    &
          (tab1%vae(itab+1,jtab) - tab1%vae (itab,jtab))/temp_1%tab_inc

          tab1w%td(itab,jtab) =     &
         (tab1w%vae(itab+1,jtab) - tab1w%vae(itab,jtab))/temp_1%tab_inc

          tab2%td  (itab,jtab) =    &
       (tab2%vae  (itab+1,jtab) - tab2%vae  (itab,jtab))/temp_1%tab_inc

          tab3%td (itab,jtab) =     &
         (tab3%vae (itab+1,jtab) - tab3%vae (itab,jtab))/temp_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for mass derivatives.
!---------------------------------------------------------------------
      do jtab=2,NUTABH2O-1
        do itab=1,NTTABH2O
          tab1%md (itab,jtab) =   &
     (tab1%vae (itab,jtab+1) - tab1%vae (itab,jtab))/mass_1%tab_inc

          tab1w%md(itab,jtab) =   &
    (tab1w%vae(itab,jtab+1) - tab1w%vae(itab,jtab))/mass_1%tab_inc

          tab2%md  (itab,jtab) =    &
   (tab2%vae  (itab,jtab+1) - tab2%vae  (itab,jtab))/mass_1%tab_inc

          tab3%md (itab,jtab) =    &
   (tab3%vae (itab,jtab+1) - tab3%vae (itab,jtab))/mass_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!    compute table entries for cross derivatives.
!---------------------------------------------------------------------
      do jtab=2,NUTABH2O-1
        do itab=1,NTTABH2O-1
      tab1%cd (itab,jtab) =    &
             (tab1%vae (itab+1,jtab+1) - tab1%vae (itab+1,jtab) -   &
              tab1%vae (itab  ,jtab+1) + tab1%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

      tab1w%cd(itab,jtab) =    &
             (tab1w%vae(itab+1,jtab+1) - tab1w%vae(itab+1,jtab) -   &
              tab1w%vae(itab  ,jtab+1) + tab1w%vae(itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)


!!  THIS NEVER USED :
!     tab2%cd  (itab,jtab) =     &
!            (tab2%vae  (itab+1,jtab+1) - tab2%vae  (itab+1,jtab) -   &
!             tab2%vae  (itab  ,jtab+1) + tab2%vae  (itab  ,jtab))/   &
!            (DTTABH2O*DUTABH2O)
!            (temp_1%tab_inc*mass_1%tab_inc)


      tab3%cd (itab,jtab) =     &
             (tab3%vae (itab+1,jtab+1) - tab3%vae (itab+1,jtab) -   &
              tab3%vae (itab  ,jtab+1) + tab3%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

        enddo
      enddo
      if (NBTRGE > 0) then
        do m=1,NBTRGE
          do jtab=1,NUTABH2O
            tab1a%vae      (:,jtab,m) = suma(:,jtab,m)/tfour(:)
            tab2a%vae (:,jtab,m) = sumdbea(:,jtab,m)/fortcu(:) 
          enddo
          do jtab=122,NUTABH2O
            tab3a%vae(:,jtab,m) = sum3a(:,jtab,m)/fortcu(:)
          enddo
          do jtab=1,3
            tab1a%vae      (:,jtab,m) = r1a(:,m)
          enddo
          do jtab=1,121
            tab3a%vae(:,jtab,m) = r2a(:,m)/2.0E+00 -     &
                                  s2a(:,m)*zroot(jtab)/3.0E+00 +   &
                                  t3a(:,m)*zmass(jtab)/8.0E+00
          enddo
        enddo

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        tab1a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%td  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab1a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%md  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab1a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab2a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
        tab3a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00

        tab1a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                 (tab1a%vae(2:NTTABH2O,1:NUTABH2O,:) -   &
                  tab1a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab2a%td (1:NTTABH2O-1,1:NUTABH2O,:) =   &
                (tab2a%vae (2:NTTABH2O,1:NUTABH2O,:) -    &
                 tab2a%vae (1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab3a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                  (tab3a%vae(2:NTTABH2O,1:NUTABH2O,:) -  &
                   tab3a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

        tab1a%md(1:NTTABH2O,2:NUTABH2O-1,:) =     &
                  (tab1a%vae(1:NTTABH2O,3:NUTABH2O,:) -   &
                   tab1a%vae(1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab2a%md (1:NTTABH2O,2:NUTABH2O-1,:) =   &
                 (tab2a%vae (1:NTTABH2O,3:NUTABH2O,:) -   &
                  tab2a%vae (1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab3a%md(1:NTTABH2O,2:NUTABH2O-1,:) =     &
                  (tab3a%vae(1:NTTABH2O,3:NUTABH2O,:) -   &
                   tab3a%vae(1:NTTABH2O,2:NUTABH2O-1,:))/mass_1%tab_inc

        tab1a%cd(1:NTTABH2O-1,2:NUTABH2O-1,:) =     &
                        (tab1a%vae(2:NTTABH2O,3:NUTABH2O,:) -    &
                         tab1a%vae(2:NTTABH2O,2:NUTABH2O-1,:)   -  &
                         tab1a%vae(1:NTTABH2O-1,3:NUTABH2O,:)   +  &
                         tab1a%vae(1:NTTABH2O-1,2:NUTABH2O-1,:))/  &
                                         (temp_1%tab_inc*mass_1%tab_inc)

        tab3a%cd(1:NTTABH2O-1,2:NUTABH2O-1,:) =    &
                        (tab3a%vae(2:NTTABH2O,3:NUTABH2O,:) -    &
                         tab3a%vae(2:NTTABH2O,2:NUTABH2O-1,:)   -  &
                         tab3a%vae(1:NTTABH2O-1,3:NUTABH2O,:)   +  &
                         tab3a%vae(1:NTTABH2O-1,2:NUTABH2O-1,:))/  &
                                        (temp_1%tab_inc*mass_1%tab_inc)
     
!---------------------------------------------------------------------
!    deallocate local arrays.
!---------------------------------------------------------------------
        deallocate ( r1a    )
        deallocate (r2a     )
        deallocate (s2a     )
        deallocate ( t3a     )
        deallocate ( suma    )
        deallocate ( sumdbea )
        deallocate ( sum3a   )
        deallocate ( sum4a   )
        deallocate ( sum6a   )
        deallocate ( sum7a   )
        deallocate (sum8a  )
      endif

!-------------------------------------------------------------------


end subroutine table



!####################################################################


                  end module longwave_tables_mod


