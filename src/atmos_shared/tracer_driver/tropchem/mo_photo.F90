      module MO_PHOTO_MOD
!----------------------------------------------------------------------
!        ... Photolysis interp table and related arrays
!----------------------------------------------------------------------
      use fms_mod,          only : mpp_clock_begin, mpp_clock_id, &
                                   mpp_clock_end, CLOCK_MODULE
      use mpp_mod,          only : mpp_error, FATAL
      use mpp_io_mod,       only : mpp_open, MPP_RDONLY, MPP_ASCII,MPP_MULTI, &
                                   MPP_SINGLE, mpp_close
      use time_manager_mod, only : time_type, get_date
      use constants_mod,    only : PI

      implicit none

      private
      public :: prate_init, photo, set_ub_col, setcol, sundis

      save

      integer, parameter :: jdim     = 40
      integer, parameter :: altdim   = 46
      integer, parameter :: zangdim  = 11
      integer, parameter :: o3ratdim = 7
      integer, parameter :: albdim   = 4
      integer, parameter :: t500dim  = 3
      integer, parameter :: t200dim  = 2
      integer, parameter :: tabdim   = jdim*altdim*zangdim*o3ratdim*albdim*t500dim*t200dim

      integer ::  offset(7)
      integer ::  indexer(jdim)
      integer ::  jno_ndx, jpooh_ndx, jc2h5ooh_ndx, jc3h7ooh_ndx, jrooh_ndx, &
                  jch3co3h_ndx, jmpan_ndx, jmacr_a_ndx, jmacr_b_ndx, jonitr_ndx, &
                  jxooh_ndx, jisopooh_ndx, jglyald_ndx, jhyac_ndx, jch3ooh_ndx, &
                  jh2o2_ndx, jpan_ndx, jch3cho_ndx, jho2no2_ndx, &
                  jn2o5_ndx, jo3p_ndx, jno2_ndx, jno3_ndx, &
                  jclono2_ndx, jhocl_ndx, jcl2o2_ndx, jbrono2_ndx, jhobr_ndx, &
                  jbrcl_ndx, jbro_ndx, jcl2_ndx, jh2o_ndx, jn2o_ndx, jhno3_ndx
      integer ::  ox_ndx, o3_ndx
      real    ::  ajl(jdim,altdim,zangdim,o3ratdim,albdim,t500dim,t200dim) = 0., &
                  ajl_solarmin(jdim,altdim,zangdim,o3ratdim,albdim,t500dim,t200dim) = 0.
!RSH  real    ::  ajl2(2,2,2,2,2,2)
      real    ::  vo3(0:80), vo3_solarmin(0:80)
      real    ::  delvo3(0:79)
      real    ::  delz(altdim-1)
      real    ::  delang(zangdim-1)
      real    ::  delv(o3ratdim-1)
      real    ::  delalb(albdim-1)
      real    ::  delt500(t500dim-1)
      real    ::  delt200(t200dim-1)

      real, parameter :: zz(altdim) = &
        (/  0.,  1.,  2.,  3.,  4.,  5.,  6.,  8., 10., 12., 14., 16., 18., 20., 22., &
           24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 52., &
           54., 56., 58., 60., 62., 64., 66., 68., 70., 72., 74., 76., 78., 80., 82., &
           85. /)
      real, parameter :: vcos(zangdim) = &
        (/ -0.07, -0.05, -0.01, 0.01, 0.05, 0.1, 0.2,   0.4,   0.6,  0.8,  1.0 /)
      real, parameter :: xv3(o3ratdim) = (/ .5, .75, 1., 1.25, 1.5, 2., 5. /)
      real, parameter :: albev(albdim) = (/ .05, .2, .5, 1. /)
      real, parameter :: t500(t500dim) = (/ 228., 248., 268. /)
      real, parameter :: t200(t200dim) = (/ 205., 225. /)
      real, parameter :: coszen_min = vcos(1)


      integer, parameter :: &
         TAB_NDX_JO2        = 1,  TAB_NDX_JO1D       = 2, &
         TAB_NDX_JO3P       = 3,  TAB_NDX_JNO2       = 4, &
         TAB_NDX_JNO3       = 5,  TAB_NDX_JN2O5      = 6, &
         TAB_NDX_JN2O5_225  = 7,  TAB_NDX_JN2O5_250  = 8, &
         TAB_NDX_JN2O5_300  = 9,  TAB_NDX_JN2O       = 10, &
         TAB_NDX_JN2O_200   = 11, TAB_NDX_JN2O_250   = 12, &
         TAB_NDX_JN2O_300   = 13, TAB_NDX_JH2O2      = 14, &
         TAB_NDX_JHNO3      = 15, TAB_NDX_JHNO3_200  = 16, &
         TAB_NDX_JHNO3_250  = 17, TAB_NDX_JHNO3_300  = 18, &
         TAB_NDX_JHO2NO2    = 19, TAB_NDX_JCH2Oa     = 20, &
         TAB_NDX_JCH2Ob     = 21, TAB_NDX_JCH3CHO    = 22, &
         TAB_NDX_JMGLY      = 23, TAB_NDX_JACET      = 24, &
         TAB_NDX_JCH3OOH    = 25, TAB_NDX_JPAN       = 26, &
         TAB_NDX_JCLONO2    = 27, TAB_NDX_JCLONO2_200= 28, &
         TAB_NDX_JCLONO2_250= 29, TAB_NDX_JCLONO2_300= 30, &
         TAB_NDX_JBRONO2    = 31, TAB_NDX_JCL2       = 32, &
         TAB_NDX_JMVK       = 33, TAB_NDX_JMACRa     = 34, &
         TAB_NDX_JCL2O2     = 35, TAB_NDX_JHYAC      = 36, &
         TAB_NDX_JHOBR      = 37, TAB_NDX_JBR2       = 38, &
         TAB_NDX_JHOCL      = 39, TAB_NDX_JBRCL      = 40

      logical :: use_tdep_jvals, use_solar_cycle
      real    :: o3_column_top, jno_scale_factor

character(len=128), parameter :: version     = '$Id: mo_photo.F90,v 19.0 2012/01/06 20:34:00 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.
integer                       :: photo_clock

      CONTAINS
      
! <SUBROUTINE NAME="prate_init">
!   <OVERVIEW>
!     Initialize photolysis rate lookup table calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the calculation of photolysis rates
!     from the TUV lookup table
!   </DESCRIPTION>
!   <TEMPLATE>
!     call prate_init( filename, filename_solarmin, lpath, mspath, use_tdep_jvals_in )
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character(len=*)">
!     Filename for ASCII lookup table file (if solar cycle used, this is the filename
!     for solar maximum conditions).
!   </IN>
!   <IN NAME="filename_solarmin" TYPE="character(len=*)">
!     Filename for ASCII lookup table file for solar minimum conditions
!     (if solar cycle used)
!   </IN>
!   <IN NAME="lpath" TYPE="character(len=*)">
!     Local directory path for input files
!   </IN>
!   <IN NAME="mspath" TYPE="character(len=*)">
!     Remote directory path for input files (not used)
!   </IN>
!   <IN NAME="use_tdep_jvals_in" TYPE="logical">
!     Does the j-value lookup table contain temperature-dependent photolysis rates?
!   </IN>
!   <IN NAME="o3_column_top_in" TYPE="real">
!     Ozone column above model top (DU)
!   </IN>
!   <IN NAME="jno_scale_factor_in" TYPE="real">
!     Scale factor for NO photolysis rate (jNO)
!   </IN>
      subroutine prate_init( filename, filename_solarmin, lpath, mspath, &
                             use_tdep_jvals_in, o3_column_top_in,  &
                             jno_scale_factor_in, retain_cm3_bugs )
!----------------------------------------------------------------------
!     ... Read in the photorate tables and arrays
!         Results are "returned" via the common block photo_tables
!          This is for the new, expanded chemistry (11/21/94)
!----------------------------------------------------------------------
        
      use mo_chem_utls_mod,  only : get_spc_ndx, get_rxt_ndx
      implicit none

!----------------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------------
      character(len=*), intent(in) :: filename, filename_solarmin, lpath, mspath
      logical,          intent(in) :: use_tdep_jvals_in
      real,             intent(in) :: o3_column_top_in, &
                                      jno_scale_factor_in
      logical,          intent(in) :: retain_cm3_bugs

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer    :: it500, it200, izen, ialb, idob
      integer    :: ios
      integer    :: unit
!     integer    :: retval
!     logical    :: cosb
!     real       :: temp(tabdim)
      character(len=128) :: msg

!     unit = navu()
!----------------------------------------------------------------------
!        ... Only masternode gets photorate table
!----------------------------------------------------------------------
!     if( masternode ) then
!        retval = ATTACH( TRIM(lpath) // TRIM( filename ), &
!                         TRIM( mspath ) // TRIM( filename ), &
!                         unit, &
!                         .false., &                 ! non binary dataset
!                         cosb )
!        if( retval /= 0 ) then
!           write(*,*) 'PRATE_INIT: Failure opening file ',TRIM(lpath) // TRIM( filename )
!           write(*,*) 'Error code = ',retval
!           call ENDRUN
!        end if
!        CLOSE( unit )
!     end if
#ifdef USE_MPI
!----------------------------------------------------------------------
!        ... All compute nodes wait for masternode to acquire file
!----------------------------------------------------------------------
!     call MPI_BARRIER( mpi_comm_comp, ios )
!     if( ios /= MPI_SUCCESS ) then
!        write(*,*) 'PRATE_INIT: Mpi barrier failed; error = ',ios
!        call ENDRUN
!     end if
#endif
!----------------------------------------------------------------------
!        ... all compute nodes open file
!----------------------------------------------------------------------
!     OPEN( unit   = unit, &
!           file   = TRIM(lpath) // TRIM(filename), &
!           status = 'old', &
!           form   = 'formatted', &
!           recl   = 4500, &
!          iostat = ios )
!     if( ios /= 0 ) then
!----------------------------------------------------------------------
!        ... Open error exit
!----------------------------------------------------------------------
!        write(6,'('' PRATE_INIT : Error ('',i5,'') opening file '',a)') &
!           ios, TRIM(lpath) // TRIM(filename)
!        call ENDRUN
!     end if

!----------------------------------------------------------------------
!            ... open file using mpp_open
!----------------------------------------------------------------------

      call mpp_open( unit, trim(lpath)//trim(filename), MPP_RDONLY, MPP_ASCII, &
                     threading = MPP_MULTI, fileset = MPP_SINGLE, &
                     recl = 4500)

!----------------------------------------------------------------------
!        ... Readin the reference o3 column and photorate table
!----------------------------------------------------------------------
      read(unit,*,iostat=ios) vo3
      if( ios /= 0 ) then
         msg = ' PRATE_INIT: Failed to read o3 column'
         call ENDRUN(msg)
      end if

      do it500 = 1,t500dim
         do it200 = 1,t200dim
            do izen = 1,zangdim
               do ialb = 1,albdim
                  do idob = 1,o3ratdim
                     read(unit,*,iostat=ios) ajl(:,:,izen,idob,ialb,it500,it200)
                     if( ios /= 0 ) then
                        msg = ' PRATE_INIT: Failed to read photo table; error = '//char(ios)
                        call ENDRUN(msg)
                     end if
                  end do
               end do
            end do
         end do
      end do

!----------------------------------------------------------------------
!        ... Set module variables
!----------------------------------------------------------------------
      delz(:altdim-1) = 1. / (zz(2:altdim) - zz(:altdim-1))
      delvo3(0:79) = vo3(1:80) - vo3(0:79)
      delang(:zangdim-1)  = 1. / (vcos(2:zangdim) - vcos(:zangdim-1))
      delv(:o3ratdim-1)   = 1. / (xv3(2:o3ratdim) - xv3(:o3ratdim-1))
      delalb(:albdim-1)   = 1. / (albev(2:albdim) - albev(:albdim-1))
      delt500(:t500dim-1) = 1. / (t500(2:t500dim) - t500(:t500dim-1))
      delt200(:t200dim-1) = 1. / (t200(2:t200dim) - t200(:t200dim-1))

      offset(1) = jdim
      offset(2) = offset(1)*altdim
      offset(3) = offset(2)*zangdim
      offset(4) = offset(3)*o3ratdim
      offset(5) = offset(4)*albdim
      offset(6) = offset(5)*t500dim
      offset(7) = SUM( offset(1:6) )

!     close( unit )
      call mpp_close( unit )

!-----------------------------------------------------------------
!           ... check whether using solar cycle
!-----------------------------------------------------------------
      if (filename_solarmin == '' .or. filename_solarmin == filename) then
         use_solar_cycle = .false.
      else
         use_solar_cycle = .true.

!----------------------------------------------------------------------
!            ... open file using mpp_open
!----------------------------------------------------------------------
         call mpp_open( unit, trim(lpath)//trim(filename_solarmin), MPP_RDONLY, MPP_ASCII, &
                        threading = MPP_MULTI, fileset = MPP_SINGLE, &
                        recl = 4500)

!----------------------------------------------------------------------
!        ... Readin the reference o3 column and photorate table 
!            for solar minimum
!----------------------------------------------------------------------
         read(unit,*,iostat=ios) vo3_solarmin
         if( ios /= 0 ) then
            msg = ' PRATE_INIT: Failed to read solarmin o3 column'
            call ENDRUN(msg)
         end if

         do it500 = 1,t500dim
         do it200 = 1,t200dim
         do izen = 1,zangdim
         do ialb = 1,albdim
         do idob = 1,o3ratdim
            read(unit,*,iostat=ios) ajl_solarmin(:,:,izen,idob,ialb,it500,it200)
            if( ios /= 0 ) then
               msg = ' PRATE_INIT: Failed to read solarmin photo table; error = '//char(ios)
               call ENDRUN(msg)
            end if
         end do
         end do
         end do
         end do
         end do
         call mpp_close( unit )
      end if

!-----------------------------------------------------------------
!           ... setup mapping array, indexer, from table to model
!-----------------------------------------------------------------
      indexer(TAB_NDX_JO2)      = get_rxt_ndx( 'jo2' )
      indexer(TAB_NDX_JO1D)     = get_rxt_ndx( 'jo1d' )
      indexer(TAB_NDX_JO3P)     = get_rxt_ndx( 'jo3p' )
      indexer(TAB_NDX_JNO2)     = get_rxt_ndx( 'jno2' )
      indexer(TAB_NDX_JNO3)     = get_rxt_ndx( 'jno3' )
      indexer(TAB_NDX_JN2O5)    = get_rxt_ndx( 'jn2o5' )
      indexer(TAB_NDX_JN2O5_225)= 0
      indexer(TAB_NDX_JN2O5_250)= 0
      indexer(TAB_NDX_JN2O5_300)= 0
      indexer(TAB_NDX_JN2O)     = get_rxt_ndx( 'jn2o' )
      indexer(TAB_NDX_JN2O_200) = 0
      indexer(TAB_NDX_JN2O_250) = 0
      indexer(TAB_NDX_JN2O_300) = 0
      indexer(TAB_NDX_JH2O2)    = get_rxt_ndx( 'jh2o2' )
      indexer(TAB_NDX_JHNO3)    = get_rxt_ndx( 'jhno3' )
      indexer(TAB_NDX_JHNO3_200)= 0
      indexer(TAB_NDX_JHNO3_250)= 0
      indexer(TAB_NDX_JHNO3_300)= 0
      indexer(TAB_NDX_JHO2NO2)  = get_rxt_ndx( 'jho2no2' )
      indexer(TAB_NDX_JCH2Oa)   = get_rxt_ndx( 'jch2o_a' )
      indexer(TAB_NDX_JCH2Ob)   = get_rxt_ndx( 'jch2o_b' )
      indexer(TAB_NDX_JCH3CHO)  = get_rxt_ndx( 'jch3cho' )
      indexer(TAB_NDX_JMGLY)    = get_rxt_ndx( 'jmgly' )
      indexer(TAB_NDX_JACET)    = get_rxt_ndx( 'jacet' )
      indexer(TAB_NDX_JCH3OOH)  = get_rxt_ndx( 'jch3ooh' )
      indexer(TAB_NDX_JPAN)     = get_rxt_ndx( 'jpan' )
      indexer(TAB_NDX_JCLONO2)  = get_rxt_ndx( 'jclono2' )
      indexer(TAB_NDX_JCLONO2_200)=0
      indexer(TAB_NDX_JCLONO2_250)=0
      indexer(TAB_NDX_JCLONO2_300)=0
      indexer(TAB_NDX_JBRONO2)  = get_rxt_ndx( 'jbrono2' )
      indexer(TAB_NDX_JCL2)     = get_rxt_ndx( 'jcl2' )
      indexer(TAB_NDX_JMVK)     = get_rxt_ndx( 'jmvk' )
      indexer(TAB_NDX_JMACRa)   = get_rxt_ndx( 'jmacr_a' )
      indexer(TAB_NDX_JCL2O2)   = get_rxt_ndx( 'jcl2o2' )
      indexer(TAB_NDX_JHYAC)    = get_rxt_ndx( 'jhyac' )
      indexer(TAB_NDX_JHOBR)    = get_rxt_ndx( 'jhobr' )
      indexer(TAB_NDX_JBR2)     = get_rxt_ndx( 'jbr2' )
      indexer(TAB_NDX_JHOCL)    = get_rxt_ndx( 'jhocl' )
      indexer(TAB_NDX_JBRCL)    = get_rxt_ndx( 'jbrcl' )

      jno_ndx      = get_rxt_ndx( 'jno' )
      jpooh_ndx    = get_rxt_ndx( 'jpooh' )
      jc2h5ooh_ndx = get_rxt_ndx( 'jc2h5ooh' )
      jc3h7ooh_ndx = get_rxt_ndx( 'jc3h7ooh' )
      jrooh_ndx    = get_rxt_ndx( 'jrooh' )
      jch3co3h_ndx = get_rxt_ndx( 'jch3co3h' )
      jmpan_ndx    = get_rxt_ndx( 'jmpan' )
      jmacr_a_ndx  = get_rxt_ndx( 'jmacr_a' )
      jmacr_b_ndx  = get_rxt_ndx( 'jmacr_b' )
      jonitr_ndx   = get_rxt_ndx( 'jonitr' )
      jxooh_ndx    = get_rxt_ndx( 'jxooh' )
      jisopooh_ndx = get_rxt_ndx( 'jisopooh' )
      jglyald_ndx  = get_rxt_ndx( 'jglyald' )
      jhyac_ndx    = get_rxt_ndx( 'jhyac' )
      jch3ooh_ndx  = get_rxt_ndx( 'jch3ooh' )
      jh2o2_ndx    = get_rxt_ndx( 'jh2o2' )
      jpan_ndx     = get_rxt_ndx( 'jpan' )
      jch3cho_ndx  = get_rxt_ndx( 'jch3cho' )
      jho2no2_ndx  = get_rxt_ndx( 'jho2no2' )
      jn2o5_ndx    = get_rxt_ndx( 'jn2o5' )
      jo3p_ndx     = get_rxt_ndx( 'jo3p' )
      jno2_ndx     = get_rxt_ndx( 'jno2' )
      jno3_ndx     = get_rxt_ndx( 'jno3' )
      jclono2_ndx  = get_rxt_ndx( 'jclono2' )
      jhocl_ndx    = get_rxt_ndx( 'jhocl' )
      jcl2o2_ndx   = get_rxt_ndx( 'jcl2o2' )
      jbrono2_ndx  = get_rxt_ndx( 'jbrono2' )
      jhobr_ndx    = get_rxt_ndx( 'jhobr' )
      jbrcl_ndx    = get_rxt_ndx( 'jbrcl' )
      jbro_ndx     = get_rxt_ndx( 'jbro' )
      jcl2_ndx     = get_rxt_ndx( 'jcl2' )
      jh2o_ndx     = get_rxt_ndx( 'jh2o' )
      jn2o_ndx     = get_rxt_ndx( 'jn2o' )
      jhno3_ndx    = get_rxt_ndx( 'jhno3' )

!      ox_ndx = get_spc_ndx( 'OX' )
!     for ox budget (jmao,1/7/2011)
   if (retain_cm3_bugs) then
      ox_ndx = get_spc_ndx( 'OX' )
      o3_ndx = get_spc_ndx( 'O3' )
   else
      ox_ndx = get_spc_ndx( 'O3' )
      o3_ndx = get_spc_ndx( 'O3' )
   endif

      use_tdep_jvals   = use_tdep_jvals_in
      o3_column_top    = o3_column_top_in
      jno_scale_factor = jno_scale_factor_in

      photo_clock = mpp_clock_id ('Tropchem:Photo', grain=CLOCK_MODULE)

      end subroutine prate_init
! </SUBROUTINE>


! <SUBROUTINE NAME="PHOTO">
!   <OVERVIEW>
!     Calculate photolysis rates
!   </OVERVIEW>
!   <DESCRIPTION>
!     Calculate photolysis rates from the TUV lookup table
!   </DESCRIPTION>
!   <TEMPLATE>
!     call PHOTO( photos, pmid, pdel, temper, zmid, col_dens, coszen,  & 
!                 srf_alb, lwc, clouds, esfact, solar_phase, plonl )
!   </TEMPLATE>
!   <IN NAME="photos" TYPE="real" DIM="(:,:,:)">
!     Photodissociation rates (s^-1)
!   </IN>
!   <IN NAME="pmid" TYPE="real" DIM="(:,:)">
!     Full level pressures (Pa)
!   </IN>
!   <IN NAME="pdel" TYPE="real" DIM="(:,:)">
!     Half level (interface) pressures (Pa)
!   </IN>
!   <IN NAME="temper" TYPE="real" DIM="(:,:)">
!     Full level temperatures (K)
!   </IN>
!   <IN NAME="zmid" TYPE="real" DIM="(:,:)">
!     Full level absolute geopotential altitudes (km)
!   </IN>
!   <IN NAME="col_dens" TYPE="real" DIM="(:,:,:)">
!     Column densities
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:)">
!     Cosine of solar zenith angle
!   </IN>
!   <IN NAME="srf_alb" TYPE="real" DIM="(:)">
!     Surface albedo
!   </IN>
!   <IN NAME="lwc" TYPE="real" DIM="(:,:)">
!     Cloud liquid water content (kg/kg)
!   </IN>
!   <IN NAME="clouds" TYPE="real" DIM="(:,:)">
!     Cloud fraction
!   </IN>
!   <IN NAME="esfact" TYPE="real">
!     Earth-sun distance factor
!   </IN>
!   <IN NAME="solar_phase" TYPE="real">
!     Solar cycle phase (1=max, 0=min)
!   </IN>
!   <IN NAME="plonl" TYPE="integer">
!     Size of longitude dimension
!   </IN>
      subroutine PHOTO( photos, pmid, pdel, temper, zmid, &
                        col_dens, &
                        coszen,  & 
                        srf_alb, lwc, clouds, &
                        esfact, solar_phase, &
                        plonl )

      use CHEM_MODS_MOD, only : ncol_abs, phtcnt
!     use M_RXT_ID_MOD

      implicit none

!-----------------------------------------------------------------
!           ... Dummy arguments
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) ::   esfact, &                 ! earth sun distance factor
                            solar_phase               ! solar cycle phase (1=max, 0=min)
      real, intent(in) ::   col_dens(:,:,:), &        ! column densities
                            coszen(:), &              ! solar zenith angle
                            srf_alb(:), &             ! surface albedo
                            lwc(:,:), &               ! liquid water content (mass mr)
                            clouds(:,:), &            ! cloud fraction
                            pmid(:,:), &              ! midpoint pressure in pascals
                            pdel(:,:), &              ! del pressure about midpoint in pascals
                            zmid(:,:), &              ! midpoint height
                            temper(:,:)               ! midpoint temperature
      real, intent(out) ::  photos(:,:,:)             ! photodissociation rates

!-----------------------------------------------------------------
!            ... Local variables
!-----------------------------------------------------------------
      integer  ::  i, k, m                 ! indicies
      integer  ::  plev
      logical  ::  zagtz(size(coszen))     ! zenith angle > 0 flag array
      real     ::  t500, t200              ! 500 & 200 mb temperatures
      real, dimension(size(zmid,2)) :: &
                   fac1, &                ! work space for J(no) calc
                   fac2, &                ! work space for J(no) calc
                   colo3, &               ! vertical o3 column density
                   zarg, &                ! vertical height array
                   pline, &               ! vertical pressure array
                   tline, &               ! vertical temperature array
                   cld_line, &            ! vertical cloud array
                   lwc_line, &            ! vertical lwc array
                   eff_alb, &             ! effective albedo from cloud modifications
                   cld_mult               ! clould multiplier
      real, dimension(plonl,size(zmid,2)) :: &
                   tmp, &                        ! wrk array
                   tmp_jch3ooh, &                ! wrk array
                   tmp_jpan, &                   ! wrk array
                   tmp_jh2o2, &                  ! wrk array
                   tmp_jch3cho, &                ! wrk array
                   tmp_jmacr_a, &                ! wrk array
                   tmp_jn2o_200, &               ! wrk array
                   tmp_jn2o_250, &               ! wrk array
                   tmp_jn2o_300, &               ! wrk array
                   tmp_jn2o5_225, &              ! wrk array
                   tmp_jn2o5_250, &              ! wrk array
                   tmp_jn2o5_300, &              ! wrk array
                   tmp_jhno3_200, &              ! wrk array
                   tmp_jhno3_250, &              ! wrk array
                   tmp_jhno3_300, &              ! wrk array
                   tmp_jclono2_200, &            ! wrk array
                   tmp_jclono2_250, &            ! wrk array
                   tmp_jclono2_300, &            ! wrk array
                   wgt200, wgt225, wgt250, wgt300, &     ! wrk array
                   tmp_jno
      real    ::   prates(jdim,size(zmid,2))        ! photorates matrix

      call mpp_clock_begin (photo_clock)
      plev = SIZE(zmid,2)
!-----------------------------------------------------------------
!        ... Zero all photorates
!-----------------------------------------------------------------
      do m = 1,max(1,phtcnt)
         do k = 1,plev
            photos(:,k,m) = 0.
         end do
      end do
      do k = 1,plev
         tmp_jch3ooh(:,k)     = 0.
         tmp_jpan(:,k)        = 0.
         tmp_jh2o2(:,k)       = 0.
         tmp_jch3cho(:,k)     = 0.
         tmp_jmacr_a(:,k)     = 0.
         tmp_jno(:,k)         = 0.
         tmp_jn2o_200(:,k)    = 0.
         tmp_jn2o_250(:,k)    = 0.
         tmp_jn2o_300(:,k)    = 0.
         tmp_jn2o5_225(:,k)   = 0.
         tmp_jn2o5_250(:,k)   = 0.
         tmp_jn2o5_300(:,k)   = 0.
         tmp_jhno3_200(:,k)   = 0.
         tmp_jhno3_250(:,k)   = 0.
         tmp_jhno3_300(:,k)   = 0.
         tmp_jclono2_200(:,k) = 0.
         tmp_jclono2_250(:,k) = 0.
         tmp_jclono2_300(:,k) = 0.
      end do
      zagtz(:) = coszen(:) >= coszen_min

      do i = 1,plonl
         if( zagtz(i) ) then
!           secant = 1. / coszen(i)
!           if( secant <= 50. ) then
!           if( coszen(i) >= coszen_min ) then
            zarg(:)     = zmid(i,:)
            colo3(:)    = col_dens(i,:,1)
            pline(:)    = pmid(i,:)
            fac1(:)     = pdel(i,:)
            tline(:)    = temper(i,:)
            lwc_line(:) = lwc(i,:)
            cld_line(:) = clouds(i,:)
            call cloud_mod( coszen(i), cld_line, lwc_line, fac1, srf_alb(i), &
                            eff_alb, cld_mult )
            call T_INT( pline, tline, t500, t200 )
            call PHOTO_INTERP( zarg, coszen(i), colo3, eff_alb, t500, &
                               t200, solar_phase, prates )
            do m = 1,jdim
               if( indexer(m) > 0 ) then
                  photos(i,:,indexer(m)) = esfact *prates(m,:) * cld_mult(:)
               else
                  select case( m )
                     case( TAB_NDX_JCH3OOH )
                        tmp_jch3ooh(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JH2O2 )
                        tmp_jh2o2(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JCH3CHO )
                        tmp_jch3cho(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JPAN )
                        tmp_jpan(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JMACRa )
                        tmp_jmacr_a(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O_200 )
                        tmp_jn2o_200(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O_250 )
                        tmp_jn2o_250(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O_300 )
                        tmp_jn2o_300(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O5_225 )
                        tmp_jn2o5_225(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O5_250 )
                        tmp_jn2o5_250(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JN2O5_300 )
                        tmp_jn2o5_300(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JHNO3_200 )
                        tmp_jhno3_200(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JHNO3_250 )
                        tmp_jhno3_250(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JHNO3_300 )
                        tmp_jhno3_300(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JCLONO2_200 )
                        tmp_jclono2_200(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JCLONO2_250 )
                        tmp_jclono2_250(i,:) = esfact *prates(m,:) * cld_mult(:)
                     case( TAB_NDX_JCLONO2_300 )
                        tmp_jclono2_300(i,:) = esfact *prates(m,:) * cld_mult(:)
                  end select
               end if
            end do
!-----------------------------------------------------------------
!        ... Calculate J(no) from formula
!-----------------------------------------------------------------
            if (coszen(i) > 0. ) then
               fac1(:) = 1.e-8  * (col_dens(i,:,2)/coszen(i))**.38
               fac2(:) = 5.e-19 * col_dens(i,:,1) / coszen(i)
               if( jno_ndx > 0 ) then
                  photos(i,:,jno_ndx) = 4.5e-6 * esfact * exp( -(fac1(:) + fac2(:)) ) &
                                               * cld_mult(:) * jno_scale_factor
               else
                  tmp_jno(i,:) = 4.5e-6 * esfact * exp( -(fac1(:) + fac2(:)) ) &
                                        * cld_mult(:) * jno_scale_factor
               end if
            end if
!-----------------------------------------------------------------
!        ... ho2no2 near-IR photolysis
!-----------------------------------------------------------------
            if (jho2no2_ndx > 0) then
               photos(i,:,jho2no2_ndx) = photos(i, :, jho2no2_ndx) + &
                    1.e-5 * esfact * cld_mult(:)
            endif
!        end if
         end if
      end do

!-----------------------------------------------------------------
!        ... Set J(pooh) from J(ch3ooh)
!                J(c2h5ooh) from J(ch3ooh)
!                J(c3h7ooh) from J(ch3ooh)
!                J(rooh) from J(ch3ooh)
!                J(ch3coooh) = .28 * J(h2o2)
!                J(mpan) from J(pan)
!                J(macr_a) and J(macr_b) = 1/2 * J(macr_total)
!               J(onitr) from j(ch3cho)
!               J(xooh) from J(ch3ooh)
!               J(isopooh) from J(ch3ooh)
!                J(glyald) = 3 * J(ch3cho)
!               J(hyac) from J(ch3cho)
!-----------------------------------------------------------------
      if( jch3ooh_ndx > 0 ) then
         tmp(:,:) = photos(:,:,jch3ooh_ndx)
      else
         tmp(:,:) = tmp_jch3ooh(:,:)
      end if
      if( jpooh_ndx > 0 ) then
         photos(:,:,jpooh_ndx)    = tmp(:,:)
      end if
      if( jc2h5ooh_ndx > 0 ) then
         photos(:,:,jc2h5ooh_ndx) = tmp(:,:)
      end if
      if( jc3h7ooh_ndx > 0 ) then
         photos(:,:,jc3h7ooh_ndx) = tmp(:,:)
      end if
      if( jrooh_ndx > 0 ) then
         photos(:,:,jrooh_ndx)    = tmp(:,:)
      end if
      if( jxooh_ndx > 0 ) then
        photos(:,:,jxooh_ndx)    = tmp(:,:)
      end if
      if( jisopooh_ndx > 0 ) then
         photos(:,:,jisopooh_ndx) = tmp(:,:)
      end if
      if( jch3co3h_ndx > 0 ) then
         if( jh2o2_ndx > 0 ) then
            photos(:,:,jch3co3h_ndx) = .28 * photos(:,:,jh2o2_ndx)
         else
            photos(:,:,jch3co3h_ndx) = .28 * tmp_jh2o2(:,:)
         end if
      end if
      if( jmpan_ndx > 0 ) then
         if( jpan_ndx > 0 ) then
            photos(:,:,jmpan_ndx)    = photos(:,:,jpan_ndx)
         else
            photos(:,:,jmpan_ndx)    = tmp_jpan(:,:)
         end if
      end if
      if( jmacr_a_ndx > 0 ) then
         photos(:,:,jmacr_a_ndx)  = .5 * photos(:,:,jmacr_a_ndx)
      end if
      if( jmacr_b_ndx > 0 ) then
         if( jmacr_a_ndx > 0 ) then
            photos(:,:,jmacr_b_ndx)  = photos(:,:,jmacr_a_ndx)
         else
            photos(:,:,jmacr_b_ndx)  = .5 * tmp_jmacr_a(:,:)
         end if
      end if
      if( jonitr_ndx > 0 ) then
         if( jch3cho_ndx > 0 ) then
            photos(:,:,jonitr_ndx)   = photos(:,:,jch3cho_ndx)
         else
            photos(:,:,jonitr_ndx)   = tmp_jch3cho(:,:)
         end if
      end if
      if( jglyald_ndx > 0 ) then
         if( jch3cho_ndx > 0 ) then
            photos(:,:,jglyald_ndx)  = 3. * photos(:,:,jch3cho_ndx)
         else
            photos(:,:,jglyald_ndx)   = 3. *tmp_jch3cho(:,:)
         end if
      end if
      if( jh2o_ndx > 0 ) then
         if( jno_ndx > 0 ) then
            photos(:,:,jh2o_ndx) = 0.1*photos(:,:,jno_ndx)
         else
            photos(:,:,jh2o_ndx) = 0.1*tmp_jno(:,:)
         end if
      end if
      if( jn2o_ndx > 0 .and. use_tdep_jvals ) then
         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
         photos(:,:,jn2o_ndx)    = wgt200(:,:)*tmp_jn2o_200(:,:) + &
                                   wgt250(:,:)*tmp_jn2o_250(:,:) + &
                                   wgt300(:,:)*tmp_jn2o_300(:,:)
      end if
      if( jn2o5_ndx > 0 .and. use_tdep_jvals ) then
         wgt225(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/25. ) )
         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
         wgt250(:,:)  = 1. - wgt225(:,:) - wgt300(:,:)
         photos(:,:,jn2o5_ndx)   = wgt225(:,:)*tmp_jn2o5_225(:,:) + &
                                   wgt250(:,:)*tmp_jn2o5_250(:,:) + &
                                   wgt300(:,:)*tmp_jn2o5_300(:,:)
      end if
      if( jhno3_ndx > 0 .and. use_tdep_jvals ) then
         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
         photos(:,:,jhno3_ndx)   = wgt200(:,:)*tmp_jhno3_200(:,:) + &
                                   wgt250(:,:)*tmp_jhno3_250(:,:) + &
                                   wgt300(:,:)*tmp_jhno3_300(:,:)
      end if
      if( jclono2_ndx > 0 .and. use_tdep_jvals ) then
         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
         photos(:,:,jclono2_ndx) = wgt200(:,:)*tmp_jclono2_200(:,:) + &
                                   wgt250(:,:)*tmp_jclono2_250(:,:) + &
                                   wgt300(:,:)*tmp_jclono2_300(:,:)
      end if

      call mpp_clock_end   (photo_clock)
      end subroutine PHOTO
! </SUBROUTINE>

      subroutine cloud_mod( coszen, clouds, lwc, delp, srf_alb, &
                            eff_alb, cld_mult )
!-----------------------------------------------------------------------
!         ... Cloud alteration factors for photorates and albedo
!-----------------------------------------------------------------------


      implicit none

      real, parameter :: gi = 1./9.80616

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      real, intent(in)    ::  coszen             ! cosine of zenith angle
      real, intent(in)    ::  srf_alb            ! surface albedo
      real, intent(in)    ::  clouds(:)          ! cloud fraction
      real, intent(in)    ::  lwc(:)             ! liquid water content (mass mr)
      real, intent(in)    ::  delp(:)            ! del press about midpoint in pascals
      real, intent(out)   ::  eff_alb(:)         ! effective albedo
      real, intent(out)   ::  cld_mult(:)        ! photolysis mult factor

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      integer :: k
      integer :: plev, plevm
      real    :: coschi
      real    :: del_lwp(SIZE(clouds,1))
      real    :: del_tau(SIZE(clouds,1))
      real    :: above_tau(SIZE(clouds,1))
      real    :: below_tau(SIZE(clouds,1))
      real    :: above_cld(SIZE(clouds,1))
      real    :: below_cld(SIZE(clouds,1))
      real    :: above_tra(SIZE(clouds,1))
      real    :: below_tra(SIZE(clouds,1))
      real    :: fac1(SIZE(clouds,1))
      real    :: fac2(SIZE(clouds,1))
      
      plev = SIZE(clouds,1)
      plevm = plev-1
!---------------------------------------------------------
!        ... Modify lwc for cloud fraction and form
!            liquid water path for each layer
!---------------------------------------------------------
      where( clouds(:) /= 0. )
         del_lwp(:) = gi * lwc(:) * delp(:) * 1.e3 / clouds(:)
      elsewhere
         del_lwp(:) = 0.
      endwhere
!---------------------------------------------------------
!            ... Form tau for each model layer
!---------------------------------------------------------
      where( clouds(:) /= 0. )
         del_tau(:) = del_lwp(:) *.155 * clouds(:)**1.5
      elsewhere
         del_tau(:) = 0.
      end where
!---------------------------------------------------------
!            ... Form integrated tau from top down
!---------------------------------------------------------
      above_tau(1) = 0.
      do k = 1,plevm
         above_tau(k+1) = del_tau(k) + above_tau(k)
      end do
!---------------------------------------------------------
!            ... Form integrated tau from bottom up
!---------------------------------------------------------
      below_tau(plev) = 0.
      do k = plevm,1,-1
         below_tau(k) = del_tau(k+1) + below_tau(k+1)
      end do
!---------------------------------------------------------
!        ... Form vertically averaged cloud cover above and below
!---------------------------------------------------------
      above_cld(1) = 0.
      do k = 1,plevm
         above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
      end do
      do k = 2,plev
         if( above_tau(k) /= 0. ) then
            above_cld(k) = above_cld(k) / above_tau(k)
         else
            above_cld(k) = above_cld(k-1)
         end if
      end do
      below_cld(plev) = 0.
      do k = plevm,1,-1
         below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
      end do
      do k = plevm,1,-1
         if( below_tau(k) /= 0. ) then
            below_cld(k) = below_cld(k) / below_tau(k)
         else
            below_cld(k) = below_cld(k+1)
         end if
      end do
!---------------------------------------------------------
!        ... Modify above_tau and below_tau via jfm
!---------------------------------------------------------
      where( above_cld(2:plev) /= 0. )
         above_tau(2:plev) = above_tau(2:plev) / above_cld(2:plev)
      end where
      where( below_cld(:plevm) /= 0. )
         below_tau(:plevm) = below_tau(:plevm) / below_cld(:plevm)
      end where
      where( above_tau(2:plev) < 5. )
            above_cld(2:plev) = 0.
      end where
      where( below_tau(:plevm) < 5. )
         below_cld(:plevm) = 0.
      end where
!---------------------------------------------------------
!        ... Form transmission factors
!---------------------------------------------------------
      above_tra(:) = 11.905 / (9.524 + above_tau(:))
      below_tra(:) = 11.905 / (9.524 + below_tau(:))
!---------------------------------------------------------
!        ... Form effective albedo
!---------------------------------------------------------
      where( below_cld(:) /= 0. )
         eff_alb(:) = srf_alb + below_cld(:) * (1. - below_tra(:)) &
                                             * (1. - srf_alb)
      elsewhere
         eff_alb(:) = srf_alb
      end where
      coschi = max( coszen,.5 )
      where( del_lwp(:)*.155 < 5. )
         fac1(:) = 0.
      elsewhere
         fac1(:) = 1.4 * coschi - 1.
      end where
      fac2(:)     = MIN( 0.,1.6*coschi*above_tra(:) - 1. )
      cld_mult(:) = 1. + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
      cld_mult(:) = MAX( .05,cld_mult(:) )

      end subroutine cloud_mod

      subroutine T_INT( p, t, t500, t200 )
!----------------------------------------------------------------
!        ... Interpolate for temperature on 500 and 200mb surfaces
!----------------------------------------------------------------


      implicit none

!----------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------
      real, intent(in)  ::  p(:)              ! pressure in pascals
      real, intent(in)  ::  t(:)              ! temperature on grid
      real, intent(out) ::  t500, t200        ! temp at 500 and 200mb

!----------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------
      integer :: k, k1
      real    :: delp
      integer :: plev, plevp, plevm
      
      plev = SIZE(p)
      plevp = plev+1
      plevm = plev-1

      if( p(plev) < 500.e2 ) then
         t500 = t(plev)
         k1 = plevp
      else
         do k = plevm,1,-1
            if( p(k) < 500.e2 ) then
               k1 = k
               exit
            end if
         end do
         delp = LOG( 500.e2/p(k) ) / LOG( p(k+1)/p(k) )
         t500 = t(k) + delp * (t(k+1) - t(k))
      end if
      do k = k1-1,1,-1
         if( p(k) < 200.e2 ) then
            exit
         end if
      end do
      delp = LOG( 200.e2/p(k) ) / LOG( p(k+1)/p(k) )
      t200 = t(k) + delp * (t(k+1) - t(k))

      end subroutine T_INT

      subroutine PHOTO_INTERP( zin, cin, vin, albin, t500in, &
                               t200in, solar_phase, ajout )
!----------------------------------------------------------------------
!           ... Loglinear interpolation for the photodissociation rates
!            Note: this subroutine computes photorates for a vertical
!                  column at a given longitude and latitude
!           This routine uses a six parameter table via a Taylor
!           series expansion. 
!           Stacy Walters, Sep 30, 1996.  Changed code to strictly limit
!           the 200mb and 500mb temperature interpolation to the table
!           endpoints; i.e. no extrapolation beyond the table is allowed.
!----------------------------------------------------------------------


      implicit none

!----------------------------------------------------------------------
!        ... Dummy arguments
!----------------------------------------------------------------------
      real, intent(in)  ::   zin(:), &              ! geo height of midpoints
                             cin, &                 ! cosine solar zenith angle
                             vin(:), &              ! o3 column density
                             albin(:), &            ! surface albedo
                             t500in, &              ! temp on 500mb surface
                             t200in, &              ! temp on 200mb surface
                             solar_phase            ! phase of solar cycle (1=max, 0=min)
      real, intent(out) ::   ajout(:,:)             ! photodissociation rates

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer  ::  plev
      integer  ::  iz, is, iv, ial, nn, it500, it200
      integer  ::  izp1, isp1, ivp1, ialp1, it500p1, it200p1
      integer  ::  i, k
      integer  ::  izl
      integer, dimension(SIZE(zin)) :: altind, ratind, albind
      real     ::  wght0
      real     ::  v3std
      real     ::  dels(6)
      real, dimension(SIZE(zin)) :: v3rat
      real     :: ajout_tmp
!RSH ADD HERE:
      real    ::  ajl2(2,2,2,2,2,2)
      
      plev = SIZE(zin)

!----------------------------------------------------------------------
!        ... Find the zenith angle index ( same for all levels )
!----------------------------------------------------------------------
      do is = 1,zangdim
         if( vcos(is) > cin ) then
            exit
         end if
      end do
      is       = MAX( MIN( is,zangdim ) - 1,1 )
      isp1     = is + 1
      dels(2)  = MIN( 1.,MAX( 0.,(cin - vcos(is)) * delang(is) ) )
      if (dels(2) > 0.5) then
         dels(2) = 1. - dels(2)
         isp1 = is
         is   = isp1 + 1
      end if

!----------------------------------------------------------------------
!        ... Find the 500mb temp index ( same for all levels )
!----------------------------------------------------------------------
      do it500 = 1,t500dim
         if( t500(it500) > t500in ) then
            exit
         end if
      end do
      it500    = MAX( MIN( it500,t500dim ) - 1,1 )
      it500p1  = it500 + 1
      dels(5)  = MIN( 1.,MAX( 0.,(t500in - t500(it500)) * delt500(it500) ) )
      if (dels(5) > 0.5) then
         dels(5) = 1. - dels(5)
         it500p1  = it500
         it500    = it500p1 + 1
      end if

!----------------------------------------------------------------------
!        ... Find the 200mb temp index ( same for all levels )
!----------------------------------------------------------------------
      it200 = 1
      it200p1 = 2
      dels(6)  = MIN( 1.,MAX( 0.,(t200in - t200(it200)) * delt200(it200) )) 
      if (dels(6) > 0.5) then
         dels(6) = 1. - dels(6)
         it200    = 2
         it200p1  = 1
      end if

      izl = 1
      do k = plev,1,-1
!----------------------------------------------------------------------
!        ... Find albedo indicies
!----------------------------------------------------------------------
         do ial = 1,albdim
            if( albev(ial) > albin(k) ) then
               exit
            end if
         end do
         albind(k) = MAX( MIN( ial,albdim ) - 1,1 )
!----------------------------------------------------------------------
!        ... Find level indicies
!----------------------------------------------------------------------
         do iz = izl,altdim
            if( zz(iz) > zin(k) ) then
               izl = iz
               exit
            end if
         end do
         altind(k) = MAX( MIN( iz,altdim ) - 1,1 )
!----------------------------------------------------------------------
!        ... Find "o3 ratio" indicies
!----------------------------------------------------------------------
         i        = MAX( MIN( 79,INT( zin(k) ) ),0 )
         v3std    = vo3(i) + (zin(k) - REAL(i)) * delvo3(i)
         v3rat(k) = vin(k) / v3std
         do iv = 1,o3ratdim
            if( xv3(iv) > v3rat(k) ) then
               exit
            end if
         end do
         ratind(k) = MAX( MIN( iv,o3ratdim ) - 1,1 )
      end do
Vert_loop : &
      do k = 1,plev
!----------------------------------------------------------------------
!        ... Interval deltas and primary weight
!----------------------------------------------------------------------
         iz    = altind(k)
         izp1  = iz + 1
         dels(1)  = MIN( 1.,MAX( 0.,(zin(k) - zz(iz)) * delz(iz) ) )
         if (dels(1) > 0.5) then
            dels(1) = 1. - dels(1)
            izp1 = iz
            iz   = izp1 + 1
         end if
         iv    = ratind(k)
         ivp1  = iv + 1
         dels(3)  = MIN( 1.,MAX( 0.,(v3rat(k) - xv3(iv)) * delv(iv) ) )
         if (dels(3) > 0.5) then
            dels(3) = 1. - dels(3)
            ivp1 = iv
            iv   = ivp1 + 1
         end if
         ial   = albind(k)
         ialp1 = ial + 1
         dels(4)  = MIN( 1.,MAX( 0.,(albin(k) - albev(ial)) * delalb(ial) ) )
         if (dels(4) > 0.5) then
            dels(4) = 1. - dels(4)
            ialp1 = ial
            ial   = ialp1 + 1
         end if
         wght0    = 1. - SUM( dels )
Rate_loop : &
         do nn = 1,jdim
            ajout(nn,k) = wght0     * ajl(nn,iz,is,iv,ial,it500,it200) &
                          + dels(1) * ajl(nn,izp1,is,iv,ial,it500,it200) &
                          + dels(2) * ajl(nn,iz,isp1,iv,ial,it500,it200) &
                          + dels(3) * ajl(nn,iz,is,ivp1,ial,it500,it200) &
                          + dels(4) * ajl(nn,iz,is,iv,ialp1,it500,it200) &
                          + dels(5) * ajl(nn,iz,is,iv,ial,it500p1,it200) &
                          + dels(6) * ajl(nn,iz,is,iv,ial,it500,it200p1)
            ajl2(:,:,:,:,:,:) = ajl(nn,(/iz,izp1/),(/is,isp1/),(/iv,ivp1/),(/ial,ialp1/), &
                                    (/it500,it500p1/),(/it200,it200p1/))
            ajout(nn,k) = MAX( MIN( ajout(nn,k), MAXVAL(ajl2) ), MINVAL(ajl2) )
            ajout(nn,k) = EXP( ajout(nn,k) )
         end do Rate_loop
         if (use_solar_cycle .and. solar_phase /= 1.) then
            do nn = 1,jdim
               ajout_tmp = wght0     * ajl_solarmin(nn,iz,is,iv,ial,it500,it200) &
                           + dels(1) * ajl_solarmin(nn,izp1,is,iv,ial,it500,it200) &
                           + dels(2) * ajl_solarmin(nn,iz,isp1,iv,ial,it500,it200) &
                           + dels(3) * ajl_solarmin(nn,iz,is,ivp1,ial,it500,it200) &
                           + dels(4) * ajl_solarmin(nn,iz,is,iv,ialp1,it500,it200) &
                           + dels(5) * ajl_solarmin(nn,iz,is,iv,ial,it500p1,it200) &
                           + dels(6) * ajl_solarmin(nn,iz,is,iv,ial,it500,it200p1)
               ajl2(:,:,:,:,:,:) = ajl_solarmin(nn,(/iz,izp1/),(/is,isp1/),(/iv,ivp1/),(/ial,ialp1/), &
                                                (/it500,it500p1/),(/it200,it200p1/))
               ajout_tmp = MAX( MIN( ajout_tmp, MAXVAL(ajl2) ), MINVAL(ajl2) )
               ajout_tmp = EXP( ajout_tmp )
               ajout(nn,k) = ajout(nn,k) + (ajout_tmp-ajout(nn,k)) * (1.-solar_phase)
               ajout(nn,k) = MAX(ajout(nn,k),0.)
            end do
         end if
      end do Vert_loop

      end subroutine PHOTO_INTERP

      subroutine set_ub_col( col_delta, vmr, invariants, pdel, ptop, plonl )
!---------------------------------------------------------------
!        ... Set the column densities at the upper boundary
!---------------------------------------------------------------

      use CHEM_MODS_MOD, only : ncol_abs

      implicit none

!---------------------------------------------------------------
!        ... Dummy args
!---------------------------------------------------------------
      integer, intent(in) ::  plonl
      real, intent(out)   ::  col_delta(:,0:,:)  ! /cm**2
      real, intent(in)    ::  vmr(:,:,:), &               ! xported species vmr
                              invariants(:,:,:), &        ! invariant species
                              pdel(:,:), &                ! pressure thickness of model layers (Pa)
                              ptop(:)                     ! model top pressure (Pa)

!---------------------------------------------------------------
!        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!              The factor 10. is to convert pdel
!              from pascals to dyne/cm**2.
!---------------------------------------------------------------
      real, parameter :: pa_to_dyncm2 = 10. !unit = (dyn/cm2)/pa
      real, parameter :: mw_air = 28.9644   !g/mole
      real, parameter :: grav = 981.         !cm/s2
      real, parameter :: navo = 6.023e23   ! molec/mole
      real, parameter :: xfactor = pa_to_dyncm2 * navo /(grav * mw_air)
      integer :: k, spc_ndx
      integer :: plev
      
      plev = SIZE(invariants,2)

!---------------------------------------------------------------
!        ... Assign column density at the upper boundary
!            The first column is O3 and the second is O2.
!            Add O3 column above top of model.
!---------------------------------------------------------------
      spc_ndx = ox_ndx
      if( spc_ndx < 1 ) then
         spc_ndx = o3_ndx
      end if
      if( spc_ndx > 0 ) then
         col_delta(:,0,1) = 2.687e16*o3_column_top
         do k = 1,plev
            col_delta(:,k,1) = xfactor * pdel(:,k) * vmr(:,k,spc_ndx) ! O3
         end do
      end if
!     col_delta(:,0,2) = 2.8e22
      col_delta(:,0,2) = xfactor * ptop(:) * invariants(:,plev,3)/invariants(:,plev,1)
      do k = 1,plev
         col_delta(:,k,2) = xfactor * pdel(:,k) * invariants(:,k,3)/invariants(:,k,1) ! O2
      end do

      end subroutine set_ub_col

      subroutine setcol( col_delta, col_dens, pdel, plonl )
!---------------------------------------------------------------
!             ... Set the column densities
!---------------------------------------------------------------

      use CHEM_MODS_MOD, only : ncol_abs

      implicit none

!---------------------------------------------------------------
!             ... Dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: plonl
!     real, intent(in)  ::   vmr(plonl,plev,pcnstm1)           ! xported species vmr
      real, intent(in)  ::   pdel(:,:)                         ! delta about midpoints
      real, intent(in)  ::   col_delta(:,0:,:)                 ! layer column densities (molecules/cm^2)
      real, intent(out) ::   col_dens(:,:,:)                   ! column densities ( /cm**2 )

!---------------------------------------------------------------
!        The local variables
!---------------------------------------------------------------
      integer  ::   k, km1      ! alt indicies
      integer  ::   spc_ndx
      integer  ::   plev
      
!---------------------------------------------------------------
!        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!              The factor 10. is to convert pdel
!              from pascals to dyne/cm**2.
!---------------------------------------------------------------
!     real, parameter :: xfactor = 2.8704e21/(9.80616*1.38044)

      plev = SIZE(pdel,2)

!---------------------------------------------------------------
!           ... Compute column densities down to the
!           current eta index in the calling routine.
!           The first column is O3 and the second is O2.
!---------------------------------------------------------------
      spc_ndx = ox_ndx
      if( spc_ndx < 1 ) then
         spc_ndx = o3_ndx
      end if
      if( spc_ndx > 0 ) then
         col_dens(:,1,1) = col_delta(:,0,1) + .5 * col_delta(:,1,1)
         do k = 2,plev
            km1 = k - 1
            col_dens(:,k,1) = col_dens(:,km1,1) + .5 * (col_delta(:,km1,1) + col_delta(:,k,1))
         end do
      end if
      col_dens(:,1,2) = col_delta(:,0,2) + .5 * col_delta(:,1,2)
      do k = 2,plev
         km1 = k - 1
         col_dens(:,k,2) = col_dens(:,km1,2) + .5 * (col_delta(:,km1,2) + col_delta(:,k,2))
      end do

      end subroutine SETCOL


      real function SUNDIS( Time )
!-----------------------------------------------------------------------------
!=  PURPOSE:                                                                 =*
!=  Calculate Earth-Sun distance variation for a given date.  Based on       =*
!=  Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier      =*
!=  series representation of the position of the sun, Search, 2:172          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  IDATE  - INTEGER, specification of the date, from YYMMDD              (I)=*
!=  ESRM2  - REAL, variation of the Earth-sun distance                    (O)=*
!=           ESRM2 = (average e/s dist)^2 / (e/s dist on day IDATE)^2        =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  01/95  Changed computation of trig function values                       =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------------
      type(time_type), intent(in) :: Time             ! time

!-----------------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------------
      integer :: iyear, imonth, iday, ihour, iminute, isecond
      integer :: mday, month, jday
      integer, save :: imn(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      real    :: dayn, thet0
      real    :: sinth, costh, sin2th, cos2th
      character(len=128) :: msg
      
!-----------------------------------------------------------------------------
!         ... Parse date to find day number (Julian day)
!-----------------------------------------------------------------------------
      call get_date( Time, iyear, imonth, iday, ihour, iminute, isecond )
      if( imonth > 12 ) then
         write(msg,*) 'Month in date exceeds 12, month = ',imonth
         call endrun(msg)
      end if

      if( MOD(iyear,4) == 0 ) then
         imn(2) = 29
      else
         imn(2) = 28
      end if

      if( iday > imn(imonth) ) then
         write(msg,*) 'Day in date exceeds days in month, day = ',iday,', month = ',imonth
         call endrun(msg)
      end if

      mday = 0
      do month = 1,imonth-1
         mday = mday + imn(month)                     
      end do
      jday = mday + iday
      dayn = REAL(jday - 1) + .5

!-----------------------------------------------------------------------------
!         ... Define angular day number and compute esrm2:
!-----------------------------------------------------------------------------
      thet0 = 2.*PI*dayn/365.

!-----------------------------------------------------------------------------
!         ... Calculate SIN(2*thet0), COS(2*thet0) 
!-----------------------------------------------------------------------------
      sinth   = SIN( thet0 )
      costh   = COS( thet0 )
      sin2th  = 2.*sinth*costh
      cos2th  = costh*costh - sinth*sinth
      SUNDIS  = 1.000110 + .034221*costh  +  .001280*sinth + .000719*cos2th +  .000077*sin2th

      end function SUNDIS


      subroutine endrun(msg)

      character(len=128), intent(in) :: msg
      call mpp_error(FATAL, msg)
      
      end subroutine endrun        

      end module MO_PHOTO_MOD
