      module MO_FPHOTO_MOD
!----------------------------------------------------------------------
!        ... Read FAST-JX Photolysis rates and pass them to chemdr
!            June. 1st, 2010
!            Junfeng.Liu
!----------------------------------------------------------------------
      use mpp_mod,              only : mpp_error
      use fms_mod,              only : write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pe,             &
                                       mpp_clock_begin, mpp_clock_end, &
                                       mpp_clock_id, CLOCK_MODULE, &
                                       error_mesg, FATAL
      use tracer_manager_mod,   only : get_tracer_index,   &
                                       get_number_tracers, &
                                       get_tracer_names,   &
                                       get_tracer_indices
      use field_manager_mod,    only : MODEL_ATMOS       
      use MO_FASTJX_MOD,        only : JVN_, fastjx_init, fastjx_end, fastjx_photo
      use sat_vapor_pres_mod, only : compute_qs      
         
      implicit none
      private
      public :: fprate_init, fphoto

      character(len=128)            :: version     = '$Id: mo_fphoto.F90,v 19.0 2012/01/06 20:33:54 fms Exp $'
      character(len=128)            :: tagname     = '$Name: tikal $'

      integer ::  fastjx_clock
      integer, parameter :: jdim     = JVN_     ! number of fastjx species 62
      integer ::  indexer(jdim)                 ! index corresponding to AM3 J order 
      
      integer ::  jno_ndx, jpooh_ndx, jc2h5ooh_ndx, jc3h7ooh_ndx, jrooh_ndx, &
                  jch3co3h_ndx, jmpan_ndx, jmacr_a_ndx, jmacr_b_ndx, jonitr_ndx, &
                  jxooh_ndx, jisopooh_ndx, jglyald_ndx, jhyac_ndx, jch3ooh_ndx, &
                  jh2o2_ndx, jpan_ndx, jch3cho_ndx, &
                  jn2o5_ndx, jo3p_ndx, jno2_ndx, jno3_ndx, &
                  jclono2_ndx, jhocl_ndx, jcl2o2_ndx, jbrono2_ndx, jhobr_ndx, &
                  jbrcl_ndx, jbro_ndx, jcl2_ndx, jh2o_ndx, jn2o_ndx, jhno3_ndx, &
!jul++
                  jmek_ndx, jbigald_ndx, jglyoxal_ndx, jalkooh_ndx, jmekooh_ndx, &
                  jtolooh_ndx, jterpooh_ndx, jacet_ndx, jmgly_ndx
      integer ::  jh2o2a_ndx, jhno2_ndx, jo3a_ndx , jo1d_ndx 
      integer ::  so4_ndx, bc1_ndx, bc2_ndx, oc1_ndx, oc2_ndx, soa_ndx, &
                  ssa_ndx(5), dust_ndx(5)   
!jul--
      integer ::  ox_ndx, o3_ndx, nqa, nqi, nql, nqq
 
      integer, parameter :: & 
         TAB_NDX_JO2            = 1, &
         TAB_NDX_JO3            = 2, & 
         TAB_NDX_JO1D           = 3, & 
         TAB_NDX_JNO            = 4, & 
         TAB_NDX_JH2COa         = 5, & 
         TAB_NDX_JH2COb         = 6, & 
         TAB_NDX_JH2O2          = 7, & 
         TAB_NDX_JCH3OOH        = 8, & 
         TAB_NDX_JNO2           = 9, & 
         TAB_NDX_JNO3           = 10, & 
         TAB_NDX_JN2O5          = 11, & 
         TAB_NDX_JHNO2          = 12, & 
         TAB_NDX_JHNO3          = 13, & 
         TAB_NDX_JHNO4          = 14, & 
         TAB_NDX_JClNO3a        = 15, & 
         TAB_NDX_JClNO3b        = 16, & 
         TAB_NDX_JCl2           = 17, & 
         TAB_NDX_JHOCl          = 18, & 
         TAB_NDX_JOClO          = 19, & 
         TAB_NDX_JCl2O2         = 20, & 
         TAB_NDX_JClO           = 21, & 
         TAB_NDX_JBrO           = 22, & 
         TAB_NDX_JBrNO3         = 23, & 
         TAB_NDX_JHOBr          = 24, & 
         TAB_NDX_JBrCl          = 25, & 
         TAB_NDX_JN2O           = 26, & 
         TAB_NDX_JCFCl3         = 27, & 
         TAB_NDX_JCF2Cl2        = 28, & 
         TAB_NDX_JF113          = 29, & 
         TAB_NDX_JF114          = 30, & 
         TAB_NDX_JF115          = 31, & 
         TAB_NDX_JCCl4          = 32, & 
         TAB_NDX_JCH3Cl         = 33, & 
         TAB_NDX_JMeCCl3        = 34, & 
         TAB_NDX_JCH2Cl2        = 35, & 
         TAB_NDX_JCHF2Cl        = 36, & 
         TAB_NDX_JF123          = 37, & 
         TAB_NDX_JF141b         = 38, & 
         TAB_NDX_JF142b         = 39, & 
         TAB_NDX_JCH3Br         = 40, & 
         TAB_NDX_JH1211         = 41, & 
         TAB_NDX_JH1301         = 42, & 
         TAB_NDX_JH2402         = 43, & 
         TAB_NDX_JCH2Br2        = 44, & 
         TAB_NDX_JCHBr3         = 45, & 
         TAB_NDX_JCH3I          = 46, & 
         TAB_NDX_JCF3I          = 47, & 
         TAB_NDX_JOCS           = 48, & 
         TAB_NDX_JPAN           = 49, & 
         TAB_NDX_JCH3NO3        = 50, &  
         TAB_NDX_JActAld        = 51, &  
         TAB_NDX_JMeVK          = 52, &  
         TAB_NDX_JMeAcr         = 53, &  
         TAB_NDX_JGlyAld        = 54, &  
         TAB_NDX_JMEKeto        = 55, &  
         TAB_NDX_JEAld          = 56, &  
         TAB_NDX_JMGlyxl        = 57, &  
         TAB_NDX_JGlyxla        = 58, &  
         TAB_NDX_JGlyxlb        = 59, &  
         TAB_NDX_JAcet_a        = 60, &  
         TAB_NDX_JAcet_b        = 61, &  
         TAB_NDX_JHYAC          = 62  
!      logical :: use_tdep_jvals, use_solar_cycle
         real                   :: o3_column_top  !
         logical                :: module_is_initialized = .false.
      CONTAINS
     
      
!--------------------------------------------------      
! <SUBROUTINE NAME="fprate_init">
!   <OVERVIEW>
!     Initialize FAST-JX photolysis rate calculation, must be called after chemini
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the calculation of FASTJX photolysis rates
!     from the TUV lookup table
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fprate_init( o3_column_top_in)
!   </TEMPLATE>
!   <IN NAME="o3_column_top_in" TYPE="real">
!     Ozone column above model top (DU)
!   </IN>

  subroutine fprate_init(o3_column_top_in)
!----------------------------------------------------------------------
!     ... Initialize FAST-JX module 
!----------------------------------------------------------------------        
      use mo_chem_utls_mod,  only : get_spc_ndx, get_rxt_ndx   

      implicit none
!----------------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------------
      real,             intent(in) :: o3_column_top_in !(DU)

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer    :: it500, it200, izen, ialb, idob
      integer    :: ios
      integer    :: unit
      character(len=128) :: msg     
      
      if (module_is_initialized) return

!-------------------------------------------------------------------------
!     write version number and tagname to stdlog.
!-------------------------------------------------------------------------
      call write_version_number(version, tagname)

      indexer(TAB_NDX_JO2)      = get_rxt_ndx( 'jo2' )
      indexer(TAB_NDX_JO3)      = get_rxt_ndx( 'jo3p' ) 
      indexer(TAB_NDX_JO1D)     = get_rxt_ndx( 'jo1d' ) 
      indexer(TAB_NDX_JNO)      = get_rxt_ndx( 'jno' )       !!
      indexer(TAB_NDX_JH2COa)   = get_rxt_ndx( 'jch2o_a' ) 
      indexer(TAB_NDX_JH2COb)   = get_rxt_ndx( 'jch2o_b' ) 
      indexer(TAB_NDX_JH2O2)    = get_rxt_ndx( 'jh2o2' ) 
      indexer(TAB_NDX_JCH3OOH)  = get_rxt_ndx( 'jch3ooh' )
      indexer(TAB_NDX_JNO2)     = get_rxt_ndx( 'jno2' ) 
      indexer(TAB_NDX_JNO3)     = get_rxt_ndx( 'jno3' ) 
      indexer(TAB_NDX_JN2O5)    = get_rxt_ndx( 'jn2o5' )       
      indexer(TAB_NDX_JHNO2)    = get_rxt_ndx( 'jhno2' )    
      indexer(TAB_NDX_JHNO3)    = get_rxt_ndx( 'jhno3' ) 
      indexer(TAB_NDX_JHNO4)    = get_rxt_ndx( 'jho2no2' ) 
      
      indexer(TAB_NDX_JClNO3a)  = get_rxt_ndx( 'jclono2' )   !!
      indexer(TAB_NDX_JClNO3b)  = 0       !get_rxt_ndx( 'j' )  
      indexer(TAB_NDX_JCl2)     = get_rxt_ndx( 'jcl2' ) 
      indexer(TAB_NDX_JHOCl)    = get_rxt_ndx( 'jhocl' ) 
      indexer(TAB_NDX_JOClO)    = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JCl2O2)   = get_rxt_ndx( 'jcl2o2' ) 
      indexer(TAB_NDX_JClO)     = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JBrO)     = get_rxt_ndx( 'jbro' )      !!
      indexer(TAB_NDX_JBrNO3)   = get_rxt_ndx( 'jbrono2' ) 
      indexer(TAB_NDX_JHOBr)    = get_rxt_ndx( 'jhobr' ) 
      indexer(TAB_NDX_JBrCl)    = get_rxt_ndx( 'jbrcl' ) 
      
      indexer(TAB_NDX_JN2O)     = get_rxt_ndx( 'jn2o' ) 
      
      indexer(TAB_NDX_JCFCl3)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCF2Cl2)  = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JF113)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JF114)    = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JF115)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCCl4)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCH3Cl)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JMeCCl3)  = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JCH2Cl2)  = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCHF2Cl)  = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JF123)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JF141b)   = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JF142b)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCH3Br)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JH1211)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JH1301)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JH2402)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCH2Br2)  = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JCHBr3)   = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JCH3I)    = 0       !get_rxt_ndx( 'j' )
      indexer(TAB_NDX_JCF3I)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JOCS)     = 0       !get_rxt_ndx( 'j' )
     
      indexer(TAB_NDX_JPAN)     = get_rxt_ndx( 'jpan' )
      indexer(TAB_NDX_JCH3NO3)  = 0       !get_rxt_ndx( 'j' )  
      indexer(TAB_NDX_JActAld)  = get_rxt_ndx( 'jch3cho' ) 
      indexer(TAB_NDX_JMeVK)    = get_rxt_ndx( 'jmvk' )  
      indexer(TAB_NDX_JMeAcr)   = get_rxt_ndx( 'jmacr_a' ) 
      indexer(TAB_NDX_JGlyAld)  = get_rxt_ndx( 'jglyald' )  
      indexer(TAB_NDX_JMEKeto)  = get_rxt_ndx( 'jmek' )     !!
      indexer(TAB_NDX_JEAld)    = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JMGlyxl)  = get_rxt_ndx( 'jmgly' )  
      indexer(TAB_NDX_JGlyxla)  = get_rxt_ndx( 'jglyoxal' )  !!??
      indexer(TAB_NDX_JGlyxlb)  = 0       !get_rxt_ndx( 'j' )  
      indexer(TAB_NDX_JAcet_a)  = get_rxt_ndx( 'jacet' )  
      indexer(TAB_NDX_JAcet_b)  = 0       !get_rxt_ndx( 'j' ) 
      indexer(TAB_NDX_JHYAC)    = get_rxt_ndx( 'jhyac' )             
      
    if (mpp_pe() == mpp_root_pe()) &
      write(*,*)'fphoto_init: indexer',indexer     
      
!-----------------------------------------------------------------
!           ... these indices are the additional ones which don't read directly from the table
!-----------------------------------------------------------------
      jh2o_ndx     = get_rxt_ndx( 'jh2o' )
      jpooh_ndx    = get_rxt_ndx( 'jpooh' )
      jch3co3h_ndx = get_rxt_ndx( 'jch3co3h' )
      jmpan_ndx    = get_rxt_ndx( 'jmpan' )
      jmacr_b_ndx  = get_rxt_ndx( 'jmacr_b' )      
      jc2h5ooh_ndx = get_rxt_ndx( 'jc2h5ooh' )
      jc3h7ooh_ndx = get_rxt_ndx( 'jc3h7ooh' )      
      jrooh_ndx    = get_rxt_ndx( 'jrooh' )      
      jxooh_ndx    = get_rxt_ndx( 'jxooh' )
      jonitr_ndx   = get_rxt_ndx( 'jonitr' )
      jisopooh_ndx = get_rxt_ndx( 'jisopooh' )
      jbigald_ndx  = get_rxt_ndx( 'jbigald' )      
      jalkooh_ndx  = get_rxt_ndx( 'jalkooh' )      
      jmekooh_ndx  = get_rxt_ndx( 'jmekooh' )
      jtolooh_ndx  = get_rxt_ndx( 'jtolooh' )
      jterpooh_ndx = get_rxt_ndx( 'jterpooh' )      
!      
      jno_ndx      = get_rxt_ndx( 'jno' )
      jmacr_a_ndx  = get_rxt_ndx( 'jmacr_a' )
      jglyald_ndx  = get_rxt_ndx( 'jglyald' )
      jhyac_ndx    = get_rxt_ndx( 'jhyac' )
      jch3ooh_ndx  = get_rxt_ndx( 'jch3ooh' )
      jh2o2_ndx    = get_rxt_ndx( 'jh2o2' )
      jpan_ndx     = get_rxt_ndx( 'jpan' )
      jch3cho_ndx  = get_rxt_ndx( 'jch3cho' )
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
!
      jn2o_ndx     = get_rxt_ndx( 'jn2o' )
      jhno3_ndx    = get_rxt_ndx( 'jhno3' )
!jul++
      jacet_ndx         = get_rxt_ndx( 'jacet' ) 
      jmgly_ndx         = get_rxt_ndx( 'jmgly' )
      jmek_ndx          = get_rxt_ndx( 'jmek' )
      jglyoxal_ndx      = get_rxt_ndx( 'jglyoxal' )
      jhno2_ndx    = get_rxt_ndx( 'jhno2' )
      jo1d_ndx     = get_rxt_ndx( 'jo1d' )
        
!jul--  

      o3_ndx = get_tracer_index(MODEL_ATMOS,'o3')
      if(  o3_ndx <1 )then
        o3_ndx = get_tracer_index(MODEL_ATMOS,'OX')      
      end if
      so4_ndx     = get_tracer_index(MODEL_ATMOS,'so4')
      bc1_ndx     = get_tracer_index(MODEL_ATMOS,'bcphob')
      bc2_ndx     = get_tracer_index(MODEL_ATMOS,'bcphil')
      oc1_ndx     = get_tracer_index(MODEL_ATMOS,'omphob')
      oc2_ndx     = get_tracer_index(MODEL_ATMOS,'omphil')
      soa_ndx     = get_tracer_index(MODEL_ATMOS,'SOA')
      ssa_ndx(1)    = get_tracer_index(MODEL_ATMOS,'ssalt1')
      ssa_ndx(2)    = get_tracer_index(MODEL_ATMOS,'ssalt2')
      ssa_ndx(3)    = get_tracer_index(MODEL_ATMOS,'ssalt3')
      ssa_ndx(4)    = get_tracer_index(MODEL_ATMOS,'ssalt4')
      ssa_ndx(5)    = get_tracer_index(MODEL_ATMOS,'ssalt5')
      dust_ndx(1)   = get_tracer_index(MODEL_ATMOS,'dust1')
      dust_ndx(2)   = get_tracer_index(MODEL_ATMOS,'dust2') 
      dust_ndx(3)   = get_tracer_index(MODEL_ATMOS,'dust3')
      dust_ndx(4)   = get_tracer_index(MODEL_ATMOS,'dust4')
      dust_ndx(5)   = get_tracer_index(MODEL_ATMOS,'dust5')
       
      nqa = get_tracer_index(MODEL_ATMOS,'cld_amt')
      nqi = get_tracer_index(MODEL_ATMOS,'ice_wat')
      nql = get_tracer_index(MODEL_ATMOS,'liq_wat')
      nqq = get_tracer_index(MODEL_ATMOS,'sphum')
      
      if(  o3_ndx <1 )   call error_mesg ('ATMOS: fphoto','Failed to find O3_ndx', FATAL)
      if(  so4_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find so4_ndx', FATAL)
      if(  bc1_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find bc1_ndx', FATAL)
      if(  bc2_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find bc2_ndx', FATAL)
      if(  oc1_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find oc1_ndx', FATAL)
      if(  oc2_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find oc2_ndx', FATAL)
      if(  soa_ndx <1 )  call error_mesg ('ATMOS: fphoto','Failed to find soa_ndx', FATAL)
      if(  ssa_ndx(1) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ssa1_ndx', FATAL)
      if(  ssa_ndx(2) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ssa2_ndx', FATAL)
      if(  ssa_ndx(3) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ssa3_ndx', FATAL)
      if(  ssa_ndx(4) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ssa4_ndx', FATAL)
      if(  ssa_ndx(5) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ssa5_ndx', FATAL)
      if(  dust_ndx(1) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find dust1_ndx', FATAL)
      if(  dust_ndx(2) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find dust2_ndx', FATAL)
      if(  dust_ndx(3) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find dust3_ndx', FATAL)
      if(  dust_ndx(4) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find dust4_ndx', FATAL)
      if(  dust_ndx(5) <1 )  call error_mesg ('ATMOS: fphoto','Failed to find dust5_ndx', FATAL)
      
      
      if(  nqa <1 )  call error_mesg ('ATMOS: fphoto','Failed to find cld_amt_ndx', FATAL)
      if(  nqi <1 )  call error_mesg ('ATMOS: fphoto','Failed to find ice_wat_ndx', FATAL)
      if(  nql <1 )  call error_mesg ('ATMOS: fphoto','Failed to find liq_wat_ndx', FATAL)
!      if(  nqq <1 )  call error_mesg ('ATMOS: fphoto','Failed to find c_ndx', FATAL)
        
      
      
      if (mpp_pe() == mpp_root_pe() ) write(*,*) 'fphoto:nqa,nqi,nql,nqq= ',nqa,nqi,nql,nqq
!      use_tdep_jvals   = use_tdep_jvals_in
      o3_column_top    = o3_column_top_in * 2.6867E+20/1.e4 !DU -> molec/cm2
!      jno_scale_factor = jno_scale_factor_in

      fastjx_clock = mpp_clock_id( 'Tropchem: Fast_JX', &
           grain=CLOCK_MODULE )

      call fastjx_init
      module_is_initialized = .true.

      end subroutine fprate_init
! </SUBROUTINE>


! <SUBROUTINE NAME="FPHOTO">
!   <OVERVIEW>
!     Calculate FAST-JX photolysis rates
!   </OVERVIEW>
!   <DESCRIPTION>
!     Calculate photolysis rates from FAST-JX codes
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
      subroutine FPHOTO( photos, &
                         pmid, &         ! pfull
                         pdel, &
                         temper, &       !tfull
                         zmid, &
                         col_dens, &
                         coszen,  & 
                         srf_alb, &
                         j_ndx,&
                         lwc, &
                         clouds, &
                         esfact, &
                         solar_phase, &
                         plonl, &
                         use_lsc_in_fastjx, &
                         phalf,&
                         zhalf,&
                         pwt , &
                         qfld, &
                         r  )

      use CHEM_MODS_MOD, only : ncol_abs, phtcnt
      use time_manager_mod, only : time_type      

      implicit none

!-----------------------------------------------------------------
!           ... Dummy arguments
!-----------------------------------------------------------------
!      type(time_type), intent(in) :: Time             ! time
      integer, intent(in) :: plonl, j_ndx              ! lat index
      logical, intent(in) :: use_lsc_in_fastjx         ! are lsc clouds
                                                       ! seen in fastjx ?
      real,    intent(in) :: esfact, &                 ! earth sun distance factor
                             solar_phase               ! solar cycle phase (1=max, 0=min)
      real,    intent(in) :: col_dens(:,:,:), &        ! column densities
                             coszen(:), &              ! solar zenith angle
                             srf_alb(:), &             ! surface albedo
                             lwc(:,:), &               ! liquid water content (mass mr)
!                            lwc_water(:,:),&          ! water cloud (Kg/Kg)
!                     lwc_ice(:,:),&                   ! ice cloud (Kg/Kg)
                             clouds(:,:), &            ! cloud fraction
                             pmid(:,:), &              ! midpoint pressure in pascals
                             pdel(:,:), &              ! del pressure about midpoint in pascals
                             zmid(:,:), &              ! midpoint height
                             temper(:,:), &            ! midpoint temperature
                             phalf(:,:),&              ! pressure at boundaries (pa)
                             zhalf(:,:),&              ! boundary hight (m)
                             pwt(:,:) , &              ! air column density (Kg/m2)
                             qfld(:,:), &              ! specific humnidity (Kg/Kg)
                             r(:,:,:)                  ! tracers' concentrtaions             
      real,   intent(out) :: photos(:,:,:)             ! photodissociation rates (s-1)
!-----------------------------------------------------------------
!            ... Local variables
!-----------------------------------------------------------------
      integer, parameter   ::  naero = 16   ! number of aerosol type
      integer, parameter   ::  ncloud = 2   ! number of aerosol type
      
      integer  ::  i,j, k, m                 ! indicies
      integer  ::  plev
!      integer, dimension(size(zmid,2)) :: &
!                   fac1, &                ! work space for J(no) calc
!                   fac2, &                ! work space for J(no) calc
!                   colo3, &               ! vertical o3 column density
!                   zarg, &                ! vertical height array
!                   pline, &               ! vertical pressure array
!                   tline, &               ! vertical temperature array
!                   cld_line, &            ! vertical cloud array
!                   lwc_line, &            ! vertical lwc array
!                   eff_alb, &             ! effective albedo from cloud modifications
!                   cld_mult               ! clould multiplier
      integer, dimension(plonl,size(zmid,2),ncloud) :: clouds_ndx ! cloud type index
      real, dimension(plonl,size(zmid,2),ncloud) :: clouds_lwc, clouds_fraction ! cloud liquid water contains (KG/KG)    
      integer, dimension(plonl,size(zmid,2),naero)    :: aeron ! aerosol type index
      real, dimension(plonl,size(zmid,2),naero) :: aerop        ! aerosol path (g/m2) 
                     

      real, dimension(plonl,size(zmid,2))       :: XO3          ! O3 (VMR)

    
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
                   tmp_jno, &
                   tmp_jglyxlb, &
                   tmp_jo1d   !jul++
   
      real*8    :: prates(size(zmid,2),jdim)        ! photorates matrix
      real, dimension(SIZE(temper,1),SIZE(temper,2)) :: relhum                   ! relative humidity

      call mpp_clock_begin (fastjx_clock)

!-------------check model initialization -------------------------------

      if (.not. module_is_initialized) &
            call error_mesg ('ATMOS: fphoto','fprate_init must be called first.', FATAL)
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
      j = j_ndx
      plev = SIZE(zmid,2)
!-----------------------------------------------------------------
!        ... Zero all photorates
!-----------------------------------------------------------------
      do m = 1,max(1,phtcnt)
         do k = 1,plev
            photos(:,k,m) = 0.
         end do
      end do
      
      XO3 = r(:,:,o3_ndx)
      
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
         tmp_jglyxlb(:,k)     = 0.
         tmp_jo1d(:,k)     = 0. !jul++
      end do

!-----------------------------------------------------------------
!        ... Calc relative humidity
!-----------------------------------------------------------------            
      do k =1, plev
         call rh_calc( pmid(:,k), temper(:,k), qfld(:,k), relhum(:,k) )
      end do
!-----------------------------------------------------------------
!        ... Assign clouds info
!-----------------------------------------------------------------
      
      clouds_lwc(:,:,:) = 0.
      clouds_fraction(:,:,:) = 0.

      if(use_lsc_in_fastjx)then
         clouds_lwc(:,:,1)      = max(0.,r(:,:, nql))  
         clouds_lwc(:,:,2)      = max(0.,r(:,:, nqi))  
         clouds_fraction(:,:,1) = r(:,:, nqa)       
      end if
      
   
      clouds_fraction(:,:,1) = max( 0.,min( 1.,clouds_fraction(:,:,1) ) )
      clouds_fraction(:,:,2) = clouds_fraction(:,:,1)
!  assign cloud index
      clouds_ndx(:,:,1) = 9          !C1, r_eff=12um
      where (temper(:,:) > 233.15)   !>-40C
         clouds_ndx(:,:,2) = 12      !Ice-Hexagonal
      elsewhere
         clouds_ndx(:,:,2) = 13      !Ice-Irregular
      endwhere 
!-----------------------------------------------------------------
!        ... Assign aerosols info
!-----------------------------------------------------------------
      call set_aerosol_mc(r(:,:,:),pwt(:,:),relhum(:,:), aerop(:,:,:),aeron(:,:,:))  ! aerop: g/m2
      do i = 1,plonl
!           if (mpp_pe() == mpp_root_pe()  ) then 
!       write(*,*)'Fphoto: ilon,jlat: ', i,j_ndx
!   endif
   call fastjx_photo(    coszen(i), &
                         esfact, &
                         phalf(i,:),&
                         zhalf(i,:),& 
                         pmid(i,:),&
                         temper(i,:),&
                         XO3(i,:),  &
                         pwt(i,:) , &
                         clouds_lwc(i,:,:),  &
                         clouds_fraction(i,:,:), &
                         clouds_ndx(i,:,:) ,&
                         qfld(i,:), &
                         srf_alb(i), &
                         aerop(i,:,:), &
                         aeron(i,:,:), &
                         o3_column_top, &
                         prates &
                         )

            prates(:,:) = max(0., prates(:,:)) 
            do m = 1,jdim         
                  if( indexer(m) > 0 ) then
                    photos(i,:,indexer(m)) = prates(:,m)
                  else
                     select case( m )
                        case( TAB_NDX_JNO )
                           tmp_jno(i,:) = prates(:,m)
                        case( TAB_NDX_JCH3OOH )
                           tmp_jch3ooh(i,:) = prates(:,m)
                        case( TAB_NDX_JH2O2 )
                           tmp_jh2o2(i,:) = prates(:,m)
                        case( TAB_NDX_JActAld )
                           tmp_jch3cho(i,:) = prates(:,m)
                        case( TAB_NDX_JPAN )
                           tmp_jpan(i,:) = prates(:,m)
                        case( TAB_NDX_JMeAcr )
                           tmp_jmacr_a(i,:) = prates(:,m)
                        case( TAB_NDX_JGlyxlb)
                           tmp_jglyxlb(i,:) = prates(:,m)
                        case( TAB_NDX_JO1D)  !jul++
                           tmp_jo1d(i,:) = prates(:,m)
   
!                        case( TAB_NDX_JN2O_200 )
!                           tmp_jn2o_200(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JN2O_250 )
!                           tmp_jn2o_250(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JN2O_300 )
!                           tmp_jn2o_300(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JN2O5_225 )
!                           tmp_jn2o5_225(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JN2O5_250 )
!                           tmp_jn2o5_250(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JN2O5_300 )
!                           tmp_jn2o5_300(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JHNO3_200 )
!                           tmp_jhno3_200(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JHNO3_250 )
!                           tmp_jhno3_250(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JHNO3_300 )
!                           tmp_jhno3_300(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JCLONO2_200 )
!                           tmp_jclono2_200(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JCLONO2_250 )
!                           tmp_jclono2_250(i,:) = esfact *prates(m,:) * cld_mult(:)
!                        case( TAB_NDX_JCLONO2_300 )
!                           tmp_jclono2_300(i,:) = esfact *prates(m,:) * cld_mult(:)
                     end select
                  end if
           end do
      end do !i
        

!      jh2o_ndx     = get_rxt_ndx( 'jh2o' )      
!      jpooh_ndx    = get_rxt_ndx( 'jpooh' )
!      jch3co3h_ndx = get_rxt_ndx( 'jch3co3h' )
!      jmpan_ndx    = get_rxt_ndx( 'jmpan' )
!      jmacr_b_ndx  = get_rxt_ndx( 'jmacr_b' )      
!      jc2h5ooh_ndx = get_rxt_ndx( 'jc2h5ooh' )
!      jc3h7ooh_ndx = get_rxt_ndx( 'jc3h7ooh' )      
!      jrooh_ndx    = get_rxt_ndx( 'jrooh' )      
!      jxooh_ndx    = get_rxt_ndx( 'jxooh' )
!      jonitr_ndx   = get_rxt_ndx( 'jonitr' )
!      jisopooh_ndx = get_rxt_ndx( 'jisopooh' )
!      jbigald_ndx  = get_rxt_ndx( 'jbigald' )      
!      jalkooh_ndx  = get_rxt_ndx( 'jalkooh' )      
!      jmekooh_ndx  = get_rxt_ndx( 'jmekooh' )
!      jtolooh_ndx  = get_rxt_ndx( 'jtolooh' )
!      jterpooh_ndx = get_rxt_ndx( 'jterpooh' )      
        
      
      
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
      
!jul++
      if( jalkooh_ndx > 0 ) then
         photos(:,:,jalkooh_ndx) = tmp(:,:)
      end if
      
      if( jmekooh_ndx > 0 ) then
         photos(:,:,jmekooh_ndx) = tmp(:,:)
      end if
      
      if( jtolooh_ndx > 0 ) then
         photos(:,:,jtolooh_ndx) = tmp(:,:)
      end if
      
      if( jterpooh_ndx > 0 ) then
         photos(:,:,jterpooh_ndx) = tmp(:,:)
      end if
      
       
!      if( jmek_ndx > 0 ) then
!         if( jacet_ndx > 0 ) then
!            photos(:,:,jmek_ndx) = photos(:,:,jacet_ndx)
!!         end if
!      end if
      if( jbigald_ndx > 0 ) then
         if( jno2_ndx > 0 ) then
            photos(:,:,jbigald_ndx) = 0.2 * photos(:,:,jno2_ndx)
         end if
      end if      
!      if( jglyoxal_ndx > 0 ) then
!         if( jmgly_ndx > 0 ) then
!            photos(:,:,jglyoxal_ndx) = photos(:,:,jmgly_ndx)
!         end if
 !     end if       
!jul--         
      
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
      
!      if( jmacr_a_ndx > 0 ) then
!         photos(:,:,jmacr_a_ndx)  = .5 * photos(:,:,jmacr_a_ndx)
!      end if
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
      
!      if( jglyald_ndx > 0 ) then
!         if( jch3cho_ndx > 0 ) then
!            photos(:,:,jglyald_ndx)  = 3. * photos(:,:,jch3cho_ndx)
!         else
!            photos(:,:,jglyald_ndx)   = 3. *tmp_jch3cho(:,:)
 !        end if
!      end if
      if( jh2o_ndx > 0 ) then
         if( jno_ndx > 0 ) then
            photos(:,:,jh2o_ndx) = 0.1*photos(:,:,jno_ndx)
         else
            photos(:,:,jh2o_ndx) = 0.1*tmp_jno(:,:)
         end if
      end if
      
! be careful here about glyoxala and b in fastjx, here add both      
!      if( jglyoxal_ndx > 0 ) then
!            photos(:,:,jglyoxal_ndx) = photos(:,:,jglyoxal_ndx)+tmp_jglyxlb(:,:)
!      end if
      
      
      
!      if (mpp_pe() == mpp_root_pe()) &
!      write(*,*)'fphoto: photoes',photos(1,1,:)
     
!      if( jn2o_ndx > 0 .and. use_tdep_jvals ) then
!         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
!         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
!         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
!         photos(:,:,jn2o_ndx)    = wgt200(:,:)*tmp_jn2o_200(:,:) + &
!                                   wgt250(:,:)*tmp_jn2o_250(:,:) + &
!                                   wgt300(:,:)*tmp_jn2o_300(:,:)
!      end if
!      if( jn2o5_ndx > 0 .and. use_tdep_jvals ) then
!         wgt225(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/25. ) )
!         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
!         wgt250(:,:)  = 1. - wgt225(:,:) - wgt300(:,:)
!         photos(:,:,jn2o5_ndx)   = wgt225(:,:)*tmp_jn2o5_225(:,:) + &
!                                   wgt250(:,:)*tmp_jn2o5_250(:,:) + &
!                                   wgt300(:,:)*tmp_jn2o5_300(:,:)
!      end if
!      if( jhno3_ndx > 0 .and. use_tdep_jvals ) then
!         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
!         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
!         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
!         photos(:,:,jhno3_ndx)   = wgt200(:,:)*tmp_jhno3_200(:,:) + &
!                                   wgt250(:,:)*tmp_jhno3_250(:,:) + &
!                                   wgt300(:,:)*tmp_jhno3_300(:,:)
!      end if
!      if( jclono2_ndx > 0 .and. use_tdep_jvals ) then
!         wgt200(:,:)  = MIN( 1.,MAX( 0.,(250.-temper(:,:))/50. ) )
!         wgt300(:,:)  = MIN( 1.,MAX( 0.,(temper(:,:)-250.)/50. ) )
!         wgt250(:,:)  = 1. - wgt200(:,:) - wgt300(:,:)
!         photos(:,:,jclono2_ndx) = wgt200(:,:)*tmp_jclono2_200(:,:) + &
!                                   wgt250(:,:)*tmp_jclono2_250(:,:) + &
!                                   wgt300(:,:)*tmp_jclono2_300(:,:)
!      end if


      call mpp_clock_end (fastjx_clock)

      end subroutine FPHOTO
! </SUBROUTINE>



subroutine set_aerosol_mc(r,pwt,rh,aerop,aeron)
!----------------------------------------------------------------
!     set aerosol information for fast-jx calculation
!----------------------------------------------------------------
      implicit none
      real, intent(in)          :: r(:,:,:)     !species
      real, intent(in)          :: pwt(:,:), &  !air column mass (kg/m2)
                                   rh(:,:)      !relative humidity
      real, intent(inout)       :: aerop(:,:,:) !aerosol mass column (g/m2)
      integer, intent(out)      :: aeron(:,:,:) !aerosol category index, see am3_scat.dat

!-----------------------------------------------------------------
!     local parameter variables
!-----------------------------------------------------------------

      integer,  parameter :: nso4bc=702, nbcphil=18, nocphob =1, nocphil=21, nssalt=130, ndust=8
      integer,  parameter :: nso4bc_vol = 27, nso4bc_rh = 26, nssalt_rh = 26
      integer,  parameter :: so4bc_rh1(nso4bc_rh)=(/30,35,40,45,50,55,60,65,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99,100/)
      integer,  parameter :: so4bc_vol1(nso4bc_vol)=(/100,98,96,94,92,90,88,86,84,82,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0/)
      integer,  parameter :: bcphil_rh1(nbcphil)=(/0,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
      integer,  parameter :: ocphil_rh1(nocphil)=(/0,55,60,65,70,75,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
      integer,  parameter :: ssalt_rh1(nssalt_rh)=(/0,47,50,55,60,65,70,72,74,76,78,80,82,84,86,88,90,91,92,93,94,95,96,97,98,99/)
      
      real,  parameter    :: denso4 = 1.74, denbc = 1.0, denocphob = 1.0, denocphil = 1.5 
      real,  parameter    :: denssalt(5) = (/0.078222,0.3580154,0.7310735, 1.125,    2.2774594/)
      real,  parameter    :: dendust(8) = (/2.5,     2.5,     2.5,     2.5,     2.6,    2.6,    2.6,     2.6  /)
      real,  parameter   ::  fractiondust(8) = (/ 0.01053,0.08421,0.25263,0.65263,1.,1.,1.,1. /)
!-----------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------

      real, dimension(size(r,1),size(r,2))      ::      so4, &          ! ammonium sulfate (VMR)
                                                        bc1, &          ! Hydrophobic carbon (MMR)
                                                        bc2, &          ! Hydrophilic carbon (MMR)
!                                                       oc1, &
!                                                       oc2, &
!                                                       soa, &
!                                                       ssa1, &
!                                                       ssa2, &
!                                                       ssa3, &
!                                                       ssa4, &
!                                                       ssa5, &
!                                                       dust11, &               ! dust1 (MMR)
!                                                       dust12, &               ! dust1 (MMR)
!                                                       dust13, &               ! dust1 (MMR)
!                                                       dust14, &               ! dust1 (MMR)
!                                                       dust1, &                ! dust1 (MMR)
!                                                       dust2, &                ! dust2 (MMR) 
!                                                       dust3, &                ! dust3 (MMR)
!                                                       dust4, &                        ! dust4 (MMR)
!                                                       dust5, &                        ! dust4 (MMR)
                                                        v1, v2, relhum

      integer   :: i, j, kd, id,   st1, naero
      integer, dimension(size(r,1),size(r,2))   :: ivol, irh
      
      
      id =size(r,1)
      kd =size(r,2)
      naero = size(aerop,3)

!----------------------------------------------------------------
!     SO4 + BC internal mixing
!----------------------------------------------------------------
      so4(:,:) = r(:,:,so4_ndx)*132./28.97     !VMR => MMR
      bc1(:,:) = r(:,:,bc1_ndx)
      bc2(:,:) = r(:,:,bc2_ndx)
      relhum(:,:) = rh(:,:) *100.
!      if (mpp_pe() == mpp_root_pe() ) then 
!             write(*,*) 'fphoto: so4=',so4(2,:) 
!             write(*,*) 'fphoto: bc1=',bc1(2,:) 
!             write(*,*) 'fphoto: bc2=',bc2(2,:)    
!      end if
      
      
      v1(:,:) = so4(:,:)/denso4
      v2(:,:) = (bc1(:,:)+bc2(:,:))/denbc
      v1(:,:) = v1(:,:)/(v1(:,:)+v2(:,:))*100.
!      if (mpp_pe() == mpp_root_pe() ) then 
!             write(*,*) 'fphoto: v1=',v1(2,:)      
!      end if     
      
      ivol(:,:) = 0
      irh(:,:)  = 0
      call find_indx2(so4bc_vol1(:),v1(:,:),ivol(:,:))
      call find_indx(so4bc_rh1(:), relhum(:,:), irh(:,:))
!      if (mpp_pe() == mpp_root_pe() ) then 
!             write(*,*) 'fphoto: ivol=',ivol(2,:)      
!             write(*,*) 'fphoto: irh=',irh(2,:)        
!      end if        
      aeron(:,:,1) = (ivol(:,:)-1)*nso4bc_rh + irh(:,:)
      aerop(:,:,1) = so4(:,:)+bc1(:,:)+bc2(:,:)
!----------------------------------------------------------------
!     OC 
!----------------------------------------------------------------
!     --phobic
      st1 = nso4bc + nbcphil 
      aeron(:,:,2) = st1 + 1
      aerop(:,:,2) = r(:,:,oc1_ndx)      
!     --philic      
      st1 = st1 + nocphob
      irh(:,:)  = 0
      call find_indx(ocphil_rh1(:), relhum(:,:), irh(:,:))
      aeron(:,:,3) = st1+ irh(:,:)
      aerop(:,:,3) = r(:,:,oc2_ndx) + r(:,:,soa_ndx)
!----------------------------------------------------------------
!     sea salt
!----------------------------------------------------------------
      st1 = st1 + nocphil
      irh(:,:)  = 0
      call find_indx(ssalt_rh1(:), relhum(:,:), irh(:,:))      
      do i= 1, 5
         aeron(:,:,3+i) = st1 + (i-1)*nssalt_rh + irh(:,:) 
         aerop(:,:,3+i) = r(:,:,ssa_ndx(i))          
      end do
!----------------------------------------------------------------
!     dust
!----------------------------------------------------------------
      st1 = st1 + nssalt
      do i= 1, 4
         aeron(:,:,8+i) = st1 + i
         aerop(:,:,8+i) = r(:,:,dust_ndx(1)) *  fractiondust(i)        
      end do      
      do i= 5, 8
         aeron(:,:,8+i) = st1 + i
         aerop(:,:,8+i) = r(:,:,dust_ndx(i-3))        
      end do      
!----------------------------------------------------------------
!     convert mmr to g/m2
!----------------------------------------------------------------    
      do i=1, naero  
         aerop(:,:,i) =aerop(:,:,i) *pwt(:,:)*1000.
      end do
end subroutine set_aerosol_mc




subroutine find_indx(t,tt, it)
  implicit none

  integer, intent(in)   :: t(:)         ! defined categories
  real, intent(in)      :: tt(:,:)      ! model simulated values
  integer, intent(out)  :: it(:,:)      ! index

  
  integer            :: N_ ! number of parameters
  integer               :: i
  real, dimension(size(t,1)-1)          :: dt 
  
  
  N_ = size(t,1)
  
  do i = 1, N_-1
      dt(i) = (t(i+1)-t(i))/2.      
  end do
  
  where (tt(:,:) .le. t(1) + dt(1)  )
      it(:,:) = 1
  endwhere
      
  where (tt(:,:) .gt. t(N_-1) + dt(N_-1)  )
      it(:,:) = N_
  endwhere
  
  do i=1,N_-2
        where(tt(:,:) .gt. t(i) + dt(i) .and. tt(:,:) .le. t(i+1) + dt(i+1)) 
         it(:,:) = i+1
        endwhere      
  end do

end subroutine find_indx



subroutine find_indx2(t,tt, it)
  implicit none

  integer, intent(in)   :: t(:)         ! defined categories
  real, intent(in)      :: tt(:,:)      ! model simulated values
  integer, intent(out)  :: it(:,:)      ! index

  
  integer       :: N_    ! number of parameters
  integer       :: i
  real, dimension(size(t,1)-1)          :: dt 
  
  
  N_ = size(t,1)
  
  do i = 1, N_-1
      dt(i) = (t(i+1)-t(i))/2.      
  end do
  
  where (tt(:,:) .ge. t(1) + dt(1)  )
      it(:,:) = 1
  endwhere
      
  where (tt(:,:) .lt. t(N_-1) + dt(N_-1)  )
      it(:,:) = N_
  endwhere
  
  do i=1,N_-2
        where(tt(:,:) .lt. t(i) + dt(i) .and. tt(:,:) .ge. t(i+1) + dt(i+1)) 
        it(:,:) = i+1
        endwhere      
  end do

end subroutine find_indx2

!      subroutine endrun(msg)

!      character(len=128), intent(in) :: msg
!      call mpp_error(FATAL, msg)
      
!      end subroutine endrun        



subroutine rh_calc(pmid, temp, sh, rh)
              
        implicit none
        
        real, intent(in), dimension(:) :: pmid, temp, sh
        real, intent(out), dimension(:) :: rh
        
        real, dimension(size(temp)) :: esat
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = Rdgas/RVgas = d622
!
!       and where 1- epsilon = d378
!
!       Note that rh does not have its proper value
!       until all of the following code has been executed.  That
!       is, rh is used to store intermediary results
!       in forming the full solution.
!-----------------------------------------------------------------------
        
!-----------------------------------------------------------------------
!calculate water saturated specific humidity
!-----------------------------------------------------------------------
        call compute_qs (temp, pmid, rh, q = sh)
        
!-----------------------------------------------------------------------
!calculate rh
!-----------------------------------------------------------------------
        rh(:)= sh(:) / rh(:)
        
end subroutine rh_calc




end module MO_FPHOTO_MOD






