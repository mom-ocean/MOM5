
      Module co2int_mod

!-----------------------------------------------------------------------
!
!       CO2INT INTERPOLATES CARBON DIOXIDE TRANSMISSION FUNCTIONS
!  FROM THE 109 LEVEL GRID,FOR WHICH THE TRANSMISSION FUNCTIONS
!  HAVE BEEN PRE-CALCULATED, TO THE GRID STRUCTURE SPECIFIED BY THE
!  USER.
!
!        METHOD: 
!
!      CO2INT IS EMPLOYABLE FOR TWO PURPOSES: 1) TO OBTAIN TRANSMIS-
!  SIVITIES BETWEEN ANY 2 OF AN ARRAY OF USER-DEFINED PRESSURES; AND
!  2) TO OBTAIN LAYER-MEAN TRANSMISSIVITIES BETWEEN ANY 2 OF AN ARRAY
!  OF USER-DEFINED PRESSURE LAYERS.TO CLARIFY THESE TWO PURPOSES,SEE
!  THE DIAGRAM AND DISCUSSION BELOW.
!      CO2INT MAY BE USED TO EXECUTE ONLY ONE PURPOSE AT ONE TIME.
!
!     LET P BE AN ARRAY OF USER-DEFINED PRESSURES
!     AND PD BE USER-DEFINED PRESSURE LAYERS.
!
!       - - - - - - - - -   PD(I-1) ---
!                                     !
!       -----------------   P(I)      !  PRESSURE LAYER I  (PLM(I))
!                                     !
!       - - - - - - - - -   PD(I)  ---
!                                     !
!       -----------------   P(I+1)    !  PRESSURE LAYER I+1 (PLM(I+1))
!                                     !
!       - - - - - - - - -   PD(I+1)---
!            ...                          (THE NOTATION USED IS
!            ...                          CONSISTENT WITH THE CODE)
!            ...
!      - - - - - - - - -    PD(J-1)
!
!      -----------------    P(J)
!
!      - - - - - - - - -    PD(J)
!
!      PURPOSE 1:   THE TRANSMISSIVITY BETWEEN SPECIFIC PRESSURES
!      P(I) AND P(J) ,TAU(P(I),P(J))  IS COMPUTED BY THIS PROGRAM.
!      IN THIS MODE,THERE IS NO REFERENCE TO LAYER PRESSURES PD
!      (PD,PLM ARE NOT INPUTTED).
!
!      PURPOSE 2:   THE LAYER-MEAN TRANSMISSIVITY BETWEEN A LAYER-
!      MEAN PRESSURE PLM(J) AND PRESSURE LAYER I IS GIVEN BY
!         TAULM(PLM(I),PLM(J)). IT IS COMPUTED BY THE INTEGRAL
!
!                           PD(I)
!                           ----
!             1             !
!        -------------  *   !   TAU ( P',PLM(J) )  DP'
!        PD(I)-PD(I-1)      !
!                        ----
!                        PD(I-1)
!
!           THE LAYER-MEAN PRESSURE PLM(I) IS SPECIFIED BY THE USER.
!        FOR MANY PURPOSES,PLM WILL BE CHOSEN TO BE THE AVERAGE
!        PRESSURE IN THE LAYER-IE,PLM(I)=0.5*(PD(I-1)+PD(I)).
!           FOR LAYER-MEAN TRANSMISSIVITIES,THE USER THUS INPUTS
!        A PRESSURE ARRAY (PD) DEFINING THE PRESSURE LAYERS AND AN
!        ARRAY (PLM) DEFINING THE LAYER-MEAN PRESSURES.THE CALCULATION
!        DOES NOT DEPEND ON THE P ARRAY USED FOR PURPOSE 1 (P IS NOT
!        INPUTTED).
!
!            THE FOLLOWING PARAGRAPHS DEPICT THE UTILIZATION OF THIS
!       CODE WHEN USED TO COMPUTE TRANSMISSIVITIES BETWEEN SPECIFIC
!       PRESSURES. LATER PARAGRAPHS DESCRIBE ADDITIONAL FEATURES NEEDED
!       FOR LAYER-MEAN TRANSMISSIVITIES.
!
!          FOR A GIVEN CO2 MIXING RATIO AND STANDARD TEMPERATURE
!      PROFILE,A TABLE OF TRANSMISSION FUNCTIONS FOR A FIXED GRID
!     OF ATMOSPHERIC PRESSURES HAS BEEN PRE-CALCULATED.
!      THE STANDARD TEMPERATURE PROFILE IS COMPUTED FROM THE US
!     STANDARD ATMOSPHERE (1977) TABLE.ADDITIONALLY, THE
!     SAME TRANSMISSION FUNCTIONS HAVE BEEN PRE-CALCULATED FOR A
!     TEMPERATURE PROFILE INCREASED AND DECREASED (AT ALL LEVELS)
!     BY 25 DEGREES.
!         THIS PROGRAM READS IN THE PRESPECIFIED TRANSMISSION FUNCTIONS
!     AND A USER-SUPPLIED PRESSURE GRID (P(I)) AND CALCULATES TRANS-
!     MISSION FUNCTIONS ,TAU(P(I),P(J)), FOR ALL P(I)'S AND P(J)'S.
!     A LOGARITHMIC INTERPOLATION SCHEME IS USED.
!         THIS METHOD IS REPEATED FOR THE THREE TEMPERATURE PROFILES
!     GIVEN ABOVE .THEREFORE OUTPUTS FROM THE PROGRAM ARE THREE TABLES
!     OF TRANSMISSION FUNCTIONS FOR THE USER-SUPPLIED PRESSURE GRID.
!     THE EXISTENCE OF THE THREE TABLES PERMITS SUBSEQUENT INTERPO-
!     LATION TO A USER-SUPPLIED TEMPERATURE PROFILE USING THE METHOD
!     DESCRIBED IN THE REFERENCE.SEE LIMITATIONS SECTION IF THE
!     USER DESIRES TO OBTAIN ONLY 1 TABLE OF TRANSMISSIVITIES.
!
!     MODIFICATIONS FOR LAYER-MEAN TRANSMISSIVITIES: 
!          THE PRESSURES INPUTTED ARE THE LAYER-MEAN PRESSURES,PD,
!     AND THE LAYER-MEAN PRESSURES ,PLM. A SERIES OF TRANSMISSIVITIES
!     (TAU(P',PLM(J)) ARE COMPUTED AND THE INTEGRAL GIVEN IN THE
!     DISCUSSION OF PURPOSE 2 IS COMPUTED.FOR PLM(I) NOT EQUAL TO
!     PLM(J) SIMPSON'S RULE IS USED WITH 5 POINTS. IF PLM(I)=PLM(J)
!     (THE "NEARBY LAYER" CASE) A 49-POINT QUADRATURE IS USED FOR
!     GREATER ACCURACY.THE OUTPUT IS IN TAULM(PLM(I),PLM(J)).
!        NOTE: 
!     TAULM IS NOT A SYMMETRICAL MATRIX. FOR THE ARRAY ELEMENT
!     TAULM(PLM(I),PLM(J)),THE INNER(FIRST,MOST RAPIDLY VARYING)
!     DIMENSION IS THE VARYING LAYER-MEAN PRESSURE,PLM(I);THE OUTER
!     (SECOND) DIMENSION IS THE FIXED LAYER-MEAN PRESSURE PLM(J).
!     THUS THE ELEMENT TAULM(2,3) IS THE TRANSMISSION FUNCTION BETWEEN
!     THE FIXED PRESSURE PLM(3) AND THE PRESSURE LAYER HAVING AN AVERAGE
!     PRESSURE OF PLM(2).
!         ALSO NOTE THAT NO QUADRATURE IS PERFORMED OVER THE LAYER
!     BETWEEN THE SMALLEST NONZERO PRESSURE AND ZERO PRESSURE;
!     TAULM IS TAULM(0,PLM(J)) IN THIS CASE,AND TAULM(0,0)=1.
!
!
!             REFERENCE: 
!         S.B.FELS AND M.D.SCHWARZKOPF,"AN EFFICIENT,ACCURATE
!     ALGORITHM FOR CALCULATING CO2 15 UM BAND COOLING RATES",JOURNAL
!     OF GEOPHYSICAL RESEARCH,VOL.86,NO. C2, PP.1205-1232,1981.
!        MODIFICATIONS TO THE ALGORITHM HAVE BEEN MADE BY THE AUTHORS;
!     CONTACT S.B.F.OR M.D.S. FOR FURTHER DETAILS.A NOTE TO J.G.R.
!     IS PLANNED TO DOCUMENT THESE CHANGES.
!
!            AUTHOR:    M.DANIEL SCHWARZKOPF
!
!            DATE:      14 JULY 1983
!
!            ADDRESS: 
!
!                      G.F.D.L.
!                      P.O.BOX 308
!                      PRINCETON,N.J.08540
!                      U.S.A.
!            TELEPHONE:  (609) 452-6521
!
!-----------------------------------------------------------------------

      Use fs_profile_mod, ONLY:  pd1013,plm1013,pd810,plm810
      Use        fms_mod, ONLY:  ERROR_MESG, FATAL, WARNING, &
                                 mpp_pe, mpp_root_pe, write_version_number

implicit none
private
!      Integer,Parameter :: kind_type = selected_real_kind(15,307)

      Real, Allocatable, Dimension(:,:,:,:) :: TRNS

!-----------------------------------------------------------------------

      Public   co2int, co2int_init, co2int_end, TRNS

!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: co2int.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
 character(len=128) :: tagname = '$Name: siena_201207 $'
 logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------

      CONTAINS

!#######################################################################

      Subroutine co2int (nlev,ir,npurp,nkkk,unit1,unit2,ratio)

      Implicit None
!-----------------------------------------------------------------------
!
!      ------------   FUNCTION INTERPOLATER ROUTINE  ------------
!
!-----------------------------------------------------------------------
      Integer, Intent(IN) :: nlev,ir,npurp,nkkk,unit1,unit2
      Real,    Intent(IN) :: ratio
!-----------------------------------------------------------------------

      call error_mesg('co2int', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine co2int

!#######################################################################

      Subroutine co2int_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('co2int_init', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine co2int_init

!#######################################################################

      Subroutine co2int_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('co2int_end', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine co2int_end

!#######################################################################


      End Module co2int_mod

