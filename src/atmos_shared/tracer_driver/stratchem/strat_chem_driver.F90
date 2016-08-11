     module strat_chem_driver_mod


use mpp_mod, only: input_nml_file 
use              fms_mod, only : file_exist, &
                                 check_nml_error,  &
                                 close_file, open_namelist_file, &
                                 stdlog, write_version_number, &
                                 error_mesg, FATAL
use              mpp_io_mod, only: mpp_open, mpp_close, &
                       MPP_NATIVE, MPP_RDONLY, MPP_DELETE

use   tracer_manager_mod, only : get_tracer_index, NO_TRACER
use    field_manager_mod, only : MODEL_ATMOS

use                mpp_mod, only: mpp_pe, mpp_root_pe, stdout
use constants_mod, only : PI, TFREEZE, GRAV, PSTD_MKS, RDGAS

use STRAT_CHEM_MOD, only : chemistry, zen2, dcly_dt, sediment
     implicit none

private
!----------- ****** VERSION NUMBER ******* ---------------------------

character(len=128)  :: version =  '$Id: strat_chem_driver.F90,v 19.0 2012/01/06 20:32:12 fms Exp $'
character(len=128)  :: tagname =  '$Name: tikal $'
logical             :: module_is_initialized = .FALSE.

!-------  interfaces --------

public   strat_chem_driver_init,   strat_chem, strat_chem_driver_end


!----------- namelist -------------------


integer     :: fixed_chem_year= 0   ! 

real        :: n2o_foto_accel= 1.0     !      

logical     :: do_coupled_stratozone = .false.  ! Do I want to use this routine?

namelist / strat_chem_nml /         &
                   fixed_chem_year, n2o_foto_accel, do_coupled_stratozone 





     logical               :: run_startup = .true.
     real  ::   ozon(11,48),cosp(14),cosphc(48),photo(132,14,11,48),    &
           solardata(1801),chlb(90,15),ozb(144,90,12),tropc(151,9),    &
           dfdage(90,48,8),anoy(90,48)
     integer :: mype

!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. 
!  This should be initialized to zero. 
!
!-----------------------------------------------------------------------

integer :: nsphum = 0
integer :: nliq_wat = 0
integer :: nice_wat = 0
integer :: ncld_amt = 0
integer :: nhno3 = 0
integer :: nn2o5 = 0
integer :: nh2o2 = 0
integer :: nhcl = 0
integer :: nhocl = 0
integer :: nclono2 = 0
integer :: nh2co = 0
integer :: noy = 0
integer :: nhobr = 0
integer :: nhno4 = 0
integer :: nhbr = 0
integer :: nbrono2 = 0
integer :: nch3ooh = -1
integer :: nco = 0
integer :: nnoy = 0
integer :: ncly = 0
integer :: nbry = 0
integer :: nch4 = 0
integer :: nstrath2o = 0
integer :: nn2o = 0
integer :: nage = 0
integer :: no3 = 0
integer :: no3ch = 0
integer :: nextinct = 0
integer :: naerosol = 0

     CONTAINS

function strat_chem_driver_init()
 logical :: strat_chem_driver_init

      integer                 :: unit, ierr, io, logunit
!---------------------------------------------------------------------
!    read strat_chem namelist.
!---------------------------------------------------------------------


         if (file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
           read (input_nml_file, nml=strat_chem_nml, iostat=io)
           ierr = check_nml_error(io,'strat_chem_nml')
#else
           unit =  open_namelist_file ( )
           ierr=1; do while (ierr /= 0)
           read (unit, nml=strat_chem_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'strat_chem_nml')
           enddo
 10        call close_file (unit)
#endif
         endif

     strat_chem_driver_init = do_coupled_stratozone
     
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number(version, tagname)
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe()) write (logunit, nml=strat_chem_nml)
 
      module_is_initialized = .true.

      if (.not. do_coupled_stratozone) return
!---------------------------------------------------------------------

      call chem_startup
      
end function strat_chem_driver_init

subroutine strat_chem_driver_end

      module_is_initialized = .false.

end subroutine strat_chem_driver_end

     subroutine chem_startup()
!
! Obtains chemical lower boundary condition and arrays for photolysis 
! rate calculation
!



      real age(90,48,12)
      integer lev,ipz,icz,lv,nc,jl,ir
!
!   local variables: 

      integer                 :: unit, outunit

      if(run_startup) then
         run_startup = .false.


      nsphum    = get_tracer_index(MODEL_ATMOS,'sphum')    !Tracer #1
      nliq_wat  = get_tracer_index(MODEL_ATMOS,'liq_wat')   !Tracer #2
      nice_wat  = get_tracer_index(MODEL_ATMOS,'ice_wat')   !Tracer #3
      ncld_amt  = get_tracer_index(MODEL_ATMOS,'cld_amt')   !Tracer #4
      nhno3     = get_tracer_index(MODEL_ATMOS,'HNO3')      !Tracer #5
      nn2o5     = get_tracer_index(MODEL_ATMOS,'N2O5')      !Tracer #6
      nh2o2     = get_tracer_index(MODEL_ATMOS,'H2O2')      !Tracer #7
      nhcl      = get_tracer_index(MODEL_ATMOS,'HCl')       !Tracer #8
      nhocl     = get_tracer_index(MODEL_ATMOS,'HOCl')      !Tracer #9
      nclono2   = get_tracer_index(MODEL_ATMOS,'ClONO2')    !Tracer #10
      nh2co     = get_tracer_index(MODEL_ATMOS,'H2CO')      !Tracer #11
      noy       = get_tracer_index(MODEL_ATMOS,'Oy')        !Tracer #12
      nhobr     = get_tracer_index(MODEL_ATMOS,'HOBr')      !Tracer #13
      nhno4     = get_tracer_index(MODEL_ATMOS,'HNO4')      !Tracer #14
      nhbr      = get_tracer_index(MODEL_ATMOS,'HBr')       !Tracer #15
      nbrono2   = get_tracer_index(MODEL_ATMOS,'BrONO2')    !Tracer #16
      nch3ooh   = get_tracer_index(MODEL_ATMOS,'CH3OOH')    !Tracer #17
      nco       = get_tracer_index(MODEL_ATMOS,'CO')        !Tracer #18
      nnoy      = get_tracer_index(MODEL_ATMOS,'NOy')       !Tracer #19
      ncly      = get_tracer_index(MODEL_ATMOS,'Cly')       !Tracer #20
      nbry      = get_tracer_index(MODEL_ATMOS,'Bry')       !Tracer #21
      nch4      = get_tracer_index(MODEL_ATMOS,'CH4')       !Tracer #22
      nstrath2o = get_tracer_index(MODEL_ATMOS,'StratH2O')  !Tracer #23
      nn2o      = get_tracer_index(MODEL_ATMOS,'N2O')       !Tracer #24
      nage      = get_tracer_index(MODEL_ATMOS,'Age')       !Tracer #25
      no3       = get_tracer_index(MODEL_ATMOS,'O3')        !Tracer #26
      no3ch     = get_tracer_index(MODEL_ATMOS,'O3_chem')   !Tracer #27
      nextinct  = get_tracer_index(MODEL_ATMOS,'Extinction')!Tracer #28
      naerosol  = get_tracer_index(MODEL_ATMOS,'Aerosol')   !Tracer #29

      if ( nliq_wat  == NO_TRACER .or. nice_wat  == NO_TRACER .or. &
           ncld_amt  == NO_TRACER .or. nhno3     == NO_TRACER .or. &
           nn2o5     == NO_TRACER .or. nh2o2     == NO_TRACER .or. &
           nhcl      == NO_TRACER .or. nhocl     == NO_TRACER .or. &
           nclono2   == NO_TRACER .or. nh2co     == NO_TRACER .or. &
           noy       == NO_TRACER .or. nhobr     == NO_TRACER .or. &
           nhno4     == NO_TRACER .or. nhbr      == NO_TRACER .or. &
           nbrono2   == NO_TRACER .or. nch3ooh   == NO_TRACER .or. &
           nco       == NO_TRACER .or. nnoy      == NO_TRACER .or. &
           ncly      == NO_TRACER .or. nbry      == NO_TRACER .or. &
           nch4      == NO_TRACER .or. nstrath2o == NO_TRACER .or. &
           nn2o      == NO_TRACER .or. nage      == NO_TRACER .or. &
           no3       == NO_TRACER .or. no3ch     == NO_TRACER ) &
         call error_mesg ('Strat_chem_driver', &
                          'A necessary tracer is missing from the field_table &
                          &and thus will not allow strat_chem to run correctly.', FATAL)  

!  read in chemical lower boundary 
!  
         outunit = stdout()
         call mpp_open( unit, 'INPUT/chemlbf',action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/chemlbf'
         DO NC = 1,15                                           
           READ(unit,'(6E13.6)') (CHLB(JL,NC),JL=1,90) 
         ENDDO                                                          
         READ(unit,'(6E13.6)') OZB  
         read(unit,'(6e13.6)') tropc
         call mpp_close(unit)
!
!  read in photolysis files
!         
         call mpp_open( unit, 'INPUT/photolsmax', action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/photolsmax'
         DO LEV = 1,48                                            
         DO IPZ = 1,11                                               
         DO ICZ = 1,14                                               
         READ(unit,'(I4,E12.4,2F10.4,5(/6E12.4),/3E12.4)') &
           LV,OZON(IPZ,LEV),COSP(ICZ),COSPHC(LEV), &
           (PHOTO(IR,ICZ,IPZ,LEV),IR=1,33)
         READ(unit,'(5(6E12.4/),3E12.4)')   &
           (PHOTO(IR,ICZ,IPZ,LEV),IR=34,66)
         enddo
         enddo
         enddo
         call mpp_close(unit)
         call mpp_open( unit, 'INPUT/photolsmin', action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/photolsmin'
         DO LEV = 1,48                                            
         DO IPZ = 1,11                                               
         DO ICZ = 1,14                                               
         READ(unit,'(I4,E12.4,2F10.4,5(/6E12.4),/3E12.4)') &
           LV,OZON(IPZ,LEV),COSP(ICZ),COSPHC(LEV), &
           (PHOTO(IR,ICZ,IPZ,LEV),IR=67,99)
         READ(unit,'(5(6E12.4/),3E12.4)')   &
           (PHOTO(IR,ICZ,IPZ,LEV),IR=100,132)
         enddo
         enddo
         enddo
         call mpp_close(unit)
         call mpp_open( unit, 'INPUT/solar_f107.dat', action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/solar_f107.dat'
         read(unit,'(f6.0,5f7.0)') solardata
         call mpp_close(unit)
         DO LEV = 1,48                                              
         DO IPZ = 1,11                                               
         DO ICZ = 1,14                                               
         DO IR = 1,132                                                
         PHOTO(IR,ICZ,IPZ,LEV) = ALOG(PHOTO(IR,ICZ,IPZ,LEV)+1.E-30)     
         enddo
         enddo
         enddo
         enddo
!
!  read in data for Cly and Bry computation
!         
         call mpp_open( unit, 'INPUT/dfdage.dat', action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/dfdage.dat'
         read(unit,'(6e13.6)') age
         read(unit,'(6e13.6)') dfdage
         call mpp_close(unit)
!
!  read in data for NOy tropospheric relaxation
!         
         call mpp_open( unit, 'INPUT/noy_annual.dat', action=MPP_RDONLY )
         if (mpp_pe() == mpp_root_pe()) WRITE(outunit,*) 'INPUT/noy_annual.dat'
         read(unit,'(6e13.6)') anoy
         call mpp_close(unit)
     endif



     return
     end   subroutine chem_startup

     subroutine strat_chem(alon,alat,chems,dchems,pfull,temp,itime,    &
       is,ie,js,je,dt,coszen,&
       nsphum,chem_tend,ozone,o3_prod,aerosol,mype)
     real, intent(in), dimension(:,:,:,:)  :: chems
     real, intent(in), dimension(:,:,:,:)  :: dchems
     real, intent(in), dimension(:,:,:)    :: pfull, temp 
     real, intent(in), dimension(:,:)      :: alon,alat,coszen
     real, intent(in)                      :: dt
     real, intent(out), dimension(:,:,:,:) :: chem_tend
     real, intent(out), dimension(:,:,:)   :: ozone, o3_prod, aerosol
     integer, intent(in)                   :: is,ie,js,je,nsphum,mype
     integer, intent(inout), dimension(6)  :: itime
!
     integer :: nc,il,jl,kl,merids,lats,levs
     integer :: ipts           
!
     real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) ::     &
                               ozcol,anat,aice,alats,cly,bry,age
     real, dimension(size(pfull,1)) :: h2o,h2o_tend,dagesq

     real  ::   cozen(ie-is+1,je-js+1),vtemp(size(pfull,1)),           &
                rho(size(pfull,1)),dy(size(pfull,1),21),               &
                ch_tend(size(pfull,1),21),                             &
                drad,darg,delta,cdx,cozi,cozj1,cozj2,const,const2,xc,  &
                vtau10,fact2,cond,extinct
!
     if(.not. do_coupled_stratozone) return ! You do not want to do stratospheric chemistry
     merids = size(pfull,1)
     lats = size(pfull,2)
     levs = size(pfull,3)
     drad = 3.141592653589793/180.        
     const =  287.056*273.16/(1.01325*9.81)                 
     const2 = 2.693e19 *273.16/101325.

     chem_tend(:,:,:,:) = 0.0
!wfc Comment out for now
!     drad = PI/180.        
!     const =  RDGAS*TFREEZE/(PSTD_MKS*1.e-5*GRAV)                 
!     const2 = 2.693e19 *TFREEZE/PSTD_MKS
!

!  change model time to fit in with chemical time
      
    if( fixed_chem_year > 0 )  itime(1) = fixed_chem_year
!    if(mpp_pe()==mpp_root_pe()) print *,' time', itime

!
! Compute cosine of solar zenith angle
!
     call zen2 (itime,dt,drad,darg,delta,cdx)
     do jl = 1,je-js+1     
     do il = 1,ie-is+1
     cozi = cos(darg+alon(il,jl))
     cozj1 = cos(alat(il,jl))*cos(delta)                                 
     cozj2 = sin(alat(il,jl))*sin(delta)             
     cozen(il,jl) = cozj2 + cozj1*cozi
     enddo
     enddo
!
! Calculate overhead ozone columns for photolysis
! Allow for transport problems with a 1 e-15 minimum for ozone.
!
     do jl = 1,je-js+1
     do il = 1,ie-is+1
     xc = chems(il,jl,1,no3) 
     if(xc.lt.1.0e-15) xc = 1.0e-15
     ozcol(il,jl,1) = xc*const*pfull(il,jl,1)
     do kl = 2,levs
     xc =  sqrt(chems(il,jl,kl,no3)*chems(il,jl,kl-1,no3)) 
     if(xc.lt.1.0e-15) xc = 1.0e-15
     ozcol(il,jl,kl) = ozcol(il,jl,kl-1) +                    &
        const*(pfull(il,jl,kl) - pfull(il,jl,kl-1))*xc         
     enddo
     enddo
     enddo
!=========================================================================
! Do main loop over model levels
!
!==========================================================================

     do 1000 kl = 1,levs
     do jl = js,je    
!
!  Set up inverse temperature and density for use in the chemical model
!  CHEMISTRY computes the chemical change for one model timestep, one latitude
!  row at a time and converts this into a rate of change chem_tend
!  chems isn't necessarily positive definite, depending on the transport 
!  scheme, so the values are set to zero within the chemistry scheme
!  and the concentration is relaxed to zero with a 1-day timescale 
!
     do il = is,ie
     vtemp(il) = 1.0/temp(il,jl,kl)
     rho(il) = const2*pfull(il,jl,kl)*vtemp(il)
     dy(il,1 ) = chems(il,jl,kl,nhno3    ) + dchems(il,jl,kl,nhno3    )*dt
     dy(il,2 ) = chems(il,jl,kl,nn2o5    ) + dchems(il,jl,kl,nn2o5    )*dt
     dy(il,3 ) = chems(il,jl,kl,nh2o2    ) + dchems(il,jl,kl,nh2o2    )*dt
     dy(il,4 ) = chems(il,jl,kl,nhcl     ) + dchems(il,jl,kl,nhcl     )*dt
     dy(il,5 ) = chems(il,jl,kl,nhocl    ) + dchems(il,jl,kl,nhocl    )*dt
     dy(il,6 ) = chems(il,jl,kl,nclono2  ) + dchems(il,jl,kl,nclono2  )*dt
     dy(il,7 ) = chems(il,jl,kl,nh2co    ) + dchems(il,jl,kl,nh2co    )*dt
     dy(il,8 ) = chems(il,jl,kl,noy      ) + dchems(il,jl,kl,noy      )*dt
     dy(il,9 ) = chems(il,jl,kl,nhobr    ) + dchems(il,jl,kl,nhobr    )*dt
     dy(il,10) = chems(il,jl,kl,nhno4    ) + dchems(il,jl,kl,nhno4    )*dt
     dy(il,11) = chems(il,jl,kl,nhbr     ) + dchems(il,jl,kl,nhbr     )*dt
     dy(il,12) = chems(il,jl,kl,nbrono2  ) + dchems(il,jl,kl,nbrono2  )*dt
     dy(il,13) = chems(il,jl,kl,nch3ooh  ) + dchems(il,jl,kl,nch3ooh  )*dt
     dy(il,14) = chems(il,jl,kl,nco      ) + dchems(il,jl,kl,nco      )*dt
     dy(il,15) = chems(il,jl,kl,nnoy     ) + dchems(il,jl,kl,nnoy     )*dt
     dy(il,16) = chems(il,jl,kl,ncly     ) + dchems(il,jl,kl,ncly     )*dt
     dy(il,17) = chems(il,jl,kl,nbry     ) + dchems(il,jl,kl,nbry     )*dt
     dy(il,18) = chems(il,jl,kl,nch4     ) + dchems(il,jl,kl,nch4     )*dt
     dy(il,19) = chems(il,jl,kl,nstrath2o) + dchems(il,jl,kl,nstrath2o)*dt
     dy(il,20) = chems(il,jl,kl,nn2o     ) + dchems(il,jl,kl,nn2o     )*dt
     dy(il,21) = chems(il,jl,kl,nage     ) + dchems(il,jl,kl,nage     )*dt
     h2o_tend(il) = 0.0
     do nc = 1,21
     ch_tend(il,nc) = 0.0
     if(dy(il,nc).lt.0.0) then
        ch_tend(il,nc) = -dy(il,nc)/86400.0
        dy(il,nc) = 0.0
     endif
     enddo
     if(dy(il,8).lt.1.0e-15) dy(il,8) = 1.0e-15
     if(dy(il,15).lt.1.0e-15) dy(il,15) = 1.0e-15
!
! Set water vapour values from the main model after converting to vmr or
! use stratospheric h2o: chems(...,23), which includes sedimented terms. 
! 
!  Note: values are separated into gas and solid phases in the heterogeneous 
!  chemistry subroutine
!
     h2o(il) = chems(il,jl,kl,nsphum)*1.61
!     cond = chems(il,jl,kl,2) + chems(il,jl,kl,3)
      cond = chems(il,jl,kl,nliq_wat) + chems(il,jl,kl,nice_wat)
    h2o(il) = h2o(il) + cond*1.61
     if(h2o(il) < 1.0e-7) h2o(il) = 1.0e-7
!     dagesq(il) = (chems(il,jl,levs,25) - chems(il,jl,kl,25))**2
     dagesq(il) = (chems(il,jl,levs,nage) - chems(il,jl,kl,nage))**2
!
!  No more than 5 ppmv water in the extra-tropical lower stratosphere
!
     if(h2o(il) > 5.0e-6.and.dagesq(il) > 0.01.and.pfull(il,jl,kl) > 1.0e4) &
         h2o(il) = 5.0e-6
     ozone(il,jl,kl) = chems(il,jl,kl,no3)!26)
     o3_prod(il,jl,kl) = chems(il,jl,kl,no3ch)!27)
     extinct = 0.0
     if (nextinct > 0 ) extinct = chems(il,jl,kl,nextinct)*1000.0
     if(extinct.eq.0.0) then
       aerosol(il,jl,kl) = 0.0
     elseif(extinct <= 4.0e-3) then
       aerosol(il,jl,kl) = 4.25e-6*extinct**0.68
     elseif(extinct > 4.0e-3.and.extinct <= 2.0e-2) then
       aerosol(il,jl,kl) = 1.223e-5*extinct**0.875
     elseif(extinct > 2.0e-2) then
       aerosol(il,jl,kl) = 2.0e-5*extinct 
     endif
     enddo
     call chemistry (alon(:,jl),alat(:,jl),jl,kl,dy,h2o,dagesq,ozcol(:,jl,kl),  &
       pfull(:,jl,kl),rho,temp(:,jl,kl),vtemp,cozen(:,jl),cdx,chlb,       &
       ozb(:,:,itime(2)),anat(:,jl,kl),aice(:,jl,kl),photo,solardata,     &
       ozon,cosp,cosphc,anoy,aerosol(:,jl,kl),dt,merids,ch_tend,   &
       ozone(:,jl,kl), o3_prod(:,jl,kl),h2o_tend,mype,itime)

     do il = 1,merids

!     do nc = 1,21
!     chem_tend(il,jl,kl,nc+4     ) = ch_tend(il,nc)
     chem_tend(il,jl,kl,nhno3    ) = ch_tend(il,1 )
     chem_tend(il,jl,kl,nn2o5    ) = ch_tend(il,2 )
     chem_tend(il,jl,kl,nh2o2    ) = ch_tend(il,3 )
     chem_tend(il,jl,kl,nhcl     ) = ch_tend(il,4 )
     chem_tend(il,jl,kl,nhocl    ) = ch_tend(il,5 )
     chem_tend(il,jl,kl,nclono2  ) = ch_tend(il,6 )
     chem_tend(il,jl,kl,nh2co    ) = ch_tend(il,7 )
     chem_tend(il,jl,kl,noy      ) = ch_tend(il,8 )
     chem_tend(il,jl,kl,nhobr    ) = ch_tend(il,9 )
     chem_tend(il,jl,kl,nhno4    ) = ch_tend(il,10)
     chem_tend(il,jl,kl,nhbr     ) = ch_tend(il,11)
     chem_tend(il,jl,kl,nbrono2  ) = ch_tend(il,12)
     chem_tend(il,jl,kl,nch3ooh  ) = ch_tend(il,13)
     chem_tend(il,jl,kl,nco      ) = ch_tend(il,14)
     chem_tend(il,jl,kl,nnoy     ) = ch_tend(il,15)
     chem_tend(il,jl,kl,ncly     ) = ch_tend(il,16)
     chem_tend(il,jl,kl,nbry     ) = ch_tend(il,17)
     chem_tend(il,jl,kl,nch4     ) = ch_tend(il,18)
     chem_tend(il,jl,kl,nstrath2o) = ch_tend(il,19)
     chem_tend(il,jl,kl,nn2o     ) = ch_tend(il,20)
     chem_tend(il,jl,kl,nage     ) = ch_tend(il,21)
!     enddo
     chem_tend(il,jl,kl,nsphum) = h2o_tend(il)/1.61
!     chem_tend(il,jl,kl,2:4) = 0.0
     chem_tend(il,jl,kl,nliq_wat) = 0.0
     chem_tend(il,jl,kl,nice_wat) = 0.0
     chem_tend(il,jl,kl,ncld_amt) = 0.0
!
!  Relax tracer 23 to the total vmr of PSC
!  Relax Oy to 2e-6 if it exceeds this value in the upper mesosphere
!
     chem_tend(il,jl,kl,nstrath2o) = (aice(il,jl,kl) + 3.0*anat(il,jl,kl)  -  &
                 chems(il,jl,kl,nstrath2o))/21600.0 !23 <- nstrath2o
 !    if(pfull(il,jl,kl) < 10.0.and.chems(il,jl,kl,noy) > 2.0e-6)    &
 !        chem_tend(il,jl,kl,noy) = (2.0e-6  - chems(il,jl,kl,noy))/21600.0
     enddo
     enddo
1000 continue
! 
! compute rate of change of Cly and Bry
!
     fact2 = 1.00
     do kl = 1,levs
     alats(:,:,kl) = 1.0 + (90.0 - alat(:,:))*89.0/180.0
     enddo
     age(:,:,:) = chems(:,:,:,nage)!25)
     cly(:,:,:) = chems(:,:,:,ncly)!20)
     bry(:,:,:) = chems(:,:,:,nbry)!21)
     CALL DCLY_DT(age,dfdage,tropc,alats,cly,bry,chem_tend(:,:,:,ncly),  &
        chem_tend(:,:,:,nbry),fact2,merids,lats,levs,itime)
!
!  set rates of change of Cly and Bry to zero at night, and double during
!  the daytime
!
     do jl = 1,je-js+1
     do il = 1,merids
     if(cozen(il,jl).gt.0.0) then
         chem_tend(il,jl,:,ncly) = 2.*chem_tend(il,jl,:,ncly)
         chem_tend(il,jl,:,nbry) = 2.*chem_tend(il,jl,:,nbry)
     else
         chem_tend(il,jl,:,ncly) = 0.0
         chem_tend(il,jl,:,nbry) = 0.0
     endif
     enddo
     enddo
!
!  Rainout of HNO3, Cly and Bry in the troposphere; relax age of air to zero 
!  also in troposphere, defined by delta Age < 0.1
!
     vtau10 = 1.0/(10.0*86400.0)
     do kl = 1,levs
     do jl = 1,je-js+1
     do il = 1,merids
     dagesq(il) = (chems(il,jl,levs,nage) - chems(il,jl,kl,nage))**2
     if (dagesq(il).lt.1.0e-2)  then
          chem_tend(il,jl,kl,ncly) = (1.0e-13 - chems(il,jl,kl,ncly))*vtau10
          chem_tend(il,jl,kl,nbry) = (1.0e-15 - chems(il,jl,kl,nbry))*vtau10
          chem_tend(il,jl,kl,nage) = - chems(il,jl,kl,nage)*vtau10
     endif
     enddo
     enddo
     enddo
!
! sediment nat and ice
!
     ipts = merids*lats
     CALL SEDIMENT(anat,aice,chem_tend(:,:,:,nnoy),chem_tend(:,:,:,nstrath2o),   &
         chem_tend(:,:,:,nhno3),ipts,levs,pfull,dt,mype)
     return
     end subroutine strat_chem

     end module strat_chem_driver_mod
