      module mo_chem_utls_mod

implicit none
      private
      public :: adjh2o, inti_mr_xform, mmr2vmr, vmr2mmr, negtrc, &
                get_spc_ndx, get_het_ndx, get_extfrc_ndx, &
                has_drydep, has_srfems, get_rxt_ndx, get_grp_ndx, &
                get_grp_mem_ndx, chem_utls_init

!     save

      integer :: ox_ndx, o3_ndx, o1d_ndx, o_ndx
      logical :: do_ox

character(len=128), parameter :: version     = '$Id: mo_chem_utls.F90,v 19.0 2012/01/06 20:32:46 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine chem_utls_init (retain_cm3_bugs)
!-----------------------------------------------------------------------
!     ... Initialize the chem utils module
!-----------------------------------------------------------------------

      implicit none

      logical, intent(in) :: retain_cm3_bugs

!      ox_ndx = get_spc_ndx( 'OX' )
! for ox budget (jmao,1/7/2011)
   if (retain_cm3_bugs) then
      ox_ndx = get_spc_ndx( 'OX' )
   else
      ox_ndx = get_spc_ndx( 'O3' )
   endif
      if( ox_ndx > 0 ) then
         o3_ndx  = get_grp_mem_ndx( 'O3' )
         o1d_ndx = get_grp_mem_ndx( 'O1D' )
         o_ndx   = get_grp_mem_ndx( 'O' )
         do_ox   = o3_ndx > 0 .and. o1d_ndx > 0 .and. o_ndx > 0
      else
         o3_ndx  = 1
         o1d_ndx = 1
         o_ndx   = 1
         do_ox = .false.
      end if

      end subroutine chem_utls_init

      subroutine adjh2o( h2o, sh, mbar, vmr, do_interactive_h2o, plonl )
!-----------------------------------------------------------------------
!     ... transform water vapor from mass to volumetric mixing ratio
!-----------------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, dimension(:,:,:), intent(in)  :: vmr         ! xported species vmr
      real, dimension(:,:),   intent(in)  :: sh          ! specific humidity ( mmr )
      real, dimension(:,:),   intent(in)  :: mbar        ! atmos mean mass
      logical,                intent(in)  :: do_interactive_h2o ! include h2o sources/sinks?
      real, dimension(:,:),   intent(out) :: h2o         ! water vapor vmr

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      real, parameter :: mh2o = 1. /18.01528

      integer ::   k, ndx_ch4
      real    ::   t_value(plonl)
      integer ::   plev

      plev = SIZE(vmr,2)
!-----------------------------------------------------------------------
!       ... if not using interactive water vapor, adjust model
!           water vapor in stratosphere for source from CH4 oxidation
!-----------------------------------------------------------------------
      ndx_ch4 = get_spc_ndx( 'CH4' )
!++lwh
         do k = 1,plev
            h2o(:,k)   = mbar(:,k) * sh(:plonl,k) * mh2o
            if( .not. do_interactive_h2o .and. ndx_ch4 > 0 ) then
               t_value(:) = 6.e-6 - 2.*vmr(:,k,ndx_ch4)
!              where( t_value(:) > h2o(:,k) )
!                 h2o(:,k) = t_value(:)
!              end where
               h2o(:,k) = MAX(h2o(:,k),t_value(:))
            end if
         end do
!--lwh

      end subroutine adjh2o      

      subroutine inti_mr_xform( sh, mbar, plonl )
!-----------------------------------------------------------------
!       ... initialize mean atmospheric "wet" mass
!-----------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: sh(:,:)     ! specific humidity (kg/kg)
      real, intent(out) :: mbar(:,:)   ! mean wet atm mass ( amu )

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      real, parameter :: dry_mass = 28.966    ! amu
      real, parameter :: mfac = 1. / .622

      integer :: k
      integer :: plev
      
      plev = size(sh,2)

      do k = 1,plev
         mbar(:,k) = dry_mass
      end do

      end subroutine inti_mr_xform

      subroutine mmr2vmr( vmr, mmr, mbar, plonl )
!-----------------------------------------------------------------
!       ... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use chem_mods_mod, only : adv_mass
      use mo_grid_mod,   only : pcnstm1, pcnst

      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: mbar(:,:)
      real, intent(in)  :: mmr(:,:,:)
      real, intent(out) :: vmr(:,:,:)

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      integer :: k, m
      integer :: plev
      
      plev = size(mbar,2)

      do m = 1,pcnstm1
         if( adv_mass(m) /= 0. ) then
            do k = 1,plev
               vmr(:,k,m) = mbar(:,k) * mmr(:,k,m) / adv_mass(m)
            end do
         end if
      end do

      end subroutine mmr2vmr

      subroutine vmr2mmr( vmr, mmr, nas, grp_ratios, mbar, plonl )
!-----------------------------------------------------------------
!       ... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use chem_mods_mod, only : adv_mass, nadv_mass, grpcnt
      use mo_grid_mod,   only : pcnstm1, pcnst

      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: mbar(:,:)
      real, intent(in)  :: vmr(:,:,:)
      real, intent(out) :: mmr(:,:,:)
      real, intent(in)  :: grp_ratios(:,:,:)
      real, intent(out) :: nas(:,:,:)

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      integer :: k, m
      integer :: plev
      real    :: grp_mass(plonl)            ! weighted group mass

      plev = size(mbar,2)

!-----------------------------------------------------------------
!       ... the non-group species
!-----------------------------------------------------------------
      do m = 1,pcnstm1
         if( adv_mass(m) /= 0. ) then
            do k = 1,plev
               mmr(:,k,m) = adv_mass(m) * vmr(:,k,m) / mbar(:,k)
            end do
         end if
      end do
!-----------------------------------------------------------------
!       ... the "group" species
!-----------------------------------------------------------------
      if( do_ox ) then
         do k = 1,plev
            grp_mass(:)     = grp_ratios(:,k,o3_ndx) * nadv_mass(o3_ndx) &
                              + grp_ratios(:,k,o_ndx) * nadv_mass(o_ndx) &
                              + grp_ratios(:,k,o1d_ndx) * nadv_mass(o1d_ndx)      
            mmr(:,k,ox_ndx)  = grp_mass(:) * vmr(:,k,ox_ndx) / mbar(:,k)
            grp_mass(:)     = mmr(:,k,ox_ndx) / grp_mass(:)
            nas(:,k,o3_ndx)  = nadv_mass(o3_ndx) * grp_ratios(:,k,o3_ndx) * grp_mass(:)
            nas(:,k,o_ndx)   = nadv_mass(o_ndx) * grp_ratios(:,k,o_ndx) * grp_mass(:)
            nas(:,k,o1d_ndx) = nadv_mass(o1d_ndx) * grp_ratios(:,k,o1d_ndx) * grp_mass(:)
         end do
      end if

      end subroutine vmr2mmr

      subroutine negtrc( lat, header, fld, plonl )
!-----------------------------------------------------------------------
!       ... check for negative constituent values and
!           replace with zero value
!-----------------------------------------------------------------------

      use mo_grid_mod,    only : pcnstm1
      use m_tracname_mod, only : tracnam

      implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: lat                      ! current latitude
      integer, intent(in)          :: plonl
      character(len=*), intent(in) :: header                   ! caller tag
      real, intent(inout)          :: fld(:,:,:)               ! field to check

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter

      do m  = 1,pcnstm1
         nneg = count( fld(:,:,m) < 0. )
         if( nneg > 0 ) then
            where( fld(:,:,m) < 0. )
               fld(:,:,m) = 0.
            endwhere
!           if( pdiags%negtrc ) then
!              worst     = minval( fld(:,:,m) )
!              windex(:) = minloc( fld(:,:,m) )
!              iw        = windex(1)
!              kw        = windex(2)
!           end if
         end if
!        if( pdiags%negtrc .and. nneg > 0 ) then
!           write(*,*) header(:len(header)), tracnam(m), ' has ',nneg,' neg values'
!           write(*,*) ' worst =',worst,' @ long = ',iw,' lat = ',lat,' eta = ',kw
!        end if
      end do

      end subroutine negtrc

      integer function get_spc_ndx( spc_name )
!-----------------------------------------------------------------------
!     ... return overall species index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : pcnstm1
      use m_tracname_mod, only : tracnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_spc_ndx = -1
      do m = 1,pcnstm1
         if( trim( spc_name ) == trim( tracnam(m) ) ) then
            get_spc_ndx = m
            exit
         end if
      end do

      end function get_spc_ndx

      integer function get_grp_ndx( grp_name )
!-----------------------------------------------------------------------
!     ... return group index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : ngrp, grp_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: grp_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_grp_ndx = -1
      do m = 1,ngrp
         if( trim( grp_name ) == trim( grp_lst(m) ) ) then
            get_grp_ndx = m
            exit
         end if
      end do

      end function get_grp_ndx

      integer function get_grp_mem_ndx( mem_name )
!-----------------------------------------------------------------------
!     ... return group member index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : grpcnt
      use m_tracname_mod, only : natsnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: mem_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_grp_mem_ndx = -1
      if( grpcnt > 0 ) then
         do m = 1,max(1,grpcnt)
            if( trim( mem_name ) == trim( natsnam(m) ) ) then
               get_grp_mem_ndx = m
               exit
            end if
         end do
      end if

      end function get_grp_mem_ndx

      integer function get_het_ndx( het_name )
!-----------------------------------------------------------------------
!     ... return overall het process index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : hetcnt, het_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: het_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_het_ndx = -1
      do m = 1,max(1,hetcnt)
         if( trim( het_name ) == trim( het_lst(m) ) ) then
            get_het_ndx = m
            exit
         end if
      end do

      end function get_het_ndx

      integer function get_extfrc_ndx( frc_name )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : extcnt, extfrc_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: frc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_extfrc_ndx = -1
      if( extcnt > 0 ) then
         do m = 1,max(1,extcnt)
            if( trim( frc_name ) == trim( extfrc_lst(m) ) ) then
               get_extfrc_ndx = m
               exit
            end if
         end do
      end if

      end function get_extfrc_ndx

      integer function get_rxt_ndx( rxt_alias )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : rxt_alias_cnt, rxt_alias_lst, rxt_alias_map

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: rxt_alias

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_rxt_ndx = -1
      do m = 1,rxt_alias_cnt
         if( trim( rxt_alias ) == trim( rxt_alias_lst(m) ) ) then
            get_rxt_ndx = rxt_alias_map(m)
            exit
         end if
      end do

      end function get_rxt_ndx

      logical function has_drydep( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species dry deposition
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : drydep_cnt, drydep_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_drydep = .false.
      do m = 1,drydep_cnt
         if( trim( spc_name ) == trim( drydep_lst(m) ) ) then
            has_drydep = .true.
            exit
         end if
      end do

      end function has_drydep

      logical function has_srfems( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species surface emission
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : srfems_cnt, srfems_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_srfems = .false.
      do m = 1,srfems_cnt
         if( trim( spc_name ) == trim( srfems_lst(m) ) ) then
            has_srfems = .true.
            exit
         end if
      end do

      end function has_srfems

      end module mo_chem_utls_mod
