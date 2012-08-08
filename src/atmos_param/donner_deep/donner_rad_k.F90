
!VERSION NUMBER:
!  $Id: donner_rad_k.F90,v 16.0 2008/07/30 22:07:00 fms Exp $

!module donner_rad_inter_mod

!#include "donner_rad_interfaces.h"

!end module donner_rad_inter_mod


!###################################################################

subroutine don_r_donner_rad_driver_k   &
         (isize, jsize, nlev_lsm, Param, Col_diag, Initialized,   &
          pfull, temp, land, exit_flag, Don_conv, Don_rad, Nml, ermesg, error)

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_initialized_type, donner_conv_type,&
                             donner_column_diag_type, donner_rad_type

implicit none

!----------------------------------------------------------------------
integer,                               intent(in)    :: isize, jsize,   &
                                                        nlev_lsm
type(donner_param_type),               intent(in)    :: Param
type(donner_column_diag_type),         intent(in)    :: Col_diag
type(donner_initialized_type),         intent(in)    :: Initialized
real, dimension(isize,jsize,nlev_lsm), intent(in)    :: pfull, temp
real, dimension(isize,jsize),          intent(in)    :: land
logical, dimension(isize,jsize),       intent(in)    :: exit_flag
type(donner_conv_type),                intent(inout) :: Don_conv
type(donner_rad_type),                 intent(inout) :: Don_rad
type(donner_nml_type),                 intent(inout) :: Nml       
character(len=*),                      intent(out)   :: ermesg
integer,                               intent(out)   :: error

      ermesg= ' ' ; error = 0

!---------------------------------------------------------------------
!    call define_ice_size to define the ice particle size distribution 
!    within the anvil (Don_conv%dgeice).
!---------------------------------------------------------------------
      call don_r_define_ice_size_k    &
           (isize, jsize, nlev_lsm, Param, Col_diag, pfull, &
            Don_conv%xice, Don_conv%przm, Don_conv%prztm,    &
            Don_conv%dgeice, ermesg, error)

!---------------------------------------------------------------------
!    call define_cell_liquid_size to compute the cell liquid effective
!    droplet diameter (Don_conv%cell_liquid_eff_diam) and the number of
!    cell droplets (Don_rad%cell_droplet_number).
!---------------------------------------------------------------------
      call don_r_define_cell_liquid_size_k   &
           (isize, jsize, nlev_lsm, Param, Nml, Col_diag, pfull, temp, &
            Don_conv%cuql, Don_conv%cual, land, exit_flag,    &
            Don_conv%cell_liquid_eff_diam, Don_rad%cell_droplet_number,&
            ermesg, error)
 
!---------------------------------------------------------------------
!    since liquid is not allowed in anvil currently, set the droplet
!    number there to be 0.0.
!---------------------------------------------------------------------
      Don_rad%meso_droplet_number = 0.0

!---------------------------------------------------------------------
!    save  variables needed by the radiation package.
!---------------------------------------------------------------------
      call don_r_donner_deep_sum_k  &
           (isize, jsize, nlev_lsm, Param, Nml, Initialized, &
            Don_conv%xliq, Don_conv%xice, Don_conv%cual,  &
            Don_conv%ampta1, Don_conv%cuql, Don_conv%cuqi, &
            Don_conv%dgeice, Don_conv%cell_liquid_eff_diam, &
            Don_rad, ermesg, error)

!--------------------------------------------------------------------


end subroutine don_r_donner_rad_driver_k

!######################################################################

subroutine don_r_define_ice_size_k    &
         (isize, jsize, nlev_lsm, Param, Col_diag, pfull, xice, przm, &
          prztm, dgeice, ermesg, error)
  
!----------------------------------------------------------------------
!    subroutine define_ice_size obtains the effective ice crystal size 
!    profile (dgeice).
!----------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_column_diag_type

implicit none

!----------------------------------------------------------------------
integer,                               intent(in)  :: isize, jsize,  &
                                                      nlev_lsm
type(donner_param_type),               intent(in)  :: Param
type(donner_column_diag_type),         intent(in)  :: Col_diag
real, dimension(isize,jsize,nlev_lsm), intent(in)  :: pfull, xice
real, dimension(isize,jsize),          intent(in)  :: przm, prztm
real, dimension(isize,jsize,nlev_lsm), intent(out) :: dgeice
character(len=*),                      intent(out) :: ermesg
integer,                               intent(out) :: error
            
!---------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure at model full levels [ Pa ]
!     xice           mesoscale ice content [ kg(ice) / kg(air ]
!     przm           pressure at anvil base [ Pa ]
!     prztm          pressure at anvil top [ Pa ]
!
!   intent(out) variables:
!
!     dgeice         effective ice crystal size [ microns ]
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:
  
      integer :: i, j, k        ! do-loop indices

      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    call subroutine andge to assign effective sizes (dgeice) to the 
!    ice in the model layers of the anvil.
!    make sure the ice size is within the limits acceptable to the 
!    radiative properties parameterizations.
!---------------------------------------------------------------------
      do k=1,nlev_lsm                           
        do j=1,jsize
          do i=1,isize
            if (xice(i,j,k) > 0.0) then
              call don_r_andge_k   &
                   (i, j, Param, Col_diag, pfull(i,j,k), przm(i,j),  &
                    prztm(i,j), dgeice(i,j,k), ermesg, error)
            else
              dgeice(i,j,k) = 0.
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------


end subroutine don_r_define_ice_size_k


!######################################################################

subroutine don_r_andge_k   &
         (i, j, Param, Col_diag, press, pzm, pztm, dgeicer, ermesg, error)
  
!---------------------------------------------------------------------
!    subroutine andge defines the generalized effective ice crystal size
!    for the mesoscale anvil layers.
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_column_diag_type

implicit none

!---------------------------------------------------------------------
integer,                       intent(in)  :: i,j 
type(donner_param_type),       intent(in)  :: Param
type(donner_column_diag_type), intent(in)  :: Col_diag
real,                          intent(in)  :: press, pzm, pztm       
real,                          intent(out) :: dgeicer    
character(len=*),              intent(out) :: ermesg
integer,                       intent(out) :: error
 
!---------------------------------------------------------------------
!   intent(in) variables:
!
!     i, j      physics window indices of the current column
!     press     pressure at model full levels [ Pa ]
!     pzm       pressure at base of mesoscale anvil [ Pa ]
!     pztm      pressure at top of mesoscale anvil [ Pa ]
!
!   intent(out) variables:
! 
!     dgeicer   generalized effective size of ice crystals in anvil 
!               defined as in Fu (1996, J. Clim.). [ microns ]
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real     :: znor    ! normalized distance from anvil base
                          ! [ dimensionless ]
      integer  :: k, n    ! do-loop indices

      ermesg = ' ' ; error = 0

!-------------------------------------------------------------------
!    be sure that anvil base has higher pressure than anvil top. 
!-------------------------------------------------------------------
      if (pzm < pztm) then
        ermesg = ' andge: pzm is < pztm'
        error = 1
        return
      endif

!-------------------------------------------------------------------
!    define the relative displacement of the current pressure between
!    the anvil top and bottom. avoid calculation if anvil depth is 0.0
!    (implying a one-layer thick anvil).
!-------------------------------------------------------------------
      if (pzm == pztm) then
        znor = 0.5
      else
        znor = (pzm - press)/(pzm - pztm)
      endif

!--------------------------------------------------------------------
!    define the value of dgeice at the appropriate relative displacement
!    from anvil base (dgeicer).
!--------------------------------------------------------------------
      do k=2,Param%anvil_levels
        if ((znor >= Param%relht(k-1)) .and.   &
            (znor <= Param%relht(k))) then
          dgeicer = Param%dgeice(k-1) + ((znor - Param%relht(k-1))*  &
                   (Param%dgeice(k) - Param%dgeice(k-1))/  &
                   (Param%relht(k) - Param%relht(k-1)))
          exit
        endif
      end do

!---------------------------------------------------------------------
!   if in diagnostics column, output relevant variables related to the 
!   effective ice size calculation.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          if (j == Col_diag%j_dc(n) .and. i == Col_diag%i_dc(n)) then
            write (Col_diag%unit_dc(n), '(a, e22.12)')  'znor= ',znor
            write (Col_diag%unit_dc(n), '(a, 2e22.12)') &
                            'relhts= ',Param%relht(k-1),Param%relht(k)
            write (Col_diag%unit_dc(n), '(a, e22.12)') &
                             'dgeicer= ',dgeicer
          endif
        end do
      endif

!--------------------------------------------------------------------


end subroutine don_r_andge_k


!###################################################################

subroutine don_r_define_cell_liquid_size_k   &
         (isize, jsize, nlev_lsm, Param, Nml, Col_diag, pfull, temp, &
          cuql, cual, land, exit_flag, cell_liquid_eff_diam,  &
          cell_droplet_number, ermesg, error)

!--------------------------------------------------------------------
!    subroutine define_cell_liquid_size calculates the effective radii 
!    of liquid cloud drops in cumulus clouds following the prescription
!    of Bower et al (JAS, 1994).
!---------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_column_diag_type

implicit none

!---------------------------------------------------------------------
integer,                         intent(in)  :: isize, jsize, nlev_lsm
type(donner_param_type),         intent(in)  :: Param
type(donner_nml_type),           intent(in)  :: Nml

type(donner_column_diag_type),   intent(in)  :: Col_diag
real, dimension(isize,jsize,nlev_lsm),                 &
                                 intent(in)  :: pfull, temp, cuql, cual
real, dimension(isize,jsize),    intent(in)  :: land
logical,dimension(isize,jsize),  intent(in)  :: exit_flag
real, dimension(isize,jsize,nlev_lsm),                  &
                                 intent(out) :: cell_liquid_eff_diam, &
                                                cell_droplet_number
character(len=*),                intent(out) :: ermesg
integer,                         intent(out) :: error
   
!--------------------------------------------------------------------
!   intent(in) variables:
!
!     pfull          pressure at model full levels [ Pa ]
!     temp           temperature at model full levels [ deg K ]
!     cuql           cell liquid content [ kg(h2o) / kg(air) ]
!     land           fraction of land in grid box [ fraction ]
!     cell_liquid_eff_diam
!                    effective diameter of liquid cloud drops 
!                    [ microns ]
!     cell_droplet_number
!                    droplet number in cells [ # / kg(air) ]
!
!------------------------------------------------------------------
! local variables
  
      real, dimension (isize,jsize) ::   &
                                cell_pbase, temp_cell_pbase, &
                                cell_land_ref_delp, cell_ocean_ref_delp
                                       
      real, dimension (isize,jsize,nlev_lsm)          ::   &
                               cell_delp, cell_liquid_eff_diam_land,  &
                               cell_liquid_eff_diam_ocean, &
                               cell_droplet_number_land,   &
                               cell_droplet_number_ocean
      integer   :: i, j, k, n

!------------------------------------------------------------------
! local variables
!
!        cell_pbase                      pressure at cloud base [ Pa ]
!        temp_cell_pbase                 temperature at cloud base
!                                        [ deg K ]
!        cell_land_ref_delp              pressure difference between
!                                        cloud base and a point 
!                                        delz_land meters above cloud 
!                                        base
!        cell_ocean_ref_delp             pressure difference between
!                                        cloud base and a point 
!                                        delz_ocean meters above cloud 
!                                        base
!        cell_delp                       pressure difference between 
!                                        level k and the pressure at
!                                        cloud base (lowest level with
!                                        liquid condensate) [ Pa ]
!        cell_liquid_eff_diam_land       droplet effective diameter when
!                                        over a land surface
!                                        [ microns ]
!        cell_liquid_eff_diam_ocean      droplet effective diameter when
!                                        over an ocean surface 
!                                        [ microns ]
!        cell_droplet_number_land
!        cell_droplet_number_ocean
!        i,j,k,n                         do-loop indices
!
!-------------------------------------------------------------------

      ermesg= ' ' ; error = 0

!--------------------------------------------------------------------
!    define the pressure (cell_pbase) and temperature (temp_cell_pbase)
!    at the cell cloud base (lowest level with liquid water present).
!--------------------------------------------------------------------
      cell_pbase = pfull(:,:,1)
      temp_cell_pbase = temp(:,:,1)
      do j=1,jsize
        do i=1,isize
          do k=nlev_lsm,1,-1
            if (cuql(i,j,k) >= 1.0e-11 )  then
              cell_pbase(i,j) = pfull(i,j,k)
              temp_cell_pbase(i,j) = temp(i,j,k)
              exit
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------
!    define the pressure distance between the cell cloud base and each
!    model pressure level.
!---------------------------------------------------------------------
      do k=1,nlev_lsm             
        cell_delp(:,:,k) = cell_pbase(:,:) - pfull(:,:,k)
      end do

!--------------------------------------------------------------------
!    define the pressure distance between the cell cloud base and points
!    delz_land meters and delz_ocean meters above cloud base. between
!    cloud base and this distance above cloud base, cloud drop size 
!    will be computed; for greater distances above cloud base, the 
!    appropriate drop radii existing as parameters in this module will
!    be used.
!--------------------------------------------------------------------
      cell_land_ref_delp =      &
               cell_pbase*(1.0 - EXP( -(Param%delz_land*Param%grav/   &
                                       (Param%rdgas*temp_cell_pbase))))
      cell_ocean_ref_delp =      &
               cell_pbase*(1.0 - EXP( -(Param%delz_ocean*Param%grav/  &
                                       (Param%rdgas*temp_cell_pbase))))

!---------------------------------------------------------------------
!    if in a diagnostics column, output the vertical profiles of cloud 
!    liquid (cuql), cloud base separation (cell_delp) and pressure 
!    (pfull), along with the column values of land and ocean reference
!    pressure deltas, cell cloud base and land fraction.
!---------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          do j=1,jsize
            do i=1,isize
              if (j == Col_diag%j_dc(n) .and. i == Col_diag%i_dc(n)) then
                if (.not. exit_flag(i,j)) then
                  do k=1,nlev_lsm             
                    if (cuql(i,j,k) > 0.0) then
                      write (Col_diag%unit_dc(n), '(a, e22.12)') &
                           ' Don_conv%cuql', cuql(i,j,k)
                      write (Col_diag%unit_dc(n), '(a, e22.12)') & 
                           ' cell_delp', cell_delp(i,j,k)
                      write (Col_diag%unit_dc(n), '(a, e22.12)') & 
                              ' pfull',  pfull(i,j,k)
                    endif
                  end do
                  write (Col_diag%unit_dc(n), '(a, e22.12)') &
                      ' cell_land_ref_delp',  cell_land_ref_delp(i,j)
                  write (Col_diag%unit_dc(n), '(a, e22.12)') &
                      ' cell_ocean_ref_delp', cell_ocean_ref_delp(i,j)
                  write (Col_diag%unit_dc(n), '(a, e22.12)') &
                                  ' land', land(i,j)
                  write (Col_diag%unit_dc(n), '(a, e22.12)') &
                                 ' cell_pbase', cell_pbase(i,j)
                endif
              endif
            end do
          end do
        end do
      endif

!--------------------------------------------------------------------
!    compute the drop diameters for the land and ocean portions of the 
!    grid box, when liquid water is present.
!--------------------------------------------------------------------
      do k=1,nlev_lsm             
        do j=1,jsize
          do i=1,isize
            if (cuql(i,j,k) >= 1.0e-11 .and. &
                cual(i,j,k) > 0.0) then

!---------------------------------------------------------------------
!    if land is present in the box and the box is more than delz_land 
!    meters from cloud base, define the cloud drop diameter as the
!    value of r_conv_land.
!---------------------------------------------------------------------
              if (land(i,j) > 0.0) then
                if (cell_delp(i,j,k) >= cell_land_ref_delp(i,j)) then 
                  cell_liquid_eff_diam_land(i,j,k) =      &
                                                  2.0*Param%r_conv_land
                  cell_droplet_number_land(i,j,k) =  3.0*cuql(i,j,k)/  &
                                     (4.0*Param%pie*Param%dens_h2o*    &
                                         (Param%r_conv_land*1.0e-06)**3)
                  cell_droplet_number_land(i,j,k) =  &
                       cell_droplet_number_land(i,j,k)*pfull(i,j,k)/&
                          (Param%rdgas*temp(i,j,k))

!---------------------------------------------------------------------
!    if the box is less than delz_land meters from cloud base, calculate
!    the cloud drop diameter.
!---------------------------------------------------------------------
                else
                  cell_liquid_eff_diam_land(i,j,k) = 2.0*(1.0e6)*  &
                                 (3.0*(pfull(i,j,k)/    &
                            (Param%rdgas*temp(i,j,k)))*cuql(i,j,k)/   &
                     (4*Param%pie*Param%dens_h2o*Param%n_land))**(1./3.) 
                  cell_droplet_number_land(i,j,k) = Param%n_land/  &
                                (pfull(i,j,k)/(Param%rdgas*temp(i,j,k)))
                endif
              else
                cell_liquid_eff_diam_land(i,j,k) = 0.0
                cell_droplet_number_land(i,j,k) = 0.0
              endif

!---------------------------------------------------------------------
!    if any fraction of the grid box is over the ocean and the grid
!    point is more than delz_ocean above cloud base, define the 
!    effective cloud drop diameter for that portion of the box in terms
!    of r_conv_ocean.
!---------------------------------------------------------------------
              if (land(i,j) < 1.0) then
                if (cell_delp(i,j,k) >= cell_ocean_ref_delp(i,j)) then 
                  cell_liquid_eff_diam_ocean(i,j,k) =     &
                                                 2.0*Param%r_conv_ocean
                  cell_droplet_number_ocean(i,j,k) =  3.0*cuql(i,j,k)/ &
                                   (4.0*Param%pie*Param%dens_h2o*    &
                                       (Param%r_conv_ocean*1.0e-06)**3)
                  cell_droplet_number_ocean(i,j,k) =  &
                        cell_droplet_number_ocean(i,j,k)*pfull(i,j,k)/&
                              (Param%rdgas*temp(i,j,k))
                else
                  cell_liquid_eff_diam_ocean(i,j,k) = 2.0*(1.0e6)*  &
                         (3.0*(pfull(i,j,k)/     &
                               (Param%rdgas*temp(i,j,k)))*cuql(i,j,k)/  &
                    (4*Param%pie*Param%DENS_H2O*Param%n_ocean))**(1./3.)
                  cell_droplet_number_ocean(i,j,k) = Param%n_ocean/  &
                               (pfull(i,j,k)/(Param%rdgas*temp(i,j,k)))
                endif
              else
                cell_liquid_eff_diam_ocean(i,j,k) = 0.0
                cell_droplet_number_ocean(i,j,k) = 0.0
              endif

!---------------------------------------------------------------------
!    define the effective diameter for the grid box as the weighted av-
!    erage of the diameters over the land and ocean portions of the box.
!---------------------------------------------------------------------
              cell_liquid_eff_diam(i,j,k) =                    &
                       land(i,j) *cell_liquid_eff_diam_land(i,j,k)   + &
                (1.0 - land(i,j))*cell_liquid_eff_diam_ocean(i,j,k) 
              cell_droplet_number (i,j,k) =                    &
                       land(i,j) *cell_droplet_number_land(i,j,k)   + &
                (1.0 - land(i,j))*cell_droplet_number_ocean(i,j,k) 

!---------------------------------------------------------------------
!    when there is no liquid in the box, set the effective diameter to
!    10 microns.
!---------------------------------------------------------------------
            else
              cell_liquid_eff_diam(i,j,k) = 10.0
              cell_droplet_number(i,j,k) = 0.
            endif  ! (cuql > 1.0e-11)
          end do
        end do
      end do

!--------------------------------------------------------------------
!    limit the liquid droplet sizes to be between 8.401 and 33.199
!    microns, the limits for which the slingo parameterization for
!    radiative properties of cloud drops is applicable.
!--------------------------------------------------------------------
      if (Nml%use_memphis_size_limits) then
        cell_liquid_eff_diam = MAX(8.401, cell_liquid_eff_diam)
        cell_liquid_eff_diam = MIN(33.199, cell_liquid_eff_diam)
      endif

!--------------------------------------------------------------------
!    if this is a diagnostics column, output the grid box effective
!    size, as well as the over land and over ocean values which were
!    used to define it.
!--------------------------------------------------------------------
      if (Col_diag%in_diagnostics_window) then
        do n=1,Col_diag%ncols_in_window
          do j=1,jsize
            do i=1,isize
              if (j == Col_diag%j_dc(n) .and. i == Col_diag%i_dc(n)) then
                do k=1,nlev_lsm            
                  if (cuql(i,j,k) > 0.0) then
                    write (Col_diag%unit_dc(n), '(a, i5    )') &
                         'k', k
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                            ' cell_liquid_eff_diam',  &
                                 cell_liquid_eff_diam(i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                            ' cell_liquid_eff_diam_land',  &
                                 cell_liquid_eff_diam_land(i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                           ' cell_liquid_eff_diam_ocean',  &
                                 cell_liquid_eff_diam_ocean(i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                            ' cell_droplet_number',  &
                                 cell_droplet_number (i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                            ' cell_droplet_number_land',  &
                                 cell_droplet_number_land(i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                           ' cell_droplet_number_ocean',  &
                                 cell_droplet_number_ocean(i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                           ' cuql     ',  &
                                 cuql     (i,j,k)
                    write (Col_diag%unit_dc(n), '(a, e22.12)') &
                           ' cual     ',  &
                                    cual     (i,j,k)
                  endif
                end do
              endif
            end do
          end do
        end do
      endif

!--------------------------------------------------------------------


end subroutine don_r_define_cell_liquid_size_k


!#####################################################################

subroutine don_r_donner_deep_sum_k  &
         (isize, jsize, nlev_lsm, Param, Nml, Initialized, xliq, &
          xice, cual, ampta1, cuql, cuqi, dgeice,   &
          cell_liquid_eff_diam, Don_rad, ermesg, error)

!------------------------------------------------------------------
!    subroutine donner_deep_sum stores the cloud amount and particle 
!    size fields for the convective cells and the mesoscale anvil clouds
!    associated with the donner_deep convection parameterization. 
!------------------------------------------------------------------

use donner_types_mod, only : donner_param_type, donner_nml_type, &
                             donner_initialized_type, donner_rad_type

implicit none

!------------------------------------------------------------------
integer,                      intent(in)    :: isize, jsize, nlev_lsm
type(donner_param_type),      intent(in)    :: Param
type(donner_nml_type),        intent(in)    :: Nml      
type(donner_initialized_type),                     &
                              intent(in)    :: Initialized
real, dimension(isize,jsize,nlev_lsm),              &
                              intent(in)    :: xliq, xice, cual, cuql, &
                                               cuqi, dgeice,  &
                                               cell_liquid_eff_diam
real, dimension(isize,jsize), intent(in)    :: ampta1
type(donner_rad_type),        intent(inout) :: Don_rad
character(len=*),             intent(out)   :: ermesg
integer,                      intent(out)   :: error
                                                            
!------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie        first and last values of i index values of points 
!                   in this physics window (processor coordinates)
!     js, je        first and last values of j index values of points 
!                   in this physics window (processor coordinates)
!     xliq          mesoscale liquid water content [ kg(h2o) / kg(air) ]
!     xice          mesoscale ice content [ kg(h2o) / kg(air) ]
!     cual          cloud fractional area of convective system
!                   (vclouds plus anvil) [ fraction ]
!     ampta1        mesoscale cloud fraction (anvil) [ fraction ]
!     cuql          cell liquid water content [ kg(h2o) / kg(air) ]
!     cuqi          cell ice content [ kg(h2o) / kg(air) ]
!     cell_liquid_eff_diam
!                   cell liquid droplet effective diameter [ microns ]
!     dgeice        ice crystal generalized effective size [ microns ]
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   local variables:

      real, dimension(isize, jsize, nlev_lsm)  ::   meso_area
      integer :: i, j, k

!--------------------------------------------------------------------
!   local variables:
!
!        meso_area      fractional area of anvil [ fraction ]
!        i,j,k          do-loop indices
!
!--------------------------------------------------------------------

       ermesg= ' ' ; error = 0

!--------------------------------------------------------------------
!    if the cloud data from donner_deep_mod that is to be passed to
!    the radiation package is to be time-averaged, increment the counter
!    of the number of time levels that are included in the sum. if the
!    data is not time averaged, set the counter to 1.
!--------------------------------------------------------------------
      if (Nml%do_average) then
        Don_rad%nsum(:,:) = Don_rad%nsum(:,:) + 1
      else
        Don_rad%nsum(:,:) = 1
      endif
 
!--------------------------------------------------------------------
!    define mesoscale anvil area at each level.
!--------------------------------------------------------------------
      do k=1,nlev_lsm            
        do j=1,jsize
          do i=1,isize
            if (xice(i,j,k) == 0.0) then
              meso_area(i,j,k) = 0.0
            else
             meso_area(i,j,k) = MAX (0.0, ampta1(i,j))
            endif
            if (xliq(i,j,k) /= 0.0) then
              ermesg = ' liquid water present in anvil -- not&
                                           & currently allowed'
              error = 1
              return
            endif
          end do
        end do
      end do

!----------------------------------------------------------------------
!    define the mesoscale anvil properties needed by the radiation 
!    package (cloud fraction, liquid amount, ice amount, liquid droplet 
!    size, ice effective size) for the case where the fields are to 
!    be time-averaged.
!----------------------------------------------------------------------
      if (Nml%do_average) then
        Don_rad%meso_cloud_frac(:,:,:) =         &
                        Don_rad%meso_cloud_frac(:,:,:) + meso_area(:,:,:)
        Don_rad%meso_liquid_amt(:,:,:) = &
                             Don_rad%meso_liquid_amt(:,:,:) + xliq(:,:,:)
        Don_rad%meso_ice_amt(:,:,:)         =   &
                                Don_rad%meso_ice_amt(:,:,:) + xice(:,:,:)
        Don_rad%meso_liquid_size(:,:,:) = &
                                    Don_rad%meso_liquid_size(:,:,:) + &
                                          Nml%meso_liquid_eff_diam_input
        Don_rad%meso_ice_size(:,:,:) = &
                         Don_rad%meso_ice_size(:,:,:) + dgeice(:,:,:)

!----------------------------------------------------------------------
!    define the mesoscale anvil properties needed by the radiation 
!    package (cloud fraction, liquid amount, ice amount, liquid droplet
!    size, ice effective size) for the case where the fields are 
!    not time-averaged.
!----------------------------------------------------------------------
      else
        Don_rad%meso_cloud_frac(:,:,:)  = meso_area(:,:,:)
        Don_rad%meso_liquid_amt(:,:,:)  = xliq(:,:,:)
        Don_rad%meso_ice_amt(:,:,:)     = xice(:,:,:)
        Don_rad%meso_liquid_size(:,:,:) = Nml%meso_liquid_eff_diam_input
        Don_rad%meso_ice_size(:,:,:)    = dgeice(:,:,:)
      endif

!----------------------------------------------------------------------
!    define the cell properties needed by the radiation package (cloud 
!    fraction, liquid amount, ice amount, liquid droplet size, ice eff- 
!    ective size) for the case where the fields are to be time-averaged.
!----------------------------------------------------------------------
      if (Nml%do_average) then
        Don_rad%cell_cloud_frac(:,:,:) = &
                       Don_rad%cell_cloud_frac(:,:,:) + &
                             MAX (0.0, cual(:,:,:) - meso_area(:,:,:) )
        Don_rad%cell_liquid_amt(:,:,:) = &
                   Don_rad%cell_liquid_amt(:,:,:) + cuql(:,:,:)
        Don_rad%cell_ice_amt(:,:,:) = &
                            Don_rad%cell_ice_amt(:,:,:) + cuqi(:,:,:)

!----------------------------------------------------------------------
!     cell liquid size may be either specified via a namelist variable
!     or defined based on the bower parameterization.
!----------------------------------------------------------------------
        if (Initialized%do_input_cell_liquid_size) then
          Don_rad%cell_liquid_size(:,:,:) = &
                              Don_rad%cell_liquid_size(:,:,:) +      &
                                           Nml%cell_liquid_eff_diam_input
        else if (Initialized%do_bower_cell_liquid_size) then
          Don_rad%cell_liquid_size(:,:,:) = &
                 Don_rad%cell_liquid_size(:,:,:) + cell_liquid_eff_diam
        endif

!----------------------------------------------------------------------
!     cell ice size may be either specified using a default value or
!     an input value supplied via the namelist.
!----------------------------------------------------------------------
        if (Initialized%do_default_cell_ice_size) then
          Don_rad%cell_ice_size(:,:,:) = &
                                Don_rad%cell_ice_size(:,:,:) +    &
                                          Param%cell_ice_geneff_diam_def
        else if (Initialized%do_input_cell_ice_size) then
          Don_rad%cell_ice_size(:,:,:) = &
                              Don_rad%cell_ice_size(:,:,:) +      &
                                           Nml%cell_ice_geneff_diam_input
        endif

!----------------------------------------------------------------------
!    define the cell properties needed by the radiation package (cloud 
!    fraction, liquid amount, ice amount, liquid droplet size, ice eff- 
!    ective size) for the case where the fields are not time-averaged.
!----------------------------------------------------------------------
      else
        Don_rad%cell_cloud_frac(:,:,:) = &
                             MAX (0.0, cual(:,:,:) - meso_area(:,:,:) )
        Don_rad%cell_liquid_amt(:,:,:) = cuql(:,:,:)
        Don_rad%cell_ice_amt(:,:,:)    = cuqi(:,:,:)

!----------------------------------------------------------------------
!     cell liquid size may be either specified via a namelist variable
!     or defined based on the bower parameterization.
!----------------------------------------------------------------------
        if (Initialized%do_input_cell_liquid_size) then
          Don_rad%cell_liquid_size(:,:,:) =       &
                                          Nml%cell_liquid_eff_diam_input
        else if (Initialized%do_bower_cell_liquid_size) then
          Don_rad%cell_liquid_size(:,:,:) = cell_liquid_eff_diam
        endif

!----------------------------------------------------------------------
!     cell ice size may be either specified using a default value or
!     an input value supplied via the namelist.
!----------------------------------------------------------------------
        if (Initialized%do_default_cell_ice_size) then
          Don_rad%cell_ice_size(:,:,:)  = Param%cell_ice_geneff_diam_def
        else if (Initialized%do_input_cell_ice_size) then
          Don_rad%cell_ice_size(:,:,:) = Nml%cell_ice_geneff_diam_input
        endif
      endif

!--------------------------------------------------------------------


end subroutine don_r_donner_deep_sum_k



