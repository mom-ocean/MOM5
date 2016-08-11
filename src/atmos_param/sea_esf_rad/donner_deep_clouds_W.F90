!FDOC_TAG_GFDL

                 module donner_deep_clouds_W_mod
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <OVERVIEW>
!          donner deep cloud radiative properties module
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use time_manager_mod,       only: time_type
use mpp_mod,                only: input_nml_file
use fms_mod,                only: open_namelist_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  close_file, FATAL,  &
                                  mpp_pe, mpp_root_pe, &
                                  write_version_number, stdlog
use rad_utilities_mod,      only: microphysics_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          donner deep cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: donner_deep_clouds_W.F90,v 19.0 2012/01/06 20:15:19 fms Exp $'
   character(len=128)  :: tagname =  '$Name: tikal $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          donner_deep_clouds_W_init,   &
          donner_deep_clouds_W_end , donner_deep_clouds_amt

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: using_dge_lw = .true.
logical   :: using_dge_sw = .true.



namelist /donner_deep_clouds_W_nml /     &
       using_dge_sw, using_dge_lw


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


  logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





! <SUBROUTINE NAME="donner_deep_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)
!
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
!  <IN NAME="axes" TYPE="integer">
! 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time

      integer            :: unit, ierr, io, logunit

     if (module_is_initialized) return
!---------------------------------------------------------------------
!-----  read namelist  ------
  
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=donner_deep_clouds_W_nml, iostat=io)
      ierr = check_nml_error(io,"donner_deep_clouds_W_nml")
#else
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif
#endif

      if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         logunit = stdlog()
         write (logunit,nml=donner_deep_clouds_W_nml)
      endif

!---------------------------------------------------------------------

       module_is_initialized = .true.


end subroutine donner_deep_clouds_W_init

! <SUBROUTINE NAME="donner_deep_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_W_end
       
!----------------------------------------------------------------------
!    diag_clouds_W_end is the destructor for diag_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------


end subroutine donner_deep_clouds_W_end


!#################################################################


!---------------------------------------------------------------------

! <SUBROUTINE NAME="donner_deep_clouds_amt">
!  <OVERVIEW>
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call donner_deep_clouds_amt (is, ie, js, je, Cell_microphys,  &
!                Meso_microphys)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine donner_deep_clouds_amt (is, ie, js, je,   &
                   cell_cloud_frac, cell_liquid_amt, cell_liquid_size, &
                   cell_ice_amt, cell_ice_size, &
                   cell_droplet_number, &
                   meso_cloud_frac, meso_liquid_amt, meso_liquid_size, &
                   meso_ice_amt, meso_ice_size, &
                   meso_droplet_number,  nsum_out, &
                   Cell_microphys,  Meso_microphys)

!---------------------------------------------------------------------
!    donner_deep_clouds_amt defines the distribution of cloud water and
!    cloud ice concentration and particle size and total cloud fraction
!    in both the mesoscale and convective cell-scale components of the
!    clouds associated with donner_deep convection. these values will
!    be combined with the large-scale cloud fields to produce the dist-
!    ribution of cloud radiative properties that will be seen by the
!    radiation package.
!----------------------------------------------------------------------

integer,                 intent(in)    :: is,ie,js,je
real, dimension(:,:,:), intent(inout) ::   &
                   cell_cloud_frac, cell_liquid_amt, cell_liquid_size, &
                   cell_ice_amt, cell_ice_size, &
                   cell_droplet_number, &
                   meso_cloud_frac, meso_liquid_amt, meso_liquid_size, &
                   meso_ice_amt, meso_ice_size, &
                   meso_droplet_number
integer, dimension(:,:), intent(inout) ::  nsum_out
type(microphysics_type), intent(inout) :: Cell_microphys, Meso_microphys

!---------------------------------------------------------------------
!     call donner_deep_avg to obtain the specification fields for both
!     the mesoscale and convective cellscale clouds assocated with 
!     donner_deep convection.
!---------------------------------------------------------------------
      call donner_deep_avg (                           &
                      is, ie, js, je,           &
                      cell_cloud_frac,  Cell_microphys%cldamt, &
                      cell_liquid_amt, Cell_microphys%conc_drop, &
                      cell_liquid_size, Cell_microphys%size_drop,&
                      cell_ice_amt, Cell_microphys%conc_ice, &
                      cell_ice_size, Cell_microphys%size_ice, &
                 cell_droplet_number, Cell_microphys%droplet_number, &
                      meso_cloud_frac, Meso_microphys%cldamt,   &
                      meso_liquid_amt, Meso_microphys%conc_drop, &
                      meso_liquid_size, Meso_microphys%size_drop,&
                      meso_ice_amt, Meso_microphys%conc_ice, &
                      meso_ice_size, Meso_microphys%size_ice,   &
                 meso_droplet_number, Meso_microphys%droplet_number, &
                      nsum_out)

!---------------------------------------------------------------------



end subroutine donner_deep_clouds_amt  



!####################################################################

subroutine donner_deep_avg    &
               (is, ie, js, je,           &
                cell_cloud_frac,  cell_cloud_frac_out,   &
                cell_liquid_amt, cell_liquid_amt_out,      &
                cell_liquid_size, cell_liquid_size_out    ,&
                cell_ice_amt, cell_ice_amt_out, &
                cell_ice_size, cell_ice_size_out, &
                cell_droplet_number, cell_droplet_number_out, &
                meso_cloud_frac, meso_cloud_frac_out, &
                meso_liquid_amt, meso_liquid_amt_out, &
                meso_liquid_size, meso_liquid_size_out,&
                meso_ice_amt, meso_ice_amt_out, &
                meso_ice_size, meso_ice_size_out,   &
                meso_droplet_number, meso_droplet_number_out, &
                nsum)

!------------------------------------------------------------------
!    subroutine donner_deep_avg outputs the cloud microphysical quant-
!    ities associated with donner_deep convection for use by the rad-
!    iation package. these fields are the cloud liquid and ice amounts,
!    liquid and ice sizes, and fractional coverage, for both the con-
!    vective cell and mesoscale components of the convective system.
!--------------------------------------------------------------------

integer,                   intent(in)  :: is, ie, js, je
real,    dimension(:,:,:), intent(inout) :: cell_cloud_frac,        &
                                          cell_liquid_amt,   &
                                          cell_liquid_size,  &
                                          cell_ice_amt, &
                                          cell_ice_size,   &
                                       cell_droplet_number, &
                                          meso_cloud_frac,   &
                                          meso_liquid_amt, &
                                          meso_liquid_size, &
                                          meso_ice_amt,  &
                                          meso_ice_size, &
                                        meso_droplet_number
integer, dimension(:,:), intent(inout) :: nsum                        
real,    dimension(:,:,:), intent(out) :: cell_cloud_frac_out,        &
                                          cell_liquid_amt_out,   &
                                          cell_liquid_size_out,  &
                                          cell_ice_amt_out, &
                                          cell_ice_size_out,   &
                                        cell_droplet_number_out, &
                                          meso_cloud_frac_out,   &
                                          meso_liquid_amt_out, &
                                          meso_liquid_size_out, &
                                          meso_ice_amt_out,  &
                                          meso_ice_size_out, &
                                       meso_droplet_number_out

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!
!  intent(out) variables:
!
!     cell_cloud_frac_out  fractional coverage of convective cells in 
!                          grid box [ dimensionless ]
!     cell_liquid_amt_out  liquid water content of convective cells 
!                          [ kg(h2o) / kg(air) ]
!     cell_liquid_size_out assumed effective size of cell liquid drops
!                          [ microns ]
!     cell_ice_amt_out     ice water content of cells 
!                          [ kg(h2o) / kg(air) ]
!     cell_ice_size_out    generalized effective diameter for ice in
!                          convective cells [ microns ]
!     meso_cloud_frac_out  fractional area of mesoscale clouds in grid 
!                          box [ dimensionless ]
!     meso_liquid_amt_out  liquid water content in mesoscale clouds
!                          [ kg(h2o) / kg(air) ]
!     meso_liquid_size_out assumed effective size of mesoscale drops
!                          [ microns ]
!     meso_ice_amt_out     ice water content of mesoscale elements 
!                          [ kg(h2o) / kg(air) ]
!     meso_ice_size_out    generalized ice effective size for anvil ice
!                          [ microns ]
!
!---------------------------------------------------------------------

   
!------------------------------------------------------------------
!   local variables

      real, dimension(size(cell_cloud_frac_out,1), &
                      size(cell_cloud_frac_out,2))   :: inv_nsum
      integer                                        :: num
      integer                                        ::   k
   
!---------------------------------------------------------------------
!   local variables:
!
!       inv_sum     inverse of number of elements in the time averaged
!                   output fields
!       num         number of grid columns which have not been given
!                   values for the output variables by donner_deep_mod
!       i,j,k       do-loop indices
!
!--------------------------------------------------------------------- 

!---------------------------------------------------------------------
!    check to make sure dimensions of arguments match the module
!    variable dimensions.
!---------------------------------------------------------------------
      if (size(cell_cloud_frac_out,3) /= size(cell_cloud_frac,3)) &
        call error_mesg (  'donner_rad_mod',  &
                         'input argument has the wrong size',FATAL)

!---------------------------------------------------------------------
!    check to see that all columns have been given values for the module
!    variables that are going to be averaged. 
!----------------------------------------------------------------------
      num = count(        nsum(:,:) == 0)

!----------------------------------------------------------------------
!    if all points have values of 0, then the scheme has not yet been
!    called. use the initialized values that are present.
!----------------------------------------------------------------------
      if (num == (ie-is+1)*(je-js+1) ) then
        cell_cloud_frac_out             = cell_cloud_frac
        cell_liquid_amt_out             = cell_liquid_amt
        cell_liquid_size_out            = cell_liquid_size
        cell_ice_amt_out                = cell_ice_amt
        cell_ice_size_out               = cell_ice_size
        cell_droplet_number_out         = cell_droplet_number
        meso_cloud_frac_out             = meso_cloud_frac
        meso_liquid_amt_out             = meso_liquid_amt
        meso_liquid_size_out            = meso_liquid_size
        meso_ice_amt_out                = meso_ice_amt
        meso_ice_size_out               = meso_ice_size
        meso_droplet_number_out         = meso_droplet_number
        
!----------------------------------------------------------------------
!    if any columns have not been given values, stop execution with an 
!    error message. 
!----------------------------------------------------------------------
      else if (num > 0) then
        call error_mesg ( 'donner_rad_mod', &
                         'nsum has some zero entries', FATAL)

!----------------------------------------------------------------------
!    if all columns have valid data, produce time averaged values of
!    the desired output fields.
!----------------------------------------------------------------------
      else
        inv_nsum(:,:) = 1.0/float(nsum(  :  , :   ))
        do k=1,size(cell_cloud_frac_out,3)
          cell_cloud_frac_out(:,:,k) =   &
                              cell_cloud_frac(:,:,        k)*inv_nsum
          cell_liquid_amt_out(:,:,k) =   &
                              cell_liquid_amt(:,:,        k)*inv_nsum
          cell_liquid_size_out(:,:,k) =  &
                              cell_liquid_size(:,:,k        )*inv_nsum
          cell_ice_amt_out(:,:,k) =      &
                              cell_ice_amt(:,:,k        )*inv_nsum
          cell_ice_size_out(:,:,k) =     &
                              cell_ice_size(:,:,k        )*inv_nsum
          cell_droplet_number_out(:,:,k) =     &
                              cell_droplet_number(:,:,k      )*inv_nsum
          meso_cloud_frac_out(:,:,k) =   &
                              meso_cloud_frac(:,:,        k)*inv_nsum
          meso_liquid_amt_out(:,:,k) =   &
                              meso_liquid_amt(:,:,        k)*inv_nsum
          meso_liquid_size_out(:,:,k) =  &
                              meso_liquid_size(:,:,        k)*inv_nsum
          meso_ice_amt_out(:,:,k) =      &
                              meso_ice_amt(:,:,        k)*inv_nsum
          meso_ice_size_out(:,:,k) =     & 
                              meso_ice_size(:,:,        k)*inv_nsum
          meso_droplet_number_out(:,:,k) =     & 
                              meso_droplet_number(:,:,     k)*inv_nsum
        end do
     
!---------------------------------------------------------------------
!     prevent the occurrence of cloud area with no condensate and cond-
!     ensate with no area, for both the convective cell clouds and the
!     mesoscale cloud.
!---------------------------------------------------------------------
        where (cell_cloud_frac_out  > 0.0 .and.  &
               cell_liquid_amt_out == 0.0 .and.  &
               cell_ice_amt_out    == 0.0)
          cell_cloud_frac_out  = 0.0
        end where
        where (cell_cloud_frac_out == 0.0 .and.   &
               cell_liquid_amt_out  > 0.0)
          cell_liquid_amt_out  = 0.0
        end where
        where (cell_cloud_frac_out == 0.0 .and.   &
               cell_ice_amt_out     > 0.0)
          cell_ice_amt_out     = 0.0
        end where

        where (meso_cloud_frac_out  > 0.0 .and.  &
               meso_liquid_amt_out == 0.0 .and.  &
               meso_ice_amt_out    == 0.0)
          meso_cloud_frac_out  = 0.0
        end where
        where (meso_cloud_frac_out == 0.0 .and.   &
               meso_liquid_amt_out  > 0.0)
          meso_liquid_amt_out  = 0.0
        end where
        where (meso_liquid_amt_out == 0.0)
          meso_liquid_size_out = 0.0
        end where
        where (meso_cloud_frac_out == 0.0 .and.   &
               meso_ice_amt_out     > 0.0)
          meso_ice_amt_out     = 0.0
        end where
      endif
 
!----------------------------------------------------------------------
!    reset the variables just processed so that new sums may be begun 
!    when donner_deep is called again.
!----------------------------------------------------------------------
      cell_cloud_frac                 = 0.0
      cell_liquid_amt                 = 0.0
      cell_liquid_size                = 0.0
      cell_ice_amt                    = 0.0
      cell_ice_size                   = 0.0
      cell_droplet_number             = 0.0
      meso_cloud_frac                 = 0.0
      meso_liquid_amt                 = 0.0
      meso_liquid_size                = 0.0
      meso_ice_amt                    = 0.0
      meso_ice_size                   = 0.0
      meso_droplet_number             = 0.0
      nsum                            = 0
       
!----------------------------------------------------------------------



end subroutine donner_deep_avg


!#####################################################################

       end module donner_deep_clouds_W_mod



