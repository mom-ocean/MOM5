module nf_utils_mod

use nfu_mod ! netcdf utilities
use nfc_mod ! netcdf utilities for compressed files

implicit none
private

! ==== public interfaces =====================================================
! export stuff from nfu_mod
public :: nfu_inq_dim, nfu_inq_var, nfu_inq_att
public :: nfu_def_dim, nfu_def_var
public :: nfu_put_att
public :: nfu_get_dim, nfu_get_dim_bounds
public :: nfu_put_var, nfu_put_rec
public :: nfu_get_var, nfu_get_rec
public :: nfu_get_valid_range, nfu_is_valid, nfu_validtype, nfu_validtype2ascii
! export stuff from nfc_mod
public :: nfu_inq_compressed_dim, nfu_inq_compressed_var
public :: nfu_get_compressed_var
public :: nfu_put_compressed_var
public :: nfu_get_compressed_rec
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: nf_utils.F90,v 17.0 2009/07/21 03:02:54 fms Exp $', &
     tagname = '$Name: tikal $'

end module nf_utils_mod
