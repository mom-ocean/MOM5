module cpl_netcdf_setup_mod
!
!=========================================================================

implicit none

include 'netcdf.inc'
integer :: pLonDimId, pLatDimId, timeDimId, pDepDimId

contains 

subroutine ncheck(status, error_str) 

    implicit none

    integer, intent(in) :: status
    character(len=*), intent(in), optional :: error_str

    if (status /= nf_noerr) then
      write(*,'(/a)')   'error - from NetCDF library'
      if (present(error_str)) then
          write(*,'(a)')   error_str
      endif
      write(*,'(a/)')   trim(nf_strerror(status))
      stop
    end if

end subroutine ncheck

subroutine create_ncfile(ncfile, ncid, ii, jj, kk, ll, ilout)
    ! Create 2, 3,or 4D ncfile, depending on optional args (kk,ll)

    implicit none

    integer, intent(in) :: ii,jj     	!x, y dimension size
    integer, optional :: kk, ll !z, t dimension size
    integer, optional :: ilout  !format io file id
    character(len=*), intent(in) :: ncfile  
    integer, intent(out) :: ncid

#if defined(DEBUG)
    if (present(ilout)) then 
      write(ilout,*) 'creating a new netcdf file: ',ncfile
      write(ilout,*) '    with nx, ny = ', ii, jj
    endif
#endif

    !create a new NetCDF and define the grid:
    call ncheck(nf_create(trim(ncfile),nf_write,ncid), "Creating: "//trim(ncfile))

    !define the dimensions
    if (present(ll)) then
        call ncheck(nf_def_dim(ncid,"time", nf_unlimited,  timeDimId), 'Defining time in '//trim(ncfile)//' in create_ncfile()')
    endif
    if (present(kk)) then
        call ncheck(nf_def_dim(ncid,"nz", kk,  pDepDimId), 'Defining nz in '//trim(ncfile)//' in create_ncfile()')
    endif
    call ncheck(nf_def_dim(ncid, "ny", jj,  pLatDimId), 'Defining ny in '//trim(ncfile)//' in create_ncfile()')
    call ncheck(nf_def_dim(ncid, "nx", ii,  pLonDimId), 'Defining nx in '//trim(ncfile)//' in create_ncfile()')

    call ncheck(nf_enddef(ncid))

    return
end subroutine create_ncfile

subroutine write_nc_1Dtime(vin, nt, vname, ncid)

    implicit none

    integer, intent(in) :: ncid,nt
    integer :: varid, ncstatus
    integer, dimension(1:6) :: adate
    real, intent(in) :: vin
    ! NOTE here real is default real*8 (which is actually the same as dbl_kind!)
    ! somehow the netcdf lib used here takes 'real' as real*4. therefore we need: 
    real*4 :: vtmp
    character(len=*), intent(in) :: vname
    character*80 ctimeatt

    vtmp = real(vin)
    ncstatus=nf_inq_varid(ncid,vname,varid)

    if (ncstatus/=nf_noerr) then
      call ncheck(nf_redef(ncid))
      call ncheck(nf_def_var(ncid,trim(vname),nf_real, 1, timeDimId, varid), 'Defining: '//trim(vname)//' in write_nc_1Dtime()')
      write(ctimeatt, &
          '("seconds since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
      call ncheck(nf_put_att_text(ncid,varid,"units",len_trim(ctimeatt),trim(ctimeatt)), 'Adding attribute units to '//trim(vname)//' in write_nc_1Dtime()')
      call ncheck(nf_enddef(ncid))
    end if

    !write values into the 1D array
    call ncheck(nf_put_vara_real(ncid,varid,nt,1,vtmp), 'Writing values to '//trim(vname)//' in write_nc_1Dtime()')

    return
end subroutine write_nc_1Dtime

subroutine write_nc2D(ncid, vname, vin, prcn, nx, ny, istep, ilout)
    !
    !to output a 2D array into a 3D field (with time dimension) 
    !with either single or double precisioin depending on argumnet 'prcn'!
    !

    implicit none

    integer, intent(in) :: ncid
    integer, intent(in) :: prcn	!precision choice (1/2: signle/double)
    character(len=*), intent(in) :: vname
    integer, intent(in) :: nx, ny
    integer, intent(in) :: istep	!position in the time dim (No of record) 
    integer, optional :: ilout
    real, dimension(nx,ny), intent(in) :: vin

    integer :: varid, ncstatus 
    real*4, dimension(nx,ny) :: vtmp   !single precision

#if defined(DEBUG)
    if (present(ilout)) write(ilout,*) 'write_nc2D: handling var *** ',vname, ' rec: ', istep
#endif

    ncstatus=nf_inq_varid(ncid,vname,varid)
    if (ncstatus/=nf_noerr) then
      call ncheck(nf_redef(ncid))
      if (prcn == 1) then
        call ncheck(nf_def_var(ncid,trim(vname),nf_real, 3, &
                (/pLonDimId, pLatDimId, timeDimId/),varid), 'Defining real '//trim(vname)//' in write_nc2D')
      else
        call ncheck(nf_def_var(ncid,trim(vname),nf_double, 3, &
                (/pLonDimId, pLatDimId, timeDimId/),varid), 'Defining double '//trim(vname)//' in write_nc2D')
      endif
      call ncheck(nf_enddef(ncid))
#if defined(DEBUG)
      if (present(ilout)) write(ilout,*) 'write_nc2D: defined new var ***', vname 
#endif
    else
#if defined(DEBUG)
      if (present(ilout)) write(ilout,*) 'write_nc2D: found   old var ***', vname
#endif
    end if

    select case(prcn)
      case (1)
        vtmp = real(vin) !dbl precision to single precision
        call ncheck(nf_put_vara_real(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vtmp), 'Writing real values to'//trim(vname)//' in write_nc2D')
      case default    !case (2)
        call ncheck(nf_put_vara_double(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vin), 'Writing double values to'//trim(vname)//' in write_nc2D')
    end select

    return
end subroutine write_nc2D

end module cpl_netcdf_setup_mod
