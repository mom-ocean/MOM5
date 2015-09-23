module fv_restart_mod

  use fms_mod,    only: file_exist, error_mesg,  lowercase, FATAL, &
       NOTE, read_data, uppercase, stdout, mpp_pe, mpp_root_pe, mpp_chksum
  use fms_io_mod,   only: open_restart_file, close_file, write_data, get_instance_filename
  use mpp_io_mod, only: mpp_open, mpp_close, axistype, fieldtype,  &
       mpp_read_meta, mpp_get_info, mpp_get_fields, mpp_read, &
       MPP_NETCDF, MPP_SINGLE, MPP_MULTI, MPP_SEQUENTIAL, MPP_RDONLY, &
       MPP_WRONLY, MPP_NATIVE, MPP_OVERWR, mpp_write_meta, mpp_write, &
       atttype, mpp_get_atts, mpp_get_att_name, mpp_get_att_real_scalar, &
       mpp_get_field_name
  use fv_pack,    only: ak, bk, ks, ptop, beglat, endlat,       &
       master, nlon, mlat, nt_phys, nqrst, ng_d, ng_s, &
       ighost, mountain, coslon, sinlon, lon, lat, lonb, latb,     &
       adjust_dry_mass, use_set_eta, p_var, d2a3d,          &
       age_tracer, age_time, get_eta_level, &
       do_fms_tracer_manager, fms_tracers_file, fv_domain, &
       restart_format, use_tendency, cold_start, gid
#ifdef MARS_GCM
  use fv_pack,    only: p_ref
#endif MARS_GCM

  use tracer_manager_mod, only : tr_get_tracer_names=>get_tracer_names, &
       get_tracer_names, get_number_tracers, &
       set_tracer_profile, &
       get_tracer_index, NO_TRACER
  use field_manager_mod, only  : MODEL_ATMOS      
  use         platform_mod, only: I8_KIND

!  use shr_kind_mod, only : r8 => shr_kind_r8
#ifdef SPMD
  use mod_comm, only: mp_bcst_int, mp_bcst_r2d, mp_scatter3d, mp_bcst_real, &
       mp_send4d_ns, mp_recv4d_ns, mp_gather3d
#endif

  use fv_arrays_mod, only: fv_array_check, fv_stack_push, fv_array_sync, &
       fv_print_chksums, ksp, kep
  use mpp_mod, only: mpp_sync, mpp_root_pe, mpp_broadcast
  use mpp_domains_mod, only: mpp_get_data_domain, mpp_get_compute_domain
  use constants_mod,    only: grav, kappa, rdgas

  use pmaxmin_mod, only: pmaxming, pmaxmin

  implicit none
  private

  public :: fv_restart, read_fv_rst, write_fv_rst, add_tracers

contains

  subroutine read_fv_rst(im, jm, km, nq, dflnm, days, seconds, oform)

#include "fv_arrays.h"
    character*(*), intent(in)       :: dflnm                ! dynamics filename
    integer, intent(in) :: nq
    integer, intent(out):: im, jm, km
    integer, intent(out):: days
    integer, intent(out):: seconds
    character(len=*), intent(in) :: oform

! Local variables
    integer iuic                       ! Unit number
    integer i, j, k, n, jn, it, ntracers
    logical is_tracer_in_restart_file
    character(len=128)    :: tname
    real zsurf(nlon, beglat: endlat)
    real ginv, qmax, qmin
    real, allocatable :: w2d(:,:), w3d(:,:,:)

    type(fieldtype), allocatable :: tracer_fields(:), fields(:)
    type(atttype), allocatable :: attributes(:)
    integer :: ndim,natt,nvar,ntime
    real p0

    character(len=*), parameter :: fms_tracers_restart_file = &
         'INPUT/'//fms_tracers_file

    character(len=128) :: directory, restart_name, tr_name, dflnm_mod
    logical :: use_native_format
    integer :: ids, ide, jds, jde, ics,ice,jcs,jce

#ifdef MARS_GCM
!      ------Dont require specific water vapor fields in netcdf version 
!         of tracer restart file:   Hence this code is not necessary 
#else
    integer :: nsphum, nliq_wat, nice_wat, ncld_amt
    character(len=20) :: fv_tr_name(4)
    data fv_tr_name / "sphum","liq_wat","ice_wat","cld_amt" /
#endif


#ifdef use_shared_pointers
    pointer( p_zsurf, zsurf )
#include "fv_point.inc"
    call fv_stack_push( p_zsurf, nlon*(endlat-beglat+1) )
#endif

    ginv = 1./ GRAV

    !Modify the filename if it has an ensemble instance appendix
    !The following call does not modify the filename unless there 
    !is an ensemble instance appendix and fms_io knows about it.
    call get_instance_filename(dflnm, dflnm_mod)

! output format
    if (trim(uppercase(oform)) == 'NATIVE') then
        use_native_format = .true.
        restart_name = trim(dflnm_mod)
    else if (trim(uppercase(oform)) == 'NETCDF') then
        use_native_format = .false.
        if (trim(dflnm_mod(LEN_TRIM(dflnm_mod)-2:LEN_TRIM(dflnm_mod))) /= ".nc" ) &
             restart_name = trim(dflnm_mod)//".nc"
        if (file_exist(dflnm_mod) .and. .not. file_exist(trim(dflnm_mod)//".nc")) then
            call error_mesg( 'read_fv_rst', &
                 'NetCDF file does not exist. Switching to Native format.', &
                 NOTE)
            restart_name = trim(dflnm_mod)
            use_native_format = .true.
        endif

    else
        call error_mesg ('read_fv_rst', &
             'invalid value for argument oform='//trim(oform), FATAL)
    endif
    if (.not. file_exist(dflnm_mod) .and. .not. file_exist(trim(dflnm_mod)//".nc")) &
         call error_mesg ('read_fv_rst', &
         'Neither a Native or NetCDF restart file is present.', FATAL)

    if( master ) then
        if (use_native_format) then
            if (file_exist(dflnm_mod)) then
!!$         call mpp_open (iuic, dflnm_mod, action=MPP_RDONLY, form=MPP_NATIVE, access=MPP_SEQUENTIAL)
                iuic = open_restart_file(dflnm_mod, 'read')
                read(iuic) im, jm, km, nqrst
            endif
        else
            if (file_exist(trim(dflnm_mod)//'.nc')) then
                call mpp_open(iuic, dflnm_mod, &
                     action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_SINGLE )
                call mpp_read_meta(iuic)
                call mpp_get_info(iuic,ndim,nvar,natt,ntime)

                allocate(attributes(natt))
                call mpp_get_atts(iuic,attributes)

                do n= 1,natt
!            select case (attributes(n)%name)
                   select case ( mpp_get_att_name(attributes(n)))
                   case ('nlon')
                       im= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('mlat')
                       jm= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('nlev')
                       km= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('nt_phys')
                       nqrst= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('days')
                       days= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('seconds')
                       seconds= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('ks')
                       ks = int(mpp_get_att_real_scalar(attributes(n)))
                   end select
                enddo
                deallocate(attributes)


                allocate(fields(nvar))
                call mpp_get_fields(iuic, fields)
                do n=1,nvar
!            select case (fields(n)%name)
                   select case ( mpp_get_field_name(fields(n)))
                   case ('ak')
                       call mpp_read(iuic,fields(n),ak)
                   case ('bk')
                       call mpp_read(iuic,fields(n),bk)
                   end select
                enddo
                call mpp_close(iuic)
                deallocate(fields)
            endif
        endif
        if ( im /= nlon .or. jm /= mlat .or. km /= nlev ) then
            write(6,*) 'Dataset resolution:',im, jm, km, nqrst
            call error_mesg('read_fv_rst:','resolution inconsistent', FATAL)
        endif

        if ( nqrst > nt_phys .or. nqrst > nq) then
            write(6,*) ' '
            write(6,*) '---- Warning -------'
            write(6,*) 'Total number of tracers in restart file is greater than either'
            write(6,*) 'NT_PHYS or NQ.  The number NT_PHYS from the namelist is used.'
            write(6,*) 'NQRST=', nqrst, 'NT_PHYS=', nt_phys, ' NQ=', nq
            nqrst = nt_phys
        endif
        if (use_native_format) then
            read(iuic) days, seconds
            read(iuic) ak, bk, ks
        endif
        write(6,*) 'Restart days=', days,' sec=',seconds
        write(6,*) '----------------------------------'
        write(6,*) 'Model Top Pressure (pa)= ', ak(1)
        write(6,*) '----------------------------------'
        write(6,*) 'Checking Initial condition ...'
    endif

#ifdef SPMD           
!    call mp_bcst_int(days)
!    call mp_bcst_int(seconds)
!    call mp_bcst_int(im)
!    call mp_bcst_int(jm)
!    call mp_bcst_int(km)
!    call mp_bcst_int(nqrst)
!    call mp_bcst_int(ks)
!
!    call mp_bcst_r2d(km+1, 1, 1, 1, ak, 0)
!    call mp_bcst_r2d(km+1, 1, 1, 1, bk, 0)
    call mpp_broadcast( days, mpp_root_pe() )
    call mpp_broadcast( seconds, mpp_root_pe() )
    call mpp_broadcast( im, mpp_root_pe() )
    call mpp_broadcast( jm, mpp_root_pe() )
    call mpp_broadcast( km, mpp_root_pe() )
    call mpp_broadcast( nqrst, mpp_root_pe() )
    call mpp_broadcast( ks, mpp_root_pe() )

    call mpp_broadcast( ak, km+1, mpp_root_pe() )
    call mpp_broadcast( bk, km+1, mpp_root_pe() )
#endif

    ptop = ak(1)

    if (use_native_format) then   
!---------------------------
! read surface geopotential
!---------------------------
        allocate ( w2d(im,jm) )
#ifdef SPMD
        if ( master ) read(iuic) w2d
        call mp_scatter3d(w2d, phis, im,  jm, 1, beglat, endlat,  &
             1, 1, ighost, ighost, 0)
        call fv_array_sync
#else
        read(iuic) w2d
        do j=beglat,endlat
           do i=1,im
              phis(i,j) = w2d(i,j)
           enddo
        enddo
#endif
!       call pmaxming('PHIS', phis, im, jm, 1, beglat, endlat,   &
!                     ighost, ighost, ginv)

        allocate ( w3d(im, jm, km) )

!---------------------------
! Read u-wind
!---------------------------

#ifdef SPMD
        if ( master ) read(iuic) w3d
        call mp_scatter3d(w3d, u, im,  jm, km, beglat, endlat,  &
             1, km, ng_d, ng_s, 0)
        call fv_array_sync
#else
        read(iuic) u
#endif

!---------------------------
! Read v-wind
!---------------------------

#ifdef SPMD
        if ( master ) read(iuic) w3d
        call mp_scatter3d(w3d, v, im,  jm, km, beglat, endlat,  &
             1, km, ng_d, ng_d, 0)
        call fv_array_sync
#else
        read(iuic) v
#endif

!-----------------------------------
! Read virtual potential temperature
!-----------------------------------

#ifdef SPMD
        if ( master ) read(iuic) w3d
        call mp_scatter3d(w3d, pt, im,  jm, km, beglat, endlat,  &
             1, km, ng_d, ng_d, 0)
        call fv_array_sync
#else
        read(iuic) pt
#endif

!-----------------------------------
! Read pressure thickness (pascal)
!-----------------------------------

#ifdef SPMD
        if ( master ) read(iuic) w3d
        call mp_scatter3d(w3d, delp, im,  jm, km, beglat, endlat,  &
             1, km, 0, 0, 0)
        call fv_array_sync
#else
        read(iuic) delp
#endif

    else  
!this code may not work? all these fields do not have the same data domain
        call mpp_get_data_domain( fv_domain, ids, ide, jds, jde )
        call mpp_get_compute_domain( fv_domain, ics, ice, jcs, jce )
        call read_data ( restart_name, 'Surface_geopotential', phis(ids:ide,jcs:jce), domain = fv_domain)
        call read_data ( restart_name, 'U', u(ids:ide,jds:jde,:) , domain = fv_domain)
        call read_data ( restart_name, 'V', v(ids:ide,jds:jde,:), domain = fv_domain)
        call read_data ( restart_name, 'T', pt(ids:ide,jds:jde,:), domain = fv_domain)
        call read_data ( restart_name, 'DELP', delp(ids:ide,jcs:jce,:), domain = fv_domain) !RASF Bug fix for indices
    endif ! not use_native_format

#ifdef SHAVE_P
    call set_eta(km, ks, ptop, ak, bk)
    do j=beglat,endlat
       do i=1,im
          delp(i,j,1) = ak(2) - ak(1)
       enddo
    enddo
#endif

!------------------------------------
! Read tracers one at a time
!------------------------------------

    if (use_native_format) then 
        do n=1, nqrst

#ifdef SPMD
           if ( master ) read(iuic) w3d
           call mp_scatter3d(w3d, q(1,beglat-ng_d,1,n), im,  jm, km, beglat, endlat,  &
                1, km, ng_d, ng_d , 0)
           call fv_array_sync
#else
           read(iuic) w3d
           do k=1,nlev
              do j=beglat,endlat
                 do i=1,im
                    q(i,j,k,n) = w3d(i,j,k)
                 enddo
              enddo
           enddo
#endif
        enddo
    else        ! ----------read "dynamics" tracers  from tracer restart file 
#ifdef MARS_GCM
!rjw      ------No longer require specific water vapor fields in netcdf version of tracer restart file
!rjw          It now makes more sense to read all the tracers from the tracer file
!rjw           in a single loop, that appears below in the processing of field manager tracers
#else
        do n =1,nqrst
           call read_data(fms_tracers_restart_file,fv_tr_name(n),q(ids:ide,jds:jde,:,n) , domain = fv_domain)
        enddo
#endif MARS_GCM 
    endif ! --------------- Finished reading dynamcics tracers ------------


    if( age_tracer ) call error_mesg('read_fv_rst', 'age_tracer == .T. not yet working', FATAL)

    if (use_native_format) then 
        if ( age_tracer .and. nqrst == nq ) then  ! this check will not work properly (pletzer)
!           q = 0. ; age_time = 0.
            if(master) read(iuic) age_time
#ifdef SPMD
!            call mp_bcst_real(age_time)
            call mpp_broadcast( age_time, mpp_root_pe() )
#endif
        endif

        if ( master ) call close_file(iuic)
    endif

!-----------------------------------------
! Read additional tracers from netCDF file
!-----------------------------------------

#ifdef MARS_GCM
    if (do_fms_tracer_manager ) then !{
#else
    if (do_fms_tracer_manager .and. nq > nt_phys) then !{
#endif MARS_GCM

       !Modify the filename if it has an ensemble instance appendix
       call get_instance_filename(fms_tracers_restart_file,restart_name )
       restart_name=trim(restart_name)
        if( .not. file_exist(restart_name//'.nc') ) then !{

            if(master) then
                write(*,*)'---- Warning -------'
                write(*,*)'do_fms_tracer_manager=',do_fms_tracer_manager, &
                     ' no of tracers nq=', nq
                write(*,*)'but netcdf restart file ', restart_name//'.nc', ' does not exist'
                write(*,*)'--------------------'               
            endif

        else !}{

            if(master) then
                write(*,*)'Additional tracers present in field_table. '
                write(*,*)'These will be read from netcdf file ',restart_name
            endif

            call mpp_open(iuic, restart_name, &
                 action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_SINGLE )

            call mpp_read_meta(iuic)
            call mpp_get_info(iuic,ndim,nvar,natt,ntime)
            if ( nt_phys+nvar /= nq ) then
                if (master) then 
                    write(6,*) ' '
                    write(6,*) '---- Warning -------'
                    write(6,*) 'Total number of tracers nq from netcdf and binary restart files does not match'
                    write(6,*) 'nt_phys=', nt_phys,' nvar=', nvar,' nq=', nq, ' nt_phys+nvar /= nq'
                    write(*,*)'---------------------'        
                endif
            endif
            allocate(tracer_fields(nvar))
!if(master) 
            call mpp_get_fields(iuic, tracer_fields)

! read only the tracers that are registered in the tracer manager
            call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
            n = 0 
            call mpp_get_data_domain( fv_domain, ids, ide, jds, jde )
            do it = 1, ntracers !{
! loop over tracers
               call get_tracer_names(MODEL_ATMOS, it, tname)
#ifdef MARS_GCM

#else
               if (use_native_format) then
! exclude tracers read from native format file 
                   if(tname == "sphum"  ) cycle
                   if(tname == "liq_wat") cycle
                   if(tname == "ice_wat") cycle
                   if(tname == "cld_amt") cycle
               endif
#endif MARS_GCM
! this tracer will be added
               n = n + 1
               is_tracer_in_restart_file = .FALSE.
               do jn = 1, nvar !{
! loop over netcdf file tracers                  
!                  if(  lowercase( trim(tracer_fields(jn)%name) ) == lowercase( trim(tname) )  ) then !{
                  if(  lowercase( trim(mpp_get_field_name(tracer_fields(jn))) ) == lowercase( trim(tname) )  ) then !{
                      if(master) print *,'tracer ', trim(tname), ' found!'
! ok found the tracer in the restart file
                      is_tracer_in_restart_file = .TRUE.
! read data and store in array w3d
!                     call mpp_read(iuic, tracer_fields(jn), w3d)
                      call read_data( restart_name, tname, &
                           q(ids:ide,jds:jde,:,it), domain = fv_domain)
                  endif !}
                  if(is_tracer_in_restart_file) exit
               enddo !}
               if(.NOT. is_tracer_in_restart_file) then !{
! missing tracer from restart file
                   if(master) then
                       write(6,*) '---- Warning -------'
                       write(6,*) 'Tracer ', tname, ' was defined in field table but could not be found in restart file.'
                       write(*,*)'---------------------'        
                   endif
! set default profile
                   call set_tracer_profile (MODEL_ATMOS, it, &
                        q(ids:ide,jds:jde,:,it))
               endif !}
            enddo !}


            call mpp_close(iuic)
            deallocate(tracer_fields)

        endif !}

    endif !}

    if (use_native_format) deallocate ( w3d )
!..............................................................................................

    do j=beglat,endlat
       do i=1,im
          zsurf(i,j) = phis(i,j)
       enddo
    enddo

    call pmaxmin('ZS', zsurf, qmin, qmax, im*(endlat-beglat+1), 1, ginv)
    if ( qmax > 1. ) then
        mountain = .true.
    else
        mountain = .false.
    endif

    call pmaxming('U ',       u, im, jm, km, beglat, endlat, ng_d, ng_s, 1.)
    call pmaxming('V ',       v, im, jm, km, beglat, endlat, ng_d, ng_d, 1.)
!    call pmaxming('DELP ', delp, im, jm, km, beglat, endlat, 0,       0, 1.)
    call pmaxming('T ',      pt, im, jm, km, beglat, endlat, ng_d, ng_d, 1.)
    do n=1, nqrst
       call pmaxming('Q ', q(1,beglat-ng_d,1,n), im, jm, km, beglat, endlat, ng_d, ng_d, 1.)
    enddo

!-----------------------
! Add additional tracers:
!-----------------------

#ifdef SET_Q_TOP
! Set water vapor at top layer to climatology value (for stratosphere)
    do j=beglat, endlat
       do i=1,im
          q(i,j,1,1) = 3.E-6
       enddo
    enddo
#endif

    if( use_set_eta ) call set_eta(km, ks, ptop, ak,  bk)

    call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,           &
         pe, peln,  pk,  pkz,  kappa, q, ng_d, ncnst, adjust_dry_mass )

!---------------------------------------------------
! Read surface pressure & winds for coupling purpose
!---------------------------------------------------

    !Modify the filename if it has an ensemble instance appendix
    call get_instance_filename('INPUT/fv_srf_wnd.res',restart_name )
    restart_name=trim(restart_name)

    if (file_exist(restart_name) .or. file_exist(trim(restart_name)//'.nc')) then
        if (use_native_format) then 
            if( master ) then
!!$         call mpp_open (iuic, restart_name, action=MPP_RDONLY,      &
!!$                        form=MPP_NATIVE, access=MPP_SEQUENTIAL)
                iuic = open_restart_file(restart_name, 'read')
                read(iuic) im, jm
                if ( im /= nlon .or. jm /= mlat ) then
                    write(6,*) 'Dataset resolution:',im, jm
                    call error_mesg('read_fv_rst:','resolution inconsistent in fv_srf_wnd', FATAL)
                endif
            endif

! Read ps_bp (surface pressure before physics update):
            if(master) read(iuic) w2d
#ifdef SPMD
            call mp_scatter3d(w2d, ps_bp, im,  jm, 1, beglat, endlat,  &
                 1, 1, 0, 0, 0)
            call fv_array_sync
#endif
! Read u-comp
            if(master) read(iuic) w2d
#ifdef SPMD
            call mp_scatter3d(w2d, u_srf, im,  jm, 1, beglat, endlat,  &
                 1, 1, 0, 0, 0)
            call fv_array_sync
#endif
! Read v-comp
            if(master)  read(iuic) w2d
#ifdef SPMD
            call mp_scatter3d(w2d, v_srf, im,  jm, 1, beglat, endlat,  &
                 1, 1, 0, 0, 0)
            call fv_array_sync
#endif
            if(master)  call close_file(iuic)


        else ! not native format

            call mpp_open(iuic, restart_name, &
                 action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_SINGLE )
            if (master ) then
                call mpp_read_meta(iuic)
                call mpp_get_info(iuic,ndim,nvar,natt,ntime)

                allocate(attributes(natt))
                call mpp_get_atts(iuic,attributes)

                do n= 1,natt
!            select case (attributes(n)%name)
                   select case ( mpp_get_att_name(attributes(n)))
                   case ('nlon')
                       im= int(mpp_get_att_real_scalar(attributes(n)))
                   case ('mlat')
                       jm= int(mpp_get_att_real_scalar(attributes(n)))
                   end select
                enddo
                deallocate(attributes)
            endif

            call mpp_close(iuic)

            if ( im /= nlon .or. jm /= mlat ) then
                write(6,*) 'Dataset resolution:',im, jm
                call error_mesg('read_fv_rst:','resolution inconsistent in fv_srf_wnd', FATAL)
            endif
            call read_data ( restart_name, 'ps_bp', ps_bp(ids:ide,jcs:jce), domain = fv_domain)
            call read_data ( restart_name, 'U_srf', u_srf(ids:ide,jcs:jce), domain = fv_domain)
            call read_data ( restart_name, 'V_srf', v_srf(ids:ide,jcs:jce), domain = fv_domain)
        endif

    else

        if(master) then
            write(*,*) ' '
            write(*,*) '--------------------------------------------'
            write(*,*) 'File not found: ',restart_name
            write(*,*) 'Surface winds will be generated from D grid'
            write(*,*) '--------------------------------------------'
        endif


        do j=beglat,endlat
           do i=1,im
              ps_bp(i,j) = ps(i,j)
           enddo
        enddo

        call d2a3d(u(1,beglat-ng_d,km), v(1,beglat-ng_d,km), u_srf, v_srf, &
             nlon, mlat, 1, beglat, endlat, ng_d, ng_s, coslon, sinlon)
    endif

    call pmaxming('u_srf ', u_srf, im, jm, 1, beglat, endlat, 0,  0, 1.)
    call pmaxming('v_srf ', v_srf, im, jm, 1, beglat, endlat, 0,  0, 1.)

    if (use_native_format) deallocate ( w2d )


! check sum input data
    call print_check_sum ('Check sums for FV input data:')
    if( master ) then 
        call error_mesg ('read_fv_rst', 'FV netcdf restart file read.', NOTE)
    endif


#ifdef PLET
! Do Pure Lagrangian Error Test
    if ( nq >= 7 ) then
        p0 = (1.E5) ** kappa
!$omp  parallel do private(i,j,k)
        do k=ksp,kep
           do j=beglat, endlat
              do i=1,im
                 q(i,j,k,nq-2) = pt(i,j,k)                   
                 q(i,j,k,nq-1) = pt(i,j,k)/pkz(i,j,k) * p0
                 q(i,j,k,nq) = pkz(i,j,k) * p0
              enddo
           enddo
        enddo
        call fv_array_sync()
    endif
#endif

  end subroutine read_fv_rst


  subroutine add_tracers( q, im, jm, km, nq, nqrst, beglat, endlat, ng,  &
       age_tracer, age_time  )

    integer, intent(in):: im,jm,km,nq
    integer, intent(in):: nqrst
    integer, intent(in):: beglat, endlat
    integer, intent(in):: ng
    logical, intent(in):: age_tracer

    real, intent(inout):: age_time
    real, intent(inout):: q(im, beglat-ng:endlat+ng, km, nq)
! Local:
    integer i,j,k,iq

    call fv_array_check( LOC(q) )
    do iq=nqrst+1, nq
!$omp  parallel do private(i,j,k)
       do k=ksp,kep
          do j=beglat, endlat
             do i=1,im
                q(i,j,k,iq) = 0.
             enddo
          enddo
       enddo
    enddo
    call fv_array_sync()

    if ( age_tracer ) age_time = 0.

  end subroutine add_tracers

  subroutine fv_restart(days, seconds)

#include "fv_arrays.h"

    integer, intent(out):: days
    integer, intent(out):: seconds

    real, allocatable :: g_phis(:,:)
    real, allocatable :: qtmp(:,:,:)
    integer :: iuhs          ! Unit to read in terrain data for cold start
    integer :: i,j,k
    integer :: im, jm, km
    character(len=128) :: restart_name

#ifdef MARS_GCM 
    character (len=128) :: filename, fieldname
#endif MARS_GCM

#include "fv_point.inc"

    call fv_print_chksums( 'Entering  fv_restart' )
    if (cold_start) then

#ifdef SW_DYN

        if(nlev /= 1)  call error_mesg('FV_init:','Set nlev=1 for Shallow Water Dynamics', FATAL)

!-----------------------
! Shallow water dynamics
!-----------------------
        ptop  = 0.
        ak(1) = 0.
        ak(2) = 0.
        bk(1) = 0.
        bk(2) = 1.
        KS    = 0
#ifndef USE_LIMA
        n_spong = 0
#endif
        call init_sw_ic(icase, nlon, mlat, nlev, ncnst, beglat, endlat, ng_d, ng_s,  &
             phis, u, v, pt, delp, q, grav, radius, omega, f_d, ighost )
#else

! Start from a "dry & cold" isothermal atmosphere at rest

#ifdef MARS_GCM
        filename= 'INPUT/surf.res.nc' 
        fieldname= 'topo'

        if (file_exist( trim(filename) )) then
            mountain = .true.

           if(master) then
              allocate ( g_phis(nlon,mlat) )
              call read_data( trim(filename), trim(fieldname), g_phis, no_domain=.true. )
           endif

           call mp_scatter3d(g_phis, phis,  nlon, mlat,  1,  &
                             beglat, endlat, 1, 1, ighost, ighost, 0)

           call fv_array_sync
           if(master) deallocate ( g_phis )

            phis = phis * grav
        else
            mountain = .false.
            if(master) write(6,*) 'No surface dataset found; setting phis = 0'
            phis = 0.
        endif
#else

        !Modify the filename if it has an ensemble instance appendix
        call get_instance_filename('INPUT/surf.res',restart_name )
        restart_name=trim(restart_name)

        if (file_exist(restart_name)) then
            mountain = .true.
            if(master) call mpp_open(iuhs, restart_name, action=MPP_RDONLY,        &
                 form=MPP_NATIVE, access=MPP_SEQUENTIAL)
#ifdef SPMD
            if(master) then
                allocate ( g_phis(nlon,mlat) )
                read(iuhs) g_phis
            endif
            call mp_scatter3d(g_phis, phis,  nlon, mlat,  1,  &
                 beglat, endlat, 1, 1, ighost, ighost, 0)
            call fv_array_sync
            if(master) deallocate ( g_phis )
#else
            read(iuhs) phis
#endif
            if(master) close (iuhs)
            phis = phis * grav
        else
            mountain = .false.
            if(master) write(6,*) 'No surface dataset found; setting phis = 0'
            phis = 0.
        endif

#endif MARS_GCM

        call set_eta(nlev, ks, ptop, ak, bk)
        call init_dry_atm(mountain, kappa, grav, rdgas)

        if ( age_tracer ) then
            age_time = 0.
!$omp parallel do private (i, j, k)
!         do k=1,nlev
            do k=ksp,kep
               do j=beglat, endlat
                  do i=1,nlon
                     q(i,j,k,ncnst) = 0.
                  enddo
               enddo
            enddo
            call fv_array_sync()
        endif

#endif
    else
!---------------------------
! Read from the restart file
!---------------------------
        restart_name ='INPUT/fv_rst.res'
        !The call get_instance_filename() is inside the read_fv_rst subroutine 
        call read_fv_rst(im, jm, km, ncnst, restart_name, days, seconds, restart_format)

!   call p_var(nlon, mlat, nlev, beglat, endlat, ptop, delp, ps,           &
!              pe, peln,  pk,  pkz,  kappa, q, ng_d, ncnst, adjust_dry_mass )
    endif

#ifdef SPMD
! Ghost phis boundary N/S
    if ( ighost /= 0 ) then
        call mp_send4d_ns(nlon, mlat, 1, 1, beglat, endlat,       &
             1, 1, ighost, ighost, phis)
        call mp_recv4d_ns(nlon, mlat, 1, 1, beglat, endlat,       &
             1, 1, ighost, ighost, phis)
        call fv_array_sync()
    endif
#endif

!------------------------------------------------
! Safety check of the (a,b) of the eta coordinate
!------------------------------------------------

    call check_eta(nlev, ks, ptop, ak, bk)

! Generate winds on the mass grid for the physics/coupler:
    if ( use_tendency ) then
        call d2a3d(u, v, ua, va, nlon, mlat, nlev, beglat, endlat,  &
             ng_d, ng_s, coslon, sinlon)
    endif
    call fv_print_chksums( 'Exiting  fv_restart' )
  end subroutine fv_restart
! $Id: fv_restart.F90,v 17.0 2009/07/21 02:53:25 fms Exp $

  subroutine write_fv_rst( dflnm, days, seconds, grav, oform)

#include "fv_arrays.h"

! !INPUT PARAMETERS:
    integer,  intent(in):: days
    integer,  intent(in):: seconds
    real, intent(in):: grav
    character*(*), intent(in):: dflnm           ! dynamics restart filename
    character(len=*), intent(in) :: oform

! local
    integer :: iuic                 ! Unit number
    real, allocatable :: w2d(:,:), w3d(:,:,:), latu(:), phis_alloc(:,:), delp_alloc(:,:,:)
    real :: ginv, qmax, qmin
    real :: pfull(nlev), phalf(nlev+1)
    integer i, j, k, n

    character(len=128)           :: tracer_name, tracer_longname, tracer_units
    type(axistype)               :: axis_lon, axis_lat, axis_lev, axis_levh, &
                                    axis_lonv, axis_latu, axis_time
    type(fieldtype), allocatable :: tracer_fields(:)
    type(fieldtype)              :: u_field, v_field, t_field, phis_field, delp_field

    type(fieldtype)              :: ak_field, bk_field

    character(len=*), parameter :: fms_tracers_restart_file = &
         'RESTART/'//fms_tracers_file

    character(len=128) :: directory, restart_name, tr_name
    logical :: use_native_format
    integer :: ids,ide, jds,jde, ics,ice,jcs,jce
    integer :: nsphum, nliq_wat, nice_wat, ncld_amt, ntracers
    real :: dummy(1) =1.0

#ifdef MARS_GCM
!rjw     fv_tr_name is not used in this subroutine      
#else
    character(len=20) :: fv_tr_name(4)
    data fv_tr_name / "sphum","liq_wat","ice_wat","cld_amt" /
#endif MARS_GCM

#include "fv_point.inc"

! output format
    if (trim(uppercase(oform)) == 'NATIVE') then
        use_native_format = .true.
    else if (trim(uppercase(oform)) == 'NETCDF') then
        use_native_format = .false.
    else
        call error_mesg ('write_fv_rst', &
             'invalid value for argument oform='//trim(oform), FATAL)
    endif

! check sum output data
  call print_check_sum ('Check sums for FV output data:')

    if (use_native_format) then
        if ( master ) then
!!$         if (file_exist(dflnm)) then
!!$             call mpp_open ( iuic, dflnm, action=MPP_OVERWR,        &
!!$                             form=MPP_NATIVE, access=MPP_SEQUENTIAL )
!!$         else
!!$             call mpp_open ( iuic, dflnm, action=MPP_WRONLY,        &
!!$                             form=MPP_NATIVE, access=MPP_SEQUENTIAL )
!!$            
!!$         endif
           !Modify the filename if it has an ensemble instance appendix
            call get_instance_filename(dflnm,restart_name )
            restart_name=trim(restart_name)
            iuic = open_restart_file(restart_name, 'write')
            write(iuic) nlon, mlat, nlev, nt_phys
            write(iuic) days, seconds
            write(iuic) ak, bk, ks
            write(6,*) ' '
            write(6,*) 'FV restart file header written for days=', days, 'sec=',seconds
        endif

        call pmaxmin('PS', ps(1,beglat), qmin, qmax, nlon, endlat-beglat+1, 0.01)

!---------------------------
! Write surface geopotential
!---------------------------
#ifdef SPMD
        allocate ( w2d(nlon,mlat) )
        call mp_gather3d(phis, w2d, nlon, mlat, 1, beglat, endlat, 1, 1,  &
             ighost, ighost, 0)
        if ( master ) then
            write(iuic) w2d
            ginv = 1./ grav
            call maxmin_global('Surface height', w2d, nlon, 1, mlat, ginv)
        endif
#else
        write(iuic) phis
#endif

        allocate ( w3d(nlon,mlat,nlev) )

!---------------------------
! Write u-wind
!---------------------------

#ifdef SPMD
        call mp_gather3d(u, w3d, nlon, mlat, nlev, beglat, endlat, 1, nlev, ng_d, ng_s, 0)
        if ( master ) then
            write(iuic) w3d
            call maxmin_global('U', w3d, nlon, mlat, nlev, 1.)
        endif
#else
        write(iuic) u
#endif

!---------------------------
! Write v-wind
!---------------------------

#ifdef SPMD
        call mp_gather3d(v, w3d, nlon, mlat, nlev, beglat, endlat, 1, nlev, ng_d, ng_d, 0)
        if ( master ) then
            write(iuic) w3d
            call maxmin_global('V', w3d, nlon, mlat, nlev, 1.)
        endif
#else
        write(iuic) v
#endif

!------------------------------------
! Write virtual potential temperature
!------------------------------------

#ifdef SPMD
        call mp_gather3d(pt, w3d, nlon, mlat, nlev, beglat, endlat, 1, nlev, ng_d, ng_d, 0)
        if ( master ) then
            write(iuic) w3d
            call maxmin_global('T', w3d, nlon, mlat, nlev, 1.)
        endif
#else
        write(iuic) pt
#endif

!------------------------------------
! Write pressure thickness (pascal)
!------------------------------------

#ifdef SPMD
        call mp_gather3d(delp, w3d, nlon, mlat, nlev, beglat, endlat, 1, nlev, 0, 0, 0)
        if ( master ) then
            write(iuic) w3d
!            call maxmin_global('DELP', w3d, nlon, mlat, nlev, 1.)
        endif
#else
        write(iuic) delp
#endif


        if(ncnst /=  0) then 
!------------------------------------
! Write tracers one at a time
!------------------------------------

            do n=1,nt_phys 

#ifdef SPMD
               call mp_gather3d(q(1,beglat-ng_d,1,n), w3d, nlon, mlat, nlev,    &
                    beglat, endlat, 1, nlev, ng_d, ng_d, 0)
#else
               do k=1,nlev
                  do j=beglat,endlat
                     do i=1,nlon
                        w3d(i,j,k) = q(i,j,k,n)
                     enddo
                  enddo
               enddo
#endif
               if ( master ) then
                   write(iuic) w3d
                   call maxmin_global('Q', w3d, nlon, mlat, nlev, 1.)
               endif
            enddo
            if ( master .and. age_tracer ) write(iuic) age_time
        endif

        if( master ) then 
            call close_file(iuic)
            write(*,*) 'FV restart file written and closed'
        endif

!--------------------------------------------------------------
! Write surface winds on A grid for coupler to restart properly
!--------------------------------------------------------------

        if ( master ) then
           !Modify the filename if it has an ensemble instance appendix
            call get_instance_filename('RESTART/fv_srf_wnd.res',restart_name )
            restart_name=trim(restart_name)
            iuic = open_restart_file(restart_name, 'write')
            write(iuic) nlon, mlat
        endif

#ifdef SPMD
        call mp_gather3d(ps_bp, w2d, nlon, mlat, 1, beglat, endlat, 1, 1,  &
             0, 0, 0)
#endif
        if (master) write(iuic) w2d

#ifdef SPMD
        call mp_gather3d(u_srf, w2d, nlon, mlat, 1, beglat, endlat, 1, 1,  &
             0, 0, 0)
#endif
        if (master) write(iuic) w2d

#ifdef SPMD
        call mp_gather3d(v_srf, w2d, nlon, mlat, 1, beglat, endlat, 1, 1,  &
             0, 0, 0)
#endif
        if (master) then
            write(iuic) w2d
            call close_file(iuic)
        endif

!----------------------------------------
! Write additional tracers to netCDF file
!----------------------------------------

        if (do_fms_tracer_manager .and. ncnst > nt_phys) then !{

            if(master) &
                 & write(*,*)'Additional tracers present in field_table. These will be saved in a netcdf file.'

#ifdef MARS_GCM
            call get_eta_level(nlev, p_ref, pfull, phalf, 0.01)
#else
            call get_eta_level(nlev, 1.E5, pfull, phalf, 0.01)
#endif MARS_GCM

            !Modify the filename if it has an ensemble instance appendix
            call get_instance_filename(fms_tracers_restart_file,restart_name )
            restart_name=trim(restart_name)
            call mpp_open(iuic, restart_name, &
                 action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

! define axes
            call mpp_write_meta(iuic, axis_lon, 'lon', 'degree_E', 'longitude', &
                 data=lon )
            call mpp_write_meta(iuic, axis_lat, 'lat', 'degree_N', 'latitude', &
                 data=lat )
            call mpp_write_meta(iuic, axis_lev, 'pfull', 'mb', 'ref full pressure level', &
                 data=pfull )

! define fields
            allocate(tracer_fields(ncnst))

            do n = nt_phys+1, ncnst
               call tr_get_tracer_names(MODEL_ATMOS, n, &
                    tracer_name, tracer_longname, tracer_units)
               call mpp_write_meta(iuic, tracer_fields(n), &
                    (/axis_lon, axis_lat, axis_lev/), &
                    tracer_name, tracer_units, tracer_longname, pack=1) ! want doubles
            enddo

! write axes
            call mpp_write(iuic, axis_lon)
            call mpp_write(iuic, axis_lat)
            call mpp_write(iuic, axis_lev)

! write tracers
            do n = nt_phys+1, ncnst
#ifdef SPMD
               call mp_gather3d(q(1,beglat-ng_d,1,n), w3d, nlon, mlat, nlev,    &
                    beglat, endlat, 1, nlev, ng_d, ng_d, 0)
#else
               do k=1,nlev
                  do j=beglat,endlat
                     do i=1,nlon
                        w3d(i,j,k) = q(i,j,k,n)
                     enddo
                  enddo
               enddo
#endif
               if(master) call mpp_write(iuic, tracer_fields(n), data=w3d)

            enddo

            call mpp_close(iuic)
            deallocate(tracer_fields)

        endif !}


        deallocate ( w3d )
        deallocate ( w2d )

    else ! ! ------------------NetCDF restart files---------
        call mpp_get_data_domain   ( fv_domain, ids, ide, jds, jde )
        call mpp_get_compute_domain( fv_domain, ics, ice, jcs, jce )
        restart_name = trim(dflnm) ! add netcdf file suffix

        !Modify the filename if it has an ensemble instance appendix
        call get_instance_filename(restart_name,restart_name )
        restart_name=trim(restart_name)
            
        call mpp_open(iuic, restart_name, &
             action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

! define axes
#ifdef MARS_GCM
        call get_eta_level(nlev, p_ref, pfull, phalf, 0.01)
#else
        call get_eta_level(nlev, 1.E5, pfull, phalf, 0.01)
#endif MARS_GCM
!        if(master) write(*,*) "phalf values are ",phalf    
! define longitude axis for u and tracer fields
        call mpp_write_meta(iuic, axis_lon, 'lon', 'degree_E', 'longitude', &
             data=lon )

! define longitude axis for v field
        call mpp_write_meta(iuic, axis_lonv, 'lonv', 'degree_E', 'longitude', &
             data=lonb(1:size(lonb)-1) )


! define latitude axis for v and tracer fields
        call mpp_write_meta(iuic, axis_lat, 'lat', 'degree_N', 'latitude', &
             data=lat )

! define latitude axis for u field
        allocate(latu(size(lat)))
        ! The first latitude corresponds to meaningless data for the South Pole. 
        ! For interpolation reasons, it has been decided to set this latitude to 
        ! the second latitude with a small offset for Ferret compatibility.
        latu(1)  = latb(2) - 1e-5
        latu(2:) = latb(2:size(lat))
        call mpp_write_meta(iuic, axis_latu, 'latu', 'degree_N', 'latitude', &
             data=latu )

! define pressure axis
        call mpp_write_meta(iuic, axis_lev, 'pfull', 'mb', 'ref full pressure level', &
             data=pfull )
        call mpp_write_meta(iuic, axis_levh, 'phalf', 'mb', 'ref half pressure level', &
             data=phalf )

        call mpp_write_meta(iuic, axis_time, 'Time', 'time level', 'Time', &
                            cartesian='T', data=dummy)

        call mpp_write_meta(iuic, ak_field, (/axis_levh/), &
             'ak', 'Pa', 'ak-field', pack=1) ! want doubles
        call mpp_write_meta(iuic, bk_field, (/axis_levh/), &
             'bk', ' ', 'bk-field', pack=1)

        call mpp_write_meta(iuic, 0, 'nlon', ival=nlon)     
        call mpp_write_meta(iuic, 0, 'mlat', ival=mlat)     
        call mpp_write_meta(iuic, 0, 'nlev', ival=nlev)     

        call mpp_write_meta(iuic, 0, 'nt_phys', nt_phys)     
        call mpp_write_meta(iuic, 0, 'days', days)     
        call mpp_write_meta(iuic, 0, 'seconds', seconds)     
        call mpp_write_meta(iuic, 0, 'ks', ks)     

        call mpp_write_meta(iuic, u_field, (/axis_lon , axis_latu, axis_lev, axis_time/), &
             'U', 'm/s', 'u-field', pack=1)
        call mpp_write_meta(iuic, v_field, (/axis_lonv, axis_lat , axis_lev, axis_time/), &
             'V', 'm/s', 'v-field', pack=1)
        call mpp_write_meta(iuic, t_field, (/axis_lon , axis_lat , axis_lev, axis_time/), &
             'T', 'K', 'Temperature', pack=1)
        call mpp_write_meta(iuic, phis_field, (/axis_lon, axis_lat, axis_time/), &
             'Surface_geopotential', 'm**2/s**2', 'Surface Height', pack=1)
        call mpp_write_meta(iuic, delp_field, (/axis_lon, axis_lat, axis_lev, axis_time/), &
             'DELP', 'pa', 'delp-field', pack=1)

! write axes
        call mpp_write(iuic, axis_time)
        call mpp_write(iuic, axis_lon)
        call mpp_write(iuic, axis_lonv)
        call mpp_write(iuic, axis_lat)
        call mpp_write(iuic, axis_latu)
        call mpp_write(iuic, axis_lev)

        call mpp_write(iuic, bk_field, bk)     
        call mpp_write(iuic, ak_field, ak)     

        call mpp_write(iuic, u_field, fv_domain, u(ids:ide,jds:jde,:))
        call mpp_write(iuic, v_field, fv_domain, v(ids:ide,jds:jde,:))
        call mpp_write(iuic, t_field, fv_domain, pt(ids:ide,jds:jde,:))
! The mpp_io code uses Cray pointers which doesn't work well with the PSET 
! 2D arrays of the FV code. Therefore a copy of the 2D array is needed before 
! passing the array to the mpp_write routine.
        allocate (phis_alloc(ids:ide,jds:jde))
        phis_alloc = -999.9
        phis_alloc(ids:ide,jcs:jce) = phis(ids:ide,jcs:jce)
        call mpp_write(iuic, phis_field, fv_domain, phis_alloc)
        allocate (delp_alloc(ids:ide,jds:jde,nlev))
        delp_alloc(ids:ide,jcs:jce,:) = delp(ids:ide,jcs:jce,:)   !RASF Bug fix for array bounds
        call mpp_write(iuic, delp_field, fv_domain, delp_alloc)   !RASF Bug fix for array bounds
        deallocate(delp_alloc)


        call mpp_close(iuic)


        call error_mesg ('write_fv_rst', 'Writing NetCDF formatted restart file: '//trim(restart_name), NOTE)

        call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)


        restart_name = trim(fms_tracers_restart_file)!//'1.nc' ! add netcdf file suffix
        !Modify the filename if it has an ensemble instance appendix
        call get_instance_filename(restart_name,restart_name )
        restart_name=trim(restart_name)

        call mpp_open(iuic, restart_name, &
             action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

! define axes

        call mpp_write_meta(iuic, axis_time, 'Time', 'time level', &
                            'Time',cartesian='T', data=dummy)

! define longitude axis for tracer fields
        call mpp_write_meta(iuic, axis_lon, 'lon', 'degree_E', 'longitude', &
             data=lon )

! define latitude axis for tracer fields
        call mpp_write_meta(iuic, axis_lat, 'lat', 'degree_N', 'latitude', &
             data=lat )

! define pressure axis
        call mpp_write_meta(iuic, axis_lev, 'pfull', 'mb', 'ref full pressure level', &
             data=pfull )

! define fields
        allocate(tracer_fields(ntracers))

        do n = 1,ntracers
           call tr_get_tracer_names(MODEL_ATMOS, n, &
                tracer_name, tracer_longname, tracer_units)
           call mpp_write_meta(iuic, tracer_fields(n), &
                (/axis_lon, axis_lat, axis_lev, axis_time/), &
                tracer_name, tracer_units, tracer_longname, pack=1)
        enddo

! write axes
        call mpp_write(iuic, axis_lon)
        call mpp_write(iuic, axis_lat)
        call mpp_write(iuic, axis_lev)
        call mpp_write(iuic, axis_time)

        do n = 1, ntracers
           call mpp_write(iuic, tracer_fields(n), fv_domain, q(ids:ide,jds:jde,:,n))
        enddo

        deallocate(tracer_fields)
        call mpp_close(iuic)


        !Modify the filename if it has an ensemble instance appendix
        call get_instance_filename('RESTART/fv_srf_wnd.res',restart_name )
        restart_name=trim(restart_name)
        call mpp_open(iuic, restart_name, &
             action=MPP_OVERWR, form=MPP_NETCDF, threading=MPP_SINGLE )

! define axes

! define longitude axis for u and tracer fields
        call mpp_write_meta(iuic, axis_lon, 'lon', 'degree_E', 'longitude', &
             data=lon )

! define longitude axis for v field
        call mpp_write_meta(iuic, axis_lonv, 'lonv', 'degree_E', 'longitude', &
             data=lonb(1:size(lonb)-1) )


! define latitude axis for v and tracer fields
        call mpp_write_meta(iuic, axis_lat, 'lat', 'degree_N', 'latitude', &
             data=lat )

! define latitude axis for u field
        latu(1)  = lat(1)
        latu(2:) = latb(2:size(lat))
        call mpp_write_meta(iuic, axis_latu, 'latu', 'degree_N', 'latitude', &
             data=latu )

        call mpp_write_meta(iuic, 'nlon', ival=nlon)     
        call mpp_write_meta(iuic, 'mlat', ival=mlat)     

        call mpp_write_meta(iuic, u_field, (/axis_lon , axis_latu/), &
             'U_srf', 'm/s', 'surface u-field', pack=1)
        call mpp_write_meta(iuic, v_field, (/axis_lonv, axis_lat /), &
             'V_srf', 'm/s', 'surface v-field', pack=1)
        call mpp_write_meta(iuic, phis_field, (/axis_lon, axis_lat/), &
             'ps_bp', 'pa', 'surface pressure', pack=1)

! write axes
        call mpp_write(iuic, axis_lon)
        call mpp_write(iuic, axis_lonv)
        call mpp_write(iuic, axis_lat)
        call mpp_write(iuic, axis_latu)
! The mpp_io code uses Cray pointers which doesn't work well with the PSET 
! 2D arrays of the FV code. Therefore a copy of the 2D array is needed before 
! passing the array to the mpp_write routine.
        phis_alloc = -999.9
        phis_alloc(ids:ide,jcs:jce) = ps_bp(ids:ide,jcs:jce)
        call mpp_write(iuic, phis_field, fv_domain, phis_alloc)
        phis_alloc(ids:ide,jcs:jce) = u_srf(ids:ide,jcs:jce)
        call mpp_write(iuic, u_field, fv_domain, phis_alloc)
        phis_alloc(ids:ide,jcs:jce) = v_srf(ids:ide,jcs:jce)
        call mpp_write(iuic, v_field, fv_domain, phis_alloc)




    endif ! end IF NETCDF
    call fv_array_sync()

  end subroutine write_fv_rst
!---------------------------------------------------------------------

  subroutine maxmin_global(qname, a, im, jm, km, fac)

    implicit none

    character*(*)  qname
    integer im, jm, km
    real a(im,jm, km)
    real fac
    real pmax(km), pmin(km)
    real qmax, qmin
    integer i, j, k

!$omp parallel do private(i, j, k)
!Balaji: not parallelized, called with unshared arrays
    do k=1, km
       pmax(k) = a(1,1,k)
       pmin(k) = a(1,1,k)
       do j=1,jm
          do i=1, im
             pmax(k) = max (pmax(k), a(i,j,k))
             pmin(k) = min (pmin(k), a(i,j,k))
          enddo
       enddo
    enddo

    qmax = pmax(1)
    qmin = pmin(1)
    do k=2,km
       qmax = max (qmax, pmax(k))
       qmin = min (qmin, pmin(k))
    enddo

    write(6,*) qname, ' max = ', qmax * fac, ' min = ', qmin * fac
  end subroutine maxmin_global


subroutine print_check_sum (label)

#include "fv_arrays.h"

 character(len=*),      intent(in) :: label
 integer :: n, chksum_unit
 integer :: isl, iel, jsl, jel
 integer(I8_KIND) :: zsum
#ifdef use_shared_pointers
#include "fv_point.inc"
#endif

  call mpp_get_compute_domain( fv_domain, isl, iel, jsl, jel )
  chksum_unit = stdout()

  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,'(/a)') trim(label)
! sum all restart variables (except eta,peta)
  zsum = mpp_chksum(pt(isl:iel,jsl:jel,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(t)',zsum
  zsum = mpp_chksum(u(isl:iel,jsl:jel,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(u)',zsum
  zsum = mpp_chksum(v(isl:iel,jsl:jel,:))
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(v)',zsum
  zsum = 0
! sum of all tracers
  do n = 1, size(q,4)
    zsum = zsum + mpp_chksum(q(isl:iel,jsl:jel,:,n))
  enddo
  if (mpp_pe() == mpp_root_pe()) write (chksum_unit,10) '(r)',zsum
10 format ('chksum',a7,' = ',z16)

end subroutine print_check_sum

end module fv_restart_mod
