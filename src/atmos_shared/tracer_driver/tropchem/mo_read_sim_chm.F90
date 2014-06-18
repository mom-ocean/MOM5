      module MO_READ_SIM_CHM_MOD

      use mpp_mod,    only : mpp_error, FATAL
      use mpp_io_mod, only : mpp_open, MPP_RDONLY, MPP_ASCII,MPP_MULTI, &
                             MPP_SINGLE, mpp_close

implicit none
character(len=128), parameter :: version     = '$Id: mo_read_sim_chm.F90,v 19.0 2012/01/06 20:34:02 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      CONTAINS
        
      subroutine READ_SIM_CHM( sim_data_flsp, &
                               sim_file_cnt )
!--------------------------------------------------------
!            ... Initialize chemistry modules
!--------------------------------------------------------

      use CHEM_MODS_MOD,     only : explicit, implicit, rodas, grpcnt, &
                                nadv_mass, adv_mass, pcnstm1, &
                                drydep_cnt, drydep_lst, &
                                srfems_cnt, srfems_lst, &
                                hetcnt, het_lst, extcnt, extfrc_lst, &
                                rxt_alias_cnt, rxt_alias_lst, rxt_alias_map, &
                                ngrp, grp_mem_cnt, grp_lst
      use M_TRACNAME_MOD,    only : tracnam, natsnam

      implicit none

!--------------------------------------------------------
!            ... Dummy args
!--------------------------------------------------------
      integer, intent(out) :: sim_file_cnt
      character(len=32), intent(in) :: sim_data_flsp

!--------------------------------------------------------
!            ... Local variables
!--------------------------------------------------------
      integer, parameter :: inst = 1, avrg = 2
      integer, parameter :: max_hst_ind = 17
      integer  ::  ios, funit
      character(len=128) :: msg

!     funit = NAVU()
!--------------------------------------------------------
!            ... Open chem input unit
!--------------------------------------------------------
!     OPEN( unit = funit, &
!           file = TRIM( sim_data_flsp ), &
!           status = 'old', &
!           recl   = 2048, &
!           iostat = ios )
!     if( ios /= 0 ) then
!        write(*,*) ' READ_SIM_CHM: Failed to open file ',TRIM( sim_data_flsp )
!        write(*,*) ' Error code = ',ios
!        call ENDRUN
!     end if
      
      call mpp_open( funit, trim(sim_data_flsp), MPP_RDONLY,MPP_ASCII, threading = MPP_MULTI, &
                     fileset=MPP_SINGLE, recl = 2048)

        
!--------------------------------------------------------
!        ... Read map info from data file
!--------------------------------------------------------
      if( explicit%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) explicit%cls_rxt_cnt
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read explicit cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) explicit%clsmap
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read explicit clscnt; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( implicit%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) implicit%cls_rxt_cnt
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%clsmap
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit clscnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%permute
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit permute; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) implicit%diag_map
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read implicit diag_map; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( rodas%clscnt > 0 ) then
         read(funit,'(4i4)',iostat=ios) rodas%cls_rxt_cnt
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas cls_rxt_cnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%clsmap
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas clscnt; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%permute
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas permute; error = ', ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rodas%diag_map
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rodas diag_map; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( pcnstm1 > 0 ) then
         read(funit,*,iostat=ios) adv_mass(:pcnstm1)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read adv_mass; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( grpcnt > 0 ) then
         read(funit,*,iostat=ios) nadv_mass(:grpcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read nadv_mass; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( pcnstm1 > 0 ) then
         read(funit,'(10a8)',iostat=ios) tracnam(:pcnstm1)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read tracnam; error = ', ios
            call ENDRUN(msg)
         end if
      end if
      if( grpcnt > 0 ) then
         read(funit,'(i4)',iostat=ios) ngrp
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read ngrp; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( grp_mem_cnt(ngrp),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate grp_mem_cnt; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( grp_lst(ngrp),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate grp_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) grp_mem_cnt(:ngrp)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read grp_mem_cnt; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) grp_lst(:ngrp)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read grp_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) natsnam(1:grpcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read natsnam; error = ',ios
            call ENDRUN(msg)
         end if
      end if
      read(funit,'(i4)',iostat=ios) srfems_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read srfems_cnt; error = ',ios
            call ENDRUN(msg)
         end if
      if( srfems_cnt > 0 ) then
         allocate( srfems_lst(srfems_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate srfems_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) srfems_lst(1:srfems_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read srfems_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if
      read(funit,'(i4)',iostat=ios) drydep_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read drydep_cnt; error = ',ios
            call ENDRUN(msg)
      end if
      if( drydep_cnt > 0 ) then
         allocate( drydep_lst(drydep_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate drydep_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(10a8)',iostat=ios) drydep_lst(1:drydep_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read drydep_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if
      if( hetcnt > 0 ) then
         read(funit,'(10a8)',iostat=ios) het_lst(1:hetcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read het_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if
      if( extcnt > 0 ) then
         read(funit,'(10a8)',iostat=ios) extfrc_lst(1:extcnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read extfrc_lst; error = ',ios
            call ENDRUN(msg)
         end if
      end if
      read(funit,'(i4)',iostat=ios) rxt_alias_cnt
      if( ios /= 0 ) then
         write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_cnt; error = ',ios
            call ENDRUN(msg)
      end if
      if( rxt_alias_cnt > 0 ) then
         allocate( rxt_alias_lst(rxt_alias_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate rxt_alias_lst; error = ',ios
            call ENDRUN(msg)
         end if
         allocate( rxt_alias_map(rxt_alias_cnt),stat=ios )
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to allocate rxt_alias_map; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(5a16)',iostat=ios) rxt_alias_lst(1:rxt_alias_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_lst; error = ',ios
            call ENDRUN(msg)
         end if
         read(funit,'(20i4)',iostat=ios) rxt_alias_map(1:rxt_alias_cnt)
         if( ios /= 0 ) then
            write(msg,*) 'READ_SIM_CHM: Failed to read rxt_alias_map; error = ',ios
            call ENDRUN(msg)
      end if
      end if

!     read(funit,'(i4)',iostat=ios) moz_file_cnt
!     if( ios /= 0 ) then
!        write(msg,*) 'READ_SIM_CHM: Failed to read moz_file_cnt; error = ',ios
!            call ENDRUN(msg)
!     end if
!     sim_file_cnt = MAX( moz_file_cnt,match_file_cnt )
!     do file = 1,moz_file_cnt
!        read(funit,'(10i4)',iostat=ios) hfile(file)%histout_cnt(:,:)
!        if( ios /= 0 ) then
!           write(msg,*) 'READ_SIM_CHM: Failed to read histout_cnt for file ',file,'; error = ',ios
!            call ENDRUN(msg)
!        end if
!     end do

!     do file = 1,sim_file_cnt
!        do k = 1,hstdim
!            if( hstinst(file)%list(k) == ' ' ) then
!               exit
!            end if
!        end do
!        hfile(file)%match_cnt(1) = k - 1
!        do k = 1,hstdim
!            if( hsttimav(file)%list(k) == ' ' ) then
!               exit
!            end if
!        end do
!        hfile(file)%match_cnt(2) = k - 1
!        do i = inst,avrg
!           moz_cnt(i) = SUM( hfile(file)%histout_cnt(:,i) )
!        end do
!        hfile(file)%mxoutflds = MAX( hfile(file)%match_cnt(1)+moz_cnt(1), &
!                                     hfile(file)%match_cnt(2)+moz_cnt(2) )
!        ALLOCATE( minst(hfile(file)%mxoutflds), mtimav(hfile(file)%mxoutflds), stat=astat )
!        if( astat /= 0 ) then
!           write(msg,*) 'READ_SIM_CHM: Failed to allocate minst,mtimav; error = ',astat
!            call ENDRUN(msg)
!        end if
!        ALLOCATE( hfile(file)%outinst(hfile(file)%mxoutflds), &
!                  hfile(file)%outtimav(hfile(file)%mxoutflds), stat=astat )
!        if( astat /= 0 ) then
!           write(msg,*) 'READ_SIM_CHM: Failed to allocate outinst,outtimav for file ',file,'; error = ',astat
!            call ENDRUN(msg)
!        end if
!        minst(:)  = 0
!        mtimav(:) = 0
!        if( hfile(file)%match_cnt(1) > 0 ) then
!           hfile(file)%outinst(:hfile(file)%match_cnt(1)) = hstinst(file)%list(:hfile(file)%match_cnt(1))
!        end if
!        if( hfile(file)%match_cnt(2) > 0 ) then
!           hfile(file)%outtimav(:hfile(file)%match_cnt(2)) = hsttimav(file)%list(:hfile(file)%match_cnt(2))
!        end if

!        do i = inst,avrg
!           end = hfile(file)%match_cnt(i)
!            endi = 0
!           do k = 1,max_hst_ind
!               if( hfile(file)%histout_cnt(k,i) /= 0 ) then
!                  start  = end + 1
!                  end    = start + hfile(file)%histout_cnt(k,i) - 1
!                  starti = endi + 1
!                  endi   = starti + hfile(file)%histout_cnt(k,i) - 1
!                  hfile(file)%histout_ind(k,i) = starti
!                  if( i == inst ) then
!                    read(funit,'(4a32)',iostat=ios) hfile(file)%outinst(start:end)
!                    read(funit,'(20i4)',iostat=ios) minst(starti:endi)
!                  else if( i == avrg ) then
!                    read(funit,'(4a32)',iostat=ios) hfile(file)%outtimav(start:end)
!                    read(funit,'(20i4)',iostat=ios) mtimav(starti:endi)
!                  end if
!               end if
!           end do
!            if( endi > 0 ) then
!               if( i == inst ) then
!                  ALLOCATE( hfile(file)%hist_inst(end), stat=astat )
!                  if( astat /= 0 ) then
!                     write(msg,*) ' READ_SIM_CHM: Failed to allocate hist_ind for file ',file,'; error = ',astat
!                     call ENDRUN(msg)
!                  end if
!                  hfile(file)%hist_inst(:end-hfile(file)%match_cnt(i)) = &
!                           hfile(file)%outinst(hfile(file)%match_cnt(i)+1:end)
!                  ALLOCATE( hfile(file)%inst_map(endi), stat=astat )
!                  if( astat /= 0 ) then
!                     write(msg,*) ' READ_SIM_CHM: Failed to allocate inst_map for file ',file,'; error = ',astat
!                    call ENDRUN(msg)
!                  end if
!                  hfile(file)%inst_map(:endi) = minst(:endi)
!               else if( i == avrg ) then
!                  ALLOCATE( hfile(file)%hist_timav(end), stat=astat )
!                  if( astat /= 0 ) then
!                     write(msg,*) ' READ_SIM_CHM: Failed to allocate hist_timav for file ',file,'; error = ',astat
!                     call ENDRUN(msg)
!                  end if
!                  hfile(file)%hist_timav(:end-hfile(file)%match_cnt(i)) = &
!                         hfile(file)%outtimav(hfile(file)%match_cnt(i)+1:end)
!                  ALLOCATE( hfile(file)%timav_map(endi), stat=astat )
!                  if( astat /= 0 ) then
!                     write(msg,*) ' READ_SIM_CHM: Failed to allocate timav_map for file ',file,'; error = ',astat
!                     call ENDRUN(msg)
!                  end if
!                  hfile(file)%timav_map(:endi) = mtimav(:endi)
!               end if
!            end if
!        end do
!         DEALLOCATE( minst, mtimav )
!     end do
!     read(funit,'(i3)') ndiags

      call mpp_close( funit )

!      write(*,*) '---------------------------------------------------------------------------------'
!      write(*,*) ' '

      end subroutine READ_SIM_CHM

      subroutine ENDRUN(msg)
         character(len=128), intent(in) :: msg
         call mpp_error(FATAL, msg)
      end subroutine ENDRUN        

      end module MO_READ_SIM_CHM_MOD
