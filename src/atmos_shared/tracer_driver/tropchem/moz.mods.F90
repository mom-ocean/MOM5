      module mo_grid_mod
!---------------------------------------------------------------------
!       ... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none

      save

character(len=128), parameter :: version     = '$Id: moz.mods.F90,v 19.0 2012/01/06 20:34:14 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: tikal $'
logical                       :: module_is_initialized = .false.

      integer, parameter :: &
                pcnst    =    84+1, &     ! number of advected constituents including cloud water
                pcnstm1  =    84, &     ! number of advected constituents excluding cloud water
                plev     =   1, &         ! number of vertical levels
                plevp    = plev+1, &      ! plev plus 1
                plevm    = plev-1, &      ! plev minus 1
                plon     =   1, &         ! number of longitudes
                plat     =   1            ! number of latitudes

      integer, parameter :: &
                pnats    =     0    ! number of non-advected trace species



      integer :: nodes                ! mpi task count
      integer :: plonl                ! longitude tile dimension
      integer :: pplon                ! longitude tile count
      integer :: plnplv               ! plonl * plev

      end module mo_grid_mod

      module chem_mods_mod
!--------------------------------------------------------------
!       ... basic chemistry array parameters
!--------------------------------------------------------------

      use mo_grid_mod, only : pcnstm1
!++lwh
      use mpp_mod,     only : mpp_error, FATAL
!--lwh

      implicit none

      save

      integer, parameter :: hetcnt     =     0, &    ! number of heterogeneous processes
                            phtcnt     =    41, &    ! number of photo processes
                            rxntot     =   224, &    ! number of total reactions
                            gascnt     =   183, &    ! number of gas phase reactions
                            nfs        =     3, &       ! number of "fixed" species
                            relcnt     =     0, &    ! number of relationship species
                            grpcnt     =     0, &    ! number of group members
                            imp_nzcnt  =   829, &     ! number of non-zero implicit matrix entries
                            rod_nzcnt  =     0, &     ! number of non-zero rodas matrix entries
                            extcnt     =     0, &    ! number of species with external forcing
                            clscnt1    =     8, &  ! number of species in explicit class
                            clscnt2    =     0, &  ! number of species in hov class
                            clscnt3    =     0, &  ! number of species in ebi class
                            clscnt4    =    76, &  ! number of species in implicit class
                            clscnt5    =     0, &  ! number of species in rodas class
                            indexm     =     1, &    ! index of total atm density in invariant array
                            ncol_abs   =     2, &    ! number of column densities
                            indexh2o   =     0, &    ! index of water vapor density
                            clsze      = 1       ! loop length for implicit chemistry

      integer ::            ngrp       = 0
      integer ::            drydep_cnt = 0
      integer ::            srfems_cnt = 0
      integer ::            rxt_alias_cnt = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_alias_map(:)
      real      :: adv_mass(pcnstm1)
      real      :: nadv_mass(grpcnt)
      character(len=16), allocatable :: rxt_alias_lst(:)
      character(len=8), allocatable  :: drydep_lst(:)
      character(len=8), allocatable  :: srfems_lst(:)
      character(len=8), allocatable  :: grp_lst(:)
      character(len=8)               :: het_lst(max(1,hetcnt))
      character(len=8)               :: extfrc_lst(max(1,extcnt))
      character(len=8)               :: inv_lst(max(1,nfs))

      type solver_class
         integer :: clscnt
         integer :: lin_rxt_cnt
         integer :: nln_rxt_cnt
         integer :: indprd_cnt
         integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class

      type(solver_class) :: explicit, implicit, rodas

      contains

      subroutine endrun(msg)

      implicit none

      character(len=128), intent(in), optional  :: msg
      call mpp_error(FATAL, msg)

      end subroutine endrun 

      subroutine chem_mods_init
!--------------------------------------------------------------
!       ... intialize the class derived type
!--------------------------------------------------------------

      implicit none

      integer :: astat

      explicit%clscnt       =     8
      explicit%indprd_cnt   =    52

      implicit%clscnt       =    76
      implicit%lin_rxt_cnt  =    70
      implicit%nln_rxt_cnt  =   151
      implicit%indprd_cnt   =     3
      implicit%iter_max     =    11

      rodas%clscnt          =     0
      rodas%lin_rxt_cnt     =     0
      rodas%nln_rxt_cnt     =     0
      rodas%indprd_cnt      =     0

      if( explicit%clscnt > 0 ) then
         allocate( explicit%clsmap(explicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate explicit%clsmap ; error = ',astat
            call endrun
         end if
         explicit%clsmap(:)  = 0
      end if
      if( implicit%clscnt > 0 ) then
         allocate( implicit%permute(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate implicit%permute ; error = ',astat
            call endrun
         end if
         implicit%permute(:)  = 0
         allocate( implicit%diag_map(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate implicit%diag_map ; error = ',astat
            call endrun
         end if
         implicit%diag_map(:)  = 0
         allocate( implicit%clsmap(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate implicit%clsmap ; error = ',astat
            call endrun
         end if
         implicit%clsmap(:)  = 0
      end if
      if( rodas%clscnt > 0 ) then
         allocate( rodas%permute(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate rodas%permute ; error = ',astat
            call endrun
         end if
         rodas%permute(:)  = 0
         allocate( rodas%diag_map(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate rodas%diag_map ; error = ',astat
            call endrun
         end if
         rodas%diag_map(:)  = 0
         allocate( rodas%clsmap(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(*,*) 'chem_mods_init: failed to allocate rodas%clsmap ; error = ',astat
            call endrun
         end if
         rodas%clsmap(:)  = 0
      end if

      end subroutine chem_mods_init

      end module chem_mods_mod
  
      module M_SPC_ID_MOD
  
      implicit none                                                             
  
      integer, parameter :: id_O3 =   1
      integer, parameter :: id_O =   2
      integer, parameter :: id_O1D =   3
      integer, parameter :: id_N2O =   4
      integer, parameter :: id_N =   5
      integer, parameter :: id_NO =   6
      integer, parameter :: id_NO2 =   7
      integer, parameter :: id_NO3 =   8
      integer, parameter :: id_HNO3 =   9
      integer, parameter :: id_HO2NO2 =  10
      integer, parameter :: id_N2O5 =  11
      integer, parameter :: id_CH4 =  12
      integer, parameter :: id_CH3O2 =  13
      integer, parameter :: id_CH3OOH =  14
      integer, parameter :: id_CH2O =  15
      integer, parameter :: id_CO =  16
      integer, parameter :: id_OH =  17
      integer, parameter :: id_HO2 =  18
      integer, parameter :: id_H2O2 =  19
      integer, parameter :: id_C3H6 =  20
      integer, parameter :: id_ISOP =  21
      integer, parameter :: id_PO2 =  22
      integer, parameter :: id_CH3CHO =  23
      integer, parameter :: id_POOH =  24
      integer, parameter :: id_CH3CO3 =  25
      integer, parameter :: id_CH3COOOH =  26
      integer, parameter :: id_PAN =  27
      integer, parameter :: id_ONIT =  28
      integer, parameter :: id_C2H6 =  29
      integer, parameter :: id_C2H4 =  30
      integer, parameter :: id_C4H10 =  31
      integer, parameter :: id_MPAN =  32
      integer, parameter :: id_ISOPO2 =  33
      integer, parameter :: id_MVK =  34
      integer, parameter :: id_MACR =  35
      integer, parameter :: id_MACRO2 =  36
      integer, parameter :: id_MACROOH =  37
      integer, parameter :: id_MCO3 =  38
      integer, parameter :: id_C2H5O2 =  39
      integer, parameter :: id_C2H5OOH =  40
      integer, parameter :: id_C10H16 =  41
      integer, parameter :: id_C3H8 =  42
      integer, parameter :: id_C3H7O2 =  43
      integer, parameter :: id_C3H7OOH =  44
      integer, parameter :: id_CH3COCH3 =  45
      integer, parameter :: id_ROOH =  46
      integer, parameter :: id_CH3OH =  47
      integer, parameter :: id_C2H5OH =  48
      integer, parameter :: id_GLYALD =  49
      integer, parameter :: id_HYAC =  50
      integer, parameter :: id_EO2 =  51
      integer, parameter :: id_EO =  52
      integer, parameter :: id_HYDRALD =  53
      integer, parameter :: id_RO2 =  54
      integer, parameter :: id_CH3COCHO =  55
      integer, parameter :: id_ISOPNO3 =  56
      integer, parameter :: id_ONITR =  57
      integer, parameter :: id_XO2 =  58
      integer, parameter :: id_XOOH =  59
      integer, parameter :: id_ISOPOOH =  60
      integer, parameter :: id_H2 =  61
      integer, parameter :: id_O3S = 62
      integer, parameter :: id_SO2 = 63
      integer, parameter :: id_SO4 = 64
      integer, parameter :: id_DMS = 65
      integer, parameter :: id_NH3 = 66
      integer, parameter :: id_NH4NO3 = 67
      integer, parameter :: id_NH4 = 68
      integer, parameter :: id_HCl = 69
      integer, parameter :: id_HOCl = 70
      integer, parameter :: id_ClONO2 = 71
      integer, parameter :: id_Cl = 72
      integer, parameter :: id_ClO = 73
      integer, parameter :: id_Cl2O2 = 74
      integer, parameter :: id_Cl2 = 75
      integer, parameter :: id_HOBr = 76
      integer, parameter :: id_HBr = 77
      integer, parameter :: id_BrONO2 = 78
      integer, parameter :: id_Br = 79
      integer, parameter :: id_BrO = 80
      integer, parameter :: id_BrCl = 81
      integer, parameter :: id_LCH4 = 82
      integer, parameter :: id_H = 83
      integer, parameter :: id_H2O = 84
  
  
      end module M_SPC_ID_MOD
                                                                                
      module M_RXT_ID_MOD
                                                                                
      implicit none                                                             
                                                                                
      integer, parameter :: rid_jo2 =    1                                      
      integer, parameter :: rid_jo1d =    2                                     
      integer, parameter :: rid_jo3p =    3                                     
      integer, parameter :: rid_jn2o =    4                                     
      integer, parameter :: rid_jno =    5                                      
      integer, parameter :: rid_jno2 =    6                                     
      integer, parameter :: rid_jn2o5 =    7                                    
      integer, parameter :: rid_jhno3 =    8                                    
      integer, parameter :: rid_jno3 =    9                                     
      integer, parameter :: rid_jho2no2 =   10                                  
      integer, parameter :: rid_jch3ooh =   11                                  
      integer, parameter :: rid_jch2o_a =   12                                  
      integer, parameter :: rid_jch2o_b =   13                                  
      integer, parameter :: rid_jh2o =   14                                     
      integer, parameter :: rid_jh2o2 =   15                                    
      integer, parameter :: rid_jch3cho =   16                                  
      integer, parameter :: rid_jpooh =   17                                    
      integer, parameter :: rid_jch3co3h =   18                                 
      integer, parameter :: rid_jpan =   19                                     
      integer, parameter :: rid_jmpan =   20                                    
      integer, parameter :: rid_jmacr_a =   21                                  
      integer, parameter :: rid_jmacr_b =   22                                  
      integer, parameter :: rid_jmvk =   23                                     
      integer, parameter :: rid_jc2h5ooh =   24                                 
      integer, parameter :: rid_jc3h7ooh =   25                                 
      integer, parameter :: rid_jrooh =   26                                    
      integer, parameter :: rid_jacet =   27                                    
      integer, parameter :: rid_jmgly =   28                                    
      integer, parameter :: rid_jxooh =   29                                    
      integer, parameter :: rid_jonitr =   30                                   
      integer, parameter :: rid_jisopooh =   31                                 
      integer, parameter :: rid_jhyac =   32                                    
      integer, parameter :: rid_jglyald =   33                                  
      integer, parameter :: rid_jclono2 =   34                                  
      integer, parameter :: rid_jhocl =   35                                    
      integer, parameter :: rid_jcl2o2 =   36                                   
      integer, parameter :: rid_jbrono2 =   37                                  
      integer, parameter :: rid_jhobr =   38                                    
      integer, parameter :: rid_jbrcl =   39                                    
      integer, parameter :: rid_jbro =   40                                     
      integer, parameter :: rid_jcl2 =   41                                     
      integer, parameter :: rid_usr1 =   42                                     
      integer, parameter :: rid_o1d_n2 =   44                                   
      integer, parameter :: rid_o1d_o2 =   45                                   
      integer, parameter :: rid_ox_l1 =   46                                    
      integer, parameter :: rid_ox_p1 =   49                                    
      integer, parameter :: rid_usr2 =   54                                     
      integer, parameter :: rid_usr3 =   55                                     
      integer, parameter :: rid_usr4 =   57                                     
      integer, parameter :: rid_usr5 =   58                                     
      integer, parameter :: rid_usr6 =   60                                     
      integer, parameter :: rid_usr7 =   62                                     
      integer, parameter :: rid_ox_p2 =   65                                    
      integer, parameter :: rid_usr8 =   72                                     
      integer, parameter :: rid_usr8a =   73                                    
      integer, parameter :: rid_ox_l2 =   77                                    
      integer, parameter :: rid_ox_l3 =   78                                    
      integer, parameter :: rid_usr9 =   79                                     
      integer, parameter :: rid_usr10 =   84                                    
      integer, parameter :: rid_ox_l4 =   85                                    
      integer, parameter :: rid_ox_p3 =   87                                    
      integer, parameter :: rid_ox_p4 =   92                                    
      integer, parameter :: rid_usr11 =   93                                    
      integer, parameter :: rid_usr12 =   97                                    
      integer, parameter :: rid_ox_l5 =   99                                    
      integer, parameter :: rid_ox_p5 =  101                                    
      integer, parameter :: rid_usr13 =  106                                    
      integer, parameter :: rid_ox_l6 =  110                                    
      integer, parameter :: rid_ox_p6 =  113                                    
      integer, parameter :: rid_ox_l7 =  119                                    
      integer, parameter :: rid_ox_l8 =  121                                    
      integer, parameter :: rid_ox_p7 =  122                                    
      integer, parameter :: rid_ox_p8 =  129                                    
      integer, parameter :: rid_usr14 =  135                                    
      integer, parameter :: rid_usr15 =  136                                    
      integer, parameter :: rid_ox_l9 =  138                                    
      integer, parameter :: rid_usr16 =  140                                    
      integer, parameter :: rid_usr17 =  141                                    
      integer, parameter :: rid_ox_p9 =  145                                    
      integer, parameter :: rid_usr22 =  149                                    
      integer, parameter :: rid_ox_p10 =  150                                   
      integer, parameter :: rid_ox_p11 =  164                                   
      integer, parameter :: rid_usr21 =  170                                    
      integer, parameter :: rid_usr24 =  180                                    
      integer, parameter :: rid_usr25 =  182                                    
      integer, parameter :: rid_strat13 =  184                                  
      integer, parameter :: rid_strat14 =  185                                  
      integer, parameter :: rid_strat20 =  186                                  
      integer, parameter :: rid_strat21 =  187                                  
      integer, parameter :: rid_strat22 =  188                                  
      integer, parameter :: rid_strat23 =  189                                  
      integer, parameter :: rid_strat24 =  190                                  
      integer, parameter :: rid_strat25 =  191                                  
      integer, parameter :: rid_strat26 =  192                                  
      integer, parameter :: rid_strat27 =  193                                  
      integer, parameter :: rid_strat28 =  194                                  
      integer, parameter :: rid_strat29 =  195                                  
      integer, parameter :: rid_strat33 =  196                                  
      integer, parameter :: rid_strat35 =  197                                  
      integer, parameter :: rid_strat37 =  198                                  
      integer, parameter :: rid_strat38 =  199                                  
      integer, parameter :: rid_strat39 =  200                                  
      integer, parameter :: rid_strat40 =  201                                  
      integer, parameter :: rid_strat41 =  202                                  
      integer, parameter :: rid_strat42 =  203                                  
      integer, parameter :: rid_strat43 =  204                                  
      integer, parameter :: rid_strat44 =  205                                  
      integer, parameter :: rid_strat45 =  206                                  
      integer, parameter :: rid_strat46 =  207                                  
      integer, parameter :: rid_strat47 =  208                                  
      integer, parameter :: rid_strat48 =  209                                  
      integer, parameter :: rid_strat69 =  210                                  
      integer, parameter :: rid_strat58 =  211                                  
      integer, parameter :: rid_strat59 =  212                                  
      integer, parameter :: rid_strat64 =  213                                  
      integer, parameter :: rid_strat71 =  214                                  
      integer, parameter :: rid_strat72 =  215                                  
      integer, parameter :: rid_strat73 =  216                                  
      integer, parameter :: rid_strat74 =  217                                  
      integer, parameter :: rid_strat75 =  218                                  
      integer, parameter :: rid_strat76 =  219                                  
      integer, parameter :: rid_strat77 =  220                                  
      integer, parameter :: rid_strat78 =  221                                  
      integer, parameter :: rid_strat79 =  222                                  
      integer, parameter :: rid_strat80 =  223                                  
                                                                                
      integer, parameter :: rid_r0043 =   43                                    
      integer, parameter :: rid_r0047 =   47                                    
      integer, parameter :: rid_r0048 =   48                                    
      integer, parameter :: rid_r0050 =   50                                    
      integer, parameter :: rid_r0051 =   51                                    
      integer, parameter :: rid_r0052 =   52                                    
      integer, parameter :: rid_r0053 =   53                                    
      integer, parameter :: rid_r0056 =   56                                    
      integer, parameter :: rid_r0059 =   59                                    
      integer, parameter :: rid_r0061 =   61                                    
      integer, parameter :: rid_r0063 =   63                                    
      integer, parameter :: rid_r0064 =   64                                    
      integer, parameter :: rid_r0066 =   66                                    
      integer, parameter :: rid_r0067 =   67                                    
      integer, parameter :: rid_r0068 =   68                                    
      integer, parameter :: rid_r0069 =   69                                    
      integer, parameter :: rid_r0070 =   70                                    
      integer, parameter :: rid_r0071 =   71                                    
      integer, parameter :: rid_r0074 =   74                                    
      integer, parameter :: rid_r0075 =   75                                    
      integer, parameter :: rid_r0076 =   76                                    
      integer, parameter :: rid_r0080 =   80                                    
      integer, parameter :: rid_r0081 =   81                                    
      integer, parameter :: rid_r0082 =   82                                    
      integer, parameter :: rid_r0083 =   83                                    
      integer, parameter :: rid_r0086 =   86                                    
      integer, parameter :: rid_r0088 =   88                                    
      integer, parameter :: rid_r0089 =   89                                    
      integer, parameter :: rid_r0090 =   90                                    
      integer, parameter :: rid_r0091 =   91                                    
      integer, parameter :: rid_r0094 =   94                                    
      integer, parameter :: rid_r0095 =   95                                    
      integer, parameter :: rid_r0096 =   96                                    
      integer, parameter :: rid_r0098 =   98                                    
      integer, parameter :: rid_r0100 =  100                                    
      integer, parameter :: rid_r0102 =  102                                    
      integer, parameter :: rid_r0103 =  103                                    
      integer, parameter :: rid_r0104 =  104                                    
      integer, parameter :: rid_r0105 =  105                                    
      integer, parameter :: rid_r0107 =  107                                    
      integer, parameter :: rid_r0108 =  108                                    
      integer, parameter :: rid_r0109 =  109                                    
      integer, parameter :: rid_r0111 =  111                                    
      integer, parameter :: rid_r0112 =  112                                    
      integer, parameter :: rid_r0114 =  114                                    
      integer, parameter :: rid_r0115 =  115                                    
      integer, parameter :: rid_r0116 =  116                                    
      integer, parameter :: rid_r0117 =  117                                    
      integer, parameter :: rid_r0118 =  118                                    
      integer, parameter :: rid_r0120 =  120                                    
      integer, parameter :: rid_r0123 =  123                                    
      integer, parameter :: rid_r0124 =  124                                    
      integer, parameter :: rid_r0125 =  125                                    
      integer, parameter :: rid_r0126 =  126                                    
      integer, parameter :: rid_r0127 =  127                                    
      integer, parameter :: rid_r0128 =  128                                    
      integer, parameter :: rid_r0130 =  130                                    
      integer, parameter :: rid_r0131 =  131                                    
      integer, parameter :: rid_r0132 =  132                                    
      integer, parameter :: rid_r0133 =  133                                    
      integer, parameter :: rid_r0134 =  134                                    
      integer, parameter :: rid_r0137 =  137                                    
      integer, parameter :: rid_r0139 =  139                                    
      integer, parameter :: rid_r0142 =  142                                    
      integer, parameter :: rid_r0143 =  143                                    
      integer, parameter :: rid_r0144 =  144                                    
      integer, parameter :: rid_r0146 =  146                                    
      integer, parameter :: rid_r0147 =  147                                    
      integer, parameter :: rid_r0148 =  148                                    
      integer, parameter :: rid_r0151 =  151                                    
      integer, parameter :: rid_r0152 =  152                                    
      integer, parameter :: rid_r0153 =  153                                    
      integer, parameter :: rid_r0154 =  154                                    
      integer, parameter :: rid_r0155 =  155                                    
      integer, parameter :: rid_r0156 =  156                                    
      integer, parameter :: rid_r0157 =  157                                    
      integer, parameter :: rid_r0158 =  158                                    
      integer, parameter :: rid_r0159 =  159                                    
      integer, parameter :: rid_r0160 =  160                                    
      integer, parameter :: rid_r0161 =  161                                    
      integer, parameter :: rid_r0162 =  162                                    
      integer, parameter :: rid_r0163 =  163                                    
      integer, parameter :: rid_r0165 =  165                                    
      integer, parameter :: rid_r0166 =  166                                    
      integer, parameter :: rid_r0167 =  167                                    
      integer, parameter :: rid_r0168 =  168                                    
      integer, parameter :: rid_r0169 =  169                                    
      integer, parameter :: rid_r0171 =  171                                    
      integer, parameter :: rid_r0172 =  172                                    
      integer, parameter :: rid_r0173 =  173                                    
      integer, parameter :: rid_r0174 =  174                                    
      integer, parameter :: rid_r0175 =  175                                    
      integer, parameter :: rid_r0176 =  176                                    
      integer, parameter :: rid_r0177 =  177                                    
      integer, parameter :: rid_r0178 =  178                                    
      integer, parameter :: rid_r0179 =  179                                    
      integer, parameter :: rid_r0181 =  181                                    
      integer, parameter :: rid_r0183 =  183                                    
      integer, parameter :: rid_r0224 =  224                                    
                                                                                
      end module M_RXT_ID_MOD
                                                                                
      module M_HET_ID_MOD
                                                                                
      implicit none                                                             
                                                                                
                                                                                
      end module M_HET_ID_MOD
