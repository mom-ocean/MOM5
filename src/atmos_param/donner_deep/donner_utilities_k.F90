
!VERSION NUMBER:
!  $Id: donner_utilities_k.F90,v 19.0 2012/01/06 20:08:36 fms Exp $

!module donner_utilities_inter_mod

!#include "donner_utilities_interfaces.h"

!end module donner_utilities_inter_mod



!module donner_types_mod

!#include "donner_types.h"

!end module donner_types_inter_mod


!#####################################################################

subroutine don_u_set_column_integral_k   &
         (nk, x_in, ptop, pbot, int_value, p_in, intgl_in, intgl_out, &
          x_out, ermesg, error)

!--------------------------------------------------------------------
!    subroutine don_u_set_column_integral_k modifies the input
!    profile x_in of size nk with column integral intgl_in to produce an
!    output profile x_out, also of size nk, whose column integral 
!    intgl_out has the specified value int_value. 
!    the integral constraint is satisfied by defining a constant cont-
!    ribution to the column integral between pbot and ptop, which 
!    balances the contribution from the rest of the profile. currently
!    the x_in profile has non-zero values only above ptop; if x_in has
!    non-zero values below ptop, the algorithm must be changed in order
!    to work properly.
!    NOTE: currently pbot is by default the surface pressure. this needs
!    to be generalized.
!    any error message is returned in ermesg.
!    "Verav notes 1/7/04" (available from Leo Donner) explain this 
!    routine in more detail, especially the procedures used to enforce 
!    the integral constraint.
!--------------------------------------------------------------------

implicit none

integer,               intent(in)    ::  nk
real, dimension(nk),   intent(in)    ::  x_in
real,                  intent(in)    ::  ptop, pbot, int_value        
real, dimension(nk+1), intent(in)    ::  p_in
real,                  intent(out)   ::  intgl_in, intgl_out
real, dimension(nk),   intent(out)   ::  x_out
character(len=*),      intent(out)   ::  ermesg
integer,               intent(out)   ::  error

!--------------------------------------------------------------------
!  intent(in) variables:
!  
!     nk        size of input profile
!     x_in      profile on lo-res model grid whose column integral 
!               is adjusted
!     ptop      pressure at top of region where integrand is 
!               adjusted [ Pa ] 
!     pbot      pressure at bottom of region where integrand is
!               adjusted[ Pa ]
!     int_value value desired for column integral of output field
!     p_in      lo-res model interface pressure levels [ Pa ]
!
!  intent(out) variables:
!
!     intgl_in  column integral of input field x_in [ units of x_in ]
!     intgl_out column integral of output field x_out [ units of x_in ]
!     x_out     vertical profile of x_in, after adjustment of values
!               at levels between pbot and ptop to produce column 
!               integral with value of int_value 
!     ermesg    error message, if error is encountered
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

     real    :: int_above !  contribution to column integral from layers
                          !  fully above adjustment region
     real    :: int_needed_below 
                          !  integrand contribution from the adjustment
                          !  layer that is required to balance the 
                          !  remainder of the profile
     real    :: int_below !  contribution to column integral after 
                          !  adjustment that comes from layers fully 
                          !  within the adjustment layer
     integer :: k         !  do-loop index


!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!----------------------------------------------------------------------
!    obtain the total column integrand (column) and the partial integ-
!    rand from layers completely above the adjustment layer (int_above) 
!    of the input quantity (x_in). 
!----------------------------------------------------------------------
      intgl_in = 0.
      int_above = 0.
      do k=1,nk              
        intgl_in = intgl_in + x_in(k)*(p_in(k) - p_in(k+1))
        if (p_in(k) <= ptop)  then
          int_above = int_above + x_in(k)*(p_in(k) - p_in(k+1))
        endif
      end do

!----------------------------------------------------------------------
!    define the value of the integrand needed from the adjustment layer
!    (int_needed_below). it is the negative of the column sum divided 
!    by the total delta p between pbot and ptop.
!----------------------------------------------------------------------
      int_needed_below = int_value - intgl_in/(pbot - ptop)

!----------------------------------------------------------------------
!    begin loop assigning new value to output variable in layers
!    completely or partially in the adjustment layer.
!----------------------------------------------------------------------
      do k=1,nk              

!---------------------------------------------------------------------
!    case of layer being completely in adjustment layer. define the 
!    output variable value as int_needed_below.
!---------------------------------------------------------------------
        if (p_in(k+1) >= ptop) then
          x_out(k) = int_needed_below
         
!---------------------------------------------------------------------
!    case of layer straddling top of adjustment layer.
!---------------------------------------------------------------------
        else if (p_in(k+1) < ptop )  then

!----------------------------------------------------------------------
!    define the amount of the needed adjustment layer value that has 
!    been assigned to the layers fully in the adjustment layer 
!    (int_below).
!----------------------------------------------------------------------
          int_below = int_needed_below*(pbot - p_in(k))

!----------------------------------------------------------------------
!    define the portion of the column sum which must be assigned to 
!    this layer straddling the top of the adjustment layer in order
!    to obtain a value of int_value for the column integral; ie, that 
!    which is needed to balance the contributions from above (int_above)
!    and from below (int_below) this layer.
!----------------------------------------------------------------------
          x_out(k) = (int_value -int_above - int_below)/  &
                                                  (p_in(k) - p_in(k+1)) 

!---------------------------------------------------------------------
!    case of layer completely above adjustment layer. define the output
!    field as equal to the input field at this and all higher levels; 
!    then exit the loop. 
!---------------------------------------------------------------------
          x_out(k+1:) = x_in(k+1:)
          exit
        endif
      end do

!---------------------------------------------------------------------
!    recalculate the column sum of the conservative quantity. it should
!    now be of order machine roundoff. return value to calling routine.
!---------------------------------------------------------------------
      intgl_out = 0.
      do k=1,nk              
        intgl_out = intgl_out + x_out(k)*(p_in(k) - p_in(k+1))
      end do

!----------------------------------------------------------------------


end subroutine don_u_set_column_integral_k

!#####################################################################

subroutine don_u_map_hires_c_to_lores_c_k   &
         (n_gcm, n_clo, x_clo, p_clo, ptop, p_gcm, x_gcm, intgl_clo,  &
          intgl_gcm, ermesg, error)

!--------------------------------------------------------------------
!    subroutine don_u_map_hires_c_to_lores_c_k maps the vertical
!    profile between cloud base and and a specified upper pressure level
!    (ptop) of a variable (x_clo) of size n_clo on the cloud-model 
!    vertical grid (p_clo) to a GCM vertical grid (p_gcm) of size 
!    n_gcm. the output profile is x_gcm.
!    vertical integrals on the cloud grid (intgl_clo) and on the 
!    GCM grid (intgl_gcm) are also returned, so that the profile 
!    mapping may be shown to be conservative.
!    any error message is returned in ermesg.
!    "Verav notes 1/7/04" (available from Leo Donner) explain this 
!    routine in more detail. 
!    The routine can handle grid thicknesses of arbitrary size for
!    both the cloud grid and GCM grid. there is no restriction as
!    to the relative sizes of the grids.
!--------------------------------------------------------------------
!    GCM grid:
!     
!     -------------- p_gcm(n_gcm+1)
!     -------------- x_gcm(n_gcm)
!     -------------- p_gcm(n_gcm-1)
!           ...
!     -------------- p_gcm(k+1)
!     -------------- p_gcm(k)
!     -------------- p_gcm(k-1)
!           ...
!     -------------- p_gcm(2)
!     -------------- x_gcm(1)
!     -------------- p_gcm(1)
!
!     Cloud-Model grid:
!
!     -------------- p_clo(n_clo)=pcloha(n_clo) (cloud top)
!     -------------- p_cloha(n_clo-1)
!     -------------- p_clo(n_clo-1)
!           ...   
!     -------------- p_clo(k2)
!     -------------- p_cloha(k2-1)
!     -------------- p_clo(k2-1)
!           ...
!     -------------- p_clo(2)
!     -------------- p_cloha(1)
!     -------------- p_clo(1)  (cloud base)
!--------------------------------------------------------------------   


implicit none

integer,                  intent(in)  :: n_gcm, n_clo
real, dimension(n_clo),   intent(in)  :: x_clo, p_clo
real,                     intent(in)  :: ptop
real, dimension(n_gcm+1), intent(in)  :: p_gcm
real, dimension(n_gcm),   intent(out) :: x_gcm
real,                     intent(out) :: intgl_clo, intgl_gcm
character(len=*),         intent(out) :: ermesg
integer,                  intent(out) :: error

!--------------------------------------------------------------------
!   intent(in) variables:
!
!        n_gcm        number of layers in GCM 
!        n_clo        number of layers in cloud model
!        x_clo        vertical profile of variable on cloud-model grid
!        p_clo        cloud-model pressure profile
!        ptop         pressure denoting top of region of interest [ Pa ]
!        p_gcm        pressure half levels in GCM [ Pa ]
!
!    intent(out) variables:
!
!        x_gcm        vertical profile of variable on GCM vertical grid
!        intgl_clo    vertical integral of variable on cloud-model grid
!        intgl_gcm    vertical integral of variable on GCM grid
!        ermesg       error message, if error is encountered
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      real, dimension(n_clo)  ::  pcloha
      real                    ::  rint, p_upper, p_lower, ptopa 
      integer                 ::  k2start
      integer                 ::  k, k2

!----------------------------------------------------------------------
!  local variables:
!
!     pcloha           pressures at mid-layers in cloud grid [ Pa ]
!     rint             accumulates pressure-weighted sum of input 
!                      variable over the cloud-model layers which over-
!                      lap a given GCM grid layer
!     p_upper          pressure at upper limit (farthest from ground) 
!                      of the current sub-layer [ Pa ]
!     p_lower          pressure at lower limit (closest to ground) of
!                      the current sub-layer [ Pa ]
!     ptopa            uppermost (farthest from ground) pressure value 
!                      on cloud grid
!     k2start          cloud-model vertical index of lowest cloud-model
!                      layer overlapping the current or next gcm layer
!     k, k2            do-loop indices
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define the lowest k index of the cloud model which overlaps the
!    first gcm layer (k2start). initialize the integral over pressure 
!    of the output field on the gcm grid (intgl_gcm). define the lowest
!    pressure at which the fields will have non-zero values as the 
!    pressure at the top of the cloud model (ptopa). define the 
!    bottommost pressure at which cloud is present (p_lower) as the 
!    pressure at cloud base.
!--------------------------------------------------------------------

      k2start  = 1
      intgl_gcm = 0.
      ptopa = MAX (p_clo(n_clo), ptop)
      p_lower = p_clo(1) 

!--------------------------------------------------------------------
!    define the interface pressures on the cloud grid. cloud base is 
!    defined to be at the midpoint of a cloud layer.
!--------------------------------------------------------------------
      do k2=1,n_clo-1
        pcloha(k2) = 0.5*(p_clo(k2) + p_clo(k2+1))
      end do
      pcloha(n_clo) = p_clo(n_clo)

!--------------------------------------------------------------------
!    march upward through gcm model layers until cloud base is found. 
!    calculate the output profile on the gcm grid. 
!--------------------------------------------------------------------
      do k=1,n_gcm

!--------------------------------------------------------------------
!    check if this gcm layer has any overlap with the cloud. if it does,
!    initialize the output field in this layer to be 0.0.
!--------------------------------------------------------------------
        if (p_gcm(k) > ptopa .and. p_lower > p_gcm(k+1)) then
          rint = 0.

!--------------------------------------------------------------------
!    march through the cloud model layers which overlap this gcm
!    layer.
!--------------------------------------------------------------------
          do k2=k2start,n_clo-1

!---------------------------------------------------------------------
!    define the topmost pressure in this gcm layer for which the 
!    current cloud model layer value will apply (p_upper). it will be 
!    the larger of the cloud top pressure, the pressure at the top of
!    the gcm layer, or the pressure at the top interface of the current
!    cloud layer.
!---------------------------------------------------------------------
            p_upper = MAX (ptopa, p_gcm(k+1), pcloha(k2))

!---------------------------------------------------------------------
!    define the contribution to the current gcm layer from this cloud 
!    layer, weighted by the pressure depth over which it applies.     
!---------------------------------------------------------------------
            rint = rint + x_clo(k2)*(p_lower - p_upper)

!--------------------------------------------------------------------
!    define the pressure at the bottom of the next layer to be 
!    processed as the pressure at the top of the current layer. if 
!    either the next gcm layer or cloud top has been reached, save the 
!    cloud model layer index and exit this loop. otherwise, there are
!    additional cloud model layer(s) contributing to the integral at
!    this GCM level. Cycle through the loop, adding the remaining
!    contributions.
!--------------------------------------------------------------------
            p_lower = p_upper
            if (p_upper /= pcloha(k2)) then
              k2start = k2
              exit
            endif
          end do

!--------------------------------------------------------------------
!    define the output variable at this gcm level and add the contrib-
!    ution to the integral from this layer to the integrand.
!--------------------------------------------------------------------
          x_gcm(k) = rint/(p_gcm(k) - p_gcm(k+1))
          intgl_gcm = intgl_gcm + rint

!---------------------------------------------------------------------
!    if this gcm layer does not overlap the cloud, define its value to
!    be 0.0.
!---------------------------------------------------------------------
        else
          x_gcm(k) = 0.
        endif   
      end do

!--------------------------------------------------------------------
!    evaluate integral over pressure at cloud-model resolution.
!--------------------------------------------------------------------
      p_upper = max (pcloha(1), ptopa)
      intgl_clo = x_clo(1)*(p_clo(1) - p_upper)
      do k2=1,n_clo 
        p_upper = max (pcloha(k2+1), ptopa)
        intgl_clo = intgl_clo + x_clo(k2+1)*(pcloha(k2) - p_upper)
        if (pcloha(k2+1) <= ptopa) exit
      end do

!--------------------------------------------------------------------


end subroutine don_u_map_hires_c_to_lores_c_k



!#####################################################################

subroutine don_u_compare_integrals_k    &
         (hi_res_intgl, lo_res_intgl, diag_unit, ermesg, error)

!---------------------------------------------------------------------
!    subroutine don_u_compare_integrals_k determines if two 
!    input vertical integrals, computed on different grids may be 
!    considered equal, after allowing for the roundoff differences 
!    inherent in their calculation.
!    any error message is returned in ermesg.
!---------------------------------------------------------------------

implicit none

real,               intent(in)  :: hi_res_intgl, lo_res_intgl
integer,            intent(in)  :: diag_unit 
character(len=*),   intent(out) :: ermesg
integer,            intent(out) :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        hi_res_intgl    integral calculated on higher resolution 
!                        vertical grid
!        lo_res_intgl    integral calculated on lower resolution 
!                        vertical grid
!        diag_unit       unit number of column diagnostics file to 
!                        which output message is written
!
!   intent(out) variables:
!
!        ermesg          error message, if error is encountered
!
!----------------------------------------------------------------------
      
!---------------------------------------------------------------------
!   local variables:

      integer    :: inteq  ! flag indicating status of integral equality

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!---------------------------------------------------------------------
!    call don_u_integrals_are_equal_k to determine the 
!    equality of hi_res_intgl and lo_res_intgl. the return flag inteq
!    will be 0 if the integrals may be considered equal. if it is 
!    non-zero, then it may indicate that the mapping of vertical 
!    integrals between the hi- and lo-resolution grids should be exam-
!    ined, or that the roundoff in the offending integral may simply be 
!    slightly larger than the current tolerance, perhaps due to the 
!    nature of the integral.
!---------------------------------------------------------------------
      call don_u_integrals_are_equal_k    &
                      (hi_res_intgl, lo_res_intgl, ermesg, error, inteq)
 
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, return to calling program where it will be processed.
!----------------------------------------------------------------------
      if (error /= 0 ) return

!-----------------------------------------------------------------------
!    output a warning message to the column diagnostics output file if
!    the integral equality test produces suspicious results.
!-----------------------------------------------------------------------
      if (inteq /= 0) then
        write (diag_unit, '(a)')  &
           'WARNING: cloud model and ls model intgls differ &
                     &non-trivially  -- perhaps significantally ?'
      endif

!---------------------------------------------------------------------


end subroutine don_u_compare_integrals_k 


!######################################################################

subroutine don_u_apply_integral_source_k      &
         (nk, x_in, ptop, pbot, src, p_in, i_in, i_out, x_out, ermesg, error)

!----------------------------------------------------------------------
!    subroutine don_u_apply_integral_source_k adds a specified
!    value src to the input field x_in of size nk on a pressure grid 
!    p_in between pressure levels pbot and ptop, resulting in a change 
!    of column integral of the field from i_in to i_out, and producing 
!    the output field x_out.
!    NOTE: currently pbot is by default the surface pressure. this needs
!    to be generalized.
!    any error message is returned in ermesg.
!----------------------------------------------------------------------

implicit none

integer,               intent(in)    :: nk
real, dimension(nk),   intent(in)    :: x_in
real,                  intent(in)    :: ptop, pbot, src
real, dimension(nk+1), intent(in)    :: p_in
real,                  intent(out)   :: i_in, i_out
real, dimension(nk),   intent(out)   :: x_out
character(len=*),      intent(out)   :: ermesg
integer,               intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      nk        size of input profile
!      x_in      variable to which the integral source src is to
!                be applied [ units of x_in ]
!      ptop      topmost pressure at which src is applied [ Pa ] 
!      pbot      bottommost pressure at which src is applied [ Pa ] 
!      src       integral source to be applied to variable x_in
!      p_in      interface pressure levels of grid of x_in [ Pa ]
!
!   intent(out) variables:
!
!      i_in      column integral of input variable x_in 
!                [ units of x_in * Pa ]
!      i_out     column integral of output variable x_out
!                [ units of x_in * Pa ]
!      x_out     variable field after the source integral is applied
!                [ units of x_in ]
!      ermesg    error message, if error is encountered
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
      integer   :: k      ! do_loop index

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    compute the column integral (i_in) of the input field x_in.
!---------------------------------------------------------------------
      i_in = 0.
      do k=1,nk              
        i_in = i_in + x_in(k)*(p_in(k) - p_in(k+1))
      end do

!----------------------------------------------------------------------
!    apply the integral source src to each layer in the specified 
!    region. for layers fully included, add src to the existing value. 
!    for the layer containing the top of the specified region, add the 
!    value of src appropriately pressure-weighted. NOTE: this source WAS
!    applied ONLY to the layer straddling the top of the region; it was
!    not applied below cloud base. I HAVE MODIFIED THIS CODE SO THAT THE
!    VALUE IS APPLIED FROM SFC TO TOP OF SPECIFIED REGION (CLOUD BASE)
!    IS THIS CORRECT AND WHAT WAS INTENDED ?? 
!----------------------------------------------------------------------
      do k=1,nk              
        if (p_in(k+1) >= ptop ) then                   
          x_out(k) = x_in(k) + src
        else if (p_in(k+1) < ptop .and. p_in(k) > ptop)  then
          x_out(k) = x_in(k) + (src/(p_in(k) - p_in(k+1))*  &
                                           (p_in(k) - ptop))
          x_out(k+1:) = x_in(k+1:)
          exit
        endif
      end do

!---------------------------------------------------------------------
!    compute the column integral (i_out) of the output field x_out.
!---------------------------------------------------------------------
      i_out = 0.
      do k=1,nk              
        i_out = i_out + x_out(k)*(p_in(k) - p_in(k+1))
      end do

!--------------------------------------------------------------------



end subroutine don_u_apply_integral_source_k


!#####################################################################

subroutine don_u_integrals_are_equal_k    &
         (x, y, ermesg, error, inteq)

!--------------------------------------------------------------------
!    subroutine don_u_integrals_are_equal_k determines if two 
!    integrals x and y are within a roundoff tolerance eps_i of one 
!    another. it computes the difference x - y, and returns inteq, which
!    is given a value of -10 if (x - y) < -eps_i, a value of 10 if 
!    (x - y) > eps_i, and a value of 0 if ABS (x - y) < eps_i. 
!    any error message is returned in ermesg.
!--------------------------------------------------------------------

implicit none

real,               intent(in)   :: x, y
character(len=*),   intent(out)  :: ermesg
integer,            intent(out)  :: error, inteq

!---------------------------------------------------------------------
!   intent (in) variables:
!
!          x        first variable
!          y        second variable
!
!   intent(out) variables:
!
!          ermesg    error message, if error is encountered
!          inteq      ==   0 if x = y, within a tolerance of eps_i
!                     ==  10 if x > y by at least eps_i
!                     == -10 if x < y by at least eps_i
!
!---------------------------------------------------------------------
      
real, PARAMETER  :: eps_i = 1.0e-12 ! roundoff tolerance for equality of
                                    ! two vertical integrals; if the 
                                    ! magnitude of the numerical differ-
                                    ! ence between two integrals is 
                                    ! smaller than eps_i, they are 
                                    ! assumed equal by function 
                                    ! integrals_are_equal.

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    if the integral values are "large", then compare the ratio of
!    integral difference to integral value to the tolerance eps_i, 
!    since the difference itself could be much larger than the toler-
!    ance and still be roundoff. 
!--------------------------------------------------------------------
      if (ABS(x) > 1.0) then

!--------------------------------------------------------------------
!    define integrals_are_equal dependent on the relationship between 
!    (x - y) and eps_i.
!--------------------------------------------------------------------
        if ( (x - y)/ ABS(x) > eps_i) then
          inteq = 10
        else if ( (x - y)/ABS(x) < -eps_i) then
          inteq = -10
        else
          inteq = 0
        endif

!--------------------------------------------------------------------
!    if the integral values are "small", simply compare their difference
!    to the tolerance, since as the integral magnitudes approach zero,
!    the ratio of difference / magnitude could become quite large and 
!    still only represent roundoff error. 
!--------------------------------------------------------------------
      else

!--------------------------------------------------------------------
!    define integrals_are_equal dependent on the relationship between 
!    (x - y) and eps_i.
!--------------------------------------------------------------------
        if ( (x - y) > eps_i) then
          inteq = 10
        else if ( (x - y) < -eps_i) then
          inteq = -10
        else
          inteq = 0
        endif
      endif

!-------------------------------------------------------------------


end subroutine don_u_integrals_are_equal_k





!#####################################################################



!######################################################################

subroutine don_u_map_hires_i_to_lores_c_k    &
         (n_lo, intgl_hi, pbot, ptop, p_lo, x_lo, ermesg, error)

!------------------------------------------------------------------
!    subroutine don_u_map_hires_i_to_lores_c_k assigns
!    an integral quantity defined on a high-res grid (intgl_hi), valid 
!    over a specified pressure depth (pbot, ptop), to the elements of an 
!    array (x_lo) of size n_lo on a lower resolution grid (p_lo).
!    any error message is returned in ermesg.
!------------------------------------------------------------------

implicit none

integer,                  intent(in)     :: n_lo
real,                     intent(in)     :: intgl_hi, pbot, ptop
real, dimension(n_lo+1),  intent(in)     :: p_lo
real, dimension(n_lo),    intent(out)    :: x_lo
character(len=*),         intent(out)    :: ermesg
integer,                  intent(out)    :: error

!------------------------------------------------------------------
!   intent(in) variables:
!
!       n_lo       size of output profile x_lo
!       intgl_hi   input integral from hi-res model
!                  [ units of intgl_hi ]
!       pbot       pressure at bottom of integral's extent [ Pa ]
!       ptop       pressure at top of integral's extent  [ Pa ]
!       p_lo       lo-res model interface pressure levels [ Pa ]
!
!   intent(out) variables:
!
!       x_lo       values of integrand at lo-res model levels
!                  [ units of intgl_hi ]
!       ermesg     error message, if error is encountered
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!   local variables:
!
      integer :: k    ! do-loop index

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!---------------------------------------------------------------------
!    verify that the specified pressure limits are reasonable. if not,
!    define an error message and return.
!---------------------------------------------------------------------
      if ( pbot < ptop ) then
        ermesg = ' input pressure pbot is less than input pressure ptop'
        error = 1
        return
      endif

!---------------------------------------------------------------------
!    assign the proper value to each large-scale model level.
!---------------------------------------------------------------------
      do k=1,n_lo            

!---------------------------------------------------------------------
!    case of bottom of lo-res model layer being within the integral's 
!    specified region.
!---------------------------------------------------------------------
        if (pbot > p_lo(k))  then 

!---------------------------------------------------------------------
!    case of top of lo-res model layer being above the integral's
!    specfied topmost pressure. in such case all lo-res model levels 
!    above are assigned values of 0.0 and the loop exited.
!---------------------------------------------------------------------
          if (ptop >= p_lo(k))  then
            x_lo(k:) = 0.
            exit

!---------------------------------------------------------------------
!    case of lo-res model layer being completely within the integral's
!    specified region.
!---------------------------------------------------------------------
          else if ( ptop <  p_lo(k+1) )  then 
            x_lo(k) = intgl_hi

!---------------------------------------------------------------------
!    case of lo-res model layer extending above the top of the integ-
!    ral's specified region. in such case, only the portion of the
!    integral contained within the appropriate pressure fraction of
!    the lo-res model layer is assigned to the output field.
!---------------------------------------------------------------------
          else
            x_lo(k) = intgl_hi*(p_lo(k) - ptop)/(p_lo(k) - p_lo(k+1))
          endif

!---------------------------------------------------------------------
!    case of bottom of lo-res model layer being below the integral's
!    specified region.
!---------------------------------------------------------------------
        else

!---------------------------------------------------------------------
!    case of top of lo-res model layer being below the integral's
!    specified region. in this case, lo-res model layer is completely
!    outside the integral's range and is assigned a value of 0.0.
!---------------------------------------------------------------------
          if (pbot <= p_lo(k+1)) then
            x_lo(k) = 0.

!---------------------------------------------------------------------
!    case of bottom of lo-res model layer being below and top of lo-res
!    model layer being above integral's specified region. in such case,
!    only the portion of the integral contained within the appropriate 
!    pressure fraction of the lo-res model layer is assigned to the 
!    output field. values at levels above this are assigned values of
!    0.0 and the loop exited.
!---------------------------------------------------------------------
          else if (ptop >  p_lo(k+1))  then
            x_lo(k) = intgl_hi*(pbot - ptop)/(p_lo(k) - p_lo(k+1))
            x_lo(k+1:) = 0.
            exit

!---------------------------------------------------------------------
!    case of bottom of lo-res model layer being below bottom of 
!    integral's specified region and top of lo-res model layer being
!    below top of integral's specified region. in such case, only the 
!    portion of the integral contained within the appropriate pressure 
!    fraction of the lo-res model layer is assigned to the output field.
!---------------------------------------------------------------------
          else if (ptop <= p_lo(k+1) )  then
            x_lo(k) = intgl_hi*(pbot - p_lo(k+1))/(p_lo(k) - p_lo(k+1))
          endif
        endif
      end do

!--------------------------------------------------------------------



end subroutine don_u_map_hires_i_to_lores_c_k



!######################################################################


subroutine don_u_numbers_are_equal_k                    &
         (x, y, ermesg, error, numeq)

!--------------------------------------------------------------------
!    subroutine don_u_numbers_are_equal_k determines if two 
!    numbers x and y are within a roundoff tolerance eps_n of one 
!    another. it computes the difference x - y, and returns a value numeq
!    which is -10 if (x - y) < -eps_n, 10 if (x - y) > eps_n, and 0 if 
!    ABS (x - y) < eps_n. 
!    any error message is returned in ermesg.
!--------------------------------------------------------------------

implicit none

real,               intent(in)  :: x, y
character(len=*),   intent(out) :: ermesg
integer,            intent(out) :: error, numeq

!---------------------------------------------------------------------
!   intent (in) variables:
!
!          x       first variable
!          y       second variable
!
!   intent (out) variables:
!
!          ermesg  error message, if error is encountered
!          numeq   ==   0 if x = y, within a tolerance of eps_n
!                  ==  10 if x > y by at least eps_n
!                  == -10 if x < y by at least eps_n
!
!---------------------------------------------------------------------
      
real, PARAMETER  :: eps_n = 1.0e-13 ! roundoff tolerance for equality 
                                    ! of two real numbers; if the mag-
                                    ! nitude of the numerical difference
                                    ! between two numbers is smaller than
                                    ! eps_n, they are assumed equal by 
                                    ! function numbers_are_equal.


!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = '  ' ; error = 0

!--------------------------------------------------------------------
!    define numbers_are_equal based on the relationship between 
!    (x - y) and eps_n.
!--------------------------------------------------------------------
      if ( (x - y) > eps_n) then
        numeq = 10
      else if ( (x - y) < -eps_n) then
        numeq = -10
      else
        numeq = 0
      endif

!-------------------------------------------------------------------


end subroutine don_u_numbers_are_equal_k


!#####################################################################




!interface donner_utilities_map_lo_res_col_to_hi_res_col_k

!     subroutine donner_utilities_lo1d_to_hi1d_k     
!     subroutine donner_utilities_lo1d_to_hi0d_linear_k 
!     subroutine donner_utilities_lo1d_to_hi0d_log_k    

!end interface donner_utilities_map_lo_res_col_to_hi_res_col_k

subroutine don_u_lo1d_to_hi1d_k     &
         (n_lo, n_hi, x_lo, p_lo, p_hi, x_hi, ermesg, error)                

!------------------------------------------------------------------
!    subroutine don_u_lo1d_to_hi1d_k interpolates the input 
!    field x_lo of size n_lo on pressure grid p_lo to produce an output 
!    field x_hi of size n_hi on a pressure grid p_hi. 
!    NOTE: vertical index 1 is closest to the ground in all arrays
!    used here.
!    any error message is returned in ermesg.
!------------------------------------------------------------------
 
implicit none

integer,                  intent(in)    :: n_lo, n_hi
real, dimension(n_lo),    intent(in)    :: x_lo
real, dimension(n_lo+1),  intent(in)    :: p_lo
real, dimension(n_hi),    intent(in)    :: p_hi
real, dimension(n_hi),    intent(out)   :: x_hi
character(len=*),         intent(out)   :: ermesg
integer,                  intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     n_lo          size of lo-res profile
!     n_hi          size of hi-res profile
!     x_lo          field to be interpolated on lo-res pressure grid
!                   [ units of x_lo ]
!     p_lo          lo-res pressure grid [ Pa ]
!     p_hi          hi-res pressure grid [ Pa ]
!     
!  intent(out) variables:
!
!     x_hi          value of field at pressure p_hi
!                   [ units of x_lo ]
!     ermesg        error message, if error is encountered
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer :: kkstart            !  k index at which to start search
                                    !  for input field value
      integer :: k, kk              !  do-loop indices

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define the k index in the low resolution grid (index 1 nearest
!    surface) at which to begin searching for the cape model pressure 
!    value.
!--------------------------------------------------------------------
      kkstart = 1

!-------------------------------------------------------------------
!    for each pressure level of the high resolution grid, find the low 
!    resolution pressure levels that bracket it. when found, use linear
!    interpolation to define the field value at the high resolution grid
!    pressure level. if it is beyond either end of the low resolution
!    grid, obtain a value at the high resolution pressure value by 
!    extrapolating the gradient at the low level grid boundary.
!----------------------------------------------------------------------
      do k=1,n_hi           

!--------------------------------------------------------------------
!    if the requested hi-res model pressure is greater than the lowest 
!    lo-res model pressure, obtain the field value at the hi-res pres-
!    sure by extrapolating the lo-res model gradient.
!---------------------------------------------------------------------
        if (p_hi(k) > p_lo(1)) then
          x_hi(k) = x_lo(1) + (p_hi(k) - p_lo(1))* &
                          ((x_lo(2) - x_lo(1))/(p_lo(2) - p_lo(1)))

!---------------------------------------------------------------------
!    if the requested hi-res model pressure is less than the highest 
!    lo-res model pressure, obtain the field value at the hi-res pres-
!    sure by extrapolating the lo-res model gradient.
!---------------------------------------------------------------------
        else if (p_hi(k) < p_lo(n_lo)) then
            x_hi(k) =  x_lo(n_lo) + (p_hi(k) - p_lo(n_lo))*   &
                            ((x_lo(n_lo) - x_lo(n_lo-1))/  &
                             (p_lo(n_lo) - p_lo(n_lo-1)))

!---------------------------------------------------------------------
!    if the requested hi-res model pressure level lies within the bounds
!    of the lo-res model pressure profile, march through the  lo-res
!    profile until the desired hi-res model pressure is reached. define
!    the hi-res field value by interpolating to the hi-res model pres-
!    sure. define the lo-res model starting level index to be used in
!    the search for the next desired hi-res pressure (kkstart), and 
!    exit the vertical loop.
!---------------------------------------------------------------------
        else
          do kk=kkstart,n_lo-1
            if (p_hi(k) >= p_lo(kk+1)) then
              x_hi(k) = x_lo(kk+1) + (p_hi(k) - p_lo(kk+1))* &
                             ((x_lo(kk+1) - x_lo(kk))/  &
                              (p_lo(kk+1) - p_lo(kk)))
              kkstart = kk
              exit
            endif
          end do
        endif
      end do

!---------------------------------------------------------------------


end subroutine don_u_lo1d_to_hi1d_k  


!#####################################################################

subroutine don_u_lo1d_to_hi0d_linear_k  &
         (n_lo, x_lo, p_lo, p_hi, x_hi, ermesg, error)


!--------------------------------------------------------------------
!    subroutine lo1d_to_hi0d linearly interpolates within the 1d input 
!    field x_lo on 1d pressure grid p_lo to determine a scalar output 
!    variable x_hi at pressure p_hi. 
!    NOTE: vertical index 1 is closest to the ground in all arrays
!    used here.
!------------------------------------------------------------------

integer,                 intent(in)    :: n_lo
real, dimension(n_lo),   intent(in)    :: x_lo
real, dimension(n_lo+1), intent(in)    ::  p_lo
real,                    intent(in)    :: p_hi
real,                    intent(out)   :: x_hi
character(len=*),        intent(out)   :: ermesg
integer,                 intent(out)   :: error

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     x_lo          field to be interpolated on lo-res pressure grid
!                   [ units of x_lo ]
!     p_lo          lo-res pressure grid [ Pa ]
!     p_hi          hi-res pressure grid [ Pa ]
!     
!  intent(out) variables:
!
!     x_hi          value of field at pressure p_hi
!                   [ units of x_lo ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (1)          :: x_hi_1d, p_hi_1d
      integer                      :: n_hi


!--------------------------------------------------------------------
!   local variables:
!    
!          x_hi_1d    1d array containing the output field x_hi
!          p_hi_1d    1d array containing the input field p_hi
!          ermesg     character string containing any error message
!                     generated in the kernel subroutines accessed from
!                     here
!          n_lo       vertical size of lo_res model profile           
!          n_hi       vertical size of array containing output profile  
!
!---------------------------------------------------------------------

      ermesg = ' ' ; error = 0

!----------------------------------------------------------------------
!    define the dimensions of input and output arrays. 
!----------------------------------------------------------------------
      n_hi = 1

!---------------------------------------------------------------------
!    define an array containing the hi-res pressure at which a variable
!    value is desired.
!---------------------------------------------------------------------
      p_hi_1d(1) = p_hi

!--------------------------------------------------------------------
!    call don_u_lo1d_to_hi1d_k to obtain the desired field 
!    value.
!--------------------------------------------------------------------
      call don_u_lo1d_to_hi1d_k     & 
           (n_lo, n_hi, x_lo, p_lo, p_hi_1d, x_hi_1d, ermesg, error)

!---------------------------------------------------------------------
!    check to be sure no errors were encountered in the kernel routine. 
!---------------------------------------------------------------------
      if (error /= 0 ) return

!--------------------------------------------------------------------
!    move the returned value from the output array to the scalar output 
!    argument to be returned to the calling routine.
!--------------------------------------------------------------------
       x_hi = x_hi_1d(1)        

!--------------------------------------------------------------------


end subroutine don_u_lo1d_to_hi0d_linear_k


!#####################################################################

subroutine don_u_lo1d_to_hi0d_log_k     &
         (n_lo, x_lo, sig_lo, ps, p_hi, x_hi, ermesg, error)

!--------------------------------------------------------------------
!    subroutine don_u_lo1d_to_hi0d_log_k uses logarithmic 
!    interpolation to define the scalar output variable x_hi at pressure
!    p_hi from the input profile x_lo of size n_lo on sigma grid sig_lo 
!    (with associated surface pressure ps).
!    NOTE: vertical index 1 is closest to the ground in all arrays
!    used here.
!    any error message is returned in ermesg.
!------------------------------------------------------------------

implicit none

integer,               intent(in)   :: n_lo
real, dimension(n_lo), intent(in)   :: x_lo, sig_lo 
real,                  intent(in)   :: ps
real,                  intent(in)   :: p_hi
real,                  intent(out)  :: x_hi     
character(len=*),      intent(out)  :: ermesg
integer,               intent(out)  :: error

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     n_lo     size of input profile on lo-res grid
!     x_lo     field to be interpolated [ units of x_lo ]
!     sig_lo   sigma profile defining the specified grid [ fraction ]
!     ps       surface pressure defining the specified grid [ Pa ]
!     p_hi     pressure whose location on the specified grid is desired 
!              [ Pa ]
!
!   intent(out) variables:
!
!     x_hi     output variable value at pressure p_hi [ units of x_lo ] 
!     ermesg   error message, if error is encountered
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer  :: indx  ! nearest vertical grid index on specified grid 
                        ! surfaceward of the desired pressure
      real     :: displ ! logarithmic displacement of desired pressure 
                        ! from the indx grid level
      integer  :: k     ! do-loop index

!----------------------------------------------------------------------
!    initialize the error message character string.
!----------------------------------------------------------------------
      ermesg = ' ' ; error = 0

!--------------------------------------------------------------------
!    define an initial value for indx which can be checked to verify
!    if the search was successful.
!--------------------------------------------------------------------
      indx = 0

!--------------------------------------------------------------------
!    march through the vertical grid until the desired pressure is 
!    bracketed. define the lower index (indx) and calculate the logar-
!    ithmic displacement of the desired pressure between the two sur-
!    rounding grid levels.
!--------------------------------------------------------------------
      do k=1,n_lo-1
        if ((ps*sig_lo(k) >= p_hi) .and.   &
            (ps*sig_lo(k+1) <= p_hi) ) then 
          indx = k
          displ = alog(p_hi/(ps*sig_lo(k)))/alog(sig_lo(k+1)/sig_lo(k)) 
          x_hi = x_lo(indx) + (x_lo(indx+1) - x_lo(indx))*displ
          exit
        endif
      end do

!---------------------------------------------------------------------
!    if pressure was outside the limits of the input profile, write 
!    error message.
!---------------------------------------------------------------------
      if (indx == 0) then
        ermesg = 'unable to bracket the input pressure within the input &
                                           &pressure profile'
        error = 1
        return
      endif

!--------------------------------------------------------------------


end subroutine don_u_lo1d_to_hi0d_log_k

!#####################################################################

subroutine don_u_process_monitor_k (variable, i,j,nlev_lsm, Monitor)
 
use donner_types_mod, only : donner_monitor_type, MAXVAL, MAXMAG, &
                             MINMAG, MINVAL

implicit none

real, dimension(nlev_lsm),  intent(in)    :: variable
integer,                    intent(in)    :: i, j, nlev_lsm
type(donner_monitor_type),  intent(inout) :: Monitor

 
      integer :: k

      select case (Monitor%limit_type)
        case (MAXMAG)              
          do k=1,nlev_lsm 
            if (abs(variable(k)) > Monitor%threshold) then
              Monitor%hits(i,j,k) = Monitor%hits(i,j,k) + 1.0
            endif 
            if (abs(variable(k)) > Monitor%extrema(i,j,k))  then
              Monitor%extrema(i,j,k) = abs(variable(k)) 
            endif 
          end do
        case (MINMAG)              
          do k=1,nlev_lsm 
            if (abs(variable(k)) < Monitor%threshold) then
              Monitor%hits(i,j,k) = Monitor%hits(i,j,k) + 1.0
            endif 
            if (abs(variable(k)) < Monitor%extrema(i,j,k))  then
              Monitor%extrema(i,j,k) = abs(variable(k)) 
            endif 
          end do
        case (MAXVAL)              
          do k=1,nlev_lsm 
            if (variable(k) > Monitor%threshold) then
              Monitor%hits(i,j,k) = Monitor%hits(i,j,k) + 1.0
            endif 
            if (variable(k) > Monitor%extrema(i,j,k))  then
              Monitor%extrema(i,j,k) = variable(k) 
            endif 
          end do
        case (MINVAL)              
          do k=1,nlev_lsm 
            if (variable(k) < Monitor%threshold) then
              Monitor%hits(i,j,k) = Monitor%hits(i,j,k) + 1.0
            endif 
            if (variable(k) < Monitor%extrema(i,j,k))  then
              Monitor%extrema(i,j,k) = variable(k)
            endif 
          end do
      end select

end subroutine don_u_process_monitor_k


!######################################################################

subroutine find_tropopause (nlev, temp, pfull, ptrop, itrop)

integer,                intent(in) :: nlev
real, dimension (nlev), intent(in) :: temp, pfull
real,                   intent(out) :: ptrop
integer,                intent(out) :: itrop

!---------------------------------------------------------------------
!    find lowest temperature in the column between 50 and 400 hPa. the 
!    corresponding pressure is returned as the tropopause pressure.
!    temp and pfull have index 1 closest to the ground.
!---------------------------------------------------------------------

     real      :: ttrop
     integer   :: k


     ttrop = 400.0
     itrop = 1
     ptrop = 1.0

     do k=nlev,1,-1  
       if (pfull(k) < 400.e+02) then
         if (pfull(k) > 50.e+02 .and. temp(k) < ttrop) then 
           ttrop = temp(k)
           itrop = k
           ptrop = pfull(k)
         endif
       else
         exit
       endif
     end do

!---------------------------------------------------------------------
     
end subroutine find_tropopause 


!######################################################################

subroutine define_arat_erat (option, kpar, eratb, erat0, erat_min, &
                             erat_max, erat, arat)

integer, intent(in) :: option, kpar
real, dimension(kpar), intent(in) :: eratb
real, intent(in) :: erat0, erat_min, erat_max
real, dimension(kpar), intent(out) :: erat, arat

     integer :: i
     real    :: x 
     real, dimension(kpar) :: eratb_loc

     if (option == 1) then
 
       do i=1,kpar
! eq (5):
         arat(i) = (exp(-eratb(i)/erat0) - exp(-eratb(i+1)/erat0))  /  &
                   (exp(-eratb(1)/erat0) - exp(-eratb(kpar+1)/erat0)) 
! eq (6):
         erat(i) = 0.5*(eratb(i+1) + eratb(i))   
       end do

     else if (option == 2) then
 
       eratb_loc(1) = erat_min
       do i=2,kpar+1
! eq (7):
         x = exp(-eratb_loc(i-1)/erat0) -     &
              (exp(-erat_min/erat0) - exp(-erat_max/erat0))/kpar
         eratb_loc(i) = -erat0*log(x)    
       end do
       do i=1,kpar
! eq (5):
         arat(i) =    &
            (exp(-eratb_loc(i)/erat0) - exp(-eratb_loc(i+1)/erat0))/ &
            (exp(-eratb_loc(1)/erat0) - exp(-eratb_loc(kpar+1)/erat0)) 
! eq (6):
         erat(i) = 0.5*(eratb_loc(i+1) + eratb_loc(i)) 
       end do
 
     end if



end subroutine define_arat_erat 


!######################################################################
