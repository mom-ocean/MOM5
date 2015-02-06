module tracer_type_mod

implicit none
private

public :: tracer_type
public :: tracer_type_version, tracer_type_tagname

character(len=128) :: tracer_type_version = '$Id: tracer_type.F90,v 11.0 2004/09/28 19:30:05 fms Exp $'
character(len=128) :: tracer_type_tagname = '$Name: tikal $'

type tracer_type
  character(len=32) :: name, numerical_representation, advect_horiz, advect_vert, hole_filling
  real :: robert_coeff
end type

end module tracer_type_mod
