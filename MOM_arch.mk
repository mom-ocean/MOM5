#
# System dependent FLAGS for FVdycoreCubed.
#

  FCNAME = $(word 1,$(shell $(FC) --version))

# Intel Fortran Compiler
# ----------------------
ifeq ($(FCNAME),ifort)

  USER_FFLAGS = -stack_temps -safe_cray_ptr -i_dynamic -assume byterecl \
                -vec-report0 -fp-model precise -fp-model source \
                -ftz -w95 -align all -fno-alias -align dcommons

endif

ifeq ($(FCNAME), GNU)  # gfortran

        USER_FFLAGS = -DNO_R16 -DNO_CRAY_POINTERS

endif

ifeq ($(ESMA_FC), ftn)
      USER_FFLAGS = -DNO_R16
endif

