#
# Makefile for ESMA components.
#
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
# 07Aug2003  Sawyer    Modified for FVdycore
# 27Aug2003  Sawyer    Added a2d3d and d2a3d
# 22Sep2004  Sawyer    Modifications for merge with current code
# 06Oct2004  Sawyer    Removed spmd_dyn
# 17Feb2005  Sawyer    Added epvd and ppme
# 17May2005  Sawyer    Added FVdycore_wrapper
# 18oct2005  da Silva  Moved ALTIX specific flags to FVdycore.mk
# 18Jan2006  Putman    Added mfz_comp
# 17May2006  Sawyer    Added diag_module.F90
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies
  include                    MOM_arch.mk  # arch dependent flags 

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS := $(shell basename `pwd`)
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
#	$(CP) -p *.rc   $(ESMAETC)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] 

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

#MSRCS := $(shell cat path_names)
#SRCS := $(addprefix src/, $(MSRCS))
SRCS := $(shell cat list_of_files)
OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS)))) 
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

DEPS := $(filter-out nsclock.d threadloc.d, $(DEPS)) # filter these out

SRC_DIRS = $(sort $(dir $(SRCS)))
INC_DIRS = . src/mom5/ocean_param/gotm-4.0/include $(SRC_DIRS) $(INC_GFDL_FMS) \
            $(INC_GMAO_SHARED) $(INC_SDF) $(INC_ESMF) $(INC_MPI)
MOD_DIRS = . $(INC_DIRS) 

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS) /usr/include

USER_FFLAGS += $(BIG_ENDIAN) 
USER_FDEFS  = $(D)MAPL_MODE $(D)EIGHT_BYTE $(D)SPMD $(D)TIMING \
              $(D)use_libMPI $(D)use_netCDF $(D)USE_OCEAN_BGC
USER_CDEFS  = $(USER_FDEFS)
USER_FFLAGS += $(FREAL8) 

FREAL = $(FREAL8)
THIS_GFDL_FMS = GFDL_fms_r8

ifneq ( $(wildcard ../mom4_odas_iau.F90), $(null) ) 
    FOPT += -fpe0
endif


$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


#                  --------------------
#                      Dependencies
#                  --------------------

ocmip2_biotic.o: ocmip2_biotic.F90
	$(FC) -c $(F90FLAGS) -O1 $<

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros



