# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies
  include       ../../../../MOM_arch.mk  # arch dependent flags 

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

esma_install install: ocean_grid_generator

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] 

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := grids_util.f90 grids_type.f90 hgrid.f90 check_mask.F90 topog.f90 vgrid.f90 ocean_grid_generator.f90
OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS)))) 
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

SRC_DIRS = $(sort $(dir $(SRCS)))
INC_DIRS = $(SRC_DIRS) $(INC_GFDL_FMS) \
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

ocean_grid_generator : $(DEPS) $(OBJS)
	$(FC) -o ocean_grid_generator $(F90FLAGS) $(OBJS) $(ESMALIB)/libGFDL_fms.a $(LIB_NETCDF)


# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros



