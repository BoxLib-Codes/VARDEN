NDEBUG := t
MPI    := t
HYPRE  := 

OMP    :=

COMP := gfortran

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90

include $(FPARALLEL)/varden/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/source/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg/source

include $(FPARALLEL)/extern/LAPACK/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/LAPACK

ifdef HYPRE

# on manda
  HYPRE_DIR = ../../extern/hypre

# on hopper
# HYPRE_DIR =  /usr/common/acts/Hypre/2.6.0b/src/cray-xt5_O

  LDFLAGS += -L $(HYPRE_DIR)/lib
  FINCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  hypre_libraries = -lHYPRE

endif

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

