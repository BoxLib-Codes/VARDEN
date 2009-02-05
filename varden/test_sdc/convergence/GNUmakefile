NDEBUG := t
MPI    := 
OMP    :=

COMP := Intel

FPARALLEL := ../../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/varden/C_GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden

include $(FPARALLEL)/varden/test_sdc/convergence/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden/test_sdc/convergence


main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

