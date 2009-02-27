NDEBUG := t
MPI    := 
OMP    :=

COMP := Intel

FPARALLEL := ../../..

include $(FPARALLEL)/mk/GMakedefs.mak

#f90sources += strang_main.f90
#f90sources += subt_main.f90
f90sources += main.f90
f90sources += checkpoint.f90

include $(FPARALLEL)/varden/test_sdc/convergence/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden/test_sdc/convergence

#include $(FPARALLEL)/varden/checkpoint.f90
VPATH_LOCATIONS += $(FPARALLEL)/varden

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

