# putting anything here sets it to true.  leave blank for false
NDEBUG := t 
MPI    := t
OMP    :=

COMP := Intel

FPARALLEL := ../../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90

include $(FPARALLEL)/varden/test_sdc/radau/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden/test_sdc/radau

include $(FPARALLEL)/varden/test_sdc/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden/test_sdc

include $(FPARALLEL)/varden/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

