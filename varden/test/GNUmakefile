NDEBUG :=
MPI    :=
OMP    :=

COMP := Intel

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/varden/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

