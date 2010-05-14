NDEBUG := t
MPI    := t
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

#LDFLAGS += -L ../../extern/hypre/lib

#FINCLUDE_LOCATIONS += ../../extern/hypre/include

#hypre_libraries = -lHYPRE

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

