NDEBUG :=
NDEBUG :=
MPI    :=
OMP    :=

COMP := Intel

FPARALLEL := ../../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90
fsources += t_dsyev.f
fsources += t_dsygv.f

include $(FPARALLEL)/extern/LAPACK/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/LAPACK

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

