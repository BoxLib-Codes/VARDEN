
NDEBUG := t
HYPRE  := 

MPI    := t
OMP    :=

COMP   := gfortran

TOP    := ../..

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

f90sources += main.f90

include $(TOP)/src/GPackage.mak
VPATH_LOCATIONS += $(TOP)/src

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

