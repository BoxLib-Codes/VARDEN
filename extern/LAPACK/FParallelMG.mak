# FParallelMG.mak
#
# Use this to access fParallel's fortran solvers in Parallel/BoxLib code.
# See iamrlib/run2d/GNUmakefile for an example.
# Note: that you can do a parallel make if you are using this file if you
# are using the MODDEP dependency evaluater that is in Parallel/scripts 
# and is used in the default Parallel/mk/Make.{rules,defs}. Otherwise,
# you can't since the dependencies for the object/module files can not be
# inferred from the names or ordering of thes files.
#
# ------------- in your GNUmakefile --------------
# 
# FBOXLIB_HOME =
# FBOXLIB_HOME = ../../../fParallel 
# ifdef FBOXLIB_HOME
#   include FParallelMG.mak
#   DEFINES += -DMG_USE_FBOXLIB
#   Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
#   Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
# endif
#
# VPATH_LOCATIONS   += . $(Blocs) $(Flocs)
# INCLUDE_LOCATIONS += . $(Blocs) $(Flocs)
#
# ------------------------------------------------

f90EXE_sources += daxpy.f90
