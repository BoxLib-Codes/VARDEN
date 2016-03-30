# A set of useful macros for putting together a VARDEN application.

# check the version number -- this comes from the GNU Make cookbook
NEED := 3.81
OK := $(filter $(NEED),$(firstword $(sort $(MAKE_VERSION) $(NEED))))

ifndef OK
  $(error your version of GNU make is too old.  You need atleast version $(NEED))
endif

ifeq ($(findstring ~, $(BOXLIB_HOME)), ~)
   $(error you cannot include the ~ character in your BOXLIB_HOME variable)
endif

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


ifdef TILING
VPATH_LOCATIONS += $(VARDEN_TOP_DIR)/src_tiled
endif


PROBIN_TEMPLATE := $(VARDEN_TOP_DIR)/src/probin.template

# list of the directories to search for _parameters files
PROBIN_PARAMETER_DIRS = ./ 
PROBIN_PARAMETER_DIRS += $(VARDEN_TOP_DIR)/src

# list of all valid _parameters files for probin
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" 
	@echo " "



include $(VARDEN_TOP_DIR)/src/GPackage.mak
VPATH_LOCATIONS += $(VARDEN_TOP_DIR)/src

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG


main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(BOXLIB_HOME)/Tools/F_scripts/makebuildinfo.py \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --boxlib_home "$(BOXLIB_HOME)" \
           --source_home "$(VARDEN_TOP_DIR)" \
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove 
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean:: 
	$(RM) probin.f90 
	$(RM) build_info.f90

