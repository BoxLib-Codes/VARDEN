TOP = ..\..

!INCLUDE $(TOP)\mk\NMakedefs.mak

boxlib_dir=$(TOP)\boxlib
!INCLUDE $(boxlib_dir)\NPackage.mak
mg_dir=$(TOP)\mg
!INCLUDE $(mg_dir)\NPackage.mak
sparskit_dir=$(TOP)\extern\SPARSKIT
!INCLUDE $(sparskit_dir)\NPackage.mak
lapack_dir=$(TOP)\extern\LAPACK
!INCLUDE $(lapack_dir)\NPackage.mak
varden_dir=$(TOP)\varden
!INCLUDE $(varden_dir)\NPackage.mak

f90sources = $(f90sources) main.f90
f90objects = $(f90objects) main.obj

main.exe: $(objects)
	$(FOR) $(FFLAGS) $(LDFLAGS) /out:main.exe $(objects)

!INCLUDE $(TOP)\mk\NMakerules.mak
