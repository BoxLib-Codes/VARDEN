f90sources = $(f90sources) $(varden_dir)\advance_premac.f90
f90sources = $(f90sources) $(varden_dir)\checkpoint.f90
f90sources = $(f90sources) $(varden_dir)\define_bc_tower.f90
f90sources = $(f90sources) $(varden_dir)\estdt.f90
f90sources = $(f90sources) $(varden_dir)\hgproject.f90
f90sources = $(f90sources) $(varden_dir)\initdata.f90
f90sources = $(f90sources) $(varden_dir)\inlet_bc.f90
f90sources = $(f90sources) $(varden_dir)\laplac.f90
f90sources = $(f90sources) $(varden_dir)\macproject.f90
f90sources = $(f90sources) $(varden_dir)\mkdivu.f90
f90sources = $(f90sources) $(varden_dir)\mkflux.f90
f90sources = $(f90sources) $(varden_dir)\mkforce.f90
f90sources = $(f90sources) $(varden_dir)\mkutrans.f90
f90sources = $(f90sources) $(varden_dir)\ml_solve.f90
f90sources = $(f90sources) $(varden_dir)\restart.f90
f90sources = $(f90sources) $(varden_dir)\scalar_advance.f90
f90sources = $(f90sources) $(varden_dir)\setbc.f90
f90sources = $(f90sources) $(varden_dir)\slope.f90
f90sources = $(f90sources) $(varden_dir)\update.f90
f90sources = $(f90sources) $(varden_dir)\varden.f90
f90sources = $(f90sources) $(varden_dir)\velocity_advance.f90
f90sources = $(f90sources) $(varden_dir)\viscsolve.f90

f90objects = $(f90objects) $(obj_dir)\advance_premac.obj
f90objects = $(f90objects) $(obj_dir)\checkpoint.obj
f90objects = $(f90objects) $(obj_dir)\define_bc_tower.obj
f90objects = $(f90objects) $(obj_dir)\estdt.obj
f90objects = $(f90objects) $(obj_dir)\hgproject.obj
f90objects = $(f90objects) $(obj_dir)\initdata.obj
f90objects = $(f90objects) $(obj_dir)\inlet_bc.obj
f90objects = $(f90objects) $(obj_dir)\laplac.obj
f90objects = $(f90objects) $(obj_dir)\macproject.obj
f90objects = $(f90objects) $(obj_dir)\mkdivu.obj
f90objects = $(f90objects) $(obj_dir)\mkflux.obj
f90objects = $(f90objects) $(obj_dir)\mkforce.obj
f90objects = $(f90objects) $(obj_dir)\mkutrans.obj
f90objects = $(f90objects) $(obj_dir)\ml_solve.obj
f90objects = $(f90objects) $(obj_dir)\restart.obj
f90objects = $(f90objects) $(obj_dir)\scalar_advance.obj
f90objects = $(f90objects) $(obj_dir)\setbc.obj
f90objects = $(f90objects) $(obj_dir)\slope.obj
f90objects = $(f90objects) $(obj_dir)\update.obj
f90objects = $(f90objects) $(obj_dir)\varden.obj
f90objects = $(f90objects) $(obj_dir)\velocity_advance.obj
f90objects = $(f90objects) $(obj_dir)\viscsolve.obj

{$(varden_dir)}.f90{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)\"
	$(FOR) /c $(FFLAGS) $< /object:$(obj_dir)
