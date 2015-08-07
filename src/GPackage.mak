f90sources += advance_premac.f90  
f90sources += advance_timestep.f90  
f90sources += checkpoint.f90
f90sources += estdt.f90
f90sources += explicit_diffusive_term.f90
f90sources += hgproject.f90
f90sources += hg_multigrid.f90
f90sources += initdata.f90
f90sources += initialize.f90
ifdef HYPRE
f90sources +=  hg_hypre.f90
f90sources += mac_hypre.f90
else
f90sources += mac_hypre_stub.f90
f90sources +=  hg_hypre_stub.f90
endif
f90sources += mac_multigrid.f90
f90sources += macproject.f90
f90sources += make_at_halftime.f90
f90sources += make_new_grids.f90
f90sources += makevort.f90
f90sources += mkflux.f90
f90sources += mkforce.f90
f90sources += probin.f90
f90sources += proj_parameters.f90
f90sources += restart.f90
f90sources += scalar_advance.f90  
f90sources += slope.f90
f90sources += tag_boxes.f90
f90sources += update.f90  
f90sources += varden.f90
f90sources += velocity_advance.f90  
f90sources += velpred.f90
f90sources += viscsolve.f90
f90sources += main.f90
