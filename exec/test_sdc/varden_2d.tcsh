#!/bin/tcsh

# calc does integer calculations, calcfx does double precision
alias calc 'awk "BEGIN{ print \!* }" '
alias calcfx ' awk -v CONVFMT="%12.2f" -v OFMT="%.9g" "BEGIN{  print \!* }" '

#alias calc='echo "scale=4; $1" | bc'

foreach DT ( 32 64 128 256)
foreach SDCITER ( 0 1 2 )


set FIXEDDT=`calcfx 0.5/$DT`

echo "dt is " $DT
echo "sdciter is " $SDCITER
echo "fixed_dt is " $FIXEDDT

cat > input << EOF
&PROBIN

 dim_in = 2

! number of scalar quantities = density + # tracers 
 nscal = 4
! number of traces/species
 nspec = 3

! inputs for reactions
  reactions = T
! sdc_iters < 0 means use strang splitting
  sdc_iters = $SDCITER

! 1 = Crank-Nicolson or SDC, 2 = Backward Euler
! SDC uses source terms in the diffusion solve, like CN does, 
! so need to use type 1
 diffusion_type = 1

 stop_time = 0.0625
! stop_time = 1.15
! fixed_dt = 0.00390625 !128
! fixed_dt = 0.001953125 !256
 fixed_dt = $FIXEDDT
! fixed_dt = 0.015625 !32
! fixed_dt = 0.0078125 !64
! fixed_dt = 0.000078125 !64
! fixed_dt = 0.03125 !16

 cflfac = 0.5
 init_shrink = 0.1

! visc_coef = 0.001
 visc_coef = 0.00
 diff_coef = 0.01

  n_rxn_steps = 100
!  k_rxn1 = 0.0
!  k_rxn2 = 0.0
  k_rxn1 = 1050.0d0
  k_rxn2 = 132.0d0

! default = T
 mass_fractions = F


!!!! Input for fixed grids !!!!
! fixed_grids = "gr0_2d_2levels"
 fixed_grids = "gr0_2d.$DT"

max_grid_size = 1024

!!!! Inputs for adaptive grids !!!!
! nlevs = 3
! n_cellx = 32
! n_celly = 32    
!   n_cellx = 64
!   n_celly = 64
! ref_ratio = 2
! regrid_int = 2

! nothing is ever done with this data in varden
! prob_lo_x = 0
! prob_lo_y = -0.5

 prob_hi_x = 1.0
 prob_hi_y = 1.0

 max_step  = 100
 init_iter = 1

 plot_int  = 10
 chk_int   = 10

 grav = 0

! bc flags (definintions in boxlib/bc.f90: 
!	    -1: periodic/interior
!           11: inflow
!           12: outflow
!           13: symmetry
!           14: slip wall
!           15: no-slip wall
 bcx_lo = -1
 bcx_hi = -1
 bcy_lo = -1
 bcy_hi = -1

! pmask_x = T : periodic in x direction
! pmask_y = T : periodic in y direction
 pmask_x = T
 pmask_y = T

 verbose = 0
 mg_verbose = 0
 cg_verbose = 0

! restart = 1200
/
EOF

mpiexec -n 2 LRmain.Linux.Intel.mpi.exe input > output_${SDCITER}_${DT}

if ( $DT == 32 ) then
 mv plt0004 convergence/plt${SDCITER}032
 mv chk0004 convergence/chk${SDCITER}032
endif
if ( $DT == 64 ) then
 mv plt0008 convergence/plt${SDCITER}064
 mv chk0008 convergence/chk${SDCITER}064
endif
if ( $DT == 128 ) then
 mv plt0016 convergence/plt${SDCITER}128
 mv chk0016 convergence/chk${SDCITER}128
endif
if ( $DT == 256 ) then
 mv plt0032 convergence/plt${SDCITER}256
 mv chk0032 convergence/chk${SDCITER}256
endif

#mv plt0

end
end
