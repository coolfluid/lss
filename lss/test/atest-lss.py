#!/usr/bin/python

import coolfluid as cf
import math
import time

cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug



print('x1')
g = cf.root.create_component('g','cf3.lss.GaussianEliminationXXX')
print('x2')
g.execute()
print('x3')



exit()



## -----------------------------------------------------------------------------
## SETUP
## -----------------------------------------------------------------------------

#"""
#In this test we calculate laminar flow in a 3d pipe with a 90 degree bend.
#The flow enters the pipe with a uniform velocity w_in (z-direction). At the
#outlet the gauge pressure is imposed to zero. At the walls, the no-slip
#condition dictates that the velocity is zero.

#No flow goes through the symmetry plane (YZ plane), hence an additional
#boundary condition imposes u to zero on the symmetry plane.

#The pressure p and velocity v=(u,v) are calculated by solving the stationary
#incompressible Navier-Stokes equations for laminar flow:

#mass conservation:       div v = 0                                 (scalar eq.)
#momentum conservation:   1/rho grad p - nu lapl v + (v grad)v = 0  (vector eq.)
#"""

## inlet velocity
#u_in  = 1. # [m/s]
#v_in  = 0. # [m/s]
#w_in  = 0. # [m/s]

## fluid properties
#rho   = 1.   # [kg/m3] - density
#nu    = 2e-2 # [m2/s]  - kinematic viscosity

## load the mesh and make aliases to simplify the test case
#mesh   = cf.root.create_component('mesh', 'cf3.mesh.Domain').load_mesh(file=cf.URI('atest-muphys-laminarduct-3d-3kn.msh'), name='fluid')
#fluid     = [ mesh.topology.InnerCells ]
#symm_xz   = [ mesh.topology.ymax       ]
#symm_xy   = [ mesh.topology.zmax,      \
#              mesh.topology.zmin       ]
#inlet     = [ mesh.topology.xmin       ]
#outlet    = [ mesh.topology.xmax       ]
#walls     = [ mesh.topology.ymin       ]


## -----------------------------------------------------------------------------
## SIMULATOR
## -----------------------------------------------------------------------------

#simulator = cf.root.create_component('simulator', 'cf3.muphys.Stationary')
#space = simulator.set_space(dictionary=mesh.geometry, coordsystem="Cartesian")

## reference velocity
#vref = math.sqrt(u_in*u_in+v_in*v_in+w_in*w_in)

## create and initialize solution field
#solution = space.create_field(name="solution", variables="solution[4]")
#for i in range(len(solution)) :
#    solution[i][0] = 0.;        # p
#    solution[i][1] = vref/1000.; # u
#    solution[i][2] = vref/1000.; # v
#    solution[i][2] = vref/1000.; # w

## define some functions that will be used in the equations
#nu     = simulator.create_quantity(name="nu",      type="Constant").set(nu)
#invrho = simulator.create_quantity(name="inv_rho", type="Constant").set(1./rho)
#vref   = simulator.create_quantity(name="vref",    type="Constant").set(vref)

## define the system of equations and unknowns
#NS = simulator.add_system(name="NavierStokes")
#mass = NS.add_equation(name="massconservation")
#momx = NS.add_equation(name="momentumx")
#momy = NS.add_equation(name="momentumy")
#momz = NS.add_equation(name="momentumz")
#p = mass.add_unknown(name="p", field=solution, index=0)
#u = momx.add_unknown(name="u", field=solution, index=1)
#v = momy.add_unknown(name="v", field=solution, index=2)
#w = momz.add_unknown(name="w", field=solution, index=3)

## convert momentumx and momentumy into a vector equation
#mom = NS.add_equation(name="momentum")
#momx.move_component(mom.uri())
#momy.move_component(mom.uri())
#momz.move_component(mom.uri())


## -----------------------------------------------------------------------------
## NAVIER-STOKES EQUATIONS AND BOUNDARY CONDITIONS
## -----------------------------------------------------------------------------

#divv        = mass.add_term(name="divv",        type="DivvA",                   regions=fluid)
#PSPG        = mass.add_term(name="PSPG",        type="PSPGLaminarNavierStokes", regions=fluid, discretisation="WeakFEMDefault")
#invrhoGradp = mom .add_term(name="invrhoGradp", type="sAGradsB",                regions=fluid)
#nuLaplu     = momx.add_term(name="-nuLaplu",    type="-DivsAGradsB",            regions=fluid, discretisation="WeakFEMDefault")
#nuLaplv     = momy.add_term(name="-nuLaplv",    type="-DivsAGradsB",            regions=fluid, discretisation="WeakFEMDefault")
#nuLaplw     = momz.add_term(name="-nuLaplw",    type="-DivsAGradsB",            regions=fluid, discretisation="WeakFEMDefault")
#vGradv      = mom .add_term(name="vGradv",      type="vAGradvB",                regions=fluid)
#SUPG        = mom .add_term(name="SUPG",        type="SUPGLaminarNavierStokes", regions=fluid, discretisation="WeakFEMDefault")

#PSPG.p  = SUPG.p  = invrhoGradp.B = p
#divv.Ax = PSPG.vx = vGradv.Ax = vGradv.Bx = SUPG.vx = nuLaplu.B = u
#divv.Ay = PSPG.vy = vGradv.Ay = vGradv.By = SUPG.vy = nuLaplv.B = v
#divv.Az = PSPG.vz = vGradv.Az = vGradv.Bz = SUPG.vz = nuLaplw.B = w
#PSPG.nu = SUPG.nu = nuLaplu.A = nuLaplv.A = nuLaplw.A = nu
#PSPG.invRho = SUPG.invRho = invrhoGradp.A = invrho
#PSPG.vRef   = vref

## Dirichlet boundary conditions
#simulator.add_bc(name="u_in",    unknown=u, regions=inlet,   value=u_in)
#simulator.add_bc(name="v_in",    unknown=v, regions=inlet,   value=v_in)
#simulator.add_bc(name="w_in",    unknown=w, regions=inlet,   value=w_in)
#simulator.add_bc(name="u_walls", unknown=u, regions=walls,   value=0.)
#simulator.add_bc(name="v_walls", unknown=v, regions=walls,   value=0.)
#simulator.add_bc(name="w_walls", unknown=w, regions=walls,   value=0.)
#simulator.add_bc(name="p_out",   unknown=p, regions=outlet,  value=0.)
#simulator.add_bc(name="v_symm",  unknown=v, regions=symm_xz, value=0.)
#simulator.add_bc(name="w_symm",  unknown=w, regions=symm_xy, value=0.)


## -----------------------------------------------------------------------------
## COMPUTE SOLUTION
## -----------------------------------------------------------------------------

#simulator.solver = "Pardiso_mkl"
#simulator.initialize()

#cd = simulator.create_component('cd', 'cf3.muphys.tools.ConvergenceDetector')
#cd.convergence_limit = 1e-14

#iteration = 0
#while True :
#  iteration += 1
#  print "Newton iteration %d" % (iteration)
#  tic = time.clock()
#  simulator.assemble()
#  print "  assembly : %.2f s" % (time.clock()-tic)
#  simulator.calculate_L2norm()
#  print "  residual : %.2e" % (simulator.L2norm)
#  cd.add_residual(simulator.L2norm)
#  if cd.converged : break
#  tic = time.clock()
#  simulator.solve()
#  print "  solving  : %.2f s" % (time.clock()-tic)


## -----------------------------------------------------------------------------
## POSTPROCESSING
## -----------------------------------------------------------------------------

## output (disabled for acceptance tests)
## mesh.write_mesh(file=cf.URI('file:atest-muphys-laminarduct-3d.plt'), fields=[solution.uri()])

