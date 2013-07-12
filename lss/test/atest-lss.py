#!/usr/bin/python

import coolfluid as cf
import math
import time

cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
  'GaussianEliminationA',
  'GaussianEliminationV',
  'GaussianEliminationAA',
  'GaussianEliminationVV',
  'GMRES' ]

for solver in list_of_solvers:
  g1 = cf.root.create_component('MySolver','cf3.lss.' + solver)

  g1.resize(r=3,c=3)
  g1.A = [ 1, 2, 1,
           1, 1, 1,
           2, 1, 1 ]
  g1.b = [ 8, 6, 7 ]
  g1.solve()
  g1.output()

  g1.resize(r=4,c=4)
  g1.A = [  2, -2,  0,  0,
           -1,  2, -1,  0,
            0, -1,  2, -1,
            1,  0, -1,  2 ]
  g1.b = [  2,  4,  6,  8 ]
  g1.solve()
  g1.output()
