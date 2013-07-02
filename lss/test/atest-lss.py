#!/usr/bin/python

import coolfluid as cf
import math
import time

cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


for solver in [ 'GaussianEliminationA',
                'GaussianEliminationV',
                'GaussianEliminationAA',
                'GaussianEliminationVV' ]:
  g1 = cf.root.create_component('MySolver','cf3.lss.' + solver)
  g1.resize(r=3,c=3)
  g1.A = [ 1, 2, 1,
           1, 1, 1,
           2, 1, 1 ]
  g1.b = [ 8, 6, 7 ]
  g1.solve()
  g1.output()

