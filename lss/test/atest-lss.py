#!/usr/bin/python


import coolfluid as cf
import math
import time

cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
  'LAPACK',
  'GaussianElimination']

for solver in list_of_solvers:
  lss = cf.root.create_component('MySolver_' + solver,'cf3.lss.' + solver)

  lss.initialize(b="xaxcsc")
  lss.initialize(i=3,j=3)

  lss.A = [ 1, 2, 1,
            1, 1, 1,
            2, 1, 1 ]
  lss.A(i=1,j=2,value=5)
  lss.zerorow(i=1)
  lss.output()

  lss.initialize(i=4,j=4)
  lss.A = [  2, -2,  0,  0,
            -1,  2, -1,  0,
             0, -1,  2, -1,
             1,  0, -1,  2 ]
  lss.b = [  2,  4,  6,  8 ]
  lss.solve()
  lss.output()

