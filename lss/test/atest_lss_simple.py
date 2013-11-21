#!/usr/bin/python


import coolfluid as cf
import math
import time


cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
  'LAPACK',
  'LAPACK_SinglePrecision',
  'GaussianElimination',
  'GaussianElimination_SinglePrecision',
  'GMRES',
  'Dlib',
# 'mkl.pardiso',
  ]
for solver in list_of_solvers:
  lss = cf.root.create_component('MySolver_' + solver,'cf3.lss.' + solver)

  print solver
  lss.initialize(i=4,j=4,k=2)
  lss.A = [  2, -2,  0,  0,
            -1,  2, -1,  0,
             0, -1,  2, -1,
             1,  0, -1,  2 ]
  lss.b = [  2,  4,
             4,  8,
             6, 12,
             8, 16 ]
  lss.output(A=3,b=3,x=3)
  lss.solve()

