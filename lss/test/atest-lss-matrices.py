#!/usr/bin/python


import coolfluid as cf
import math
import time


cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
  'LAPACK',
# 'LAPACK_SinglePrecision',
# 'GaussianElimination',
# 'GaussianElimination_SinglePrecision',
  'GMRES',
  ]

#                                                                            'fidapm13',
#                                                                            'fidapm13',
list_of_matrices = [ 'matrices/fidap/'  +a+'.mtx'      for a in [ 'fidapm03',           'fidapm33' ]] #+ [ 'matrices/drivcav/'+a+'.mtx' for a in [ 'e05r0100','e05r0500','e40r0100','e40r0500' ]]
list_of_vectors  = [ 'matrices/fidap/'  +a+'_rhs1.mtx' for a in [ 'fidapm03',           'fidapm33' ]] #+ [ 'matrices/drivcav/'+a+'.rua' for a in [ 'e05r0100','e05r0500','e40r0100','e40r0500' ]]


for solver in list_of_solvers:
  lss = cf.root.create_component('Solver_' + solver,'cf3.lss.' + solver)
  for system in zip(list_of_matrices,list_of_vectors):
    print 'test: (solver:'+solver + ' matrix:'+system[0] + ' rhs:'+system[1]+')...'

    print 'test initialize...'; d=time.time()
    lss.initialize(A=system[0],b=system[1])
    print 'test initialize. time: %.2fs.'%(time.time()-d)

    print 'test solve...'; d=time.time()
    lss.solve()
    print 'test solve. time: %.2fs.'%(time.time()-d)

    print 'test.'
  lss.delete_component()

