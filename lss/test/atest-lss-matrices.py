#!/usr/bin/python


import coolfluid as cf
import math
import time


cf.env.log_level = 4 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
  #'LAPACK',
# 'LAPACK_SinglePrecision',
# 'GaussianElimination',
# 'GaussianElimination_SinglePrecision',
  'GMRES',
  'pardiso',
  ]

list_of_matrices = []
list_of_vectors  = []
#list_of_matrices.extend([ 'matrices/'+a for a in [ 'simple_sparse_matrix_base_0.csr','simple_sparse_matrix_base_1.csr','simple_dense_matrix.mtx' ]])
#list_of_vectors .extend([ None, None, None ])
#list_of_matrices.extend([ 'matrices/fidap/'  +a+'.mtx'      for a in [ 'fidapm03','fidapm13','fidapm33' ]])
#list_of_vectors .extend([ 'matrices/fidap/'  +a+'_rhs1.mtx' for a in [ 'fidapm03','fidapm13','fidapm33' ]])
list_of_matrices.extend([ 'matrices/drivcav/'+a+'.mtx'      for a in [ 'e05r0100','e05r0500' ]]) #,'e40r0100','e40r0500' ]])
list_of_vectors .extend([ 'matrices/drivcav/'+a+'_rhs1.mtx' for a in [ 'e05r0100','e05r0500' ]]) #,'e40r0100','e40r0500' ]])


for solver in list_of_solvers:
  lss = cf.root.create_component('Solver_'+solver,'cf3.lss.'+solver)
  for system in zip(list_of_matrices,list_of_vectors):
    mat = system[0] if system[0] else ''
    vec = system[1] if system[1] else ''
    print 'test: (solver:'+solver + ' matrix:'+mat + ' rhs:'+vec+')...'

    print 'test initialize...'; d=time.time()
    lss.initialize(A=mat,b=vec)
    print 'test initialize. time: %.2fs.'%(time.time()-d)

    print 'test solve...'; d=time.time()
    lss.b = [1]
    lss.solve()
    print 'test solve. time: %.2fs.'%(time.time()-d)
 
    lss.output(A=1,b=1,x=1)
    lss.A = [ 0 ]
    print 'test.'
  lss.delete_component()

