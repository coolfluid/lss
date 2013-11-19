#!/usr/bin/python


import coolfluid as cf
import math
import time


cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


solvers = [
# 'GaussianElimination',
# 'GaussianElimination_SinglePrecision',
  'GMRES',
# 'LAPACK',
# 'LAPACK_SinglePrecision',
  'pardiso.pardiso',
# 'mkl.pardiso',
# 'mkl.dss',
  'wsmp.wsmp',
  ]

systems=[
# ('matrices/intel_mkl_simple_sparse_matrix.csr', ''),
  ('matrices/samg_demo_matrix.csr', 'matrices/samg_demo_rhs.mtx'),
# ('matrices/drivcav/e05r0100.mtx', 'matrices/drivcav/e05r0100_rhs1.mtx'),
# ('matrices/drivcav/e05r0500.mtx', 'matrices/drivcav/e05r0500_rhs1.mtx'),
# ('matrices/drivcav/e40r0100.mtx', 'matrices/drivcav/e40r0100_rhs1.mtx'),
# ('matrices/drivcav/e40r0500.mtx', 'matrices/drivcav/e40r0500_rhs1.mtx'),
# ('matrices/fidap/fidapm03.mtx',   'matrices/fidap/fidapm03_rhs1.mtx'),
# ('matrices/fidap/fidapm13.mtx',   'matrices/fidap/fidapm13_rhs1.mtx'),
# ('matrices/fidap/fidapm33.mtx',   'matrices/fidap/fidapm33_rhs1.mtx'),
  ]
for t in solvers:
  lss = cf.root.create_component('Solver_'+t,'cf3.lss.'+t)
  for s in systems:
    print 'test: (solver:'+t + ' matrix:'+s[0] + ' rhs:'+s[1]+')...'

    d=time.time()
    lss.initialize(A=s[0],b=s[1])
    if not len(s[1]): lss.b=[2]
    print '  initialize time: %.2fs.'%(time.time()-d)
    #lss.output(A=3)

    d=time.time()
    lss.solve()
    lss.output(file='asd_'+t,x=1)
    print '  solve time: %.2fs.'%(time.time()-d)

  lss.delete_component()

