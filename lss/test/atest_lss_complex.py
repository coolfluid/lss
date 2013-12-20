#!/usr/bin/python


import coolfluid as cf


cf.env.log_level = 3 #  1=error, 2=warning, 3=info, 4=debug


list_of_solvers = [
    'LAPACK_LongPrecisionComplex',
    'LAPACK_ShortPrecisionComplex',
    ]
for solver in list_of_solvers:
  lss = cf.root.create_component('linsys','cf3.lss.' + solver)

  lss.initialize(
      A='matrices/ibm_simple_dense_complex_A.mtx',
      b='matrices/ibm_simple_dense_complex_b.mtx' )

  lss.solve()
  lss.output(A=1,b=1)
  
  lss.delete_component()
