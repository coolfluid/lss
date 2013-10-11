
TODO:
=====


  - work in progress:
      + (core plugin, dense) GaussianElimination
      + (core plugin, dense) LAPACK[2]
      + (core plugin, sparse) GMRES
      + (separate plugin, sparse) PETSc[5]
      + (separate plugin, sparse) Pardiso/Basel[7] (version 4)
      + (separate plugin, sparse) Pardiso/MKL (version 3)


  - to consider:
      + (core plugin, dense) Eigen[1]
      + (core plugin, sparse) Eigen[1] (extending the sparse matrices pruning)
      + (separate plugin, sparse) DSS/MKL (version 3)
      + (separate plugin, sparse) MUMPS
      + (separate plugin, sparse) SPARSKIT[8] (neglected uncle of pARMS?)
      + (separate plugin, sparse) SuperLU[9]
      + (separate plugin, sparse) Trilinos/AztecOO
      + (separate plugin, sparse) Trilinos/ShyLU
      + (separate plugin, sparse) WSMP


  - probably not reasonable:
      + (separate plugin, hypergraph partitioning) Mondriaan[11]
      + (separate plugin, hypergraph partitioning) PaToH (used via Zoltan or Mondriaan)[12]
      + (separate plugin, sparse) Aztec
      + (separate plugin, sparse) HSL Mathematical Software Library[4]
      + (separate plugin, sparse) PSPIKE[6] (abandoned?)
      + (separate plugin, sparse) SAMG
      + (separate plugin, sparse) SAMGp
      + (separate plugin, sparse) Trilinos (other solvers?)
      + (separate plugin, sparse) pARMS[10]


[1]: http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
[2]: http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
[3]: http://www.netlib.org/lapack/lawns/
[4]: http://www.hsl.rl.ac.uk/
[5]: http://lib.bioinfo.pl/files/courses/materials/537/Matrices_in_PETSc.pdf
[6]: http://www.pspike-project.org/
[7]: http://www.pardiso-project.org/
[8]: http://www-users.cs.umn.edu/~saad/software/SPARSKIT/
[9]: http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
[10]: http://www-users.cs.umn.edu/~saad/software/pARMS/
[11]: http://www.staff.science.uu.nl/~bisse101/Mondriaan/
[12]: http://bmi.osu.edu/~umit/software.html#patoh


