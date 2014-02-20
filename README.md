
# Coolfluid3 plugin: lss

This is a collection of plugins to solve linear systems.


## Installation

* Create a plugins directory, and clone the lss sources inside:

```
mkdir -p $PLUGIN_DIR
cd $PLUGIN_DIR
git clone https://github.com/coolfluid/lss.git $PLUGIN_DIR/lss
```

* Rerun cmake in the build directory:

```
cd $CF3_BUILD_DIR
cmake .  -DCF3_PLUGIN_DIRS=$PLUGIN_DIR -DCF3_PLUGIN_LSS=ON
```


## Introduction

In the particular case of the Finite Element Method, assembling a linear system of equations results in a well-known system matrix structure: almost of the matrix entries are zero, and a *symmetric pattern* of non-zero values. This matrix is known as a structurally symmetric *sparse* matrix. When solving the assembled linear system, there is great advantage in exploiting the *sparsity* of the matrix, because for a given useful simulation *size* (read, resolving the physical phenomena layers for many variables simultaneously) it actually becomes impossible to store every single entry in the system matrix, as well as unnecessarily multiplying anything by zero as you know.

There are some useful PDE's that, discretized with FEM, actually result in a *symmetric matrix* in the entries values themselves in addition to their position in the matrix (the structure). While more than a few linear solvers exploit this, most of the derived discretizations of interesting PDE's don't have this property so we don't try to exploit this here; non-linear solvers anyway generally handle these cases very well [citation needed].

If you don't wish to read further, the short list of linear system solvers most useful for FEM users is (the category of each solver refers to the system matrix):

* cf3.lss.mkl.dss: Intel MKL generic direct sparse solver
* cf3.lss.mkl.pardiso: Intel MKL Pardiso solver, sparse matrix direct s. (more sophisticated than the above)
* cf3.lss.pardiso.pardiso: Universita della Svizzera Italiana Pardiso solver, sparse matrix direct s.
* cf3.lss.wsmp.wsmp: IBM sparse matrix direct s.
* cf3.lss.mkl.iss_fgmres: sparse matrix iterative s.
* cf3.lss.GMRES: generic sparse matrix iterative s.
* cf3.lss.LAPACK: generic dense matrix dense s.

This list is not comprehensive but for the time being there is little reason to explore further if you are not a developer. For the linear systems considered here, the left and right-hand side vectors are actually dense matrices for simultaneous solutions, and these solver categories are explained further below.


## Direct solvers

These are slow and memory-consuming solvers, but are very accurate and typically require less non-linear iterations (less system matrix assembly steps) because the resulting solution is better than with iterative solvers for your linear problem. However, you might be far form the real solution to the problem so you can't evaluate your Jacobian matrix very well, so you might keep in mind this insight from Neils Bohr: never express yourself more clearly than you are able to think. So there is no advantage in being very precise when in fact you can only be approximate, you can just be optimistic -- but optimism can pay off though.

A rather large company that goes by the name of Intel is quite concerned that we think their processors are very fast. So they give us nice software libraries and show big computations that make us think the processors are really good, but actually the software is actually the responsible -- we have some AMD processors and they are quite good, too. This library is the Math Kernel Library **MKL** and it implements a wide range of solvers: direct, iterative, sparse or dense or any combination that makes sense here. For the direct sparse solvers they provide a nice wrapper which they creativelly called Direct Sparse Solvers **DSS**.

The only detail to remark is MKL's direct sparse matrix solvers includes the particularly useful **Pardiso** solver because it's somewhat faster than DSS. Pardiso is developed by the Universita della Svizzera Italiana and MKL includes version 3. The university is currently on version 5, and so is our provided interface, and it should be better at solving any problem - indeed... except when it isn't for unexpected reasons (upgrades don't always improve quality either), and in addition you have to be registered with the university to use it.

A rather even larger company is IBM and it also has got it's hand at trying to solve humanity big problems (they play chess quite well). In this effort is the Watson Sparse Matrix Package **WSMP** which implements some of the bright ideas of Pardiso and then some, at the trade-off of some computational effort. It actually works better than Pardiso for the most problematic matrices because it can tackle matrices with even more exotic features (citation needed). It's the most used direct solver in large-scale computers around the world and its a reference in the field.

Sometimes however issues come not from the matrix itself but from your method of calculating it, and then the linear solver doesn't help. However, a wise option might just be what you need, so choose wisely. The described solvers are all very performant in their own magical realm, the one of direct *sparse* solvers, these are their component names:

* cf3.lss.mkl.dss
* cf3.lss.mkl.pardiso
* cf3.lss.pardiso.pardiso
* cf3.lss.wsmp.wsmp

The above are sparse solvers and  of sparse direct solvers. In fact they are the most performant and accurate solvers out there, so you can generally not complain about the solver if you don't get a solution out of it: you are probably doing something you shouldn't be, the first thing to do is reevaluate your numerical strategy before the second (or *n*-th) attempt.

In addition to sparse solvers there are also dense solvers. Some PDE solving methods do assemble dense matrices, such as the Boundary Element method, or if you can't get rid of the hyperbolic behaviour of you equations then every unknown in your system depends of the value of every other unknown (but you can generally avoid this with some clever tricks).

A particularlly known direct matrix inversion method, and thus linear system solving, is the Gaussian elimination. It is also particular in its spectacular computational inneficiency, and because you actually never really look into the inverse anyway, you just want $x$ out of $A x = B$. So there are better ways to solve dense linear systems, and that's the reason for the existence of **LAPACK**. A lot of engineers have put their effort ever since computers push electrons around, its performance it without rival for dense systems and it even has its roots in Fortran - imagine that!

This applies to complex matrices too (in addition to real), which occur in many electronics physics PDEs, equations describing waves and their interactions, or your method involves this and you need it anyway (BEM stands out quite well here). BEM also traditionally results in a dense system matrix, so the following direct *dense* solver components are available:

* cf3.lss.LAPACK
* cf3.lss.LAPACK_LongPrecisionComplex


## Iterative solvers

Iterative solvers typically have lighter computational loads, at the expense of accuracy. With the employed numerical approaches, parallel architectures and quality of implementation varying so wildly (really!) there are situations where direct solvers might outperform iterative solvers.

(to be continued)

* cf3.lss.mkl.iss_fgmres
* cf3.lss.GMRES


## Only for the curious, seriously

* cf3.lss.petsc.petsc_seq
* cf3.lss.GaussianElimination
* cf3.lss.GaussianElimination_SinglePrecision
* cf3.lss.LAPACK_ShortPrecisionComplex
* cf3.lss.LAPACK_ShortPrecisionReal


