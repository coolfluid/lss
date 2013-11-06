// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_pardiso_h
#define cf3_lss_pardiso_h


#include <cstdio>  // for sscanf

#include "LibLSS_PARDISO.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * implementation of ...
 */
class lss_pardiso_API pardiso : public
  linearsystem< double,
    detail::sparse_matrix_csr< double, 1 >,
    detail::dense_matrix_v< double > >
{
  // utility definitions
  typedef detail::sparse_matrix_csr< double, 1 > matrix_t;
  typedef detail::dense_matrix_v< double > vector_t;
  typedef linearsystem< double, matrix_t, vector_t > linearsystem_t;


 public:
  // framework interfacing
  static std::string type_name() { return "pardiso"; }


  /// Construction
  pardiso(const std::string& name,
        const size_t& _size_i=size_t(),
        const size_t& _size_j=size_t(),
        const size_t& _size_k=1,
        const double& _value=double() ) : linearsystem_t(name) {

    for (int i=0; i<64; ++i)  iparm[i] = 0;
    for (int i=0; i<64; ++i)  dparm[i] = 0.;

    char* nthreads = getenv("OMP_NUM_THREADS");
    sscanf(nthreads? nthreads:"1","%d",&iparm[2]);
    std::cout << "info: number of threads: " << iparm[2] << " (OMP_NUM_THREADS: " << (nthreads? "":"not ") << "set)" << std::endl;

    iparm[ 7] = 0;  // max numbers of iterative refinement steps
    iparm[31] = 0;  // [0|1] sparse direct solver or multi-recursive iterative solver
    maxfct    = 1,  // maximum number of numerical factorizations
    mnum      = 1,  // which factorization to use
    mtype     = 1,  // real structurally symmetric matrix
    nrhs      = 1;  // number of right hand sides

    linearsystem_t::initialize(_size_i,_size_j,_size_k,_value);
  }

  /// Solve
  pardiso& solve();


  // internal functions
 private:
  int call_pardisoinit();
  int call_pardiso(int phase, int msglvl);


  // linear system components access
 public:
        matrix_t& A()       { return m_A; }
        vector_t& b()       { return m_b; }
        vector_t& x()       { return m_x; }
  const matrix_t& A() const { return m_A; }
  const vector_t& b() const { return m_b; }
  const vector_t& x() const { return m_x; }


  // members
 protected:
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;
std::vector< int > ia, ja;
std::vector< double > va, vb, vx;

  void*  pt[64];  // internal memory pointer (void* for both 32/64-bit)
  int    iparm[64];
  double dparm[64];
  int    maxfct,
         mnum,
         mtype,
         nrhs;

};


}  // namespace lss
}  // namespace cf3


#endif

