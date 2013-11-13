// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_pardiso_hpp
#define cf3_lss_mkl_pardiso_hpp


#include "LibLSS_MKL.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief Interface to Pardiso linear system solver (Intel MKL version).
 * @note the matrix structure is expected as:
 * - including a diagonal entry for each row
 * - row indices sorted in increasing order
 * @author Pedro Maciel
 */
class lss_API pardiso : public
  linearsystem< double,
    detail::sparse_matrix< double, 1>,
    detail::dense_matrix_v< double > >
{
  // utility definitions
  typedef detail::sparse_matrix< double, 1> matrix_t;
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
    const double& _value=double() );

  /// Destruction
  ~pardiso();

  /// Solve
  pardiso& solve();


  // internal functions
 private:
  int call_pardiso(int _phase);


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

  void*  pt[64];  // internal memory pointer (void* for both 32/64-bit)
  int    iparm[64],
         err,
         maxfct,
         mnum,
         msglvl,
         mtype,
         nrhs,
         phase;
  std::vector< int > perm;

};


}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

