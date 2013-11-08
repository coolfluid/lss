// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_WSMP_hpp
#define cf3_lss_WSMP_hpp


#include "LibLSS_WSMP.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Interface to WSMP linear system solver (serial version).
 * @author Pedro Maciel
 */
class lss_API WSMP : public
  linearsystem< double,
    detail::sparse_matrix_csr< double, 0 >,
    detail::dense_matrix_v< double, detail::column_oriented > >
{
  // utility definitions
  typedef detail::sparse_matrix_csr< double, 0 > matrix_t;
  typedef detail::dense_matrix_v< double, detail::column_oriented > vector_t;
  typedef linearsystem< double, matrix_t, vector_t > linearsystem_t;


 public:
  // framework interfacing
  static std::string type_name() { return "wsmp"; }


  /// Construction
  WSMP(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1,
    const double& _value=double() );

  /// Solve
  WSMP& solve();


  // internal functions
 private:
  int call_wsmp(int _task);


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

  double dparm[64];
  int    iparm[64],
         task;

};


}  // namespace lss
}  // namespace cf3


#endif

