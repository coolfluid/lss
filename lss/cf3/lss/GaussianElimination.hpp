// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GaussianElimination_h
#define cf3_lss_GaussianElimination_h


#include <cmath>
#include "boost/progress.hpp"

#include "LibLSS.hpp"
#include "linearsystem.h"
#include "detail/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * example linear system solver, using Gaussian elimination
 * (precision is configurable)
 */
template< typename T >
class lss_API GaussianElimination : public
  linearsystem,
  detail::linearsystem< T,
    detail::dense_matrix_v< T >,
    detail::dense_matrix_v< T > >
{
  // utility definitions
  typedef detail::dense_matrix_v< T > matrix_t;
  typedef detail::dense_matrix_v< T > vector_t;
  typedef detail::linearsystem< T, matrix_t, vector_t > linearsystem_t;

  // framework interfacing
 public:
  static std::string type_name();
  GaussianElimination(const std::string& name,
                      const size_t& _size_i=size_t(),
                      const size_t& _size_j=size_t(),
                      const size_t& _size_k=1,
                      const double& _value=T() ) : linearsystem(name) { linearsystem_t::resize(_size_i,_size_j,_size_k,_value); }

  /// Linear system resizing (consistently)
  GaussianElimination& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const double& _value=double()) { linearsystem_t::resize(_size_i,_size_j,_size_k,static_cast< T >(_value)); return *this; }

  /// Linear system initialization from file(s)
  GaussianElimination& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) { linearsystem_t::initialize(_Afname,_bfname,_xfname); return *this; }

  /// Linear system initialization from vectors of values (lists, in right context)
  GaussianElimination& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< double >(),
      const std::vector< double >& vx=std::vector< double >()) { linearsystem_t::initialize(vA,vb,vx); return *this; }

  /// Linear system solving
  GaussianElimination& solve() {
    const size_t N(linearsystem_t::size(0));
    double C;
    m_x = m_b;

    boost::progress_display pbar(N-1);
    for (size_t m=0; m<N-1; ++m, ++pbar) {

      // put row with highest diagonal element on top
      C = m_A(m,m);
      for (size_t n=m+1; n<N; ++n) {
        if (std::abs(m_A(n,m)) > std::abs(C)) {
          for (size_t p=m; p<N; ++p) {
            C = m_A(m,p);
            m_A(m,p) = m_A(n,p);
            m_A(n,p) = C;
          }
          std::swap( m_x(m), m_x(n) );
          C = m_A(m,m);
        }
      }

      // check if diagonal element is (close to) zero
      if (std::abs(C)<1.e-32) {
        std::cerr
          << std::endl
          << "error: matrix is singular (line:" << m << ",C:" << std::abs(C) << ")"
          << std::endl;
        throw 42;
      }

      // normalize row m
      for (size_t n=m+1; n<N; ++n)
        m_A(m,n) /= C;
      m_x(m) /= C;

      // subtract row m from subsequent rows
      for (size_t n=m+1; n<N; ++n) {
        C = m_A(n,m);
        for (size_t p=m+1; p<N; ++p)
          m_A(n,p) -= C*m_A(m,p);
        m_x(n) -= C * m_x(m);
      }
    }

    // solve by back substitution
    m_x(N-1) /= m_A(N-1,N-1);
    for (size_t p=0; p<N-1; ++p) {
      size_t m = N-p-2;
      for (size_t n=m+1; n<N; ++n)
        m_x(m) -= m_A(m,n) * m_x(n);
    }

    return *this;
  }


 public:
        matrix_t& A()       { return m_A; }
        vector_t& b()       { return m_b; }
        vector_t& x()       { return m_x; }
  const matrix_t& A() const { return m_A; }
  const vector_t& b() const { return m_b; }
  const vector_t& x() const { return m_x; }


 protected:
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;

};


}  // namespace lss
}  // namespace cf3


#endif
