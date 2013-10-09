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
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * example linear system solver, using Gaussian elimination
 * (precision is configurable)
 */
template<
    typename T,
    typename INDEX=index_hierarchy_t< index_hierarchy_t_end > >
struct lss_API GaussianElimination :
  public linearsystem< T, INDEX >
{
  // utility definitions
  typedef lss_matrix::dense_matrix_v< T > matrix_t;
  typedef lss_matrix::dense_matrix_v< T > vector_t;
  typedef linearsystem< T, INDEX > linearsystem_base_t;

  // framework interfacing
  static std::string type_name();
  GaussianElimination(const std::string& name) :
    linearsystem_base_t(name) {}

  // linear system solver addressing
  const double& A(const size_t i) const { return m_A(i); }
        double& A(const size_t i)       { return m_A(i); }
  const double& A(const size_t r, const size_t c) const { return m_A(r,c); }
        double& A(const size_t r, const size_t c)       { return m_A(r,c); }

  // linear system solver interfacing
  size_t size()               const { return m_A.size(); }
  size_t size(const size_t d) const { return m_A.size(d); }

  GaussianElimination& resize(size_t Nequations, size_t Nvariables, const double& v=double()) {
    if (!Nequations || !Nvariables)
      return clear();
    m_A.resize(Nequations,Nvariables);
    m_b.assign(Nequations,v);
    m_x.assign(Nvariables,v);
    return *this;
  }

  GaussianElimination& zerorow(const size_t r) {
    m_A.zerorow(r);
    m_b(r) = 0;
    m_x(r) = 0;
    return *this;
  }

  GaussianElimination& clear() {
    m_A.clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }

  GaussianElimination& solve() {
    const size_t N(size(0));
    double C;
    m_x = m_b;

    boost::progress_display pbar(N-1);
    for (size_t m=0; m<N-1; ++m, ++pbar) {

      // put row with highest diagonal element on top
      C = A(m,m);
      for (size_t n=m+1; n<N; ++n) {
        if (std::abs(A(n,m)) > std::abs(C)) {
          for (size_t p=m; p<N; ++p) {
            C = A(m,p);
            A(m,p) = A(n,p);
            A(n,p) = C;
          }
          std::swap( m_x(m), m_x(n) );
          C = A(m,m);
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
        A(m,n) /= C;
      m_x(m) /= C;

      // subtract row m from subsequent rows
      for (size_t n=m+1; n<N; ++n) {
        C = A(n,m);
        for (size_t p=m+1; p<N; ++p)
          A(n,p) -= C*A(m,p);
        m_x(n) -= C * m_x(m);
      }
    }

    // solve by back substitution
    m_x(N-1) /= A(N-1,N-1);
    for (size_t p=0; p<N-1; ++p) {
      size_t m = N-p-2;
      for (size_t n=m+1; n<N; ++n)
        m_x(m) -= A(m,n) * m_x(n);
    }

    return *this;
  }

  void output_A(std::ostream& out) const { m_A.output(out); }

  // variables
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;
};


}  // namespace lss
}  // namespace cf3


#endif
