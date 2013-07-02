// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GaussianElimination_h
#define cf3_lss_GaussianElimination_h


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * example linear system solver, using Gaussian elimination
 * precision and matrix types are arbitrary, though it is only expected to work
 * for dense matrices
 */
template< typename T, typename MATRIX >
struct lss_API GaussianElimination :
  public linearsystem< T >
{
  typedef linearsystem< T > linearsystem_t;

  // framework interfacing
  static std::string type_name() ;//{ return "GaussianElimination"; } // (mandatory!)
  GaussianElimination(const std::string& name) :
    linearsystem_t(name) {}

  // linear system solver addressing
  const T& A(const size_t i) const { return m_A(i); }
        T& A(const size_t i)       { return m_A(i); }
  const T& A(const size_t r, const size_t c) const { return m_A(r,c); }
        T& A(const size_t r, const size_t c)       { return m_A(r,c); }

  // linear system solver interfacing
  size_t size()               const { return m_A.size(); }
  size_t size(const size_t d) const { return m_A.size(d); }

  GaussianElimination& resize(size_t Nequations, size_t Nvariables, const T& v=T()) {
    if (!Nequations || !Nvariables)
      return clear();
    m_A.resize(Nequations,Nvariables);
    linearsystem_t::m_b.assign(Nequations,v);
    linearsystem_t::m_x.assign(Nvariables,v);
    return *this;
  }

  GaussianElimination& zerorow(const size_t r) {
    m_A.zerorow(r);
    linearsystem_t::b(r) = 0;
    return *this;
  }

  GaussianElimination& clear() {
    m_A.clear();
    linearsystem_t::m_b.clear();
    linearsystem_t::m_x.clear();
    return *this;
  }

  GaussianElimination& solve() {
    const size_t N(size(0));
    T C;
    linearsystem_t::m_x = linearsystem_t::m_b;

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
          std::swap( linearsystem_t::x(m), linearsystem_t::x(n) );
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
      linearsystem_t::x(m) /= C;

      // subtract row m from subsequent rows
      for (size_t n=m+1; n<N; ++n) {
        C = A(n,m);
        for (size_t p=m+1; p<N; ++p)
          A(n,p) -= C*A(m,p);
        linearsystem_t::x(n) -= C * linearsystem_t::x(m);
      }
    }

    // solve by back substitution
    linearsystem_t::x(N-1) /= A(N-1,N-1);
    for (size_t p=0; p<N-1; ++p) {
      size_t m = N-p-2;
      for (size_t n=m+1; n<N; ++n)
        linearsystem_t::x(m) -= A(m,n) * linearsystem_t::x(n);
    }

    return *this;
  }





  void output_A(std::ostream& out) const { m_A.output(out); }

  // variables
  MATRIX m_A;
};


}  // namespace lss
}  // namespace cf3


#endif
