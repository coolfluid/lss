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


// implementation of a linear system solver, using gaussian elimination (double p.)
struct lss_API GaussianElimination :
  public linearsystem< double >
{
  typedef linearsystem< double > linearsystem_t;

  // framework interfacing
  static std::string type_name() { return "GaussianElimination"; } // (mandatory!)
  GaussianElimination(const std::string& name) :
    linearsystem_t(name) {}

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
    b(r) = 0;
    return *this;
  }

  GaussianElimination& clear() {
    m_A.clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }

  GaussianElimination& solve();

  void output_A(std::ostream& out) const { m_A.output(out); }

  // variables
  matrix_dense_vv< double > m_A;
};


}  // namespace lss
}  // namespace cf3


#endif
