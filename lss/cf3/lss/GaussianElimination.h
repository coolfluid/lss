// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_lss_GaussianElimination_h
#define cf3_lss_GaussianElimination_h


#include "mlinearsystem.hpp"


namespace cf3 {
namespace lss {


// implementation of a linear system solver, using gaussian elimination (double p.)
class GaussianElimination : public mlinearsystem< double > {
 public:
  // constructor
  GaussianElimination() : mlinearsystem< double >() { issparse=false; }
  // interfacing functions
  void reset(const double& v=0.) { mlinearsystem< double >::reset(v); m_A.reset(v); }
  void solve();
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb) {
    mlinearsystem< double >::initialize(_Ne,_Nv,_Nb);
    m_A.initialize(Ne,Nv,Nb);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {}
  // indexing functions (absolute indexing)
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
  // indexing functions (block indexing)
  const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R*Nb+r,C*Nb+c); }
        double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R*Nb+r,C*Nb+c); }
  // members
 private:
  mmatrix_aa< double > m_A;
};


}  // namespace lss
}  // namespace cf3


#endif

