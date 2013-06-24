// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_lss_mlinearsystem_h
#define cf3_lss_mlinearsystem_h

#include <vector>
#include "mmatrix.h"


namespace cf3 {
namespace lss {


/*
 * description of a linear system solver
 * (solution and right-hand side vectors are included, matrix should be
 * included in the derived class)
 */
template< typename T >
class mlinearsystem {
 public:
  // cons/destructor
  mlinearsystem() : issparse(true), Ne(0), Nv(0), Nb(1) {}
  virtual ~mlinearsystem() {}
  // interfacing methods
  virtual void reset(const T& v=T()) { m_X.assign(Nb*Ne,v); m_B.assign(Nb*Nv,v); }
  virtual void solve()                   = 0;
  virtual void zerorow(const unsigned r) = 0;
  virtual void zerorow(const unsigned R, const unsigned r) { zerorow(R*Nb+r); }
  virtual void print(std::ostream& o, bool pmatrix=false) {
    o << "m::mlinearsystem::X(Nv:" << Nv << ",Nb:" << Nb << "):" << std::endl;
    for (unsigned J=0; J<Nv; ++J) for (unsigned j=0; j<Nb; ++j)
      o << ' ' << X(J,j) << std::endl;
    o << "m::mlinearsystem::B(Ne:" << Nv << ",Nb:" << Nb << "):" << std::endl;
    for (unsigned I=0; I<Ne; ++I) for (unsigned i=0; i<Nb; ++i)
      o << ' ' << B(I,i) << std::endl;
    o << "m::mlinearsystem::A(issparse?" << (issparse? "true":"false") << "):" << std::endl;
    for (unsigned I=0; I<Ne && pmatrix; ++I) for (unsigned i=0; i<Nb; ++i) {
      for (unsigned J=0; J<Nv; ++J) for (unsigned j=0; j<Nb; ++j)
        o << ' ' << A(I,J,i,j);
      o << std::endl;
    }
  }
  // initialize methods for dense/sparse variations
  virtual void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb=1) {
    Ne = _Ne;
    Nv = _Nv;
    Nb = _Nb;
    mlinearsystem< T >::reset(T());
  }
  virtual void initialize(const std::vector< std::vector< unsigned > >& nz) = 0;
  // indexing functions (absolute indexing)
  virtual const T& A(const unsigned r, const unsigned c) const = 0;
  virtual       T& A(const unsigned r, const unsigned c)       = 0;
          const T& X(const unsigned r) const { return m_X[r]; }
                T& X(const unsigned r)       { return m_X[r]; }
          const T& B(const unsigned c) const { return m_B[c]; }
                T& B(const unsigned c)       { return m_B[c]; }
  // indexing functions (block indexing)
  virtual const T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const = 0;
  virtual       T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       = 0;
          const T& X(const unsigned R, const unsigned r) const { return X(R*Nb+r); }
                T& X(const unsigned R, const unsigned r)       { return X(R*Nb+r); }
          const T& B(const unsigned C, const unsigned c) const { return B(C*Nb+c); }
                T& B(const unsigned C, const unsigned c)       { return B(C*Nb+c); }
 public:
  // members
  bool issparse;
 protected:
 public:  // FIXME bypass, should be protected
  // members
  unsigned Ne;  // number of (block) equations
  unsigned Nv;  // ... (block) variables
  unsigned Nb;  // ... block size
  std::vector< T > m_X;
  std::vector< T > m_B;
};


}  // namespace lss
}  // namespace cf3


#endif

