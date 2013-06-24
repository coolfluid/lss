#ifndef mlinearsystem_h
#define mlinearsystem_h

#include <vector>
#include "ext/xmlParser.h"
#include "mmatrix.h"


namespace m {


/*
 * description of a linear system solver
 * (solution and right-hand side vectors are included, matrix should be
 * included in the derived class)
 */
template< typename T >
class mlinearsystem {
 public:
  // cons/destructor
  mlinearsystem() :
    issparse(true),
    xml(XMLNode::createXMLTopNode("mlinearsystem")),
    Ne(0), Nv(0), Nb(1)    {}
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
  XMLNode xml;
 protected:
 public:  // FIXME bypass, should be protected
  // members
  unsigned Ne;  // number of (block) equations
  unsigned Nv;  // ... (block) variables
  unsigned Nb;  // ... block size
  std::vector< T > m_X;
  std::vector< T > m_B;
};


// implementation of a linear system solver, using gaussian elimination (double p.)
class ls_gauss : public mlinearsystem< double > {
 public:
  // constructor
  ls_gauss() : mlinearsystem< double >() { issparse=false; }
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


}  // namespace m


#endif

