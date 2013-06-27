// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_h
#define cf3_lss_linearsystem_h


#include <vector>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/keyword.hpp>
#include "common/Action.hpp"


namespace cf3 {
namespace lss {


/**
 * description of a matrix row/column indexer
 * the idea is to return a single size_t index to a (very long and boring)
 * vector of entries that populate a sparse, dense, or block-addressable matrix.
 */
struct index {
  virtual void print(std::ostream& o) = 0;
#if 0
unsigned Ne;  // number of (block) equations
unsigned Nv;  // ... (block) variables
unsigned Nb;  // ... block size
virtual const T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const = 0;
virtual       T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       = 0;
const T& X(const unsigned R, const unsigned r) const { return X(R*Nb+r); }
T& X(const unsigned R, const unsigned r)       { return X(R*Nb+r); }
const T& B(const unsigned C, const unsigned c) const { return B(C*Nb+c); }
T& B(const unsigned C, const unsigned c)       { return B(C*Nb+c); }
virtual void resize(const std::vector< std::vector< unsigned > >& nz) = 0;
virtual void zerorow(const unsigned R, const unsigned r) { zerorow(R*Nb+r); }
void print(std::ostream& o) {     o << "m::linearsystem::B(Ne:" << Nv << ",Nb:" << Nb << "):" << std::endl;     }
// indexing functions (block indexing)
const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R*Nb+r,C*Nb+c); }
double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R*Nb+r,C*Nb+c); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return P::operator()(R,C,r,c); };
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return P::operator()(R,C,r,c); };

void zerorow(const unsigned R, const unsigned r) { P::zerorow(R,r); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }

void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
void zerorow(const unsigned R, const unsigned r) {
const int b = (int) P::Nb;
const int i = getindex(R,R,r,0);
for (int j=0; j<b*b*(bpntr[R+1]-bpntr[R]); j+=b)
val[i+j] = T();
}
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }


  BOOST_PARAMETER_FUNCTION(
    (void),                // 1. parenthesized return type
    depth_first_search,    // 2. name of the function template
    tag,                   // 3. namespace of tag types

    (required
      (graph, *) ) // 4. one required parameter, and

    (optional              //    four optional parameters, with defaults
      (visitor,           *, boost::dfs_visitor<>())
      (root_vertex,       *, *vertices(graph).first)
      (index_map,         *, get(boost::vertex_index,graph))
      (in_out(color_map), *, default_color_map(num_vertices(graph), index_map)) )
    )
    {
      // ... body of function goes here...
      // use graph, visitor, index_map, and color_map
    }

#endif
};


/**
 * description of a matrix
 * (T: storage type, P: parent class)
 * it is based on parent class template extension (CRTP) to simulate virtual
 * methods (not a true polymorphic type) so it can't be used with factories.
 */
template< typename T, class P >
struct matrix {

  // constructor
  matrix() : Nr(0), Nc(0) {}

  // indexing
  const T& operator()(const size_t r, const size_t c) const { return P::operator()(r,c); }
        T& operator()(const size_t r, const size_t c)       { return P::operator()(r,c); }

  // interfacing
  matrix& clear(const T& v=T())     { P::clear(v);   return *this; }
  matrix& zerorow(const size_t r) { P::zerorow(r); return *this; }

  // size/resizing
  size_t size()                { return Nr*Nc; }
  size_t size(const size_t& d) { return (d==0? Nr : (d==1? Nc : 0)); }
  size_t resize(size_t _Nr, size_t _Nc) { Nr=_Nr; Nc=_Nc; return size(); }

  // members
  size_t Nr;  // number of rows
  size_t Nc;  // ... columns

};


/**
 * implementation of a dense matrix, std::vector< std::vector< T > > based
 */
template< typename T >
struct matrix_dense_vv : matrix< T,matrix_dense_vv< T > > {
  typedef matrix< T,matrix_dense_vv< T > > P;

  // indexing functions
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }

  // interfacing functions
  matrix_dense_vv& operator=(const T& v) { return clear(v); }
  matrix_dense_vv& clear(const T& v=T())   { if (size()) a.assign(P::Nr,std::vector< T >(P::Nc,v)); return *this; }
  matrix_dense_vv& zerorow(const size_t r) { if (size()) a[r].assign(P::Nc,T()); return *this; }

  // size/resizing
  size_t size()                { return P::size(); }
  size_t size(const size_t& d) { return P::size(d); }
  size_t resize(size_t _Nr, size_t _Nc) { P::resize(_Nr,_Nc); clear(); return size(); }

  // members
  std::vector< std::vector< T > > a;
};


/**
 * implementation of a dense matrix, T[][] array based
 */
template< typename T >
struct matrix_dense_aa : matrix< T,matrix_dense_aa< T > > {
  typedef matrix< T,matrix_dense_aa< T > > P;

  // destructor
  ~matrix_dense_aa() {
    if (P::Nr || P::Nc) {
      delete[] a[0];
      delete[] a;
    }
  }

  // indexing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }

  // interfacing
  matrix_dense_aa& operator=(const T& v) { return clear(v); }
  matrix_dense_aa& clear(const T& v=T())   { for (size_t r=0; r<P::Nr; ++r) for (size_t c=0; c<P::Nc; ++c) a[r][c] = v; return *this; }
  matrix_dense_aa& zerorow(const size_t r) { for (size_t c=0; c<P::Nc; ++c) a[r][c] = T(); return *this; }

  // size/resizing
  size_t size()                { return P::size(); }
  size_t size(const size_t& d) { return P::size(d); }
  void resize(size_t _Nr, size_t _Nc) {
    if (_Nr && _Nc) {
      P::resize(_Nr,_Nc);
      a    = new T*[ P::Nr ];
      a[0] = new T [ P::Nr * P::Nc ];
      for (size_t r=1; r<P::Nr; ++r)
        a[r] = a[r-1] + P::Nc;
      clear();
    }
  }

  // members
  T **a;

};


/**
 * implementation of a sparse matrix, compressed sparse rows (CSR) format
 * note: BASE={0,1}: {0,1}-based indexing (other values probably don't work)
 */
template< typename T, int BASE >
struct matrix_sparse_csr : matrix< T,matrix_sparse_csr< T,BASE > > {
  typedef matrix< T,matrix_sparse_csr< T,BASE > > P;

  // cons/destructor
  matrix_sparse_csr() : P(), zero(T()), nnz(0), nnu(0) {}
  ~matrix_sparse_csr() {
    if (nnz) {
      delete[] a;
      delete[] ja;
      delete[] ia;
    }
  }

  // indexing
  const T& operator()(const size_t r, const size_t c) const { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }
        T& operator()(const size_t r, const size_t c)       { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }

  // interfacing
  void clear(const T& v=T())                       { for (int i=0; i<nnz; ++i) a[i] = v; }
  void zerorow(const size_t r)                   { for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k) a[k] = T(); }

  // resizing
  void resize(size_t _Nr, size_t _Nc) { P::resize(_Nr,_Nc); }
  void resize(const std::vector< std::vector< size_t > >& nz) {

    // set row indices
    nnu = (int) nz.size();
    ia = new int[nnu+1];
    ia[0] = BASE;
    int k = 0;
    for (size_t r=0; r<nz.size(); ++r)
      ia[k+1] = ia[k] + (int) nz[r].size();
    nnz = ia[nnu]-BASE;

    // set column indices
    ja = new int[nnz];
    for (size_t r=0; r<nz.size(); ++r) {
      k = ia[r]-BASE;
      for (size_t i=0; i<nz[r].size(); ++i)
        ja[k++] = (int) (nz[r][i]) + BASE;
    }

    // set entries
    a = new T[nnz];
    clear();
  }
  int getindex(const size_t r, const size_t c) const {
    for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k)
      if (ja[k]-BASE==(int) c)
        return k;
    return -1;
  }

  // members
  T zero;
  int *ia;
  int *ja;
  T   *a;
  int nnz;
  int nnu;

};


// NOTE: missing Aztec MSR and VBR (maybe this will never happen)


/**
 * description of a linear system
 * solution and right-hand side vectors are included, matrix should be
 * included (by aggregation) in derived linearsystem classes
 */
template< typename T, typename M >
class linearsystem : public common::Action {

 public:
  // framework interfacing
  linearsystem(const std::string& name) : common::Action(name) {}
  virtual ~linearsystem() {}
  void execute() { solve(); }

 public:
  // linear system solver interfacing
  virtual linearsystem& solve()                 = 0;
  virtual linearsystem& zerorow(const size_t r) = 0;
  linearsystem& resize(size_t Nequations, size_t Nvariables, const T& v=T()) {
    (Nequations? m_b.resize(Nequations) : m_b.clear());
    (Nvariables? m_x.resize(Nvariables) : m_x.clear());
    (Nvariables && Nequations? m_A.resize(Nequations,Nvariables) : (void) m_A.clear() );
    return clear(v);
  }
  linearsystem& clear(const T& v=T()) {
    if (m_x.size()) m_x.assign(m_x.size(),v);
    if (m_b.size()) m_b.assign(m_b.size(),v);
    if (m_A.size()) m_A = v;
    return *this;
  }

  void print_A(std::ostream& o) { m_A.print(o); o << std::endl; }
  void print_x(std::ostream& o) { std::copy(m_x.begin(), m_x.end(), std::ostream_iterator< double >(o,'\n'));  o << std::endl; }
  void print_b(std::ostream& o) { std::copy(m_b.begin(), m_b.end(), std::ostream_iterator< double >(o,'\n'));  o << std::endl; }

  // indexing
  virtual const T& A(const size_t r, const size_t c) const = 0;
  virtual       T& A(const size_t r, const size_t c)       = 0;
          const T& x(const size_t r) const { return m_x[r]; }
                T& x(const size_t r)       { return m_x[r]; }
          const T& b(const size_t c) const { return m_b[c]; }
                T& b(const size_t c)       { return m_b[c]; }

 protected:
  // member variables
  std::vector< T > m_x;
  std::vector< T > m_b;
  M                m_A;

};


}  // namespace lss
}  // namespace cf3


#endif

