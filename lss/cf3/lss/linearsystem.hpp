// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_h
#define cf3_lss_linearsystem_h


#include <vector>
//include <boost/foreach.hpp>
//include <boost/parameter/preprocessor.hpp>
//include <boost/parameter/keyword.hpp>
#include "common/Action.hpp"
#include "common/Log.hpp"


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
 * only describes an interface, storage is provided by the derived structs.
 * (T: storage type, P: parent class)
 * it is based on parent class template extension (CRTP) to simulate virtual
 * methods (not a true polymorphic type) so it can't be used with factories.
 */
template< typename T, class P >
struct matrix {

  // constructor
  matrix() : Nr(0), Nc(0) {}

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return P::operator()(r,c); }
        T& operator()(const size_t r, const size_t c)       { return P::operator()(r,c); }

  matrix& operator=(const T& v)   { return P::assign(v); }
  matrix& assign(const T& v=T())  { return P::assign(v); }
  matrix& zerorow(const size_t r) { return P::zerorow(r); }

  size_t size()                const { return Nr*Nc; }
  size_t size(const size_t& d) const { return (d==0? Nr : (d==1? Nc : 0)); }
  matrix& clear() { Nr = Nc = 0; }
  matrix& resize(size_t _Nr, size_t _Nc) { Nr=_Nr; Nc=_Nc; return *this; }

  // members
  size_t Nr;  // number of rows
  size_t Nc;  // ... columns

};


/**
 * dense matrix implementation, std::vector< std::vector< T > > based
 */
template< typename T >
struct matrix_dense_vv : matrix< T,matrix_dense_vv< T > > {
  typedef matrix< T,matrix_dense_vv< T > > P;

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }

  matrix_dense_vv& assign(const T& v=T())  { if (P::size()) a.assign(P::Nr,std::vector< T >(P::Nc,v)); return *this; }
  matrix_dense_vv& zerorow(const size_t r) { if (P::size()) a[r].assign(P::Nc,T()); return *this; }
  void print(std::ostream& out) {
    bool firstline = false;
    BOOST_FOREACH(const std::vector< T >& r, a) {
      out << (firstline++? "\n  ":"[ ");
      std::copy(r.begin(),r.end(), std::ostream_iterator< T >(out," "));
    }
    out << ']' << std::endl;
  }

  //size_t size()                { return P::size(); }
  //size_t size(const size_t& d) { return P::size(d); }
  matrix_dense_vv& clear() { a.clear(); return *this; }
  matrix_dense_vv& resize(size_t _Nr, size_t _Nc) { P::resize(_Nr,_Nc); return assign(T()); }

  // members
  std::vector< std::vector< T > > a;
};


/**
 * dense matrix implementation, T[][] array based
 */
template< typename T >
struct matrix_dense_aa : matrix< T,matrix_dense_aa< T > > {
  typedef matrix< T,matrix_dense_aa< T > > P;

  // destructor
  ~matrix_dense_aa() { matrix_dense_aa::clear(); }

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }

  matrix_dense_aa& assign(const T& v=T())  {
    for (size_t r=0; r<P::Nr; ++r)
      for (size_t c=0; c<P::Nc; ++c)
        a[r][c] = v;
    return *this;
  }
  matrix_dense_aa& zerorow(const size_t r) { for (size_t c=0; c<P::Nc; ++c) a[r][c] = T(); return *this; }
  void print(std::ostream& out) {
    for (size_t r=0; r<P::Nr; ++r) {
      out << (r? "\n  ":"[ ");
      std::copy(&(a[r][0]),&(a[r][0])+P::Nc, std::ostream_iterator< T >(out," "));
    }
    out << ']' << std::endl;
  }

  matrix_dense_aa& clear() {
    if (P::size()) {
      delete[] a[0];
      delete[] a;
    }
    P::clear();
  }
  matrix_dense_aa& resize(size_t Nrows, size_t Ncolumns) {
    if (Nrows && Ncolumns) {
      P::resize(Nrows,Ncolumns);
      a    = new T*[ P::Nr ];
      a[0] = new T [ P::Nr * P::Nc ];
      for (size_t r=1; r<P::Nr; ++r)
        a[r] = a[r-1] + P::Nc;
      assign();
    }
  }

  // members
  T **a;

};


/**
 * sparse matrix implementation, in compressed sparse rows (CSR) format
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

  // interfacing
  const T& operator()(const size_t r, const size_t c) const { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }
        T& operator()(const size_t r, const size_t c)       { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }

  void clear(const T& v=T())   { for (int i=0; i<nnz; ++i) a[i] = v; }
  void zerorow(const size_t r) { for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k) a[k] = T(); }

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
template< typename T, typename MATRIX >
class linearsystem :
    public common::Action {

 public:
  // framework interfacing
  linearsystem(const std::string& name) : common::Action(name) {}
  virtual ~linearsystem() { clear(); }
  void execute() { solve(); }

 public:
  // linear system solver interfacing
  const T& A(const size_t r, const size_t c) const { return m_A(r,c); }
        T& A(const size_t r, const size_t c)       { return m_A(r,c); }
  const T& b(const size_t c) const { return m_b[c]; }
        T& b(const size_t c)       { return m_b[c]; }
  const T& x(const size_t r) const { return m_x[r]; }
        T& x(const size_t r)       { return m_x[r]; }

  linearsystem& zerorow(const size_t r) { m_A.zerorow(r); m_b(r) = 0; }
  linearsystem& assign(const T& v=T()) {
    if (size()) {
      m_A.assign(v);
      m_b.assign(m_b.size(),v);
      m_x.assign(m_x.size(),v);
    }
    return *this;
  }
  virtual linearsystem& solve() = 0;

  size_t size()                const { return m_A.size(0)*m_A.size(1); }
  size_t size(const size_t& d) const { return m_A.size(d); }
  linearsystem& clear() {
    m_A.clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }
  linearsystem& resize(size_t Nequations, size_t Nvariables, const T& v=T()) {
    if (!Nequations || !Nvariables)
      return clear();
    m_A.resize(Nequations,Nvariables);
    m_b.resize(Nequations);
    m_x.resize(Nvariables);
    return assign(v);
  }

  template< typename T2, typename MATRIX2 >
  friend std::ostream& operator<<(std::ostream&, linearsystem< T2, MATRIX2>&);

 protected:
  // member variables
  MATRIX           m_A;
  std::vector< T > m_b;
  std::vector< T > m_x;

};


template< typename T, typename MATRIX >
std::ostream& operator<<(std::ostream& out, linearsystem< T, MATRIX >& lss) {
  out << "linearsystem::A("
    << "Nequations=" << lss.m_A.size(0) << ","
    << "Nvariables=" << lss.m_A.size(1) << ","
    << "Nentries="   << lss.m_A.size()  << "):\n";
  lss.m_A.print(out);

  out << "linearsystem::b(Nequations=" << lss.m_b.size() << "):\n[ ";
  std::copy(lss.m_b.begin(), lss.m_b.end(), std::ostream_iterator< T >(out," "));
  out << "]\n"
      << "linearsystem::x(Nvariables=" << lss.m_x.size() << "):\n[ ";
  std::copy(lss.m_x.begin(), lss.m_x.end(), std::ostream_iterator< T >(out," "));
  return out << "]";
}


}  // namespace lss
}  // namespace cf3


#endif

