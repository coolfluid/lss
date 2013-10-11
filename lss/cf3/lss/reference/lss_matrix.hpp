// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_lss_matrix_h
#define cf3_lss_lss_matrix_h


#include <vector>
#include <stack>
#include "boost/lexical_cast.hpp"

#include "common/Signal.hpp"
#include "common/Action.hpp"
#include "lss_index.hpp"


namespace cf3 {
namespace lss {
namespace lss_matrix {


/**
 * description of a matrix
 * only describes an interface, storage is provided by the derived structs.
 * (T: storage type, SPECIAL: specialization class)
 * it is based on parent class template extension (CRTP) to simulate virtual
 * methods (not a true polymorphic type) so it can't be used with factories.
 */
template< typename T, class SPECIAL >
struct matrix
{
  matrix() : m_zero(T()) { clear(); }

  // matrix interfacing (derived matrices should have these defined)
  const T& operator()(const size_t r, const size_t c) const { return m_zero; }
        T& operator()(const size_t r, const size_t c)       { return m_zero; }

  size_t size()                const { return Nr*Nc; }
  size_t size(const size_t& d) const { return (d==0? Nr : (d==1? Nc : 0)); }

  matrix& assign    (const T& v=T()) { return *this; }
  matrix& operator= (const T& v) { return SPECIAL::assign(v); }
  matrix& resize    (size_t r, size_t c, const T& v=T()) { Nr=r; Nc=c; return *this; }
  matrix& zerorow   (const size_t r) { return *this; }
  matrix& clear() { Nr=Nc=0; return *this; }

  void output(std::ostream& out) const { SPECIAL::output(out); }

  T m_zero;
  size_t Nr;  // number of rows
  size_t Nc;  // ... columns
};


#if 0
/**
 * dense matrix implementation, T[] array based
 */
template< typename T >
struct dense_matrix_a :
    matrix< T,dense_matrix_a< T > > {
  typedef matrix< T,dense_matrix_a< T > > P;

  // destructor
  ~dense_matrix_a() { dense_matrix_a::clear(); }

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[P::Nc*r+c]; }
        T& operator()(const size_t r, const size_t c)       { return a[P::Nc*r+c]; }
  const T& operator()(const size_t i) const { return a[i]; }
        T& operator()(const size_t i)       { return a[i]; }

  dense_matrix_a& resize(const size_t r, const size_t c, const T& v=T()) {
    if (r && c) {
      P::resize(r,c);
      a = new T [ P::Nr * P::Nc ];
      assign(v);
    }
  }

  dense_matrix_a& zerorow(const size_t r) {
    for (size_t i=P::Nc*r; i<P::Nc*(r+1); ++i)
      a[i] = T();
    return *this;
  }

  dense_matrix_a& assign(const T& v=T())  {
    for (size_t i=0; i<P::Nr*P::Nc; ++i)
      a[i] = v;
    return *this;
  }

  dense_matrix_a& clear() {
    if (P::size())
      delete[] a;
    P::clear();
  }

  void output(std::ostream& out) const {
    out << "[ ";
    for (size_t r=0; r<P::Nr; ++r) {
      std::copy(&(a[P::Nc*r]),&(a[P::Nc*(r+1)]), std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  T *a;
};
#endif


#if 0
/**
 * dense matrix implementation, T[][] array based
 * (internally, it is stored row-oriented)
 */
template< typename T >
struct dense_matrix_aa :
    matrix< T,dense_matrix_aa< T > > {
  typedef matrix< T,dense_matrix_aa< T > > P;

  // destructor
  ~dense_matrix_aa() { dense_matrix_aa::clear(); }

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }
  const T& operator()(const size_t i) const { const size_t b(P::size(1)); return a[i/b][i%b]; }
        T& operator()(const size_t i)       { const size_t b(P::size(1)); return a[i/b][i%b]; }

  dense_matrix_aa& resize(const size_t r, const size_t c, const T& v=T()) {
    if (r && c) {
      P::resize(r,c);
      a    = new T*[ P::Nr ];
      a[0] = new T [ P::Nr * P::Nc ];
      for (size_t r=1; r<P::Nr; ++r)
        a[r] = a[r-1] + P::Nc;
      assign(v);
    }
  }

  dense_matrix_aa& zerorow(const size_t r) {
    for (size_t c=0; c<P::Nc; ++c)
      a[r][c] = T();
    return *this;
  }

  dense_matrix_aa& assign(const T& v=T())  {
    for (size_t r=0; r<P::Nr; ++r)
      for (size_t c=0; c<P::Nc; ++c)
        a[r][c] = v;
    return *this;
  }

  dense_matrix_aa& clear() {
    if (P::size()) {
      delete[] a[0];
      delete[] a;
    }
    P::clear();
  }

  void output(std::ostream& out) const {
    out << "[ ";
    for (size_t r=0; r<P::Nr; ++r) {
      std::copy(&(a[r][0]),&(a[r][0])+P::Nc, std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  T **a;
};
#endif


#if 1
/**
 * dense matrix implementation, std::vector< T > based
 */
template< typename T >
struct dense_matrix_v :
    matrix< T,dense_matrix_v< T > > {
  typedef matrix< T,dense_matrix_v< T > > P;

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[P::Nc*r+c]; }
        T& operator()(const size_t r, const size_t c)       { return a[P::Nc*r+c]; }
  const T& operator()(const size_t i) const { return a[i]; }
        T& operator()(const size_t i)       { return a[i]; }

  dense_matrix_v& resize  (const size_t r, const size_t c, const T& v=T()) { P::resize(r,c); return assign(v); }
  dense_matrix_v& zerorow (const size_t r) { if (P::size()) std::fill_n(&a[P::Nc*r],P::Nc,T()); return *this; }
  dense_matrix_v& assign  (const T& v=T()) { if (P::size()) a.assign(P::Nr*P::Nc,v); return *this; }
  dense_matrix_v& clear() { a.clear(); P::clear(); return *this; }

  void output(std::ostream& out) const {
    out << "[ ";
    for (size_t r=0; r<P::Nr; ++r) {
      std::copy(&a[P::Nc*r],&a[P::Nc*(r+1)], std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  std::vector< T > a;
};
#endif


#if 0
/**
 * dense matrix implementation, std::vector< std::vector< T > > based
 * (internally, storage is row-oriented)
 */
template< typename T >
struct dense_matrix_vv :
    matrix< T,dense_matrix_vv< T > > {
  typedef matrix< T,dense_matrix_vv< T > > P;

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }
  const T& operator()(const size_t i) const { const size_t b(P::size(0)); return a[i/b][i%b]; }
        T& operator()(const size_t i)       { const size_t b(P::size(0)); return a[i/b][i%b]; }

  dense_matrix_vv& resize  (const size_t r, const size_t c, const T& v=T()) { P::resize(r,c); return assign(v); }
  dense_matrix_vv& zerorow (const size_t r) { if (P::size()) a[r].assign(P::Nc,T());                    return *this; }
  dense_matrix_vv& assign  (const T& v=T()) { if (P::size()) a.assign(P::Nr,std::vector< T >(P::Nc,v)); return *this; }
  dense_matrix_vv& clear() { a.clear(); P::clear(); return *this; }

  void output(std::ostream& out) const {
    out << "[ ";
    BOOST_FOREACH(const std::vector< T >& r, a) {
      std::copy(r.begin(),r.end(), std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  std::vector< std::vector< T > > a;
};
#endif


#if 1
/**
 * sparse matrix implementation, in compressed sparse row (CSR) format
 * note: BASE={0,1}: {0,1}-based indexing (other values probably don't work)
 */
template< typename T, int BASE >
struct sparse_matrix_csr : matrix< T,sparse_matrix_csr< T,BASE > > {
  typedef matrix< T,sparse_matrix_csr< T,BASE > > P;

  // cons/destructor
  sparse_matrix_csr() : P(), zero(T()), nnz(0), nnu(0) {}
  ~sparse_matrix_csr() {
    if (nnz) {
      delete[] a;
      delete[] ja;
      delete[] ia;
    }
  }

  // interfacing
  const T& operator()(const size_t r, const size_t c) const { const int i=getindex(r,c); return (i<0? zero : a[i]); }
        T& operator()(const size_t r, const size_t c)       { const int i=getindex(r,c); return (i<0? zero : a[i]); }
  const T& operator()(const size_t i) const { return a[i]; }
        T& operator()(const size_t i)       { return a[i]; }

  void clear(const T& v=T())   { for (int i=0; i<nnz; ++i) a[i] = v; }
  void zerorow(const size_t r) { for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k) a[k] = T(); }

  void resize(size_t r, size_t c) { P::resize(r,c); }
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

  void output(std::ostream& out) const {
    out << "[ ";
    out << "(CSR sparse matrix) ";
    //FIXME: implement sparse matrix output
//    BOOST_FOREACH(const std::vector< T >& r, a) {
//      std::copy(r.begin(),r.end(), std::ostream_iterator< T >(out,", "));
//      out << "\n  ";
//    }
    out << ']';
  }

  // members
  T zero;
  int *ia;
  int *ja;
  T   *a;
  int nnz;
  int nnu;

};
#endif


// NOTE: missing MSR and VBR (not very useful?)


}  // namespace lss_matrix
}  // namespace lss
}  // namespace cf3


#endif

