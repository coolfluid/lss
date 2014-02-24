// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_matrix_hpp
#define cf3_lss_matrix_hpp


#include <cmath>
#include <iterator>
#include <numeric>

#include "common/Log.hpp"

#include "utilities.hpp"


namespace cf3 {
namespace lss {


/* -- matrix helper definitions --------------------------------------------- */

/// @brief Matrix print level
enum print_t { print_auto=0, print_size, print_signs, print_full, print_file };


/// @brief Matrix orientations
enum orientation_t { sort_by_column=0, sort_by_row=1 };


/// @brief Sparse/coordinate matrix entry
template< typename T >
struct coord_t : std::pair< idx_t, T > {
  coord_t(const idx_t& _idx, const T& _value) : std::pair< idx_t, T >(_idx,_value) {}
};


/// @brief Sparse/coordinate matrix sorter (functor, to be specialized)
template< typename _Key, int ORIENT >
struct sort_t {
  bool operator()(const _Key& a, const _Key& b) const {
    throw std::runtime_error("sort_t: not implemented.");
    return true;
  }
};


/// @brief Sparse/coordinate matrix sorter, specialized to column-oriented
template< typename _Key >
struct sort_t< _Key, sort_by_column >
{
  bool operator()(const _Key& a, const _Key& b) const {
    return (a.first.j<b.first.j? true  :
           (a.first.j>b.first.j? false :
           (a.first.i<b.first.i)));
  }
};


/// @brief Sparse/coordinate matrix sorter, specialized to row-oriented
template< typename _Key >
struct sort_t< _Key, sort_by_row >
{
  bool operator()(const _Key& a, const _Key& b) const {
    return (a.first.i<b.first.i? true  :
           (a.first.i>b.first.i? false :
           (a.first.j<b.first.j) ));
  }
};


/* -- matrix interface  ----------------------------------------------------- */

/**
 * @brief Generic matrix basic interface, based on class template implementation
 * (like profile pattern) to implement functionality (this is NOT polymorphism).
 * T: storage type
 * Impl: implementation type
 */
template< typename T, class IMPL >
struct matrix
{
  // constructor
  matrix() :
    m_zero(std::numeric_limits< T >::signaling_NaN()),
    m_size()
  {}
  virtual ~matrix() {}

  // -- functionality in remote implementation

  virtual matrix& initialize(
      const size_t& i,
      const size_t& j,
      const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >() ) = 0;

  virtual matrix& initialize(const std::vector< double >& _vector) {
    if (_vector.size()==1)
      return operator=(_vector[0]);
    else if (m_size.i*m_size.j!=_vector.size())
      throw std::runtime_error("matrix: assignment not consistent with current size.");
    for (size_t i=0, k=0; i<size(0); ++i)
      for (size_t j=0; j<size(1); ++j, ++k)
        operator()(i,j) = (type_is_equal< T, double >()? (const T) _vector[k] : static_cast< const T >(_vector[k]) );
    return *this;
  }

  virtual matrix& initialize(const std::string& _fname) {
    using namespace std;
    clear();
    m_size.invalidate();
    try {

      ifstream f(_fname.c_str());
      if (!f) throw runtime_error("matrix: cannot open file.");
      const bool hasdot(string("."+_fname).find_last_of("."));


      // read format: MatrixMarket (*.mtx)
      if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {


        // read matrix size and properties, and allocate
        MatrixMarket::typecode_t t;
        int nnz(0);
        if (!MatrixMarket::read_banner(f,t))                    throw runtime_error("matrix: MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
        if (!MatrixMarket::read_size(f,m_size.i,m_size.j,nnz))  throw runtime_error("matrix: MatrixMarket: invalid matrix/array size.");
        if (!t.is_general())                                    throw runtime_error("matrix: MatrixMarket: only \"general\" matrices supported.");

        initialize(m_size.i,m_size.j);

        // read into row/column-oriented dense matrix, line by line
        idx_t ij(1,1);
        string line;
        while (getline(f,line)) {
          if (line.find_first_of("%")==0) {}
          else {
            T& v(operator()(ij.i-1,ij.j-1));
            istringstream strm(line);

            double
                a(0.),
                b(0.);
            t.is_sparse()?  strm >> ij.i >> ij.j : strm;
            t.is_complex()? strm >> a >> b :
            t.is_real()?    strm >> a :
            t.is_integer()? strm >> a :
                            strm;

            if (t.is_pattern()) {
              v = T();
            }
            else if (type_is_equal< T, zdouble >()) {
              ((zdouble&) v).real(a);
              ((zdouble&) v).imag(b);
            }
            else if (type_is_equal< T, zfloat >()) {
              ((zfloat&) v).real(static_cast< float >(a));
              ((zfloat&) v).imag(static_cast< float >(b));
            }
            else {
              v = static_cast< T >(a);
            }

            if (t.is_dense() && (++ij.i>m_size.i)) {
              ij.i = 1;
              ij.j++;
            }
          }
        }


      }
      // read format: CSR (*.csr)
      else if (hasdot && _fname.substr(_fname.find_last_of("."))==".csr") {


        // read matrix size, and allocate
        string line;
        while (!m_size.is_valid_size()) {
          if (getline(f,line) && line.find_first_of("%")!=0)
            istringstream(line) >> m_size.i >> m_size.j;
        }

        initialize(m_size.i,m_size.j);

        // read rows and column indices (base 0), then fill in non-zeros
        vector< int > ia;
        ia.reserve(m_size.i+1);
        for (int i; ia.size()<m_size.i+1 && f >> i;)
          ia.push_back(i);

        vector< int > ja;
        ja.reserve(ia.back()-ia.front());
        for (int i; ja.size()<ia.back()-ia.front() && f >> i;)
          ja.push_back(i);

        for_each(ja.begin(),ja.end(),base_conversion_t(-ia.front()));
        for_each(ia.begin(),ia.end(),base_conversion_t(-ia.front()));

        for (size_t i=0; i<m_size.i; ++i)
          for (int k=ia[i]; k<ia[i+1] && f; ++k)
            f >> operator()(i,ja[k]);


      }


      // read format: file format not detected
      else {
        throw runtime_error("matrix: file format not detected.");
      }
      f.close();


    }
    catch (const runtime_error& e) {
      throw runtime_error("matrix: cannot read file. " + string(e.what()));
    }
    return *this;
  }

  matrix& clear() {
    m_zero = std::numeric_limits< T >::signaling_NaN();
    m_size.clear();
    return *this;
  }

  virtual matrix& operator=(const double& _value) = 0;
  matrix& operator=(const matrix& _other) { return IMPL::operator=(_other); }
  matrix& zerorow(const size_t& i)        { return IMPL::zerorow(i); }
  matrix& sumrows(const size_t& i, const size_t& isrc) { return IMPL::sumrows(i,isrc); }
  T sumrows(const size_t& j=0) const { return IMPL::sumrows(j); }
  virtual T norm1(const size_t& j=0) { T n(0); for (size_t i=0; i<m_size.i; ++i) n += std::abs(operator()(i,j)); return n; }
  virtual T norm2(const size_t& j=0) { T n(0); for (size_t i=0; i<m_size.i; ++i) { const T a(operator()(i,j)); n += std::abs(a*a); } return std::sqrt(n); }
  virtual T normi(const size_t& j=0) { T n(0); for (size_t i=0; i<m_size.i; ++i) n = std::max(std::abs(n),std::abs(operator()(i,j))); return n; }

  // -- intrinsic functionality

  void swap(matrix& _other) {
    std::swap(_other.m_size,m_size);
    std::swap(_other.m_zero,m_zero);
  }

  virtual void print(std::ostream& o, const print_t& l=print_auto) const {

    const double eps = 1.e3*static_cast< double >(std::abs(std::numeric_limits< T >::epsilon()));
    const print_t lvl(l? std::max(print_size,std::min(l,print_file)) :
                     (m_size.i>100 || m_size.j>100? print_size  :
                     (m_size.i> 10 || m_size.j> 10? print_signs :
                                                    print_full )));
    std::string row_s;
    std::ostringstream ss;

    switch (lvl) {
      case print_size:
        o << "(" << size(0) << 'x' << size(1) << ") [ ... ]";
        break;

      case print_signs:
        o << "(" << size(0) << 'x' << size(1) << ") [";
        for (size_t i=0; i<size(0); ++i) {
          row_s.resize(size(1));
          for (size_t j=0; j<size(1); ++j) {
            const T& v(operator()(i,j));
            const double a(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).real()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).real() : (const double&) v );
            row_s[j] = type_is_complex< T >()? (static_cast< double >(std::abs(v))>eps? '*':'0') :
                       a> eps? '+' :
                       a<-eps? '-' :
                               '0' ;
          }
          o << "\n  " << row_s;
        }
        o << " ]";
        break;

      case print_full:
        o << "(" << size(0) << 'x' << size(1) << ") [";
        for (size_t i=0; i<size(0); ++i) {
          ss.str("");
          for (size_t j=0; j<size(1); ++j)
            ss << operator()(i,j) << ", ";
          o << "\n  " << ss.str();
        }
        o << " ]";
        break;

      case print_file:
        CFinfo << (typeid(T).name()) << CFendl;
        o << "%%MatrixMarket matrix array"
          << (type_is_complex< T >()? " complex":" real")
          << " general\n"
          << size(0) << ' ' << size(1) << '\n';
        for (size_t j=0; j<size(1); ++j)
          for (size_t i=0; i<size(0); ++i) {
            const T& v(operator()(i,j));
            const double
                a(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).real()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).real() : (const double&) v ),
                b(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).imag()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).imag() : double() );
            type_is_complex< T >()? o << a << ' ' << b << '\n' :
                                    o << a << '\n';
          }
        break;

      case print_auto:
      default:
        break;
    }
  }

  size_t size(const size_t& _d) const {
    return (_d==0? m_size.i : (_d==1? m_size.j : 0));
  }

  // indexing (defer to implementation)
  virtual const T& operator()(const size_t& i, const size_t& j=0) const = 0;
  virtual       T& operator()(const size_t& i, const size_t& j=0)       = 0;

  // members (implementation should maintain this)
  T m_zero;
  idx_t m_size;
};


/* -- matrix implementations ------------------------------------------------ */

/**
 * @brief Dense matrix, stored in row or column oriented vector< vector< T > >
 * (can work as column vectors if that is the intention, and sometimes it is...)
 * T: storage type
 * ORIENT: if storage is row or column (default) oriented
 */
template< typename T, int ORIENT=sort_by_column >
struct dense_matrix_vv :
    matrix< T, dense_matrix_vv< T > >
{
  typedef matrix< T, dense_matrix_vv< T > > matrix_base_t;
  using matrix_base_t::size;

  // initializations

  dense_matrix_vv& initialize(
      const size_t& i,
      const size_t& j,
      const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >() ) {
    if (idx_t(i,j).is_valid_size()) {
      if (i!=size(0) && j!=size(1))
        a.clear();
      if (i*j)
        a.assign(ORIENT? i:j, std::vector< T >(ORIENT? j:i, T()));
      matrix_base_t::m_size = idx_t(i,j);
    }
    return *this;
  }

  dense_matrix_vv& initialize(const std::vector< double >& _vector) {
    if (_vector.size()==1)
      return operator=(_vector[0]);
    matrix_base_t::initialize(_vector);
    return *this;
  }

  dense_matrix_vv& initialize(const std::string& _fname) {
    matrix_base_t::initialize(_fname);
    return *this;
  }

  // assignments

  dense_matrix_vv& operator=(const dense_matrix_vv& _other) {
    a = _other.a;
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
    return *this;
  }

  dense_matrix_vv& operator=(const double& _value) {
    if (size(0)*size(1))
      a.assign(ORIENT? size(0):size(1),std::vector< T >(
               ORIENT? size(1):size(0),static_cast< T >(_value) ));
    return *this;
  }

  dense_matrix_vv& clear() {
    matrix_base_t::clear();
    a.clear();
    return *this;
  }

  dense_matrix_vv& zerorow(const size_t& i) {
    if (i>=size(0))
      throw std::runtime_error("dense_matrix_vv: row index outside bounds.");
    if (ORIENT) a[i].assign(size(1),T());
    else {
      for (size_t j=0; j<size(1); ++j)
        a[j][i] = T();
    }
    return *this;
  }

  dense_matrix_vv& sumrows(const size_t& i, const size_t& isrc) {
    if (std::max(i,isrc)>=size(0))
      throw std::runtime_error("dense_matrix_vv: row index(es) outside bounds.");
    for (size_t j=0; j<size(1); ++j)
      operator()(i,j) += operator()(isrc,j);
    return *this;
  }

  T sumrows(const size_t& j=0) const {
    if (j>=size(1))
      throw std::runtime_error("dense_matrix_vv: column index outside bounds.");
    T s(0);
    if (ORIENT) for (size_t i=0; i<size(0); ++i) s+=a[i][j];
    else        s = std::accumulate(a[j].begin(),a[j].end(),T());
    return s;
  }

  dense_matrix_vv& swap(dense_matrix_vv& other) {
    other.a.swap(a);
    matrix_base_t::swap(other);
    return *this;
  }

  // indexing
  const T& operator()(const size_t& i, const size_t& j=0) const { return ORIENT? a[i][j]:a[j][i]; }
        T& operator()(const size_t& i, const size_t& j=0)       { return ORIENT? a[i][j]:a[j][i]; }

  // storage
  std::vector< std::vector< T > > a;
};


/**
 * @brief Dense matrix, stored in row or column oriented vector< vector< T > >
 * (can work as column vectors if that is the intention, and sometimes it is...)
 * T: storage type
 * ORIENT: if storage is row or column (default) oriented
 */
template< typename T, int ORIENT=sort_by_column >
struct dense_matrix_v :
    matrix< T, dense_matrix_v< T > >
{
  typedef matrix< T, dense_matrix_v< T > > matrix_base_t;
  using matrix_base_t::size;

  // initializations

  dense_matrix_v& initialize(
      const size_t& i,
      const size_t& j,
      const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >() ) {
    if (idx_t(i,j).is_valid_size()) {
      if (i!=size(0) && j!=size(1))
        a.clear();
      if (i*j)
        a.assign(i*j,T());
      matrix_base_t::m_size = idx_t(i,j);
    }
    return *this;
  }

  dense_matrix_v& initialize(const std::vector< double >& _vector) {
    if (_vector.size()==1)
      return operator =(_vector[0]);
    matrix_base_t::initialize(_vector);
    return *this;
  }

  dense_matrix_v& initialize(const std::string& _fname) {
    matrix_base_t::initialize(_fname);
    return *this;
  }

  // assignments

  dense_matrix_v& operator=(const dense_matrix_v& _other) {
    a = _other.a;
    matrix_base_t::m_size = _other.dense_matrix_v::matrix_base_t::m_size;
    return *this;
  }

  dense_matrix_v& operator=(const double& _value) {
    if (size(0)*size(1))
      a.assign(size(0)*size(1),static_cast< T >(_value));
    return *this;
  }

  void swap(dense_matrix_v& other) {
    other.a.swap(a);
    matrix_base_t::swap(other);
  }

  // clearing

  dense_matrix_v& clear() {
    matrix_base_t::clear();
    a.clear();
    return *this;
  }

  dense_matrix_v& zerorow(const size_t& i) {
    if (i>=size(0))
      throw std::runtime_error("dense_matrix_v: row index outside bounds.");
    if (ORIENT) std::fill_n(a.begin()+i*size(1),size(1),T());
    else {
      for (size_t j=0, k=i; j<size(1); ++j, k+=size(0))
        a[k] = T();
    }
    return *this;
  }

  dense_matrix_v& sumrows(const size_t& i, const size_t& isrc) {
    if (std::max(i,isrc)>=size(0))
      throw std::runtime_error("dense_matrix_v: row index(es) outside bounds.");
    for (size_t j=0; j<size(1); ++j)
      operator()(i,j) += operator()(isrc,j);
    return *this;
  }

  T sumrows(const size_t& j=0) const {
    if (j>=size(1))
      throw std::runtime_error("dense_matrix_v: column index outside bounds.");
    if (!ORIENT)
      return std::accumulate(a.begin()+size(0)*(j),a.begin()+size(0)*(j+1),T());
    T s(0);
      for (size_t i=0; i<size(0); ++i)
        s += operator()(i,j);
    return s;
  }

  // indexing
  const T& operator()(const size_t& i, const size_t& j=0) const { return ORIENT? a[i*matrix_base_t::m_size.j+j]:a[j*matrix_base_t::m_size.i+i]; }
        T& operator()(const size_t& i, const size_t& j=0)       { return ORIENT? a[i*matrix_base_t::m_size.j+j]:a[j*matrix_base_t::m_size.i+i]; }

  // storage
  std::vector< T > a;
};


/**
 * @brief Sparse matrix: coordinate matrix with particular ordering, to allow
 * compression of rows/columns as necessary. It suffers from schizophrenia,
 * alternating between compressed/uncompressed personalities to provide dynamic
 * entries insertion/retrieval.
 * T: storage type
 * ORIENT: if storage is row (default) or column oriented
 * BASE: column & row numbering base (0 or 1, other values won't work)
 */
template< typename T, int ORIENT=sort_by_row, int BASE=0 >
struct sparse_matrix :
  matrix< T,sparse_matrix< T,ORIENT,BASE > >
{
  // utility definitions
  typedef matrix< T,sparse_matrix< T,ORIENT,BASE > > matrix_base_t;

  // uncompressed/compressed matrix structures definitions
  typedef std::set< coord_t<T>, sort_t< coord_t<T>, ORIENT > > matrix_uncompressed_t;
  struct matrix_compressed_t {
    void clear() {
      nnu = nnz = 0;
      ia.clear();
      ja.clear();
      a.clear();
    }
    matrix_compressed_t& swap(matrix_compressed_t& _other) {
      swap(nnu,_other.nnu);
      swap(nnu,_other.nnz);
      ia.swap(_other.ia);
      ja.swap(_other.ja);
      a .swap(_other.a);
    }
    int nnu, nnz;               // number of rows/nonzeros
    std::vector< int > ia, ja;  // rows/column indices
    std::vector< T > a;         // values
  };

  // constructor
  sparse_matrix() : matrix_base_t() {
    if (BASE!=0 && BASE!=1)
      throw std::logic_error("sparse_matrix: indexing base should be 0 or 1.");
  }

  // initializations

  sparse_matrix& initialize(
      const size_t& i,
      const size_t& j,
      const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >() ) {
    if (idx_t(i,j).is_valid_size()) {
      matu.clear();
      matc.clear();
      matrix_base_t::m_size = idx_t(i,j);
      if (_nnz.size() && ORIENT) {

        // build (already compressed) row and column indices, and allocate values
        matc.nnu = static_cast< int >(_nnz.size());
        matc.ia.reserve(matc.nnu+1);
        matc.ia.push_back(BASE);
        for (size_t r=0; r<_nnz.size(); ++r)
          matc.ia.push_back(matc.ia.back() + static_cast< int >(_nnz[r].size()));

        matc.nnz = matc.ia.back()-BASE;
        matc.ja.resize(matc.nnz);
        for (size_t r=0, k=0; r<_nnz.size(); ++r)
          for (size_t c=0; c<_nnz[r].size(); ++c)
            matc.ja[k++] = static_cast< int >(_nnz[r][c])+BASE;

        matc.a.assign(matc.nnz,T());

      }
      else if (_nnz.size()) {

        for (size_t r=0; r<_nnz.size(); ++r)
          for (std::vector< size_t >::const_iterator c=_nnz[r].begin(); c!=_nnz[r].end(); ++c)
            operator()(r,*c) = T();
        compress();

      }
    }
    else {
      CFwarn << "sparse_matrix: invalid size: (" << i << ',' << j << ')' << CFendl;
    }
    return *this;
  }

  sparse_matrix& initialize(const std::vector< double >& _vector) {
    if (_vector.size()==1)
      return operator =(_vector[0]);
    matrix_base_t::initialize(_vector);
    compress();
    return *this;
  }

  sparse_matrix& initialize(const std::string& _fname) {
    matrix_base_t::initialize(_fname);
    compress();
    return *this;
  }

  sparse_matrix& clear() {
    matrix_base_t::clear();
    matu.clear();
    matc.clear();
    return *this;
  }

  sparse_matrix& operator=(const double& _value) {
    CFdebug << "sparse_matrix: assigning a value only affects populated entries." << CFendl;
    const T value = static_cast< T >(_value);
    std::fill(matc.a.begin(),matc.a.end(),value);
    for (typename matrix_uncompressed_t::iterator it = matu.begin(); it!=matu.end(); ++it)
      const_cast< T& >(it->second) = value;
    return *this;
  }

  sparse_matrix& operator=(const sparse_matrix& _other) {
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
    matu = _other.matu;
    matc = _other.matc;
    return *this;
  }

  sparse_matrix& zerorow(const size_t& i) {
    if (i>=matrix_base_t::m_size.i)
      throw std::runtime_error("sparse_matrix: row index out of bounds.");
    if (is_compressed()) {
      for (int k=matc.ia[i]-BASE; ORIENT && k<matc.ia[i+1]-BASE; ++k)
        matc.a[k] = T();
      for (int k=0; !ORIENT && k<matc.nnz; ++k)
        if (matc.ia[k]-BASE==(int) i)
            matc.a[k] = T();
    }
    else {
      for (typename matrix_uncompressed_t::iterator it = matu.begin(); it!=matu.end(); ++it)
        if (it->first.i==(size_t) i)
          const_cast< T& >(it->second) = T();
    }
    return *this;
  }

  sparse_matrix& sumrows(const size_t& i, const size_t& isrc) {
    if (std::max(i,isrc)>=matrix_base_t::m_size.i)
      throw std::runtime_error("sparse_matrix: row index(es) outside bounds.");

    matrix_uncompressed_t rowu;  // (row buffer, in case matrix is compressed)
    if (is_compressed()) {
      for (int k=matc.ia[isrc]-BASE; ORIENT && k<matc.ia[isrc+1]-BASE; ++k)
        rowu.insert(coord_t<T>(idx_t(isrc,matc.ja[k]-BASE), matc.a[k] ));
      for (int j=0; !ORIENT && j<matc.nnu; ++j)
        for (int k=matc.ja[j]-BASE; k<matc.ja[j+1]-BASE; ++k)
          if (matc.ia[k]-BASE==isrc)
            rowu.insert(coord_t<T>(idx_t(isrc,j), matc.a[k] ));
    }

    matrix_uncompressed_t& mat = is_compressed()? rowu : matu;
    for (typename matrix_uncompressed_t::const_iterator it = mat.begin(); it!=mat.end(); ++it)
      if (it->first.i==isrc)
        operator()(i,it->first.j) += it->second;
    return *this;
  }

  T sumrows(const size_t& j=0) const {
    if (j>=this->size(1))
      throw std::runtime_error("sparse_matrix: column index outside bounds.");

    T s(0);
    if (!is_compressed()) {
      // uncompressed
      for (typename matrix_uncompressed_t::const_iterator it=matu.begin(); it!=matu.end(); ++it)
        if (it->first.j==j)
          s += it->second;
    }
    else if (ORIENT) {
      // compressed & sorted by row
      for (int i=0; i<matc.nnu; ++i)
        for (int k=matc.ia[i]-BASE; k<matc.ia[i+1]-BASE; ++k)
          if (matc.ja[k]-BASE==j) {
            s += matc.a[k];
            break;
          }
    }
    else {
      // compressed & sorted by column
      s = std::accumulate( matc.a.begin()+(matc.ja[j  ]-BASE),
                           matc.a.begin()+(matc.ja[j+1]-BASE), T() );
    }
    return s;
  }

  T norm1(const size_t& j=0) { T n(0); for (size_t i=0; i<this->size(0); ++i) n += std::abs(operator()(i,j)); return n; }
  T norm2(const size_t& j=0) { T n(0); for (size_t i=0; i<this->size(0); ++i) { const T a(operator()(i,j)); n += a*a; } return std::sqrt(n); }
  T normi(const size_t& j=0) { T n(0); for (size_t i=0; i<this->size(0); ++i) n = std::max(n,std::abs(operator()(i,j))); return n; }

  sparse_matrix& swap(sparse_matrix& _other) {
    matrix_base_t::swap(_other);
    std::swap(matu,_other.matu);
    std::swap(matc,_other.matc);
    return *this;
  }

  void print(std::ostream& o, const print_t& l=print_auto) const {
    using namespace std;
    const double eps = 1.e3*static_cast< double >(abs(numeric_limits< T >::epsilon()));
    const idx_t&  size = matrix_base_t::m_size;
    const print_t lvl(l? max(print_size,min(l,print_file)) :
                     (size.i>100 || size.j>100? print_size  :
                     (size.i> 10 || size.j> 10? print_signs :
                                                print_full )));

    if (lvl==print_size)  {
      o << "(" << size.i << 'x' << size.j << ">=" << (matc.nnz+matu.size()) << ") [ ... ]";
    }
    else {

      // (requires an uncompressed structure)
      matrix_uncompressed_t tmp;
      if (is_compressed())
        uncompress(size,tmp,matc);
      const matrix_uncompressed_t& mat(is_compressed()? tmp:matu);
      switch (lvl) {

        case print_signs:
          o << "(" << size.i << 'x' << size.j << ">=" << mat.size() << ") [ ";
          for (size_t i=0; i<size.i; ++i) {
            string row(size.j,' ');
            for (typename matrix_uncompressed_t::const_iterator it = mat.begin(); it!=mat.end(); ++it) {
              const T& v(it->second);
              const double a(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).real()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).real() : (const double&) v );
              if (it->first.i==i)
                row[ it->first.j ] = type_is_complex< T >()? (static_cast< double >(abs(v))>eps? '*':'.') :
                                       a> eps? '+' :
                                       a<-eps? '-' :
                                               '.' ;
            }
            o << "\n  " << row;
          }
          o << " ]";
          break;

        case print_full:
          o << "(" << size.i << 'x' << size.j << ">=" << mat.size() << ") [ ";
          for (size_t i=0; i<size.i; ++i) {
            vector< T >row(size.j,T());
            for (typename matrix_uncompressed_t::const_iterator it = mat.begin(); it!=mat.end(); ++it)
              if (it->first.i==i)
                row[ it->first.j ] = it->second;
            o << "\n  ";
            copy(row.begin(),row.end(),ostream_iterator< T >(o,", "));
          }
          o << " ]";
          break;

        case print_file:
          o << "%%MatrixMarket matrix coordinate"
            << (type_is_complex< T >()? " complex":" real")
            << " general\n"
            << size.i << ' ' << size.j << ' ' << mat.size() << '\n';
          for (typename matrix_uncompressed_t::const_iterator i = mat.begin(); i!=mat.end(); ++i) {
            o << (i->first.i+1) << ' '<< (i->first.j+1) << ' ';
            const T& v(i->second);
            const double
                a(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).real()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).real() : (const double&) v ),
                b(type_is_equal< T, zfloat >()? static_cast< double >(((const zfloat&) v).imag()) : type_is_equal< T, zdouble >()? ((const zdouble&) v).imag() : double() );
            type_is_complex< T >()? o << a << ' ' << b << '\n' :
                                    o << a << '\n';
          }
          break;

        case print_auto:
        case print_size:
        default:
          break;
      }

    }
  }

  // indexing
  const T& operator()(const size_t& i, const size_t& j) const {
    if (i>=matrix_base_t::m_size.i || j>=matrix_base_t::m_size.j) {
      CFwarn << "sparse_matrix: index not available, (" << i << ',' << j << ") >= (" << matrix_base_t::m_size.i << ',' << matrix_base_t::m_size.j << ")." << CFendl;
      return matrix_base_t::m_zero;
    }
    if (is_compressed()) {
      for (int k=matc.ia[i]-BASE; is_compressed() && ORIENT && k<matc.ia[i+1]-BASE; ++k)
        if (matc.ja[k]-BASE==(int) j)
          return matc.a[k];
      for (int k=matc.ja[j]-BASE; is_compressed() && !ORIENT && k<matc.ja[j+1]-BASE; ++k)
        if (matc.ia[k]-BASE==(int) i)
          return matc.a[k];
    }
    else {
      typename matrix_uncompressed_t::const_iterator it = matu.find(coord_t<T>(idx_t(i,j),T()));
      if (it!=matu.end())
        return it->second;
    }
    CFwarn << "sparse_matrix: index not found: (" << i << ',' << j << ")." << CFendl;
    return matrix_base_t::m_zero;
  }

  T& operator()(const size_t& i, const size_t& j) {
    if (i>=matrix_base_t::m_size.i || j>=matrix_base_t::m_size.j) {
      CFwarn << "sparse_matrix: index not available, (" << i << ',' << j << ") >= (" << matrix_base_t::m_size.i << ',' << matrix_base_t::m_size.j << ")." << CFendl;
      return matrix_base_t::m_zero;
    }
    if (is_compressed()) {
      for (int k=matc.ia[i]-BASE; ORIENT && k<matc.ia[i+1]-BASE; ++k)
        if (matc.ja[k]-BASE==(int) j)
          return matc.a[k];
      for (int k=matc.ja[j]-BASE; !ORIENT && k<matc.ja[j+1]-BASE; ++k)
        if (matc.ia[k]-BASE==(int) i)
          return matc.a[k];
    }
    // find/insert new entry, and structurally symmetric pair. the constness
    // removal is safe because the entry value does not change matrix ordering
    uncompress();
    std::pair< typename matrix_uncompressed_t::iterator, bool > p =
      matu.insert( coord_t<T>(idx_t(i,j),T()) );
    if (p.second && i!=j)
      matu.insert( coord_t<T>(idx_t(j,i),T()) );
    return const_cast< T& >((p.first)->second);
  }


  // compression/uncompression

  matrix_compressed_t& compress() {
    if (!is_compressed()) {
      CFinfo << "sparse_matrix::compress..." << CFendl;
      const size_t nmodif = ensure_structural_symmetry(matrix_base_t::m_size,matu);
      if (nmodif)
        CFinfo << "sparse_matrix: symmetry preserving additional entries: " << nmodif << CFendl;
      compress(matrix_base_t::m_size,matu,matc);
      matu.clear();
      CFinfo << "sparse_matrix::compress." << CFendl;
    }
    return matc;
  }

  matrix_uncompressed_t& uncompress() {
    if (is_compressed()) {
      CFinfo << "sparse_matrix::uncompress..." << CFendl;
      uncompress(matrix_base_t::m_size,matu,matc);
      matc.clear();
      const size_t nmodif = ensure_structural_symmetry(matrix_base_t::m_size,matu);
      if (nmodif)
        CFinfo << "sparse_matrix: symmetry preserving additional entries: " << nmodif << CFendl;
      CFinfo << "sparse_matrix::uncompress." << CFendl;
    }
    return matu;
  }


 private:
  // compression utilities (not to use outside this context)

  inline bool is_compressed() const { return matc.nnz; }

  static void compress(
    const idx_t& _size,
    const matrix_uncompressed_t & _u,
    matrix_compressed_t& _c)
  {
    _c.clear();
    typename matrix_uncompressed_t::const_iterator it=_u.begin();
    if (it==_u.end())
      return;

    _c.nnu = (ORIENT? _size.i:_size.j);
    _c.nnz = _u.size();
    _c.ia.reserve(ORIENT? _c.nnu+1 : _c.nnz  );
    _c.ja.reserve(ORIENT? _c.nnz   : _c.nnu+1);
    _c.a .reserve(_c.nnz);

    // row-oriented matrix
    if (ORIENT) {
      _c.ia.push_back(0);
      for (size_t count, r=0; r<_size.i; ++r) {
        for (count=0; r==(it->first.i) && it!=_u.end(); ++it, ++count) {
          _c.ja.push_back(it->first.j);
          _c.a .push_back(it->second);
        }
        _c.ia.push_back(_c.ia.back() + count);
      }
    }

    // column-oriented matrix
    else {
      _c.ja.push_back(0);
      for (size_t count, c=0; c<_size.j; ++c) {
        for (count=0; c==(it->first.j) && it!=_u.end(); ++it, ++count) {
          _c.ia.push_back(it->first.i);
          _c.a .push_back(it->second);
        }
        _c.ja.push_back(_c.ja.back() + count);
      }
    }

    if (BASE) {
      std::for_each(_c.ja.begin(),_c.ja.end(),base_conversion_t(BASE));
      std::for_each(_c.ia.begin(),_c.ia.end(),base_conversion_t(BASE));
    }
  }


  static void uncompress(
      const idx_t& _size,
      matrix_uncompressed_t & _u,
      const matrix_compressed_t& _c)
  {
    _u.clear();
    for (int i=0; ORIENT && i<_c.ia.size()-1; ++i)
      for (int k=_c.ia[i]-BASE; k<_c.ia[i+1]-BASE; ++k)
        _u.insert(coord_t<T>(idx_t(i,_c.ja[k]-BASE),_c.a[k]));
    for (int j=0; !ORIENT && j<_c.ja.size()-1; ++j)
      for (int k=_c.ja[j]-BASE; k<_c.ja[j+1]-BASE; ++k)
        _u.insert(coord_t<T>(idx_t(_c.ia[k]-BASE,j),_c.a[k]));
  }


  static size_t ensure_structural_symmetry(
    const idx_t& _size,
    matrix_uncompressed_t & _u)
  {
    size_t nmodif = 0;
    for (size_t i=0; i<std::min(_size.i,_size.j); ++i)
      _u.insert(coord_t<T>(idx_t(i,i),T())).second? ++nmodif:nmodif;
#if 0
    for (bool modif=true; modif;) {
      modif = false;
      for (typename matrix_uncompressed_t::const_reverse_iterator c=_u.rbegin(); !modif && c!=_u.rend(); ++c)
        if (c->first.j != c->first.i) {
          coord_t<T> p(idx_t(c->first.j,c->first.i),T());
          (modif = _u.insert(p).second)? ++nmodif:nmodif;
        }
    }
#else
    for (typename matrix_uncompressed_t::const_iterator c=_u.begin(); c!=_u.end(); ++c)
      if (c->first.j != c->first.i)
        _u.insert( coord_t<T>(idx_t(c->first.j,c->first.i),T()) ).second? ++nmodif:nmodif;
#endif
    return nmodif;
  }


 private:
  // storage
  matrix_uncompressed_t matu;  // (uncompressed, in 0-based indexing)
  matrix_compressed_t   matc;  // (compressed, in BASE indexing)

};


}  // namespace lss
}  // namespace cf3


#endif

