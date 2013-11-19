// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_matrix_hpp
#define cf3_lss_matrix_hpp


#include <iterator>

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
    m_zero(std::numeric_limits< T >::quiet_NaN()),
    m_size()
  {}
  virtual ~matrix() {}

  // -- functionality in remote implementation

  matrix& initialize(const size_t& i, const size_t& j, const double& _value=double()) { return IMPL::initialize(i,j,_value); }
  matrix& initialize(const std::vector< double >& _vector) { return IMPL::initialize(_vector); }
  matrix& initialize(const std::string& _fname)            { return IMPL::initialize(_fname); }

  matrix& clear()                         { return IMPL::initialize(m_size.i,m_size.j,double()); }
  matrix& operator=(const double& _value) { return IMPL::initialize(m_size.i,m_size.j,_value); }
  matrix& operator=(const matrix& _other) { return IMPL::operator=(_other); }
  matrix& zerorow(const size_t& i)        { return IMPL::zerorow(i); }

  // -- intrinsic functionality

  void swap(matrix& other) {
    std::swap(other.m_size,m_size);
    std::swap(other.m_zero,m_zero);
  }

  virtual void print(std::ostream& o, const print_t& l=print_auto) const {

    const T eps = 1.e3*std::numeric_limits< T >::epsilon();
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
          for (size_t j=0; j<size(1); ++j)
            row_s[j] = (operator()(i,j)> eps? '+'
                     : (operator()(i,j)<-eps? '-'
                     :                        '0' ));
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
        o << "%%MatrixMarket matrix array real general\n"
          << size(0) << ' ' << size(1) << '\n';
        for (size_t j=0; j<size(1); ++j)
          for (size_t i=0; i<size(0); ++i)
            o << operator()(i,j) << '\n';
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

  dense_matrix_vv& initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    if (idx_t(i,j).is_valid_size()) {
      clear();
      matrix_base_t::m_size = idx_t(i,j);
      if (size(0)*size(1))
        a.assign(ORIENT? size(0):size(1),std::vector< T >(
                 ORIENT? size(1):size(0),static_cast< T >(_value) ));
    }
    return *this;
  }

  dense_matrix_vv& initialize(const std::vector< double >& _vector) {
    if (size(0)*size(1)!=_vector.size())
      throw std::runtime_error("dense_matrix_vv: assignment not consistent with current size.");
    for (size_t i=0, k=0; i<size(0); ++i)
      for (size_t j=0; j<size(1); ++j, ++k)
        operator()(i,j) = (type_is_equal< T, double >()? (const T) _vector[k]
                                                       : static_cast< const T >(_vector[k]) );
    return *this;
  }

  dense_matrix_vv& initialize(const std::string& _fname) {
    using namespace std;
    clear();
    idx_t& size = matrix_base_t::m_size;
    size.invalidate();
    try {

      ifstream f(_fname.c_str());
      if (!f) throw runtime_error("sparse_matrix: cannot open file.");
      const bool hasdot(string("."+_fname).find_last_of("."));


      // read format: MatrixMarket (*.mtx)
      if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {


        // read matrix size and properties
        MatrixMarket::typecode_t t;
        int nnz(0);
        if (!MatrixMarket::read_banner(f,t))                throw runtime_error("dense_matrix_vv: MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
        if (!MatrixMarket::read_size(f,size.i,size.j,nnz))  throw runtime_error("dense_matrix_vv: MatrixMarket: invalid matrix/array size.");
        if (!t.is_real() || !t.is_general())                throw runtime_error("dense_matrix_vv: MatrixMarket: only \"(coordinate|array) real general\" supported.");

        // read into row/column-oriented dense matrix, line by line
        a.assign(ORIENT? size.i:size.j,vector< T >(
                 ORIENT? size.j:size.i,T() ));
        idx_t ij(1,1);
        string line;
        T v;
        while (getline(f,line)) {
          if (line.find_first_of("%")==0) {}
          else if (t.is_dense()) {
            istringstream(line) >> operator()(ij.i-1,ij.j-1);
            if (++ij.i > size.i) {
              ij.i = 1;
              ij.j++;
            }
          }
          else {
            istringstream(line) >> ij.i >> ij.j >> v;
            operator()(ij.i-1,ij.j-1) = v;
          }
        }


      }
      // read format: CSR (*.csr)
      else if (hasdot && _fname.substr(_fname.find_last_of("."))==".csr") {


        // read matrix size
        string line;
        while (!size.is_valid_size()) {
          if (getline(f,line) && line.find_first_of("%")!=0)
            istringstream(line) >> size.i >> size.j;
        }

        // read rows and column indices, converting to base 0 (easier)
        vector< int > ia;
        ia.reserve(size.i+1);
        for (int i; ia.size()<size.i+1 && f >> i;)
          ia.push_back(i);

        vector< int > ja;
        ja.reserve(ia.back()-ia.front());
        for (int i; ja.size()<ia.back()-ia.front() && f >> i;)
          ja.push_back(i);

        for_each(ja.begin(),ja.end(),base_conversion_t(-ia.front()));
        for_each(ia.begin(),ia.end(),base_conversion_t(-ia.front()));

        // populate per row
        a.assign(ORIENT? size.i:size.j,vector< T >(
                 ORIENT? size.j:size.i,T() ));
        for (size_t i=0; i<size.i; ++i)
          for (int k=ia[i]; k<ia[i+1] && f; ++k)
            f >> operator()(i,ja[k]);


      }


      // read format: file format not detected
      else {
        throw runtime_error("dense_matrix_vv: file format not detected.");
      }
      f.close();


    }
    catch (const runtime_error& e) {
      throw runtime_error("dense_matrix_vv: cannot read file. " + string(e.what()));
    }
    return *this;
  }

  dense_matrix_vv& clear() {
    a.clear();
    matrix_base_t::m_size.clear();
    return *this;
  }

  dense_matrix_vv& operator=(const dense_matrix_vv& _other) {
    a = _other.a;
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
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

  dense_matrix_v& initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    if (idx_t(i,j).is_valid_size()) {
      clear();
      matrix_base_t::m_size = idx_t(i,j);
      if (size(0)*size(1))
        a.assign(size(0)*size(1),static_cast< T >(_value));
    }
    return *this;
  }

  dense_matrix_v& initialize(const std::vector< double >& _vector) {
    if (a.size()!=_vector.size())
      throw std::runtime_error("dense_matrix_v: assignment not consistent with current size.");

    std::vector< T > w;
    if (!type_is_equal< T, double >()) {
      w.resize(_vector.size());
      std::transform(_vector.begin(),_vector.end(),w.begin(),type_conversion_t< double, T >() );
    }
    const std::vector< T >& v(type_is_equal< T, double >()? (const std::vector< T >&) _vector : w);

    if (ORIENT) { a = v; }
    else {
      for (size_t i=0, k=0; i<size(0); ++i)
        for (size_t j=0; j<size(1); ++j, ++k)
          operator()(i,j) = v[k];
    }
    return *this;
  }

  dense_matrix_v& initialize(const std::string& _fname) {
    using namespace std;
    clear();
    idx_t& size = matrix_base_t::m_size;
    size.invalidate();
    try {

      ifstream f(_fname.c_str());
      if (!f) throw runtime_error("sparse_matrix: cannot open file.");
      const bool hasdot(string("."+_fname).find_last_of("."));


      // read format: MatrixMarket (*.mtx)
      if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {


        // read matrix size and properties
        MatrixMarket::typecode_t t;
        int nnz(0);
        if (!MatrixMarket::read_banner(f,t))                throw runtime_error("dense_matrix_v: MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
        if (!MatrixMarket::read_size(f,size.i,size.j,nnz))  throw runtime_error("dense_matrix_v: MatrixMarket: invalid matrix/array size.");
        if (!t.is_real() || !t.is_general())                throw runtime_error("dense_matrix_v: MatrixMarket: only \"(coordinate|array) real general\" supported.");

        // read into row/column-oriented dense matrix, line by line
        a.assign(size.i*size.j,T());
        idx_t ij(1,1);
        string line;
        T v;
        while (getline(f,line)) {
          if (line.find_first_of("%")==0) {}
          else if (t.is_dense()) {
            istringstream(line) >> operator()(ij.i-1,ij.j-1);
            if (++ij.i > size.i) {
              ij.i = 1;
              ij.j++;
            }
          }
          else {
            istringstream(line) >> ij.i >> ij.j >> v;
            operator()(ij.i-1,ij.j-1) = v;
          }
        }


      }
      // read format: CSR (*.csr)
      else if (hasdot && _fname.substr(_fname.find_last_of("."))==".csr") {


        // read matrix size
        string line;
        while (!size.is_valid_size()) {
          if (getline(f,line) && line.find_first_of("%")!=0)
            istringstream(line) >> size.i >> size.j;
        }

        // read rows and column indices, converting to base 0 (easier)
        vector< int > ia;
        ia.reserve(size.i+1);
        for (int i; ia.size()<size.i+1 && f >> i;)
          ia.push_back(i);

        vector< int > ja;
        ja.reserve(ia.back()-ia.front());
        for (int i; ja.size()<ia.back()-ia.front() && f >> i;)
          ja.push_back(i);

        for_each(ja.begin(),ja.end(),base_conversion_t(-ia.front()));
        for_each(ia.begin(),ia.end(),base_conversion_t(-ia.front()));

        // populate per row
        a.assign(size.i*size.j,T());
        for (size_t i=0; i<size.i; ++i)
          for (int k=ia[i]; k<ia[i+1] && f; ++k)
            f >> operator()(i,ja[k]);


      }


      // read format: file format not detected
      else {
        throw runtime_error("dense_matrix_v: file format not detected.");
      }
      f.close();


    }
    catch (const runtime_error& e) {
      throw runtime_error("dense_matrix_v: cannot read file. " + string(e.what()));
    }
    return *this;
  }

  // assignments

  dense_matrix_v& operator=(const dense_matrix_v& _other) {
    a = _other.a;
    matrix_base_t::m_size = _other.dense_matrix_v::matrix_base_t::m_size;
    return *this;
  }

  dense_matrix_v& operator=(const double& _value) {
    return initialize(size(0),size(1),_value);
  }

  void swap(dense_matrix_v& other) {
    other.a.swap(a);
    matrix_base_t::swap(other);
  }

  // clearing

  dense_matrix_v& clear() {
    a.clear();
    matrix_base_t::m_size.clear();
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

  // cons/destructor
  sparse_matrix() : matrix_base_t() { clear(); }
  ~sparse_matrix() { clear(); }

  // initializations

  sparse_matrix& initialize(const size_t& i, const size_t& j, const double&) {
    if (idx_t(i,j).is_valid_size()) {
      // initialization value ignored (after all, this is a sparse matrix)
      clear();
      matrix_base_t::m_size = idx_t(i,j);
      ensure_structural_symmetry(matrix_base_t::m_size,matu);
    }
    else {
      CFwarn << "sparse_matrix: invalid size: (" << i << ',' << j << ')' << CFendl;
    }
    return *this;
  }

  sparse_matrix& initialize(const std::vector< double >& _vector) {
    idx_t& size = matrix_base_t::m_size;
    if (size.i*size.j!=_vector.size())
      throw std::runtime_error("sparse_matrix: assignment not consistent with current size.");
    CFwarn << "sparse_matrix: initialize from vector fully populates the (otherwise sparse) matrix." << CFendl;
    matu.clear();
    matc.clear();
    for (size_t i=0, k=0; i<size.i; ++i)
      for (size_t j=0; j<size.j; ++j, ++k)
        matu.insert(coord_t<T>(idx_t(i+BASE,j+BASE),static_cast< T >(_vector[k])));
    compress();  // since matrix is fully populated, this won't hurt.
    return *this;
  }

  sparse_matrix& initialize(const std::string& _fname) {
    using namespace std;
    clear();
    idx_t& size = matrix_base_t::m_size;
    size.invalidate();
    try {
      ifstream f(_fname.c_str());
      if (!f)
        throw runtime_error("sparse_matrix: cannot open file.");
      const bool hasdot(string("."+_fname).find_last_of("."));


      // read format: MatrixMarket (*.mtx)
      if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {
        MatrixMarket::typecode_t t;
        if (!MatrixMarket::read_banner(f,t))  throw runtime_error("sparse_matrix: MatrixMarket invalid header (\"%%MatrixMarket ...\" not found).");
        if (!t.is_real() || !t.is_general())  throw runtime_error("sparse_matrix: MatrixMarket only \"(coordinate|array) real general\" supported.");
        if (!MatrixMarket::read_size(f,size.i,size.j,matc.nnz))
          throw runtime_error("sparse_matrix: MatrixMarket invalid matrix/array size.");
      
        // read into set, line by line
        coord_t<T> p(idx_t(BASE,BASE),T());
        string line;        
        while (getline(f,line)) {
          if (line.find_first_of("%")==0) {}
          else if (t.is_dense()) {
            istringstream(line) >> p.second;
            matu.insert(p);
            if (++p.first.i-BASE >= size.i) {
              p.first.i = BASE;
              p.first.j++;
            }
          }
          else {
            istringstream(line) >> p.first.i >> p.first.j >> p.second;
            p.first.i += (BASE-1);
            p.first.j += (BASE-1);
            matu.insert(p);
          }
        }
      }


      // read format: CSR (*.csr)
      else if (hasdot && _fname.substr(_fname.find_last_of("."))==".csr") {


        // read matrix size
        string line;
        while (!size.is_valid_size()) {
          if (getline(f,line) && line.find_first_of("%")!=0)
            istringstream(line) >> size.i >> size.j;
        }

        /*
         * read into temporary compressed structure and convert to intended base
         * - ia has indices in increasing order
         * - ja first entry is read while building ia
         */
        vector< T > &a = matc.a;
        vector< int >
            &ia = matc.ia,
            &ja = matc.ja;
        int &nnu = matc.nnu,
            &nnz = matc.nnz;

        ia.reserve(size.i+1);
        ja.assign(1,0);
        for (int &i=ja[0], j=-1; f>>i && j<i;)
          ia.push_back(j=i);
        nnu = static_cast< int >(size.i);
        nnz = static_cast< int >(ia.back()-ia.front());

        ja.reserve(nnz);
        a .reserve(nnz);
        for (int i; ja.size()<nnz && f>>i;)  ja.push_back(i);
        for (T   i; a .size()<nnz && f>>i;)  a .push_back(i);

        // convrt to intended base, and uncompress
        for_each(ja.begin(),ja.end(),base_conversion_t(BASE-ia.front()));
        for_each(ia.begin(),ia.end(),base_conversion_t(BASE-ia.front()));
        uncompress();

      }


      // read format: file format not detected
      else {
        throw runtime_error("sparse_matrix: file format not detected.");
      }
      f.close();


    }
    catch (const runtime_error& e) {
      throw runtime_error("sparse_matrix: cannot read file. " + string(e.what()));
    }


    // compress and return
    compress();
    if (size.i!=matc.nnu || size.j<=*max_element(matc.ja.begin(),matc.ja.end())-BASE)
      throw runtime_error("sparse_matrix: after reading file, indexing not correct.");
    return *this;
  }

  sparse_matrix& clear() {
    matrix_base_t::m_size.clear();
    matu.clear();
    matc.clear();
    return *this;
  }

  sparse_matrix& operator=(const double& _value) {
    CFwarn << "sparse_matrix: assigning a value to a sparse matrix only affects the populated entries." << CFendl;
    if (is_compressed())
      matc.a.assign(matc.a.size(),static_cast< const T >(_value));
    else
      for (typename matrix_uncompressed_t::iterator i = matu.begin(); i!=matu.end(); ++i)
        const_cast< T& >(i->second) = _value;
    return *this;
  }

  sparse_matrix& operator=(const sparse_matrix& _other) {
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
    matu = _other.matu;
    matc = _other.matc;
    return *this;
  }

  sparse_matrix& zerorow(const size_t& _i) {
    if ( _i<matrix_base_t::m_size.i) {
      if (is_compressed() && ORIENT) {
        // compressed, row-oriented matrix
        std::fill(&matc.ia[_i],&matc.ia[_i+1],T());
      }
      else if (is_compressed() && !ORIENT) {
        // compressed, column-oriented matrix
        for (int j=0; j<matc.nnu; ++j)
          for (int k=matc.ja[j]-BASE; k<matc.ja[j+1]-BASE; ++k)
            if (matc.ia[k]-BASE==(int) _i)
              matc.a[k] = T();
      }
      else {
        // uncompressed, coordinate matrix
        for (typename matrix_uncompressed_t::iterator it = matu.begin(); it!=matu.end(); ++it)
          if (it->first.i-BASE==(size_t) _i)
            const_cast< T& >(it->second) = T();
      }
    }
    return *this;
  }

  sparse_matrix& swap(sparse_matrix& _other) {
    matrix_base_t::swap(_other);
    std::swap(matu,_other.matu);
    std::swap(matc,_other.matc);
    return *this;
  }

  void print(std::ostream& o, const print_t& l=print_auto) const {
    using namespace std;
    const idx_t&  size = matrix_base_t::m_size;
    const T       eps  = 1.e3*numeric_limits< T >::epsilon();
    const print_t lvl(l? max(print_size,min(l,print_file)) :
                     (size.i>100 || size.j>100? print_size  :
                     (size.i> 10 || size.j> 10? print_signs :
                                                print_full )));

    if (lvl==print_size)  {
      o << "(" << size.i << 'x' << size.j << ">=" << (matc.nnu+matu.size()) << ") [ ... ]";
    }
    else {

      // (requires an uncompressed structure)
      matrix_uncompressed_t tmp;
      if (is_compressed())
        uncompress(size,tmp,matc);
      const matrix_uncompressed_t& mat(is_compressed()? tmp:matu);

      string row_s;
      vector< T > row_v;
      switch (lvl) {

        case print_signs:
          o << "(" << size.i << 'x' << size.j << ">=" << mat.size() << ") [ ";
          for (size_t i=0; i<size.i; ++i) {
            row_s.assign(size.j,' ');
            for (typename matrix_uncompressed_t::const_iterator it = mat.begin(); it!=mat.end(); ++it)
              if (it->first.i==i+BASE)
                row_s[ it->first.j-BASE ] = (it->second> eps? '+'
                                          : (it->second<-eps? '-'
                                          :                   '.' ));
            o << "\n  " << row_s;
          }
          o << " ]";
          break;

        case print_full:
          o << "(" << size.i << 'x' << size.j << ">=" << mat.size() << ") [ ";
          for (size_t i=0; i<size.i; ++i) {
            row_v.assign(size.j,T());
            for (typename matrix_uncompressed_t::const_iterator it = mat.begin(); it!=mat.end(); ++it)
              if (it->first.i==i+BASE)
                row_v[ it->first.j-BASE ] = it->second;
            o << "\n  ";
            copy(row_v.begin(),row_v.end(),ostream_iterator< T >(o,", "));
          }
          o << " ]";
          break;

        case print_file:
          o << "%%MatrixMarket matrix coordinate real general\n"
            << size.i << ' ' << size.j << ' ' << mat.size() << '\n';
          for (typename matrix_uncompressed_t::const_iterator i = mat.begin(); i!=mat.end(); ++i)
            o << (i->first.i-BASE+1) << ' '<< (i->first.j-BASE+1) << ' '<< i->second << '\n';
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
      typename matrix_uncompressed_t::const_iterator it = matu.find(coord_t<T>(idx_t(i+BASE,j+BASE),T()));
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
    // find/insert new entry, and a structurally symmetric pair. the constness
    // removal is safe because the entry value does not change matrix ordering
    uncompress();
    if (i!=j)
      matu.insert( coord_t<T>(idx_t(j+BASE,i+BASE),T()) );
    typename matrix_uncompressed_t::iterator it =
      matu.insert( coord_t<T>(idx_t(i+BASE,j+BASE),T()) ).first;
    return const_cast< T& >(it->second);
  }


  // compression/uncompression

  matrix_compressed_t& compress() {
    if (!is_compressed()) {
      compress(matrix_base_t::m_size,matu,matc);
      matu.clear();
    }
    return matc;
  }

  matrix_uncompressed_t& uncompress() {
    if (is_compressed()) {
      uncompress(matrix_base_t::m_size,matu,matc);
      matc.clear();
    }
    return matu;
  }


 private:
  // compression utilities (not to use outside this context)

  inline bool is_compressed() const { return matc.nnu; }

  static void compress(
    const idx_t& _size,
    const matrix_uncompressed_t & _u,
    matrix_compressed_t& _c)
  {
    _c.clear();

    typename matrix_uncompressed_t::const_iterator it=_u.begin();
    if (it==_u.end())
      return;

    // row-oriented matrix
    if (ORIENT) {
      _c.ia.push_back(it->first.i);
      _c.ja.reserve(_u.size());
      _c.a .reserve(_u.size());
      for (size_t count,
           r =  _u.begin() ->first.i;
           r <= _u.rbegin()->first.i;
           ++r) {
        for (count=0; r==(it->first.i) && it!=_u.end(); ++it, ++count) {
          _c.ja.push_back(it->first.j);
          _c.a .push_back(it->second);
        }
        _c.ia.push_back(_c.ia.back() + count);
      }
    }

    // column-oriented matrix
    else {
      _c.ja.push_back(it->first.j);
      _c.ia.reserve(_u.size());
      _c.a .reserve(_u.size());
      for (size_t count,
           c =  _u.begin() ->first.j;
           c <= _u.rbegin()->first.j;
           ++c) {
        for (count=0; c==(it->first.j) && it!=_u.end(); ++it, ++count) {
          _c.ia.push_back(it->first.i);
          _c.a .push_back(it->second);
        }
        _c.ja.push_back(_c.ja.back() + count);
      }
    }

    _c.nnu = (ORIENT? _size.i:_size.j);
    _c.nnz = _u.size();
  }


  static void uncompress(
      const idx_t& _size,
      matrix_uncompressed_t & _u,
      const matrix_compressed_t& _c)
  {
    _u.clear();
    for (int i=0; ORIENT && i<_c.ia.size()-1; ++i)
      for (int k=_c.ia[i]-BASE; k<_c.ia[i+1]-BASE; ++k)
        _u.insert(coord_t<T>(idx_t(i+BASE,_c.ja[k]),_c.a[k]));
    for (int j=0; !ORIENT && j<_c.ja.size()-1; ++j)
      for (int k=_c.ja[j]-BASE; k<_c.ja[j+1]-BASE; ++k)
        _u.insert(coord_t<T>(idx_t(_c.ia[k]+BASE,j+BASE),_c.a[k]));
    ensure_structural_symmetry(_size,_u);
  }


  static size_t ensure_structural_symmetry(
    const idx_t& _size,
    matrix_uncompressed_t & _u)
  {
    size_t nmodif = 0;
    for (size_t i=0; i<std::min(_size.i,_size.j); ++i)
      _u.insert(coord_t<T>(idx_t(i+BASE,i+BASE),T())).second? ++nmodif:nmodif;
    for (bool modif=true; modif;) {
      modif = false;
      for (typename matrix_uncompressed_t::const_reverse_iterator c=_u.rbegin(); !modif && c!=_u.rend(); ++c)
        if (c->first.j != c->first.i) {
          coord_t<T> p(idx_t(c->first.j,c->first.i),T());
          (modif = _u.insert(p).second)? ++nmodif:nmodif;
        }
    }
    return nmodif;
  }


 private:
  // storage
  matrix_uncompressed_t matu;  // (uncompressed)
  matrix_compressed_t   matc;  // (compressed)

};


}  // namespace lss
}  // namespace cf3


#endif

