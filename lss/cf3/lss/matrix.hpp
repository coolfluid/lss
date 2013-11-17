// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_matrix_hpp
#define cf3_lss_matrix_hpp


#include <iterator>

#include "common/Log.hpp"

#include "index.hpp"
#include "utilities.hpp"


namespace cf3 {
namespace lss {


/* -- matrix helper definitions --------------------------------------------- */

/// @brief Matrix print level
enum print_t { print_auto=0, print_size, print_signs, print_full, print_file };


/// @brief Matrix print level converter utility
print_t print_level(const int& i);


/// @brief Matrix structure compressed
template< typename T >
struct matrix_compressed_t {
  void clear() {
    nnu = nnz = 0;
    ia.clear();
    ja.clear();
    a.clear();
  }
  matrix_compressed_t& swap(matrix_compressed_t& _other) {
    std::swap(nnu,_other.nnu);
    std::swap(nnu,_other.nnz);
    ia.swap(_other.ia);
    ja.swap(_other.ja);
    a .swap(_other.a);
  }
  int nnu, nnz;               // number of rows/nonzeros
  std::vector< T > a;         // values
  std::vector< int > ia, ja;  // rows/column indices
};


/// @brief Dense matrix orientations
enum orientation_t { sort_by_column=0, sort_by_row=1 };


/// @brief Sparse/coordinate matrix entry
template< typename T >
struct coord_t : std::pair< idx_t, T > {
  coord_t(const idx_t& _idx, const T& _value) : std::pair< idx_t, T >(_idx,_value) {}
};


/// @brief Sparse/coordinate matrix sorting and compression, by row (functor)
template< typename T, typename _Key >
struct sort_by_row_t
{
  static const orientation_t orient = sort_by_row;

  virtual bool operator()(const _Key& a, const _Key& b) const {
    return (a.first.i<b.first.i? true  :
           (a.first.i>b.first.i? false :
           (a.first.j<b.first.j) ));
  }

  static void compress(
    const std::set< _Key, sort_by_row_t< T, _Key > >& _entries,
    matrix_compressed_t< T >& _comp)
  {
    typename std::set< _Key, sort_by_row_t< T, _Key > >::const_iterator it=_entries.begin();
    _comp.clear();
    _comp.ia.push_back(it->first.i);
    _comp.ja.reserve(_entries.size());
    _comp.a .reserve(_entries.size());
    for (size_t count,
         r =  _entries.begin() ->first.i;
         r <= _entries.rbegin()->first.i;
         ++r) {
      for (count=0; r==(it->first.i) && it!=_entries.end(); ++it, ++count) {
        _comp.ja.push_back(it->first.j);
        _comp.a .push_back(it->second);
      }
      _comp.ia.push_back(_comp.ia.back() + count);
    }
    _comp.nnu = _comp.ia.size() - 1;
    _comp.nnz = _entries.size();
  }

  static void uncompress(
    std::set< _Key, sort_by_row_t< T, _Key > >& _entries,
    const matrix_compressed_t< T >& _comp)
  {
    _entries.clear();
    const int base = _comp.ia.front();
    for (int i=0; i<_comp.ia.size()-1; ++i)
      for (int k=_comp.ia[i]-base; k<_comp.ia[i+1]-base; ++k)
        _entries.insert(coord_t<T>(idx_t(i+base,_comp.ja[k]),_comp.a[k]));
  }

};


/// @brief Sparse/coordinate matrix sorting and compression, by column (functor)
template< typename T, typename _Key >
struct sort_by_column_t
{
  static const orientation_t orient = sort_by_column;

  virtual bool operator()(const _Key& a, const _Key& b) const {
    return (a.first.j<b.first.j? true  :
           (a.first.j>b.first.j? false :
           (a.first.i<b.first.i)));
  }

  static void compress(
    const std::set< _Key, sort_by_column_t< T, _Key > >& _entries,
    matrix_compressed_t< T >& _comp)
  {
    typename std::set< _Key, sort_by_column_t< T, _Key > >::const_iterator it=_entries.begin();
    _comp.clear();
    _comp.ja.push_back(it->first.j);
    _comp.ia.reserve(_entries.size());
    _comp.a .reserve(_entries.size());
    for (size_t count,
         c =  _entries.begin() ->first.j;
         c <= _entries.rbegin()->first.j;
         ++c) {
      for (count=0; c==(it->first.j) && it!=_entries.end(); ++it, ++count) {
        _comp.ia.push_back(it->first.i);
        _comp.a .push_back(it->second);
      }
      _comp.ja.push_back(_comp.ja.back() + count);
    }
    _comp.nnu = _comp.ja.size() - 1;
    _comp.nnz = _entries.size();
  }

  static void uncompress(
      std::set< _Key, sort_by_column_t< T, _Key > >& _entries,
      const matrix_compressed_t< T >& _comp)
  {
    _entries.clear();
    const int base = _comp.ja.front();
    for (int j=0; j<_comp.ja.size()-1; ++j)
      for (int i=_comp.ja[j]-base; i<_comp.ja[j+1]-base; ++i)
        _entries.insert(coord_t<T>(idx_t(_comp.ia[i]+base,j+base),_comp.a[i]));
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
    m_size(),
    m_print(print_auto)
  {}
  virtual ~matrix() {}

  // -- functionality in remote implementation

  matrix& initialize(const size_t& i, const size_t& j, const double& _value=double()) { return IMPL::initialize(i,j,_value); }
  matrix& initialize(const std::vector< double >& _vector) { return IMPL::initialize(_vector); }
  matrix& initialize(const std::string& _fname)            { return IMPL::initialize(_fname); }
  matrix& initialize(const index_t& _index)                { return IMPL::initialize(_index); }

  matrix& clear()                         { return IMPL::initialize(m_size.i,m_size.j,double()); }
  matrix& operator=(const double& _value) { return IMPL::initialize(m_size.i,m_size.j,_value); }
  matrix& operator=(const matrix& _other) { return IMPL::operator=(_other); }
  matrix& zerorow(const size_t& i)        { return IMPL::zerorow(i); }

  // -- intrinsic functionality

  void swap(matrix& other) {
    std::swap(other.m_size,m_size);
    std::swap(other.m_zero,m_zero);
    std::swap(other.m_print,m_print);
  }

  virtual std::ostream& print(std::ostream& o) const {

    const T eps = 1.e3*std::numeric_limits< T >::epsilon();
    const print_t print_level(m_print? std::max(print_auto,std::min(m_print,print_file)) :
                             (m_size.i>100 || m_size.j>100? print_size  :
                             (m_size.i> 10 || m_size.j> 10? print_signs :
                                                            print_full )));

    if (print_level==print_size)  {
      o << "(" << size(0) << 'x' << size(1) << ") [ ... ]";
    }
    else if (print_level==print_signs) {
      o << "(" << size(0) << 'x' << size(1) << ") [";
      for (size_t i=0; i<size(0); ++i) {
        std::string str(size(1),'0');
        for (size_t j=0; j<size(1); ++j)
          str[j] = (operator()(i,j)>eps? '+' : (operator()(i,j)<-eps? '-':'0'));
        o << "\n  " << str;
      }
      o << " ]";
    }
    else if (print_level==print_full) {
      o << "(" << size(0) << 'x' << size(1) << ") [";
      for (size_t i=0; i<size(0); ++i) {
        std::ostringstream ss;
        for (size_t j=0; j<size(1); ++j)
          ss << operator()(i,j) << ", ";
        o << "\n  " << ss.str();
      }
      o << " ]";
    }
    else if (print_level==print_file) {
      o << "%%MatrixMarket matrix array real general\n"
        << size(0) << ' ' << size(1) << '\n';
      for (size_t j=0; j<size(1); ++j)
        for (size_t i=0; i<size(0); ++i)
          o << operator()(i,j) << '\n';
    }
    return o;
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
  print_t m_print;
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

  dense_matrix_vv& initialize(const index_t& _index) {
    clear();
    //FIXME implement
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

  dense_matrix_v& initialize(const index_t& _index) {
    clear();
    //FIXME implement
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
 * compression of rows/columns as necessary
 * T: storage type
 * BASE: column & row numbering base (0 or 1, other values won't work)
 * ORIENT: ordering of coordinate entries, as (i,j,v) tuples
 */
template< typename T, int BASE=0, typename ORIENT=sort_by_row_t< T, coord_t<T> > >
struct sparse_matrix :
  matrix< T,sparse_matrix< T,BASE,ORIENT > >
{
  typedef matrix< T,sparse_matrix< T,BASE,ORIENT > > matrix_base_t;
  typedef std::set< coord_t<T>, ORIENT > set_t;

  // cons/destructor
  sparse_matrix() : matrix_base_t() { clear(); }
  ~sparse_matrix() { clear(); }

  // initializations

  sparse_matrix& initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    if (idx_t(i,j)==matrix_base_t::m_size)
      return operator=(_value);
    CFwarn << "sparse_matrix: clearing matrix." << CFendl;
    clear();
    return *this;
  }

  sparse_matrix& initialize(const std::vector< double >& _vector) {
    idx_t& size = matrix_base_t::m_size;
    if (size.i*size.j!=_vector.size())
      throw std::runtime_error("sparse_matrix: assignment not consistent with current size.");
    CFwarn << "sparse_matrix: initialize from vector fully populates the otherwise sparse matrix." << CFendl;
    clear();
    for (size_t i=0, k=0; i<size.i; ++i)
      for (size_t j=0; j<size.j; ++j, ++k)
        entries.insert(coord_t<T>(idx_t(i+BASE,j+BASE),static_cast< T >(_vector[k])));
    return *this;
  }

  sparse_matrix& initialize(const std::string& _fname) {
    clear();
    idx_t& size = matrix_base_t::m_size;
    size.invalidate();
    try {

      using namespace std;
      ifstream f(_fname.c_str());
      if (!f) throw runtime_error("sparse_matrix: cannot open file.");
      const bool hasdot(string("."+_fname).find_last_of("."));


      // read format: MatrixMarket (*.mtx)
      if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {
        MatrixMarket::typecode_t t;
        if (!MatrixMarket::read_banner(f,t))  throw runtime_error("sparse_matrix: MatrixMarket invalid header (\"%%MatrixMarket ...\" not found).");
        if (!t.is_real() || !t.is_general())  throw runtime_error("sparse_matrix: MatrixMarket only \"(coordinate|array) real general\" supported.");
        if (!MatrixMarket::read_size(f,size.i,size.j,compressed.nnz))
          throw runtime_error("sparse_matrix: MatrixMarket invalid matrix/array size.");
      
        // read into set, line by line
        coord_t<T> p(idx_t(BASE,BASE),T());
        string line;        
        while (getline(f,line)) {
          if (line.find_first_of("%")==0) {}
          else if (t.is_dense()) {
            istringstream(line) >> p.second;
            entries.insert(p);
            if (++p.first.i-BASE >= size.i) {
              p.first.i = BASE;
              p.first.j++;
            }
          }
          else {
            istringstream(line) >> p.first.i >> p.first.j >> p.second;
            p.first.i += (BASE-1);
            p.first.j += (BASE-1);
            entries.insert(p);
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
        vector< T > &a = compressed.a;
        vector< int >
            &ia = compressed.ia,
            &ja = compressed.ja;
        int &nnu = compressed.nnu,
            &nnz = compressed.nnz;

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
    catch (const std::runtime_error& e) {
      throw std::runtime_error("sparse_matrix: cannot read file. " + std::string(e.what()));
    }


    // compress and return
    compress();
    if (size.i!=compressed.nnu || size.j<=*std::max_element(compressed.ja.begin(),compressed.ja.end())-BASE)
      throw std::runtime_error("sparse_matrix: after reading file, indexing not correct.");
    return *this;
  }

  sparse_matrix& initialize(const index_t& _index) {
    clear();
    //FIXME implement

    //    void initialize(const std::vector< std::vector< size_t > >& nz) {
    //      const int Nb=1;

    //      // set row indices
    //      nnu = Nb * (int) nz.size();
    //      ia.assign(1,BASE);
    //      ia.reserve(nnu+1);
    //      for (int R=0; R<(int) nz.size(); ++R)
    //        for (int i=0; i<Nb; ++i)
    //          ia.push_back(ia.back() + Nb * (int) nz[R].size());
    //      nnz = ia.back() - ia.front();

    //      // set column indices
    //      ja.assign(nnz,0);
    //      for (size_t R=0; R<(size_t) nz.size(); ++R) {
    //        for (int r=0; r<Nb; ++r) {
    //          int k = ia[R*Nb+r]-BASE;
    //          for (size_t I=0; I<(size_t) nz[R].size(); ++I)
    //            for (int c=0; c<Nb; ++c)
    //              ja[k++] = (int) (Nb*nz[R][I]) + c + BASE;
    //        }
    //      }

    //      // set entries
    //      a.assign(nnz,T());
    //    }

    return *this;
  }

  sparse_matrix& clear() {
    // NOTE: clearing a matrix puts it into uncompressed state (nnu==0)
    matrix_base_t::m_size.clear();
    entries.clear();
    compressed.clear();
    return *this;
  }

  sparse_matrix& operator=(const double& _value) {
    CFwarn << "sparse_matrix: assigning a value to a sparse matrix only affects the populated entries." << CFendl;
    if (is_compressed())
      compressed.a.assign(compressed.a.size(),static_cast< const T >(_value));
    else
      for (typename set_t::iterator i = entries.begin(); i!=entries.end(); ++i)
        const_cast< T& >(i->second) = _value;
    return *this;
  }

  sparse_matrix& operator=(const sparse_matrix& _other) {
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
    entries = _other.entries;
    compressed = _other.compressed;
    return *this;
  }

  sparse_matrix& zerorow(const size_t& _i) {
    //FIXME: this is specific to CSR structure! this container should also work for CSC!
    std::fill_n(compressed.a.begin()+compressed.ia[_i]-BASE,compressed.ia[_i+1]-compressed.ia[_i],T());
    return *this;
  }

  sparse_matrix& swap(sparse_matrix& _other) {
    matrix_base_t::swap(_other);
    std::swap(entries,_other.entries);
    std::swap(compressed,_other.compressed);
    return *this;
  }

  std::ostream& print(std::ostream& o) const {
    using namespace std;
    const idx_t&   size  = matrix_base_t::m_size;
    const print_t& print = matrix_base_t::m_print;

    const T eps = 1.e3*numeric_limits< T >::epsilon();
    const print_t print_level(print? max(print_auto,min(print,print_file)) :
                             (size.i>100 || size.j>100? print_size  :
                             (size.i> 10 || size.j> 10? print_signs :
                                                        print_full )));

    if (print_level==print_size)  {
      o << "(" << size.i << 'x' << size.j << ">=" << compressed.a.size() << ") [ ... ]";
    }
    else if (print_level==print_signs || print_level==print_full) {
      o << "(" << size.i << 'x' << size.j << ">=" << compressed.a.size() << ") [ ";

      matrix_compressed_t< T > tmp_comp;
      if (!is_compressed())
        ORIENT::compress(entries,tmp_comp);
      const matrix_compressed_t< T >& _comp(is_compressed()? compressed : tmp_comp);

      if (print_level==print_signs) {
        string str;
        for (size_t i=0; i<size.i; ++i) {
          str.assign(size.j,' ');
          for (int k=_comp.ia[i]-BASE; k<_comp.ia[i+1]-BASE; ++k)
            str[ _comp.ja[k]-BASE ] = (_comp.a[k]>eps? '+' : (_comp.a[k]<-eps? '-' : '.' ));
          o << "\n  " << str;
        }
      }
      else {  // print_level==print_full
        for (size_t i=0; i<size.i; ++i) {
          vector< T > row(size.j,T());
          for (int k=_comp.ia[i]-BASE; k<_comp.ia[i+1]-BASE; ++k)
            row[ _comp.ja[k]-BASE ] = _comp.a[k];
          o << "\n  ";
          copy(row.begin(),row.end(),ostream_iterator< T >(o,", "));
        }
      }

      o << " ]";
    }
    else if (print_level==print_file) {

      set_t tmp_entries;
      if (is_compressed())
        ORIENT::uncompress(tmp_entries,compressed);
      const set_t& _entries(is_compressed()? tmp_entries:entries);

      // write format: MatrixMarket (*.mtx)
      o << "%%MatrixMarket matrix coordinate real general\n"
        << size.i << ' ' << size.j << ' ' << _entries.size() << '\n';
      for (typename set_t::const_iterator i = _entries.begin(); i!=_entries.end(); ++i)
        o << (i->first.i-BASE+1) << ' '<< (i->first.j-BASE+1) << ' '<< i->second << '\n';

    }
    return o;
  }

  // indexing
  const T& operator()(const size_t& _i, const size_t& _j) const {
    if (is_compressed()) {
      const int i=getindex(_i,_j);
      return (i<0? matrix_base_t::m_zero:compressed.a[i]);
    }
    const coord_t<T> p(idx_t(_i+BASE,_j+BASE),T());
    typename set_t::const_iterator i = entries.find(p);
    return i!=entries.end()? i->second : matrix_base_t::m_zero;
  }

  T& operator()(const size_t& _i, const size_t& _j) {
    if (is_compressed()) {
      const int i=getindex(_i,_j);
      return i>=BASE? compressed.a[i] : matrix_base_t::m_zero;
    }
    const coord_t<T> p(idx_t(_i+BASE,_j+BASE),T());
    typename set_t::iterator i = entries.insert(p).first;
    /*
     * note: removing constness is safe because changing the coordinate value
     * does not affect the comparison operator in the matrix ordering
     */
    return const_cast< T& >(i->second);
  }

  int getindex(const size_t& i, const size_t& j) const {
    if (ORIENT::orient==sort_by_row)
      for (int k=compressed.ia[i]-BASE; k<compressed.ia[i+1]-BASE; ++k)
        if (compressed.ja[k]-BASE==(int) j)
          return k;
    else if (ORIENT::orient==sort_by_column)
      for (int k=compressed.ja[j]-BASE; k<compressed.ja[j+1]-BASE; ++k)
        if (compressed.ia[k]-BASE==(int) i)
          return k;
    return -1;
  }


  // compression/uncompression

  inline bool is_compressed() const { return compressed.nnu; }

  matrix_compressed_t< T >& compress() {
    if (!is_compressed()) {
      CFdebug << "sparse_matrix: compress..." << CFendl;
      CFdebug << "sparse_matrix: additional entries to preserving symmetry: " << ensure_structural_symmetry() << '.' << CFendl;
      ORIENT::compress(entries,compressed);
      CFdebug << "sparse_matrix: compress." << CFendl;
    }
    entries.clear();
    return compressed;
  }

  set_t& uncompress() {
    if (is_compressed()) {
      CFdebug << "sparse_matrix: uncompress..." << CFendl;
      ORIENT::uncompress(entries,compressed);
      CFdebug << "sparse_matrix: additional entries to preserving symmetry: " << ensure_structural_symmetry() << '.' << CFendl;
      CFdebug << "sparse_matrix: uncompress." << CFendl;
    }
    compressed.clear();
    return entries;
  }

  size_t ensure_structural_symmetry() {
    if (is_compressed())
      uncompress();
    size_t nmodif = 0;
    for (size_t i=0; i<this->size(0); ++i)
      entries.insert(coord_t<T>(idx_t(i+BASE,i+BASE),T())).second? ++nmodif:nmodif;
    for (bool modif=true; modif;) {
      modif = false;
      for (typename set_t::const_reverse_iterator c=entries.rbegin(); !modif && c!=entries.rend(); ++c) {
        coord_t<T> p(idx_t(c->first.j,c->first.i),T());
        (modif = entries.insert(p).second)? ++nmodif:nmodif;
      }
    }
    return nmodif;
  }

 // storage
 private:
  set_t entries;                        // (uncompressed)
  matrix_compressed_t< T > compressed;  // (compressed)

};


}  // namespace lss
}  // namespace cf3


#endif

