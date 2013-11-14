// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_detail_matrix_hpp
#define cf3_lss_detail_matrix_hpp


#include <iterator>

#include "common/Log.hpp"

#include "index.hpp"
#include "utilities.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- matrix interface and implementations ---------------------------------- */


// utilities
enum orientation_t { sort_by_column=0, sort_by_row=1 };
enum print_t { print_auto, print_size, print_signs, print_full };
print_t print_level(const int& i);


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
    const print_t print_level(m_print? std::min(m_print,print_full) :
                             (m_size.i>100 || m_size.j>100? print_size  :
                             (m_size.i> 10 || m_size.j> 10? print_signs :
                                                            print_full )));
    o << "(" << size(0) << 'x' << size(1) << ") [ ";
    if      (print_level==print_size)  { o << "..."; }
    else if (print_level==print_signs) {
      for (size_t i=0; i<size(0); ++i) {
        std::string str(size(1),'0');
        for (size_t j=0; j<size(1); ++j)
          str[j] = (operator()(i,j)>eps? '+' : (operator()(i,j)<-eps? '-':'0'));
        o << "\n  " << str;
      }
    }
    else if (print_level==print_full) {
      for (size_t i=0; i<size(0); ++i) {
        std::ostringstream ss;
        for (size_t j=0; j<size(1); ++j)
          ss << operator()(i,j) << ", ";
        o << "\n  " << ss.str();
      }
    }
    return o << " ]";
  }

  virtual void print(std::string& _fname) const {
    using namespace std;
    ofstream f(_fname.c_str());
    if (!f)
      throw runtime_error("matrix: cannot open file.");

    const bool hasdot(string("."+_fname).find_last_of("."));
    if (hasdot && _fname.substr(_fname.find_last_of("."))==".mtx") {
      f << "%%MatrixMarket matrix array real general\n"
        << size(0) << ' ' << size(1) << '\n';
      for (size_t j=0; j<size(1); ++j)
        for (size_t i=0; i<size(0); ++i)
          f << operator()(i,j) << '\n';
    }

    else {
      ostringstream msg;
      msg << "matrix: cannot write file \"" << _fname << "\" (only \"*.mtx\" is supported).";
      throw runtime_error(msg.str());
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
    clear();
    try { read_dense< T >(_fname,ORIENT,matrix_base_t::m_size,a); }
    catch (const std::runtime_error& e) {
      throw std::runtime_error("dense_matrix_vv: " + std::string(e.what()) + " cannot read file.");
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
    clear();
    std::vector< std::vector< T > > another_a;
    try { read_dense< T >(_fname,ORIENT,matrix_base_t::m_size,another_a); }
    catch (const std::runtime_error& e) {
      throw std::runtime_error("dense_matrix_v: " + std::string(e.what()) + " cannot read file.");
    }
    a.resize(size(0)*size(1));
    for (size_t i=0, k=0; i<(ORIENT? size(0):size(1)); ++i, k+=(ORIENT? size(1):size(0)))
      std::copy(another_a[i].begin(),another_a[i].end(),a.begin()+k);
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
 * SORT: ordering of coordinate entries, as (i,j,v) tuples
 */
template< typename T, int BASE=0, typename SORT=sort_by_row_t< coord_t<T> > >
struct sparse_matrix :
  matrix< T,sparse_matrix< T,BASE,SORT > >
{
  typedef matrix< T,sparse_matrix< T,BASE,SORT > > matrix_base_t;
  typedef std::set< coord_t<T>, SORT > set_t;

  // cons/destructor
  sparse_matrix() : matrix_base_t() { clear(); }
  ~sparse_matrix() { clear(); }

  // initializations

  sparse_matrix& initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    if (idx_t(i,j)==matrix_base_t::m_size)
      return operator=(_value);
    throw std::runtime_error("sparse_matrix: resizing not available.");
    return *this;
  }

  sparse_matrix& initialize(const std::vector< double >& _vector) {
    throw std::runtime_error("sparse_matrix_csr: initialize from vector is not possible.");
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
        if (!MatrixMarket::read_banner(f,t))                throw runtime_error("sparse_matrix: MatrixMarket invalid header (\"%%MatrixMarket ...\" not found).");
        if (!MatrixMarket::read_size(f,size.i,size.j,nnz))  throw runtime_error("sparse_matrix: MatrixMarket invalid matrix/array size.");
        if (!t.is_real() || !t.is_general())                throw runtime_error("sparse_matrix: MatrixMarket only \"(coordinate|array) real general\" supported.");
      
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
         * temporarily read into ia and ja arrays and convert to intended base
         * - ia has indices in increasing order
         * - ja first entry is read while building ia
         */
        ia.reserve(size.i+1);
        ja.assign(1,0);
        for (int &i=ja[0], j=-1; f>>i && j<i;)
          ia.push_back(j=i);
        ja.reserve(ia.back()-ia.front());
        for (int i; ja.size()<ja.capacity() && f>>i;)
          ja.push_back(i);

        for_each(ja.begin(),ja.end(),base_conversion_t(BASE-ia.front()));
        for_each(ia.begin(),ia.end(),base_conversion_t(BASE-ia.front()));

        // convert into set
        for (int i=0; i<ia.size()-1; ++i) {
          coord_t<T> p(idx_t(i+BASE,BASE),T());
          for (int j=ia[i]; j<ia[i+1] && f>>p.second; ++j) {
            p.first.j = ja[j-BASE];
            entries.insert(p);
          }
        }
        ia.clear();
        ja.clear();
      }


      // read format: Harwell-Boeing exchange format (*.rua)
      else if (hasdot && _fname.substr(_fname.find_last_of("."))==".rua") {
        throw runtime_error("sparse_matrix: Harwell-Boeing file format not implemented.");
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

    // force diagonal entries and structural symmetry
    {
      size_t nmodif = 0;
      for (size_t i=0; i<size.i; ++i)
        entries.insert(coord_t<T>(idx_t(i+BASE,i+BASE),T())).second? ++nmodif:nmodif;
      for (bool modif=true; modif;) {
        modif = false;
        for (typename set_t::const_reverse_iterator c=entries.rbegin(); !modif && c!=entries.rend(); ++c) {
          coord_t<T> p(idx_t(c->first.j,c->first.i),T());
          (modif = entries.insert(p).second)? ++nmodif:nmodif;
        }
      }
      CFdebug << "sparse_matrix: additional entries to preserve symmetry: " << nmodif << CFendl;
    }

    // compress and return
    compress();
    if (matrix_base_t::size(0)!=nnu ||
        matrix_base_t::size(1)<=*std::max_element(ja.begin(),ja.end())-BASE)
      throw std::runtime_error("sparse_matrix: after reading file, indexing not correct.");
    return *this;
  }

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

  sparse_matrix& clear() {
    matrix_base_t::m_size.clear();
    entries.clear();
    nnu = nnz = 0;
    ia.clear();
    ja.clear();
    a.clear();
    return *this;
  }

  sparse_matrix& operator=(const double& _value) {
    a.assign(a.size(),static_cast< const T >(_value));
    return *this;
  }

  sparse_matrix& operator=(const sparse_matrix& _other) {
    matrix_base_t::m_size = _other.matrix_base_t::m_size;
    nnu = _other.nnu;
    nnz = _other.nnz;
    ia = _other.ia;
    ja = _other.ja;
    a  = _other.a;
    return *this;
  }

  sparse_matrix& zerorow(const size_t& _i) {
    //FIXME
    std::fill_n(a.begin()+ia[_i]-BASE,ia[_i+1]-ia[_i],T());
    return *this;
  }

  sparse_matrix& swap(sparse_matrix& _other) {
    matrix_base_t::swap(_other);
    std::swap(nnu,_other.nnu);
    std::swap(nnz,_other.nnz);
    ia.swap(_other.ia);
    ja.swap(_other.ja);
    a .swap(_other.a);
    return *this;
  }

  std::ostream& print(std::ostream& o) const {
    using namespace std;
    const idx_t& size = matrix_base_t::m_size;

    const T eps = 1.e3*numeric_limits< T >::epsilon();
    const print_t print_level(matrix_base_t::m_print?   min(matrix_base_t::m_print,print_full) :
                             (size.i>100 || size.j>100? print_size  :
                             (size.i> 10 || size.j> 10? print_signs :
                                                        print_full )));
    o << "(" << size.i << 'x' << size.j << ">=" << a.size() << ") [ ";
    if      (print_level==print_size)  { o << "..."; }
    else if (print_level==print_signs) {
      string str;
      for (size_t i=0; i<size.i; ++i) {
        str.assign(size.j,' ');
        for (int k=ia[i]-BASE; k<ia[i+1]-BASE; ++k)
          str[ ja[k]-BASE ] = (a[k]>eps? '+' : (a[k]<-eps? '-' : '.' ));
        o << "\n  " << str;
      }
    }
    else if (print_level==print_full) {
      for (size_t i=0; i<size.i; ++i) {
        vector< T > row(size.j,T());
        for (int k=ia[i]-BASE; k<ia[i+1]-BASE; ++k)
          row[ ja[k]-BASE ] = a[k];
        o << "\n  ";
        copy(row.begin(),row.end(),ostream_iterator< T >(o,", "));
      }
    }
    return o << " ]";
  }

  void print(std::string& _fname) const {
    using namespace std;
    ofstream f(_fname.c_str());
    if (!f)
      throw runtime_error("sparse_matrix: cannot open file.");

    const bool hasdot(string("."+_fname).find_last_of("."));
    if (hasdot && _fname.substr(_fname.find_last_of("."))==".csr") {

      // compress to temporary memory if necessary (to go around the "const" qualifier)
      vector< int > another_ia;
      vector< int > another_ja;
      vector< T   > another_a;
      if (!is_compressed()) {
        SORT::compress(entries,another_ia,another_ja);
        for (typename set_t::const_iterator c=entries.begin(); c!=entries.end(); ++c)
          another_a.push_back(c->second);
      }
      const vector< int >& my_ia = is_compressed()? ia : another_ia;
      const vector< int >& my_ja = is_compressed()? ja : another_ja;
      const vector< T   >& my_a  = is_compressed()? a  : another_a;

      // output
      f << matrix_base_t::size(0) << ' ' << matrix_base_t::size(1) << '\n';
      copy(my_ia.begin(),my_ia.end(),ostream_iterator< int >(f," ")); f << '\n';
      copy(my_ja.begin(),my_ja.end(),ostream_iterator< int >(f," ")); f << '\n';
      copy(my_a .begin(),my_a .end(),ostream_iterator< T   >(f," ")); f << '\n';
    }

    else {
      ostringstream msg;
      msg << "sparse_matrix: cannot write file \"" << _fname << "\" (only \"*.csr\" is supported).";
      throw runtime_error(msg.str());
    }
  }


  // indexing
  const T& operator()(const size_t& _i, const size_t& _j) const {
    if (is_compressed()) {
      const int i=getindex(_i,_j);
      return (i<0? matrix_base_t::m_zero:a[i]);
    }
    const coord_t<T> p(idx_t(_i+BASE,_j+BASE),T());
    typename set_t::const_iterator i = entries.find(p);
    return i!=entries.end()? i->second : matrix_base_t::m_zero;
  }

  T& operator()(const size_t& _i, const size_t& _j) {
    if (is_compressed()) {
      const int i=getindex(_i,_j);
      return (i<0? matrix_base_t::m_zero:a[i]);
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
    //FIXME: this is specific to CSR structure! this container should also work for CSC!
    for (int k=ia[i]-BASE; k<ia[i+1]-BASE; ++k)
      if (ja[k]-BASE==(int) j)
        return k;
    return -1;
  }


  // compression/uncompression

  inline bool is_compressed() const { return nnu; }

  void compress() {
    if (is_compressed())
      return;
    ia.clear();
    ja.clear();
    a .clear();  a.reserve(nnz);

    nnu = SORT::compress(entries,ia,ja);
    nnz = entries.size();
    for (typename set_t::const_iterator c=entries.begin(); c!=entries.end(); ++c)
      a.push_back(c->second);

    entries.clear();
  }

  void uncompress() {
    if (!is_compressed())
      return;
    entries.clear();
    for (int i=0; i<ia.size()-1; ++i)
      for (int j=ia[i]; j<ia[i+1]; ++j)
        entries.insert(coord_t<T>(idx_t(i+BASE,ja[j-BASE]),a[j-BASE]));
    ia.clear();
    ja.clear();
    a .clear();
    nnu = nnz = 0;
  }


  // storage (uncompressed)
  set_t entries;

  // storage (compressed)
  std::vector< T > a;
  std::vector< int > ia;
  std::vector< int > ja;
  int nnu;
  int nnz;
};


}  // namespace detail
}  // namespace lss
}  // namespace cf3


#endif

