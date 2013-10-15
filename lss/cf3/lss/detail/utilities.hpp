// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_detail_utilities_hpp
#define cf3_lss_detail_utilities_hpp


#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>


namespace cf3 {
namespace lss {
namespace detail {


/* -- basic types for index manipulation ------------------------------------ */


/// @brief Index pair, fundamental dereferencing type
struct idx_t
{
  size_t i, j;
  idx_t(const size_t& _i=std::numeric_limits< size_t >::max(),
        const size_t& _j=std::numeric_limits< size_t >::max()) : i(_i), j(_j) {}

  bool operator<  (const idx_t& _other) const { return (i<_other.i? true : i>_other.i? false : (j<_other.j)); }
  bool operator>  (const idx_t& _other) const { return idx_t(_other)<*this; }
  bool operator== (const idx_t& _other) const { return i==_other.i && j==_other.j; }
  bool operator!= (const idx_t& _other) const { return i!=_other.i || j!=_other.j; }

  idx_t& invalidate() { return (*this = idx_t()); }
  bool is_valid_size()  const { return operator>(idx_t(0,0)) && operator<(idx_t()); }
  bool is_square_size() const { return i==j; }
  bool is_diagonal()    const { return is_square_size(); }
};


/* -- indexing conversion/application types --------------------------------- */

/// @brief Indexing base conversion tool (functor)
struct base_conversion_t
{
  base_conversion_t(const int& _diff) : diff (_diff) {}
  int& operator()(int& v) { return (v+=diff); }
  const int& diff;
};


/// @brief Storage type conversion (POD to POD, no fancy stuff) (functor)
template< typename Tin, typename Tout >
struct storage_conversion_t
{
  Tout operator()(const Tin& in) { return static_cast< Tout >(in); }
};


/// @brief Index conversion base class (pure abstract)
struct index_conversion_t
{
  virtual ~index_conversion_t() {}
  virtual void clear() = 0;
  virtual size_t size()                   const = 0;
  virtual size_t size(const size_t& d)    const = 0;
  virtual idx_t& dereference(idx_t& _idx) const = 0;
};


/// @brief Composite indexing: index hierarchy type list termination type
struct index_hierarchy_t_end
{
  idx_t& dereference(idx_t& _idx) const { return _idx; }
  void initialize() {}
};


/// @brief Composite indexing: index hierarchy definition
template<
    typename Idx,
    typename IdxNested=index_hierarchy_t_end >
struct index_hierarchy_t
{
  idx_t& dereference(idx_t& _idx) const { return IdxNested::dereference(Idx::dereference(_idx)); }
  void initialize() { Idx::initialize(); IdxNested::initialize(); }
};


/* -- vector transformations by type ---------------------------------------- */
/* (common operations for building row or column-oriented sparsity patterns)  */

/// @brief Type for recursive application of vector transformations
struct transform_list_t_end
{
  static void apply(std::vector< size_t >&, size_t) {}
};


/// @brief Type for recursive application of vector transformations
template<
    typename Transf,
    typename TransfNested=transform_list_t_end >
struct transform_list_t
{
  static void apply(std::vector< size_t >& v, int e) {
    TransfNested::apply(v,e);
    Transf::apply(v,e);
  }
};


/// @brief Vector transformation: sort, removing duplicate entries
struct vector_sort_unique
{
  static void apply(std::vector< size_t >& v, int) {
    std::sort(v.begin(),v.end());
    v.erase(std::unique(v.begin(),v.end()),v.end());
  }
};


/// @brief Vector transformation: add vector entry at start
struct vector_element_push_front
{
  static void apply(std::vector< size_t >& v, size_t e) { v.insert(v.begin(),e); }
};


/// @brief Vector transformation: add vector entry at end
struct vector_element_push_back
{
  static void apply(std::vector< size_t >& v, int e) { v.push_back(e); }
};


/// @brief Vector transformation: remove entries of specific value, if existing
struct vector_element_remove
{
  struct equal_to
  {
    equal_to(const size_t& _elem) : elem(_elem) {}
    bool operator()(const size_t& _elem) const { return elem==_elem; }
    const size_t elem;
  };
  static void apply(std::vector< size_t >& v, int e) {
    v.erase(std::remove_if(v.begin(),v.end(),
      equal_to(static_cast< size_t >(e))),v.end());
  }
};


/// @brief Vector transformation: sort, removing duplicate entries
struct vector_add_value
{
  struct add_value
  {
    add_value(const int& _diff) : diff(_diff) {}
    size_t operator()(size_t &v) const {
      return static_cast< size_t >(static_cast< int >(v)+diff);
    }
    const int diff;
  };
  static void apply(std::vector< size_t >& v, int e) {
    std::for_each(v.begin(),v.end(),add_value(e));
  }
};


/**
 * @brief Vector transformation: sorted indices vector
 * (useful for specific CSR/CSC matrix linear solvers)
 */
typedef transform_list_t< vector_sort_unique,
        transform_list_t< vector_element_push_back > >
  vector_sorted_t;

/**
 * @brief Vector transformation: sorted indices vector with specific entry at start
 * (useful for specific CSR/CSC matrix linear solvers)
 */
typedef transform_list_t< vector_element_push_front,
        transform_list_t< vector_sort_unique,
        transform_list_t< vector_element_remove > > >
  vector_sorted_diagonal_first_t;


/* -- Matrix Market I/O (or, say, just I) ----------------------------------- */

namespace MatrixMarket
{


/**
 * @brief read_dense: read MatrixMarket file to dense structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param size: output matrix/array size
 * @param a: output dense matrix/array
 * @return if reading is successful
 */
bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > >& a );


/**
 * @brief read_sparse: read MatrixMarket file to sparse structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param base: sparse structure index base
 * @param size: output matrix/array size
 * @param a: output sparse matrix/array
 * @param ia: output i-coordinate indices
 * @param ja: output j-coordinate indices
 * @return if reading is successful
 */
bool read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< double >& a,
    std::vector< int >& ia,
    std::vector< int >& ja );


}


/* -- Harwell-Boeing I/O (or, say again, just I) ---------------------------- */

namespace HarwellBoeing
{


/**
 * @brief read_dense: read Harwell-Boeing file format to dense structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param size: output matrix/array size
 * @param a: output dense matrix/array
 * @return if reading is successful
 */
bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > >& a );


/**
 * @brief read_sparse: read Harwell-Boeing file format to sparse structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param base: sparse structure index base
 * @param size: output matrix/array size
 * @param a: output sparse matrix/array
 * @param ia: output i-coordinate indices
 * @param ja: output j-coordinate indices
 * @return if reading is successful
 */
bool read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< double >& a,
    std::vector< int >& ia,
    std::vector< int >& ja );


}


/* -- CSR I/O (just I) ------------------------------------------------------ */

namespace CSR
{


/**
 * @brief read_dense: read CSR file format (a hack on MM) to dense structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param size: output matrix/array size
 * @param a: output dense matrix/array
 * @return if reading is successful
 */
bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > >& a );


/**
 * @brief read_sparse: read CSR file format (a hack on MM) to sparse structure
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param base: sparse structure index base
 * @param size: output matrix/array size
 * @param a: output sparse matrix/array
 * @param ia: output i-coordinate indices
 * @param ja: output j-coordinate indices
 * @return if reading is successful
 */
bool read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< double >& a,
    std::vector< int >& ia,
    std::vector< int >& ja );


}


/* -- Generic I/O (interfacing the above) ----------------------------------- */


/**
 * @brief read_dense: interface reading of files in different formats to dense
 * data structure, templasized with the container storage type
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param size: output matrix/array size
 * @param a: output dense matrix/array
 * @return if reading is successful
 */
template< typename T >
bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< T > >& a )
{
  // if storage is of different type than double, it miraculously converts
  const bool storage_conversion_needed(typeid(T)!=typeid(double));
  std::vector< std::vector< double > > another_a;
  std::vector< std::vector< double > >& storage(
        storage_conversion_needed? another_a
                                 : (std::vector< std::vector< double > >&) a);

  // read file contents
  if      (fname.substr(fname.find_last_of("."))==".mtx") { MatrixMarket ::read_dense(fname,roworiented,size,storage); }
  else if (fname.substr(fname.find_last_of("."))==".rua") { HarwellBoeing::read_dense(fname,roworiented,size,storage); }
  else if (fname.substr(fname.find_last_of("."))==".csr") { CSR          ::read_dense(fname,roworiented,size,storage); }
  else {
    throw std::runtime_error("file format not detected (\""+fname+"\").");
  }

  // perform storage conversion if necessary, and return
  if (storage_conversion_needed) {
    a.assign(roworiented? size.i:size.j,std::vector< T >(
             roworiented? size.j:size.i,T() ));
    for (size_t i=0; i<another_a.size(); ++i)
      transform( another_a[i].begin(),another_a[i].end(),a[i].begin(),
                 storage_conversion_t< double, T >() );
  }
  return true;
}


/**
 * @brief read_sparse: interface reading of files in different formats to sparse
 * data structure, templasized with the container storage type
 * @param fname: input filename
 * @param roworiented: if result should be row (most common) or column oriented
 * @param base: sparse structure index base
 * @param size: output matrix/array size
 * @param a: output sparse matrix/array
 * @param ia: output i-coordinate indices
 * @param ja: output j-coordinate indices
 * @return if reading is successful
 */
template< typename T >
bool read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< T >& a,
    std::vector< int >& ia,
    std::vector< int >& ja )
{
  // if storage is of different type than double, it miraculously converts
  const bool storage_conversion_needed(typeid(T)!=typeid(double));
  std::vector< double > another_a;
  std::vector< double >& storage(
        storage_conversion_needed? another_a : (std::vector< double >&) a);

  // read file contents
  if      (fname.substr(fname.find_last_of("."))==".mtx") { MatrixMarket ::read_sparse(fname,roworiented,base,size,storage,ia,ja); }
  else if (fname.substr(fname.find_last_of("."))==".rua") { HarwellBoeing::read_sparse(fname,roworiented,base,size,storage,ia,ja); }
  else if (fname.substr(fname.find_last_of("."))==".csr") { CSR          ::read_sparse(fname,roworiented,base,size,storage,ia,ja); }
  else {
    throw std::runtime_error("file format not detected (\""+fname+"\").");
  }

  // perform storage conversion if necessary, and return
  if (storage_conversion_needed) {
    a.assign(another_a.size(),T());
    transform( another_a.begin(),another_a.end(),a.begin(),
               storage_conversion_t< double, T >() );
  }
  return true;
}


}  // namespace detail
}  // namespace lss
}  // namespae cf3


#endif
