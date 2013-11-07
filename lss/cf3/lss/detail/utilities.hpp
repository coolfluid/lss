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

#include "index.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- generic utilities ----------------------------------------------------- */

/// @brief Type comparison at compilation time
template< typename A, typename B >
bool type_is_equal() { return (typeid(A)==typeid(B)); }


/// @brief Type conversion functor (POD to POD, no fancy stuff)
template< typename Tin, typename Tout >
struct type_conversion_t
{
  Tout operator()(const Tin& in) { return static_cast< Tout >(in); }
};


/* -- indexing conversion/application types --------------------------------- */

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

/// @brief Indexing base conversion tool (functor)
struct base_conversion_t
{
  base_conversion_t(const int& _diff) : diff (_diff) {}
  int& operator()(int& v) { return (v+=diff); }
  const int& diff;
};


/// @brief Elementary vector transformation: sort, removing duplicate entries
struct vector_sort_unique
{
  void apply(std::vector< int >& v, int) {
    std::sort(v.begin(),v.end());
    v.erase(std::unique(v.begin(),v.end()),v.end());
  }
};


/// @brief Elementary vector transformation: add vector entry at start
struct vector_element_push_front
{
  void apply(std::vector< int >& v, int e) { v.insert(v.begin(),e); }
};


/// @brief Elementary vector transformation: add vector entry at end
struct vector_element_push_back
{
  void apply(std::vector< int >& v, int e) { v.push_back(e); }
};


/// @brief Elementary vector transformation: remove entries of specific value, if existing
struct vector_element_remove
{
  struct equal_to
  {
    equal_to(const int& _elem) : elem(_elem) {}
    bool operator()(const int& _elem) const { return elem==_elem; }
    const int elem;
  };
  void apply(std::vector< int >& v, int e) {
    v.erase(std::remove_if(v.begin(),v.end(),equal_to(e)),v.end());
  }
};


/// @brief Elementary vector transformation: apply supplied difference to vector
struct vector_apply_diff
{
  struct apply_diff
  {
    apply_diff(const int& _diff) : diff(_diff) {}
    int operator()(int &v) const { return (v+diff); }
    const int diff;
  };
  void apply(std::vector< int >& v, int e) {
    std::for_each(v.begin(),v.end(),apply_diff(e));
  }
};


/// @brief Type for recursive application of vector transformations
struct vector_transform_t
{
  virtual void apply(std::vector< int >&, int) {}
};


/// @brief Type for recursive application of vector transformations
template<
    typename Transf,
    typename TransfNested=vector_transform_t >
struct vector_transform_list_t : vector_transform_t
{
  void apply(std::vector< int >& v, int e) {
    TransfNested::apply(v,e);
    Transf::apply(v,e);
  }
};


/**
 * @brief Vector transformation: sorted indices vector
 * (useful for specific CSR/CSC matrix linear solvers)
 */
typedef vector_transform_list_t< vector_sort_unique,
        vector_transform_list_t< vector_element_push_back > >
  vector_sorted_with_diagonal_t;

/**
 * @brief Vector transformation: sorted indices vector with specific entry at start
 * (useful for specific CSR/CSC matrix linear solvers)
 */
typedef vector_transform_list_t< vector_element_push_front,
        vector_transform_list_t< vector_sort_unique,
        vector_transform_list_t< vector_element_remove > > >
  vector_sorted_with_diagonal_first_t;


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
 * @param base: sparse structure index base
 * @param size: output matrix/array size
 * @param a: output sparse matrix/array
 * @param ia: output i-coordinate indices
 * @param ja: output j-coordinate indices
 * @return if reading is successful
 */
bool read_sparse(const std::string& fname,
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
void read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< T > >& a )
{
  // if storage is of different type than double, it miraculously converts
  std::vector< std::vector< double > > another_a;
  std::vector< std::vector< double > >& storage(
    type_is_equal< T, double >()? (std::vector< std::vector< double > >&) a
                                : another_a);

  // read file contents
  const bool hasdot(std::string("."+fname).find_last_of("."));
  if      (hasdot && fname.substr(fname.find_last_of("."))==".mtx") { MatrixMarket ::read_dense(fname,roworiented,size,storage); }
  else if (hasdot && fname.substr(fname.find_last_of("."))==".csr") { CSR          ::read_dense(fname,roworiented,size,storage); }
/*else if (hasdot && fname.substr(fname.find_last_of("."))==".rua") { HarwellBoeing::read_dense(fname,roworiented,size,storage); }*/
  else
    throw std::runtime_error("file format not detected.");

  // perform storage conversion if necessary, and return
  if (!type_is_equal< T, double >()) {
    a.assign(roworiented? size.i:size.j,std::vector< T >(
             roworiented? size.j:size.i,T() ));
    for (size_t i=0; i<another_a.size(); ++i)
      transform( another_a[i].begin(),another_a[i].end(),a[i].begin(),
                 type_conversion_t< double, T >() );
  }
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
void read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< T >& a,
    std::vector< int >& ia,
    std::vector< int >& ja )
{
  // if storage is of different type than double, it miraculously converts
  std::vector< double > another_a;
  std::vector< double >& storage(
    type_is_equal< T, double >()? (std::vector< double >&) a : another_a);

  // read file contents
  const bool hasdot(std::string("."+fname).find_last_of("."));
  if      (hasdot && fname.substr(fname.find_last_of("."))==".mtx") { MatrixMarket ::read_sparse(fname,roworiented,base,size,storage,ia,ja); }
  else if (hasdot && fname.substr(fname.find_last_of("."))==".csr") { CSR          ::read_sparse(fname,            base,size,storage,ia,ja); }
/*else if (hasdot && fname.substr(fname.find_last_of("."))==".rua") { HarwellBoeing::read_sparse(fname,roworiented,base,size,storage,ia,ja); }*/
  else
    throw std::runtime_error("file format not detected.");

  // perform storage conversion if necessary, and return
  if (!type_is_equal< T, double >()) {
    a.assign(another_a.size(),T());
    transform( another_a.begin(),another_a.end(),a.begin(),
               type_conversion_t< double, T >() );
  }
}


}  // namespace detail
}  // namespace lss
}  // namespae cf3


#endif
