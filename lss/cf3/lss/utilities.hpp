// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_utilities_hpp
#define cf3_lss_utilities_hpp


#include <algorithm>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <fstream>
#include <set>
#include <sstream>

#include "index.hpp"


namespace cf3 {
namespace lss {


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


/// @brief Indexing base conversion tool (functor)
struct base_conversion_t
{
  base_conversion_t(const int& _diff) : diff (_diff) {}
  int& operator()(int& v) { return (v+=diff); }
  const int& diff;
};


/* -- Matrix Market I/O (or, say, just I) ----------------------------------- */

namespace MatrixMarket
{


// matrix typecodes
enum type     { no_type, matrix, all_types };
enum format   { no_format, coordinate, array, all_formats };
enum field    { no_field, real, integer, cmplex, pattern, all_fields };
enum symmetry { no_symmetry, general, symmetric, skew_symmetric, hermitian, all_symmetries };


// matrix typecode manipulator
struct typecode_t {

  type     m_type;
  format   m_format;
  field    m_field;
  symmetry m_symmetry;

  typecode_t() { clear(); }
  void clear() {
    m_type     = no_type;
    m_format   = no_format;
    m_field    = no_field;
    m_symmetry = general;
  }

  // -- typecode query functions

  bool is_matrix()     { return m_type==matrix; }

  bool is_sparse()     { return m_format==coordinate; }
  bool is_coordinate() { return m_format==coordinate; }
  bool is_dense()      { return m_format==array; }
  bool is_array()      { return m_format==array; }

  bool is_complex()    { return m_field==cmplex; }
  bool is_real()       { return m_field==real; }
  bool is_pattern()    { return m_field==pattern; }
  bool is_integer()    { return m_field==integer; }

  bool is_symmetric()  { return m_symmetry==symmetric; }
  bool is_general()    { return m_symmetry==general; }
  bool is_skew()       { return m_symmetry==skew_symmetric; }
  bool is_hermitian()  { return m_symmetry==hermitian; }

  bool is_valid() {
    return is_matrix()
        && !(is_complex())  // (NOTE: not supported in my implementation!)
        && !(is_dense() && is_pattern())
        && !(is_real() && is_hermitian())
        && !(is_pattern() && (is_hermitian() || is_skew()));
  }

  // -- typecode modify functions

  void set_type     (const type&     i) { m_type=i;     }
  void set_format   (const format&   i) { m_format=i;   }
  void set_field    (const field&    i) { m_field=i;    }
  void set_symmetry (const symmetry& i) { m_symmetry=i; }

};


// read file utility: process file header and matrix/array size
bool read_banner(std::ifstream& f, typecode_t& t);
bool read_size(std::ifstream& f, size_t& i, size_t& j, int& nz);


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


}  // namespace MatrixMarket


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


}  // namespace lss
}  // namespae cf3


#endif
