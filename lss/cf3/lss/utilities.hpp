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

  typecode_t() {
    m_type     = no_type;
    m_format   = no_format;
    m_field    = no_field;
    m_symmetry = general;
  }

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
        && !(is_complex())  // (not supported in this implementation)
        && !(is_dense() && is_pattern())
        && !(is_real() && is_hermitian())
        && !(is_pattern() && (is_hermitian() || is_skew()));
  }

};


// read file utility: process file header and matrix/array size
bool read_banner(std::ifstream& f, typecode_t& t);
bool read_size(std::ifstream& f, size_t& i, size_t& j, int& nz);


}  // namespace MatrixMarket


}  // namespace lss
}  // namespae cf3


#endif
