// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>


#include "common/Log.hpp"
#include "utilities.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- coordinate matrix helper structures ----------------------------------- */

/**
 * @brief Coordinate matrix entry
 */
typedef std::pair< idx_t, double > coord_t;


/**
 * @brief Coordinate matrix sorting tool, by row (functor)
 */
struct coord_ordering_by_row_t
{
  bool operator()(const coord_t& a, const coord_t& b) const {
    return (a.first.i<b.first.i? true  :
           (a.first.i>b.first.i? false : (a.first.j<b.first.j)));
  }
};


/**
 * @brief Coordinate matrix sorting tool, by column (functor)
 */
struct coord_ordering_by_column_t
{
  bool operator()(const coord_t& a, const coord_t& b) const {
    return (a.first.j<b.first.j? true  :
           (a.first.j>b.first.j? false : (a.first.i<b.first.i)));
  }
};


/**
 * @brief Coordinate matrix row compression tool (functor)
 */
struct coord_row_equal_to_t
{
  const size_t i;
  coord_row_equal_to_t(const size_t& _i) : i(_i) {}
  bool operator()(const coord_t& a) const { return i==a.first.i; }
};


/**
 * @brief Coordinate matrix column compression tool (functor)
 */
struct coord_column_equal_to_t
{
  const size_t j;
  coord_column_equal_to_t(const size_t& _j) : j(_j) {}
  bool operator()(const coord_t& a) const { return j==a.first.j; }
};


/* -- MatrixMarket I/O helper structures ------------------------------------ */

namespace MatrixMarket {


namespace mm {


  // matrix typecodes
  enum type     { no_type, matrix, all_types };
  enum format   { no_format, coordinate, array, all_formats };
  enum field    { no_field, real, integer, cmplex, pattern, all_fields };
  enum symmetry { no_symmetry, general, symmetric, skew_symmetric, hermitian, all_symmetries };


  // matrix typecode descriptions
  const char* str_type[]     = { "?", "matrix" };
  const char* str_format[]   = { "?", "coordinate", "array" };
  const char* str_field[]    = { "?", "real", "integer", "complex", "pattern" };
  const char* str_symmetry[] = { "?", "general", "symmetric", "skew_symmetric", "hermitian" };


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


  // read file utility: process file header
  bool read_banner(std::ifstream& f, mm::typecode_t& t) {
    std::string line, word;
    std::getline(f,line);
    if (line.find_first_of("%%MatrixMarket")==0) {
      t.clear();
      std::istringstream lstream(line);

      lstream >> word /*"%%MatrixMarket"*/ >> word;
      for (int i=mm::no_type; i<mm::all_types; ++i)
        if (word==mm::str_type[i])
          t.set_type((mm::type) i);

      lstream >> word;
      for (int i=mm::no_format; i<mm::all_formats; ++i)
        if (word==mm::str_format[i])
          t.set_format((mm::format) i);

      lstream >> word;
      for (int i=mm::no_field; i<mm::all_fields; ++i)
        if (word==mm::str_field[i])
          t.set_field((mm::field) i);

      lstream >> word;
      for (int i=mm::no_symmetry; i<mm::all_symmetries; ++i)
        if (word==mm::str_symmetry[i])
          t.set_symmetry((mm::symmetry) i);

    }
    return t.is_valid();
  }


  // read file utility: process matrix/array size
  bool read_size(std::ifstream& f, size_t& i, size_t& j, size_t& nz) {
    std::string line, word;
    i = j = nz = 0;
    while (std::getline(f,line)) {
      if (line.find_first_of("%")==0 || !line.length()) {}
      else {
        std::istringstream lstream(line);
        lstream >> i >> j >> nz;
        break;
      }
    }
    return idx_t(i,j).is_valid_size();
  }


}  // namespace mm


bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > > &a )
{
  using namespace std;

  size.invalidate();
  a.clear();
  mm::typecode_t t;
  size_t nnz(0);

  ifstream f(fname.c_str());
  if (!f)                                   throw runtime_error("cannot open file.");
  if (!mm::read_banner(f,t))                throw runtime_error("MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
  if (!mm::read_size(f,size.i,size.j,nnz))  throw runtime_error("MatrixMarket: invalid matrix/array size.");
  if (!t.is_real() || !t.is_general())      throw runtime_error("MatrixMarket: only \"(coordinate|array) real general\" supported.");

  // read into row/column-oriented dense matrix, line by line
  a.assign(roworiented? size.i:size.j, vector< double >(
           roworiented? size.j:size.i, double()));
  idx_t ij(1,1);
  string line;
  double v;
  while (getline(f,line)) {
    if (line.find_first_of("%")==0) {}
    else if (t.is_dense()) {
      istringstream(line) >> v;
      if (roworiented) a[ ij.i-1 ][ ij.j-1 ] = v;
      else             a[ ij.j-1 ][ ij.i-1 ] = v;
      if (++ij.i > size.i) {
        ij.i = 1;
        ij.j++;
      }
    }
    else {
      istringstream(line) >> ij.i >> ij.j >> v;
      if (roworiented) a[ ij.i-1 ][ ij.j-1 ] = v;
      else             a[ ij.j-1 ][ ij.i-1 ] = v;
    }
  }
  return true;
}


bool read_sparse(
    const std::string& fname,
    const bool& roworiented,
    const int& base,
    idx_t& size,
    std::vector< double >& a,
    std::vector< int >& ia,
    std::vector< int >& ja )
{
  using namespace std;

  size.invalidate();
  ia.clear();
  ja.clear();
  a.clear();

  mm::typecode_t t;
  size_t nnz(0);
  ifstream f(fname.c_str());
  if (!f)                                   throw runtime_error("cannot open file.");
  if (!mm::read_banner(f,t))                throw runtime_error("MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
  if (!mm::read_size(f,size.i,size.j,nnz))  throw runtime_error("MatrixMarket: invalid matrix/array size.");
  if (!t.is_real() || !t.is_general())      throw runtime_error("MatrixMarket: only \"(coordinate|array) real general\" supported.");

  coord_t p = make_pair(idx_t(1,1),0.);
  string line;
  size_t fixes(0);
  if (roworiented) {

    // read into ordered set, line by line
    typedef set< coord_t, coord_ordering_by_row_t > set_t;
    set_t entries;
    while (getline(f,line)) {
      if (line.find_first_of("%")==0) {}
      else if (t.is_dense()) {
        istringstream(line) >> p.second;
        entries.insert(p);
        if (++p.first.i > size.i) {
          p.first.i = 1;
          p.first.j++;
        }
      }
      else {
        istringstream(line) >> p.first.i >> p.first.j >> p.second;
        entries.insert(p);
      }
    }

    // force diagonal entries and structural symmetry
    for (int i=1; i<=size.i; ++i)
      entries.insert(make_pair(idx_t(i,i),0.)).second? ++fixes:fixes;
    for (bool modif=true; modif;) {
      modif = false;
      for (set_t::const_reverse_iterator c=entries.rbegin(); !modif && c!=entries.rend(); ++c)
        (modif = entries.insert(make_pair(idx_t(c->first.j,c->first.i),0.)).second)? ++fixes:fixes;
    }

    // compress each row, then set indices and values
    ia.reserve(size.i+1);
    ja.reserve(entries.size());
    a .reserve(entries.size());
    ia.push_back(1);  // (base is always 1 in this format)
#if 0
    for (size_t r=1; r<=size.i; ++r)
      ia.push_back(ia.back() + count_if(entries.begin(),entries.end(),coord_row_equal_to_t(r)));
#else
    set_t::const_iterator c=entries.begin();
    for (size_t r=1, count; r<=size.i; ++r) {
      for (count=0; r==(c->first.i) && c!=entries.end(); ++c, ++count) {}
      ia.push_back(ia.back() + count);
    }
#endif
    for (set_t::const_iterator c=entries.begin(); c!=entries.end(); ++c) {
      ja.push_back(c->first.j);
      a.push_back(c->second);
    }

  }
  else {

    // read into ordered set, line by line
    typedef set< coord_t, coord_ordering_by_column_t > set_t;
    set_t entries;
    while (getline(f,line)) {
      if (line.find_first_of("%")==0) {}
      else if (t.is_dense()) {
        istringstream(line) >> p.second;
        entries.insert(p);
        if (++p.first.i > size.i) {
          p.first.i = 1;
          p.first.j++;
        }
      }
      else {
        istringstream(line) >> p.first.i >> p.first.j >> p.second;
        entries.insert(p);
      }
    }

    // force diagonal entries and structural symmetry
    for (int i=1; i<=size.i; ++i)
      entries.insert(make_pair(idx_t(i,i),0.)).second? ++fixes:fixes;
    for (bool modif=true; modif;) {
      modif = false;
      for (set_t::const_reverse_iterator c=entries.rbegin(); !modif && c!=entries.rend(); ++c)
        (modif = entries.insert(make_pair(idx_t(c->first.j,c->first.i),0.)).second)? ++fixes:fixes;
    }

    // compress each column, then set indices and values
    ja.reserve(size.i+1);
    ia.reserve(entries.size());
    a .reserve(entries.size());
    ja.push_back(1);  // (base is always 1 in this format)
#if 0
    for (size_t c=1; c<=size.j; ++c)
      ja.push_back(ja.back() + count_if(entries.begin(),entries.end(),coord_column_equal_to_t(c)));
#else
    set_t::const_iterator c=entries.begin();
    for (size_t j=1, count; j<=size.j; ++j) {
      for (count=0; j==(c->first.j) && c!=entries.end(); ++c, ++count) {}
      ja.push_back(ja.back() + count);
    }
#endif
    for (set_t::const_iterator c=entries.begin(); c!=entries.end(); ++c) {
      ia.push_back(c->first.i);
      a.push_back(c->second);
    }

  }
  if (fixes) {
    CFdebug << "MatrixMarket::read_sparse: additional entries to preserve symmetry: " << fixes << CFendl;
  }

  // conversion to intended indexing base
  for_each(ja.begin(),ja.end(),base_conversion_t(base-1));
  for_each(ia.begin(),ia.end(),base_conversion_t(base-1));
  return true;
}


}  // namespace MatrixMarket


/* -- CSR I/O helper structures --------------------------------------------- */

namespace CSR {


bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > > &a )
{
  using namespace std;

  size.invalidate();
  a.clear();

  // open file and read matrix size
  string line;
  ifstream f(fname.c_str());
  if (!f) throw runtime_error("cannot open file.");
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
  a.assign(roworiented? size.i:size.j,vector< double >(
           roworiented? size.j:size.i,double() ));
  for (int i=0; i<static_cast< int >(size.i); ++i)
    for (int k=ia[i]; k<ia[i+1]; ++k)
      f >> a[ roworiented? i:ja[k] ]
            [ roworiented? ja[k]:i ];

  return true;
}


bool read_sparse(
    const std::string& fname,
    const int& base,
    idx_t& size,
    std::vector< double >& a,
    std::vector< int >& ia,
    std::vector< int >& ja )
{
  using namespace std;

  size.invalidate();
  a.clear();
  ia.clear();
  ja.clear();

  // open file and read matrix size
  ifstream f(fname.c_str());
  {
    string line;
    if (!f)
      return false;
    while (!size.is_valid_size()) {
      if (getline(f,line) && line.find_first_of("%")!=0)
        istringstream(line) >> size.i >> size.j;
    }
  }

  // read ia and ja arrays and convert to intended base
  // (ia has indices in increasing order, ja first entry read while building ia)
  ia.reserve(size.i+1);
  ja.push_back(0);
  for (int &i=ja[0], j=-1; f>>i && j<i;)
    ia.push_back(j=i);

  const int
      nnu(ia.size()-1),
      nnz(ia.back()-ia.front());

  ja.reserve(nnz);
  for (int i; ja.size()<nnz && f>>i;)
    ja.push_back(i);

  for_each(ja.begin(),ja.end(),base_conversion_t(base-ia.front()));
  for_each(ia.begin(),ia.end(),base_conversion_t(base-ia.front()));

  // create a coordinate matrix, forcing diagonal entry and structural symmetry
  typedef set< coord_t, coord_ordering_by_row_t > set_t;
  set_t entries;
  for (int i=0; i<nnu; ++i) {
    coord_t p = make_pair(idx_t(),0.);
    for (int j=ia[i]; j<ia[i+1] && f>>p.second; ++j) {      
      p.first.i = i + base;
      p.first.j = ja[j-base];
      entries.insert(p);
    }
  }

  // force diagonal entries and structural symmetry
  size_t fixes(0);
  for (int i=0; i<nnu; ++i)
    entries.insert(make_pair(idx_t(i+base,i+base),0.)).second? ++fixes:fixes;
  for (bool modif=true; modif;) {
    modif = false;
    for (set_t::const_reverse_iterator c=entries.rbegin(); !modif && c!=entries.rend(); ++c)
      (modif = entries.insert(make_pair(idx_t(c->first.j,c->first.i),0.)).second)? ++fixes:fixes;
  }
  if (fixes) {
    CFdebug << "CSR::read_sparse: additional entries to preserve symmetry: " << fixes << CFendl;
  }

  // recompress row indices, and rebuild column indices and matrix values
  ia.clear();  ia.reserve(nnu);
  ja.clear();  ja.reserve(entries.size());
               a .reserve(entries.size());
  ia.push_back(base);
#if 0
  for (int i=0; i<nnu; ++i)
    ia.push_back(ia.back() + count_if(entries.begin(),entries.end(),coord_row_equal_to_t(i+base)));
#else
  set_t::const_iterator c=entries.begin();
  for (int i=0, count; i<nnu; ++i) {
    for (count=0; (i+base)==(c->first.i) && c!=entries.end(); ++c, ++count) {}
    ia.push_back(ia.back() + count);
  }
#endif
  for (set_t::const_iterator c=entries.begin(); c!=entries.end(); ++c) {
    ja.push_back(c->first.j);
    a .push_back(c->second);
  }

  return true;
}


}  // namespace CSR


}  // namespace detail
}  // namespace lss
}  // namespace cf3
