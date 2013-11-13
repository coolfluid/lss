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


/* -- MatrixMarket I/O helper structures ------------------------------------ */

namespace MatrixMarket {


bool read_banner(std::ifstream& f, typecode_t& t) {

  const char
    *str_type[]     = { "?", "matrix" },
    *str_format[]   = { "?", "coordinate", "array" },
    *str_field[]    = { "?", "real", "integer", "complex", "pattern" },
    *str_symmetry[] = { "?", "general", "symmetric", "skew_symmetric", "hermitian" };

  std::string line, word;
  std::getline(f,line);
  if (line.find_first_of("%%MatrixMarket")==0) {
    t.clear();
    std::istringstream lstream(line);

    lstream >> word /*"%%MatrixMarket"*/ >> word;
    for (int i=no_type; i<all_types; ++i)
      if (word==str_type[i])
        t.set_type((type) i);

    lstream >> word;
    for (int i=no_format; i<all_formats; ++i)
      if (word==str_format[i])
        t.set_format((format) i);

    lstream >> word;
    for (int i=no_field; i<all_fields; ++i)
      if (word==str_field[i])
        t.set_field((field) i);

    lstream >> word;
    for (int i=no_symmetry; i<all_symmetries; ++i)
      if (word==str_symmetry[i])
        t.set_symmetry((symmetry) i);

  }
  return t.is_valid();
}


bool read_size(std::ifstream& f, size_t& i, size_t& j, int& nz) {
  std::string line;
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


bool read_dense(
    const std::string& fname,
    const bool &roworiented,
    idx_t& size,
    std::vector< std::vector< double > > &a )
{
  using namespace std;

  size.invalidate();
  a.clear();
  typecode_t t;
  int nnz(0);

  ifstream f(fname.c_str());
  if (!f)                               throw runtime_error("cannot open file.");
  if (!read_banner(f,t))                throw runtime_error("MatrixMarket: invalid header, \"%%MatrixMarket ...\" not found.");
  if (!read_size(f,size.i,size.j,nnz))  throw runtime_error("MatrixMarket: invalid matrix/array size.");
  if (!t.is_real() || !t.is_general())  throw runtime_error("MatrixMarket: only \"(coordinate|array) real general\" supported.");

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


}  // namespace CSR


}  // namespace detail
}  // namespace lss
}  // namespace cf3
