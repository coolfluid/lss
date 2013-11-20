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


/* -- MatrixMarket I/O helper structures ------------------------------------ */

namespace MatrixMarket {


bool read_banner(std::ifstream& f, typecode_t& t) {

  const char
    *str_type[]     = { "?", "matrix" },
    *str_format[]   = { "?", "coordinate", "array" },
    *str_field[]    = { "?", "real", "integer", "complex", "pattern" },
    *str_symmetry[] = { "?", "general", "symmetric", "skew_symmetric", "hermitian" };

  std::string line, word[5];
  std::getline(f,line);
  if (line.find_first_of("%%MatrixMarket")==0) {
    std::istringstream lstream(line);
    lstream >> word[0] /*"%%MatrixMarket"*/
            >> word[1]
            >> word[2]
            >> word[3]
            >> word[4];
    for (int i=1; i<all_types;      ++i) if (word[1]==str_type    [i]) t.m_type     = (type)     i;
    for (int i=1; i<all_formats;    ++i) if (word[2]==str_format  [i]) t.m_format   = (format)   i;
    for (int i=1; i<all_fields;     ++i) if (word[3]==str_field   [i]) t.m_field    = (field)    i;
    for (int i=1; i<all_symmetries; ++i) if (word[4]==str_symmetry[i]) t.m_symmetry = (symmetry) i;
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


}  // namespace MatrixMarket


}  // namespace lss
}  // namespace cf3

