// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "matrix.hpp"


namespace cf3 {
namespace lss {


/* -- matrix interface and implementations ---------------------------------- */


// utilities
print_t print_level(const int& i) {
  return std::max(print_auto,std::min(print_file,static_cast< print_t >(i)));
}


}  // namespace lss
}  // namespace cf3
