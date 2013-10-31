// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "matrix.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- matrix interface and implementations ---------------------------------- */


// utilities
print_t print_level(const int& i) {
  return i>=print_auto || i<=print_full? static_cast< print_t >(i) : print_auto;
}

}  // namespace detail
}  // namespace lss
}  // namespace cf3
