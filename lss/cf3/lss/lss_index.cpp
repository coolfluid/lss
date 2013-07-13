// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "lss_index.hpp"
#include <cstdarg>  // NOTE: below the above because it's a "very ugly" header


namespace cf3 {
namespace lss {
namespace lss_index {


// NOTE: I opted for the constructor registering multiple methods, thus using
// variable arguments, because source files can only contain definitions (no
// declarations) so I can't call member methods, overload operators (), << or >>
// to register all entries into the factory instance outside of run-time.

/// constructor registering built-in indexing techniques
index_factory_t::index_factory_t(const size_t n, index_creator_t* c, ...) {
  if (n==1 && c!=NULL) Register(c);
  if (n<2) return;

  va_list ap;
  va_start(ap,c);
  for (size_t i=1; i<n; ++i)
    Register( va_arg(ap,index_creator_t*) );
  va_end(ap);
}


#if 0
// NOTE: maybe this should be addressed directly because it's internal storage and we don't want users to mess around these things?
// register built-in matrix storage indexing techniques
index_factory_t idx_factory_matrix_storage(
    1,
    new index_creator< idx_sparse_matrix_csr >("sparse_matrix_csr_indexing")
    );
#endif


// register built-in linear system matrix addressing indexing techniques
index_factory_t idx_factory_matrix_indexing(
    1,
    new index_creator< idx_regular_block >("regular_block_indexing")
    );


}  // namespace lss_index
}  // namespace lss
}  // namespace cf3
