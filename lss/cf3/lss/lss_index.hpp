// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_lss_index_hpp
#define cf3_lss_lss_index_hpp


#include <vector>
#include <algorithm>


namespace cf3 {
namespace lss {
namespace lss_index {


/* -- indexing techniques --------------------------------------------------- */

/**
 * @brief Storage indexing in compressed sparse row
 * BASE: {0,1}-based indexing, other bases won't work
 */
template< int BASE >
struct index_compressed_sparse_row_t : index_conversion_t
{
  index_compressed_sparse_row_t() { clear(); }

#if 0
  void initialize(const index_compressed_sparse_row_t& x) {
    const std::vector< std::vector< size_t > >& nz;

    nz.clear();
    nz.reserve(_nz.size());
    for (size_t r=0; r<_nz.size(); ++r)
      lss_index::vector_sorted_t::apply(nz.back(),r);

    //FIXME
    const size_t Nb = 1;

    // set row indices
    nnu = (int) (Nb * nz.size());
    ia.resize(nnu+1);
    ia[0] = BASE;
    int k = 0;
    for (int R=0; R<(int) nz.size(); ++R)
      for (int i=0; i<(int) Nb; ++i, ++k)
        ia[k+1] = ia[k] + (int) Nb * (int) nz[R].size();

    // set column indices
    nnz = ia[nnu] - BASE;
    ja.resize(nnz);
    for (size_t R=0; R<nz.size(); ++R) {
      for (int r=0; r<(int) Nb; ++r) {
        k = ia[R*Nb+r] - BASE;
        for (size_t I=0; I<nz[R].size(); ++I)
          for (int c=0; c<(int) Nb; ++c)
            ja[k++] = (int) (Nb*nz[R][I]) + c + BASE;
      }
    }
  }
#endif

  void clear() {
    ia.clear();
    ja.clear();
    nnu = nnz = 0;
  }

  size_t size() const { return static_cast< size_t >(nnz); }
  size_t size(const size_t& d) const {
    return (d==0? static_cast< size_t >(nnu) :
                  std::numeric_limits< size_t >::max());
  }

  idx_t& dereference(idx_t& _idx) const {
    const int i = static_cast< const int >(_idx.i),
              j = static_cast< const int >(_idx.j);
    _idx.invalidate();
    for (int k=ia[i]-BASE; k<ia[i+1]-BASE; ++k)
      if (ja[k]-BASE==j) {
        _idx = idx_t(static_cast< size_t >(k),0);
        break;
      }
    return _idx;
  }

  // storage
  std::vector< int > ia, ja;
  int nnu, nnz;
};


}  // namespace lss_index
}  // namespace lss
}  // namespace cf3


#endif

