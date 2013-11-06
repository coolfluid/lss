#ifndef lss_lss_index_hpp
#define lss_lss_index_hpp


#include <vector>
#include <algorithm>


#include "lss_utilities.hpp"


/* -- indexing techniques --------------------------------------------------- */

namespace lss_index {


/**
 * @brief Type list for indexing hierarchy definition
 * (note: for instantiation, composition is used for nested indexing)
 */
struct index_hierarchy_t_end
{
  idx_t& dereference(idx_t& _idx) const { return _idx; }
  void initialize() {}
};

template<
    typename Idx,
    typename IdxNested=index_hierarchy_t_end >
struct index_hierarchy_t
{
  idx_t& dereference(idx_t& _idx) const { return IdxNested::dereference(Idx::dereference(_idx)); }
  void initialize() { Idx::initialize(); IdxNested::initialize(); }
};


/**
 * @brief Index conversion base class (pure abstract)
 */
struct index_t
{
  virtual ~index_t() {}
  virtual void clear() = 0;
  virtual size_t size()                const = 0;
  virtual size_t size(const size_t& d) const = 0;
  virtual idx_t& dereference(idx_t& _idx) const = 0;
};


/**
 * @brief Indexing in compressed sparse row
 * BASE: {0,1}-based indexing, other bases won't work
 */
template< int BASE >
struct index_compressed_sparse_row_t : index_t
{
  index_compressed_sparse_row_t() { clear(); }
  ~index_compressed_sparse_row_t() { clear(); }

  void initialize(const std::vector< std::vector< size_t > >& nz) {

#if 0
    nz.clear();
    nz.reserve(_nz.size());
    for (size_t r=0; r<_nz.size(); ++r)
      lss_index::vector_sorted_t::apply(nz.back(),r);
#endif
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

  // members
  std::vector< int > ia, ja;
  int nnu, nnz;
};


}  // namespace lss_index


#endif

