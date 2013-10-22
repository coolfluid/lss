// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_detail_index_hpp
#define cf3_lss_detail_index_hpp


#include <vector>
#include <algorithm>


namespace cf3 {
namespace lss {
namespace detail {


/* -- fundamental index pair (tuple?) type ---------------------------------- */


#if 0
/// @brief Basic index pair (tuple?) input/return type as (i[,j])
typedef struct ij {
  size_t idx[2];
  ij(const size_t& i=size_t(), const size_t& j=size_t()) { idx[0]=i; idx[1]=j; }
} ij;
#endif


/// @brief Basic index pair (tuple?) input/return type as (i[,j])
struct idx_t
{
//  size_t ij[2];  // TODO: check this option too?
  size_t i, j;
  idx_t(const size_t& _i=std::numeric_limits< size_t >::max(),
        const size_t& _j=std::numeric_limits< size_t >::max()) : i(_i), j(_j) {}

  bool operator<  (const idx_t& _other) const { return (i<_other.i? true : i>_other.i? false : (j<_other.j)); }
  bool operator>  (const idx_t& _other) const { return idx_t(_other)<*this; }
  bool operator== (const idx_t& _other) const { return i==_other.i && j==_other.j; }
  bool operator!= (const idx_t& _other) const { return i!=_other.i || j!=_other.j; }

  idx_t& invalidate() { return (*this = idx_t()); }
  bool is_valid_size()  const { return operator>(idx_t(0,0)) && operator<(idx_t()); }
  bool is_square_size() const { return i==j; }
  bool is_diagonal()    const { return is_square_size(); }
};


/* -- indexing techniques --------------------------------------------------- */


/// @brief Index dereferencing base class (pure abstract)
struct index_t
{
  virtual ~index_t() {}
  virtual void clear() = 0;
  virtual size_t size()                   const = 0;
  virtual size_t size(const size_t& d)    const = 0;
  virtual idx_t& dereference(idx_t& _idx) const = 0;
};


#if 0
namespace {
  static const std::string idx_name[] = { "",   "block_regular", "block_irregular", ""  };
  enum                     idx_type     { none,  block_regular,   block_irregular,  all };
  idx_type idx_name_to_type(const std::string& name) {
    for (int i=1; i<all; ++i)
      if (idx_name[i]==name)
        return (idx_type) i;
    return none;
  }
}  // namespace (unnamed)
#endif


#if 0
typename INDEX=detail::index_hierarchy_t< detail::index_hierarchy_t_end > >
#endif


/// @brief Matrix indexer assuming a irregular size block
struct index_irregular_block_t : index_t
{
  index_irregular_block_t() {}
  index_irregular_block_t(std::vector< size_t >&) {}

  // setup in construction
  index_irregular_block_t(
    const std::vector< size_t >& _block_size_i,
    const std::vector< size_t >& _block_size_j=std::vector< size_t >() )
  { setup(_block_size_i,_block_size_j); }

  // setup
  void setup(
    const std::vector< size_t >& _block_size_i,
    const std::vector< size_t >& _block_size_j=std::vector< size_t >() )
  {
    const std::vector< size_t >&
        sizei(_block_size_i),
        sizej(_block_size_j.size()? _block_size_j : _block_size_i);
  }

  // application
  idx_t& index(idx_t& idx) const { return idx; }
};


/// @brief Matrix indexer assuming a irregular size block
struct index_multi_domain_t : index_t
{
  index_multi_domain_t() {}
  index_multi_domain_t(std::vector< std::string >&) {}

  // setup in construction
  index_multi_domain_t(
    const std::vector< std::string >& _domain_names,
    const std::vector< size_t >& _domain_eqs_per_node)
  { setup(_domain_names,_domain_eqs_per_node); }

  // setup
  void setup(
    const std::vector< std::string >& _domain_names,
    const std::vector< size_t >& _domain_eqs_per_node )
  {
    std::vector< size_t > domsizes(_domain_names.size());
    // set sparsity pattern
    //const std::vector< std::vector< size_t > >&_nz
    // do it
  }

};


/// @brief Matrix indexer assuming a regular size block
struct index_regular_block_t : index_t
{
  index_regular_block_t() {}
  index_regular_block_t(std::vector< size_t >&) {}
  index_regular_block_t(size_t&) {}

  // setup in construction
  index_regular_block_t(
    size_t _block_size_i,
    size_t _block_size_j=size_t() )
  { setup(_block_size_i,_block_size_j); }

  // setup
  void setup(
    size_t _block_size_i,
    size_t _block_size_j=size_t() )
  {
    const size_t
      sizei(_block_size_i),
      sizej(_block_size_j? _block_size_j : _block_size_i);
  }

  // application
  idx_t& index(idx_t& idx) const { return idx; }
};


/**
 * @brief Matrix storage indexer
 *
 * This object wraps a method to return a single size_t index to a (very long
 * and boring) vector of entries that populate a sparse, dense, and/or
 * block-addressable matrix.
 *
 * Internally, the indexing is recursivelly nested to layer functionality as
 * desired. Addressing related to the structured storage is the bottom layer,
 * however block-addressing or multi-domain compositions can be layered over the
 * base storage indexing. Specialized indexing which combine multiple indexing
 * techniques in one layer can also be specified, at least that's the intent...
 */




/**
 * @brief Storage indexing in compressed sparse row
 * BASE: {0,1}-based indexing, other bases won't work
 */
template< int BASE >
struct index_compressed_sparse_row_t : index_t
{
  index_compressed_sparse_row_t() { clear(); }

  void setup(index_compressed_sparse_row_t& x) {
#if 0
    const std::vector< std::vector< size_t > >& nz;

    nz.clear();
    nz.reserve(_nz.size());
    for (size_t r=0; r<_nz.size(); ++r)
      index::vector_sorted_t::apply(nz.back(),r);

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
#endif
  }

  void swap(index_compressed_sparse_row_t& other) {
    other.ia.swap(ia);
    other.ia.swap(ja);
    std::swap(other.nnu,nnu);
    std::swap(other.nnz,nnz);
  }

  void clear() { ia.clear(); ja.clear(); nnu = nnz = 0; }

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


}  // namespace detail
}  // namespace lss
}  // namespace cf3


#endif

