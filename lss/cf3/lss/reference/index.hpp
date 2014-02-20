// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_index_hpp
#define cf3_lss_index_hpp


#include <limits>


namespace cf3 {
namespace lss {


/* -- indexing techniques --------------------------------------------------- */


/// @brief Index dereferencing base class (pure abstract)
struct index_t
{
  virtual ~index_t() {}
  virtual size_t size(const size_t& d)    const = 0;
  virtual idx_t& dereference(idx_t& _idx) const = 0;

  const idx_t& size() const { return m_size; }
  idx_t m_size;
};


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


/// @brief Hierarchical indexing: type list termination type
struct index_hierarchy_t_end
{
  idx_t& dereference(idx_t& _idx) const { return _idx; }
  void initialize() {}
};


/// @brief Hierarchical indexing: nested hierarchy definition
template<
    typename Idx,
    typename IdxNested=index_hierarchy_t_end >
struct index_hierarchy_t
{
  idx_t& dereference(idx_t& _idx) const { return IdxNested::dereference(Idx::dereference(_idx)); }
  void initialize() { Idx::initialize(); IdxNested::initialize(); }
};


}  // namespace lss
}  // namespace cf3


#endif
