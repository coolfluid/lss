// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_lss_index_h
#define cf3_lss_lss_index_h


#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "boost/parameter.hpp"


namespace cf3 {
namespace lss {
namespace lss_index {


// matrix index return and input type as (i[,j]) pair
typedef struct ij {
  size_t idx[2];
  ij(const size_t& i=size_t(), const size_t& j=size_t()) { idx[0]=i; idx[1]=j; }
} ij;






struct index_t
{
  void index(ij& _idx) {}
};


struct index_multi_domain_t : index_t
{
  index_multi_domain_t() {}
  index_multi_domain_t(std::vector< std::string >&) {}
};


struct index_irregular_block_t : index_t
{
  index_irregular_block_t() {}
  index_irregular_block_t(std::vector< size_t >&) {}
};


struct index_regular_block_t : index_t
{
  index_regular_block_t() {}
  index_regular_block_t(std::vector< size_t >&) {}
  index_regular_block_t(size_t&) {}
};










#if 0
/**
 * @brief Matrix indexer
 */
class index_t
{
  template< typename T, typename MATRIX >
  friend class linearsystem;

// (deliberately force right methods on the right types at runtime)
#define NO_IMPLEMENTATION throw std::runtime_error("index_t: wrong method call")

#if 0
  index_t() : type(detail::none) {}
  index_t(const std::string& _type) : type(detail::idx_name_to_type(_type)) {}
  const detail::idx_type type;
#endif

 public:
  // index configuration methods
  virtual void index(const size_t& _i, const size_t& _j) { NO_IMPLEMENTATION; }

 private:
  // index application
  virtual ij& index(ij& idx) const { NO_IMPLEMENTATION; return idx; }

#undef NO_IMPLEMENTATION
};


/**
 * @brief Matrix indexer assuming a regular size block
 */
struct index_regular_block_t : public index_t
{
  // setup in construction
  index_regular_block_t(
    size_t _block_size_i,
    size_t _block_size_j=size_t() )
  {
    const size_t
      sizei(_block_size_i),
      sizej(_block_size_j? _block_size_j : _block_size_i);
    // do it
  }

  // configuration

  // application
  ij& index(ij& idx) const {
    return idx;
  }
};


/**
 * @brief Matrix indexer assuming a irregular size block
 */
struct index_irregular_block_t : public index_t
{
  // no-setup construction (for internal handling)
  index_irregular_block_t() {}

  // setup in construction
  index_irregular_block_t(
    const std::vector< size_t >& _block_size_i,
    const std::vector< size_t >& _block_size_j=std::vector< size_t >() )
  { setup_index_irregular_block(_block_size_i,_block_size_j); }

  // setup
  void setup_index_irregular_block(
    const std::vector< size_t >& _block_size_i,
    const std::vector< size_t >& _block_size_j=std::vector< size_t >() ) {
    const std::vector< size_t >&
        sizei(_block_size_i),
        sizej(_block_size_j.size()? _block_size_j : _block_size_i);
  }

  // configuration

  // application
  ij& index(ij& idx) const {
    return idx;
  }
};


/**
 * @brief Matrix indexer assuming a regular size block
 */
struct index_multi_domain_t : index_irregular_block_t
{
  // setup in construction
  index_multi_domain_t(
    const std::vector< std::string >& _domain_names,
    const std::vector< size_t >& _domain_eqs_per_node)
  {
    std::vector< size_t > domsizes(_domain_names.size());

    // set sparsity pattern
    // const std::vector< std::vector< size_t > >&_nz

    // do it
    // return index_irregular_block_t::setup_block_irregular_size(domsizes);
  }
};
#endif


#if 0
namespace detail {
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
}  // namespace detail
#endif


#if 0
  // indexing setup
  BOOST_PARAMETER_NAME(sparsity_structure_per_row)
  BOOST_PARAMETER_NAME(block_regular_size_i)
  BOOST_PARAMETER_NAME(block_regular_size_j)
  BOOST_PARAMETER_NAME(block_irregular_size_i)
  BOOST_PARAMETER_NAME(block_irregular_size_j)
  BOOST_PARAMETER_NAME(domains)

  // indexing configuration
  BOOST_PARAMETER_NAME(block_i)
  BOOST_PARAMETER_NAME(block_j)
  BOOST_PARAMETER_NAME(domain)

  // indexing application
  BOOST_PARAMETER_NAME(idx)


  BOOST_PARAMETER_MEMBER_FUNCTION( (index_t&), setup, tag, (optional
    (sparsity_structure_per_row, *, std::vector< std::vector< size_t > >())
    (block_regular_size_i, (size_t), 0)
    (block_regular_size_j, (size_t), 0)
    (block_irregular_size_i, (std::vector< size_t >), 0)
    (block_irregular_size_j, (std::vector< size_t >), 0)
                                                            ))
  {
    switch (type) {

      case block_regular : {
        idx[0] += block_index[0]*block_size[0];
        idx[1] += block_index[1]*block_size[1];
        break;
      }

      case block_irregular : { break; }

      case none: case all: default: {}
    }
    return *this;
  }

  BOOST_PARAMETER_CONST_MEMBER_FUNCTION( (idx_t&), index, tag, (required
    (in_out(idx), (idx_t&)) ))
  {
    switch (type) {
      case block_regular : {
        idx[0] += block_index[0]*block_size[0];
        idx[1] += block_index[1]*block_size[1];
        break;
      }
      case block_irregular : { break; }
      case none:
      case all:
      default: {}
    }
    return idx;
  }

  BOOST_PARAMETER_CONST_MEMBER_FUNCTION( (idx_t&), index, tag, (required
    (in_out(idx), (idx_t&)) ))
  {
    switch (type) {
      case block_regular :   { break; }
      case block_irregular : { break; }
      case none:
      case all:
      default: {}
    }
    return idx;
  }
#endif


#if 0
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
 * @brief Auxiliary type to generic "index creator"
 */
struct index_creator_t {
  index_creator_t(const std::string& _key) : key(_key) {}
  virtual ~index_creator_t() {}
  virtual index_t* create() = 0;
  const std::string key;
};


/**
 * @brief Auxiliary type to specific "index creator"
 */
template< class INDEX >
struct index_creator : index_creator_t {
  index_creator(const std::string& _key) : index_creator_t(_key) {}
  index_t* create() { return new INDEX(); }
};


/**
 * @brief Handler of "index creators" registered against std::string keys, which
 * are used to create specific "index" objects (factory pattern).
 * The matrix addressing factory is declared below and instantiated outside.
 */
class index_factory_t
{
 private:
  // registered indexing creators
  typedef std::vector< index_creator_t* > icreators_t;
  icreators_t m_creators;

  // helper class for checking if specific creator keys are already registered
  struct creator_has_key
  {
    creator_has_key(const std::string& _key) : key(_key) {}
    bool operator() (index_creator_t*& c) const { return key==c->key; }
    const std::string key;
  };

 public:
  // indexing creation registration
  void Register(index_creator_t* creator) {
    if ( creator->key.length() && m_creators.end()==std::find_if(
          m_creators.begin(),m_creators.end(),creator_has_key(creator->key)) ) {
      m_creators.push_back(creator);
      return;
    }
    std::cerr << "index_factory_t: unable to register key \"" << creator->key << "\"." << std::endl;
  }

  // indexing creation from key
  index_t* Create(const std::string& key) {
    icreators_t::iterator c(std::find_if(m_creators.begin(),m_creators.end(),creator_has_key(key)));
    if (c==m_creators.end()) {
      std::cerr << "index_factory_t: unable to create from key \"" << key << "\"." << std::endl;
      return NULL;
    }
    return (*c)->create();
  }

  void output(std::ostream& out) {
    out << "index_factory_t: " << m_creators.size() << " key(s) registered.";
    BOOST_FOREACH(index_creator_t*& c, m_creators) {
      out << "\n  - \"" << c->key << '"';
    }
    out << std::endl;
  }

  // constructor registering (a series of) indexing techniques
  index_factory_t(const size_t n, index_creator_t* c, ...);

  // destructor to handle deletion of indexing techniques creators
  ~index_factory_t() {
    icreators_t::reverse_iterator c;
    for (c=m_creators.rbegin(); c!=m_creators.rend(); ++c)
      delete (*c);
  }

}
extern idx_factory_matrix_indexing;
#endif


}  // namespace lss_index
}  // namespace lss
}  // namespace cf3


#endif

