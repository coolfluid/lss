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
#include "boost/foreach.hpp"


namespace cf3 {
namespace lss {
namespace lss_index {


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
struct index {

  // destructor
  virtual ~index() {}

  // setup
  virtual index& setup(
    const std::vector< std::vector< size_t > >& _nz,
    const size_t _regular_block_size) = 0;
};


/**
 * @brief Matrix indexer assuming regular number of system equations
 */
class idx_regular_block : public index
{
 public:
  index& setup(
     const std::vector< std::vector< size_t > >& _nz,
     const size_t _regular_block_size) {
    return *this;
  }
};


/**
 * @brief Matrix storage in compressed sparse row format (CSR)
 */
class idx_sparse_matrix_csr : public index
{
 public:
  index& setup(
     const std::vector< std::vector< size_t > >& _nz,
     const size_t _regular_block_size) {
    return *this;
  }
};


/**
 * @brief Auxiliary type to generic "index creator"
 */
struct index_creator_t {
  index_creator_t(const std::string& _key) : key(_key) {}
  virtual ~index_creator_t() {}
  virtual index* create() = 0;
  const std::string key;
};


/**
 * @brief Auxiliary type to specific "index creator"
 */
template< class INDEX >
struct index_creator : index_creator_t {
  index_creator(const std::string& _key) : index_creator_t(_key) {}
  index* create() { return new INDEX(); }
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
  index* Create(const std::string& key) {
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


}  // namespace lss_index
}  // namespace lss
}  // namespace cf3


#endif

