// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_lss_index_h
#define cf3_lss_lss_index_h


#include <vector>


namespace cf3 {
namespace lss {
namespace lss_index {


// helper forward declaration to generic index base type
struct index;


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
 * @brief Instance of the handler of "index creators" registered against a key
 * (a string), which can then create "index" objects using such key.
 */
struct index_factory_t {

  // indexing creation registration
  void Register(const index_creator_t* creator) {
    if (creator!=NULL)
      if (creator->key.length() && std::find(m_icreators.begin(),m_icreators.end(),creator->key)==m_icreators.end()) {
        m_icreators.push_back(creator);
        return;
      }
    std::cerr << "index_factory_t: warning: cannot register key \"" << key << "\"!" << std::endl;
  }

  // indexing creation from key
  index* Create(const std::string& key) {
    if (std::find(m_icreators.begin(),m_icreators.end(),key)==m_icreators.end()) {
      std::cerr << "index_factory_t: warning: creation of key \"" << key << "\" not possible!" << std::endl;
      return NULL;
    }
    return std::find(m_icreators.begin(),m_icreators.end(),key)->create();
  }

  // destructor
  ~index_factory_t() {
    for (size_t i=0; i<m_icreators.size(); ++i)
      delete m_icreators[i];
  }

  // registered indexing creators
  std::vector< index_creator_t* > m_icreators;

}
index_factory;


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
#if 0
virtual void resize(const std::vector< std::vector< unsigned > >& nz) = 0;
#endif
};




}  // namespace lss_index
}  // namespace lss
}  // namespace cf3


#endif

