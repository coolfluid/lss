// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_hpp
#define cf3_lss_linearsystem_hpp


#include <limits>
#include <sstream>
#include <stdexcept>

#include "common/Action.hpp"

#include "detail/index.hpp"
#include "detail/matrix.hpp"


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */
/* note: the class name forces the distinction of (this) acessible plugin     */
/* and the detail::linearsystem building class (to avoid 'using namespace's)  */

template<
    typename T,
    typename INDEX=detail::index_hierarchy_t< detail::index_hierarchy_t_end > >
class linearsystem :
  public common::Action
{
  typedef INDEX index_t;

  // -- Construction and destruction
 public:

  /// Construct the linear system
  linearsystem(const std::string& name) : common::Action(name) {}

  /// Destruct the linear system
  virtual ~linearsystem() {}


  // -- Interfacing (pure)
 public:

  /// Linear system resizing (consistently)
  virtual linearsystem& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const T& _value=T()) = 0;

  /// Linear system initialization from file(s)
  virtual linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) = 0;

  /// Linear system initialization from vectors of values (lists, in right context)
  virtual linearsystem& initialize(
      const std::vector< T >& vA,
      const std::vector< T >& vb=std::vector< T >(),
      const std::vector< T >& vx=std::vector< T >()) = 0;

  /// Linear system solving
  virtual linearsystem& solve() = 0;


#if 0
  // -- Interfacing
 public:

  /// Value assignment (operator)
  linearsystem& operator=(const T& _value) { return this->operator =(_value); }

  /// Value assignment (method)
  linearsystem& assign(const T& _value=T()) { return this->operator=(_value); }

  /// Clear the contents
  linearsystem& clear() { return this->clear(); }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const { return this->size(d); }

  /// Checks whether the linear system matrix is empty
  bool empty() { return this->empty(); }
#endif

  /// Linear system solving, aliased from execute
  void execute() { solve(); }


  // -- Storage
 private:
  index_t m_idx;

};


}  // namespace lss
}  // namespace cf3


#endif
