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

#include "lss_utilities.hpp"
#include "lss_index.hpp"
#include "lss_matrix.hpp"
#include "lss_linearsystem.hpp"


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */

template<
    typename T,
    typename INDEX=index_hierarchy_t< index_hierarchy_t_end > >
class linearsystem :
  public detail::lss_linearsystem<
    T,
    detail::dense_matrix_v< T >,
    detail::dense_matrix_v< T > >,
  public common::Action
{ public:  //FIXME update permissions

  // utility definitions
  typedef INDEX index_t;

  // construction, destruction, initialization and execution

  /// Construct the linear system, by direct initialization
  linearsystem(const std::string& name) : common::Action(name) {}

  /// Destruct the linear system
  virtual ~linearsystem() {}

  /// Initialize linear system from file(s)
  linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    return this->initialize(_Afname,_bfname,_xfname);
  }

  /// Initialize linear system from vectors of values (lists, in right context)
  linearsystem& initialize(
      const std::vector< T >& vA,
      const std::vector< T >& vb=std::vector< T >(),
      const std::vector< T >& vx=std::vector< T >()) {
    return this->initialize(vA,vb,vx);
  }

  /// Linear system solve aliased to action execute
  void execute() { this->solve(); }


  // interfacing

  /// Linear system value assignment operator
  linearsystem& operator=(const T& _value) { return this->operator =(_value); }

  /// Linear system value assignment method
  linearsystem& assign(const T& _value=T()) { return this->operator=(_value); }

  /// Changes the number of elements stored
  linearsystem& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const T& _value=T()) { return this->resize(_size_i,_size_j,_size_k); }

  /// Clears the contents
  linearsystem& clear() { return this->clear(); }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const { return this->size(d); }

  /// Checks whether the linear system matrix is empty
  bool empty() { return this->empty(); }

  // -- Utilities

  // a private friend, how promiscuous
 private:
  template< typename aT, typename aINDEX >
  friend std::ostream& operator<< (std::ostream&, const linearsystem< aT, aINDEX >&);

  // storage
 protected:
  index_t m_idx;  // (only indexing is part of the base type)
};


/* -- linear system output -------------------------------------------------- */

template< typename T, typename INDEX >
std::ostream& operator<< (std::ostream& o, const linearsystem< T, INDEX >& lss)
{
  o << "linearsystem: A: "; lss.A_print(o); o << std::endl;
  o << "linearsystem: b: "; lss.b_print(o); o << std::endl;
  o << "linearsystem: x: "; lss.x_print(o); o << std::endl;
  return o;//.print(o);
}


}  // namespace lss
}  // namespace cf3


#endif
