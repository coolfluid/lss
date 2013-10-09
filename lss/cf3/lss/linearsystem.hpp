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


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */

template<
    typename T,
    typename INDEX=index_hierarchy_t< index_hierarchy_t_end > >
class linearsystem :
  public common::Action
{ public:  //FIXME update permissions

  // utility definitions
  typedef INDEX index_t;

  // construction, destruction and initialization

  /// Construct the linear system, by direct initialization
  linearsystem(const std::string& name) : common::Action(name) {}
  linearsystem(
      const size_t& _size_i=size_t(),
      const size_t& _size_j=size_t(),
      const size_t& _size_k=1,
      const T& _value=T() ) { resize(_size_i,_size_j,_size_k,_value); }

  /// Construct the linear system, by copy
  linearsystem(const linearsystem& _other) { operator=(_other); }

  /// Destructs the linear system
  virtual ~linearsystem() {}

  /// Initialize linear system from file(s)
  linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) A_initialize(_Afname);
    if (_bfname.length()) b_initialize(_bfname); else b_resize(size(0),1);
    if (_xfname.length()) x_initialize(_xfname); else x_resize(size(1),size(2));
    is_size_consistent(A_size(0),A_size(1),b_size(0),b_size(1),x_size(0),x_size(1));
    return *this;
  }

  /// Initialize linear system from vectors of values (lists, in right context)
  linearsystem& initialize(
      const std::vector< T >& vA,
      const std::vector< T >& vb=std::vector< T >(),
      const std::vector< T >& vx=std::vector< T >()) {
    if (vA.size()) A_initialize(vA);
    if (vb.size()) b_initialize(vb); else b_resize(size(0),1);
    if (vx.size()) x_initialize(vx); else x_resize(size(1),size(2));
    is_size_consistent(A_size(0),A_size(1),b_size(0),b_size(1),x_size(0),x_size(1));
    return *this;
  }

  /// Execute redirects to solve
  void execute () { solve(); }


  // (pure methods)
  /// Accessors/mutators
  virtual       T& A(const size_t& i, const size_t& j)         = 0;
  virtual       T& b(const size_t& i, const size_t& j=0)       = 0;
  virtual       T& x(const size_t& i, const size_t& j=0)       = 0;
  virtual const T& A(const size_t& i, const size_t& j)   const = 0;
  virtual const T& b(const size_t& i, const size_t& j=0) const = 0;
  virtual const T& x(const size_t& i, const size_t& j=0) const = 0;


  // (pure methods)
  /// Matrix/vectors resizing
  virtual linearsystem& A_resize(const size_t& _size_i, const size_t& _size_j,   const T& _value=T()) = 0;
  virtual linearsystem& b_resize(const size_t& _size_i, const size_t& _size_j=1, const T& _value=T()) = 0;
  virtual linearsystem& x_resize(const size_t& _size_i, const size_t& _size_j=1, const T& _value=T()) = 0;


  // (pure methods)
  /// Matrix/vectors initialization from given vector or filename
  virtual linearsystem& A_initialize(const std::vector< T >&) = 0;
  virtual linearsystem& b_initialize(const std::vector< T >&) = 0;
  virtual linearsystem& x_initialize(const std::vector< T >&) = 0;
  virtual linearsystem& A_initialize(const std::string&)      = 0;
  virtual linearsystem& b_initialize(const std::string&)      = 0;
  virtual linearsystem& x_initialize(const std::string&)      = 0;


  // (pure methods)
  /// Matrix/vectors clearing
  virtual void A_clear() = 0;
  virtual void b_clear() = 0;
  virtual void x_clear() = 0;


  // (pure methods)
  /// Matrix/vectors output
  virtual void A_print(std::ostream&) = 0;
  virtual void b_print(std::ostream&) = 0;
  virtual void x_print(std::ostream&) = 0;


  // (pure methods)
  /// Matrix/vectors size
  virtual size_t A_size(const size_t&) const = 0;
  virtual size_t b_size(const size_t&) const = 0;
  virtual size_t x_size(const size_t&) const = 0;


/*
  // (pure method)
  /// Linear system copy assignment operator
  virtual linearsystem& operator=(const linearsystem& _other) = 0;
*/


  // (pure method)
  /// Linear system solve (what everyone is waiting for!)
  virtual linearsystem& solve() = 0;

  // interfacing

  /// Linear system value assignment operator
  linearsystem& operator=(const T& _value) { return resize(size(0),size(1),size(2),_value); }

  /// Linear system value assignment method
  linearsystem& assign(const T& _value=T()) { return operator=(_value); }

  /// Changes the number of elements stored
  linearsystem& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const T& _value=T()) {
    is_size_consistent(_size_i,_size_j,_size_i,_size_k,_size_j,_size_k);
    A_resize(_size_i,_size_j,_value);
    b_resize(_size_i,_size_k,_value);
    x_resize(_size_j,_size_k,_value);
    return *this;
  }

  /// Clears the contents
  linearsystem& clear() {
    A_clear();
    b_clear();
    x_clear();
    return *this;
  }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const {
    return (d< 2? A_size(d) :
           (d==2? b_size(1) : std::numeric_limits< size_t >::max()));
  }

  /// Checks whether the linear system matrix is empty
  bool empty() { return !(size(0)*size(1)*size(2)); }

  // -- Utilities

 protected:
  bool is_size_consistent(const size_t& Ai, const size_t& Aj,
                          const size_t& bi, const size_t& bj,
                          const size_t& xi, const size_t& xj) const {
    if (Ai!=bi || Aj!=xi || bj!=xj || !(Ai*Aj*bi*bj*xi*xj)) {
      std::ostringstream msg;
      msg << "linearsystem: size is not consistent: "
          << "A(" << Ai << 'x' << Aj << ") "
          << "x(" << xi << 'x' << xj << ") = "
          << "b(" << bi << 'x' << bj << ") ";
      throw std::runtime_error(msg.str());
    }
    return true;
  }

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
