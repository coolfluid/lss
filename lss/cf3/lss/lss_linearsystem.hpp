// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_lss_linearsystem_hpp
#define cf3_lss_lss_linearsystem_hpp


#include <limits>
#include <sstream>
#include <stdexcept>

#include "lss_utilities.hpp"
#include "lss_index.hpp"
#include "lss_matrix.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- linear system --------------------------------------------------------- */

template<
    typename T,
    typename MATRIX=detail::dense_matrix_v< T >,
    typename VECTOR=detail::dense_matrix_v< T > >
class lss_linearsystem
{ public:  //FIXME update permissions

  // construction, destruction and initialization

  /// Construct the linear system, by direct initialization
  lss_linearsystem(const std::string& name) : common::Action(name) {}
  lss_linearsystem(
      const size_t& _size_i=size_t(),
      const size_t& _size_j=size_t(),
      const size_t& _size_k=1,
      const T& _value=T() ) { resize(_size_i,_size_j,_size_k,_value); }

  /// Destruct the linear system
  virtual ~lss_linearsystem() {}

  /// Initialize linear system from file(s)
  lss_linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) A().initialize(_Afname);
    if (_bfname.length()) b().initialize(_bfname); else b().resize(size(0),1);
    if (_xfname.length()) x().initialize(_xfname); else x().resize(size(1),size(2));
    is_size_consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Initialize linear system from vectors of values (lists, in right context)
  lss_linearsystem& initialize(
      const std::vector< T >& vA,
      const std::vector< T >& vb=std::vector< T >(),
      const std::vector< T >& vx=std::vector< T >()) {
    if (vA.size()) A().initialize(vA);
    if (vb.size()) b().initialize(vb); else b().resize(size(0),1);
    if (vx.size()) x().initialize(vx); else x().resize(size(1),size(2));
    is_size_consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Linear system solve aliased to action execute
  void execute() { solve(); }


  // (pure method)
  /// Linear system solve (what everyone is waiting for!)
  virtual lss_linearsystem& solve() = 0;


  // interfacing

  /// Linear system value assignment operator
  lss_linearsystem& operator=(const T& _value) { return resize(size(0),size(1),size(2),_value); }

  /// Linear system value assignment method
  lss_linearsystem& assign(const T& _value=T()) { return operator=(_value); }

  /// Changes the number of elements stored
  lss_linearsystem& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const T& _value=T()) {
    is_size_consistent(_size_i,_size_j,_size_i,_size_k,_size_j,_size_k);
    A().resize(_size_i,_size_j,_value);
    b().resize(_size_i,_size_k,_value);
    x().resize(_size_j,_size_k,_value);
    return *this;
  }

  /// Clears the contents
  lss_linearsystem& clear() {
    A().clear();
    b().clear();
    x().clear();
    return *this;
  }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const {
    return (d< 2? A().size(d) :
           (d==2? b().size(1) : std::numeric_limits< size_t >::max()));
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

  // (pure methods)
  /// Linear system components access
 private:
  virtual       MATRIX& A()       = 0;
  virtual       VECTOR& b()       = 0;
  virtual       VECTOR& x()       = 0;
  virtual const MATRIX& A() const = 0;
  virtual const VECTOR& b() const = 0;
  virtual const VECTOR& x() const = 0;

  // a private friend, how promiscuous
 private:
  template< typename aT, typename aMATRIX, typename aVECTOR >
  friend std::ostream& operator<< (std::ostream&, const lss_linearsystem< aT, aMATRIX, aVECTOR >&);

};


/* -- linear system output -------------------------------------------------- */

template< typename T, typename MATRIX, typename VECTOR >
std::ostream& operator<< (std::ostream& o, const lss_linearsystem< T, MATRIX, VECTOR >& lss)
{
  o << "linearsystem: A: "; lss.A_print(o); o << std::endl;
  o << "linearsystem: b: "; lss.b_print(o); o << std::endl;
  o << "linearsystem: x: "; lss.x_print(o); o << std::endl;
  return o;//.print(o);
}


}  // namespace detail
}  // namespace lss
}  // namespace cf3


#endif
