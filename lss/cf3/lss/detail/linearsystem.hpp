// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_detail_linearsystem_hpp
#define cf3_lss_detail_linearsystem_hpp


#include <limits>
#include <sstream>
#include <stdexcept>

#include "utilities.hpp"
#include "index.hpp"
#include "matrix.hpp"


namespace cf3 {
namespace lss {
namespace detail {


/* -- linear system --------------------------------------------------------- */

template<
    typename T,
    typename MATRIX=detail::dense_matrix_v< T >,
    typename VECTOR=detail::dense_matrix_v< T > >
class linearsystem
{

  // -- Construction and destruction
 public:

  /// Construct the linear system
  linearsystem(
      const size_t& _size_i=size_t(),
      const size_t& _size_j=size_t(),
      const size_t& _size_k=1,
      const T& _value=T() ) { resize(_size_i,_size_j,_size_k,_value); }

  /// Destruct the linear system
  virtual ~linearsystem() {}


  // -- Interfacing
 public:

  /// Resize the linear system (consistently)
  linearsystem& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const double& _value=T()) {
    consistent(_size_i,_size_j,_size_i,_size_k,_size_j,_size_k);
    const T value(static_cast< T >(_value));
    A().resize(_size_i,_size_j,value);
    b().resize(_size_i,_size_k,value);
    x().resize(_size_j,_size_k,value);
    return *this;
  }

  /// Initialize linear system from file(s)
  linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) A().initialize(_Afname);
    if (_bfname.length()) b().initialize(_bfname); else b().resize(size(0),1);
    if (_xfname.length()) x().initialize(_xfname); else x().resize(size(1),size(2));
    consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Initialize linear system from vectors of values (lists, in the right context)
  linearsystem& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< T >(),
      const std::vector< double >& vx=std::vector< T >()) {

    // input parameter type insulation from the templated base class
    const bool conversion_needed(typeid(T)!=typeid(double));
    std::vector< T >
      another_A(vA.size()),
      another_b(vb.size()),
      another_x(vx.size());
    if (conversion_needed) {
      std::transform( vA.begin(),vA.end(),another_A.begin(),storage_conversion_t< double, T >() );
      std::transform( vb.begin(),vb.end(),another_b.begin(),storage_conversion_t< double, T >() );
      std::transform( vx.begin(),vx.end(),another_x.begin(),storage_conversion_t< double, T >() );
    }

    // propper initialization of linear system components
    std::vector< T >&
      correct_A(conversion_needed? another_A : (std::vector< T >&) vA),
      correct_b(conversion_needed? another_b : (std::vector< T >&) vb),
      correct_x(conversion_needed? another_x : (std::vector< T >&) vx);
    if (vA.size()) A().initialize(correct_A);
    if (vb.size()) b().initialize(correct_b); else b().resize(size(0),1);
    if (vx.size()) x().initialize(correct_x); else x().resize(size(1),size(2));

    consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Value assignment (operator)
  linearsystem& operator=(const T& _value) { return resize(size(0),size(1),size(2),_value); }

  /// Value assignment (method)
  linearsystem& assign(const T& _value=T()) { return operator=(_value); }

  /// Clear the contents
  linearsystem& clear() {
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


  /// Checks whether the matrix/vectors sizes are consistent in the system
  bool consistent(const size_t& Ai, const size_t& Aj,
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

 private:

  /// Output to given std::ostream
  template< typename aT, typename aMATRIX, typename aVECTOR >
  friend std::ostream& operator<< (std::ostream&, const linearsystem< aT, aMATRIX, aVECTOR >&);
  // (a private generic friend? how promiscuous!)


  // -- Pure methods
 public:

  /// Linear system solving and components access
  virtual linearsystem& solve() = 0;
  virtual       MATRIX& A()       = 0;
  virtual       VECTOR& b()       = 0;
  virtual       VECTOR& x()       = 0;
  virtual const MATRIX& A() const = 0;
  virtual const VECTOR& b() const = 0;
  virtual const VECTOR& x() const = 0;

};


/// Output to given std::ostream (non-member version)
template< typename T, typename MATRIX, typename VECTOR >
std::ostream& operator<< (std::ostream& o, const linearsystem< T, MATRIX, VECTOR >& lss)
{
  o << "linearsystem: A: "; lss.A_print(o); o << std::endl;
  o << "linearsystem: b: "; lss.b_print(o); o << std::endl;
  o << "linearsystem: x: "; lss.x_print(o); o << std::endl;
  return o;
}


}  // namespace detail
}  // namespace lss
}  // namespace cf3


#endif
