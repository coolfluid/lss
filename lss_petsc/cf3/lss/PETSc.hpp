// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_PETSc_hpp
#define cf3_lss_PETSc_hpp


#include "LibLSS_PETSC.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


namespace detail {


/// @brief PETSc matrix wrapper (consistent with cf3::lss::detail::matrix<...>)
struct petsc_matrix_wrapper  :
  matrix< double, petsc_matrix_wrapper >
{
  petsc_matrix_wrapper() : v(double()) {}

  petsc_matrix_wrapper& initialize(const size_t& i, const size_t& j, const double& _value=double()) { return *this; }
  petsc_matrix_wrapper& initialize(const std::vector< double >& _vector) { return *this; }
  petsc_matrix_wrapper& initialize(const std::string& _fname)            { return *this; }

  petsc_matrix_wrapper& clear()                  { return *this; }
  petsc_matrix_wrapper& zerorow(const size_t& i) { return *this; }

  petsc_matrix_wrapper& operator=(const petsc_matrix_wrapper& _other) { return *this; }
  petsc_matrix_wrapper& operator=(const double& _value)               { return *this; }

  const double& operator()(const size_t& i, const size_t& j=0) const { return v; }
        double& operator()(const size_t& i, const size_t& j=0)       { return v; }

  double v;
};


/// @brief PETSc vector wrapper (consistent with cf3::lss::detail::matrix<...>)
struct petsc_vector_wrapper  :
  matrix< double, petsc_vector_wrapper >
{
  petsc_vector_wrapper() : v(double()) {}

  petsc_vector_wrapper& initialize(const size_t& i, const size_t& j, const double& _value=double()) { return *this; }
  petsc_vector_wrapper& initialize(const std::vector< double >& _vector) { return *this; }
  petsc_vector_wrapper& initialize(const std::string& _fname)            { return *this; }

  petsc_vector_wrapper& clear()                  { return *this; }
  petsc_vector_wrapper& zerorow(const size_t& i) { return *this; }

  petsc_vector_wrapper& operator=(const petsc_vector_wrapper& _other) { return *this; }
  petsc_vector_wrapper& operator=(const double& _value)               { return *this; }

  const double& operator()(const size_t& i, const size_t& j=0) const { return v; }
        double& operator()(const size_t& i, const size_t& j=0)       { return v; }

  double v;
};


}  // namespace detail


/**
 * @brief Interface to PETSc linear system solver.
 * @author Pedro Maciel
 */
class lss_API PETSc : public
  linearsystem< double,
    detail::petsc_matrix_wrapper,
    detail::petsc_vector_wrapper >
{
  // utility definitions
  typedef detail::petsc_matrix_wrapper matrix_t;
  typedef detail::petsc_vector_wrapper vector_t;
  typedef linearsystem< double, matrix_t, vector_t > linearsystem_t;


 public:
  // framework interfacing
  static std::string type_name() { return "petsc"; }


  /// Construction
  PETSc(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1,
    const double& _value=double() );

  /// Solve
  PETSc& solve();


  // linear system components access
 public:
        matrix_t& A()       { return m_A; }
        vector_t& b()       { return m_b; }
        vector_t& x()       { return m_x; }
  const matrix_t& A() const { return m_A; }
  const vector_t& b() const { return m_b; }
  const vector_t& x() const { return m_x; }


  // members
 protected:
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;

};


}  // namespace lss
}  // namespace cf3


#endif

