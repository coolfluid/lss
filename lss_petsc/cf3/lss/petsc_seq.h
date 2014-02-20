// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_petsc_seq_hpp
#define cf3_lss_petsc_seq_hpp


#include "petscksp.h"
#include "petscpc.h"

#include "LibLSS_PETSC.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Interface to PETSc linear system solver, sequential (serial) version
 * @author Pedro Maciel
 */
class lss_API petsc_seq : public linearsystem< double >
{
  // utility definitions
//typedef petsc::matrix_wrapper matrix_t;
  typedef sparse_matrix< double, sort_by_row, 0 > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "petsc_seq"; }

  /// Construction
  petsc_seq(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1 );

  /// Destruction
  ~petsc_seq();

  /// Linear system solving
  petsc_seq& solve();

  /// Linear system copy
  petsc_seq& copy(const petsc_seq& _other);

  /// Linear system forward multiplication
  petsc_seq& multi(const double& _alpha=1., const double& _beta=0.);


 private:
  // internal functions

  /// Verbose error message
  static const std::string err_message(const int& err, const char* basemsg);

  /// Verbose converged/diverged message
  static const std::string converged_message(const KSPConvergedReason& rsn);


 protected:
  // linear system matrix interfacing

  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) { m_A.initialize(i,j,_nnz); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___assign(const double& _value) { m_A = _value;   }
  void A___clear()                      { m_A.clear();    }
  void A___zerorow(const size_t& i)     { m_A.zerorow(i); }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc); }

  void   A___print(std::ostream& o, const print_t& l=print_auto) const { m_A.print(o,l); }
  size_t A___size(const size_t& d)  const { return m_A.size(d);  }


 protected:
  // storage

  matrix_t m_A;  // sparse system matrix
  KSP      ksp;  // (PETSc) Krylov method context
  PC       pc;   // (PETSc) preconditioner context
  struct {
    std::string ksptype;
    std::string pctype;
    bool        monitor;
    PetscInt    restart;
    PetscInt    ovl;
    PetscInt    maxits;
    PetscReal   rtol;
    PetscReal   abstol;
    PetscReal   dtol;
  } opt;    // (PETSc) options

};


}  // namespace lss
}  // namespace cf3


#endif

