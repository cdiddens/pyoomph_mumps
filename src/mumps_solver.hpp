// A thin, arithmetic-generic C++ wrapper around the MUMPS C interface.
//
// Deliberately free of nanobind: everything here works on raw pointers and std::vector, so the same
// class can be driven from a C++ test as well as from Python, and so that the binding layer in
// bindings.cpp has exactly one job (converting numpy arrays and releasing the GIL).
//
// The one class is a template over the scalar, instantiated in mumps_solver.cpp for double (dmumps)
// and std::complex<double> (zmumps). Both MUMPS structs are identical apart from the scalar type, so
// a traits struct is all that is needed to share the implementation - which matters more than it
// looks: the complex arithmetic is not an exotic extra here but the ordinary case for the
// eigensolver, where J - sigma*M is complex for every complex shift even when J and M are real.

#pragma once

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <dmumps_c.h>
#include <zmumps_c.h>

namespace pyoomph_mumps
{

  /// Raised for anything MUMPS reports as a negative INFOG(1), and for the misuse this wrapper
  /// detects itself. The binding layer maps it onto a Python exception that pyoomph's solver turns
  /// into a retryable SolverError, so a Newton or arclength step shrinks and tries again rather than
  /// ending the run.
  class MumpsError : public std::runtime_error
  {
  public:
    MumpsError(const std::string &what, int infog1 = 0, int infog2 = 0)
        : std::runtime_error(what), infog1(infog1), infog2(infog2) {}
    int infog1;
    int infog2;
  };

  template <typename Scalar> struct MumpsTraits;

  template <> struct MumpsTraits<double>
  {
    using Struct = DMUMPS_STRUC_C;
    using MumpsScalar = DMUMPS_COMPLEX; // a plain double in the d arithmetic
    static void call(Struct *p) { dmumps_c(p); }
    static const char *name() { return "dmumps"; }
  };

  template <> struct MumpsTraits<std::complex<double>>
  {
    using Struct = ZMUMPS_STRUC_C;
    // struct {double r,i;} - the same object representation as std::complex<double>, which the
    // standard guarantees is a two-element array of the value type. So the value arrays are cast
    // rather than converted, and nothing is copied to change layout.
    using MumpsScalar = ZMUMPS_COMPLEX;
    static void call(Struct *p) { zmumps_c(p); }
    static const char *name() { return "zmumps"; }
  };

  /// One MUMPS instance: a matrix, its analysis and its factorisation.
  ///
  /// The lifetime rule that matters: this object OWNS its irn/jcn/a storage. The arrays pyoomph
  /// hands to a linear solver are zero-copy views into oomph-lib's CRDoubleMatrix buffers, which are
  /// freed and reallocated as soon as the Jacobian is reassembled - so keeping MUMPS's pointers
  /// aimed at them would be a use-after-free between the factorisation and the back-substitution.
  /// (pyoomph's Pardiso backend copies for exactly this reason; see the comment on
  /// PardisoSolver.get_jacobian_matrix.)
  template <typename Scalar> class MumpsSolver
  {
  public:
    using Traits = MumpsTraits<Scalar>;
    using MumpsScalar = typename Traits::MumpsScalar;

    /// @param sym          0 = unsymmetric (LU), 2 = general symmetric (LDL^T). Never 1: that is the
    ///                     positive-definite Cholesky, and positive definiteness is not something
    ///                     pyoomph can prove symbolically - only symmetry is.
    /// @param comm_fortran the Fortran handle of the communicator, i.e. MPI_Comm_c2f(comm). Ignored
    ///                     by a serial build.
    MumpsSolver(int sym, int comm_fortran, int verbosity);
    ~MumpsSolver();

    MumpsSolver(const MumpsSolver &) = delete;
    MumpsSolver &operator=(const MumpsSolver &) = delete;

    /// Install a centralized matrix from 0-based CSR (the layout pyoomph's serial solve_serial
    /// receives, under SuperLU's argument names).
    ///
    /// @return true when the sparsity pattern is bit-for-bit the one already installed, so that the
    ///         caller may skip analyse() and go straight to factorize(). The pattern is COMPARED,
    ///         never assumed: pyoomph's jacobian_structure_id promises it has not moved, but that
    ///         promise does not hold on an augmented (bifurcation-tracking) system, whose pattern is
    ///         value-filtered - and reusing an elimination tree for a pattern that has changed is a
    ///         silently wrong answer, not a crash.
    bool set_matrix_csr(int n, const int *row_start, int nrow, const int *col_index,
                        const Scalar *values);

    /// The distributed counterpart (MUMPS ICNTL(18)=3): this rank supplies the rows
    /// [first_row, first_row+nrow) of the global n x n matrix, with GLOBAL column indices - which is
    /// exactly what pyoomph's solve_distributed hands over.
    bool set_matrix_csr_distributed(int n, int first_row, const int *row_start, int nrow,
                                    const int *col_index, const Scalar *values);

    void analyse();   ///< JOB=1
    void factorize(); ///< JOB=2

    /// JOB=3, in place. @p rhs is the dense centralized right-hand side and must have n entries on
    /// the host; on the other ranks it is ignored and may be empty. Every rank must call this: the
    /// MUMPS phases are collective.
    void solve(Scalar *rhs, std::size_t rhs_len);

    void set_icntl(int i, int value); ///< 1-based, as in the MUMPS user guide
    int get_icntl(int i) const;
    int info(int i) const;
    int infog(int i) const;
    double rinfog(int i) const;

    /// Sign of the determinant, or 0 when it could not be determined. Requires ICNTL(33)=1 to have
    /// been set before the factorisation.
    int determinant_sign() const;

    int n() const { return id_.n; }
    std::size_t nnz() const { return irn_.size(); }
    bool is_distributed() const { return distributed_; }
    bool has_analysis() const { return analysed_; }
    bool has_factorisation() const { return factorised_; }
    static std::string version();

  private:
    void run_job(int job, const char *what);
    void check_error(const char *what);
    /// Expand CSR into MUMPS's 1-based COO triplets, keeping only the lower triangle when sym != 0.
    /// @return true if the resulting irn/jcn equal the ones already stored.
    bool build_pattern(int row_offset, const int *row_start, int nrow, const int *col_index,
                       const Scalar *values);

    typename Traits::Struct id_;
    std::vector<int> irn_, jcn_;
    std::vector<Scalar> a_;
    bool distributed_ = false;
    bool analysed_ = false;
    bool factorised_ = false;
    bool initialised_ = false;
    int sym_ = 0;
  };

  extern template class MumpsSolver<double>;
  extern template class MumpsSolver<std::complex<double>>;

} // namespace pyoomph_mumps
