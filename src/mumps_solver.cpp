#include "mumps_solver.hpp"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <type_traits>

#ifdef PYOOMPH_MUMPS_HAS_MPI
#include <mpi.h>
#endif

namespace pyoomph_mumps
{

  namespace
  {
    /// The MUMPS user guide's sentinel for "use MPI_COMM_WORLD", spelled the same way in every
    /// arithmetic. A serial (libmpiseq) build ignores comm_fortran entirely.
    constexpr int MUMPS_USE_COMM_WORLD = -987654;

    /// The negative INFOG(1) codes that are actually reachable from here, with what to do about
    /// them. MUMPS itself prints a line about most of these, but only when its diagnostic streams
    /// are on - and pyoomph runs it quiet, because a Newton solve that recovers from a failed
    /// factorisation would otherwise print a screenful per retry.
    std::string explain_infog(int infog1, int infog2)
    {
      switch (infog1)
      {
      case -5:
        return "out of memory during the analysis";
      case -6:
        return "the matrix is structurally singular";
      case -7:
        return "out of memory during the analysis (integer workspace)";
      case -8:
      case -9:
      case -14:
      case -15:
      case -17:
      case -20:
        return "the working space MUMPS estimated was too small; raise ICNTL(14) (the percentage "
               "increase in working space, default 20) and factorise again";
      case -10:
        return "the matrix is numerically singular";
      case -13:
        return "MUMPS could not allocate " + std::to_string(infog2) +
               (infog2 < 0 ? " million" : "") + " words of workspace - the problem is too large for "
                                                "the memory available";
      case -19:
        return "the maximum allowed size of working memory (ICNTL(23)) is too small";
      case -40:
        return "the matrix was declared symmetric positive definite but a negative pivot was found";
      case -44:
        return "the factors were discarded (ICNTL(31)); a solve cannot follow this factorisation";
      case -46:
      case -47:
        return "an unsupported combination of options for this solve phase";
      default:
        return "";
      }
    }
  } // namespace

  template <typename Scalar>
  MumpsSolver<Scalar>::MumpsSolver(int sym, int comm_fortran, int verbosity) : sym_(sym)
  {
    if (sym != 0 && sym != 2)
    {
      // sym==1 (Cholesky) is refused rather than silently mapped onto 2: pyoomph proves SYMMETRY
      // symbolically and cannot prove positive definiteness, so nothing upstream is in a position
      // to ask for it correctly, and a wrong yes here fails deep inside the factorisation.
      throw MumpsError("MUMPS sym must be 0 (unsymmetric) or 2 (general symmetric); sym=1 "
                       "(positive definite) is not offered because positive definiteness cannot be "
                       "established symbolically");
    }
#ifdef PYOOMPH_MUMPS_HAS_MPI
    {
      // MUMPS calls MPI_Comm_f2c() during JOB=-1, and MPI aborts the whole process if MPI_Init has
      // not run - "Local abort before MPI_INIT completed", with no exception to catch and no hint
      // about who was at fault. pyoomph_mumps/__init__.py imports mpi4py precisely so this cannot
      // normally happen; the check is for the case where somebody bypassed it (a bare `import
      // _pyoomph_mumps_core`, or an embedding application), and turns an abort into a sentence.
      int flag = 0;
      MPI_Initialized(&flag);
      if (!flag)
      {
        throw MumpsError("MPI has not been initialised, and this pyoomph_mumps was built with MPI "
                         "support, so MUMPS cannot start. Import pyoomph_mumps (or pyoomph, or "
                         "mpi4py.MPI) rather than _pyoomph_mumps_core directly - importing the "
                         "package initialises MPI.");
      }
    }
#endif
    std::memset(&id_, 0, sizeof(id_));
    id_.par = 1; // the host takes part in the factorisation as well as coordinating it
    id_.sym = sym;
    id_.comm_fortran = comm_fortran;
    id_.job = -1;
    Traits::call(&id_);
    initialised_ = true;
    check_error("initialising MUMPS");

    // Diagnostic streams (1-based ICNTL). All off at verbosity 0, error stream included: a failed
    // factorisation is not an error in pyoomph, it is a signal to shrink the step and try again, and
    // an arclength continuation through a fold would otherwise print a MUMPS error banner per retry.
    // Nothing is lost by it - INFOG(1) and INFOG(2), which is what those banners carry, are in the
    // exception this wrapper raises.
    id_.icntl[0] = (verbosity > 0) ? 6 : -1;        // error messages
    id_.icntl[1] = (verbosity > 1) ? 6 : -1;        // diagnostics, warnings and statistics
    id_.icntl[2] = (verbosity > 1) ? 6 : -1;        // global information
    id_.icntl[3] = (verbosity > 0) ? verbosity : 0; // 0 = say nothing at all
  }

  template <typename Scalar> MumpsSolver<Scalar>::~MumpsSolver()
  {
    if (initialised_)
    {
      id_.job = -2;
      // No check_error: a destructor that throws during stack unwinding terminates the process, and
      // there is nothing useful to do about a failed teardown anyway.
      Traits::call(&id_);
    }
  }

  template <typename Scalar> void MumpsSolver<Scalar>::check_error(const char *what)
  {
    if (id_.infog[0] >= 0) return;
    const int i1 = id_.infog[0], i2 = id_.infog[1];
    std::ostringstream oss;
    oss << "MUMPS failed while " << what << ": INFOG(1)=" << i1 << ", INFOG(2)=" << i2;
    const std::string why = explain_infog(i1, i2);
    if (!why.empty()) oss << " (" << why << ")";
    // Whatever state MUMPS is in after an error, it is not one a later phase may build on. Saying so
    // here is what forces the caller back through a full analyse+factorise rather than letting a
    // retry walk into the handle that just failed.
    analysed_ = false;
    factorised_ = false;
    throw MumpsError(oss.str(), i1, i2);
  }

  template <typename Scalar> void MumpsSolver<Scalar>::run_job(int job, const char *what)
  {
    id_.job = job;
    Traits::call(&id_);
    check_error(what);
  }

  template <typename Scalar>
  bool MumpsSolver<Scalar>::build_pattern(int row_offset, const int *row_start, int nrow,
                                          const int *col_index, const Scalar *values)
  {
    // oomph-lib hands out a local block whose row_start begins at 0, but reading the base from the
    // array rather than assuming it costs nothing and survives a caller that does not.
    const int base = (nrow > 0) ? row_start[0] : 0;
    const std::size_t nnz_in = (nrow > 0) ? static_cast<std::size_t>(row_start[nrow] - base) : 0;

    std::vector<int> irn, jcn;
    std::vector<Scalar> a;
    irn.reserve(nnz_in);
    jcn.reserve(nnz_in);
    a.reserve(nnz_in);

    for (int r = 0; r < nrow; ++r)
    {
      const int gr = row_offset + r + 1; // MUMPS indices are 1-based
      for (int k = row_start[r] - base; k < row_start[r + 1] - base; ++k)
      {
        const int gc = col_index[k] + 1;
        // For a symmetric factorisation MUMPS wants ONE triangle only; handing it both halves makes
        // it see every off-diagonal entry twice and add them together.
        if (sym_ != 0 && gc > gr) continue;
        irn.push_back(gr);
        jcn.push_back(gc);
        a.push_back(values[k]);
      }
    }

    const bool same = (irn.size() == irn_.size()) && (irn == irn_) && (jcn == jcn_);
    a_ = std::move(a);
    if (!same)
    {
      irn_ = std::move(irn);
      jcn_ = std::move(jcn);
      // A pattern change invalidates both the elimination tree and the factors built on it.
      analysed_ = false;
      factorised_ = false;
    }
    else
    {
      // Same pattern, new values: the analysis stands, the numbers in the factors do not.
      factorised_ = false;
    }
    return same;
  }

  template <typename Scalar>
  bool MumpsSolver<Scalar>::set_matrix_csr(int n, const int *row_start, int nrow,
                                           const int *col_index, const Scalar *values)
  {
    if (nrow != n)
    {
      throw MumpsError("A centralized MUMPS matrix must be square: got " + std::to_string(nrow) +
                       " rows for n=" + std::to_string(n));
    }
    const bool same = build_pattern(0, row_start, nrow, col_index, values);
    distributed_ = false;
    id_.n = n;
    id_.icntl[17] = 0; // ICNTL(18)=0: the whole matrix is supplied on the host
    id_.nz = 0;        // superseded by nnz, which is the 64-bit count MUMPS prefers
    id_.nnz = static_cast<MUMPS_INT8>(irn_.size());
    id_.irn = irn_.data();
    id_.jcn = jcn_.data();
    id_.a = reinterpret_cast<MumpsScalar *>(a_.data());
    return same;
  }

  template <typename Scalar>
  bool MumpsSolver<Scalar>::set_matrix_csr_distributed(int n, int first_row, const int *row_start,
                                                       int nrow, const int *col_index,
                                                       const Scalar *values)
  {
    const bool same = build_pattern(first_row, row_start, nrow, col_index, values);
    distributed_ = true;
    id_.n = n; // required on the host; harmless and simpler to set everywhere
    id_.icntl[17] = 3; // ICNTL(18)=3: every rank supplies its own share of the assembled matrix
    id_.nz_loc = 0;
    id_.nnz_loc = static_cast<MUMPS_INT8>(irn_.size());
    id_.irn_loc = irn_.data();
    id_.jcn_loc = jcn_.data();
    id_.a_loc = reinterpret_cast<MumpsScalar *>(a_.data());
    return same;
  }

  template <typename Scalar> void MumpsSolver<Scalar>::analyse()
  {
    run_job(1, "analysing the matrix");
    analysed_ = true;
    factorised_ = false;
  }

  template <typename Scalar> void MumpsSolver<Scalar>::factorize()
  {
    if (!analysed_)
    {
      throw MumpsError("MUMPS was asked to factorise before the matrix had been analysed. This is "
                       "the reuse path used wrongly: analyse() may only be skipped when "
                       "set_matrix_csr() reported that the sparsity pattern was unchanged.");
    }
    run_job(2, "factorising the matrix");
    factorised_ = true;
  }

  template <typename Scalar> void MumpsSolver<Scalar>::solve(Scalar *rhs, std::size_t rhs_len)
  {
    if (!factorised_)
    {
      throw MumpsError("MUMPS was asked to solve before the matrix had been factorised");
    }
    // The right-hand side is dense and centralized (ICNTL(20)=0) and so is the solution
    // (ICNTL(21)=0): only the host supplies and receives it. A rank with nothing to supply passes an
    // empty array, which is why this is a length check and not an assertion that it is n.
    if (rhs_len > 0 && rhs_len < static_cast<std::size_t>(id_.n))
    {
      throw MumpsError("The centralized right-hand side has " + std::to_string(rhs_len) +
                       " entries but the system has " + std::to_string(id_.n));
    }
    id_.nrhs = 1;
    id_.lrhs = id_.n;
    id_.rhs = (rhs_len > 0) ? reinterpret_cast<MumpsScalar *>(rhs) : nullptr;
    run_job(3, "solving with the factorised matrix");
    id_.rhs = nullptr;
  }

  template <typename Scalar> void MumpsSolver<Scalar>::set_icntl(int i, int value)
  {
    if (i < 1 || i > 60) throw MumpsError("ICNTL index out of range (1..60): " + std::to_string(i));
    id_.icntl[i - 1] = value;
  }

  template <typename Scalar> int MumpsSolver<Scalar>::get_icntl(int i) const
  {
    if (i < 1 || i > 60) throw MumpsError("ICNTL index out of range (1..60): " + std::to_string(i));
    return id_.icntl[i - 1];
  }

  template <typename Scalar> int MumpsSolver<Scalar>::info(int i) const
  {
    if (i < 1 || i > 80) throw MumpsError("INFO index out of range (1..80): " + std::to_string(i));
    return id_.info[i - 1];
  }

  template <typename Scalar> int MumpsSolver<Scalar>::infog(int i) const
  {
    if (i < 1 || i > 80) throw MumpsError("INFOG index out of range (1..80): " + std::to_string(i));
    return id_.infog[i - 1];
  }

  template <typename Scalar> double MumpsSolver<Scalar>::rinfog(int i) const
  {
    if (i < 1 || i > 40) throw MumpsError("RINFOG index out of range (1..40): " + std::to_string(i));
    return static_cast<double>(id_.rinfog[i - 1]);
  }

  template <typename Scalar> int MumpsSolver<Scalar>::determinant_sign() const
  {
    // Only meaningful for a real matrix: a complex determinant has an argument, not a sign.
    if (!std::is_same<Scalar, double>::value) return 0;
    if (id_.icntl[32] != 1) return 0; // ICNTL(33): the determinant was not computed
    const double mantissa = static_cast<double>(id_.rinfog[11]); // RINFOG(12)
    if (mantissa > 0.0) return 1;
    if (mantissa < 0.0) return -1;
    return 0;
  }

  template <typename Scalar> std::string MumpsSolver<Scalar>::version()
  {
    // MUMPS_VERSION is #defined by whichever arithmetic header was included first and guarded
    // against redefinition by the others, so d and z necessarily report the same string - they come
    // from one library.
    return std::string(MUMPS_VERSION);
  }

  template class MumpsSolver<double>;
  template class MumpsSolver<std::complex<double>>;

} // namespace pyoomph_mumps
