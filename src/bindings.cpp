// The nanobind layer: numpy arrays in, numpy arrays out, and nothing else. All the MUMPS knowledge
// lives in mumps_solver.hpp/.cpp so that it can be read - and tested - without nanobind in the way.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <complex>
#include <cstdint>

#include "mumps_solver.hpp"

#ifdef PYOOMPH_MUMPS_HAS_MPI
#include <mpi.h>
#endif

namespace nb = nanobind;
using namespace nb::literals;
using pyoomph_mumps::MumpsError;
using pyoomph_mumps::MumpsSolver;

namespace
{
  template <typename Scalar> using ConstVec = nb::ndarray<const Scalar, nb::ndim<1>, nb::c_contig>;
  template <typename Scalar> using MutVec = nb::ndarray<Scalar, nb::ndim<1>, nb::c_contig>;
  using ConstIntVec = nb::ndarray<const int, nb::ndim<1>, nb::c_contig>;

  void check_csr(const ConstIntVec &row_start, int nrow, const ConstIntVec &col_index,
                 std::size_t nvalues)
  {
    if (row_start.shape(0) != static_cast<std::size_t>(nrow) + 1)
    {
      throw MumpsError("row_start must have nrow+1 = " + std::to_string(nrow + 1) +
                       " entries, got " + std::to_string(row_start.shape(0)));
    }
    if (col_index.shape(0) != nvalues)
    {
      throw MumpsError("col_index and values must have the same length, got " +
                       std::to_string(col_index.shape(0)) + " and " + std::to_string(nvalues));
    }
  }

  template <typename Scalar> void bind_solver(nb::module_ &m, const char *name, const char *doc)
  {
    using Solver = MumpsSolver<Scalar>;
    nb::class_<Solver>(m, name, doc)
        .def(nb::init<int, int, int>(), "sym"_a = 0, "comm_fortran"_a = -987654, "verbosity"_a = 0,
             "Create a MUMPS instance (JOB=-1).\n\n"
             "sym: 0 = unsymmetric (LU), 2 = general symmetric (LDL^T). sym=1 (Cholesky) is\n"
             "     deliberately not offered - see the C++ comment.\n"
             "comm_fortran: MPI_Comm_c2f() of the communicator; the default is MUMPS's own\n"
             "     'use MPI_COMM_WORLD' sentinel. Use comm_fortran_from_address() to obtain one\n"
             "     from an mpi4py communicator. Ignored by a build without MPI.\n"
             "verbosity: 0 keeps only MUMPS's error stream, higher values switch its diagnostic\n"
             "     and statistics streams back on.")

        .def(
            "set_matrix_csr",
            [](Solver &s, int n, ConstIntVec row_start, ConstIntVec col_index,
               ConstVec<Scalar> values)
            {
              const int nrow = static_cast<int>(row_start.shape(0)) - 1;
              check_csr(row_start, nrow, col_index, values.shape(0));
              return s.set_matrix_csr(n, row_start.data(), nrow, col_index.data(), values.data());
            },
            "n"_a, "row_start"_a, "col_index"_a, "values"_a,
            "Install a centralized matrix from 0-based CSR.\n\n"
            "Returns True when the sparsity pattern is exactly the one already installed, in which\n"
            "case analyse() may be skipped and factorize() called directly. The pattern is\n"
            "compared, not assumed - reusing an elimination tree for a pattern that has moved is a\n"
            "silently wrong answer rather than a crash.\n\n"
            "The arrays are copied: MUMPS dereferences its pointers again at the factorisation and\n"
            "solve phases, by which time a caller that passed a view of a matrix it has since\n"
            "reassembled would have freed them.")

        .def(
            "set_matrix_csr_distributed",
            [](Solver &s, int n, int first_row, ConstIntVec row_start, ConstIntVec col_index,
               ConstVec<Scalar> values)
            {
              const int nrow = static_cast<int>(row_start.shape(0)) - 1;
              check_csr(row_start, nrow, col_index, values.shape(0));
              return s.set_matrix_csr_distributed(n, first_row, row_start.data(), nrow,
                                                  col_index.data(), values.data());
            },
            "n"_a, "first_row"_a, "row_start"_a, "col_index"_a, "values"_a,
            "Install this rank's row block of a distributed matrix (MUMPS ICNTL(18)=3).\n\n"
            "row_start is local and 0-based; col_index holds GLOBAL column indices. Returns the\n"
            "same pattern-unchanged flag as set_matrix_csr(). Collective: every rank must call it,\n"
            "including one that owns no rows.")

        .def(
            "analyse", [](Solver &s) { nb::gil_scoped_release r; s.analyse(); },
            "Run the analysis phase (JOB=1): ordering, symbolic factorisation, task mapping.")
        .def(
            "factorize", [](Solver &s) { nb::gil_scoped_release r; s.factorize(); },
            "Run the numerical factorisation (JOB=2). Requires a preceding analysis.")

        .def(
            "solve",
            [](Solver &s, MutVec<Scalar> rhs)
            {
              // In place, like pyoomph's own solve_serial contract, so a Newton step needs no extra
              // copy of the residual vector.
              Scalar *ptr = rhs.shape(0) ? rhs.data() : nullptr;
              const std::size_t len = rhs.shape(0);
              nb::gil_scoped_release r;
              s.solve(ptr, len);
            },
            "rhs"_a,
            "Back-substitute (JOB=3), overwriting rhs with the solution.\n\n"
            "The right-hand side is dense and centralized: in a distributed solve only the host\n"
            "passes an array of length n, the other ranks pass an empty one. Collective - every\n"
            "rank must call it.")

        .def("set_icntl", &Solver::set_icntl, "i"_a, "value"_a,
             "Set ICNTL(i), 1-based as in the MUMPS user guide.")
        .def("get_icntl", &Solver::get_icntl, "i"_a, "Read ICNTL(i), 1-based.")
        .def("info", &Solver::info, "i"_a, "Read INFO(i) (this rank), 1-based.")
        .def("infog", &Solver::infog, "i"_a, "Read INFOG(i) (global), 1-based.")
        .def("rinfog", &Solver::rinfog, "i"_a, "Read RINFOG(i) (global), 1-based.")
        .def("determinant_sign", &Solver::determinant_sign,
             "Sign of the determinant (+1/-1/0), or 0 when unavailable.\n\n"
             "Requires ICNTL(33)=1 to have been set before the factorisation, and is 0 for the\n"
             "complex arithmetic, where a determinant has an argument rather than a sign.")
        .def_prop_ro("n", &Solver::n, "Order of the system.")
        .def_prop_ro("nnz", &Solver::nnz, "Number of stored entries on this rank.")
        .def_prop_ro("is_distributed", &Solver::is_distributed)
        .def_prop_ro("has_analysis", &Solver::has_analysis)
        .def_prop_ro("has_factorisation", &Solver::has_factorisation);
  }
} // namespace

NB_MODULE(_pyoomph_mumps_core, m)
{
  m.doc() = "Direct bindings to MUMPS (http://mumps-solver.org) for pyoomph.";

  nb::exception<MumpsError>(m, "MumpsError");

  bind_solver<double>(m, "MumpsSolverReal", "A MUMPS instance over real (float64) matrices - dmumps.");
  bind_solver<std::complex<double>>(
      m, "MumpsSolverComplex",
      "A MUMPS instance over complex (complex128) matrices - zmumps.\n\n"
      "Not an exotic extra: the eigensolver's shift-and-invert factorises J - sigma*M, which is\n"
      "complex for every complex shift even when J and M are real - i.e. for every Hopf\n"
      "bifurcation and every azimuthal stability problem.");

#ifdef PYOOMPH_MUMPS_HAS_MPI
  m.attr("has_mpi") = true;
  m.def(
      "comm_fortran_from_address",
      [](std::uintptr_t address)
      {
        // The address of the caller's MPI_Comm variable, not the handle itself: mpi4py's
        // MPI._addressof(comm) gives a pointer into the object, and dereferencing it here is what
        // makes this work for both ABIs - Open MPI's MPI_Comm is a pointer, MPICH's is an int, and
        // neither can be passed through a single Python integer without knowing which.
        MPI_Comm comm = *reinterpret_cast<MPI_Comm *>(address);
        return static_cast<int>(MPI_Comm_c2f(comm));
      },
      "address"_a,
      "Turn the address of an MPI_Comm (mpi4py's MPI._addressof(comm)) into the Fortran handle\n"
      "the solver constructor wants.");
#else
  m.attr("has_mpi") = false;
  m.def(
      "comm_fortran_from_address", [](std::uintptr_t) -> int
      { throw MumpsError("this pyoomph_mumps was built without MPI support"); },
      "address"_a, "Unavailable: this build has no MPI support.");
#endif

#ifdef PYOOMPH_MUMPS_HAS_OPENMP
  m.attr("has_openmp") = true;
#else
  m.attr("has_openmp") = false;
#endif

  m.attr("mumps_version") = MumpsSolver<double>::version();
  m.attr("USE_COMM_WORLD") = -987654;
#ifdef PYOOMPH_MUMPS_BUILD_VERSION
  m.attr("__version__") = PYOOMPH_MUMPS_BUILD_VERSION;
#endif
}
