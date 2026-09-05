#  @file
#  @author Christian Diddens <c.diddens@utwente.nl>
#
#  @section LICENSE
#
#  pyoomph_mumps - direct MUMPS bindings for pyoomph
#  Copyright (C) 2021-2026  Christian Diddens
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the conditions of the BSD 3-Clause License are met.
#  See the LICENSE file distributed with this package for the full licence text.
#
# ========================================================================

"""Direct bindings to MUMPS (http://mumps-solver.org).

Deliberately knows nothing about pyoomph. Everything crossing this boundary is a numpy array or a
plain scalar, so the package builds, imports and passes its own tests on a machine where pyoomph is
not installed at all - which is what makes it possible to debug a MUMPS problem without a finite
element problem in the way. The pyoomph side of the coupling (the linear solver, the shift-invert
operator and the eigensolver) lives in ``pyoomph.solvers.mumps``.

The two solver classes differ only in their scalar type:

* :class:`MumpsSolverReal` wraps ``dmumps`` and takes float64 matrices;
* :class:`MumpsSolverComplex` wraps ``zmumps`` and takes complex128 ones.

Both follow the same three-phase MUMPS cycle, and the reason the phases are separate here rather
than bundled into one ``solve()`` is the whole point of the wrapper::

    unchanged = solver.set_matrix_csr(n, row_start, col_index, values)
    if not unchanged:
        solver.analyse()      # JOB=1: ordering and symbolic factorisation, the expensive part
    solver.factorize()        # JOB=2
    solver.solve(b)           # JOB=3, in place

``set_matrix_csr`` returns whether the sparsity pattern is bit-for-bit the one already installed, so
a transient run whose pattern never moves pays for the analysis once. It compares the pattern rather
than trusting the caller, because reusing an elimination tree for a pattern that has changed
produces a wrong answer instead of an error.
"""

from ._pyoomph_mumps_core import (  # type: ignore[import-not-found]
    MumpsSolverReal,
    MumpsSolverComplex,
    MumpsError,
    comm_fortran_from_address,
    has_mpi,
    has_openmp,
    mumps_version,
    USE_COMM_WORLD,
)

# An MPI-enabled MUMPS calls MPI_Comm_f2c() as soon as an instance is created, and MPI aborts the
# process outright if MPI_Init has not run. Importing mpi4py.MPI initialises MPI (and registers the
# finalisation), which is the same thing pyoomph's own pyoomph.generic.mpi does - and it is
# idempotent, so it does not matter which of the two is imported first.
if has_mpi:  # pragma: no cover - depends on how the extension was built
    try:
        from mpi4py import MPI as _MPI  # noqa: F401  # type: ignore[import-not-found]
    except ImportError as _e:  # pragma: no cover
        raise ImportError(
            "pyoomph_mumps was built with MPI support but mpi4py is not installed. MUMPS aborts the "
            "process rather than raising when MPI has not been initialised, so this is refused here "
            "instead. Install mpi4py, or rebuild with "
            "-DPYOOMPH_MUMPS_USE_MPI=OFF (which then requires a pyoomph built without MPI too)."
        ) from _e


__all__ = [
    "MumpsSolverReal",
    "MumpsSolverComplex",
    "MumpsError",
    "comm_fortran_from_address",
    "has_mpi",
    "has_openmp",
    "mumps_version",
    "USE_COMM_WORLD",
    "solver_for_dtype",
    "comm_fortran_for",
]

__version__ = "0.0.1"


def solver_for_dtype(dtype):
    """The solver class matching a numpy dtype: real for float64, complex for complex128.

    Callers should decide from the dtype of the matrix they are about to factorise, never from the
    matrices it was built out of. A shifted matrix ``J - sigma*M`` is complex whenever ``sigma`` is,
    while ``J`` and ``M`` stay real - which is the ordinary situation at a Hopf bifurcation, not an
    unusual one.
    """
    import numpy

    kind = numpy.dtype(dtype).kind
    if kind == "c":
        return MumpsSolverComplex
    if kind == "f":
        return MumpsSolverReal
    raise TypeError("MUMPS needs a float64 or complex128 matrix, got dtype " + str(dtype))


def comm_fortran_for(comm=None):
    """The ``comm_fortran`` handle for an mpi4py communicator, or the MPI_COMM_WORLD sentinel.

    Returns :data:`USE_COMM_WORLD` when ``comm`` is None or this build has no MPI, so a caller does
    not have to branch on either.
    """
    if comm is None or not has_mpi:
        return USE_COMM_WORLD
    from mpi4py import MPI  # type: ignore[import-not-found]

    return comm_fortran_from_address(MPI._addressof(comm))
