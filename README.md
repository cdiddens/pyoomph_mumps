# pyoomph_mumps

Direct [MUMPS](http://mumps-solver.org) bindings for [pyoomph](https://github.com/pyoomph/pyoomph).

pyoomph can already use MUMPS, but only through PETSc/SLEPc, which means building and maintaining a
full PETSc+SLEPc stack in two scalar flavours (a real one for linear solves, a complex one for
azimuthal and Floquet stability). This package talks to `dmumps` and `zmumps` directly instead. With
it installed, pyoomph gains

* the linear solver `"mumps"` - `problem.set_linear_solver("mumps")` - serial or, under MPI,
  natively distributed (MUMPS `ICNTL(18)=3`), including the sparsity-pattern reuse that lets a
  transient run pay for the ordering once;
* the eigensolver `"mumps"` - `problem.set_eigen_solver("mumps")` - which is pyoomph's existing
  Spectra Arnoldi driver with MUMPS supplying the shift-and-invert factorisation, in real *and*
  complex arithmetic.

Without it, `pyoomph.solvers.mumps` raises `ImportError` at import and neither backend is registered,
so nothing else changes.

## Installing

    pip install .

By default this downloads and builds MUMPS from source (via the CMake superbuild at
[scivision/mumps](https://github.com/scivision/mumps), which fetches the genuine upstream sources),
along with a reference BLAS/LAPACK if the machine has none. That needs a Fortran compiler.

To use a MUMPS that is already installed - much faster, and the right choice on a machine with
`libmumps-dev` or a module-loaded MUMPS:

    pip install . --config-settings=cmake.define.PYOOMPH_MUMPS_DOWNLOAD=OFF

or point at a specific prefix with `-DPYOOMPH_MUMPS_DIR=/path/to/prefix`. Either way the tree must
carry **both** the real and the double-complex arithmetic (`libdmumps` and `libzmumps`).

### Options

| CMake option | Default | Meaning |
|---|---|---|
| `PYOOMPH_MUMPS_USE_MPI` | `OFF` | **Must match how pyoomph itself was built.** See below. |
| `PYOOMPH_MUMPS_DOWNLOAD` | `ON` | Download and build MUMPS. |
| `PYOOMPH_MUMPS_DIR` | *(empty)* | Use the MUMPS installed under this prefix. |
| `PYOOMPH_MUMPS_USE_OPENMP` | `AUTO` | Build MUMPS with OpenMP, so a factorisation can use several cores per rank (`ICNTL(16)`). `ON` fails if OpenMP is missing. |
| `PYOOMPH_MUMPS_BLAS_VENDOR` | *(empty)* | `OpenBLAS`, `MKL`, `AOCL`, `Atlas`, `Netlib`. |
| `PYOOMPH_MUMPS_BLAS_ROOT` | *(empty)* | Where to look for it. |
| `PYOOMPH_MUMPS_DOWNLOAD_BLAS` | `AUTO` | Build the reference BLAS/LAPACK when none is found. |
| `PYOOMPH_MUMPS_EXTRA_LIBRARIES` | *(empty)* | Extra link libraries, for a static prebuilt tree whose dependencies are recorded nowhere. |

`build_for_develop.sh` wraps all of this for an editable install, reads the MPI setting off pyoomph's
own build, and defaults to the system MUMPS because that is the fast loop; set
`PYOOMPH_MUMPS_BUILD_MODE=download` to exercise the download route.

### The one configuration that goes wrong quietly

`PYOOMPH_MUMPS_USE_MPI` has to agree with pyoomph. A serial MUMPS links a stub library that defines
its own `MPI_Init` and friends; put that in the same process as an MPI-enabled pyoomph and the two
either initialise MPI twice or not at all, and the result is a hang or a wrong communicator rather
than an error message. `pyoomph.solvers.mumps` therefore compares the two at import time and refuses
a mismatched pair outright. `build_for_develop.sh` reads pyoomph's setting and matches it
automatically.

## Standalone use

The package does not import pyoomph and is usable on its own:

```python
import numpy, scipy.sparse, pyoomph_mumps

A = scipy.sparse.random(500, 500, density=0.01, format="csr") + scipy.sparse.eye(500)
b = numpy.ones(500)

s = pyoomph_mumps.MumpsSolverReal(sym=0)
s.set_matrix_csr(A.shape[0], A.indptr.astype(numpy.int32),
                 A.indices.astype(numpy.int32), A.data)
s.analyse()
s.factorize()
x = b.copy()
s.solve(x)          # in place
```

`set_matrix_csr` returns `True` when the sparsity pattern is exactly the one already installed, which
is the signal that `analyse()` may be skipped. It **compares** the pattern rather than trusting the
caller: reusing an elimination tree for a pattern that has moved is a wrong answer, not a crash.

## Licence

BSD-3-Clause (see [LICENSE](LICENSE)); note that pyoomph itself is GPL-3.0-or-later. MUMPS is
distributed under its own (CeCILL-C) licence and is not included here - it is downloaded, or taken
from the system, at build time.
