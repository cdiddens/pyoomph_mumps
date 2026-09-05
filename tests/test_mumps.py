"""The standalone test suite: numpy and scipy only, no pyoomph.

That independence is the point. When a run goes wrong it is rarely obvious whether MUMPS, the
binding or the finite element problem is at fault, and these tests settle the first two without the
third in the way.
"""

import numpy
import pytest
import scipy.sparse
import scipy.sparse.linalg

import pyoomph_mumps as pm


def _csr(A):
    """The three arrays in the types the binding wants: int32 indices, native-dtype values."""
    A = A.tocsr()
    return A.shape[0], A.indptr.astype(numpy.int32), A.indices.astype(numpy.int32), A.data


def _random_system(n=200, density=0.03, seed=1, complex_=False):
    # The +5*I is not cosmetic: a random sparse matrix at this density is very often structurally
    # singular, and a test that fails for that reason tests nothing.
    A = scipy.sparse.random(n, n, density=density, format="csr", random_state=seed)
    A = (A + 5.0 * scipy.sparse.eye(n)).tocsr()
    rng = numpy.random.default_rng(seed)
    b = rng.standard_normal(n)
    if complex_:
        A = A.astype(numpy.complex128)
        A.data += 0.4j * numpy.linspace(-1.0, 1.0, A.data.size)
        b = b + 1j * rng.standard_normal(n)
    return A, b


def _solve(solver, A, b):
    n, indptr, indices, data = _csr(A)
    reused = solver.set_matrix_csr(n, indptr, indices, data)
    if not reused:
        solver.analyse()
    solver.factorize()
    x = b.copy()
    solver.solve(x)
    return x, reused


@pytest.mark.parametrize("complex_", [False, True])
def test_solve_matches_scipy(complex_):
    A, b = _random_system(complex_=complex_)
    cls = pm.MumpsSolverComplex if complex_ else pm.MumpsSolverReal
    x, _ = _solve(cls(sym=0), A, b)
    ref = scipy.sparse.linalg.spsolve(A.tocsc(), b)
    assert numpy.max(numpy.abs(x - ref)) < 1e-9 * max(1.0, numpy.max(numpy.abs(ref)))


def test_solver_for_dtype_follows_the_matrix():
    assert pm.solver_for_dtype(numpy.float64) is pm.MumpsSolverReal
    assert pm.solver_for_dtype(numpy.complex128) is pm.MumpsSolverComplex
    with pytest.raises(TypeError):
        pm.solver_for_dtype(numpy.int32)


def test_pattern_reuse_skips_the_analysis():
    """New values on an unchanged pattern must reuse the analysis and still give the right answer."""
    A, b = _random_system()
    s = pm.MumpsSolverReal(sym=0)
    x1, reused1 = _solve(s, A, b)
    assert reused1 is False  # nothing was installed before

    A2 = A.copy()
    A2.data = A.data * 1.7 + 0.25
    n, indptr, indices, data = _csr(A2)
    reused2 = s.set_matrix_csr(n, indptr, indices, data)
    assert reused2 is True
    # The claim under test: the analysis survived the new values, so factorize() alone is legal.
    assert s.has_analysis is True
    s.factorize()
    x2 = b.copy()
    s.solve(x2)
    ref2 = scipy.sparse.linalg.spsolve(A2.tocsc(), b)
    assert numpy.allclose(x2, ref2, atol=1e-9)
    # ... and the first answer was not quietly reused for the second system.
    assert not numpy.allclose(x1, x2)


def test_a_changed_pattern_is_detected():
    """A pattern that has moved must be reported, even when the number of entries is the same.

    This is the case that makes trusting a caller's "the pattern is unchanged" promise unsafe: the
    sizes match, so nothing but comparing the indices themselves can tell the two apart.
    """
    A, _ = _random_system(n=50, density=0.1)
    A = A.tolil()
    B = A.copy()
    # Move one entry, keeping the count identical.
    rows, cols = A.nonzero()
    r, c = rows[0], cols[0]
    free = next(j for j in range(50) if B[r, j] == 0)
    B[r, c] = 0.0
    B[r, free] = 1.0
    s = pm.MumpsSolverReal(sym=0)
    n, indptr, indices, data = _csr(A)
    assert s.set_matrix_csr(n, indptr, indices, data) is False
    s.analyse()
    n, indptr, indices, data = _csr(B)
    assert A.nnz == B.nnz
    assert s.set_matrix_csr(n, indptr, indices, data) is False
    assert s.has_analysis is False


def test_symmetric_matches_unsymmetric():
    """sym=2 takes only the lower triangle; the answer must not depend on which path was used."""
    A, b = _random_system(n=150, seed=3)
    A = ((A + A.T) * 0.5).tocsr()
    x_gen, _ = _solve(pm.MumpsSolverReal(sym=0), A, b)
    x_sym, _ = _solve(pm.MumpsSolverReal(sym=2), A, b)
    assert numpy.allclose(x_gen, x_sym, atol=1e-9)
    ref = scipy.sparse.linalg.spsolve(A.tocsc(), b)
    assert numpy.allclose(x_sym, ref, atol=1e-9)


def test_sym_one_is_refused():
    # Not an arbitrary restriction: positive definiteness cannot be established symbolically, so
    # nothing upstream is in a position to ask for the Cholesky correctly.
    with pytest.raises(pm.MumpsError):
        pm.MumpsSolverReal(sym=1)


def test_singular_matrix_raises_rather_than_returning_garbage():
    n = 60
    A = scipy.sparse.eye(n, format="lil")
    A[3, 3] = 0.0  # a structurally present but numerically zero pivot
    A[3, 4] = 0.0
    A = A.tocsr()
    b = numpy.ones(n)
    s = pm.MumpsSolverReal(sym=0)
    with pytest.raises(pm.MumpsError):
        _solve(s, A, b)


def test_phase_order_is_enforced():
    A, b = _random_system(n=40)
    s = pm.MumpsSolverReal(sym=0)
    n, indptr, indices, data = _csr(A)
    s.set_matrix_csr(n, indptr, indices, data)
    with pytest.raises(pm.MumpsError):
        s.factorize()  # no analysis yet
    s.analyse()
    with pytest.raises(pm.MumpsError):
        s.solve(b.copy())  # no factorisation yet


def test_rhs_length_is_checked():
    A, b = _random_system(n=40)
    s = pm.MumpsSolverReal(sym=0)
    _solve(s, A, b)
    with pytest.raises(pm.MumpsError):
        s.solve(numpy.ones(10))
    # An empty right-hand side is legal, not an error: that is what a non-host rank passes in a
    # distributed solve, where the solution is centralized on the host.
    s.solve(numpy.zeros(0))


def test_non_square_centralized_matrix_is_refused():
    A, _ = _random_system(n=40)
    n, indptr, indices, data = _csr(A)
    s = pm.MumpsSolverReal(sym=0)
    with pytest.raises(pm.MumpsError):
        s.set_matrix_csr(n + 1, indptr, indices, data)


def test_determinant_sign():
    # Only the SIGN, and only its changes, are meaningful - see the docstring of
    # GenericLinearSystemSolver.get_determinant_sign. Here the two matrices differ by one row
    # negation, so the signs must differ.
    n = 30
    A = scipy.sparse.eye(n, format="lil")
    A[0, 0] = 2.0
    B = A.copy()
    B[0, 0] = -2.0
    signs = []
    for M in (A.tocsr(), B.tocsr()):
        s = pm.MumpsSolverReal(sym=0)
        s.set_icntl(33, 1)  # compute the determinant
        n_, indptr, indices, data = _csr(M)
        s.set_matrix_csr(n_, indptr, indices, data)
        s.analyse()
        s.factorize()
        signs.append(s.determinant_sign())
    assert signs[0] != 0 and signs[1] != 0
    assert signs[0] == -signs[1]


def test_determinant_sign_is_unavailable_without_icntl33():
    A, _ = _random_system(n=30)
    s = pm.MumpsSolverReal(sym=0)
    n, indptr, indices, data = _csr(A)
    s.set_matrix_csr(n, indptr, indices, data)
    s.analyse()
    s.factorize()
    assert s.determinant_sign() == 0


def test_determinant_sign_is_unavailable_for_complex():
    A, _ = _random_system(n=30, complex_=True)
    s = pm.MumpsSolverComplex(sym=0)
    s.set_icntl(33, 1)
    n, indptr, indices, data = _csr(A)
    s.set_matrix_csr(n, indptr, indices, data)
    s.analyse()
    s.factorize()
    # A complex determinant has an argument, not a sign.
    assert s.determinant_sign() == 0


def test_icntl_round_trip_and_range():
    s = pm.MumpsSolverReal(sym=0)
    s.set_icntl(14, 35)
    assert s.get_icntl(14) == 35
    with pytest.raises(pm.MumpsError):
        s.set_icntl(0, 1)
    with pytest.raises(pm.MumpsError):
        s.get_icntl(61)


def test_unsorted_column_indices_are_accepted():
    """oomph-lib does not emit the columns of a row in ascending order, and MUMPS does not need it.

    Worth asserting rather than assuming: the CSR->COO expansion here never sorts, so if MUMPS ever
    did care, every pyoomph solve would be quietly wrong rather than failing.
    """
    A, b = _random_system(n=80, seed=5)
    A = A.tocsr()
    A.has_sorted_indices = False
    shuffled = A.copy()
    rng = numpy.random.default_rng(0)
    for r in range(shuffled.shape[0]):
        lo, hi = shuffled.indptr[r], shuffled.indptr[r + 1]
        perm = rng.permutation(hi - lo)
        shuffled.indices[lo:hi] = shuffled.indices[lo:hi][perm]
        shuffled.data[lo:hi] = shuffled.data[lo:hi][perm]
    shuffled.has_sorted_indices = False
    x, _ = _solve(pm.MumpsSolverReal(sym=0), shuffled, b)
    ref = scipy.sparse.linalg.spsolve(A.tocsc(), b)
    assert numpy.allclose(x, ref, atol=1e-9)


def test_module_reports_its_build():
    assert isinstance(pm.mumps_version, str) and pm.mumps_version
    assert isinstance(pm.has_mpi, bool)
    assert isinstance(pm.has_openmp, bool)
