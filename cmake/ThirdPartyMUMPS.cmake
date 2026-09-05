# Finds (or downloads and builds) MUMPS - http://mumps-solver.org - together with the BLAS/LAPACK it
# needs.
#
# Three mutually exclusive routes, in the order they are tried:
#
#   1. PYOOMPH_MUMPS_DIR       an already installed MUMPS tree (a system one, or a prefix produced by
#                              a previous download build). Used in place, never touched.
#   2. PYOOMPH_MUMPS_DOWNLOAD  (the default) download and build MUMPS from source.
#   3. neither                 look for a system MUMPS on the default search paths, and fail with an
#                              explanatory message when there is none.
#
# Sets, for use by the top-level CMakeLists.txt:
#   PYOOMPH_MUMPS_INCLUDE_DIRS   include directories carrying dmumps_c.h / zmumps_c.h
#   PYOOMPH_MUMPS_LIBRARIES      what the extension module must link
#   PYOOMPH_MUMPS_DEPENDS        ExternalProject targets the extension must depend on (empty when
#                                nothing was downloaded)
#   PYOOMPH_MUMPS_NEEDS_FORTRAN  TRUE when the Fortran runtime has to be linked explicitly, i.e. when
#                                the MUMPS libraries are static
#
# WHY NOT THE UPSTREAM TARBALL. MUMPS itself is not a CMake project and never has been: its build is
# driven by a Makefile.inc that the user is expected to copy from Make.inc/ and edit by hand for the
# local compilers, BLAS, ScaLAPACK and ordering libraries. Driving that from ExternalProject means
# generating that file ourselves and keeping it correct across four platforms, which is a great deal
# of machinery to maintain for something somebody else already maintains. So the download route goes
# through the CMake superbuild at https://github.com/scivision/mumps, which fetches the genuine
# upstream MUMPS sources from mumps-solver.org and only supplies the build system. It is a build
# wrapper, not a fork - the solver code compiled is upstream's.

include(ExternalProject)
include(GNUInstallDirs)

set(PYOOMPH_MUMPS_DIR "" CACHE PATH
    "An already installed MUMPS tree (containing include/dmumps_c.h and lib/libdmumps) to use instead of downloading one. Used in place and left untouched.")
option(PYOOMPH_MUMPS_DOWNLOAD "Download and build MUMPS from source via the scivision/mumps CMake superbuild" ON)

# Pinned, exactly as cmake/ThirdPartySpectra.cmake pins Spectra: the option names this file passes
# below (MUMPS_parallel, MUMPS_openmp, BUILD_COMPLEX16, ...) are the pinned ref's, and a floating
# branch is free to rename them - at which point CMake silently ignores the unknown -D and builds
# something other than what was asked for, with no error anywhere.
set(PYOOMPH_MUMPS_SUPERBUILD_REF "31999e0f398d7f486c23c265ba5c40f8e1288a3c" CACHE STRING
    "Git ref of the scivision/mumps CMake superbuild to download")
set(PYOOMPH_MUMPS_SUPERBUILD_URL "" CACHE STRING
    "Full URL of a scivision/mumps source archive, overriding PYOOMPH_MUMPS_SUPERBUILD_REF")
# Left empty on purpose: the superbuild then uses its own default, which is the version its
# cmake/source.json has a SHA256 for. Pinning a different one here is only safe while
# mumps-solver.org still hosts that tarball, and it does not keep old ones indefinitely.
set(PYOOMPH_MUMPS_UPSTREAM_VERSION "" CACHE STRING
    "MUMPS version for the download route (empty = the superbuild's own default)")

set(PYOOMPH_MUMPS_EXTRA_LIBRARIES "" CACHE STRING
    "Extra libraries to link after the MUMPS ones. Needed only for a STATIC prebuilt tree supplied through PYOOMPH_MUMPS_DIR, whose ScaLAPACK/BLAS/ordering dependencies are not recorded anywhere.")

# --- BLAS/LAPACK -------------------------------------------------------------------------------
# The superbuild does `find_package(LAPACK REQUIRED)` and, unlike ScaLAPACK, has NO download fallback
# for it - a machine without a BLAS fails its configure with a message about LAPACK that says nothing
# about MUMPS. So the fallback is here.
set(PYOOMPH_MUMPS_DOWNLOAD_BLAS "AUTO" CACHE STRING
    "Reference BLAS/LAPACK for the download route: AUTO=only if none is found, ON=always, OFF=never (fail instead)")
set_property(CACHE PYOOMPH_MUMPS_DOWNLOAD_BLAS PROPERTY STRINGS AUTO ON OFF)
set(PYOOMPH_MUMPS_BLAS_VENDOR "" CACHE STRING
    "BLAS/LAPACK vendor, passed on as LAPACK_VENDOR: Netlib (default), OpenBLAS, MKL, AOCL, Atlas")
set(PYOOMPH_MUMPS_BLAS_ROOT "" CACHE PATH
    "Where to look for the BLAS/LAPACK, passed on as LAPACK_ROOT")
set(PYOOMPH_MUMPS_LAPACK_REF "v3.12.1" CACHE STRING
    "Git ref of Reference-LAPACK to download when PYOOMPH_MUMPS_DOWNLOAD_BLAS applies")

set(PYOOMPH_MUMPS_INCLUDE_DIRS "")
set(PYOOMPH_MUMPS_LIBRARIES "")
set(PYOOMPH_MUMPS_DEPENDS "")
set(PYOOMPH_MUMPS_NEEDS_FORTRAN FALSE)

# ================================================================================================
# A small FindMUMPS. CMake ships none, and neither MUMPS nor Debian's libmumps-dev installs a
# pkg-config file or a CMake package config, so there is nothing to find_package().
# ================================================================================================
function(pyoomph_find_mumps hint_dir out_found out_includes out_libraries out_is_static)
  set(_hints)
  if(hint_dir)
    set(_hints HINTS "${hint_dir}" "${hint_dir}/include" "${hint_dir}/lib" NO_DEFAULT_PATH)
  endif()

  # dmumps_c.h and zmumps_c.h are asked for SEPARATELY rather than assumed to travel together: the
  # arithmetics are independent build options upstream (BUILD_DOUBLE / BUILD_COMPLEX16), so a tree
  # with only the real one is a perfectly ordinary MUMPS install - it just cannot do the complex
  # shift-invert this package exists to provide, and saying so here beats a link error later.
  find_path(PYOOMPH_MUMPS_D_INCLUDE_DIR NAMES dmumps_c.h ${_hints})
  find_path(PYOOMPH_MUMPS_Z_INCLUDE_DIR NAMES zmumps_c.h ${_hints})
  find_library(PYOOMPH_MUMPS_D_LIBRARY NAMES dmumps ${_hints})
  find_library(PYOOMPH_MUMPS_Z_LIBRARY NAMES zmumps ${_hints})
  find_library(PYOOMPH_MUMPS_COMMON_LIBRARY NAMES mumps_common ${_hints})
  # Only present in a serial build; a parallel MUMPS has no libmpiseq and must not link one.
  find_library(PYOOMPH_MUMPS_MPISEQ_LIBRARY NAMES mpiseq ${_hints})
  # MUMPS >= 5.3 keeps its portability shims here; older versions fold them into mumps_common.
  find_library(PYOOMPH_MUMPS_PORD_LIBRARY NAMES pord ${_hints})

  if(NOT PYOOMPH_MUMPS_D_INCLUDE_DIR OR NOT PYOOMPH_MUMPS_D_LIBRARY OR NOT PYOOMPH_MUMPS_COMMON_LIBRARY)
    set(${out_found} FALSE PARENT_SCOPE)
    return()
  endif()
  if(NOT PYOOMPH_MUMPS_Z_INCLUDE_DIR OR NOT PYOOMPH_MUMPS_Z_LIBRARY)
    message(FATAL_ERROR
      "Found a MUMPS with the real (dmumps) arithmetic but not the double-complex one (zmumps). "
      "pyoomph_mumps needs both: the complex shift-and-invert of the eigensolver factorises "
      "J - sigma*M, which is complex for any complex shift even when J and M are real - i.e. for "
      "every Hopf bifurcation and every azimuthal stability problem. Rebuild MUMPS with "
      "BUILD_COMPLEX16=ON, install your distribution's complex MUMPS, or leave PYOOMPH_MUMPS_DIR "
      "empty to let this build download and build a complete one.")
  endif()

  set(_libs "${PYOOMPH_MUMPS_D_LIBRARY}" "${PYOOMPH_MUMPS_Z_LIBRARY}" "${PYOOMPH_MUMPS_COMMON_LIBRARY}")
  if(PYOOMPH_MUMPS_PORD_LIBRARY)
    list(APPEND _libs "${PYOOMPH_MUMPS_PORD_LIBRARY}")
  endif()
  if(PYOOMPH_MUMPS_MPISEQ_LIBRARY AND NOT PYOOMPH_MUMPS_USE_MPI)
    list(APPEND _libs "${PYOOMPH_MUMPS_MPISEQ_LIBRARY}")
  endif()

  # Static libraries record none of their own dependencies, so the caller has to know to link the
  # Fortran runtime (and, via PYOOMPH_MUMPS_EXTRA_LIBRARIES, ScaLAPACK and the BLAS) after them.
  get_filename_component(_ext "${PYOOMPH_MUMPS_D_LIBRARY}" EXT)
  if(_ext STREQUAL "${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(${out_is_static} TRUE PARENT_SCOPE)
  else()
    set(${out_is_static} FALSE PARENT_SCOPE)
  endif()

  set(_incs "${PYOOMPH_MUMPS_D_INCLUDE_DIR}" "${PYOOMPH_MUMPS_Z_INCLUDE_DIR}")
  list(REMOVE_DUPLICATES _incs)
  set(${out_found} TRUE PARENT_SCOPE)
  set(${out_includes} "${_incs}" PARENT_SCOPE)
  set(${out_libraries} "${_libs}" PARENT_SCOPE)
endfunction()


set(PYOOMPH_MUMPS_SCALAPACK "" CACHE FILEPATH
    "The ScaLAPACK library MUMPS was built against. Only needed for a static MUMPS in an MPI build; found automatically when it has one of the usual names.")

# CMake ships no FindSCALAPACK either. Only the parallel MUMPS needs it, and only when MUMPS is
# static - a shared libdmumps.so records the dependency itself.
function(pyoomph_find_scalapack out_lib)
  if(PYOOMPH_MUMPS_SCALAPACK)
    set(${out_lib} "${PYOOMPH_MUMPS_SCALAPACK}" PARENT_SCOPE)
    return()
  endif()
  # scalapack-openmpi / scalapack-mpich: Debian and Ubuntu build one per MPI implementation.
  # "scalapack": Netlib, Fedora/RHEL, and what the reference build installs.
  find_library(PYOOMPH_MUMPS_SCALAPACK_LIBRARY
               NAMES scalapack scalapack-openmpi scalapack-mpich mkl_scalapack_lp64)
  if(PYOOMPH_MUMPS_SCALAPACK_LIBRARY)
    set(${out_lib} "${PYOOMPH_MUMPS_SCALAPACK_LIBRARY}" PARENT_SCOPE)
  else()
    set(${out_lib} "" PARENT_SCOPE)
  endif()
endfunction()


# A static MUMPS - however it was obtained - needs BLAS/LAPACK, and ScaLAPACK when it is parallel,
# after it on the link line. Factored out so the supplied-tree route gets the same treatment as the
# download one: without it a user pointing PYOOMPH_MUMPS_DIR at a static prefix gets "undefined
# symbol: ztrtrs_" at import time, which names neither MUMPS nor LAPACK.
function(pyoomph_append_static_mumps_deps out_libs)
  set(_libs ${${out_libs}})
  if(PYOOMPH_MUMPS_USE_MPI)
    pyoomph_find_scalapack(_sca)
    if(_sca)
      list(APPEND _libs "${_sca}")
      message(STATUS "pyoomph_mumps: ScaLAPACK: ${_sca}")
    else()
      message(WARNING
        "The MUMPS libraries are static and this is an MPI build, but no ScaLAPACK was found. If the "
        "module fails to import with an undefined pblas/blacs symbol, pass "
        "-DPYOOMPH_MUMPS_SCALAPACK=/path/to/libscalapack.so or add it to "
        "PYOOMPH_MUMPS_EXTRA_LIBRARIES.")
    endif()
  endif()
  find_package(LAPACK)
  if(LAPACK_FOUND)
    list(APPEND _libs ${LAPACK_LIBRARIES})
    message(STATUS "pyoomph_mumps: BLAS/LAPACK: ${LAPACK_LIBRARIES}")
  else()
    message(WARNING
      "The MUMPS libraries are static but no BLAS/LAPACK was found. If the module fails to import "
      "with an undefined symbol such as ztrtrs_, name one through PYOOMPH_MUMPS_EXTRA_LIBRARIES.")
  endif()
  set(${out_libs} "${_libs}" PARENT_SCOPE)
endfunction()


# ================================================================================================
# Route 1: a supplied tree
# ================================================================================================
if(PYOOMPH_MUMPS_DIR)
  pyoomph_find_mumps("${PYOOMPH_MUMPS_DIR}" _found _incs _libs _is_static)
  if(NOT _found)
    message(FATAL_ERROR
      "PYOOMPH_MUMPS_DIR=${PYOOMPH_MUMPS_DIR} does not contain a usable MUMPS: expected "
      "include/dmumps_c.h and lib/libdmumps plus libmumps_common. Point it at the install PREFIX "
      "(the directory holding include/ and lib/), not at either of them.")
  endif()
  set(PYOOMPH_MUMPS_INCLUDE_DIRS "${_incs}")
  set(PYOOMPH_MUMPS_LIBRARIES "${_libs}")
  set(PYOOMPH_MUMPS_NEEDS_FORTRAN ${_is_static})
  if(_is_static)
    pyoomph_append_static_mumps_deps(PYOOMPH_MUMPS_LIBRARIES)
  endif()
  message(STATUS "pyoomph_mumps: using the MUMPS in ${PYOOMPH_MUMPS_DIR}")

# ================================================================================================
# Route 2: download and build
# ================================================================================================
elseif(PYOOMPH_MUMPS_DOWNLOAD)
  if(NOT CMAKE_Fortran_COMPILER)
    message(FATAL_ERROR
      "Building MUMPS from source needs a Fortran compiler, and none was found. Install one "
      "(gfortran on Linux/macOS, or the Intel/LLVM Fortran compiler), or point "
      "-DPYOOMPH_MUMPS_DIR=... at a MUMPS somebody else has already built - a system package "
      "(libmumps-dev on Debian/Ubuntu) needs no Fortran compiler here, only its runtime.")
  endif()

  set(_pyoomph_mumps_prefix "${CMAKE_BINARY_DIR}/mumps")
  set(_pyoomph_mumps_install "${_pyoomph_mumps_prefix}/install")

  # --- the BLAS/LAPACK the superbuild will require ---
  set(_pyoomph_lapack_args "")
  set(_pyoomph_lapack_depends "")
  set(_pyoomph_need_lapack_download FALSE)
  if(PYOOMPH_MUMPS_DOWNLOAD_BLAS STREQUAL "ON")
    set(_pyoomph_need_lapack_download TRUE)
  elseif(PYOOMPH_MUMPS_DOWNLOAD_BLAS STREQUAL "AUTO")
    # Probe with the same call the superbuild will make, so the answer here predicts the answer
    # there. Not REQUIRED: a miss is the whole point of the probe.
    if(PYOOMPH_MUMPS_BLAS_ROOT)
      set(LAPACK_ROOT "${PYOOMPH_MUMPS_BLAS_ROOT}")
    endif()
    find_package(LAPACK QUIET)
    if(NOT LAPACK_FOUND)
      set(_pyoomph_need_lapack_download TRUE)
      message(STATUS "pyoomph_mumps: no BLAS/LAPACK found, a reference one will be downloaded "
                     "(-DPYOOMPH_MUMPS_DOWNLOAD_BLAS=OFF to fail instead)")
    endif()
  endif()

  if(_pyoomph_need_lapack_download)
    set(_pyoomph_lapack_install "${CMAKE_BINARY_DIR}/lapack/install")
    # The Netlib reference implementation, which is a real CMake project and builds both the
    # reference BLAS and LAPACK. It is SLOW compared with OpenBLAS or MKL and is meant as the
    # fallback that makes an unattended install work at all, not as the recommended configuration:
    # anyone who cares about factorisation speed should point PYOOMPH_MUMPS_BLAS_VENDOR/_ROOT at a
    # tuned BLAS, and the message says so.
    message(WARNING
      "pyoomph_mumps is building the REFERENCE (Netlib) BLAS/LAPACK. It is correct but slow, and "
      "MUMPS spends most of its time in BLAS3. For real runs configure with "
      "-DPYOOMPH_MUMPS_BLAS_VENDOR=OpenBLAS (or MKL/AOCL) and, if needed, "
      "-DPYOOMPH_MUMPS_BLAS_ROOT=/path/to/it.")
    ExternalProject_Add(pyoomph_lapack_external
      URL "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/${PYOOMPH_MUMPS_LAPACK_REF}.tar.gz"
      PREFIX "${CMAKE_BINARY_DIR}/lapack"
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${_pyoomph_lapack_install}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_TESTING=OFF
        -DCBLAS=OFF
        -DLAPACKE=OFF
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      BUILD_BYPRODUCTS
        "${_pyoomph_lapack_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${_pyoomph_lapack_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    list(APPEND _pyoomph_lapack_depends pyoomph_lapack_external)
    list(APPEND _pyoomph_lapack_args
         "-DLAPACK_ROOT=${_pyoomph_lapack_install}"
         "-DCMAKE_PREFIX_PATH=${_pyoomph_lapack_install}"
         "-DMUMPS_find_static=ON")
  else()
    if(PYOOMPH_MUMPS_BLAS_VENDOR)
      list(APPEND _pyoomph_lapack_args "-DLAPACK_VENDOR=${PYOOMPH_MUMPS_BLAS_VENDOR}")
    endif()
    if(PYOOMPH_MUMPS_BLAS_ROOT)
      list(APPEND _pyoomph_lapack_args "-DLAPACK_ROOT=${PYOOMPH_MUMPS_BLAS_ROOT}")
    endif()
  endif()

  if(PYOOMPH_MUMPS_SUPERBUILD_URL)
    set(_pyoomph_mumps_url "${PYOOMPH_MUMPS_SUPERBUILD_URL}")
  else()
    set(_pyoomph_mumps_url "https://github.com/scivision/mumps/archive/${PYOOMPH_MUMPS_SUPERBUILD_REF}.tar.gz")
  endif()

  set(_pyoomph_mumps_version_arg "")
  if(PYOOMPH_MUMPS_UPSTREAM_VERSION)
    set(_pyoomph_mumps_version_arg "-DMUMPS_UPSTREAM_VERSION=${PYOOMPH_MUMPS_UPSTREAM_VERSION}")
  endif()

  # OpenMP inside MUMPS is a separate question from OpenMP in this extension's own (nonexistent)
  # loops: it is what makes a single-rank factorisation use more than one core, and it is what
  # ICNTL(16) then controls at run time.
  if(PYOOMPH_MUMPS_USE_OPENMP)
    set(_pyoomph_mumps_omp ON)
  else()
    set(_pyoomph_mumps_omp OFF)
  endif()

  # Static (BUILD_SHARED_LIBS=OFF, the superbuild's own default). A shared MUMPS in a private prefix
  # would need either an RPATH into the build tree or an LD_LIBRARY_PATH the user has to set, and
  # this prefix is inside the build directory rather than anywhere installed. Static libraries get
  # linked into the extension module and the question disappears.
  ExternalProject_Add(pyoomph_mumps_external
    URL "${_pyoomph_mumps_url}"
    PREFIX "${_pyoomph_mumps_prefix}"
    DEPENDS ${_pyoomph_lapack_depends}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${_pyoomph_mumps_install}
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
      -DBUILD_SHARED_LIBS=OFF
      # Only the two arithmetics this package binds. Building s/c as well roughly doubles the
      # compile time for code nothing here can reach.
      -DBUILD_SINGLE=OFF
      -DBUILD_DOUBLE=ON
      -DBUILD_COMPLEX=OFF
      -DBUILD_COMPLEX16=ON
      -DMUMPS_parallel=${PYOOMPH_MUMPS_USE_MPI}
      -DMUMPS_scalapack=${PYOOMPH_MUMPS_USE_MPI}
      -DMUMPS_openmp=${_pyoomph_mumps_omp}
      # 32-bit integers, matching oomph-lib's own CRDoubleMatrix index type. pyoomph hands the solver
      # int32 index arrays; a 64-bit-integer MUMPS would read them as half as many int64s and produce
      # nonsense rather than an error.
      -DMUMPS_intsize64=OFF
      -DMUMPS_BUILD_TESTING=OFF
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      ${_pyoomph_mumps_version_arg}
      ${_pyoomph_lapack_args}
    BUILD_BYPRODUCTS
      "${_pyoomph_mumps_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}dmumps${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${_pyoomph_mumps_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}zmumps${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${_pyoomph_mumps_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}mumps_common${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${_pyoomph_mumps_install}/include/dmumps_c.h"
  )
  list(APPEND PYOOMPH_MUMPS_DEPENDS pyoomph_mumps_external)

  # The libraries do not exist yet at configure time, so they are named rather than found. The order
  # is the link order a static build needs: the arithmetic libraries call into mumps_common, which
  # calls into pord.
  set(_l "${_pyoomph_mumps_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}")
  set(_s "${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(PYOOMPH_MUMPS_INCLUDE_DIRS "${_pyoomph_mumps_install}/include")
  set(PYOOMPH_MUMPS_LIBRARIES "${_l}dmumps${_s}" "${_l}zmumps${_s}" "${_l}mumps_common${_s}" "${_l}pord${_s}")
  if(NOT PYOOMPH_MUMPS_USE_MPI)
    list(APPEND PYOOMPH_MUMPS_LIBRARIES "${_l}mpiseq${_s}")
  endif()

  # What a STATIC MUMPS needs after itself. A static archive contributes only the members that
  # resolve symbols already referenced, so these have to FOLLOW the MUMPS libraries on the link
  # line - and none of them is guesswork: MUMPS's own INSTALL file names exactly this set (BLAS,
  # LAPACK, ScaLAPACK+BLACS when parallel, and the MPI bindings). It is not open-ended.
  #
  # The MPI Fortran bindings and the Fortran runtime are appended by the top-level CMakeLists.txt
  # instead, because they come from targets it owns (MPI::MPI_Fortran,
  # CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES).
  #
  # Nothing here fails the LINK when it is missing: a Python extension module is linked with
  # undefined symbols allowed, so an omission surfaces as "undefined symbol: ztrtrs_" at import
  # time. The stub-generation step imports the module and is what turns that back into a build
  # failure.
  if(PYOOMPH_MUMPS_USE_MPI)
    pyoomph_find_scalapack(_pyoomph_scalapack)
    if(NOT _pyoomph_scalapack)
      message(FATAL_ERROR
        "A parallel MUMPS is being built but no ScaLAPACK was found, and a static MUMPS cannot be "
        "linked without one. Install it (libscalapack-openmpi-dev on Debian/Ubuntu, "
        "scalapack-openmpi-devel on Fedora), point -DPYOOMPH_MUMPS_SCALAPACK=/path/to/libscalapack.so "
        "at it, or build without MPI - which then requires a pyoomph built without MPI too.")
    endif()
    list(APPEND PYOOMPH_MUMPS_LIBRARIES "${_pyoomph_scalapack}")
    message(STATUS "pyoomph_mumps: ScaLAPACK: ${_pyoomph_scalapack}")
  endif()

  if(_pyoomph_need_lapack_download)
    # Named rather than found, like the MUMPS libraries above: they do not exist yet at configure
    # time. LAPACK before BLAS - LAPACK calls into it.
    set(_ll "${_pyoomph_lapack_install}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}")
    list(APPEND PYOOMPH_MUMPS_LIBRARIES "${_ll}lapack${_s}" "${_ll}blas${_s}")
  else()
    # find_package(LAPACK) was already run by the AUTO probe above; re-running it is a no-op that
    # just makes this branch readable on its own.
    find_package(LAPACK REQUIRED)
    list(APPEND PYOOMPH_MUMPS_LIBRARIES ${LAPACK_LIBRARIES})
    message(STATUS "pyoomph_mumps: BLAS/LAPACK: ${LAPACK_LIBRARIES}")
  endif()

  set(PYOOMPH_MUMPS_NEEDS_FORTRAN TRUE)
  message(STATUS "pyoomph_mumps: downloading and building MUMPS (superbuild ${PYOOMPH_MUMPS_SUPERBUILD_REF})")

# ================================================================================================
# Route 3: a system MUMPS on the default paths
# ================================================================================================
else()
  pyoomph_find_mumps("" _found _incs _libs _is_static)
  if(NOT _found)
    message(FATAL_ERROR
      "No MUMPS found, and both PYOOMPH_MUMPS_DIR and PYOOMPH_MUMPS_DOWNLOAD are off. Either "
      "install one (libmumps-dev on Debian/Ubuntu, mumps-devel on Fedora, brew install mumps), "
      "point -DPYOOMPH_MUMPS_DIR=... at a prefix, or configure with -DPYOOMPH_MUMPS_DOWNLOAD=ON.")
  endif()
  set(PYOOMPH_MUMPS_INCLUDE_DIRS "${_incs}")
  set(PYOOMPH_MUMPS_LIBRARIES "${_libs}")
  set(PYOOMPH_MUMPS_NEEDS_FORTRAN ${_is_static})
  if(_is_static)
    pyoomph_append_static_mumps_deps(PYOOMPH_MUMPS_LIBRARIES)
  endif()
  message(STATUS "pyoomph_mumps: using the system MUMPS (${PYOOMPH_MUMPS_D_LIBRARY})")
endif()

if(PYOOMPH_MUMPS_EXTRA_LIBRARIES)
  list(APPEND PYOOMPH_MUMPS_LIBRARIES ${PYOOMPH_MUMPS_EXTRA_LIBRARIES})
endif()
