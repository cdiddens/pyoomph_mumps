cmake_minimum_required(VERSION 3.22)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Declared BEFORE project() because whether a Fortran compiler is needed depends on them, and
# enable_language() cannot be undone once project() has picked the languages. Re-declaring them in
# cmake/ThirdPartyMUMPS.cmake is harmless: a second set(... CACHE ...) without FORCE is a no-op.
set(PYOOMPH_MUMPS_DIR "" CACHE PATH
    "An already installed MUMPS tree to use instead of downloading one.")
option(PYOOMPH_MUMPS_DOWNLOAD "Download and build MUMPS from source" ON)

if(NOT DEFINED SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")
endif()

project(pyoomph_mumps LANGUAGES C CXX VERSION ${SKBUILD_PROJECT_VERSION})

# Only the download route compiles Fortran. A build against a supplied or system MUMPS needs the
# Fortran RUNTIME (which a shared libdmumps.so already records as a dependency of its own), not a
# compiler - so requiring one there would refuse perfectly good configurations.
if(PYOOMPH_MUMPS_DOWNLOAD AND NOT PYOOMPH_MUMPS_DIR)
  enable_language(Fortran)
endif()

# ----------------------------------------------------------------------------
# Build options
#
# Set on the cmake command line (-DPYOOMPH_MUMPS_...=...) or, through pip/scikit-build-core as
# ./build_for_develop.sh does, via
# `pip install --config-settings=cmake.define.PYOOMPH_MUMPS_...=...`.
# ----------------------------------------------------------------------------

# MUST match how pyoomph itself was built. A pyoomph with MPI and a serial pyoomph_mumps (or the
# reverse) is the failure mode that costs an afternoon: MUMPS's serial stub library defines its own
# MPI_Init and friends, so the two either initialise MPI twice or not at all, and what comes out is a
# hang or a wrong communicator rather than an error. pyoomph_mumps/__init__.py compares this against
# pyoomph's own answer at import time and refuses the pairing loudly - see the note there.
option(PYOOMPH_MUMPS_USE_MPI "Build with MPI support (must match how pyoomph itself was built)" OFF)

# Tri-state for the same reason pyoomph's PYOOMPH_USE_OPENMP is one: with AUTO alone, whether a build
# can use more than one core per rank is decided silently by whatever the build machine happened to
# have. ON says the intent out loud and fails at configure time when it cannot be met.
set(PYOOMPH_MUMPS_USE_OPENMP "AUTO" CACHE STRING
    "OpenMP inside MUMPS (multithreaded factorisation, controlled at run time by ICNTL(16)): AUTO=use if found, ON=require, OFF=never")
set_property(CACHE PYOOMPH_MUMPS_USE_OPENMP PROPERTY STRINGS AUTO ON OFF)

# Worth keeping on for more than the stub: nanobind.stubgen IMPORTS the freshly built module, which
# is the only step that notices an extension that linked cleanly but cannot load - the normal state
# of affairs for a Python extension, whose undefined symbols are allowed at link time. A missing
# Fortran MPI library on the static download route showed up exactly here and nowhere else.
option(PYOOMPH_MUMPS_GENERATE_STUBS "Generate the .pyi stub for the extension via nanobind.stubgen (also acts as an import check)" ON)
option(PYOOMPH_MUMPS_REQUIRE_STUBS "Fail the build when the .pyi stub cannot be generated, rather than warning" ON)
option(PYOOMPH_MUMPS_COPY_STUBS_TO_SOURCE_TREE "Also mirror the generated .pyi into the source-tree pyoomph_mumps/ directory, so Pylance/Pyright/mypy resolve it while editing without a full install" ON)

# ----------------------------------------------------------------------------
# OpenMP
#
# Resolved BEFORE the MUMPS acquisition, which reads PYOOMPH_MUMPS_USE_OPENMP to decide whether to
# build MUMPS itself with -DMUMPS_openmp=ON. Note what the flag does and does not mean: this
# extension has no parallel loops of its own, and links OpenMP only so that a MUMPS built against it
# resolves its symbols. All the threading is inside MUMPS.
# ----------------------------------------------------------------------------
if(PYOOMPH_MUMPS_USE_OPENMP)
  if(PYOOMPH_MUMPS_USE_OPENMP STREQUAL "AUTO")
    find_package(OpenMP COMPONENTS CXX)
  else()
    find_package(OpenMP REQUIRED COMPONENTS CXX)
  endif()
endif()
if(OpenMP_CXX_FOUND AND PYOOMPH_MUMPS_USE_OPENMP)
  set(PYOOMPH_MUMPS_HAS_OPENMP TRUE)
else()
  set(PYOOMPH_MUMPS_HAS_OPENMP FALSE)
endif()

# ----------------------------------------------------------------------------
# MPI
# ----------------------------------------------------------------------------
if(PYOOMPH_MUMPS_USE_MPI)
  set(_pyoomph_mumps_mpi_components C CXX)
  if(CMAKE_Fortran_COMPILER)
    # MUMPS is Fortran and calls the FORTRAN MPI bindings (mpi_type_free_, mpi_allreduce_, ...),
    # which live in a different library from the C ones - libmpi_mpifh with Open MPI, libmpifort with
    # MPICH. A shared libdmumps.so records that dependency itself, so the system route never needs
    # this; a static one, which is what the download route builds, records nothing.
    #
    # And a Python extension module is linked with undefined symbols ALLOWED (it has to be - the
    # CPython symbols are resolved by the interpreter), so leaving this out does not fail the link.
    # It fails at `import`, with "undefined symbol: mpi_type_free_" and no mention of Fortran.
    list(APPEND _pyoomph_mumps_mpi_components Fortran)
  endif()
  find_package(MPI REQUIRED COMPONENTS ${_pyoomph_mumps_mpi_components})
endif()

# ----------------------------------------------------------------------------
# MUMPS
# ----------------------------------------------------------------------------
include(cmake/ThirdPartyMUMPS.cmake)

# ----------------------------------------------------------------------------
# Python / nanobind
# ----------------------------------------------------------------------------
find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED
                        OPTIONAL_COMPONENTS Development.SABIModule)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT
  COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind CONFIG REQUIRED)

# Its OWN nanobind domain, not pyoomph's. The domain isolates the type registry, and the only things
# crossing this boundary are numpy arrays and plain scalars - never a bound C++ type - so sharing
# pyoomph's domain would buy nothing and would tie this module's ABI to pyoomph's.
nanobind_add_module(_pyoomph_mumps_core
  NB_DOMAIN pyoomph_mumps
  STABLE_ABI
  src/bindings.cpp
  src/mumps_solver.cpp)

target_compile_features(_pyoomph_mumps_core PRIVATE cxx_std_17)
target_include_directories(_pyoomph_mumps_core PRIVATE ${PYOOMPH_MUMPS_INCLUDE_DIRS})
target_link_libraries(_pyoomph_mumps_core PRIVATE ${PYOOMPH_MUMPS_LIBRARIES})

if(PYOOMPH_MUMPS_DEPENDS)
  add_dependencies(_pyoomph_mumps_core ${PYOOMPH_MUMPS_DEPENDS})
endif()

if(PYOOMPH_MUMPS_USE_MPI)
  target_compile_definitions(_pyoomph_mumps_core PRIVATE PYOOMPH_MUMPS_HAS_MPI)
  # After the MUMPS libraries, not before: a static archive only contributes the members that resolve
  # symbols already referenced, so the library that satisfies MUMPS has to follow it on the link line.
  target_link_libraries(_pyoomph_mumps_core PRIVATE MPI::MPI_C MPI::MPI_CXX)
  if(TARGET MPI::MPI_Fortran)
    target_link_libraries(_pyoomph_mumps_core PRIVATE MPI::MPI_Fortran)
  endif()
endif()

if(PYOOMPH_MUMPS_HAS_OPENMP)
  target_compile_definitions(_pyoomph_mumps_core PRIVATE PYOOMPH_MUMPS_HAS_OPENMP)
  target_link_libraries(_pyoomph_mumps_core PRIVATE OpenMP::OpenMP_CXX)
endif()

# A static MUMPS carries no record of its own dependencies, so the Fortran runtime (libgfortran,
# libquadmath, ...) has to be named explicitly. CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is what the
# Fortran compiler itself would have added, which is why enable_language(Fortran) above is not just
# about compiling: it is also how this list gets populated.
if(PYOOMPH_MUMPS_NEEDS_FORTRAN)
  if(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
    target_link_libraries(_pyoomph_mumps_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    target_link_directories(_pyoomph_mumps_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  else()
    message(WARNING
      "The MUMPS libraries are static but no Fortran runtime is known to this build "
      "(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is empty). If the link fails with undefined "
      "_gfortran_* symbols, add the runtime through -DPYOOMPH_MUMPS_EXTRA_LIBRARIES=gfortran.")
  endif()
endif()

# What the Python layer reads back to describe the build it is talking to.
target_compile_definitions(_pyoomph_mumps_core PRIVATE
  PYOOMPH_MUMPS_BUILD_VERSION="${SKBUILD_PROJECT_VERSION}")

install(TARGETS _pyoomph_mumps_core LIBRARY DESTINATION pyoomph_mumps)

# ----------------------------------------------------------------------------
# .pyi stub
#
# Same script, same three options and the same rationale as pyoomph's: the package ships py.typed, so
# a build without the stub tells type checkers the API is described and then hands them an
# unintrospectable binary.
# ----------------------------------------------------------------------------
if(PYOOMPH_MUMPS_GENERATE_STUBS)
  set(PYOOMPH_MUMPS_STUB_EXTRA_ARGS "")
  if(PYOOMPH_MUMPS_COPY_STUBS_TO_SOURCE_TREE)
    list(APPEND PYOOMPH_MUMPS_STUB_EXTRA_ARGS --extra-copy-dir "${CMAKE_SOURCE_DIR}/pyoomph_mumps")
  endif()
  if(PYOOMPH_MUMPS_REQUIRE_STUBS)
    list(APPEND PYOOMPH_MUMPS_STUB_EXTRA_ARGS --required)
  endif()
  add_custom_command(TARGET _pyoomph_mumps_core POST_BUILD
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/cmake/generate_stubs.py"
            --module-dir "$<TARGET_FILE_DIR:_pyoomph_mumps_core>"
            --module-name "_pyoomph_mumps_core"
            --stage-dir "${CMAKE_BINARY_DIR}/pyi-stage"
            --python "${Python_EXECUTABLE}"
            ${PYOOMPH_MUMPS_STUB_EXTRA_ARGS}
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}" VERBATIM)
  if(PYOOMPH_MUMPS_REQUIRE_STUBS)
    install(DIRECTORY "${CMAKE_BINARY_DIR}/pyi-stage/" DESTINATION pyoomph_mumps)
  else()
    install(DIRECTORY "${CMAKE_BINARY_DIR}/pyi-stage/" DESTINATION pyoomph_mumps OPTIONAL)
  endif()
endif()

# ----------------------------------------------------------------------------
message(STATUS "pyoomph_mumps: MPI=${PYOOMPH_MUMPS_USE_MPI}  OpenMP=${PYOOMPH_MUMPS_HAS_OPENMP}")
message(STATUS "pyoomph_mumps: MUMPS libraries: ${PYOOMPH_MUMPS_LIBRARIES}")
