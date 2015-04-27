# - Try to find LibMagma
# Once done this will define
#  LIBMAGMA_FOUND - System has LibMagma
#  LIBMAGMA_INCLUDE_DIRS - The LibMagma include directories
#  LIBMAGMA_LIBRARIES - The libraries needed to use LibMagma
#  LIBMAGMA_DEFINITIONS - Compiler switches required for using LibMagma

find_package(PkgConfig)
pkg_check_modules(PC_LIBMAGMA QUIET magma)
set(LIBMAGMA_DEFINITIONS ${PC_LIBMAGMA_CFLAGS_OTHER})

find_path(LIBMAGMA_INCLUDE_DIR magma.h
          HINTS ${PC_LIBMAGMA_INCLUDE_DIR} ${PC_LIBMAGMA_INCLUDE_DIRS})

find_library(LIBMAGMA_LIBRARY NAMES magma libmagma
             HINTS ${PC_LIBMAGMA_LIBDIR} ${PC_LIBMAGMA_LIBRARY_DIRS})

set(LIBMAGMA_LIBRARIES ${LIBMAGMA_LIBRARY})
set(LIBMAGMA_INCLUDE_DIRS ${LIBMAGMA_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBMAGMA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibMagma DEFAULT_MSG
                                  LIBMAGMA_LIBRARY LIBMAGMA_INCLUDE_DIR)

mark_as_advanced(LIBMAGMA_INCLUDE_DIR LIBMAGMA_LIBRARY)
