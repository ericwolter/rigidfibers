# - Try to find LibOpenBlas
# Once done this will define
#  LIBOPENBLAS_FOUND - System has LibOpenBlas
#  LIBOPENBLAS_INCLUDE_DIRS - The LibOpenBlas include directories
#  LIBOPENBLAS_LIBRARIES - The libraries needed to use LibOpenBlas

find_path(LIBOPENBLAS_INCLUDE_DIR openblas_config.h
          HINTS $ENV{OPENBLAS_ROOT}/include)

find_library(LIBOPENBLAS_LIBRARY NAMES openblas libopenblas
             HINTS $ENV{OPENBLAS_ROOT}/lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBOPENBLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibOpenBlas DEFAULT_MSG
                                  LIBOPENBLAS_LIBRARY LIBOPENBLAS_INCLUDE_DIR)

mark_as_advanced(LIBOPENBLAS_INCLUDE_DIR LIBOPENBLAS_LIBRARY)
