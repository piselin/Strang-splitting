# Module to find the FFTW library.
# following variables will be defined:
# FFTW_FOUND            if the library was found
# FFTW_LIBRARIES        path to all libraries
# FFTW_INCLUDE_DIRS     
# FFTW_LIB              default fftw3 library (double-precision)
# FFTW_LIB_FLOAT        path to single precision library
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PKG_FFTW QUIET "fftw3")
    find_library(FFTW_LIB NAMES "fftw3")
    find_library(FFTW_LIB_FLOAT NAMES "fftw3f")
    find_path(FFTW_INCLUDE_DIRS NAMES "fftw3.h")
endif()

set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTW_LIB_FLOAT})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIBRARIES FFTW_LIB FFTW_LIB_FLOAT)