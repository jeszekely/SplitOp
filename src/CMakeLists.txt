# CMakeLists.txt in src
project(SplitOp)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_COMPILER /opt/local/bin/g++)
set(CMAKE_C_COMPILER /opt/local/bin/gcc)
add_compile_options(-O3 -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw -std=c++11 -Wall -Wno-sign-compare -Wno-unused-function -Werror)

set(SplitOp_SOURCES
        chebyshev.cpp
        input_parser.cpp
        junction.cpp
        matrix.cpp
        splitop.cpp
        surrogate.cpp
        vector.cpp
        wvfxn.cpp
        main.cpp)

# Find and link boost libraries
find_package(boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
#link_directories(${Boost_LIBRARY_DIR})

# BLAS and LAPACK
find_package(LAPACK REQUIRED)

# FFTW
if (FFTW_INCLUDES)
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)
find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTW_LIBRARIES NAMES fftw3)
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

add_executable(SplitOp ${SplitOp_SOURCES})

#linking
target_link_libraries(SplitOp mkl_intel_lp64 mkl_intel_thread mkl_core mkl_sequential pthread m ${BOOST_LIBRARIES})	#${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# # MKL
# find_package(MKL REQUIRED)
# include_directories(${MKL_INCLUDE_DIRS})
# link_directories(${MKL_LIBRARIES})
# target_link_libraries(<module>
# mkl_intel_lp64
# mkl_intel_thread
# mkl_sequential
# mkl_core
# )
