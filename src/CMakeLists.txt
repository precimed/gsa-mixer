cmake_minimum_required(VERSION 2.8)

project(bgmg)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")
message("-- Build type: ${CMAKE_BUILD_TYPE}")
message("-- C++ compiler: ${CMAKE_CXX_COMPILER_ID}")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else (COMPILER_SUPPORTS_CXX11)
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif (COMPILER_SUPPORTS_CXX11)

if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -march=native -mfpmath=sse -ggdb")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -march=native -mfpmath=sse -ggdb")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
      ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
      ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))

#SET(BUILD_GMOCK OFF CACHE BOOL "Builds the googlemock subproject")
#SET(BUILD_GTEST ON CACHE BOOL "Builds the googletest subproject")
#add_subdirectory( googletest )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(zlib)

message("CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

set(BGMG_BOOST_COMPONENTS program_options filesystem system date_time)
find_package(Boost COMPONENTS REQUIRED ${BGMG_BOOST_COMPONENTS})
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find required boost libraries.")
  return()
endif (NOT Boost_FOUND)
message("-- Boost_LIBRARIES: ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})

set(SRC_LIST
	bgmg_calculator.cc
	bgmg_calculator_unified.cc
	bgmg_calculator_legacy.cc
	ld_matrix_csr.cc
	ld_matrix.cc
	bgmg.h
	bgmg_calculator.h
	ld_matrix_csr.h
	annotate_intervals.cc
	annotate_intervals.h
	bgmg_log.cc
	bgmg_log.h
	bgmg_parse.cc
	bgmg_parse.h
	fmath.hpp
	plink_ld.cc
	plink_common.cc
	semt/semt/VectorExpr.cpp
	snp_lookup.cc
	TurboPFor/vsimple.c
	TurboPFor/vint.c
	TurboPFor/bitutil.c
	nlopt/nldrmd.c
	nlopt/sbplx.c
	nlopt/redblack.c
	nlopt/stop.c
	nlopt/qsort_r.c
	nlopt/timer.c
	cubature/hcubature.c
	FastDifferentialCoding/fastdelta.c
	SIMDxorshift/xorshift128plus.c
	SIMDxorshift/simdxorshift128plus.c	
)

add_library(bgmg-static STATIC ${SRC_LIST})
set_target_properties(bgmg-static PROPERTIES LINK_FLAGS "-static" )
target_include_directories(bgmg-static PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(bgmg-static PUBLIC ${CMAKE_BINARY_DIR}/zlib)
target_include_directories(bgmg-static PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(bgmg-static PUBLIC ${PROJECT_SOURCE_DIR}/semt)
target_include_directories(bgmg-static PUBLIC ${PROJECT_SOURCE_DIR}/zlib)
target_compile_definitions(bgmg-static PRIVATE USE_SSE)

add_library(bgmg SHARED source.cc)
target_link_libraries(bgmg bgmg-static zlibstatic ${Boost_LIBRARIES})

set(SRC_TEST_LIST
	source.cc
	test/bgmg_main_unittest.cc
	test/bgmg_ld_test.cc
	test/semt_test.cc
	test/streamvbyte_test.cc
	test/SIMDxorshift_test.cc
	test/nlopt_test.cc
	test/cubature_test.cc
	test/bgmg_test_main.cc
	test/bgmg_math_test.cc
	test/annotate_intervals_test.cc
	${PROJECT_SOURCE_DIR}/googletest/googletest/src/gtest-all.cc
)
add_executable(bgmg-test ${SRC_TEST_LIST})
target_link_libraries(bgmg-test bgmg-static zlibstatic ${Boost_LIBRARIES}) 
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest/include)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/googletest/googletest)
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(bgmg-test PUBLIC ${PROJECT_SOURCE_DIR}/semt)
target_compile_definitions(bgmg-test PRIVATE BGMG_STATIC_DEFINE)

include(GenerateExportHeader)
generate_export_header(bgmg
		EXPORT_MACRO_NAME DLL_PUBLIC)

install( TARGETS bgmg DESTINATION bin )
