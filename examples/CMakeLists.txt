# -------------------------------------
# Dependencies
# -------------------------------------
include_directories("${PROJECT_SOURCE_DIR}/src")

# Use Boost
if(MSVC)
	add_definitions(-DBOOST_ALL_DYN_LINK)
endif(MSVC)
if(IS_ICPC)
    set(Boost_USE_STATIC_LIBS ON CACHE STRING "Use statically linked Boost libraries")
else(IS_ICPC)
    # Ensure dynamic linking with boost unit_test_framework
    add_definitions(-DBOOST_TEST_DYN_LINK)
endif(IS_ICPC)
find_package(Boost 1.55 COMPONENTS chrono date_time filesystem program_options regex serialization system thread unit_test_framework REQUIRED)
if (Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIR}")
  set(LIBS ${LIBS} ${Boost_LIBRARIES}) 
else()
  message("Could not find Boost")
endif()


# Use Eigen3 (with LGSM module)
if (NOT $ENV{EIGEN3_INCLUDE_DIR} STREQUAL "")
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
endif()
find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message("Could not find Eigen3")
endif()
# Set EIGEN_DONT_ALIGN definition for 32bits platforms
# XXX crashes if not set also under MSVC10 Win64 _ 
#if(NOT CMAKE_CL_64)
	add_definitions(-DEIGEN_DONT_ALIGN)
#endif()



macro(add_example name)
add_executable(${name} ${name}.cpp)
target_link_libraries(${name} qserl ${Boost_LIBRARIES} )
add_dependencies(${name} qserl)
add_test(${name} ${name})
endmacro(add_example)
add_example(example0)


