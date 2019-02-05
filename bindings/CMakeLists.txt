add_library(rod3d SHARED
  rod3d.cc)
set_target_properties(rod3d PROPERTIES
  PREFIX "")

# XXX macro SEARCH_FOR_EIGEN has normally included Eigen directories
# for all the project
target_include_directories(rod3d
  PRIVATE
#    ${EIGEN3_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}/unsupported
  )

if(${CMAKE_VERSION} VERSION_GREATER 3.8)
  target_compile_features(rod3d PRIVATE cxx_std_11)
endif()

target_link_libraries(rod3d
  qserl
  Boost::python
  )

PKG_CONFIG_USE_DEPENDENCY(rod3d "eigenpy")

install(TARGETS rod3d
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/python2.7/site-packages/qserl
  )
install(FILES __init__.py
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/python2.7/site-packages/qserl
  )
