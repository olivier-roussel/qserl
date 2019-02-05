FIND_PACKAGE(PythonLibs REQUIRED 2.7)

ADD_LIBRARY(rod3d SHARED rod3d.cc)
SET_TARGET_PROPERTIES(rod3d PROPERTIES
  PREFIX "")

include_directories ( ${EIGEN3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}/unsupported
  )

if(${CMAKE_VERSION} VERSION_GREATER 3.8)
  target_compile_features(rod3d PRIVATE cxx_std_11)
endif()

target_link_libraries(rod3d
  qserl Boost::boost Boost::python ${PYTHON_LIBRARIES})
PKG_CONFIG_USE_DEPENDENCY(rod3d "eigenpy")

install(TARGETS rod3d
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/python2.7/site-packages/qserl
  )
install(FILES __init__.py
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/python2.7/site-packages/qserl
  )
