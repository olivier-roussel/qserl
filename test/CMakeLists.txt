# Copyright (c) 2012-2018 CNRS
# Author: Olivier Roussel
#
# This file is part of the qserl package.
# qserl is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# qserl is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# qserl.  If not, see
# <http://www.gnu.org/licenses/>.

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------

find_package(Boost 1.55 REQUIRED MODULE COMPONENTS unit_test_framework)

#------------------------------------------------------------------------------
# Setting up target
#------------------------------------------------------------------------------

add_executable(qserl-tests
  main.cc
  rod2d_analytic_vs_numeric_dqda.cc
  rod2d_analytic_vs_numeric_mu.cc
  rod2d_analytic_vs_numeric_q.cc
  rod2d_integrated_tests.cc
  rod3d_integrated_tests.cc
  )

target_include_directories(qserl-tests
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src
  )

target_compile_features(qserl-tests PRIVATE cxx_auto_type)

target_link_libraries(qserl-tests
    PRIVATE
      qserl
      Boost::unit_test_framework
    )

include(CMakeUnitTests.txt)

add_custom_command(TARGET qserl-tests POST_BUILD
  COMMAND qserl-tests
    --log_level=test_suite
    --result_code=no
    --report_level=confirm
    --detect_memory_leak=1
  )

