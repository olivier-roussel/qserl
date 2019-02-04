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
# Setting up target
#------------------------------------------------------------------------------

macro(add_example name)
  add_executable(${name} ${name}.cpp)
  target_compile_features(${name} PRIVATE cxx_auto_type)
  target_link_libraries(${name}
    PRIVATE
    qserl)
endmacro(add_example)

add_example(example_rod_2D)
add_example(example_rod_3D)
add_example(generate_inverse_geometry_benchmark_dataset_2D)
add_example(generate_stability_dataset_2D)