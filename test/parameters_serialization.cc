/**
* Copyright (c) 2012-2018 CNRS
* Author: Olivier Roussel
*
* This file is part of the qserl package.
* qserl is free software: you can redistribute it
* and/or modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation, either version
* 3 of the License, or (at your option) any later version.
*
* qserl is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Lesser Public License for more details.  You should have
* received a copy of the GNU Lesser General Public License along with
* qserl.  If not, see
* <http://www.gnu.org/licenses/>.
**/

#include <boost/test/unit_test.hpp>

#include "qserl/rod2d/parameters.h"
#include "qserl/rod3d/parameters.h"

#include <fstream>

/* ------------------------------------------------------------------------- */
/* Parameters_serialization     																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(Parameters_serialization)

BOOST_AUTO_TEST_CASE(Rod2d_parameters_serialization)
{
  auto params = qserl::rod2d::Parameters{};
  params.integrationTime = 0.5;
  params.radius = 0.2;
  params.length = 3.;
  params.delta_t = 1.e-4;

  auto params_serialized = qserl::rod2d::Parameters{};
  BOOST_CHECK(params.integrationTime != params_serialized.integrationTime);

  const auto filename = "test_rod2d_parameters_serialization.xml";
  {
    auto ofs = std::ofstream{filename};
    cereal::XMLOutputArchive ar(ofs);
    ar(cereal::make_nvp("params", params));
  }
  {
    auto ifs = std::ifstream{filename};
    cereal::XMLInputArchive ar(ifs);
    ar( cereal::make_nvp("params", params_serialized) );
  }
  BOOST_CHECK(params.radius == params_serialized.radius);
  BOOST_CHECK(params.length == params_serialized.length);
  BOOST_CHECK(params.integrationTime == params_serialized.integrationTime);
  BOOST_CHECK(params.delta_t == params_serialized.delta_t);
  BOOST_CHECK(params.rodModel == params_serialized.rodModel);
}

BOOST_AUTO_TEST_CASE(Rod3d_parameters_serialization)
{
  auto params = qserl::rod3d::Parameters{};
  params.radius = 0.2;
  params.length = 3.;
  params.stiffnessCoefficients[0] = 10.;
  params.stiffnessCoefficients[1] = 11.;
  params.stiffnessCoefficients[2] = 12.;
  params.stiffnessCoefficients[3] = 13.;
  params.stiffnessCoefficients[4] = 14.;
  params.stiffnessCoefficients[5] = 15.;
  params.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  params.numNodes = 100;
  params.integrationTime = 0.5;

  auto params_serialized = qserl::rod3d::Parameters{};
  BOOST_CHECK(params.integrationTime != params_serialized.integrationTime);

  const auto filename = "test_rod3d_parameters_serialization.xml";
  {
    auto ofs = std::ofstream{filename};
    cereal::XMLOutputArchive ar(ofs);
    ar(cereal::make_nvp("params", params));
  }
  {
    auto ifs = std::ifstream{filename};
    cereal::XMLInputArchive ar(ifs);
    ar( cereal::make_nvp("params", params_serialized) );
  }
  BOOST_CHECK(params.radius == params_serialized.radius);
  BOOST_CHECK(params.length == params_serialized.length);
  BOOST_CHECK((params.stiffnessCoefficients - params_serialized.stiffnessCoefficients).norm() < 1.e-6);
  BOOST_CHECK(params.rodModel == params_serialized.rodModel);
  BOOST_CHECK(params.numNodes == params_serialized.numNodes);
  BOOST_CHECK(params.integrationTime == params_serialized.integrationTime);
}

BOOST_AUTO_TEST_SUITE_END();
