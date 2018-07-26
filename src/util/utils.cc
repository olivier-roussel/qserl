/**
* Copyright (c) 2012-2014 CNRS
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

#include "utils.h"

namespace qserl {
namespace util {

/** Do not use this one directly. */
std::vector<std::string> &split(const std::string &s,
                                char delim,
                                std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s,
                               char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

double rand_interval(double rmin,
                     double rmax)
{
  double t = rand() / ((double)RAND_MAX + 1);
  return (t * (rmax - rmin) + rmin);
}

Eigen::Vector4f getJetColor(double v,
                            double vmin,
                            double vmax)
{
   Eigen::Vector4f c(1.0f, 1.0f, 1.0f, 1.f); 
   double dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      c[0] = 0.;
      c[1] = static_cast<float>(4. * (v - vmin) / dv);
   } else if (v < (vmin + 0.5 * dv)) {
      c[0] = 0.;
      c[2] = static_cast<float>(1. + 4. * (vmin + 0.25 * dv - v) / dv);
   } else if (v < (vmin + 0.75 * dv)) {
      c[0] = static_cast<float>(4. * (v - vmin - 0.5 * dv) / dv);
      c[2] = 0.;
   } else {
      c[1] = static_cast<float>(1. + 4. * (vmin + 0.75 * dv - v) / dv);
      c[2] = 0.;
   }

   return c;
}


} // namespace util
} // namespace qserl
