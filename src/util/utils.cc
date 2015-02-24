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

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

namespace qserl {
namespace util {

std::string formatTime(boost::posix_time::ptime now)
{
  using namespace boost::posix_time;
  static std::locale loc(std::wcout.getloc(),
                         new boost::posix_time::time_facet("%Y%m%d_%H%M%S"));

  //std::basic_stringstream<wchar_t> wss;
  std::basic_stringstream<char> wss;
  wss.imbue(loc);
  wss << now;
  return wss.str();
}

/** Do not use this one directly. */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool createDirAndParents(const std::string& i_path)
{
	boost::filesystem::path dir(i_path);

	if (!boost::filesystem::exists(dir))
	{
		if (!boost::filesystem::create_directories(dir))
			return false;
	}
	return true;
}

double rand_interval(double rmin, double rmax)
{
  double t = rand() / ((double)RAND_MAX + 1);
  return (t * (rmax - rmin) + rmin);
}

Eigen::Vector4f getJetColor(double v, double vmin, double vmax)
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

double peaks(double i_x, double i_y)
{
	const double x = 0.3 * i_x;
	const double y = 0.3 * i_y;
	return 3*sqr(1-x) * exp(-sqr(x) - sqr(y+1)) -10*(0.2*x - pow(x,3) - pow(y,5)) * 
		exp(-sqr(x)-sqr(y)) - (1/3)*exp(-sqr(x+1) - sqr(y));
}

} // namespace util
} // namespace qserl
