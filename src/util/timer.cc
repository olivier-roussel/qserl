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

#include "qserl/util/timer.h"

namespace qserl {
namespace util {

TimePoint
getTimePoint()
{
  return boost::chrono::high_resolution_clock::now();
}

boost::chrono::nanoseconds
getElapsedTimeNsec(const TimePoint& start)
{
  return boost::chrono::duration_cast<boost::chrono::nanoseconds>(
      boost::chrono::high_resolution_clock::now() - start);//.count();
}

boost::chrono::microseconds
getElapsedTimeUsec(const TimePoint& start)
{
  return boost::chrono::duration_cast<boost::chrono::microseconds>(
      boost::chrono::high_resolution_clock::now() - start);//.count();
}

boost::chrono::milliseconds
getElapsedTimeMsec(const TimePoint& start)
{
  return boost::chrono::duration_cast<boost::chrono::milliseconds>(
      boost::chrono::high_resolution_clock::now() - start);//.count();
}

boost::chrono::seconds
getElapsedTimeSec(const TimePoint& start)
{
  return boost::chrono::duration_cast<boost::chrono::seconds>(
      boost::chrono::high_resolution_clock::now() - start);//.count();
}

} // namespace util
} // namespace qserl
