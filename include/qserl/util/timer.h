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

#ifndef QSERL_UTIL_TIMER_H_
#define QSERL_UTIL_TIMER_H_

#include <boost/chrono.hpp>

namespace qserl {
namespace util {

typedef boost::chrono::high_resolution_clock::time_point TimePoint;

TimePoint getTimePoint();

boost::chrono::nanoseconds getElapsedTimeNsec(const TimePoint& start);
boost::chrono::microseconds getElapsedTimeUsec(const TimePoint& start);
boost::chrono::milliseconds getElapsedTimeMsec(const TimePoint& start);
boost::chrono::seconds getElapsedTimeSec(const TimePoint& start);

} // namespace util
} // namespace qserl

#endif // QSERL_UTIL_TIMER_H_