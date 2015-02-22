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

/** Common helper functions. */

#ifndef QSERL_UTIL_UTILS_H_
#define QSERL_UTIL_UTILS_H_

#include <boost/array.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#pragma warning( push, 0 )
#include <Eigen/Lgsm>
#pragma warning( pop )

namespace util {

/**
* Square function.
*/
template<class T>
inline T sqr(T a) { return a*a; }

/**
* Clamps the value v to the bounds [mn, mx].
*/
template<class T>
inline T clamp(T v, T mn, T mx) { return v < mn ? mn : (v > mx ? mx : v); }

/**
* Static power function.
* \return x^y
*/
template<unsigned int X, unsigned int Y>
struct power { enum { value = X * power<X, Y - 1>::value, }; };
/**
* Specialization of the static power function as stopping criterion.
*/
template<unsigned int X>
struct power<X, 0> { enum { value = 1, }; };

/**
* Rounding function.
*/
template<class Real, int N>
inline Real round(Real val)
{
	const int p = power<10, N>::value;
	return static_cast<Real>(static_cast<int>(val * p + 0.5) / static_cast<Real>(p));
}

/**
* \brief Implementation of the signum function. Returns:
* 1 if a > 0
* 0 if a = 0
* -1 if a < 0
*/
template<typename T>
inline signed char sign(T a) { return (static_cast<signed char>(a > static_cast<T>(0)) - 
  static_cast<signed char>(a < static_cast<T>(0))); }

/**
* \brief Non-null version of the signum function. Returns:
* 1 if a >= 0
* -1 if a < 0
*/
template<typename T>
inline signed char sigpos(T a) { return (a >= static_cast<T>(0) ? 1 : -1); }

/** rand utils */

/** TODO doc. */
double rand_interval(double rmin, double rmax);

/** Stats utils */

/**
* \brief Computes mean and standard deviation of a set of values
*/
template<class T>
void computeMeanAndStdDev(const std::vector<T>& i_values, T& o_mean, T& o_stdDev)
{
	//o_mean = static_cast<T>(0);
	for (unsigned int i = 0 ; i < i_values.size() ; ++i)
		o_mean += i_values[i];
	o_mean /= i_values.size();
	//o_stdDev = static_cast<T>(0);
	for (unsigned int i = 0 ; i < i_values.size() ; ++i)
		o_stdDev += sqr(i_values[i] - o_mean);
	o_stdDev /= i_values.size();
	o_stdDev = sqrt(static_cast<T>(o_stdDev));
}

/** Eigen utils */
/** TODO doc. */

inline Eigen::Matrix3d hat(const Eigen::Vector3d& i_v)
{
	Eigen::Matrix3d out;
	out << 0, -i_v[2], i_v[1],
		i_v[2], 0, -i_v[0],
		-i_v[1], i_v[0], 0;
	return out;
}

Eigen::Matrix4d GetHomogenousMatrix(const Eigen::Displacementd & disp);

Eigen::Matrix4d GetTransformationMatrix(const Eigen::Matrix4d & H, const Eigen::Vector3d & scale);

Eigen::Matrix4d GetTransformationMatrix(const Eigen::Displacementd & disp, const Eigen::Vector3d & scale);

Eigen::Vector3d GetScaleFromMatrix(const Eigen::Matrix4d &TM);

Eigen::Vector3d GetTranslationFromMatrix(const Eigen::Matrix4d &TM);

Eigen::Rotation3d GetRotationFromMatrix(const Eigen::Matrix4d & TM);

Eigen::Matrix3d rotationMatrixFromBryantAngles(double rx, double ry, double rz);

Eigen::Matrix4d transformationMatrixFromByrantAngles(double x, double y, double z, double rx, double ry, double rz);

void transformationMatrixToByrantAngles(const Eigen::Matrix4d& trans, double& x, double& y, double& z, double& rx, double& ry, double& rz);

void rotationMatrixToByrantAngles(const Eigen::Matrix3d& rot, double& rx, double& ry, double& rz);

/** Date/time tools. *

/**
 * Returns a string from given ptime with format "CCYYMMDD_hhmmss"
 */
std::string formatTime(boost::posix_time::ptime now);


/** Split strings wrt delimiter.
*/

/** Do not use this one directly. */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
/** But use this one instead. */
std::vector<std::string> split(const std::string &s, char delim);

/**
* cfg Files contains planning problem configurations.
* Comment lines must start with a #.
* Line containing start configuration dof values must start with start keyword
* Second line is goal configuration dof values must start with goal keyword
* Dof values must be ordered as it:
* tx ty tz fx fy fz x y z rx ry rz
* (separator are spaces)
* where where wrench is given by the first six coefficient and
* displacement by the six last ones.
*/
bool readStartAngGoalConfigsFromCfgFile(const std::string& i_cfgPath, Eigen::Wrenchd& o_startWrench, Eigen::Displacementd& o_startDisp, Eigen::Wrenchd& o_goalWrench, Eigen::Displacementd& o_goalDisp);

/**
* Equivalent to the Unix mkdip -p <i_path>
*/
bool createDirAndParents(const std::string& i_path);

/**
* \brief Performs coefficient wise equality comparison for a given precision between
* two Eigen matrices.
*/
template<typename DerivedA, typename DerivedB>
bool allclose(const Eigen::DenseBase<DerivedA>& a,
              const Eigen::DenseBase<DerivedB>& b,
              const typename DerivedA::RealScalar& rtol
                  = Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
              const typename DerivedA::RealScalar& atol
                  = Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
{
  return ((a.derived() - b.derived()).array().abs()
          <= (atol + rtol * b.derived().array().abs())).all();
}

/** \brief Search insertion index of given element in sorted randomly accessible container.
* Returns a pair where first is a boolean set to true if element was found in the array,
* and second is the index where should be inserted e_ in sorted array vec_ (i.e. e_ < vec_[index] )
*/
template<typename V, typename T, typename FuncPred>
inline std::pair<bool, int> searchInsertionIndexInSortedContainer(const V& i_cont, const T& i_e, const FuncPred& i_isLess)
{
  int binf = 0, bsup = static_cast<int>(i_cont.size()) - 1;
  while (binf <= bsup)
  {
    const int mid = (binf + bsup) / 2;
    if (i_isLess(i_e, i_cont[mid]))
      bsup = mid - 1;
    else if (i_isLess(i_cont[mid], i_e))
      binf = mid + 1;
    else
      return std::pair<bool, int>(true, mid);
  }
  return std::pair<bool, int>(false, binf);
}

/** Returns color based on the Jet color map for the given 1-d value v
* which goes from vmin to vmax.
*/
Eigen::Vector4f getJetColor(double v, double vmin, double vmax);

/** Similar to the function used in matlab peaks. */
double peaks(double x, double y);

} // namespace util

#endif // QSERL_UTIL_UTILS_H_
