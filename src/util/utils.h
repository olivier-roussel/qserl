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

/** Common helper functions. */

#ifndef QSERL_UTIL_UTILS_H_
#define QSERL_UTIL_UTILS_H_

#include <Eigen/Core>
#include <string>
#include <vector>

namespace qserl {
namespace util {

/**
* Square function.
*/
template<class T>
inline T
sqr(T a)
{
  return a * a;
}

/**
* Clamps the value v to the bounds [mn, mx].
*/
template<class T>
inline T
clamp(T v,
      T mn,
      T mx)
{
  return v < mn ? mn : (v > mx ? mx : v);
}

/**
* Static power function.
* \return x^y
*/
template<unsigned int X, unsigned int Y>
struct power
{
  enum
  {
    value = X * power<X, Y - 1>::value,
  };
};

/**
* Specialization of the static power function as stopping criterion.
*/
template<unsigned int X>
struct power<X, 0>
{
  enum
  {
    value = 1,
  };
};

/**
* Rounding function.
*/
template<class Real, int N>
inline Real
round(Real val)
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
inline signed char
sign(T a)
{
  return (static_cast<signed char>(a > static_cast<T>(0)) -
          static_cast<signed char>(a < static_cast<T>(0)));
}

/**
* \brief Non-null version of the signum function. Returns:
* 1 if a >= 0
* -1 if a < 0
*/
template<typename T>
inline signed char
sigpos(T a)
{
  return (a >= static_cast<T>(0) ? 1 : -1);
}

/** rand utils */

/** TODO doc. */
double
rand_interval(double rmin,
              double rmax);

/** Stats utils */

/**
* \brief Computes mean and standard deviation of a set of values
*/
template<class T>
void
computeMeanAndStdDev(const std::vector<T>& i_values,
                     T& o_mean,
                     T& o_stdDev)
{
  //o_mean = static_cast<T>(0);
  for(unsigned int i = 0; i < i_values.size(); ++i)
  {
    o_mean += i_values[i];
  }
  o_mean /= i_values.size();
  //o_stdDev = static_cast<T>(0);
  for(unsigned int i = 0; i < i_values.size(); ++i)
  {
    o_stdDev += sqr(i_values[i] - o_mean);
  }
  o_stdDev /= i_values.size();
  o_stdDev = sqrt(static_cast<T>(o_stdDev));
}

/** Split strings wrt delimiter.
*/

/** Do not use this one directly. */
std::vector<std::string>&
split(const std::string& s,
      char delim,
      std::vector<std::string>& elems);

/** But use this one instead. */
std::vector<std::string>
split(const std::string& s,
      char delim);

///**
//* \brief Performs coefficient wise equality comparison for a given precision between
//* two Eigen matrices.
//*/
//template<typename DerivedA, typename DerivedB>
//bool
//allclose(const Eigen::DenseBase<DerivedA>& a,
//         const Eigen::DenseBase<DerivedB>& b,
//         const typename DerivedA::RealScalar& rtol
//         = Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
//         const typename DerivedA::RealScalar& atol
//         = Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
//{
//  return ((a.derived() - b.derived()).array().abs()
//          <= (atol + rtol * b.derived().array().abs())).all();
//}

///** \brief Search insertion index of given element in sorted randomly accessible container.
//* Returns a pair where first is a boolean set to true if element was found in the array,
//* and second is the index where should be inserted e_ in sorted array vec_ (i.e. e_ < vec_[index] )
//*/
//template<typename V, typename T, typename FuncPred>
//inline std::pair<bool, int>
//searchInsertionIndexInSortedContainer(const V& i_cont,
//                                      const T& i_e,
//                                      const FuncPred& i_isLess)
//{
//  int binf = 0, bsup = static_cast<int>(i_cont.size()) - 1;
//  while(binf <= bsup)
//  {
//    const int mid = (binf + bsup) / 2;
//    if(i_isLess(i_e, i_cont[mid]))
//    {
//      bsup = mid - 1;
//    }
//    else if(i_isLess(i_cont[mid], i_e))
//    {
//      binf = mid + 1;
//    }
//    else
//    {
//      return std::pair<bool, int>(true, mid);
//    }
//  }
//  return std::pair<bool, int>(false, binf);
//}

/**
 * \brief Returns color based on the Jet color map for the given 1-d value v
*  which goes from vmin to vmax.
*/
Eigen::Vector4f
getJetColor(double v,
            double vmin,
            double vmax);

} // namespace util
} // namespace qserl

#endif // QSERL_UTIL_UTILS_H_
