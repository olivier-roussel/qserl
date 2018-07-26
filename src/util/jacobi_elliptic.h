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

/** Wrapper arround inner boost jacobi elliptic internal function, as
* boost impl does not allow to access to amplitude function while computing the 
* three poles sn, cn and dn.
*/
#ifndef QSERL_UTIL_JACOBI_ELLIPTIC_H_
#define QSERL_UTIL_JACOBI_ELLIPTIC_H_

#include <boost/math/special_functions/jacobi_elliptic.hpp>


namespace boost {
namespace math {

namespace detail {
  
template <class T, class Policy>
T jacobi_imp(const T& x, const T& k, T* cn, T* dn, T* am, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   if(k < 0)
   {
      *cn = policies::raise_domain_error<T>(function, "Modulus k must be positive but got %1%.", k, pol);
      *dn = *cn;
      *am = *cn;
      return *cn;
   }
   // TODO generalize this for amplitude am
   if(k > 1)
   {
      *cn = policies::raise_domain_error<T>(function, "Modulus k must be less than one but got %1%.", k, pol);
      *dn = *cn;
      *am = *cn;
      return *cn;
   }
   //if(k > 1)
   //{
   //   T xp = x * k;
   //   T kp = 1 / k;
   //   T snp, cnp, dnp;
   //   snp = jacobi_imp(xp, kp, &cnp, &dnp, pol, function);
   //   *cn = dnp;
   //   *dn = cnp;
   //   return snp * kp;
   //}
   //
   // Special cases first:
   //
   if(x == 0)
   {
      *cn = *dn = 1;
      *am = 0;
      return 0;
   }
   if(k == 0)
   {
      *cn = cos(x);
      *dn = 1;
      *am = x;
      return sin(x);
   }
   if(k == 1)
   {
      *cn = *dn = 1 / cosh(x);
      T snu = tanh(x);
      *am = asin(snu);
      return snu;
   }
   //
   // Asymptotic forms from A&S 16.13:
   //
   if(k < tools::forth_root_epsilon<T>())
   {
      T su = sin(x);
      T cu = cos(x);
      T m = k * k;
      *dn = 1 - m * su * su / 2;
      *cn = cu + m * (x - su * cu) * su / 4;
      *am = x - m * (x - su * cu) / 4;
      return su - m * (x - su * cu) * cu / 4;
   }
   /*  Can't get this to work to adequate precision - disabled for now...
   //
   // Asymptotic forms from A&S 16.15:
   //
   if(k > 1 - tools::root_epsilon<T>())
   {
      T tu = tanh(x);
      T su = sinh(x);
      T cu = cosh(x);
      T sec = 1 / cu;
      T kp = 1 - k;
      T m1 = 2 * kp - kp * kp;
      *dn = sec + m1 * (su * cu + x) * tu * sec / 4;
      *cn = sec - m1 * (su * cu - x) * tu * sec / 4;
      T sn = tu;
      T sn2 = m1 * (x * sec * sec - tu) / 4;
      T sn3 = (72 * x * cu + 4 * (8 * x * x - 5) * su - 19 * sinh(3 * x) + sinh(5 * x)) * sec * sec * sec * m1 * m1 / 512;
      return sn + sn2 - sn3;
   }*/
   T T1;
   T kc = 1 - k;
   T k_prime = k < 0.5 ? T(sqrt(1 - k * k)) : T(sqrt(2 * kc - kc * kc));
   T T0 = detail::jacobi_recurse(x, k, T(1), k_prime, 0, &T1, pol);
   *cn = cos(T0);
   *dn = cos(T0) / cos(T1 - T0);
   *am = T0;
   return sin(T0);
}

} // namespace detail

template <class T, class U, class V, class Policy>
inline typename tools::promote_args<T, U, V>::type jacobi_elliptic(T k, U theta, V* pcn, V* pdn, V* pam, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy, 
      policies::promote_float<false>, 
      policies::promote_double<false>, 
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_elliptic<%1%>(%1%)";

   value_type sn, cn, dn, am;
   sn = detail::jacobi_imp<value_type>(static_cast<value_type>(theta), static_cast<value_type>(k), &cn, &dn, &am, forwarding_policy(), function);
   if(pcn)
      *pcn = policies::checked_narrowing_cast<result_type, Policy>(cn, function);
   if(pdn)
      *pdn = policies::checked_narrowing_cast<result_type, Policy>(dn, function);
   if(pam)
      *pam = policies::checked_narrowing_cast<result_type, Policy>(am, function);
   return policies::checked_narrowing_cast<result_type, Policy>(sn, function);;
}

template <class T, class U, class V>
inline typename tools::promote_args<T, U, V>::type jacobi_elliptic(T k, U theta, V* pcn, V* pdn, V* pam)
{
   return jacobi_elliptic(k, theta, pcn, pdn, pam, policies::policy<>());
}

} // namespace math
} // namespace boost

#endif // QSERL_UTIL_JACOBI_ELLIPTIC_H_
