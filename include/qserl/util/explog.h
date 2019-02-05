/**
* Copyright (c) 2019 CNRS
* Author: Joseph Mirabel
* Inspired by file explog.hpp of Pinocchio library.
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

#ifndef QSERL_EXPLOG_H_
#define QSERL_EXPLOG_H_

#define QSERL_ASSERT_MATRIX_SPECIFIC_SIZE(type,M,nrows,ncols)              \
  EIGEN_STATIC_ASSERT(   (type::RowsAtCompileTime == Eigen::Dynamic || type::RowsAtCompileTime == nrows) \
                      && (type::ColsAtCompileTime == Eigen::Dynamic || type::ColsAtCompileTime == ncols),\
                      THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);    \
  assert(M.rows()==nrows && M.cols()==ncols);

#ifdef __linux__
# define SINCOS(a,sa,ca) sincos(a,sa,ca);
#elif __APPLE__
# define SINCOS(a,sa,ca) __sincos(a,sa,ca);
#else // if sincos specialization does not exist
# define SINCOS(a,sa,ca) (*sa) = std::sin(a); (*ca) = std::cos(a);
#endif

#include <cmath>
#include <Eigen/Geometry>
#include <qserl/util/constants.h>

namespace qserl
{
  template<typename Scalar>
  struct TaylorSeriesExpansion
  {
    ///
    /// \brief Computes the expected tolerance of the argument of a Taylor series expansion for a certain degree
    ///        according to the machine precision of the given input Scalar.
    ///
    /// \tparam degree the degree of the Taylor series expansion.
    ///
    template<int degree>
    static Scalar precision()
    {
      static Scalar value = std::pow(std::numeric_limits<Scalar>::epsilon(),Scalar(1)/Scalar(degree+1));
      return value;
    }
  }; // struct TaylorSeriesExpansion

  template<typename Scalar>
  inline typename Eigen::Matrix<Scalar,4,4>
  inv(const Eigen::Matrix<Scalar,4,4> & M)
  {
    typedef Eigen::Matrix<Scalar,4,4> Matrix4;
    typedef Eigen::Block<const Matrix4,3,3> Block33;
    typedef Eigen::Block<const Matrix4,3,1> Block31;

    Matrix4 Minv;
    const Block33 R = M.template topLeftCorner <3,3>();
    const Block31 p = M.template topRightCorner<3,1>();

    Minv <<
      R.transpose(), -R.transpose() * p,
      0,0,0        , 1;

    return Minv;
  }

  /// \brief Exp: so3 -> SO3.
  ///
  /// Return the integral of the input angular velocity during time 1.
  ///
  /// \param[in] v The angular velocity vector.
  ///
  /// \return The rotational matrix associated to the integration of the angular velocity during time 1.
  ///
  template<typename Vector3Like>
  typename Eigen::Matrix<typename Vector3Like::Scalar,3,3>
  exp3(const Eigen::MatrixBase<Vector3Like> & v)
  {
    QSERL_ASSERT_MATRIX_SPECIFIC_SIZE (Vector3Like, v, 3, 1);

    typedef typename Vector3Like::Scalar Scalar;
    typedef Eigen::Matrix<Scalar,3,3> Matrix3;
    
    const Scalar t2 = v.squaredNorm();
    
    const Scalar t = std::sqrt(t2);
    if(t > TaylorSeriesExpansion<Scalar>::template precision<3>())
    {
      Scalar ct,st; SINCOS(t,&st,&ct);
      const Scalar alpha_vxvx = (1 - ct)/t2;
      const Scalar alpha_vx = (st)/t;
      Matrix3 res(alpha_vxvx * v * v.transpose());
      res.coeffRef(0,1) -= alpha_vx * v[2]; res.coeffRef(1,0) += alpha_vx * v[2];
      res.coeffRef(0,2) += alpha_vx * v[1]; res.coeffRef(2,0) -= alpha_vx * v[1];
      res.coeffRef(1,2) -= alpha_vx * v[0]; res.coeffRef(2,1) += alpha_vx * v[0];
      res.diagonal().array() += ct;
      
      return res;
    }
    else
    {
      const Scalar alpha_vxvx = Scalar(1)/Scalar(2) - t2/24;
      const Scalar alpha_vx = Scalar(1) - t2/6;
      Matrix3 res(alpha_vxvx * v * v.transpose());
      res.coeffRef(0,1) -= alpha_vx * v[2]; res.coeffRef(1,0) += alpha_vx * v[2];
      res.coeffRef(0,2) += alpha_vx * v[1]; res.coeffRef(2,0) -= alpha_vx * v[1];
      res.coeffRef(1,2) -= alpha_vx * v[0]; res.coeffRef(2,1) += alpha_vx * v[0];
      res.diagonal().array() += Scalar(1) - t2/2;
      
      return res;
    }
  }
  
  /// \brief Same as \ref log3
  ///
  /// \param[in] R the rotation matrix.
  /// \param[out] theta the angle value.
  ///
  /// \return The angular velocity vector associated to the rotation matrix.
  ///
  template<typename Matrix3Like>
  Eigen::Matrix<typename Matrix3Like::Scalar,3,1>
  log3(const Eigen::MatrixBase<Matrix3Like> & R,
       typename Matrix3Like::Scalar & theta)
  {
    QSERL_ASSERT_MATRIX_SPECIFIC_SIZE (Matrix3Like, R, 3, 3);

    typedef typename Matrix3Like::Scalar Scalar;
    typedef Eigen::Matrix<Scalar,3,1> Vector3;
    
    Vector3 res;
    const Scalar tr = R.trace();
    if (tr > Scalar(3))       theta = 0; // acos((3-1)/2)
    else if (tr < Scalar(-1)) theta = constants::pi; // acos((-1-1)/2)
    else                      theta = acos((tr - Scalar(1))/Scalar(2));
    assert(theta == theta && "theta contains some NaN"); // theta != NaN
    
    // From runs of hpp-constraints/tests/logarithm.cc: 1e-6 is too small.
    if (theta < constants::pi - 1e-2)
    {
      const Scalar t = ((theta > TaylorSeriesExpansion<Scalar>::template precision<3>())
                        ? theta / sin(theta)
                        : Scalar(1)) / Scalar(2);
      res(0) = t * (R (2, 1) - R (1, 2));
      res(1) = t * (R (0, 2) - R (2, 0));
      res(2) = t * (R (1, 0) - R (0, 1));
    }
    else
    {
      // 1e-2: A low value is not required since the computation is
      // using explicit formula. However, the precision of this method
      // is the square root of the precision with the antisymmetric
      // method (Nominal case).
      const Scalar cphi = cos(theta - constants::pi);
      const Scalar beta  = theta*theta / ( Scalar(1) + cphi );
      Vector3 tmp((R.diagonal().array() + cphi) * beta);
      res(0) = (R (2, 1) > R (1, 2) ? Scalar(1) : Scalar(-1)) * (tmp[0] > Scalar(0) ? sqrt(tmp[0]) : Scalar(0));
      res(1) = (R (0, 2) > R (2, 0) ? Scalar(1) : Scalar(-1)) * (tmp[1] > Scalar(0) ? sqrt(tmp[1]) : Scalar(0));
      res(2) = (R (1, 0) > R (0, 1) ? Scalar(1) : Scalar(-1)) * (tmp[2] > Scalar(0) ? sqrt(tmp[2]) : Scalar(0));
    }
    
    return res;
  }
  
  /// \brief Log: SO3 -> so3.
  ///
  /// Pseudo-inverse of log from \f$ SO3 -> { v \in so3, ||v|| \le pi } \f$.
  ///
  /// \param[in] R The rotation matrix.
  ///
  /// \return The angular velocity vector associated to the rotation matrix.
  ///
  template<typename Matrix3Like>
  Eigen::Matrix<typename Matrix3Like::Scalar,3,1>
  log3(const Eigen::MatrixBase<Matrix3Like> & R)
  {
    QSERL_ASSERT_MATRIX_SPECIFIC_SIZE (Matrix3Like, R, 3, 3);

    typename Matrix3Like::Scalar theta;
    return log3(R.derived(),theta);
  }

  /// \brief Exp: se3 -> SE3.
  ///
  /// Return the integral of the input spatial velocity during time 1.
  ///
  /// \param[in] v The twist represented by a vector.
  ///
  /// \return The rigid transformation associated to the integration of the twist vector during time 1..
  ///
  template<typename Vector6Like>
  Eigen::Matrix<typename Vector6Like::Scalar,4,4>
  exp6(const Eigen::MatrixBase<Vector6Like> & nu)
  {
    QSERL_ASSERT_MATRIX_SPECIFIC_SIZE (Vector6Like, nu, 6, 1);

    typedef typename Vector6Like::Scalar Scalar;
    typedef Eigen::Matrix<Scalar,4,4> Matrix4;
    typedef Eigen::VectorBlock<const Vector6Like,3> VectorBlock3;
    typedef Eigen::Block<Matrix4,3,3> Block33;
    typedef Eigen::Block<Matrix4,3,1> Block31;

    Matrix4 res; res.row(3) << 0,0,0,1;
    Block33 R = res.template topLeftCorner <3,3>();
    Block31 p = res.template topRightCorner<3,1>();
    
    const VectorBlock3 w (nu.template head<3>()); // Angular
    const VectorBlock3 v (nu.template tail<3>()); // Linear
    
    const Scalar t2 = w.squaredNorm();
    
    const Scalar t = std::sqrt(t2);
    if(t < TaylorSeriesExpansion<Scalar>::template precision<3>())
    {
      // Taylor expansion
      const Scalar alpha_wxv = Scalar(1)/Scalar(2) - t2/24;
      const Scalar alpha_v = Scalar(1) - t2/6;
      const Scalar alpha_w = (Scalar(1)/Scalar(6) - t2/120)*w.dot(v);
      
      // Linear
      p.noalias() = (alpha_v*v + alpha_w*w + alpha_wxv*w.cross(v));
      
      // Rotational
      R.noalias() = alpha_wxv * w * w.transpose();
      R.coeffRef(0,1) -= alpha_v * w[2]; R.coeffRef(1,0) += alpha_v * w[2];
      R.coeffRef(0,2) += alpha_v * w[1]; R.coeffRef(2,0) -= alpha_v * w[1];
      R.coeffRef(1,2) -= alpha_v * w[0]; R.coeffRef(2,1) += alpha_v * w[0];
      R.diagonal().array() += Scalar(1) - t2/2;
    }
    else
    {
      Scalar ct,st; SINCOS(t,&st,&ct);
      
      const Scalar inv_t2 = Scalar(1)/t2;
      const Scalar alpha_wxv = (Scalar(1) - ct)*inv_t2;
      const Scalar alpha_v = (st)/t;
      const Scalar alpha_w = (Scalar(1) - alpha_v)*inv_t2 * w.dot(v);
      
      // Linear
      p.noalias() = (alpha_v*v + alpha_w*w + alpha_wxv*w.cross(v));
      
      // Rotational
      R.noalias() = alpha_wxv * w * w.transpose();
      R.coeffRef(0,1) -= alpha_v * w[2]; R.coeffRef(1,0) += alpha_v * w[2];
      R.coeffRef(0,2) += alpha_v * w[1]; R.coeffRef(2,0) -= alpha_v * w[1];
      R.coeffRef(1,2) -= alpha_v * w[0]; R.coeffRef(2,1) += alpha_v * w[0];
      R.diagonal().array() += ct;
    }
    
    return res;
  }

  /// \brief Log: SE3 -> se3.
  ///
  /// Pseudo-inverse of exp from SE3 -> { v,w \in se3, ||w|| < 2pi }.
  ///
  /// \param[in] R The rigid transformation represented as an homogenous matrix.
  ///
  /// \return The twist associated to the rigid transformation during time 1.
  ///
  template<typename Matrix4Like>
  Eigen::Matrix<typename Matrix4Like::Scalar,6,1>
  log6(const Eigen::MatrixBase<Matrix4Like> & M)
  {
    QSERL_ASSERT_MATRIX_SPECIFIC_SIZE (Matrix4Like, M, 4, 4);

    typedef typename Matrix4Like::Scalar Scalar;
    typedef Eigen::Matrix<Scalar,6,1> Vector6;
    typedef Eigen::VectorBlock<Vector6,3> VectorBlock3;
    typedef Eigen::Block<const Matrix4Like,3,3> Block33;
    typedef Eigen::Block<const Matrix4Like,3,1> Block31;

    const Block33 R = M.template topLeftCorner <3,3>();
    const Block31 p = M.template topRightCorner<3,1>();
    
    Scalar t;
    Vector6 res;
    VectorBlock3 w (res.template head<3>()); // Angular
    w = log3(R,t); // t in [0,π]
    const Scalar t2 = t*t;
    Scalar alpha, beta;
    if (t < TaylorSeriesExpansion<Scalar>::template precision<3>())
    {
      alpha = Scalar(1) - t2/Scalar(12) - t2*t2/Scalar(720);
      beta = Scalar(1)/Scalar(12) + t2/Scalar(720);
    }
    else
    {
      Scalar st,ct; SINCOS(t,&st,&ct);
      alpha = t*st/(Scalar(2)*(Scalar(1)-ct));
      beta = Scalar(1)/t2 - st/(Scalar(2)*t*(Scalar(1)-ct));
    }
    
    // Linear
    res.template tail<3>() = alpha * p - 0.5 * w.cross(p) + beta * w.dot(p) * w;
    return res;
  }

} // namespace pinocchio

#endif // QSERL_EXPLOG_H_
