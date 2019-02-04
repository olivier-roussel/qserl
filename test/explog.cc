//
// Copyright (c) 2016,2018 CNRS
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include <qserl/util/explog.h>

using Eigen::Vector3d;
using Eigen::Matrix3d;

typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,4,4> Matrix4d;

BOOST_AUTO_TEST_SUITE ( explog )

BOOST_AUTO_TEST_CASE(exp)
{
  Matrix4d M(Matrix4d::Random());

  Vector6d v(Vector6d::Random());

  Matrix3d R = qserl::exp3(v.tail<3>());
  BOOST_CHECK(R.isApprox(Eigen::AngleAxis<double>(v.tail<3>().norm(), v.tail<3>().normalized()).matrix()));
  BOOST_CHECK(v.tail<3>().isApprox (qserl::log3(R)));

  Matrix3d R0 = qserl::exp3(Vector3d::Zero());
  BOOST_CHECK(R0.isIdentity());
  
  M = qserl::exp6(v);
  
  BOOST_CHECK(R.isApprox(M.topLeftCorner<3,3>()));
  BOOST_CHECK(v.isApprox (qserl::log6(M)));
}

/*
BOOST_AUTO_TEST_CASE(log)
{
  SE3 M(SE3::Identity());
  Motion v(Motion::Random()); v.linear().setZero();
  
  SE3::Vector3 omega = log3(M.rotation());
  BOOST_CHECK(omega.isZero());
  
  M.setRandom();
  M.translation().setZero();
  
  v = log6(M);
  omega = log3(M.rotation());
  BOOST_CHECK(omega.isApprox(v.angular()));
  
  // Quaternion
  Eigen::Quaterniond quat(SE3::Matrix3::Identity());
  omega = quaternion::log3(quat);
  BOOST_CHECK(omega.isZero());

  for(int k = 0; k < 1e3; ++k)
  {
    SE3::Vector3 w; w.setRandom();
    quaternion::exp3(w,quat);
    SE3::Matrix3 rot = exp3(w);
    
    BOOST_CHECK(quat.toRotationMatrix().isApprox(rot));
    double theta;
    omega = quaternion::log3(quat,theta);
    const double PI_value = PI<double>();
    BOOST_CHECK(omega.norm() <= PI_value);
    double theta_ref;
    SE3::Vector3 omega_ref = log3(quat.toRotationMatrix(),theta_ref);
    
    BOOST_CHECK(omega.isApprox(omega_ref));
  }


  // Check QuaternionMap
  Eigen::Vector4d vec4;
  Eigen::QuaternionMapd quat_map(vec4.data());
  quat_map = SE3::Random().rotation();
  BOOST_CHECK(quaternion::log3(quat_map).isApprox(log3(quat_map.toRotationMatrix())));
}

BOOST_AUTO_TEST_CASE(explog3)
{
  SE3 M(SE3::Random());
  SE3::Matrix3 M_res = exp3(log3(M.rotation()));
  BOOST_CHECK(M_res.isApprox(M.rotation()));
  
  Motion::Vector3 v; v.setRandom();
  Motion::Vector3 v_res = log3(exp3(v));
  BOOST_CHECK(v_res.isApprox(v));
}

BOOST_AUTO_TEST_CASE(explog3_quaternion)
{
  SE3 M(SE3::Random());
  Eigen::Quaterniond quat;
  quat = M.rotation();
  Eigen::Quaterniond quat_res;
  quaternion::exp3(quaternion::log3(quat),quat_res);
  BOOST_CHECK(quat_res.isApprox(quat) || quat_res.coeffs().isApprox(-quat.coeffs()));
  
  Motion::Vector3 v; v.setRandom();
  quaternion::exp3(v,quat);
  BOOST_CHECK(quaternion::log3(quat).isApprox(v));
  
  SE3::Matrix3 R_next = M.rotation() * exp3(v);
  Motion::Vector3 v_est = log3(M.rotation().transpose() * R_next);
  BOOST_CHECK(v_est.isApprox(v));
  
  SE3::Quaternion quat_v;
  quaternion::exp3(v,quat_v);
  SE3::Quaternion quat_next = quat * quat_v;
  v_est = quaternion::log3(quat.conjugate() * quat_next);
  BOOST_CHECK(v_est.isApprox(v));
}

BOOST_AUTO_TEST_CASE(explog6)
{
  SE3 M(SE3::Random());
  SE3 M_res = exp6(log6(M));
  BOOST_CHECK(M_res.isApprox(M));
  
  Motion v(Motion::Random());
  Motion v_res = log6(exp6(v));
  BOOST_CHECK(v_res.toVector().isApprox(v.toVector()));
}

BOOST_AUTO_TEST_CASE (test_basic)
{
  typedef pinocchio::SE3::Vector3 Vector3;
  typedef pinocchio::SE3::Matrix3 Matrix3;
  typedef Eigen::Matrix4d Matrix4;
  typedef pinocchio::Motion::Vector6 Vector6;
  
  const double EPSILON = 1e-12;
  
  // exp3 and log3.
  Vector3 v3(Vector3::Random());
  Matrix3 R(pinocchio::exp3(v3));
  BOOST_CHECK(R.transpose().isApprox(R.inverse(), EPSILON));
  BOOST_CHECK_SMALL(R.determinant() - 1.0, EPSILON);
  Vector3 v3FromLog(pinocchio::log3(R));
  BOOST_CHECK(v3.isApprox(v3FromLog, EPSILON));
  
  // exp6 and log6.
  pinocchio::Motion nu = pinocchio::Motion::Random();
  pinocchio::SE3 m = pinocchio::exp6(nu);
  BOOST_CHECK(m.rotation().transpose().isApprox(m.rotation().inverse(),
                                                EPSILON));
  BOOST_CHECK_SMALL(m.rotation().determinant() - 1.0, EPSILON);
  pinocchio::Motion nuFromLog(pinocchio::log6(m));
  BOOST_CHECK(nu.linear().isApprox(nuFromLog.linear(), EPSILON));
  BOOST_CHECK(nu.angular().isApprox(nuFromLog.angular(), EPSILON));
  
  Vector6 v6(Vector6::Random());
  pinocchio::SE3 m2(pinocchio::exp6(v6));
  BOOST_CHECK(m2.rotation().transpose().isApprox(m2.rotation().inverse(),
                                                 EPSILON));
  BOOST_CHECK_SMALL(m2.rotation().determinant() - 1.0, EPSILON);
  Matrix4 M = m2.toHomogeneousMatrix();
  pinocchio::Motion nu2FromLog(pinocchio::log6(M));
  Vector6 v6FromLog(nu2FromLog.toVector());
  BOOST_CHECK(v6.isApprox(v6FromLog, EPSILON));
}
*/

BOOST_AUTO_TEST_SUITE_END()
