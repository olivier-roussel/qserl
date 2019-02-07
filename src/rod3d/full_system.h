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

#ifndef QSERL_3D_FULL_SYSTEM_H_
#define QSERL_3D_FULL_SYSTEM_H_

#include "qserl/exports.h"

#include <functional>
#include "qserl/rod3d/workspace_integrated_state.h"

namespace qserl {
namespace rod3d {

class QSERL_EXPORT FullSystem
{
public:
  typedef std::array<double, 94> state_type; /**< 6 first are costate mu,
                                                  16 following for state q (4x4 matrix),
                                                  72 following for M and J matrices resp (2 x (6x6) matrices).
                                                  See index helpers <x>_index() methods below. */

  /**
  * Constructors, destructors
  */
  FullSystem(const Parameters& i_params,
                double i_dt);

  virtual ~FullSystem();

  void
  operator()(const state_type& i_x,
             state_type& o_dxdt,
             double i_t);

  /** Returns default state value. */
  static state_type
  defaultState();

  /**
   * @return starting index of costate mu within a full state
   */
  static size_t
  mu_index() { return 0ul; }

  /**
   * @return starting index of state q within a full state
   */
  static size_t
  q_index() { return 6ul; }

  /**
   * @return starting index of jacobian state (M and J matrices) within a full state
   */
  static size_t
  MJ_index() { return 22ul; }

  /**
   * @brief The determinant of the jacobian, starting at 0 for t=0, is checked for zero-crossing
   * for t > t_threshold, where t_threshold is the first value of t (starting at t = 0)
   * where |det(J(t))| = stability_threshold
   * @return stability_threshold
   */
  double
  jacobianStabilityThreshold() const;

  /**
   * @brief The determinant of the jacobian, starting at 0 for t=0, is checked for zero-crossing
   * for t > t_threshold, where t_threshold is the first value of t (starting at t = 0)
   * where |det(J(t))| = stability_threshold
   * @param stability_threshold
   */
  void
  jacobianStabilityThreshold(double stability_threshold);

  /**
   * @brief The determinant of the jacobian will be assumed as null if |det(J(t))| < stability_tolerance,
   * for any t \in [t_threshold, t_max]. See jacobianStabilityThreshold() for t_threshold description.
   * @return stability_tolerance
   */
  double
  jacobianStabilityTolerance() const;

  /**
 * @brief The determinant of the jacobian will be assumed as null if |det(J(t))| < stability_tolerance,
 * for any t \in [t_threshold, t_max]. See jacobianStabilityThreshold() for t_threshold description.
 * @param stability_tolerance
 */
  void
  jacobianStabilityTolerance(double stability_tolerance);

private:

  Eigen::Matrix<double, 6, 1> m_inv_c;    /**< Inverse stiffness coefficients (already stored in parameters, but used to speedup the computation. */
  Eigen::Matrix<double, 6, 1> m_b;      /**<	Precomputed values from inverse stiffness coefficients, where:
                                          b(1) = inv_c(3) - inv_c(2)
                                          b(2) = inv_c(1) - inv_c(3)
                                          b(3) = inv_c(2) - inv_c(1)
                                          b(4) = inv_c(6) - inv_c(5)
                                          b(5) = inv_c(4) - inv_c(6)
                                          b(6) = inv_c(5) - inv_c(4)
                                          */
  Eigen::Vector4d m_w_x_0;                /** Gravity field in base frame. */
  Parameters m_rodParameters;
  double m_dt;
  double m_stability_threshold;
  double m_stability_tolerance;

  std::function<void(const state_type&,
                     state_type&,
                     double)> m_evaluationCallback;

  /**
  * Derivative evaluation at time t for the inextensible (RM_INEXTENSIBLE) rod model.
  */
  void
  evaluateInextensible(const state_type& i_x,
                       state_type& o_dxdt,
                       double i_t);

  /**
  * Derivative evaluation at time t for the extensible (RM_EXTENSIBLE_SHEARABLE) rod model.
  */
  void
  evaluateExtensibleShearable(const state_type& i_x,
                              state_type& o_dxdt,
                              double i_t);

  /**
  * Derivative evaluation at time t for the inextensible with gravity (RM_INEXTENSIBLE_WITH_GRAVITY) rod model.
  */
  void
  evaluateInextensibleWithGravity(const state_type& i_x,
                                  state_type& o_dxdt,
                                  double i_t);
};

}  // namespace rod3d
}  // namespace qserl

#endif // QSERL_3D_FULL_SYSTEM_H_
