/**
* Copyright (c) 2019 CNRS
* Author: Joseph Mirabel
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

#include <qserl/rod3d/ik.h>
#include <qserl/util/explog.h>

#include <iostream>

namespace qserl {
namespace rod3d {

  InverseKinematics::InverseKinematics (const RodConstShPtr& rod) :
    m_rod (rod),
    m_squareErrorThr (1e-6),
    m_maxIter (20),
    m_verbosity (INT_MAX)
  {}

  bool InverseKinematics::compute (const WorkspaceIntegratedStateShPtr& state,
      std::size_t iNode, Displacement oMi) const
  {
    Displacement iMo (inv(oMi)), iMt;
    Wrench w (state->wrench (0)), dw;
    typedef Eigen::Matrix<double,6,1> Vector6;
    Vector6 error;
    double scale = 1.;

    typedef Eigen::FullPivLU<Matrix6d> Decomposition;
    Decomposition decomposition (6,6);

    int iter = m_maxIter;
    while (true) {
      iMt = iMo * state->nodes()[iNode];
      error = log6 (iMt);
      double errorNorm2 = error.squaredNorm();
      if (iter % m_verbosity == 0)
        std::cout << iter << '\t' << errorNorm2 << '\t' << w.transpose() << std::endl;
      if (errorNorm2 < m_squareErrorThr) return true;
      if (iter == 0) return false;

      const Matrix6d& J (state->getJMatrix (iNode));
      decomposition.compute (J);
      if (!decomposition.isInvertible()) return false;
      dw = decomposition.solve (error);

      w -= scale * dw;

      WorkspaceIntegratedState::IntegrationResultT result
        = state->integrateFromBaseWrenchRK4 (w);

      if (result != WorkspaceIntegratedState::IR_VALID) return false;

      iter--;
    }
    return false;
  }
}  // namespace rod3d
}  // namespace qserl
