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

#ifndef QSERL_3D_INVERSE_KINEMATICS_H_
#define QSERL_3D_INVERSE_KINEMATICS_H_

#include <qserl/rod3d/rod.h>

namespace qserl {
namespace rod3d {

  class InverseKinematics {
    public:
      InverseKinematics (const RodConstShPtr& rod);

      bool compute (const WorkspaceIntegratedStateShPtr& state,
          std::size_t iNode, Displacement target) const;

      void setErrorThreshold (double thr)
      {
        m_squareErrorThr = thr*thr;
      }

      double getErrorThreshold () const
      {
        return std::sqrt(m_squareErrorThr);
      }

      void setMaxIter (int iter)
      {
        m_maxIter = iter;
      }

      int getMaxIter () const
      {
        return m_maxIter;
      }

      void setVerbosity (int level)
      {
        m_verbosity = level;
      }

      int getVerbosity () const
      {
        return m_verbosity;
      }

      void setScale (double scale)
      {
        m_scale = scale;
      }

      double getScale () const
      {
        return m_scale;
      }

    private:
      RodConstShPtr m_rod;
      double m_squareErrorThr;
      int m_maxIter;
      int m_verbosity;
      double m_scale;
  };

}  // namespace rod3d
}  // namespace qserl

#endif // QSERL_3D_INVERSE_KINEMATICS_H_
