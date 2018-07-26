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

#ifndef QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_
#define QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_

#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Eigen Matrix serialization
namespace cereal {

// Eigen matrix save helper
template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct save_helper
{
  static void save(Archive& ar,
                   const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
  {
    const auto num_rows = m.rows();
    const auto num_cols = m.cols();
    std::vector<_Scalar> data_v(m.data(), m.data() + m.size());
    ar( make_nvp("num_rows", num_rows),
        make_nvp("num_cols", num_cols),
        make_nvp("data", data_v) );
  }
};

// Partial specialization for 2-d vectors
template <class Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
struct save_helper<Archive, _Scalar, 2, 1, _Options, _MaxRows, _MaxCols>
{
  static void save(Archive& ar,
                   const Eigen::Matrix<_Scalar, 2, 1, _Options, _MaxRows, _MaxCols>& m)
  {
    ar( make_nvp("x", m.x()),
        make_nvp("y", m.y()) );
  }
};

// Partial specialization for 3-d vectors
template <class Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
struct save_helper<Archive, _Scalar, 3, 1, _Options, _MaxRows, _MaxCols>
{
  static void save(Archive& ar,
                   const Eigen::Matrix<_Scalar, 3, 1, _Options, _MaxRows, _MaxCols>& m)
  {
    ar( make_nvp("x", m.x()),
        make_nvp("y", m.y()),
        make_nvp("z", m.z()) );
  }
};

// eigen matrix load helper
template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct load_helper
{
  static void load(Archive& ar,
                   Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
  {
    using Index = typename Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::Index;  Index num_rows, num_cols;
    std::vector<_Scalar> data_v;
    ar( make_nvp("num_rows", num_rows),
        make_nvp("num_cols", num_cols),
        make_nvp("data", data_v) );
    m = Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>(data_v.data());
  }
};

// Partial specialization for 2-d vectors
template <class Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
struct load_helper<Archive, _Scalar, 2, 1, _Options, _MaxRows, _MaxCols>
{
  static void load(Archive& ar,
                   Eigen::Matrix<_Scalar, 2, 1, _Options, _MaxRows, _MaxCols>& m)
  {
    ar( make_nvp("x", m.x()),
        make_nvp("y", m.y()) );
  }
};

// Partial specialization for 3-d vectors
template <class Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
struct load_helper<Archive, _Scalar, 3, 1, _Options, _MaxRows, _MaxCols>
{
  static void load(Archive& ar,
                   const Eigen::Matrix<_Scalar, 3, 1, _Options, _MaxRows, _MaxCols>& m)
  {
    ar( make_nvp("x", m.x()),
        make_nvp("y", m.y()),
        make_nvp("z", m.z()) );
  }
};

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void
save(Archive& ar,
     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
  save_helper<Archive, _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::save(ar, m);
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void
load(Archive& ar,
     Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m)
{
  load_helper<Archive, _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::load(ar, m);
}

// Eigen triplets serialization (for sparse matrices)
template <class Archive, typename _Scalar>
void
save(Archive& ar,
     const Eigen::Triplet<_Scalar>& m)
{
  ar( make_nvp("row", m.row()),
      make_nvp("col", m.col()),
      make_nvp("value", m.value()) );
}

template <class Archive, typename _Scalar>
void
load(Archive& ar,
     Eigen::Triplet<_Scalar>& m)
{
  using Index = typename Eigen::Triplet<_Scalar>::Index;
  Index row, col;
  _Scalar value;
  ar( make_nvp("row", row),
      make_nvp("col", col),
      make_nvp("value", value) );
  m = Eigen::Triplet<_Scalar>(row, col, value);
}

// Eigen sparse matrices serialization
template <class Archive, typename _Scalar, int _Options, typename _Index>
void
save(Archive& ar,
     const Eigen::SparseMatrix<_Scalar,_Options,_Index>& m)
{
  using Index = typename Eigen::SparseMatrix<_Scalar,_Options,_Index>::Index;
  const auto inner_size = m.innerSize();
  const auto outer_size = m.outerSize();
  using Triplet = typename Eigen::Triplet<_Scalar>;
  std::vector<Triplet> triplets;
  for(Index i = 0; i < outer_size; ++i)
  {
    for(typename Eigen::SparseMatrix<_Scalar, _Options,_Index>::InnerIterator it(m,i);
        it;
        ++it)
    {
      triplets.push_back(Triplet{it.row(), it.col(), it.value()});
    }
  }
  ar( make_nvp("inner_size", inner_size),
      make_nvp("outer_size", outer_size),
      make_nvp("triplets", triplets) );
}

template <class Archive, typename _Scalar, int _Options, typename _Index>
void
load(Archive& ar,
     Eigen::SparseMatrix<_Scalar,_Options,_Index>& m)
{
  using Index = typename Eigen::SparseMatrix<_Scalar,_Options,_Index>::Index;
  Index inner_size;
  Index outer_size;
  ar( make_nvp("inner_size", inner_size), make_nvp("outer_size", outer_size) );
  const auto num_rows = m.IsRowMajor ? outer_size : inner_size;
  const auto num_cols = m.IsRowMajor ? inner_size : outer_size;
  m.resize(num_rows, num_cols);
  using Triplet = typename Eigen::Triplet<_Scalar>;
  std::vector<Triplet> triplets;
  ar( make_nvp("triplets", triplets) );
  m.setFromTriplets(triplets.begin(), triplets.end());
}

template <class Archive, typename _Scalar, int _Rows, int _Cols>
void
save(Archive& ar,
     const Eigen::Transform<_Scalar, _Rows, _Cols>& t)
{
  ar( make_nvp("transform_matrix", t.matrix()) );
}

template <class Archive, typename _Scalar, int _Rows, int _Cols>
void
load(Archive& ar,
     Eigen::Transform<_Scalar, _Rows, _Cols>& t)
{
  typename Eigen::Transform<_Scalar, _Rows, _Cols>::MatrixType mat;
  ar( make_nvp("transform_matrix", mat) );
  t = Eigen::Transform<_Scalar, _Rows, _Cols>(mat);
}

} // namespace cereal

#endif // QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_
