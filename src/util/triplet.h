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

#ifndef QSERL_UTIL_TRIPLET_H_
#define QSERL_UTIL_TRIPLET_H_

namespace util {

/** \struct Triplet
* \brief 3-Tuple generic structure which can handle different types. \n
* The behaviour is very similar to the std::pair, and it can be seen just as an
* extension of std::pair to triplets.
*/

template<class T1, class T2, class T3>
struct Triplet
{
  T1 first;   /**< First element of the triplet. */
  T2 second;  /**< Second element of the triplet. */
  T3 third;   /**< Third element of the triplet. */

  /** \fn Triplet()
  * \brief Empty constructor.
  */
  Triplet() :
      first(T1()), second(T2()), third(T3())
  {
  }

  /** \fn Triplet(const T1& x_, const T2& y_, const T3& z_)
  * \brief Parametrized constructor.
  * \param[in] x_ Value of first element of triplet.
  * \param[in] y_ Value of second element of triplet.
  * \param[in] z_ Value of third element of triplet.
  */
  Triplet(const T1& x_,
          const T2& y_,
          const T3& z_) :
      first(x_), second(y_), third(z_)
  {
  }

  /** \fn template <class U, class V, class W> Triplet (const Triplet<U,V,W> & p_)
  * \brief Copy constructor.
  * \param[in] p_ Triplet to copy.
  */
  template<class U, class V, class W>
  Triplet(const Triplet<U, V, W>& p_) :
      first(p_.first), second(p_.second), third(p_.third)
  {
  }

};

/** \struct UniformTriplet
* \brief 3-Tuple generic structure which contains uniform types. \n
*/
template<class T>
struct UniformTriplet
{

  // -------------------- Static members ------------------------------------

  static const ::size_t npos = -1;  /**< Invalid index value. */


  // -------------------- Attributes ----------------------------------------

  T e[3];                         /**< Triplet elements. */


  // -------------------- Constructors / Destructors ------------------------

  /** \fn UniformTriplet()
  * \brief Empty constructor.
  */
  UniformTriplet()
  {
  }

  /** \fn UniformTriplet(const T& x_, const T& y_, const T& z_)
  * \brief Parametrized constructor.
  * \param[in] x_ Value of first element of uniform triplet.
  * \param[in] y_ Value of second element of uniform triplet.
  * \param[in] z_ Value of third element of uniform triplet.
  */
  UniformTriplet(const T& x_,
                 const T& y_,
                 const T& z_)
  {
    e[0] = x_;
    e[1] = y_;
    e[2] = z_;
  }

  /** \fn UniformTriplet (const UniformTriplet<T> &p_)
  * \brief Copy constructor.
  * \param[in] p_ Uniform triplet to copy.
  */
  UniformTriplet(const UniformTriplet<T>& p_)
  {
    e[0] = p_.e[0];
    e[1] = p_.e[1];
    e[2] = p_.e[2];
  }

  ~UniformTriplet()
  {
  }

  // ------------------------ Methods ---------------------------------------

  /** \fn const T& operator = (const T& p_) const
  * \brief Affectation operator.
  * \param[in] p_ Uniform triplet to copy.
  */
  const T&
  operator=(const T& p_) const
  {
    e[0] = p_.e[0];
    e[1] = p_.e[1];
    e[2] = p_.e[2];
    return *this;
  }

  /** \fn const T& operator[](int idx_) const
  * \brief Const accessor to triplet element.
  * \param[in] idx_ Index of the triplet element to acces.
  * \return the element at the specified index position in the triplet.
  * \pre index must be between 0 and 2 included.
  */
  const T&
  operator[](int idx_) const
  {
    return e[idx_];
  }

  /** \fn T& operator[](int idx_)
  * \brief Accessor to triplet element.
  * \param[in] idx_ Index of the triplet element to acces.
  * \return the element at the specified index position in the triplet.
  * \pre index must be between 0 and 2 included.
  */
  T&
  operator[](int idx_)
  {
    //return const_cast<T&>(static_cast<const UniformTriplet&>(*this)[idx]);
    return e[idx_];
  }

  /** \fn bool contains(const T& p_) const
  * \brief Tests if triplet contains the given element value.
  * \param[in] p_ Value to test.
  * \return true if one of the triplet element is equal to given one.
  */
  bool
  contains(const T& p_) const
  {
    return (e[0] == p_ || e[1] == p_ || e[2] == p_);
  }

  /** \fn bool isUnique() const
  * \brief Tests each element unicity within the triplet.
  * \return true if each element of the triplet is unique.
  */
  bool
  isUnique() const
  {
    return (e[0] != e[1] && e[0] != e[2] && e[1] != e[2]);
  }

  /** \fn const T& last(const T& e1_, const T& e2_) const
  * \brief Const accessor to the other element of the triplet that is not
  * equal to the two specified elements.
  * \param[in] e1_ First element of the triplet to be different from.
  * \param[in] e2_ Second element of the triplet to be different from.
  * \return the last element of the triplet that is not equal to given elements.
  * \pre the triplet \b must contains elements e1_ and e2_/
  */
  const T&
  last(const T& e1_,
       const T& e2_) const
  {
    if(e[0] == e1_)
    {
      if(e[1] == e2_)
      {
        return e[2];
      }
      else
      {
        return e[1];
      }
    }
    else if(e[1] == e1_)
    {
      if(e[2] == e2_)
      {
        return e[0];
      }
      else
      {
        return e[2];
      }
    }
    else // if (e[2] == e1_)
    if(e[0] == e2_)
    {
      return e[1];
    }
    else
    {
      return e[0];
    }
  }

  /** \fn T& last(const T& e1_, const T& e2_)
  * \brief Accessor to the other element of the triplet that is not
  * equal to the two specified elements.
  * \param[in] e1_ First element of the triplet to be different from.
  * \param[in] e2_ Second element of the triplet to be different from.
  * \return the last element of the triplet that is not equal to given elements.
  * \pre the triplet \b must contains elements e1_ and e2_/
  */
  T&
  last(const T& e1_,
       const T& e2_)
  {
    return const_cast<T&>(static_cast<const UniformTriplet&>(*this).last(e1_, e2_));
  }

  /** \fn size_t indexOf(const T& e_) const
  * \brief Accessor to the index of given element within the triplet, if
  * triplet contains this element.
  * \param[in] e_ Element of the triplet to find index.
  * \return the index within the triplet of the element, or npos if this triplet
  * do not contains given element.
  */
  ::size_t
  indexOf(const T& e_) const
  {
    if(e[0] == e_)
    {
      return 0;
    }
    else if(e[1] == e_)
    {
      return 1;
    }
    else if(e[2] == e_)
    {
      return 2;
    }
    else
    {
      return npos;
    }
  }
};

} // namespace util

#endif // QSERL_UTIL_TRIPLET_H_