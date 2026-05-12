// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/numerics/polynomial_solvers.hpp"

#include <cassert>
#include <cmath>

namespace axom
{
namespace numerics
{
int solve_linear(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 2);
  assert(roots.size() >= 1);

  int status = -1;

  // solve ax + b = 0
  const double a = coeff[1];
  const double b = coeff[0];

  if(utilities::isNearlyEqual(a, 0.))
  {
    if(utilities::isNearlyEqual(b, 0.))
    {
      // Infinite solutions: a horizontal line on the X-axis.
      status = 0;
      numRoots = -1;
    }
    else
    {
      // No solutions: a horizontal line not on the X-axis.
      numRoots = 0;
    }
  }
  else
  {
    // One solution, where the line crosses the X-axis.
    status = 0;
    numRoots = 1;
    roots[0] = -b / a;
  }

  return status;
}

//------------------------------------------------------------------------------
#ifdef _MSC_VER
  #pragma warning(push)
  #pragma warning(disable : 4244)
#endif
int solve_linear(const double* coeff, double* roots, int& numRoots)
{
  return solve_linear(ArrayView<const double>(coeff,
                                              axom::StackArray<axom::IndexType, 1> {2},
                                              axom::StackArray<axom::IndexType, 1> {1}),
                      ArrayView<double>(roots,
                                        axom::StackArray<axom::IndexType, 1> {1},
                                        axom::StackArray<axom::IndexType, 1> {1}),
                      numRoots);
}

//------------------------------------------------------------------------------
int solve_quadratic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 3);
  assert(roots.size() >= 2);

  int status = -1;

  // solve ax^2 + bx + c = 0
  const double a = coeff[2];
  const double b = coeff[1];
  const double c = coeff[0];

  if(utilities::isNearlyEqual(a, 0.))
  {
    // If this system is nearly linear, solve it as such.
    return solve_linear(coeff, roots, numRoots);
  }

  const double discriminant = b * b - 4 * a * c;
  const double overtwoa = 1. / (2 * a);

  if(utilities::isNearlyEqual(discriminant, 0.))
  {
    // One unique real root
    status = 0;
    numRoots = 1;
    roots[0] = roots[1] = -b * overtwoa;
  }
  else if(discriminant < 0)
  {
    // No real roots
    numRoots = 0;
  }
  else
  {
    // Two real roots
    status = 0;
    numRoots = 2;
    const double sqrtdisc = std::sqrt(discriminant);
    roots[0] = (-b + sqrtdisc) * overtwoa;
    roots[1] = (-b - sqrtdisc) * overtwoa;
  }

  return status;
}

//------------------------------------------------------------------------------
int solve_quadratic(const double* coeff, double* roots, int& numRoots)
{
  return solve_quadratic(ArrayView<const double>(coeff,
                                                 axom::StackArray<axom::IndexType, 1> {3},
                                                 axom::StackArray<axom::IndexType, 1> {1}),
                         ArrayView<double>(roots,
                                           axom::StackArray<axom::IndexType, 1> {2},
                                           axom::StackArray<axom::IndexType, 1> {1}),
                         numRoots);
}

inline double cuberoot(double x)
{
  // pow(x, y) returns NaN for negative finite x and noninteger y.
  if(x < 0)
  {
    return -pow(-x, 1. / 3.);
  }
  else
  {
    return pow(x, 1. / 3.);
  }
}

//------------------------------------------------------------------------------
int solve_cubic(ArrayView<const double> coeff, ArrayView<double> roots, int& numRoots)
{
  assert(coeff.size() >= 4);
  assert(roots.size() >= 3);

  int status = -1;

  // Here I use variable names as presented in Korn & Korn:
  // x^3 + ax^2 + bx + c = 0
  double cubecoeff = coeff[3];
  double a = coeff[2];
  double b = coeff[1];
  double c = coeff[0];

  if(utilities::isNearlyEqual(cubecoeff, 0.))
  {
    // If this system is nearly quadratic, solve it as such.
    return solve_quadratic(coeff, roots, numRoots);
  }

  // We normalize by dividing all by the cubic coefficient.
  const double invcubecoeff = 1. / cubecoeff;
  a *= invcubecoeff;
  b *= invcubecoeff;
  c *= invcubecoeff;

  // Note that p and q differ by a multiplicative constant from Korn,
  // because they're always used with that multiplication.
  const double p = (-a * a + 3 * b) / 9;                      //  1/3 Korn's p
  const double q = (a * (-2 * a * a + 9 * b) - 27 * c) / 54;  // -1/2 Korn's q

  const double Q = p * p * p + q * q;  // actual discriminant == -108Q
  // the term chVar occurs because we've changed variables
  // (x = y - a/3) and we need to change back to x.
  const double onethird = 1. / 3.;
  const double chVar = -a * onethird;

  if(utilities::isNearlyEqual(Q, 0.))
  {
    // We have three real roots, and at least two are equal.
    if(utilities::isNearlyEqual(q, 0.))
    {
      numRoots = 1;
    }
    else
    {
      numRoots = 2;
    }
    status = 0;

    const double cuberootq = cuberoot(q);
    roots[0] = chVar + 2 * cuberootq;
    roots[1] = chVar - cuberootq;
    roots[2] = roots[1];
  }
  else if(Q > 0)
  {
    // We have one real root, and two complex roots.
    // Right now we're calculating the real root only, but the complex roots
    // can be easily added by un-commenting the calculation.
    numRoots = 1;
    status = 0;

    const double sqrtQ = sqrt(Q);
    const double A = cuberoot(q + sqrtQ);
    const double B = cuberoot(q - sqrtQ);

    roots[0] = chVar + A + B;
    roots[1] = 0;
    roots[2] = 0;
    // double imagpart = sqrt(3) * (A - B)/2;
    // double stavg = (A + B)/2;
    // Here note use of imaginary i:
    // roots[1] = chVar - stavg - i*imagpart;
    // roots[2] = chVar - stavg + i*imagpart;
  }
  else
  {
    // Q < 0, and we have three distinct real roots.
    // The case for Q == 0 is a special case of the case for Q > 0, and
    // both correspond to Cardano's method proper, reported in section
    // 1.8-3 of Korn, where the quadratic term ax^2 is eliminated by a
    // change of variable x = y - a/3.  For "irreducible" cubics (this
    // case, where the variable substitution didn't eliminate the quadratic),
    // we use the trigonometric solution to the cubic equation, reported
    // in section 1.8-4a of Korn.
    numRoots = 3;
    status = 0;

    const double alpha = acos(q / sqrt(-p * p * p));
    const double m = 2 * sqrt(-p);

    roots[0] = chVar + m * cos(alpha * onethird);
    roots[1] = chVar - m * cos((alpha + M_PI) * onethird);
    roots[2] = chVar - m * cos((alpha - M_PI) * onethird);
  }

  return status;
}

//------------------------------------------------------------------------------
int solve_cubic(const double* coeff, double* roots, int& numRoots)
{
  return solve_cubic(ArrayView<const double>(coeff,
                                             axom::StackArray<axom::IndexType, 1> {4},
                                             axom::StackArray<axom::IndexType, 1> {1}),
                     ArrayView<double>(roots,
                                       axom::StackArray<axom::IndexType, 1> {3},
                                       axom::StackArray<axom::IndexType, 1> {1}),
                     numRoots);
}

#ifdef _MSC_VER
  #pragma warning(pop)
#endif

} /* end namespace numerics */
} /* end namespace axom */
