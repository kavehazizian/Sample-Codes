/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2010
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */
#include <iostream>
#include <lemon/lp.h>
#include <glpk.h>
#include <lemon/cplex.h>
using namespace lemon;
int main()
{
  // Create an instance of the default MIP solver class
  // (it will represent an "empty" problem at first)
  CplexMip mip;
  // Add two columns (variables) to the problem
  CplexMip::Col x1 = mip.addCol();
  CplexMip::Col x2 = mip.addCol();
  // Add rows (constraints) to the problem
  mip.addRow(x1 - 5 <= x2);
  mip.addRow(0 <= 2 * x1 + x2 <= 25);

  // Set lower and upper bounds for the columns (variables)
  mip.colLowerBound(x1, 0);
  mip.colUpperBound(x2, 10);

  // Set the type of the columns
  mip.colType(x1, CplexMip::INTEGER);
  mip.colType(x2, CplexMip::REAL);

  // Specify the objective function
  mip.max();
  mip.obj(5 * x1 + 3 * x2);

  // Solve the problem using the underlying MIP solver
  mip.solve();
  // Print the results
  if (mip.type() == CplexMip::OPTIMAL) {
    std::cout << "Objective function value: " << mip.solValue() << std::endl;
    std::cout << "x1 = " << mip.sol(x1) << std::endl;
    std::cout << "x2 = " << mip.sol(x2) << std::endl;
  } else {
    std::cout << "Optimal solution not found." << std::endl;
  }
  return 0;
}
