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
#include <armadillo>
using namespace lemon;
int main()
{
  // Create an instance of the default MIP solver class
  // (it will represent an "empty" problem at first)
  Mip mip;
  // Add two columns (variables) to the problem
  std::vector<Mip::Col> x(2);
for(int i=0;i<2;i++){
  x[i] = mip.addCol();

}
arma::sp_mat A_ineq(3,2);
arma::vec b_ineq(3);
A_ineq(0,0)=1;
A_ineq(0,1)=-1;
b_ineq(0)=5;
A_ineq(1,0)=2;
A_ineq(1,1)=1;
b_ineq(1)=25;
A_ineq(2,0)=-2;
A_ineq(2,1)=-1;
b_ineq(2)=0;


  // Add rows (constraints) to the problem
  for(int i=0;i<3;i++){
  Mip::Expr e;
  for(int j=0;j<2;j++){
e += A_ineq(i,j) * x[j];

//mip.addRow(A_ineq(i,j)*x[j] <=b_ineq(i));

  }
mip.addRow(e <= b_ineq(i));
  }


  // Set lower and upper bounds for the columns (variables)
  mip.colLowerBound(x[0], 0);
  mip.colUpperBound(x[1], 10);

  // Set the type of the columns
  mip.colType(x[0], Mip::INTEGER);
  mip.colType(x[1], Mip::REAL);

  // Specify the objective function
  mip.max();
  mip.obj(5 * x[0] + 3 * x[1]);

  // Solve the problem using the underlying MIP solver
  mip.solve();
  mip.type();
  std::cout<<mip.type()<<std::endl;
  // Print the results
  if (mip.type() == Mip::OPTIMAL) {
    std::cout << "Objective function value: " << mip.solValue() << std::endl;
    std::cout << "x1 = " << mip.sol(x[0]) << std::endl;
    std::cout << "x2 = " << mip.sol(x[1]) << std::endl;
  } else {
    std::cout << "Optimal solution not found." << std::endl;
  }
  return 0;
}
