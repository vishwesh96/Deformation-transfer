/*#include <iostream>

#include "Eigen/Sparse"
#include "Eigen/UmfPackSupport"

using namespace Eigen;

// On Fedora make sure to install the packages suitesparse*

// Everything MUST be double (or you will get "no matching function call" error)

// To do LU on a non-square matrix, compute A.adjoint() * A, and then do the LU solve on that
// http://en.wikipedia.org/wiki/Linear_least_squares_%28mathematics%29#Inverting_the_matrix_of_the_normal_equations

int main(int argc, char *argv[])
{
  Eigen::SparseMatrix<double> A(4,2);
  A.insert(0,0) = 1;
  A.insert(0,1) = 1;
  A.insert(1,0) = 1;
  A.insert(1,1) = 2;
  A.insert(2,0) = 1;
  A.insert(2,1) = 3;
  A.insert(3,0) = 1;
  A.insert(3,1) = 4;

  std::cout << "A: " << A << std::endl;

  Eigen::VectorXd b(4);
  b(0) = 6;
  b(1) = 5;
  b(2) = 7;
  b(3) = 10;

  std::cout << "b: " << b << std::endl;

  Eigen::VectorXd x(2);

  Eigen::SparseMatrix<double> A2 = A.adjoint() * A;
  Eigen::VectorXd b2 = A.adjoint() * b;

  // solve Ax = b using UmfPack:
  Eigen::SparseLU<Eigen::SparseMatrix<double> > lu_of_A(A2);
  

  std::cout << "x: " << x << std::endl;
  // Solution should be (3.5, 1.4)

  
  return 0;
}*/