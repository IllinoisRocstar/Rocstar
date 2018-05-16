///
/// @file
/// @ingroup gridconversion_group
/// @brief Example C++ source code file for GridConversion.
///
/// This file is part of an example C++ project
/// at IllinoisRocstar.   It demonstrates doxygen
/// usage and provides example constructs and 
/// fixtures to be used in example tests.
///
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @bug No known bugs
/// @ingroup gridconversion_group
///
// Doxygen processes all comments with an extra
// "/". Comments at the top of the file as above 
// are interpreted by Doxygen as the documentation
// for the file. While two-slash comments like these
// are ignored.
#include <cmath>
#include "ExampleHeader.H"


namespace GridConversion {

  // Note that the doxygen documentation for this function is with the function
  // prototype in the ExampleHeader.H.  If this were a substantial function, then
  // the regular inline code comments can be used to comment the actual implementation
  // while doxygen documentation (with the exception of that for the *file* itself)
  // can be entirely inside the header files, or even outside of any source file.
  std::string ExampleFunction(const std::string &instring) 
  { 
    return(instring); 
  };


  // Quadrature by the trapezoid rule
  // This function integrates f from x0 to xn using n partitions of the domain
  double TrapezoidQuadrature(double (*f)(double),double x0,double xn,int n)
  {
    double h = (xn - x0)/(static_cast<double>(n));
    if(std::fabs(h) < 1e-12) throw 1;
    double sum = .5*(f(x0)+f(xn));
    for(int i = 1; i < n;i++)
      sum += f(x0 + i*h);
    return(h*sum);
  };

  // Quadrature by the midpoint rule
  // This function integrates f from x0 to xn using n partitions of the domain
  double MidPointQuadrature(double (*f)(double),double x0,double xn,int n)
  {
    double h = (xn - x0)/(static_cast<double>(n));
    if(std::fabs(h) < 1e-12) throw 1;
    double sum = 0.0;
    for(int i = 1;i <= n;i++)
      sum += f(x0+((static_cast<double>(i)-.5)*h));
    return(h*sum);
  }
}
