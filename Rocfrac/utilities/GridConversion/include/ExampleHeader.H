///
/// @file
/// @ingroup gridconversion_group
/// @brief Example C++ header file for GridConversion.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @bug No known bugs
///
/// This file is part of an example C++ project
/// at IllinoisRocstar.   It demonstrates doxygen
/// usage and provides example constructs and 
/// fixtures to be used in example tests.
///
/// At IR, we like all files to have the 
/// tags (file, brief, and ingroup). Other tags
/// are optional and as needed or desired. Such
/// as \@author and \@bug.
///
///
#ifndef __EXAMPLE_HEADER_H__
#define __EXAMPLE_HEADER_H__
#include <string>

namespace GridConversion {

  /// 
  /// Example function for GridConversion (this is a brief description).
  /// 
  /// @param instring Input string to copy
  /// @returns A copy of the input string
  /// @bug No known bugs.
  /// @note All of this documentation is generated from the comments
  /// in ExampleHeader.H
  /// @note The first sentence in a doxygen comment block before a 
  /// code construct is interpreted as the "brief" description.
  /// @warning This is an example warning.
  /// @todo Nothing to do on this example function.
  ///
  /// The detailed description for ExampleFunction actually starts
  /// above with the constructs describing the input/output, and
  /// any notes, warnings, todo, bugs, author, or other information
  /// the author would like to put in.
  ///
  /// At IR, we like for all code constructs 
  /// to have a "brief" description -first sentence in this
  /// comment block, and then a "detailed" description 
  /// as needed, but it should at least list the input params and output 
  /// description. The actual example detailed documentation for this 
  /// function now follows:
  ///
  /// This function returns a copy of the input string. If there were
  /// more to say about this, then it can be said here.  
  ///
  std::string ExampleFunction(const std::string &instring);

  ///
  /// Integrates f with composite trapezoid rule.
  ///
  /// @param f Function to integrate takes a double argument, x, and returns a double value
  /// @param x0 The lower integration domain limit.
  /// @param xn The upper integration domain limit.
  /// @param n  The number of intervals into which to break the domain.
  /// @returns The numerical approximation for \f$\int_{x_0}^{x_n}{f(x)}{dx}\f$.
  ///
  /// With \f$h = \frac{(x_n - x_0)}{n}\f$, and 
  /// \f$x_i = x_0 + {i}{h}\f$, the integral is calculated as:
  /// 
  /// \f$\int_{x_0}^{x_n}{f(x)}{dx} \approx \frac{h}{2}\sum_{i=1}^{n}(f(x_{i-1}) + f(x_i))\f$
  ///
  /// The error of this method is \f$O(h^2)\f$.
  ///
  /// @test GridConversion::TestingObject::Test__TrapezoidQuadrature tests this 
  /// function by sending linear and quadratic test functions.
  ///
  double TrapezoidQuadrature(double (*f)(double),double x0,double xn,int n);

  ///
  /// Integrates f with composite midpoint rule.
  ///
  /// @param f Function to integrate takes a double argument, x, and returns a double value
  /// @param x0 The lower integration domain limit.
  /// @param xn The upper integration domain limit.
  /// @param n  The number of intervals into which to break the domain.
  /// @returns The numerical approximation for \f$\int_{x_0}^{x_n}{f(x)}{dx}\f$.
  ///
  /// With \f$h = \frac{(x_n - x_0)}{n}\f$, and 
  /// \f$x_i = x_0 + h(i-\frac{1}{2})\f$, the integral is calculated as:
  ///
  /// \f$\int_{x_0}^{x_n}{f(x)}{dx} \approx h\sum_{i=1}^{n}f(x_i)\f$
  ///
  /// The error of this method is \f$O(h^2)\f$.
  ///
  /// @test GridConversion::TestingObject::Test__MidPointQuadrature tests this 
  /// function by sending linear and quadratic test functions.
  ///
  double MidPointQuadrature(double (*f)(double),double x0,double xn,int n);

};
#endif
