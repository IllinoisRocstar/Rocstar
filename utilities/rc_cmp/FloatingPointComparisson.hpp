/**
 * @brief This library provides routines for comparing two floating numbers within a user-supplied tolerance.
 * The algorithm is based on Knuth, The Art of Computer Programming (vol II).
 * @file FloatingPointComparisson.hpp
 * @date Aug 5, 2009
 * @author George Zagaris (gzagaris@illinois.edu)
 */

#ifndef FLOATINGPOINTCOMPARISSON_HPP_
#define FLOATINGPOINTCOMPARISSON_HPP_

#include <limits>

namespace Numerics {


			 /**
			  * @brief Returns the absolute value of the value a.
			  * @param a the value whose absolute value is queried
			  * @return abs the absolute value of a abs( a ).
			  * @post abs > 0
			  */
			 inline double AbsoluteValue( double a )
			 {
				 return( ( ( a < static_cast< double >( 0 ) )? -a : a ) );
			 }

			 /**
				* @brief Performs a safe floating point division f1/f2
				* @param f1 numerator
				* @param f2 denominator
				* @return f1/f2
				* @pre f1 >= 0
				* @pre f2 >= 0
				* @post 0 <= f1/f2 <= std::numeric_limits<double>::max( )
				*/
				 inline double ___safe_fpt_division( const double f1, const double f2 )
				 {
					#ifdef ASSERT_ON
					 assert( f1 >= static_cast< double >( 0 ) );
					 assert( f2 >= static_cast< double >( 0 ) );
					#endif

					 /* Avoid overflow */
					 if( f2 < static_cast< double >( 1 ) && f1 > f2*std::numeric_limits< double >::max( ) )
						 return std::numeric_limits<double>::max( );


					 /* Avoid underflow */
					 if( f1 == static_cast< double >( 0 ) ||
					 f2 > static_cast< double >( 1 ) && f1 < f2*std::numeric_limits< double >::min( ) )
						 return static_cast< double >( 0 );

					 return( f1/f2 );
				 }

			 /**
				* @brief Checks if two floating numbers are "nearly" equal.
				* @param a first floating number to check
				* @param b second floating number to check
				* @return s true if a =~ b, else false.
				* @note This implementation is based on Knuth, The Art of Computer Programming (vol II).
				*/
				 inline bool fpointequals( const double a, const double b, double TOL=1e-9 )
				 {
					 double diff      = a-b;
					 double adiff 	  = AbsoluteValue( diff );
					 double d1		  = ___safe_fpt_division( adiff, AbsoluteValue( a ) );
					 double d2		  = ___safe_fpt_division( adiff, AbsoluteValue( b ) );
					 bool predicate1  = (d1 <= TOL);
					 bool predicate2  = (d2 <= TOL);

					 if( (predicate1 || predicate2) )
						 return true;
					 return false;
				 }


}
#endif /* FLOATINGPOINTCOMPARISSON_HPP_ */
