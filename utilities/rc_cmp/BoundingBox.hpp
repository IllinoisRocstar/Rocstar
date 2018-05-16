/**
 * @brief A concrete object that provides the functionality for representing a bounding box
 * geometric primitive and associated operations such as bounding box intersection.
 *
 * @author George Zagaris (gzagaris@illinois.edu)
 * @date 1/30/2009
 */
#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <sstream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <vector>
#include <algorithm>

#include "FloatingPointComparisson.hpp"


class BoundingBox
{
	protected:
		double max[3];
		double min[3];


	public:

		/**
		 * @brief Computes the bounding box intersection between bounding box b1,b2.
		 * @param b1 the bounding box instance to intersect with bounding box b2.
		 * @param b2 the bounding box instance to intersect with bounding box b1.
		 * @return c the bounding box intersection of b1,b2.
		 */
		static BoundingBox intersection( BoundingBox &b1, BoundingBox &b2 )
		{
			 double min[3]; double max[3];

				std::vector< double > xaxis; xaxis.resize( 4 );
				std::vector< double > yaxis; yaxis.resize( 4 );
				std::vector< double > zaxis; zaxis.resize( 4 );

				xaxis[ 0 ] = b1.getMinX( );
				xaxis[ 1 ] = b1.getMaxX( );
				xaxis[ 2 ] = b2.getMinX( );
				xaxis[ 3 ] = b2.getMaxX( );

				yaxis[ 0 ] = b1.getMinY( );
				yaxis[ 1 ] = b1.getMaxY( );
				yaxis[ 2 ] = b2.getMinY( );
				yaxis[ 3 ] = b2.getMaxY( );

				zaxis[ 0 ] = b1.getMinZ( );
				zaxis[ 1 ] = b1.getMaxZ( );
				zaxis[ 2 ] = b2.getMinZ( );
				zaxis[ 3 ] = b2.getMaxZ( );

				std::sort( xaxis.begin( ), xaxis.end( ) );
				std::sort( yaxis.begin( ), yaxis.end( ) );
				std::sort( zaxis.begin( ), zaxis.end( ) );


				min[ 0 ] = xaxis[ 1 ];
				max[ 0 ] = xaxis[ 2 ];

				min[ 1 ] = yaxis[ 1 ];
				max[ 1 ] = yaxis[ 2 ];

				min[ 2 ] = zaxis[ 1 ];
				max[ 2 ] = zaxis[ 2 ];

				xaxis.clear( );
				yaxis.clear( );
				zaxis.clear( );

				return ( BoundingBox( min[ 0 ], min[ 1 ], min[ 2 ],
															max[ 0 ], max[ 1 ], max[ 2 ] ) );

		}

		/**
		 * @brief Same as BoundingBox::getBytesize( ) but does not require
		 * a BoundingBox instance to be allocated.
		 * @return The size of a BoundingBox instance.
		 */
		static size_t Size( )
		{
			return( 6*sizeof( double ) );
		}

		/**
		 * @brief Sets the coordinates of this BoundingBox instance.
		 * @param x min x of the bounding box
		 * @param y min y of the bounding box
		 * @param z min z of the bounding box
		 * @param X max x of the bounding box
		 * @param Y max y of the bounding box
		 * @param Z max z of the bounding box
		 */
		inline void setCoordinates(
			const double x, const double y, const double z,
			const double X, const double Y, const double Z  )
			{
					this ->min[ 0 ] = x;
					this ->min[ 1 ] = y;
					this ->min[ 2 ] = z;
					this ->max[ 0 ] = X;
					this ->max[ 1 ] = Y;
					this ->max[ 2 ] = Z;
			};

		/**
		 * @brief Default Constructor. Constructs a 1x1x1 cube
		 */
		BoundingBox( )
		{
			this ->min[ 0 ]  = 0.0;
			this ->min[ 1 ]  = 0.0;
			this ->min[ 2 ]  = 0.0;
			this ->max[ 0 ]  = 1.0;
			this ->max[ 1 ]  = 1.0;
			this ->max[ 2 ]  = 1.0;
		}

		/**
		 * @brief Custom Constructor
		 * @param x min x of the bounding box
		 * @param y min y of the bounding box
		 * @param z min z of the bounding box
		 * @param X max x of the bounding box
		 * @param Y max y of the bounding box
		 * @param Z max z of the bounding box
		 */
		BoundingBox( const double x, const double y, const double z,
				const double X, const double Y, const double Z )
		{
			this ->setCoordinates( x,y,z, X,Y,Z );
		}

		/**
		 * @brief Destructor
		 */
		~BoundingBox( )
		{

		}

		/**
		 * @return x min x
		 */
		inline double getMinX( ){ return this ->min[ 0 ]; };

		/**
		 * @return y min y
		 */
		inline double getMinY( ){ return this ->min[ 1 ]; };

		/**
		 * @return z min z
		 */
		inline double getMinZ( ){ return this ->min[ 2 ]; };

		/**
		 * @return x max x
		 */
		inline double getMaxX( ){ return this ->max[ 0 ]; };

		/**
		 * @return y max y
		 */
		inline double getMaxY( ){ return this ->max[ 1 ]; };

		/**
		 * @return z max z
		 */
		inline double getMaxZ( ){ return this ->max[ 2 ]; };

		/**
		 * @brief Determines if a point is inside the bounding box.
		 * @param x the x-coordinate of the point in query.
		 * @param y the y-coordinate of the point in query.
		 * @param z the z-coordinate of the point in query.
		 * @return status true if the point (x,y,z) is inside this instance of
		 * bounding box, else false.
		 */
		inline bool hasPoint( const double x, const double y, const double z )
		{

			return(
				( x > this ->getMinX( ) || Numerics::fpointequals( x, this ->getMinX( ) ) ) &&
				( x < this ->getMaxX( ) || Numerics::fpointequals( x, this ->getMaxX( ) ) ) &&
				( y > this ->getMinY( ) || Numerics::fpointequals( y, this ->getMinY( ) ) ) &&
				( y < this ->getMaxY( ) || Numerics::fpointequals( y, this ->getMaxY( ) ) ) &&
				( z > this ->getMinZ( ) || Numerics::fpointequals( z, this ->getMinZ( ) ) ) &&
				( z < this ->getMaxZ( ) || Numerics::fpointequals( z, this ->getMaxZ( ) ) )
			);

		}

};

#endif
