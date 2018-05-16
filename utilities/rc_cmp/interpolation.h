#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <list>
#include <vector>

#include "util.h"
#include "datatypedef.h"
#include "points.h"

using std::list;
using std::vector;


class Interpolate {

    private:

        /**
         * Get a points surrounding a target point in a mesh
         * @param pTarget Target mesh
         * @param pSource Source mesh
         * @param rTargetPoint The number of the point on the target 
         *        mesh that we'd like to find a cell for on the source mesh
         * @param rClosestPoint The number of the point on the source 
         *        mesh closest to the target point
         * @return A list of the point numbers of points surrounding 
         *         the target point.
         */
        static vector<int> getCell( points *pTarget, points *pSource, int rTargetPointNumber, int rClosestPointNumber );

        /**
         * Sort points (rGridPoints) from pSouce into a vector where
         *      0 - Bottom Left
         *      1 - Top Left
         *      2 - Top Right
         *      3 - Bottom Right
         * @param pSource Dataset of the points to be sorted
         * @param rGridPoints 4 points to be sorted as described above
         * @return Vector of point numbers sorted in the order indicated above
         */
        static vector<int> sortPoints( points *pSource, vector<int> rGridPoints );

    public:

        /**
         * Performs bilinear interpolation to find values from source mesh at
         * the target points
         */
        static points *bilinear( points *pTarget, points *pSource, int rTargetVariable, int rSourceVariable );

        //Performs bicubic interpolation to find values from source mesh at 
        //the target points
        static points *bicubic( points *target, points *source, int var_target, int var_source );

};

#endif
