#include "interpolation.h"


/*
 * Get cell (rectangle) around target point
 */
vector<int> Interpolate::getCell( points *pTarget, points *pSource, int rTargetPointNumber, int rClosestPointNumber ){

    //Retrieve the actual points
    pnt tgt_pnt = pTarget->get_point( rTargetPointNumber );
    pnt closest_pnt = pSource->get_point( rClosestPointNumber );    

    //Initialize the storage vector
    //Stores point indices of the source points that surround the target point 
    vector<int> grid_points;
    
    grid_points.push_back(-1);
    grid_points.push_back(-1);
    grid_points.push_back(-1);
    grid_points.push_back(-1);

    //Get points closest to the source point that's nearest the target point
    std::vector<int> near_pts = pSource->get_connected_points(rClosestPointNumber);

    //Find triangles using the closest points
    std::list<int*> tri_str;

    int tris = 0;
    for(int j=0; j<near_pts.size(); j++){
        for(int k=j; k<near_pts.size(); k++){
            if(k == j)
                continue;

            //Test if the target point is in the triangle specified
            //by the two current points and the closest point
            if(pnt_in_tri(tgt_pnt, pSource->get_point( near_pts[j] ),
                                   pSource->get_point( near_pts[k] ),
                                   closest_pnt ) ){
                
                //If so, store the two points
                int *t_add = new int[2];
                t_add[0] = near_pts[j];
                t_add[1] = near_pts[k];
                    
                tri_str.push_back(t_add);

                tris++;
            }
        }
    }

    //Check for two common points
    //There should be two triangles (which form a rectangle around the target point)
    //Find the two points common between the triangles
    for(list<int*>::iterator it = tri_str.begin(); it != tri_str.end(); it++){
    
        std::vector<int> near_pts_x = pSource->get_connected_points( (*it)[0] );
        std::vector<int> near_pts_y = pSource->get_connected_points( (*it)[1] );

        int matched_pts = 0;

        int mp_str[2];

        for(int j=0; j<near_pts_x.size(); j++){
            for(int k=0; k<near_pts_y.size(); k++){

                if( near_pts_x[j] == near_pts_y[k] ){                   
                    if(matched_pts < 2)
                        mp_str[matched_pts] = near_pts_x[j];
                        
                    matched_pts++;
                }
            }
        }

        if(matched_pts == 2 
                        && mp_str[0] != mp_str[1] //){
                        && mp_str[0] != (*it)[0]            //Detect Triangles
                        && mp_str[0] != (*it)[1]
                        && mp_str[1] != (*it)[0]
                        && mp_str[1] != (*it)[1] ){

            grid_points[0] = (*it)[0];
            grid_points[1] = (*it)[1];
            grid_points[2] = mp_str[0];
            grid_points[3] = mp_str[1];

            break;
        }
    }

    for(list<int*>::iterator it = tri_str.begin(); it != tri_str.end(); it++){
        delete[] (*it);
    }

    //Return final 4 points
    return grid_points;
}

/*
 * Sort the given points from the given dataset in clockwise order from bottom left
 */
vector<int> Interpolate::sortPoints( points *pSource, vector<int> rGridPoints ){
   
    vector< int > sorted_points;

    //Check input
    if(pSource == NULL || rGridPoints.size() != 4)
        return sorted_points;
 
    //Sort independent var 0 values
    while( rGridPoints.size() > 0 ){
        
        //Pop a point
        int curr_point = rGridPoints.back();
        rGridPoints.pop_back();

        pnt to_insert = pSource->get_point(curr_point);

        vector<int>::iterator it;

        for(it = sorted_points.begin(); it < sorted_points.end(); it++){
            pnt curr_sorted = pSource->get_point(*it);

            if(to_insert.vals[0] < curr_sorted.vals[0]){
                sorted_points.insert(it, curr_point);
                break;
            }
        }

        //Insert the point at the end if it hasn't already been inserted
        if(it == sorted_points.end())
            sorted_points.push_back(curr_point);
    
    }

    //Sort independent var 1 values
    if(pSource->get_point(sorted_points[0]).vals[1] > pSource->get_point(sorted_points[1]).vals[1]){
        int zero = sorted_points[0];
        sorted_points[0] = sorted_points[1];
        sorted_points[1] = zero;
    }

    if(pSource->get_point(sorted_points[2]).vals[1] < pSource->get_point(sorted_points[3]).vals[1]){
        int two = sorted_points[2];
        sorted_points[2] = sorted_points[3];
        sorted_points[3] = two;
    }

    return sorted_points;
}

/*
 * Return interpolated values of dataset "source" at the points of dataset "target"
 */
points * Interpolate::bilinear( points *pTarget, points *pSource, int rTargetVariable, int rSourceVariable ){

    //Sanity checks
    if( pTarget == NULL || pSource == NULL )
        return NULL;

    if( rTargetVariable >= pTarget->get_num_vars() || rSourceVariable >= pSource->get_num_vars() )
        return NULL;

    //Bilinear interpolation only works with 
    //plane data (2 dimensions, 3rd dimension = value)
    if( pTarget->get_num_indep_vars() != 2 ||
        pSource->get_num_indep_vars() != 2 ){
        
        return NULL;
    }

    //Create data structure for return value
    points *interp = pTarget->clone();

    int dropped_points = 0;

    //Loop over all points
    for(int i=0; i<pTarget->get_num_points(); i++){
        
        //Get point coordinates
        long double *coord = new long double[pTarget->get_num_indep_vars()];
        pnt tgt_pnt = pTarget->get_point(i);

        for(int m=0; m<pTarget->get_num_indep_vars(); m++){
            coord[m] = tgt_pnt.vals[m];
        }


        //Get nearest point in the source mesh
        int closest = pSource->get_closest(coord);
        delete[] coord;

        pnt closest_pnt = pSource->get_point(closest);

        //Get cell
        vector<int> grid_points = getCell( pTarget, pSource, i, closest );  

        //TODO
        //Handle case where less than 4 points can be established
        //IE: the point is off the grid
        if( grid_points[0] < 0 || grid_points[1] < 0 || 
            grid_points[2] < 0 || grid_points[3] < 0 ){

            dropped_points++;   
            interp->set_point_val(i, rTargetVariable, nan());
            continue;
        }

        
        //Extract specific grid information
        vector<int> sorted_points = sortPoints(pSource, grid_points);

        long double bot_right_0 = pSource->get_point(sorted_points[3]).vals[0];
        long double bot_left_0 = pSource->get_point(sorted_points[0]).vals[0];
        long double top_left_1 = pSource->get_point(sorted_points[1]).vals[1];
        long double bot_left_1 = pSource->get_point(sorted_points[0]).vals[1];

        long double source_vals[4];
        
        //0 = Bottom Left Point
        //1 = Top Left Point
        //2 = Top Right Point
        //3 = Bottom Right Point

        for(int sp_num = 0; sp_num < sorted_points.size(); sp_num++){
            source_vals[sp_num] = pSource->get_point(sorted_points[sp_num]).vals[rSourceVariable];    
        }

        //Calculate square distances
        long double dist_0 = bot_right_0 - bot_left_0;
        long double dist_1 = top_left_1  - bot_left_1;

        //Correct locations
        long double tgt_0 = tgt_pnt.vals[0];
        long double tgt_1 = tgt_pnt.vals[1];

        long double t = (tgt_0 - bot_left_0) / dist_0;
        long double u = (tgt_1 - bot_left_1) / dist_1;

        long double interp_val  = (1-t)*(1-u)*source_vals[0] + t*(1-u)*source_vals[1];
                    interp_val += t*u*source_vals[2] + (1-t)*u*source_vals[3];

        //Store interpolated value
        interp->set_point_val(i, rTargetVariable, interp_val);
    }

    return interp;
}

