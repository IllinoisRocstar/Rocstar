/*
 * Defines points classes
 * Used for storing datapoints
 */

#ifndef _POINTS_H_
#define _POINTS_H_

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
#include <limits>

#include "util.h"
#include "datatypedef.h"



/*
 * Points Storage Class
 */

class points {

    public:

        points();

        points(int num_points);

        points(const points &p);

        points & operator = (const points &p);

        virtual ~points();

        virtual points *clone();

        //Set Info
        void set_num_points(int n);
        void set_num_vars(int n);
        void set_num_indep_vars(int indep_vars);
        
        int get_num_points();
        int get_num_vars();
        int get_num_indep_vars();

        //Point manipulation
        void set_point_val(int point_number, int var, long double val);
        void clear_dep_vars(int point_number);
        pnt get_point(int n);

        //Get nearest point in specified direction
        virtual int get_closest(long double *coord);
        virtual std::vector<int> get_connected_points(int point_num);

        virtual void sort();

        //Point statistics
        long double largest_val(int dim);
        long double smallest_val(int dim);

    //Static functions

        //Delete all points in the point vector
        static void delete_all(std::vector<points*> points_vector);

        //Returns true if any of the datasets are null
        static bool contains_null(std::vector<points*> points_vector);

    protected:

        int num_points;
        int num_vars;
        int num_indep_vars;

        //POINT STORAGE
        pnt **pt_str;

        //Point Set Map
        //Used to keep track of memory usage
        bitmap pt_set;

};



class connect_points : public points {
    
    public:

        connect_points();

        virtual ~connect_points();

        connect_points & operator = (const connect_points &cp);

        connect_points(const connect_points &cp);

        virtual connect_points *clone();

        virtual int get_closest(long double *coord);

        void set_connectivity(int point_num, int c_point);

        virtual std::vector<int> get_connected_points(int point_num);

        virtual void sort();

    private:

        struct conn_points{
            
            std::vector<int> points;

            conn_points();

            ~conn_points();

            conn_points(const conn_points &cp);
            
            //TODO
            //Move implementation to cpp file
            conn_points & operator = (const conn_points &cp){
                
                if(this != &cp){
                    points = cp.points;
                }
    
                return *this;
            }

            void set_next(int n);
        };

        //Connected Points List
        conn_points *c_points;

};



class indexed_points : public points {
    
    public:
        
        indexed_points();

        virtual ~indexed_points();

        indexed_points(const indexed_points &ip);

        indexed_points & operator = (const indexed_points &ip);

        virtual indexed_points *clone();

        virtual void sort();

        pnt get_indexed_point(std::vector<int> loc);

        virtual int get_closest(long double *coord);
        virtual std::vector<int> get_connected_points(int point_num);

        virtual std::vector<int> get_dim();

    protected:

        VariableDimensionArray index;
        bool index_set;

        std::vector<int> index_layout;  
};



class manual_index_pts : public indexed_points {
    
    public:

        manual_index_pts();

        ~manual_index_pts();

        manual_index_pts(const manual_index_pts &mip);

        manual_index_pts & operator = (const manual_index_pts &mip);

        virtual manual_index_pts *clone();

        void set_layout(std::vector<int> index_layout);

        void set_index(int point_num, std::vector<int> loc);

};


#endif
