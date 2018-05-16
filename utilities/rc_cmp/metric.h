#ifndef _METRIC_H_
#define _METRIC_H_

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

#include "CImg.h"
#include "file.h"
#include "points.h"
#include "datatypedef.h"
#include "tecplotwriter.h"

#include "colormaps.h"

#include "BoundingBox.hpp"

using namespace cimg_library;
using std::cout;
using std::cerr;
using std::endl;



/* 
 * Base class for metrics
 */
class metric {

    public:

        metric(std::vector<points*> datasets1, std::vector<points*> datasets2, 
               int var_ds1, int var_ds2,
               ofstream &outfile, string outprefix);

        void set_datasets1(std::vector<points*> datasets1);
        void set_datasets2(std::vector<points*> datasets2);

        inline void setRange( BoundingBox *rangePtr )
        {

        	if( rangePtr != NULL )
        	{
        		this ->restrictToBoxRange = true;
        		this ->range						  = rangePtr;
        	}

        }

        virtual string get_res();

        virtual string get_metric_name();

    protected:

				bool restrictToBoxRange;
				BoundingBox *range;

        string outprefix;

        std::vector<points*> datasets1;
        std::vector<points*> datasets2;

        int var_ds1;
        int var_ds2;

        ofstream &outfile;

        void error_out(string err_ps);

};



/* Mean square error
 *
 *
 *//*
class mse : public metric {

    public:

        mse(std::vector<points*> datasets1, std::vector<points*> datasets2, 
            int var_ds1, int var_ds2,
            ofstream &outfile, string outprefix);

        string get_metric_name();

        string get_res();

    private:

        long double compare();

        //Make sure the files have the same number of elements
        bool check_card_match();

};
*/


/* Root mean square error
 *
 *
 *//*
class rmse : public metric {

    public:

        rmse(std::vector<points*> datasets1, points *dataset2, int var_ds1, int var_ds2,
                ofstream &outfile, string outprefix);

        string get_res();

        string get_metric_name();

    private:

        long double compare();

        //Make sure the files have the same number of elements
        bool check_card_match();

};
*/


/* Simple Cross-correlation / Template Matching
 *
 *
 *//*
class scc : public metric {

    public:

        scc(std::vector<points*> datasets1, points *dataset2, 
            int var_ds1, int var_ds2,
            ofstream &outfile, string outprefix);

        string get_res();

        string get_metric_name();

    private:

        long double compare();

        //Make sure the files have the same number of elements
        bool check_card_match();

};
*/


/* map Base Class
 *
 *
 *//*
class map : public metric {

    public:

        map(std::vector<points*> datasets1, points *dataset2, 
            int var_ds1, int var_ds2,
            ofstream &outfile, string outprefix);

    protected:

        //For datasets that have a definite start point
        //ie. 0 -> inf
        rgb_color get_color_single_sided(long double pt_val, long double range, long double min);

        //For datasets that have a definite midpoint
        //ie.  -inf <- 0 -> inf
        rgb_color get_color_double_sided(long double pt_val, long double range, long double min,
                                            long double mid_point);
};
*/


/* Difference Map
 *
 *
 *//*
class difmap : public map {

    public:

        difmap(std::vector<indexed_points*> datasets1, indexed_points *dataset2, int var_ds1, int var_ds2,
                ofstream &outfile, string outprefix);

        string get_res();

        virtual string get_metric_name();

    protected:

        virtual string compare();

        //Make sure the files have the same number of elements
        bool check_card_match();

        long double get_max_value(int ds);

        long double get_min_value(int ds);

        long double get_max_dif();

        long double get_min_dif();

        void output_info(string filename, long double max, long double min, int large_ds);
};
*/

/*
class difmap_wkey : public difmap {

    public:

        difmap_wkey(std::vector<indexed_points*> datasets1, indexed_points *dataset2, int var_ds1, int var_ds2,
                    ofstream &outfile, string outprefix);

        string get_metric_name();

    protected:

        virtual string compare();

    private:

        int border_width;
        int key_width;

        rgb_color border_color;
};
*/


/* Color Map
 *
 *
 *//*
class colmap : public map {

    public:

        colmap(std::vector<indexed_points*> datasets1, indexed_points *dataset2, int var_ds1, int var_ds2,
                    ofstream &outfile, string outprefix);

        string get_res();

        string get_metric_name();

    private:

        string compare();

        long double get_max_value(int ds, int var);

        long double get_min_value(int ds, int var);

};
*/



/**
 * Sample Correlation Coefficient
 *//*
class scorco : public metric {

    public:

        scorco(std::vector<points*> datasets1, points *dataset2, int var_ds1, int var_ds2,
                ofstream &outfile, string outprefix);

        string get_res();

        string get_metric_name();

    private:

        long double compare();
*/
        /**
         * Get the average of the var_ds1 variable in dataset 1
         * @return Average value
         *//*
        long double get_average_f1();
*/
        /**
         * Get the average of the var_ds2 variable in dataset2
         * @return Average value
         *//*
        long double get_average_f2();
*/
        //Make sure the files have the same number of elements
/*        bool check_card_match();

};
*/


/**
 * Modelling Efficiency
 *//*
class modef : public metric {

    public:

        modef(std::vector<points*> datasets1, points *dataset2, int var_ds1, int var_ds2,
                ofstream &outfile, string outprefix);

        string get_res();

        string get_metric_name();

    private:

        long double compare();

        long double get_average_f1();
        long double get_average_f2();

        //Make sure the files have the same number of elements
        bool check_card_match();

};
*/


/*
 * Metric Selection
 */

//List metrics
void print_metrics();

//Select metric from switch
metric *metric_select(int m, std::vector<points*> dss1, std::vector<points*> dss2, 
                      int var_ds1, int var_ds2,
                      ofstream &outfile, string outprefix);



/* HOW TO Write your own metric
 *
 * 1. Create a derived class from the "metric" class
 * 2. Write compare() function to do the comparison work
 * 3. Write a get_metric_name() function that returns the name of your metric
 * 4. Add your function to
 *      - print_metrics()
 *      - metric_select()
 */



#endif
