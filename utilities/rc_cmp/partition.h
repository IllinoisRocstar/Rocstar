#ifndef _TECPLOTZONES_H_
#define _TECPLOTZONES_H_

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

#include "util.h"
#include "datatypedef.h"
#include "points.h"

class partition {

    public:
    
        partition();

        virtual ~partition();   

        std::vector<int> get_layout();

        virtual void parse_layout();

        virtual bool parse_data(int num_dep_vars, int num_vars, 
                                std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                                std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

        string get_error();
        string get_status();

        int get_num_points();

        points *get_points();

        pnt get_point(int n);

    protected:

        std::vector<int> p_layout;

        //Messages
        string error;
        string status;

        //DATA POINTS
        points *data_pts;
        bool data_pts_set;

        //Helper Functions
        long double get_adj(std::vector<adj_map> adj, int var);

        //Get field mapping
        int GetFieldMapping(std::vector< std::vector<int> > rFieldMapping, int rSearchField);

};

class tpzone : public partition {
    
    public:

        tpzone(string *zheader, string *zdata);

        virtual ~tpzone();

    protected:

        string *zheader;
        string *zdata;

};

class tpz_ordered : public tpzone {
    
    public: 
        
        tpz_ordered(string *zheader, string *zdata);

        ~tpz_ordered();

        virtual void parse_layout();

        virtual bool parse_data(int num_dep_vars, int num_vars, 
                                std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                                std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

    private:

        void build_point_index(std::vector< std::vector<int> > rFieldMappings);

        int get_index_value(int i, int j, int k, index_order rIndexOrder);

        manual_index_pts *data_pts_local;

};

class tpz_fequad : public tpzone {
    
    public:

        tpz_fequad(string *zheader, string *zdata);
        void build_point_index();
        ~tpz_fequad();

        virtual void parse_layout();

        virtual bool parse_data(int num_dep_vars, int num_vars,
                                std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                                std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

    private:

        connect_points *data_pts_local;

};

#endif
