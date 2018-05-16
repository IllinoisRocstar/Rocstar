#ifndef _FILE_H_
#define _FILE_H_

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

//Helper classes
#include "partition.h"

using std::string;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::cerr;


/** 
 * Base class for file parsing.
 */
class datafile {

    public:

        /**
         * Construct data file object using specified information
         * @param infile Input file stream
         * @param outfile Log file stream
         * @param loud Log verbosity
         * @param filename Name of the input file
         * @param dim Number of dimensions in the input file
         * @param rFieldMappings Mappings of file variables to internal representation variables
         * @param rIndexOrder The order of variables for index creation (for each partition)
         * @param conv_factor Factors multiplied with corresponding variables in order to convert values
         * with difference units
         * @param norm_val Values multiplied with corresponding variables in order to reverse normalization
         */
        datafile(ifstream &infile, ofstream &outfile, bool loud, string filename, int dim,
                 std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                 std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

        /**
         * Get class name
         * @return Name of the class (datafile)
         */
        virtual string get_mod_name();

        /**
         * Free used data, delete the object
         */
        virtual ~datafile();

        /**
         * Deep copy constructor
         * @param d Datafile to copy
         */
        datafile(const datafile& d);

        /**
         * Parse the data file. Non functional, meant to be overriden by subclasses.
         */
        virtual void parse();   

        /**
         * Get the size of n dimension of the specified partition
         * @param partition Partition number
         * @param n Dimension number
         * @return Size of n in partition
         */
        int partition_layout(int partition, int n);

        /**
         * Get a partitions' points
         * @param partition Partition number
         * @return Points of the specified partition
         */
        points *get_points(int partition);

        /**
         * Get the nth point of the file
         * @param n Point to retrieve
         * @return nth point of the entire file
         */
        pnt get_point(int n);

        /**
         * Get the nth point of the specified partition
         * @param n Point number
         * @param partition Partition number
         * @return nth point of partition
         */
        pnt get_point(int n, int partition);

        /**
         * Get the number of points in the file
         * @return Number of points in file
         */
        int get_num_points();

        /**
         * Get the number of points in the specified partition
         * @param partition Partition number
         * @return Number of points in partition
         */
        int get_num_points(int partition);

        /**
         * Get the number of variables in the file
         * @return Number of variables
         */
        int get_num_vars();

        /**
         * Get file title if specified in file header
         * @return File title
         */
        string get_title();

        /**
         * Get number of partition in the file
         * @return Number of partitions
         */
        int get_num_partitions();

    protected:  

        string filename;

        //File streams for reading and loggin
        ifstream &infile;
        ofstream &outfile;
        bool loud;  

        /**
         * Log the string specified as an error event
         * @param err_ps Error message
         */
        void error_out(string err_ps);

        /**
         * Log the string specified as a status event
         * @param st_ps Status message
         */
        void status_out(string st_ps);

        /**
         * Log the field mapping data from the command line
         * and stored in mFieldMappings
         */
        void LogFieldMappings();

        //PARTITIONS
        partition **data_part;

        //File info
        int num_partitions;
        std::vector<string> var_names;
        string title;   

        /**
         * Get number of dependent variables
         * @return Number of dependent variables
         */
        int get_num_dep_vars();
        
        //Number of dimensions
        //aka independent variables
        int dimensions;

        //Conversion information
        std::vector<adj_map> conv_factor;
        std::vector<adj_map> norm_val;

        //Field mappings
        std::vector< std::vector<int> > mFieldMappings;

    private:
        
        //Dependent Variables
        int dep_vars;

    protected:

        index_order get_index_order(int partition);

        std::vector<index_order> mIndexOrder;
};



/** 
 * Base class for tecplot files
 * Implements useful functions for tecplot file parsing
 */
class tecplot_data : public datafile {

    public:

        /**
         * Get the name of the class (tecplot_data)
         * @return Name of the class
         */
        virtual string get_mod_name();

        /**
         * Constructor that initializes values
         */
        tecplot_data(ifstream &infile, ofstream &outfile, bool loud, string filename, int dim,
                     std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                     std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

        /**
         * Destructor
         */
        virtual ~tecplot_data();

        /**
         * Parse the tecplot file. Read the header using read_header(), iterate through each zone
         * and retrieve the zone header and zone body. Create a new partition object and have that
         * object parse the zone. Log and error output.
         */
        virtual void parse();

    protected:
    
        /**
         * Read the tecplot header from the file and retrieve the title, variable names,
         * and number of variables. Log the information as a status event.
         */
        void read_header(); 

        /**
         * Count the number of zones in the tecplot file, starting from the current location
         * of the file stream. For count of all zones, seek to beginning of file then call.
         * @return Number of zones from the current location.
         */
        int count_zones();

        /**
         * Seek to the next zone in the file. Starting position is the current location of
         * the infile object.
         */
        void seek_to_zone();

        /**
         * Read the zone header
         * @return String of zone header
         */
        string *get_zone_header(int zone);
        
        /**
         * Retrieve zone point data from the zone body.
         * @return String of zone data points.
         */
        string *get_zone_nodes(int zone);

        //Logging
        /**
         * Log an error message with appended zone information so the error can be located.
         * @param zone Zone number affected.
         * @param zs_err Error message.
         */
        void zone_error_out(int zone, string zs_err);
        
        /**
         * Log a status message with appended zone information.
         * @param zone Zone number affected.
         * @param zs_out Status message.
         */
        void zone_status_out(int zone, string zs_out);

        /**
         * Read the header information for the specified zone, then create a new zone object
         * based on the type of zone detected.
         * @param zone Zone number
         * @return Zone object
         */
        tpzone *zone_detect(int zone);

    private:

        static const int BUFFER_SIZE = 512;

};

#endif
