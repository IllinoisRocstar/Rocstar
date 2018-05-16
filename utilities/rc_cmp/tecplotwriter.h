#ifndef _TECPLOTWRITER_H_
#define _TECPLOTWRITER_H_

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
#include <vector>

#include "CImg.h"
#include "file.h"
#include "points.h"
#include "datatypedef.h"

#include "colormaps.h"

using namespace cimg_library;
using std::cout;
using std::cerr;
using std::endl;


/*
 * Base class for objects that write data out to 
 * tecplot data files
 */
class TecplotWriter{

    public:

        /**
         * Default constructor
         */
        TecplotWriter();

        /**
         * Constructor
         * @param rFilename Output filename
         */
        TecplotWriter(string rFilename);

        /**
         * Destructor
         */
        ~TecplotWriter();

        /**
         * Set filename
         * @param rFilename Output filename
         */
        void setOutputFilename(string rFilename);

        /**
         * Get filename
         * @return Output filename
         */
        string getOutputFilename();

        /**
         * Initialize writer
         */
        bool init();

    protected:

        string mFilename;

        ofstream mFileHandle;

};


/*
 * Handles outputing indexed_points data to tecplot 
 * Ordered file
 */
class TecplotOrderedWriter : public TecplotWriter {

    public:
        
        /**
         * Default constructor
         */
        TecplotOrderedWriter();
        
        /**
         * Constructor
         * @param rFilename Output filename
         */
        TecplotOrderedWriter(string rFilename);

        /**
         * Destructor
         */
        virtual ~TecplotOrderedWriter();
        
        /**
         * Write all information to file
         * @param rTitle Dataset Title
         * @param rVariableNames Variable names for header
         * @return True if write is successful, otherwise false
         */
        bool writeFile(string rTitle, std::vector<string> rVariableNames);        

        /**
         * Add ordered partition of points to the file write queue
         * Points are written after WriteFile call
         * Partitions must have the same number of variables
         * @param rTitle Title of partition
         * @param rpPoints Datapoints to write to file
         * @param rPartitionVariable Variable number to write (after the independent variables).
         * If -1, all variables will be written
         */
        bool addOrderedPartition(string rTitle, indexed_points *rpPoints, int rPartitionVariable);

    private:

        /**
         * Write partition to file
         * @param rPartitionPoints Partition datapoints
         * @return True if write is successful, otherwise false
         */
        bool writePartition(int rPartitionNum);

        /**
         * Write the points of the partition specified to file.
         * Points are written in block format.
         * @param rPartitionNum Index of the parition's points in local storage
         * @return True on success, otherwise false
         */
        bool writePoints(int rPartitionNum);

        std::vector<indexed_points *> mPartitionPoints;
        std::vector<string> mPartitionTitles;
        std::vector<int> mPartitionVariable;

        /** Number of points per line */
        static const int POINT_WIDTH = 5;
};

#endif
