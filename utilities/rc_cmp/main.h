#ifndef _MAIN_H_
#define _MAIN_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <time.h>
#include <getopt.h>
#include <vector>
#include <limits>

#include "file.h"
#include "fileselect.h"
#include "metric.h"
#include "interpolation.h"
#include "datatypedef.h"
#include "util.h"

#include "BoundingBox.hpp"

using std::string;
using std::ofstream;
using std::ifstream;

using namespace std;



class ComSwitch {

    public:

        /**
         * Default Constructor
         */
        ComSwitch();

        /**
         * Constructor with program arguments
         * @param rArgc Number of arguments
         * @param rpArgv Actual arguments
         */
        ComSwitch(int rArgc, char *rpArgv[]);

        /**
         * Free memory, and close open file handles
         */
        ~ComSwitch();

        /**
         * Parse command options, begin logging
         * @return True if options were successfully
         * parsed
         */
        bool parseOptions();

        /**
         * Set command line arguments
         * @param rArgc Number of arguments
         * @param rpArgv Actual arguments
         */
        void setArguments(int rArgc, char *rpArgv[]);

        /**
         * Check if required arguments are set
         * @return True if all required args are set,
         *         otherwise false
         */
        bool requiredArgsSet();

        /**
         * Get file dimensions
         * @param rFile File number (1 or 2)
         * @return Number of independent dimensions
         */
        int getFileDimensions(int rFile);

        /**
         * Get conversion factors
         * @param rFile File number (1 or 2)
         * @return Vector of conversion factors for the file
         */
        std::vector<adj_map> getConvFactor(int rFile);

        /**
         * Get normalization values
         * @param rFile File number (1 or 2)
         * @return Vector of normalization values for the file
         */
        std::vector<adj_map> getNormVal(int rFile);

        /**
         * Get logging verbosity
         * @return True if logging is verbose,
         * otherwise false
         */
        bool getVerbosity();

        /**
         * Get Input File
         * @param rFile File number (1 or 2)
         * @return Reference to ifstream
         */
        ifstream & getInputFile(int rFile);

        /**
         * Get File Name
         * @param rFile File number (1 or 2)
         * @return File name
         */
        string getFileName(int rFile);

        /**
         * Get Log File
         * @return Reference to ofstream logfile
         */
        ofstream & getLogFile();

        /**
         * Get OutPrefix
         * @return OutPrefix string
         */
        string getOutPrefix();

        /**
         * Returns a flag that indicates whether a range was supplied
         * to this instance of CommSwitch.
         * @return
         */
        inline bool isRangeSpecified( )
        {
        	return this ->restrictToRange;
        };

        /**
         * Returns a bounding box corresponding to the range of
         * the region of interest.
         * @return b a bounding box instance.
         * @see BoundingBox
         */
        inline BoundingBox getRange( )
        {
        	return this ->mRange;
        }
        /**
         * Get Comparison List
         * @return Vector comparison list
         */
        std::vector<cmp_map> getComparisonList();

        /**
         * Get the field mapping vector for the
         * specified file
         * @param rFile File number (1 or 2)
         * @return Field mapping vector
         */
        std::vector< std::vector<int> > getFieldMappings(int rFile);

        /**
         * Get the index order for the specified file
         * @param rFile File number (1 or 2)
         * @return Index order vector
         */
        std::vector<index_order> getIndexOrder(int rFile);

    private:

        /**
         * Set default member variables
         */
        void setDefaultVars();

        /**
         * Set option string and long opts struct
         * then parse arguments
         * @return True if options were successfully
         * parsed
         */
        bool readOptions();

        /**
         * Read and parse options from array
         * @param rOptString Short options string
         * @param rLongOpts Long options struct
         * @return True if options were successfully
         * parsed
         */
        bool parseValuesFromArguments(string rOptString,
                                      const struct option *rLongOpts);

        /**
         * Initialize file streams
         * @return True on successful initialization
         * false on failure
         */
        bool initializeFileStreams();

        /**
         * Write the log file header
         * Program Name - switches
         * Input files
         */
        bool writeLogHeader();

        /**
         * Help flag is set from command line options
         */
        bool mHelp;

        int mArgc;     					 /** Number of arguments. 										 															  */
        char **mpArgv; 					 /** The actual command line arguments. 			 																*/
        bool mBadSwitch;  			 /** Bad command line switch detection. 			 															  */
        int mInfile1Dimensions;  /** Infile 1 independent variable dimensions. 															  */
        int mInfile2Dimensions;  /** Infile 2 independent variable dimensions. 															  */

        bool restrictToRange;		 /** Flag used to indicate if the metrics will be restricted to a given range.*/
        BoundingBox mRange; 		 /** Range of region of interest. 						 																*/

        std::vector<adj_map> mFile1ConversionFactor;
        std::vector<adj_map> mFile2ConversionFactor;

        std::vector<adj_map> mFile1Norm;
        std::vector<adj_map> mFile2Norm;

        ofstream mLog;
        string mOutPrefix;

        string mInfileName1;
        string mInfileName2;

        ifstream mInfile1;
        ifstream mInfile2;

        std::vector<cmp_map> mComparisonList;

        std::vector< std::vector<int> > mFieldMappingsFile1;
        std::vector< std::vector<int> > mFieldMappingsFile2;

        std::vector<index_order> mIndexOrderFile1;
        std::vector<index_order> mIndexOrderFile2;

        bool mLoud;

};


string gUsage("Usage:  [--outfile | -o]      <output file prefix>\n\t"
                      "[--infile1 | -1]      <intput file 1 name>:<dimensions>\n\t"
                      "[--infile2 | -2]      <intput file 2 name>:<dimensions>\n\t"
                      "[--fieldmap | -f]     <File (1,2)>:<File Dimension>/<Mapped Dimension or \"null\">\n\t"
                      "[--index-order | -i]  <File (1,2)>:<Partition>:<Dimension 1>,<Dimension 2>,...\n\t"
                      "[--conv | -c]         <File (1,2)>:<Variable>:<Conversion Factor> \n\t"
                      "[--norm | -n]         <File (1,2)>:<Variable>:<Normalization Factor> \n\t"
                      "[--range | -r]        <xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax> \n\t"
                      "[--metric | -m]       <metric number>=<file 1 partition 1>,<file 1 partition 2>,...:<file 1 variable>/\n\t"
                      "                                      <file 2 partition>:<file 2 variable>\n\t"
                      "[--verbose | -v]\n" );

/**
 * Print the rc_cmp usage message
 * @param rErrorString The error message
 */
void printUsage(string rErrorMessage);

/**
 * Open the files specified at the command line
 * and parse them
 * @return Pointers to the files
 */
vector<datafile*> openAndParseFiles();

/**
 * Compare the two files, interpolating the points from
 * the second onto the first mesh. Then, ouput results.
 */
void compareFiles(vector<datafile*> pFiles);

/**
 * Construct a vector of pointers to datasets (point objects)
 * from a list of partition numbers
 */
std::vector<points*> retrievePoints(datafile* pFile, vector<int> rPartitions);

/**
 * Write point values to standard out
 */
void printPointValues(points* pPoints);

/**
 * Calculates the average value of the points
 * @param pPoints Datapoints
 * @param rVar Variable number to average
 * @return Returns the average
 */
long double calcAverageValue(points *pPoints, int rVar);


/**
 * Find the minimum value
 * @param pPoints Points
 * @param rVar Variable number
 */
long double minPoint(points *pPoints, int rVar);


/**
 * Find the maximum value
 * @param pPoints Points
 * @param rVar Variable number
 */
long double maxPoint(points *pPoints, int rVar);


#endif
