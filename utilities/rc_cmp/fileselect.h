#ifndef _FILESELECT_H_
#define _FILESELECT_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <sstream>
#include <iomanip>

#include "points.h"
#include "file.h"

using std::cout;
using std::cerr;
using std::endl;

class Fileselect {

    public:

        /**
         * Use test functions to determine the filetype
         * and to create a new file object for the input file
         * @param infile Open file stream for the input data file
         * @param outfile Open file stream for the log file
         * @param loud Logging verbosity
         * @param infile_name Name of the data file
         * @param dim Number of dimensions in the input data file
         * @param rFieldMappings Variables mappings that specify which file variables are mapped to which storage variables
         * when the file is read in
         * @param rIndexOrder Order of variables for index creation
         * @param conv_factor Collection of values that are multiplied with corresponding variables from the data file
         * in order to convert between units
         * @param norm_val Collection of values that are multiplied with corresponding variables from the data file
         * in order to undo normalization of the data
         */
        static datafile *detectFiletype(ifstream &infile, ofstream &outfile, bool loud, string infile_name, int dim, 
                                        std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                                        std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val);

        /**
         * Test if a file is a tecplot file
         * @param infile Open file stream
         * @param datapacking Type of datapacking used by the file
         * @return True if the file is a tecplot file, otherwise false
         */
        static bool testTecplot(ifstream &infile, string datapacking);

        /**
         * Print list valid filetypes to stdout
         */
        static void printFiletypes();

    private:

        static const int BUFFER_SIZE = 512;
};

#endif
