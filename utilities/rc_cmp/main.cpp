/* rc_cmp - compares results from simulations and experiments using multiple file
 *          formats and multiple difference metrics
 *
 *
 *  Joseph Kaczmarek, June 2009
 */

#include "main.h"

ComSwitch gCmdSwitches;

int main(int argc, char *argv[]){

    //Parse command line options
    gCmdSwitches.setArguments(argc, argv);
    gCmdSwitches.parseOptions();

    //Open files
    vector<datafile*> files = openAndParseFiles();
    if(files.size() != 2)
        return 1;

    //Perform work on the files
    compareFiles(files);

    //Free resources
    for(int file_num = 0; file_num < files.size(); file_num++){
        if( files[file_num] != NULL )
            delete files[file_num];
    }

    return 0;
}


/*
 * Open the files specified at the command line
 * and parse them
 */
vector<datafile*> openAndParseFiles(){

    vector<datafile*> file_pointers;

    //Create file objects
    datafile *file1 = Fileselect::detectFiletype(gCmdSwitches.getInputFile(1),
                                                 gCmdSwitches.getLogFile(),
                                                 gCmdSwitches.getVerbosity(),
                                                 gCmdSwitches.getFileName(1),
                                                 gCmdSwitches.getFileDimensions(1),
                                                 gCmdSwitches.getFieldMappings(1),
                                                 gCmdSwitches.getIndexOrder(1),
                                                 gCmdSwitches.getConvFactor(1),
                                                 gCmdSwitches.getNormVal(1));

    datafile *file2 = Fileselect::detectFiletype(gCmdSwitches.getInputFile(2),
                                                 gCmdSwitches.getLogFile(),
                                                 gCmdSwitches.getVerbosity(),
                                                 gCmdSwitches.getFileName(2),
                                                 gCmdSwitches.getFileDimensions(2),
                                                 gCmdSwitches.getFieldMappings(2),
                                                 gCmdSwitches.getIndexOrder(2),
                                                 gCmdSwitches.getConvFactor(2),
                                                 gCmdSwitches.getNormVal(2));

    //Check to make sure files have been properly opened
    if(file1 == NULL || file2 == NULL){

        if(file1 != NULL)
            delete file1;

        if(file2 != NULL)
            delete file2;

        file1 = NULL;
        file2 = NULL;

        return file_pointers;

    }

    //Parse file content
    file1->parse();
    file2->parse();

    file_pointers.push_back(file1);
    file_pointers.push_back(file2);

    return file_pointers;
}


/*
 * Compare two files, ouput results
 */
void compareFiles(vector<datafile*> pFiles){

    datafile *file1 = pFiles[0];
    datafile *file2 = pFiles[1];

    if(gCmdSwitches.getVerbosity())
        gCmdSwitches.getLogFile() << "\n";


    //Get comparison list
    std::vector<cmp_map> comparison_list = gCmdSwitches.getComparisonList();

    std::vector<points*> interp_pts;

    for(int i=0; i<comparison_list.size(); i++){

        cmp_map curr_comp = comparison_list[i];

        //Check if previous interp_pts satisfies current interp_pts
        if( i > 0 ){
            cmp_map prev_comp = comparison_list[(i-1)];

            if( !same_elements(curr_comp.file1_partitions, curr_comp.file1_partitions) ||
                curr_comp.file2_partition != prev_comp.file2_partition                 ||
                curr_comp.var1 != prev_comp.var1                                       ||
                curr_comp.var2 != prev_comp.var2                                       ){

                points::delete_all(interp_pts);

                for(int partition_1_num = 0; partition_1_num < curr_comp.file1_partitions.size(); partition_1_num++){
                    points *temp_points;
                    temp_points = Interpolate::bilinear( file1->get_points(curr_comp.file1_partitions[partition_1_num]),
                                                         file2->get_points(curr_comp.file2_partition),
                                                         curr_comp.var1,
                                                         curr_comp.var2 );
                    interp_pts.push_back(temp_points);
                }

            }
        }
        else{

            for(int partition_1_num = 0; partition_1_num < curr_comp.file1_partitions.size(); partition_1_num++){
                points *temp_points;
                temp_points = Interpolate::bilinear( file1->get_points(curr_comp.file1_partitions[partition_1_num]),
                                                     file2->get_points(curr_comp.file2_partition),
                                                     curr_comp.var1,
                                                     curr_comp.var2 );
                interp_pts.push_back(temp_points);
            }
        }

        if(points::contains_null(interp_pts)){

            string out_msg;
            out_msg  = "ERROR   |  Interpolation Failed, ";
            out_msg += "Possible Invalid Partition or Variable Argument: \"";
            out_msg += itoa(curr_comp.metric) + "=";
            
            for(int partition_1_num = 0; partition_1_num < curr_comp.file1_partitions.size(); partition_1_num++){
                
                if(partition_1_num != 0)
                    out_msg += ",";
                
                out_msg += itoa(curr_comp.file1_partitions[partition_1_num]);
            }
            
            out_msg += ":";
            out_msg += itoa(curr_comp.var1) + "/" + itoa(curr_comp.file2_partition) + ":";
            out_msg += itoa(curr_comp.var2) + "\"";

            gCmdSwitches.getLogFile() << out_msg << endl;
            cout << out_msg << endl;

            continue;
        }

        //Construct the outprefix string for the metric
        //Used by the metric to prefix output files
        string metric_outprefix = gCmdSwitches.getOutPrefix();
        metric_outprefix += "_part1_";

        for(int partition_1_num = 0; partition_1_num < curr_comp.file1_partitions.size(); partition_1_num++){
                
            if(partition_1_num != 0)
                metric_outprefix += "_";
                
            metric_outprefix += itoa(curr_comp.file1_partitions[partition_1_num]);
        }

        metric_outprefix += "_part2_" + itoa(curr_comp.file2_partition);

        //Create metric object
        metric *comp = metric_select(curr_comp.metric, retrievePoints(file1, curr_comp.file1_partitions),
                                     interp_pts, curr_comp.var1, 2, //Interpolation only returns the desired point
                                     gCmdSwitches.getLogFile(), metric_outprefix);

        if( gCmdSwitches.isRangeSpecified( ) )
        {
        	BoundingBox range = gCmdSwitches.getRange( );
        	comp ->setRange( &range );
        }

        //Check metric object
        if(comp == NULL){
            gCmdSwitches.getLogFile() << "ERROR   |  " << "Invalid Metric Number: \"";
            gCmdSwitches.getLogFile() <<  curr_comp.metric << "\"" << endl;

            cout << "ERROR   |  " << "Invalid Metric Number: \"";
            cout <<  curr_comp.metric << "\"" << endl;

            continue;
        }

        //Calculate metric
        string res = comp->get_res();

        //Print result
        gCmdSwitches.getLogFile().precision(10);
        gCmdSwitches.getLogFile().width(65);
        string msg_out = "RESULT  |  (" + itoa(curr_comp.file1_partitions[0]) + ":"
                                         + itoa(curr_comp.var1) + " "
                                         + itoa(curr_comp.file2_partition) + ":"
                                         + itoa(curr_comp.var2) + ") "
                                         + comp->get_metric_name() + ": ";

        gCmdSwitches.getLogFile() << std::left << msg_out;
        gCmdSwitches.getLogFile() << std::right << res << "\n";

        cout.precision(10);
        cout.width(65);
        cout << std::left << msg_out;
        cout << std::right << res << "\n";

        delete comp;
    }

    //Cleanup Memory
    points::delete_all(interp_pts);
}


/*
 * Construct a vector of pointers to datasets
 */
std::vector<points*> retrievePoints(datafile* pFile, vector<int> rPartitions){
    
    std::vector<points*> datasets;
    
    for(int partition_record = 0; partition_record < rPartitions.size(); partition_record++){
        if(rPartitions[partition_record] < pFile->get_num_partitions()){
            datasets.push_back( pFile->get_points(rPartitions[partition_record]) );
        }
        else
            datasets.push_back( NULL );
    }

    return datasets;
}


/*
 * Output point values to stdout
 */
void printPointValues(points* pPoints){

    //Print points
    for(int point_num = 0; point_num < pPoints->get_num_points(); point_num++){
        pnt curr = pPoints->get_point(point_num);

        for(int dimension = 0; dimension < curr.size; dimension++)
            cout << curr.vals[dimension] << ", ";
        cout << endl;
    }

    cout << endl << "Printed " << pPoints->get_num_points() << " points" << endl;
}

/*
 * Find the average of the var of the points specified
 */
long double calcAverageValue(points *pPoints, int rVar){

    long double sum = 0.0;

    //Sum the points
    for(int point_num = 0; point_num < pPoints->get_num_points(); point_num++){
        pnt curr = pPoints->get_point(point_num);
        sum += curr.vals[rVar];
    }

    return (sum / (long double)pPoints->get_num_points());
}


/*
 * Find the minimum point value
 */
long double minPoint(points *pPoints, int rVar){
    long double min = std::numeric_limits<long double>::infinity();

    for(int point_num = 0; point_num < pPoints->get_num_points(); point_num++){
        pnt curr = pPoints->get_point(point_num);
        if(curr.vals[rVar] < min)
            min = curr.vals[rVar];
    }

    return min;
}


/*
 * Find the maximum point value
 */
long double maxPoint(points *pPoints, int rVar){
    long double max = -std::numeric_limits<long double>::infinity();

    for(int point_num = 0; point_num < pPoints->get_num_points(); point_num++){
        pnt curr = pPoints->get_point(point_num);
        if(curr.vals[rVar] > max)
            max = curr.vals[rVar];
    }

    return max;
}



/*
 * Print the rc_cmp usage message
 */
void printUsage(string rErrorMessage){

    if(rErrorMessage.compare("") != 0)
        cout << "\n" << rErrorMessage << endl;
    cout << "\n" << gUsage << endl;

    cout << "Comparisons are made at the point locations of input file 1" << endl << endl;

    cout << "Valid Filetypes:\n";
    Fileselect::printFiletypes();
    cout << endl;

    cout << "Metrics:\n";
    print_metrics();
    cout << endl;

}



/* ComSwitch Implementation --------------------------------------------- */

ComSwitch::ComSwitch(){

    //Set default values
    setDefaultVars();
}

ComSwitch::ComSwitch(int argc, char *argv[]){

    //Set default values
    setDefaultVars();

    setArguments(argc, argv);
}

/*
 * Set Arguments
 */
void ComSwitch::setArguments(int rArgc, char *rpArgv[]){
    mArgc = rArgc;
    mpArgv = rpArgv;
}

/*
 * Set default member variables
 */
void ComSwitch::setDefaultVars(){

    //Help message flag
    mHelp = false;

    //Bad command line switches
    mBadSwitch = false;

    //File dimensions
    mInfile1Dimensions = -1;
    mInfile2Dimensions = -1;

    //Log verbosity
    mLoud = false;

    // Restrict search
    restrictToRange = false;
}

/*
 * Parse command options
 */
bool ComSwitch::parseOptions(){

    //Parse Arguments
    readOptions();

    if( mHelp ){
        printUsage(string(""));
        exit(0);
    }

    if( mBadSwitch ){
        printUsage(string(""));
        exit(1);
    }

    if( !initializeFileStreams() )
        exit(1);

    //Check completeness
    if(requiredArgsSet()){
        writeLogHeader();
    }
    else{
        string err_ps("Missing Arguments");
        printUsage(err_ps);

        exit(1);
    }

    return true;

}

ComSwitch::~ComSwitch(){

    //Close files
    mLog.close();
    mInfile1.close();
    mInfile2.close();
}

bool ComSwitch::requiredArgsSet(){

    if( mOutPrefix.compare("") == 0 )
        return false;

    if( !mLog.is_open() )
        return false;

    if( !mInfile1.is_open() || !mInfile2.is_open() )
        return false;

    if( mComparisonList.size() <= 0 )
        return false;

    if( mInfile1Dimensions <= 0 || mInfile2Dimensions <= 0 )
        return false;

    if( mBadSwitch )
        return false;

    return true;
}

/*
 * Get file dimensions
 */
int ComSwitch::getFileDimensions(int file){
    if(file == 1)
        return mInfile1Dimensions;
    else if(file == 2)
        return mInfile2Dimensions;
    else
        return -1;
}

/*
 * Get conversion factors for file
 */
std::vector<adj_map> ComSwitch::getConvFactor(int file){
    if(file == 1){
        return mFile1ConversionFactor;
    }
    else if(file == 2){
        return mFile2ConversionFactor;
    }

    std::vector<adj_map> ret;
    return ret;
}

/*
 * Get normalization values
 */
std::vector<adj_map> ComSwitch::getNormVal(int file){
    if(file == 1){
        return mFile1Norm;
    }
    else if(file == 2){
        return mFile2Norm;
    }

    std::vector<adj_map> ret;
    return ret;
}

/*
 * Get logging verbosity
 */
bool ComSwitch::getVerbosity(){
    return mLoud;
}

/*
 * Get input file
 */
ifstream & ComSwitch::getInputFile(int rFile){

    if(rFile == 1)
        return mInfile1;
    else
        return mInfile2;

}

/*
 * Get file name
 */
string ComSwitch::getFileName(int rFile){
    string ret("");
    if(rFile == 1)
        ret = mInfileName1;
    else if(rFile == 2)
        ret = mInfileName2;

    return ret;
}

/*
 * Get log file
 */
ofstream & ComSwitch::getLogFile(){
    return mLog;
}

/*
 * Get Outprefix
 */
string ComSwitch::getOutPrefix(){
    return mOutPrefix;
}

/*
 * Get comparison list
 */
std::vector<cmp_map> ComSwitch::getComparisonList(){
    return mComparisonList;
}

/*
 * Get field mapping
 */
std::vector< std::vector<int> > ComSwitch::getFieldMappings(int rFile){
    switch( rFile ){

        case 1 :
            return mFieldMappingsFile1;
            break;

        case 2 :
            return mFieldMappingsFile2;
            break;

        default :
            std::vector< std::vector<int> > empty;
            return empty;
    }
}

/*
 * Get index order
 */
std::vector<index_order> ComSwitch::getIndexOrder(int rFile){
    switch( rFile ){
        
        case 1 :
            return mIndexOrderFile1;
            break;

        case 2 :
            return mIndexOrderFile2;
            break;

        default :
            std::vector<index_order> empty;
            return empty;
    }
}

/*
 * Set option string and long opts struct
 */
bool ComSwitch::readOptions(){

    //Short options
    string opt_string("o:1:2:f:i:c:n:m:r:vh");

    //Long options
    static const struct option long_opts[] = {
        {"outfile",     required_argument, NULL, 'o'},
        {"infile1",     required_argument, NULL, '1'},
        {"infile2",     required_argument, NULL, '2'},
        {"fieldmap",    required_argument, NULL, 'f'},
        {"index-order", required_argument, NULL, 'i'},
        {"conv",        required_argument, NULL, 'c'},
        {"norm",        required_argument, NULL, 'n'},
        {"range",	    required_argument, NULL, 'r'},
        {"metric",      required_argument, NULL, 'm'},
        {"verbose",     no_argument,       NULL, 'v'},
        {"help",        no_argument,       NULL, 'h'},
        {NULL,          no_argument,       NULL,  0 }
    };

    return parseValuesFromArguments(opt_string, long_opts);
}

/*
 * Read and parse options from array
 */
bool ComSwitch::parseValuesFromArguments(string rOptString,
                                         const struct option *rLongOpts){


    //Parse options
    int opt = 0;
    opt = getopt_long( mArgc, mpArgv, rOptString.c_str(), rLongOpts, NULL );

    while( opt != -1 ){
        switch( opt ){

            //Outfile
            case 'o': {
                mOutPrefix = optarg;
            } break;

            //File 1
            case '1': {
                std::vector<string> infile_info = str_tokenize(string(optarg), string(":"));

                if(infile_info.size() != 2)
                    break;

                if( check_non_numeric_str( infile_info[1] ) )
                    break;

                mInfileName1 = infile_info[0];
                mInfile1Dimensions = atoi(infile_info[1].c_str());
            } break;

            //File 2
            case '2': {
                std::vector<string> infile_info = str_tokenize(string(optarg), string(":"));

                if(infile_info.size() != 2)
                    break;

                if( check_non_numeric_str( infile_info[1] ) )
                    break;

                mInfileName2 = infile_info[0];
                mInfile2Dimensions = atoi(infile_info[1].c_str());
            } break;

            //Field mapping
            case 'f': {
                std::vector<string> field_map = str_tokenize(string(optarg), string(":/"));

                if(field_map.size() != 3)
                    break;

                //Make sure we only have numerals or "null"
                if( check_non_numeric_str( field_map[0] ) ||
                    check_non_numeric_str( field_map[1] ) ){
                    
                    if(check_non_numeric_str( field_map[2] ) &&
                       field_map[2].compare("null") != 0     ){
                        break;
                    }
                }

                //Pack entry
                vector<int> entry;
                entry.push_back( atoi(field_map[1].c_str()) );

                if(field_map[2].compare("null") != 0)
                    entry.push_back( atoi(field_map[2].c_str()) );
                else
                    entry.push_back( -1 );

                //Add entry to list of mappings
                switch( atoi(field_map[0].c_str()) ){

                    case 1:
                        mFieldMappingsFile1.push_back(entry);
                    break;

                    case 2:
                        mFieldMappingsFile2.push_back(entry);
                        break;

                    default:
                        break;

                }

            } break;

            //Index order
            case 'i': {
                std::vector<string> order = str_tokenize(string(optarg), string(":,"));
                std::vector<int> order_int;

                if(order.size() <= 0)
                    break;

                //Get the file number
                int file_num = 0;
                if(!check_non_numeric_str(order[0])){
                    file_num = atoi(order[0].c_str());

                    if(file_num != 1 && file_num != 2){
                        break;
                    }
                }

                //Get the partition number
                int partition_num = -1;
                if(!check_non_numeric_str(order[1])){
                    partition_num = atoi(order[1].c_str());
                    if(partition_num < 0)
                        break;
                }

                bool non_numeric = false;

                //Loop through remaining values
                for(int token_num = 2; token_num < order.size(); token_num++){
                    
                    if(!check_non_numeric_str(order[token_num])){

                        //Convert to int
                        int value = atoi(order[token_num].c_str());

                        //Check for repeats
                        bool repeat = false;
                        for(int prev_num = 0; prev_num < order_int.size(); prev_num++){
                            if(value == order_int[prev_num]){
                                repeat = true;
                                break;
                            }
                        }

                        //Add value to the list
                        if(!repeat)
                            order_int.push_back(value);
                    }
                    else{
                        non_numeric = true;
                        break;
                    }
                }

                //Discard all values if any are non-numeric
                if(non_numeric)
                    break;
                
                index_order curr_index;
                curr_index.order = order_int;
                curr_index.partition_number = partition_num;

                if(file_num == 1)
                    mIndexOrderFile1.push_back(curr_index);
                else if(file_num == 2)
                    mIndexOrderFile2.push_back(curr_index);
                
            } break;

            //Conversion factor
            case 'c': {
                std::vector<string> conv_info = str_tokenize(string(optarg), string(":"));

                if(conv_info.size() != 3)
                    break;

                adj_map add;
                add.variable = atoi(conv_info[1].c_str());
                add.factor = strtold(conv_info[2].c_str(), NULL);

                int file = atoi(conv_info[0].c_str());

                if(file == 1)
                    mFile1ConversionFactor.push_back(add);
                else if(file == 2)
                    mFile2ConversionFactor.push_back(add);

            } break;

            // Range used for cropping
            case 'r': {

								// Set the flag to indicate that matrics are going to be restricted to the user-supplied range.
								restrictToRange = true;

								std::vector< std::string > conv_info = str_tokenize( std::string( optarg ), std::string( "," ) );

								//
								// The expected parameters are of the form "xmin,ymin,zmin,xmax,ymax,zmax"
								//
								if( conv_info.size( ) != 6 )
									break;

								double xmin = std::atof( conv_info[ 0 ].c_str( ) );
								double ymin = std::atof( conv_info[ 1 ].c_str( ) );
								double zmin = std::atof( conv_info[ 2 ].c_str( ) );
								double xmax = std::atof( conv_info[ 3 ].c_str( ) );
								double ymax = std::atof( conv_info[ 4 ].c_str( ) );
								double zmax = std::atof( conv_info[ 5 ].c_str( ) );

								mRange.setCoordinates( xmin, ymin, zmin, xmax, ymax, zmax );

            }

            //Normalization factor
            case 'n': {
                std::vector<string> norm_info = str_tokenize(string(optarg), string(":"));

                if(norm_info.size() != 3)
                    break;

                adj_map add;
                add.variable = atoi(norm_info[1].c_str());
                add.factor = strtold(norm_info[2].c_str(), NULL);

                int file = atoi(norm_info[0].c_str());

                if(file == 1)
                    mFile1Norm.push_back(add);
                else if(file == 2)
                    mFile2Norm.push_back(add);
            } break;

            //Comparison
            case 'm': {
                std::vector<string> comp = str_tokenize(string(optarg), string("=:/"));

                if(comp.size() != 5)
                    break;

                if( check_non_numeric_str( comp[0] ) ||
                    check_non_numeric_str( comp[2] ) ||
                    check_non_numeric_str( comp[3] ) ||
                    check_non_numeric_str( comp[4] )    ){

                    break;
                }

                bool bad_entry = false;

                //Get file 1's partition list
                std::vector<string> partitions_1 = str_tokenize(string(comp[1]), string(","));               
                std::vector<int> partitions_1_int;
                
                for(int token_num = 0; token_num < partitions_1.size(); token_num++){
                    
                    if(check_non_numeric_str( partitions_1[token_num] )){
                        bad_entry = true;
                        break;
                    }

                    partitions_1_int.push_back( atoi(partitions_1[token_num].c_str()) );
                }

                if(bad_entry)
                    break;

                cmp_map add_cmp;
                add_cmp.metric = atoi(comp[0].c_str());
                add_cmp.file1_partitions = partitions_1_int;
                add_cmp.var1 = atoi(comp[2].c_str());
                add_cmp.file2_partition = atoi(comp[3].c_str());
                add_cmp.var2 = atoi(comp[4].c_str());

                mComparisonList.push_back(add_cmp);

            } break;

            //Logging Level
            case 'v': {
                mLoud = true;
            } break;

            //Help message
            case 'h': {
                mHelp = true;
                return false;
            } break;

            case 0: {
                mBadSwitch = true;
                return false;
            } break;

            case '?': {
                mBadSwitch = true;
                return false;
            } break;
        }

        opt = getopt_long( mArgc, mpArgv, rOptString.c_str(), rLongOpts, NULL );
    }

    return true;
}

/*
 * Initialize files stream
 */
bool ComSwitch::initializeFileStreams(){

    //Intialize Log File
    if(mOutPrefix.compare("") != 0){
        mLog.open((mOutPrefix + ".log").c_str(), ofstream::out);

        if(!mLog.is_open()){
            string err_ps("Invalid Output File Prefix Specified: Directory may not exist");
            printUsage(err_ps);

            return false;
        }
    }

    //Initialize Input File 1
    if(mInfileName1.compare("") != 0){
        mInfile1.open(mInfileName1.c_str(), ifstream::in);

        if(!mInfile1.is_open()){
            string err_ps("Could Not Open Input File 1: ");
            err_ps += mInfileName1;
            printUsage(err_ps);

            return false;
        }
    }

    //Initialize Input File 2
    if(mInfileName2.compare("") != 0){
        mInfile2.open(mInfileName2.c_str(), ifstream::in);

        if(!mInfile2.is_open()){
            string err_ps("Could Not Open Input File 2: ");
            err_ps += mInfileName2;
            printUsage(err_ps);

            return false;
        }
    }

    return true;

}

/*
 * Write the log file header
 */
bool ComSwitch::writeLogHeader(){

    //Write the program name to the log file
    if(mArgc > 0)
        mLog << mpArgv[0] << endl;

    mLog << endl;

    //Write the supplied arguments to the log file
    mLog << "ARGUMENTS: " << endl;

    for(int arg_num=0; arg_num < mArgc; arg_num++){
        mLog << mpArgv[arg_num] << " ";
    }

    mLog << endl << endl;

    //Write input files to log file
    mLog << "INPUT FILES:  " << endl;
    mLog << mInfileName1;
    mLog << ", " << mInfileName2 << endl << endl;

    //Write header for run output
    mLog << "RUN OUTPUT:" << endl << flush;

    return false;
}

