#include "tecplotwriter.h"

/*
 * Base class for objects that write data out to 
 * tecplot data files
 */

/*
 * Default constructor
 */
TecplotWriter::TecplotWriter(){
    mFilename = "";
}

/*
 * Constructor
 * @param rFilename Output filename
 */
TecplotWriter::TecplotWriter(string rFilename){
    setOutputFilename(rFilename);
}

/*
 * Destructor
 */
TecplotWriter::~TecplotWriter(){

    //Close the file handle
    if(mFileHandle.is_open())
        mFileHandle.close();
    
    return;
}

/*
 * Set filename
 * @param rFilename Output filename
 */
void TecplotWriter::setOutputFilename(string rFilename){
    mFilename = rFilename;
}

/*
 * Get filename
 * @return Output filename
 */
string TecplotWriter::getOutputFilename(){
    return mFilename;
}

/*
 * Initialize writer
 */
bool TecplotWriter::init(){
    mFileHandle.open(mFilename.c_str(), ofstream::out);
    return mFileHandle.is_open();
}



/*
 * Handles outputing indexed_points data to tecplot
 * ordered file
 */

/*
 * Default constructor
 */
TecplotOrderedWriter::TecplotOrderedWriter() : TecplotWriter(){
    return;
}

/*
 * Constructor
 * @param rFilename Output filename
 */
TecplotOrderedWriter::TecplotOrderedWriter(string rFilename) : TecplotWriter(rFilename){
    return;
}

/*
 * Destructor
 */
TecplotOrderedWriter::~TecplotOrderedWriter(){
    return;
}

/*
 * Write all information to file
 * @param rTitle Dataset Title
 * @param rVariableNames Variable names for header
 * @return True if write is successful, otherwise false
 */
bool TecplotOrderedWriter::writeFile(string rTitle, std::vector<string> rVariableNames){
    
    //Check file handle
    if(!mFileHandle.is_open())
        return false;

    //Check to make sure the number of variables match
    if(mPartitionPoints.size() != 0){
        if(rVariableNames.size() != mPartitionPoints[0]->get_num_vars())
            return false;
    }

    mFileHandle << "TITLE = \"" << rTitle << "\"" << endl;
    
    //Write variable names to file
    mFileHandle << "VARIABLES = ";

    for(int var_num = 0; var_num < rVariableNames.size(); var_num++){
        mFileHandle << "\"" << rVariableNames[var_num] << "\"" << endl;
    }

    //Write partitions to file
    for(int partition_num = 0; partition_num < mPartitionPoints.size(); partition_num++){
        writePartition(partition_num);
    }

    return false;
}

/*
 * Add ordered partition of points to the file write queue
 * @param rpPoints Datapoints to write to file
 */
bool TecplotOrderedWriter::addOrderedPartition(string rTitle, indexed_points *rpPoints, int rPartitionVariable){

    if(rpPoints == NULL)
        return false;

    if(mPartitionPoints.size() != 0){

        //Calculate the number of vars in the supplied dataset
        int curr_num_vars = -1;
        if(rPartitionVariable == -1)
            curr_num_vars = rpPoints->get_num_vars();
        else
            curr_num_vars = rpPoints->get_num_indep_vars() + 1;

        //Calculate the number of vars in the previous dataset
        int prev_var = mPartitionVariable.back();
        int prev_num_vars = -1;
        if(prev_var == -1)
            prev_num_vars = mPartitionPoints.back()->get_num_vars();
        else
            prev_num_vars = mPartitionPoints.back()->get_num_indep_vars() + 1;

        //Check variable consistency
        if(curr_num_vars != prev_num_vars)
            return false;
    }

    mPartitionPoints.push_back(rpPoints);
    mPartitionTitles.push_back(rTitle);
    mPartitionVariable.push_back(rPartitionVariable);
    return true;
}

/*
 * Write partition to file
 * @param rPartitionPoints Partition datapoints
 * @return True if write is successful, otherwise false
 */
bool TecplotOrderedWriter::writePartition(int rPartitionNum){
    
    //Sanity check
    if(rPartitionNum >= mPartitionPoints.size() || rPartitionNum > mPartitionTitles.size())
        return false;
    
    indexed_points *curr_points = mPartitionPoints[rPartitionNum];

    //Zone and title
    mFileHandle << "ZONE ";
    mFileHandle << "T = \"" << mPartitionTitles[rPartitionNum] << "\"" << endl;
    
    //Don't know what this is...
    mFileHandle << "STRANDID=0, SOLUTIONTIME=0" << endl;

    if(curr_points->get_dim().size() < 3)
        return false;

    //Write i,j,k values
    int i = curr_points->get_dim()[0];
    int j = curr_points->get_dim()[1];
    int k = curr_points->get_dim()[2];

    mFileHandle << "I=" << i << ", ";
    mFileHandle << "J=" << j << ", ";
    mFileHandle << "K=" << k << ", ";

    //Write other properties
    mFileHandle << "ZONETYPE=Ordered" << endl;
    mFileHandle << "DATAPACKING=BLOCK" << endl;

    mFileHandle << "DT=(SINGLE SINGLE SINGLE SINGLE )" << endl;


    return writePoints(rPartitionNum);
}

/*
 * Write the points of the specified partition.
 * @param rPartitionNum is the number of the local storage of the partition's points
 * @return True on success, otherwise false
 */
bool TecplotOrderedWriter::writePoints(int rPartitionNum){
    
    if(mPartitionPoints.size() <= rPartitionNum)
        return false;

    indexed_points *curr_points = mPartitionPoints[rPartitionNum];
    
    int num_pts = curr_points->get_num_points();
    int vars = curr_points->get_num_vars();

    int values_printed = 0;

    //Loop through and write points
    for(int var_num = 0; var_num < vars; var_num++){
        
        if(mPartitionVariable[rPartitionNum] != -1 && var_num >= curr_points->get_num_indep_vars()){
            if(var_num != mPartitionVariable[rPartitionNum])
                continue;
        }
        
        for(int z = 0; z < curr_points->get_dim()[2]; z++){
            for(int y = 0; y < curr_points->get_dim()[1]; y++){
                for(int x = 0; x < curr_points->get_dim()[0]; x++){
                   
                    //Line break after every POINT_WIDTH values are written
                    if(values_printed == POINT_WIDTH){
                        mFileHandle << endl;
                        values_printed = 0;
                    }

                    //Calculate value index (see tecplot docs)
                    int index =   x * curr_points->get_dim()[1]
                                + y * curr_points->get_dim()[2]
                                + z;

                    mFileHandle << " " << std::uppercase << std::scientific << curr_points->get_point(index).vals[var_num];

                    values_printed++;

                }
            }
        }
    }

    mFileHandle << std::endl;

    return false;

}

