#include "file.h"

/* datafile implementation -----------------------------------------------*/

datafile::datafile(ifstream &infile, ofstream &outfile, bool loud, string filename, int dim,
                   std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                   std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val) 
    : infile(infile), outfile(outfile) {

    this->loud = loud;
    this->filename = filename;

    this->conv_factor = conv_factor;
    this->norm_val = norm_val;

    this->mFieldMappings = rFieldMappings;

    this->num_partitions = -1;
    this->dep_vars = 0;

    this->dimensions = dim;

    this->mIndexOrder = rIndexOrder;

    //Seek to beginning of file
    infile.clear();
    infile.seekg(0, std::ios::beg);
}

datafile::~datafile(){

    if(num_partitions > 0){
        for(int i=0; i<num_partitions; i++)
            delete data_part[i];
        delete[] data_part;
    }

}

datafile::datafile(const datafile& d)
    : infile(d.infile), outfile(d.outfile) {

    filename = d.filename;
    loud = d.loud;
    num_partitions = d.num_partitions;
    title = d.title;
    dep_vars = d.dep_vars;

    var_names = d.var_names;

}

string datafile::get_mod_name(){
    return string("Datafile");
}

void datafile::parse(){

    //Print status
    status_out(get_mod_name());

    return;
}

void datafile::error_out(string err_ps){
    outfile << "ERROR   |  File \"" << this->filename << "\":  " << err_ps << std::endl;
    cout << "ERROR   |  File \"" << this->filename << "\":  " << err_ps << std::endl;
}

void datafile::status_out(string st_ps){
    if(loud){
        outfile << "status  |  File \"" << this->filename << "\":  " << st_ps << std::endl;
    }
}

/*
 * Log field mapping data
 */
void datafile::LogFieldMappings(){
    for(int field=0; field<mFieldMappings.size(); field++){
        string st_out = string("Mapping VAR ");
        st_out += itoa(mFieldMappings[field][0]);
        st_out += " -> VAR ";
        st_out += itoa(mFieldMappings[field][1]);

        status_out(st_out);
    }
}

int datafile::partition_layout(int partition, int n){
    if(data_part[partition] != NULL)
        return data_part[partition]->get_layout()[n];
    else
        return 0;
}

/*
 * Get the points of the specified partition
 */
points *datafile::get_points(int partition){
    
    if(partition >= get_num_partitions())
        return NULL;
    
    return data_part[partition]->get_points();
}

/*
 * Get specified point
 */
pnt datafile::get_point(int n){

    //Calculate the correct partition
    int i;
    for(i=0; i<num_partitions; i++){
        if(data_part[i]->get_num_points() > n)
            break;

        n = n - data_part[i]->get_num_points();
    }
    
    //Retrieve point from partition
    return data_part[i]->get_point(n);
}

pnt datafile::get_point(int n, int partition){
    return data_part[partition]->get_point(n);
}

int datafile::get_num_points(){
    
    int ret = 0;
    
    //Get number of points of each partition
    for(int i=0; i<num_partitions; i++){
        if(data_part[i] != NULL)
            ret += data_part[i]->get_num_points();
    }

    return ret;
}

int datafile::get_num_points(int partition){
    
    if(data_part[partition] != NULL)
        return data_part[partition]->get_num_points();
    else
        return 0;
}

int datafile::get_num_partitions(){
    return num_partitions;
}

string datafile::get_title(){
    return title;
}

int datafile::get_num_vars(){
    return var_names.size();
}

int datafile::get_num_dep_vars(){
    return var_names.size() - dimensions;
}

index_order datafile::get_index_order(int partition){
    for(int record_num = 0; record_num < mIndexOrder.size(); record_num++){
        if(mIndexOrder[record_num].partition_number == partition){
            return mIndexOrder[record_num];
        }
    }

    index_order empty;
    return empty;
}



/* tecplot_data - tecplot file parent class --------------------------*/

tecplot_data::~tecplot_data(){  
    return;
}

tecplot_data::tecplot_data(ifstream &infile, ofstream &outfile, bool loud, string filename, int dim,
                           std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                           std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val)
    : datafile(infile, outfile, loud, filename, dim, rFieldMappings, rIndexOrder, conv_factor, norm_val){

    return;
}

string tecplot_data::get_mod_name(){
    return string("Tecplot File");
}

void tecplot_data::seek_to_zone(){

    string tmp("");

    while(!infile.eof() && tmp.compare(string("ZONE")) != 0)
        infile >> tmp;

}

void tecplot_data::read_header(){
    
    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);

    string temp;
    string header;

    infile >> temp;
    
    //Read until a Zone header is reached
    while(temp.compare(string("ZONE")) != 0){
        header += temp;
        header += " ";
        infile >> temp;
    }

    //Parse the file title
    std::vector<string> title_parse = readQuoteExpression(header, string("TITLE"));
    if(title_parse.size() == 1){
        this->title = title_parse[0];
    }

    //Parse variable names
    this->var_names = readQuoteExpression(header, string("VARIABLES"));

    status_out(itoa(this->var_names.size()) + " Variable(s) Detected");

    //Log Field Mapping
    LogFieldMappings();

    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);
    
}

void tecplot_data::zone_status_out(int zone, string zs_out){
    string st_ps("Zone ");
    st_ps = st_ps + itoa(zone) + "  |  " + zs_out;

    status_out(st_ps);
}

void tecplot_data::zone_error_out(int zone, string zs_err){ 
    string err_ps("Zone ");
    err_ps = err_ps + itoa(zone) + "  |  " + zs_err;

    error_out(err_ps);
}

/*
 * Count number of zones in the file from the current location
 */
int tecplot_data::count_zones(){

    int zones = 0;

    //Seek to the next zone until
    //end of file is reached
    seek_to_zone();
    while(!infile.eof()){
        zones++;
        seek_to_zone();
    }

    infile.clear();
    infile.seekg(0, std::ios::beg);

    return zones;

}

void tecplot_data::parse(){
    
    //Get header information
    //  - Sets num_vars
    //  - Sets title
    read_header();

    int zone_count = count_zones();

    if(zone_count == 0){
        string err_ps("No Zones Found");
        error_out(err_ps);
        return;
    }

    status_out(itoa(zone_count) + " Zone(s) Detected");

    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);

    this->num_partitions = zone_count;
    this->data_part = new partition*[this->num_partitions];

    for(int i=0; i<zone_count; i++){

        //Detect zone
        //Create proper zone object
        data_part[i] = zone_detect(i);
    
        //Parse zone header
        if(data_part[i] != NULL){   
            data_part[i]->parse_layout();
            if( data_part[i]->get_error().compare("") != 0)
                zone_error_out(i, data_part[i]->get_error());
            else
                zone_status_out(i, data_part[i]->get_status());
        }
        else{
            string err_ps("Invalid Zone");
            zone_error_out(i, err_ps);
        }
    }

    //Parse zone data
    for(int i=0; i<zone_count; i++){
        
        if(data_part[i] != NULL){
            if(data_part[i]->parse_data(get_num_dep_vars(), get_num_vars(), mFieldMappings, 
                                        get_index_order(i), conv_factor, norm_val))
                zone_status_out(i, data_part[i]->get_status());
            else
                zone_error_out(i, data_part[i]->get_error());
        }
    }

    return;
}

tpzone *tecplot_data::zone_detect(int zone){

    //Retrieve zone header and read zone type
    string *zheader = get_zone_header(zone);
    std::vector<string> ztype = readCommaExpression(*zheader, string("ZONETYPE"));

    tpzone *retzone = NULL;

    if(ztype.size() == 1){

        //Get zone data
        string *zdata = get_zone_nodes(zone);

        //If the zone is an ordered zone, create a new ordered zone object
        if(ztype[0].compare(string("Ordered")) == 0){
            retzone = new tpz_ordered( zheader, zdata );

            string zone_st("Ordered Zone Detected");
            zone_status_out(zone, zone_st);
        }
        //If the zone is FEQuad, create a new FEQuad object
        else if(ztype[0].compare(string("FEQuadrilateral")) == 0){
            retzone = new tpz_fequad( zheader, zdata );

            string zone_st("FEQuadrilateral Zone Detected");
            zone_status_out(zone, zone_st);
        }

    }   

    return retzone;
}

string *tecplot_data::get_zone_header(int zone){

    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);

    //Seek to specified zone
    for(int i=0; i<=zone; i++){
        seek_to_zone();
    }

    string temp;
    infile >> temp;

    string *header = new string;

    //TODO
    //Make finding the end of the header
    //more robust
    while( (temp.find("DT=") == string::npos || temp.find("\"") != string::npos) && !infile.eof()){
        
        *header += temp;
        *header += " ";

        infile >> temp;
    }

    *header += temp;

    //Get the last line
    char c;
    infile.get(c);
    while(c != '\n'){
        *header += c;
        infile.get(c);
    }

    return header;

}

string *tecplot_data::get_zone_nodes(int zone){

    string *ret = new string;

    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);

    //Check valid zone number
    if(zone+1 > num_partitions)
        return NULL;
    
    //Seek to specified zone
    for(int i=0; i<=zone; i++){
        seek_to_zone();
    }

    char line[BUFFER_SIZE];
    string curr_line;

    infile.getline(line, BUFFER_SIZE-1);
    curr_line = line;

    //Skip through the zone header
    while( (curr_line.find("DT=") == string::npos || curr_line.find("\"") != string::npos) 
            && !infile.eof() ){

        infile.getline(line, 511);
        curr_line = line;
    }

    //Read zone point data
    infile.getline(line, 511);
    curr_line = line;

    //Stop when the next zone is reached, or eof
    while( (curr_line.find("ZONE") == string::npos) && !infile.eof() ){

        *ret += curr_line;
        
        infile.getline(line, 511);
        curr_line = line;
    }

    if(infile.eof()){
        infile.clear();
        infile.seekg(0, std::ios::beg);
    }

    return ret;

}
