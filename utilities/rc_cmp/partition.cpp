#include "partition.h"



/* partition -- data partition base class --------------------------------*/

partition::partition(){
    data_pts_set = false;
    return;
}

partition::~partition(){
    if(data_pts_set)
        delete data_pts;
}

std::vector<int> partition::get_layout(){
    if(p_layout.size() == 0){
        parse_layout();
    }

    return p_layout;
}

void partition::parse_layout(){
    return;
}

bool partition::parse_data(int num_dep_vars, int num_vars,
                           std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                           std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val){
    return false;
}

string partition::get_error(){
    string temp = error;
    error = "";

    return temp;
}

string partition::get_status(){
    string temp = status;
    status = "";

    return temp;
}

int partition::get_num_points(){
    return data_pts->get_num_points();
}

points *partition::get_points(){
    return this->data_pts;
}

pnt partition::get_point(int n){
    return data_pts->get_point(n);
}

long double partition::get_adj(std::vector<adj_map> adj, int var){
    for(int i=0; i<adj.size(); i++){
        if(adj[i].variable == var)
            return adj[i].factor;
    }

    return 1.0;
}

int partition::GetFieldMapping(std::vector< std::vector<int> > rFieldMappings, int rSearchField){
    
    for(int map_num = 0; map_num < rFieldMappings.size(); map_num++){
        
        if(rFieldMappings[map_num][0] == rSearchField)
            return rFieldMappings[map_num][1];
    }

    return rSearchField;
}



/* tpzone -- tecplot zone base class -------------------------------------*/

tpzone::tpzone(string *zheader, string *zdata){
    
    this->zheader = zheader;
    this->zdata = zdata;
    this->error = string("");

    return;
}

tpzone::~tpzone(){

    delete zheader;
    
    if(zdata != NULL)       //Delete this after parsing
        delete zdata;
    
    return;
}



/* tpz_ordered -- ordered zone class -------------------------------------*/

tpz_ordered::tpz_ordered(string *zheader, string *zdata) : tpzone(zheader, zdata){
    data_pts_local = new manual_index_pts;
    data_pts = data_pts_local;

    data_pts_set = true;
}

tpz_ordered::~tpz_ordered(){
    return;
}

void tpz_ordered::parse_layout(){

    //Get I value
    std::vector<string> i_val = readCommaExpression( *zheader, string("I") );

    if(i_val.size() == 1){
        if( !check_non_numeric_str(string(i_val[0])) )
            p_layout.push_back(atoi(i_val[0].c_str()));
    }
    else{
        error = string("Unable to read I value from zone header");

        p_layout.erase(p_layout.begin(), p_layout.end());
        return;
    }


    //Get J value
    std::vector<string> j_val = readCommaExpression(*zheader, string("J"));

    if(j_val.size() == 1){
        if( !check_non_numeric_str(string(j_val[0])) )
            p_layout.push_back(atoi(j_val[0].c_str()));
    }
    else{
        error = string("Unable to read J value from zone header");

        p_layout.erase(p_layout.begin(), p_layout.end());
        return;
    }


    //Get K value
    std::vector<string> k_val = readCommaExpression(*zheader, string("K"));

    if(k_val.size() == 1){
        if( !check_non_numeric_str(string(k_val[0])) )
            p_layout.push_back(atoi(k_val[0].c_str()));
    }
    else{
        error = string("Unable to read K value from zone header");

        p_layout.erase(p_layout.begin(), p_layout.end());
        return;
    }

    status =  "I = " + itoa(p_layout[0]);
    status += ", J = " + itoa(p_layout[1]); 
    status += ", K = " + itoa(p_layout[2]);
}

//TODO
//Pass file stream directly to this function to decrease memory usage/copying
bool tpz_ordered::parse_data(int num_dep_vars, int num_vars,
                             std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                             std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val){

    std::stringstream data(std::stringstream::in | std::stringstream::out);
    data << *zdata; 

    //Get the length of each dimension
    int length = 1;
    for(int i=0; i < 3; i++){
        length *= p_layout[i];
    }

    //Initialize storage
    data_pts->set_num_points(length);
    data_pts->set_num_vars(num_vars);
    data_pts->set_num_indep_vars(num_vars - num_dep_vars);

    //Read Data
    string out;

    bool swi = false;

    for(int n=0; n<(num_vars); n++){

        long double conv = get_adj(conv_factor, n);
        long double norm = get_adj(norm_val, n);

        int var_location = GetFieldMapping(rFieldMappings, n);

        for(int z=0; z<p_layout[ 2 ]; z++){
            for(int y=0; y<p_layout[ 1 ]; y++){
                for(int x=0; x<p_layout[ 0 ]; x++){

                    //Read the value
                    data >> out;

                    //Discard data if dimension is mapped to null
                    if(var_location < 0)
                        continue;

                    long double read = strtold(out.c_str(), NULL);

                    //Apply conversion and normalization values
                    read = read * conv * norm;                        

                    //Calculate the storage index
                    int index_val = get_index_value(x, y, z, rIndexOrder);

/*
                    std::cout << index_val << " (" << var_location << ") ";
                    std::cout << ": " << "(" << x << ", " << y << ", " << z << ") > ";
                    std::cout << read << std::endl;
*/

                    //Store the point
                    data_pts_local->set_point_val(index_val, var_location, read);

                }
            }
        }
        //std::cout << std::endl;
    }

    //std::cout << "=====================================================" << std::endl << std::endl;

    //Write Status
    string zs_out("Extracted ");
    zs_out += itoa(length) + " datapoints";
    status = zs_out;

    //Build point index
    build_point_index(rFieldMappings);

    return true;
}

void tpz_ordered::build_point_index(std::vector< std::vector<int> > rFieldMappings){

    int dimension_0 = GetFieldMapping(rFieldMappings, 0);
    int dimension_1 = GetFieldMapping(rFieldMappings, 1);
    int dimension_2 = GetFieldMapping(rFieldMappings, 2);

    std::vector<int> point_layout;
    for(int i=0; i<p_layout.size(); i++){
        point_layout.push_back(p_layout[i]);
    }

    //Set layout
    data_pts_local->set_layout(point_layout);
    
    //Build index
    int n=0;
    for(int i=0; i < p_layout[0]; i++){
        for(int j=0; j < p_layout[1]; j++){
            for(int k=0; k < p_layout[2]; k++){

                std::vector<int> loc;

                loc.push_back(i);
                loc.push_back(j);
                loc.push_back(k);

                data_pts_local->set_index(n, loc);
                pnt curr = get_point(n);
               
                n++;
            }
        }
    }

}



int tpz_ordered::get_index_value(int i, int j, int k, index_order rIndexOrder){

    if(rIndexOrder.order.size() <= 0 || rIndexOrder.order.size() > 3){

        int index_val =   i * p_layout[ 1 ]
                        + j * p_layout[ 2 ]
                        + k;

        return index_val;
    }

    int vars[3] = {i, j, k};
    int layout[3] = {0, 1, 2};

    for(int dim = 0; dim < rIndexOrder.order.size(); dim++){
        switch( rIndexOrder.order[dim] ){
            case 0 : 
                vars[dim] = i;
                layout[dim] = 0;
                break;

            case 1 :
                vars[dim] = j;
                layout[dim] = 1;
                break;

            case 2 :
                vars[dim] = k;
                layout[dim] = 2;
                break;

            default :
                break;
        }
    }
/*
    std::cout << "------------------------" << std::endl << std::endl;

    for(int m=0; m<3; m++){
        std::cout << m << ": " << layout[m] << " | " << p_layout[ layout[m] ] << std::endl;
    }
    std::cout << "===========" << std::endl;
    std::cout << "i | " << i << ": " << vars[0] << std::endl;
    std::cout << "j | " << j << ": " << vars[1] << std::endl;
    std::cout << "k | " << k << ": " << vars[2] << std::endl;
    std::cout << std::endl;
*/
    int index_val =   vars[0] * p_layout[ layout[1] ]
                    + vars[1] * p_layout[ layout[2] ]
                    + vars[2];

    return index_val;
 
}



/* tpz_fequad -- fequad zone class ---------------------------------------*/

tpz_fequad::tpz_fequad(string *zheader, string *zdata) : tpzone(zheader, zdata){
    data_pts_local = new connect_points;
    data_pts = data_pts_local;

    data_pts_set = true;
}

tpz_fequad::~tpz_fequad(){
    return;
}

void tpz_fequad::parse_layout(){

    int nodes = 0;
    int elements = 0;

    //Parse number of nodes
    std::vector<string> nodes_str = readCommaExpression( *zheader, string("Nodes") );

    if(nodes_str.size() == 1){
        if( !check_non_numeric_str(string(nodes_str[0])) )
            nodes = atoi(nodes_str[0].c_str());
    }
    else{
        error = string("Unable to read Nodes value from zone header");

        return;
    }

    //Parse number of elements
    std::vector<string> elements_str = readCommaExpression( *zheader, string("Elements") );

    if(elements_str.size() == 1){
        if( !check_non_numeric_str(string(elements_str[0])) )
            elements = atoi(elements_str[0].c_str());
    }
    else{
        error = string("Unable to read Elements value from zone header");

        return;
    }

    //Set Status
    status = string("Nodes = ") + itoa(nodes) + string(", Elements = ") + itoa(elements);

    p_layout.erase(p_layout.begin(), p_layout.end());
    p_layout.push_back(nodes);

}

//TODO
//Pass file stream directly to this function
bool tpz_fequad::parse_data(int num_dep_vars, int num_vars,
                            std::vector< std::vector<int> > rFieldMappings, index_order rIndexOrder,
                            std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val){

    std::stringstream data(std::stringstream::in | std::stringstream::out);
    data << *zdata; 

    int length = 1;
    for(int i=0; i < (num_vars - num_dep_vars); i++){
        if(i >= p_layout.size())
            break;
        
        length *= p_layout[i];
    }

    //Initialize storage
    data_pts->set_num_points(length);
    data_pts->set_num_vars(num_vars);
    data_pts->set_num_indep_vars(num_vars - num_dep_vars);

    //Read Data
    string out;

    for(int n=0; n<num_vars; n++){

        long double conv = get_adj(conv_factor, n);
        long double norm = get_adj(norm_val, n);

        //Get variable number from mapping
        int var_location = GetFieldMapping(rFieldMappings, n);

        for(int i=0; i<length; i++){
        
            data >> out;

            //Skip saving the point if the dimension is mapped to null
            if(var_location < 0)
                continue;

            long double read = strtold(out.c_str(), NULL);

            //Apply conversion and normalization factors
            read = read * conv * norm;
            
            //Store point value
            data_pts->set_point_val(i, var_location, read);
        }
    }

    //Read Connectivity Information
    int first, sec, third, fourth;
    first = 0;
    sec = 0;
    third = 0;
    fourth = 0;

    //Read 4 values
    data >> first >> sec >> third >> fourth;

    while(!data.eof()){

        data_pts_local->set_connectivity(first, sec);
        data_pts_local->set_connectivity(first, fourth);

        data_pts_local->set_connectivity(sec, third);
        data_pts_local->set_connectivity(sec, first);

        data_pts_local->set_connectivity(third, fourth);
        data_pts_local->set_connectivity(third, sec);

        data_pts_local->set_connectivity(fourth, first);
        data_pts_local->set_connectivity(fourth, third);
        
        first = 0;
        sec = 0;
        third = 0;
        fourth = 0;
        data >> first >> sec >> third >> fourth;
    }
    

    //Set Status
    string zs_out("Extracted ");
    zs_out += itoa(length) + " datapoints";
    status = zs_out;

    return true;
}

