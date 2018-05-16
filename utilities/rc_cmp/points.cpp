#include "points.h"



/* points -- base class for point storage --------------------------------*/

points::points(){
    this->num_points = 0;
    this->num_vars = 0;
    this->num_indep_vars = 0;
    this->pt_str = NULL;
    return;
}

points::~points(){

    if(num_points > 0){

        for(int i=0; i<num_points; i++){
            if(pt_set.check_value(i)){
                delete pt_str[i];
            }
        }
        delete[] pt_str;
    }

    return;
}

points::points(const points &p){

    num_points = p.num_points;
    num_vars = p.num_vars;
    num_indep_vars = p.num_indep_vars;

    if(num_points > 0)
        pt_str = new pnt*[num_points];

    pt_set = p.pt_set;

    for(int i=0; i<num_points; i++){
        if(pt_set.check_value(i)){

            pt_str[i] = new pnt( *(p.pt_str[i]) );
            
        }
    }

}

points & points::operator = (const points &p){

    if(this != &p){

        for(int i=0; i<num_points; i++){
            if(pt_set.check_value(i))
                delete pt_str[i];
        }
    
        if(pt_str != NULL)
            delete[] pt_str;

        num_points = p.num_points;
        num_vars = p.num_vars;
        num_indep_vars = p.num_indep_vars;

        pt_set = p.pt_set;
        pt_str = new pnt*[num_points];

        for(int i=0; i<num_points; i++){
            if(pt_set.check_value(i)){
                pt_str[i] = new pnt( *(p.pt_str[i]) );  
            }
        }
    }

    return *this;
}

points *points::clone(){
    return new points(*this);
}

void points::set_num_points(int n){

    if(num_points <= 0){
        pt_str = new pnt*[n];
    }
    else{
        pnt **temp;
        temp = new pnt*[n];

        for(int i=0; i<get_num_points(); i++){
            if(pt_set.check_value(i)){
                temp[i] = pt_str[i];
                delete pt_str[i];
            }
        }

        delete[] pt_str;
        pt_str = temp;  
    }

    this->num_points = n;
    pt_set.set_size(n);

}

void points::set_num_vars(int n){
    num_vars = n;
}

void points::set_num_indep_vars(int indep_vars){
    num_indep_vars = indep_vars;
}

int points::get_num_points(){
    return num_points;
}

int points::get_num_vars(){
    return num_vars;
}

int points::get_num_indep_vars(){
    return num_indep_vars;
}

void points::set_point_val(int point_number, int var, long double val){

    //Check point index
    if(point_number > num_points)
        return;
    
    //Check if the pnt object has already been created
    if(pt_set.check_value(point_number)){

        //Check size of point object
        if(pt_str[point_number]->size != num_vars && num_vars > pt_str[point_number]->size){
            pnt *temp = new pnt;

            temp->vals = new long double[num_vars];
            temp->size = num_vars;

            for(int m=0; m<pt_str[point_number]->size; m++){
                temp->vals[m] = pt_str[point_number]->vals[m];
            }

            delete[] pt_str[point_number]->vals;
            delete[] pt_str[point_number];

            pt_str[point_number] = temp;
        }

        if(pt_str[point_number]->size > var)
            pt_str[point_number]->vals[var] = val;
    }
    else{
        pt_str[point_number] = new pnt;
        pt_str[point_number]->vals = new long double[num_vars];
        pt_str[point_number]->size = num_vars;

        pt_set.set_value(point_number);

        if(pt_str[point_number]->size > var)
            pt_str[point_number]->vals[var] = val;
        
    }
}

void points::clear_dep_vars(int point_number){

    long double *temp = new long double[num_indep_vars];

    for(int i=0; i<num_indep_vars; i++){
        temp[i] = pt_str[point_number]->vals[i];
    }

    delete[] pt_str[point_number]->vals;
    pt_str[point_number]->vals = temp;
}

pnt points::get_point(int n){   
    if(pt_set.check_value(n))
        return *(pt_str[n]);
    else{
        pnt blank;
        return blank;
    }
}

int points::get_closest(long double *coord){
    return -1;
}

std::vector<int> points::get_connected_points(int point_num){
    std::vector<int> ret;
    return ret;
}

void points::sort(){
    return;
}

long double points::largest_val(int dim){
    long double largest = -std::numeric_limits<long double>::infinity();

    for(int i=0; i<num_points; i++){
        if(pt_str[i]->vals[dim] > largest)
            largest = pt_str[i]->vals[dim];
    }

    return largest;
}

long double points::smallest_val(int dim){
    long double smallest = std::numeric_limits<long double>::infinity();

    for(int i=0; i<num_points; i++){
        if(pt_str[i]->vals[dim] < smallest)
            smallest = pt_str[i]->vals[dim];
    }

    return smallest;
}

void points::delete_all(std::vector<points*> points_vector){
    for(int set_num = 0; set_num < points_vector.size(); set_num++){
        if(points_vector[set_num] != NULL){
            delete points_vector[set_num];
        }
    }

    points_vector.erase(points_vector.begin(), points_vector.end());
}

bool points::contains_null(std::vector<points*> points_vector){
    for(int set_num = 0; set_num < points_vector.size(); set_num++){
        if(points_vector[set_num] == NULL)
            return true;
    }

    return false;
}



/* connect_points -- class for points with connectivity info -------------*/

connect_points::connect_points(){
    c_points = NULL;
}

connect_points::~connect_points(){
    delete[] c_points;
}

connect_points::connect_points(const connect_points &cp) : points(cp) {

    if(cp.num_points > 0){
    
        c_points = new conn_points[num_points];

        for(int i=0; i<num_points; i++){
            c_points[i] = cp.c_points[i];
        }
    }
}

connect_points & connect_points::operator = (const connect_points &cp){
    if(this != &cp){
        
        int curr_pts = num_points;

        points::operator = (cp);

        delete[] c_points;

        c_points = new conn_points[cp.num_points];

        for(int i=0; i<cp.num_points; i++){
            c_points[i] = cp.c_points[i];
        }
    }

    return *this;
}

connect_points *connect_points::clone(){
    return new connect_points(*this);
}

int connect_points::get_closest(long double *coord){

    long double low_dist = std::numeric_limits<long double>::infinity();
    int low_point = -1;

    //Iterate through all points and find the point that's closest
    //to the specified point
    for(int i=0; i<num_points; i++){
        
        long double *d_calc = new long double[num_indep_vars];
        
            for(int m=0; m<num_indep_vars; m++){
                d_calc[m] = pt_str[i]->vals[m];
            }

        long double distance = dist(coord, d_calc, num_indep_vars);

        if( distance < low_dist ){
            low_dist = distance;
            low_point = i;
        }

        delete[] d_calc;
    }

    return low_point;
}

void connect_points::set_connectivity(int point_num, int c_point){

    //Point numbers start at 1
    //Subtract 1 for array storage (starting at 0)
    point_num = point_num - 1;

    if(c_points == NULL && num_points > 0){
        c_points = new conn_points[num_points];
    }

    bool dup = false;

    for(int i=0; i<c_points[point_num].points.size(); i++){
        if(c_points[point_num].points[i] == c_point){
            dup = true;
            break;
        }
    }

    if(!dup)
        c_points[point_num].set_next(c_point);
}

std::vector<int> connect_points::get_connected_points(int point_num){

    std::vector<int> ret;
    
    //Point numbers start at 1
    //Subtract 1 for array storage (starting at 0)
    point_num = point_num - 1;

    if(point_num < 0)
        return ret;

    for(int i=0; i<c_points[point_num].points.size(); i++){
        ret.push_back(c_points[point_num].points[i]);
    }

    return ret;

}

void connect_points::sort(){
    return;
}

connect_points::conn_points::conn_points(){
}

connect_points::conn_points::~conn_points(){
}

connect_points::conn_points::conn_points(const conn_points &cp) {
    points = cp.points;
}

//TODO
//Move implementation here
/*
conn_points & connect_points::conn_points::operator = (const conn_points &cp){
    
}
*/

void connect_points::conn_points::set_next(int n){
    points.push_back(n);
}



/* indexed_points -- class for indexed point storage ---------------------*/

indexed_points::indexed_points() {
    index_set = false;
    return;
}

indexed_points::~indexed_points(){
    return;
}

indexed_points::indexed_points(const indexed_points &ip) : points(ip) {

    index_set = ip.index_set;
    index_layout = ip.index_layout;

    index = ip.index;
}

indexed_points & indexed_points::operator = (const indexed_points &ip) {
    
    if(this != &ip){

        points::operator = (ip);

        index_set = ip.index_set;
        index_layout = ip.index_layout;

        index = ip.index;
    }

    return *this;
}

indexed_points *indexed_points::clone(){
    return new indexed_points(*this);
}

void indexed_points::sort(){
    return;
}

pnt indexed_points::get_indexed_point(std::vector<int> loc){    

    int index_loc = index.getValue(loc);

    return *(pt_str[ index_loc ]);
}

int indexed_points::get_closest(long double *coord){

    //TODO
    //Data is ordered;
    //Improve search algorithm

    long double low_dist = std::numeric_limits<long double>::infinity();
    int low_point = -1;

    //Iterate through all points and find the point in the dataset
    //that's closest to the specified point
    for(int i=0; i<num_points; i++){
        
        long double *d_calc = new long double[num_indep_vars];
        
            for(int m=0; m<num_indep_vars; m++){
                d_calc[m] = pt_str[i]->vals[m];
            }

        long double distance = dist(coord, d_calc, num_indep_vars);

        if( distance < low_dist ){
            low_dist = distance;
            low_point = i;
        }

        delete[] d_calc;
    }

    return low_point;
}

std::vector<int> indexed_points::get_connected_points(int point_num){
    
    std::vector<int> coords = index.translate_index_to_coords(point_num);

    std::vector<int> ret;

    //Iterate through each independent variable
    for(int i=0; i<coords.size(); i++){
        std::vector<int> loc;

        //Positive Direction
        for(int j=0; j<coords.size(); j++){
            if(j == i)
                loc.push_back(coords[j] + 1);
            else
                loc.push_back(coords[j]);
        }
        int ptup = index.translate_coord_to_index(loc);

        //Erase coordinates
        loc.erase(loc.begin(), loc.end());

        //Negative Direction
        for(int j=0; j<coords.size(); j++){
            if(j == i)
                loc.push_back(coords[j] - 1);
            else
                loc.push_back(coords[j]);
        }
        int ptdown = index.translate_coord_to_index(loc);

        if(ptup > 0)
            ret.push_back(ptup);
        if(ptdown > 0)
            ret.push_back(ptdown);
    }

    return ret;
}

std::vector<int> indexed_points::get_dim(){
    return index_layout;
}



/* manual_index_pts -- class for indexed points with manual index --------*/

manual_index_pts::manual_index_pts(){
    return;
}

manual_index_pts::~manual_index_pts(){
    return;
}

manual_index_pts::manual_index_pts(const manual_index_pts &mip) : indexed_points(mip) {
    return;
}

manual_index_pts & manual_index_pts::operator = (const manual_index_pts &mip){
    indexed_points::operator = (mip);
    return *this;
}

manual_index_pts *manual_index_pts::clone(){
    return new manual_index_pts(*this);
}

void manual_index_pts::set_layout(std::vector<int> index_layout){

    this->index_layout.erase(this->index_layout.begin(), this->index_layout.end());

    for(int i=0; i<index_layout.size(); i++){
        this->index_layout.push_back(index_layout[i]);
    }
    
    index.initArray(index_layout);

    index_set = true;

}

void manual_index_pts::set_index(int point_num, std::vector<int> loc){  
    index.setValue(loc, point_num);
}
