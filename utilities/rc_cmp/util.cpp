#include "util.h"

/*
 * int to string conversion
 */
string itoa(int n){
    string s;
    std::stringstream out;
    out << n;
    return out.str();
}

bool check_non_numeric_char(char test_char){

    int num = (int)test_char;

    if(num > 57 || num < 48)
        return true;

    return false;
}

bool check_non_numeric_char_per(char test_char){

    int num = (int)test_char;

    if(num > 57 || num < 48 || test_char == '.')
        return true;

    return false;
}

bool check_non_numeric_str(string test_str){
    for(int i=0; i<test_str.length(); i++){
        if(check_non_numeric_char(test_str[i]))
            return true;
    }

    return false;
}

bool check_non_numeric_str_per(string test_str){
    for(int i=0; i<test_str.length(); i++){
        if(check_non_numeric_char_per(test_str[i]))
            return true;
    }

    return false;
}

bool check_delimiters(char *delimit, int num, char test){

    //Default delimiters if none are specified
    if(delimit == NULL){
        num = 4;
        delimit = "\n\t \0";
    }
    
    //Check for delimiters
    for(int i=0; i < num; i++){
        if(test == delimit[i])
            return true;
    }

    return false;
}



//TODO
//Rewrite this ugly function using regex

/*
 * Read expressions of the form
 * x = "test" "100" "another test"
 */
std::vector<string> readQuoteExpression(string rStack, string rNeedle){

    std::vector<string> null_res;

    //
    //Search for needle in stack
    //

    bool needle_found;

    int i;
    for(i=0; i < (rStack.length()-rNeedle.length()) + 1; i++){

        char check = '\0';
        if( i > 0 )
            check = rStack[i-1];

        if( i == 0 || check_delimiters(NULL, 0, check) ){

            needle_found = true;

            for(int j=0; j<rNeedle.length(); j++){
                if( rStack[i+j] != rNeedle[j] )
                    needle_found = false;
            }
        }
        
        if(needle_found)
            break;

    }

    //Not found
    if(!needle_found)
        return null_res;
    
    i += rNeedle.length();



    //
    //Skip Whitespace, check for '='
    //

    char c = rStack[i];
    while(c == ' '){
        i++;
        c = rStack[i];
    }

    //Check for '='
    if(c != '=')
        return null_res;
    else
        i++;

    c = rStack[i];
    while(c == ' '){
        i++;
        c = rStack[i];
    }


    std::vector<string> values;



    //
    //Read Values
    //

    while(c == '"'){

        //Check for valid parentheses
        int end = rStack.length();

        i++;
        for(end=0; end<rStack.length(); end++){
            c = rStack[i+end];
            if(c == '"')
                break;  
        }

        //Read value
        string val("");

        int limit = end+i;

        for(i; i < limit; i++){
            val += rStack[i];
        }
        
        i++;

        values.push_back(val);

        //Skip whitespace
        c = rStack[i];
        while(c == ' '){
            i++;
            c = rStack[i];
        }

    }

    return values;
}



std::vector<string> readCommaExpression(string rStack, string rNeedle){

    std::vector<string> null_res;


    //
    //Search for needle in stack
    //

    bool needle_found;

    int i;
    for(i=0; i < (rStack.length()-rNeedle.length()) + 1; i++){

        char check = '\0';
        if( i > 0 )
            check = rStack[i-1];

        if( i == 0 || check_delimiters(NULL, 0, check) ){

            needle_found = true;

            for(int j=0; j<rNeedle.length(); j++){
                if( rStack[i+j] != rNeedle[j] )
                    needle_found = false;
            }
        }
        
        if(needle_found)
            break;

    }

    //Not found
    if(!needle_found)
        return null_res;
    
    i += rNeedle.length();



    //
    //Skip Whitespace, check for '='
    //

    char c = rStack[i];
    while(c == ' '){
        i++;
        c = rStack[i];
    }

    //Check for '='
    if(c != '=')
        return null_res;
    else
        i++;

    c = rStack[i];
    while(c == ' '){
        i++;
        c = rStack[i];
    }


    std::vector<string> values;


    //
    //Read Values
    //

    //Move back one character
    // - Neccessary due to loop
    i--;

    do{

        i++;
        c = rStack[i];

        //Skip whitespace
        while(c == ' '){
            i++;
            c = rStack[i];
        }

        //Read value
        string val("");

        while(c != ' ' && c != ',' && c != '='){
            val += rStack[i];
            
            i++;
            c = rStack[i];
        }

        if( c == '=')
            break;

        values.push_back(val);
    }
    while(c == ',');


    return values;
}

bool check_bool_arr(bool test[], int size){

    for(int i=0; i<size; i++){
        if(!test[i])
            return false;
    }

    return true;

}

long double dist(long double *coord1, long double *coord2, int size){
    
    long double total = 0;
    
    for(int i=0; i<size; i++){
        total += pow( coord2[i] - coord1[i], 2 );
    }

    total = sqrt(total);

    return total;
}

bool pnt_in_tri(pnt test, pnt p1, pnt p2, pnt p3){

    if(test.size < 2 || p1.size < 2 || p2.size < 2 || p3.size < 2)
        return false;

    if( side_check(test, p1, p2, p3) &&
        side_check(test, p2, p3, p1) &&
        side_check(test, p3, p1, p2)   ){
            return true;
    }

    return false;
}

bool side_check(pnt test, pnt p1, pnt p2, pnt check){

    //Line to check
    pnt line;
    line.size = 3;
    line.vals = new long double[3];

    line.vals[0] = p2.vals[0] - p1.vals[0];
    line.vals[1] = p2.vals[1] - p1.vals[1];
    line.vals[2] = 0.0;

    //Test point vector
    pnt tpt;
    tpt.size = 3;
    tpt.vals = new long double[3];

    tpt.vals[0] = test.vals[0] - p1.vals[0];
    tpt.vals[1] = test.vals[1] - p1.vals[1];
    tpt.vals[2] = 0.0;

    //Check point vector
    pnt checkvec;
    checkvec.size = 3;
    checkvec.vals = new long double[3];

    checkvec.vals[0] = check.vals[0] - p1.vals[0];
    checkvec.vals[1] = check.vals[1] - p1.vals[1];
    checkvec.vals[2] = 0.0;

    pnt res1 = cross_product(line, tpt);
    pnt res2 = cross_product(line, checkvec);

    if( dot_product(res1, res2) >= 0.0 )
        return true;

    return false;
}

pnt cross_product(pnt vec1, pnt vec2){

    pnt ret;

    if(vec1.size < 3 || vec2.size < 3)
        return ret;

    ret.size = 3;
    ret.vals = new long double[3];

    ret.vals[0] = (vec1.vals[1]*vec2.vals[2] - vec1.vals[2]*vec2.vals[1]);
    ret.vals[1] = (vec1.vals[2]*vec2.vals[0] - vec1.vals[0]*vec2.vals[2]);
    ret.vals[2] = (vec1.vals[0]*vec2.vals[1] - vec1.vals[1]*vec2.vals[0]);

    return ret;
}

long double dot_product(pnt vec1, pnt vec2){
    if(vec1.size != vec2.size){
        return nan();
    }

    long double ret_val = 0.0;

    for(int i=0; i<vec1.size; i++){
        ret_val += vec1.vals[i] * vec2.vals[i];
    }

    return ret_val;
}

long double nan(){
    
    //Return NaN without compiler warning
    long double *zer = new long double;
    *zer = 0.0;

    long double ret = ( (*zer) / (*zer) );
    delete zer;

    return ret;
}

std::vector<string> str_tokenize(string full, string delimiters){
    
    std::vector<string> ret;
    
    string::size_type start = full.find_first_not_of(delimiters, 0);
    string::size_type end = full.find_first_of(delimiters, start);

    while(start != string::npos || end != string::npos){
        ret.push_back(full.substr(start, end-start));

        start = full.find_first_not_of(delimiters, end);
        end = full.find_first_of(delimiters, start);
    }

    return ret;
}

bool same_elements(std::vector<int> vector_1, std::vector<int> vector_2){
    for(int vec1_element_num = 0; vec1_element_num < vector_1.size(); vec1_element_num++){
        
        bool match_found = false;
        for(int vec2_element_num = 0; vec2_element_num < vector_2.size(); vec2_element_num++){
            
            //check for matching element
            if(vector_1[vec1_element_num] == vector_2[vec2_element_num]){
                match_found = true;
                break;
            }
        }

        if(!match_found)
            return false;
    }

    return true;
}



/* bitmap class -- efficient bit array -----------------------------------*/

bitmap::bitmap(){
    map_set = false;
    return;
}

bitmap::~bitmap(){
    if(map_set)
        delete[] map;
}

bitmap::bitmap(const bitmap &b){
    size = b.size;
    map_set = b.map_set;

    int arr_size = size / (sizeof(int)*8);
    arr_size++;

    if(map_set){
        map = new int[arr_size];

        for(int i=0; i<arr_size; i++){
            map[i] = b.map[i];
        }
    }
}

bitmap & bitmap::operator =(const bitmap &b){

    if(this != &b){
        
        if(map_set)
            delete[] map;

        size = b.size;
        map_set = b.map_set;

        int arr_size = size / (sizeof(int)*8);
        arr_size++;

        if(map_set){
            map = new int[arr_size];

            for(int i=0; i<arr_size; i++){
                int temp = b.map[i];
                map[i] = temp;
            }
        }
    }

    return *this;
}

void bitmap::set_size(int n){

    int new_size = n / (sizeof(int)*8);
    new_size++;

    int old_size = size / (sizeof(int)*8);
    old_size++;

    if(!map_set){

        map = new int[ new_size ];

        for(int i=0; i < new_size; i++){
            map[i] = 0;
        }
    }
    else{
        int *temp;
        temp = new int[ new_size ];

        for(int i=0; i < new_size; i++){
            temp[i] = 0;
        }

        for(int i=0; i < old_size; i++){
            temp[i] = map[i];
        }   

        delete[] map;
        map = temp;
    }

    size = n;
    map_set = true;
}

int bitmap::get_size(){
    return size;
}

void bitmap::set_value(int n){

    int index = n / (sizeof(int) * 8);
    int num   = n % (sizeof(int) * 8);

    int entry = map[index];

    int set = 0x00000001;
    set = set << num;

    entry = entry | set;

    map[index] = entry;

}

void bitmap::unset_value(int n){
    
    int index = n / (sizeof(int) * 8);
    int num   = n % (sizeof(int) * 8);

    int entry = map[index];

    int set = 0x00000001;
    set = set << num;
    set = ~set;

    entry = entry & set;

    map[index] = entry;

}

bool bitmap::check_value(int n){

    if(n >= size || n < 0)
        return false;

    int index = n / (sizeof(int) * 8);
    int num   = n % (sizeof(int) * 8);

    int entry = map[index];

    entry = entry >> num;
    entry = entry & 0x00000001;

    if(entry == 0x00000001)
        return true;

    return false;
}

