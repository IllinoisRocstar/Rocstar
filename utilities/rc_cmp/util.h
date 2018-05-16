#ifndef _UTIL_H_
#define _UTIL_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <time.h>
#include <sstream>
#include <math.h>
#include <vector>

#include "datatypedef.h"

using std::string;
using std::ofstream;
using std::ifstream;



/**
 * Convert an int to a string
 * @param n Integer to convert
 * @return String representation of the integer
 */
string itoa(int n);

/**
 * Read expressions of the form
 * x = "test" "100" "another test"
 * @param stack
 */
std::vector<string> readQuoteExpression(string rStack, string rNeedle);

//Read expressions of the form
//y = 100, test, this
std::vector<string> readCommaExpression(string rStack, string rNeedle);

//Returns true if char is non-numeric
bool check_non_numeric_char(char test_char);

//Returns true if char is non-numeric (except period)
bool check_non_numeric_char_per(char test_char);

//Returns true if any char in the string is non-numeric
bool check_non_numeric_str(string test_str);

//Returns true if any char in the string is non-numeric, except period char
bool check_non_numeric_str_per(string test_str);

bool check_delimiters(char *delimit, int num, char test);

//Checks to see if all values are true
bool check_bool_arr(bool test[], int size);

//Distance formula
long double dist(long double *coord1, long double *coord2, int size);

//Triangle Test
//Check if point is inside a triangle
bool pnt_in_tri(pnt test, pnt p1, pnt p2, pnt p3);

//Side check
//Helper function for triangle test
bool side_check(pnt test, pnt p1, pnt p2, pnt check);

//Cross product
pnt cross_product(pnt vec1, pnt vec2);

//Dot product
long double dot_product(pnt vec1, pnt vec2);

//Return's nan without compiler warning
long double nan();

//String tokenize
std::vector<string> str_tokenize(string full, string delimiters);

//Compare integer vectors
//Ignore element order
bool same_elements(std::vector<int> vector_1, std::vector<int> vector_2);



/*
 * Bitmap class
 * Used to keep track of set values
 *
 */

class bitmap {
    
    public:

        bitmap();

        ~bitmap();

        bitmap(const bitmap &b);

        bitmap & operator = (const bitmap &b);

        void set_size(int n);

        int get_size();

        void set_value(int n);

        void unset_value(int n);

        bool check_value(int n);

    private:

        //Map Data
        int *map;

        bool map_set;

        //Size
        int size;

};

#endif
