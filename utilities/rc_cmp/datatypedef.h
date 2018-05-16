/*
 * datatypedef.h
 *
 * Defines neccessary structs
 */


#ifndef _DATATYPEDEF_H_
#define _DATATYPEDEF_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using std::string;

/**
 * Used to store index order information
 */
struct index_order{
    std::vector<int> order;
    int partition_number;
};


/**
 * Used to store adjustment (normalization factor, conversion factor) information
 */
struct adj_map{
    int variable;
    long double factor;
};



/**
 * Used by cmd parser to store comparison information
 */
struct cmp_map{

    int metric;

    std::vector<int> file1_partitions;
    int var1;

    int file2_partition;
    int var2;
};



/**
 * Store color information
 */
struct rgb_color{
    int r;
    int g;
    int b;
};



/**
 * Point object that represents a single point
 */
struct pnt{
    long double *vals;
    int size;

    /**
     * Constructor that initializes values to defaults
     */
    pnt(){
        size = 0;
        vals = NULL;
    }

    /**
     * Destructor
     */
    ~pnt(){
        if(size > 0)
            delete[] vals;
    }

    /**
     * Deep copy constructor
     */
    pnt(const pnt &p){
        size = p.size;  
        
        if(size > 0){
            vals = new long double[size];

            for(int i=0; i<size; i++){
                vals[i] = p.vals[i];
            }
        }
    }

    /**
     * Equals operator, deep copy
     * @return Reference to new object
     */
    pnt & operator = (const pnt &p){
        
        if(this != &p){
            if(size > 0)
                delete[] vals;

            size = p.size;
            vals = new long double[size];

            for(int i=0; i<size; i++){
                vals[i] = p.vals[i];
            }
        }
        
        return *this;
    }
};



/**
 * Variable dimension array, an array whose dimensions can 
 * be specified at runtime. Values are stored in a flat array 
 * (an array with 1 dimension), and arithmetic is used to translate
 * between coordinates (with the array's dimensionality) to flat indices.
 */
struct VariableDimensionArray {
    
    /**
     * Default constructor that initializes values to default
     */
    VariableDimensionArray(){
        this->index = NULL;
        index_size = 0;
    }

    /**
     * Destructor
     */
    ~VariableDimensionArray(){
        delete[] index;
    }

    /**
     * Deep copy
     * @param vda Array to copy
     */
    VariableDimensionArray(const VariableDimensionArray &vda){
        this->sizes = vda.sizes;
        this->index_size = vda.index_size;

        if(index_size > 0){
            this->index = new int[index_size];

            for(int i=0; i<index_size; i++)
                index[i] = vda.index[i];
        }
    }

    /**
     * Equals operator deep copy
     * @param vda Array to copy
     * @return Reference to new object
     */
    VariableDimensionArray & operator = (const VariableDimensionArray &vda){
        if(this != &vda){
            if(index_size > 0)
                delete[] index;
        
            this->sizes = vda.sizes;
            this->index_size = vda.index_size;

            if(index_size > 0){
                this->index = new int[index_size];

                for(int i=0; i<index_size; i++)
                    index[i] = vda.index[i];
            }
        }

        return *this;
    }

    /**
     * Initialize the array with the specified dimensions
     * @param dim_sizes The size of each dimension of the array
     */
    void initArray(std::vector<int> dim_sizes){

        if(this->index != NULL)
            return;

        int i_size = 1;
        for(int i=0; i<dim_sizes.size(); i++){
            i_size *= dim_sizes[i];
        }

        index = new int[i_size];
        index_size = i_size;

        for(int i=0; i<i_size; i++)
            index[i] = -1;

        for(int i=0; i<dim_sizes.size(); i++){
            this->sizes.push_back(dim_sizes[i]);
        }
    }

    /**
     * Translate coordinates to the flat index actually used for storage
     * @param vals Coordinates to translate
     * @return Flat index location corresponding to the coordinates
     */
    int translate_coord_to_index(std::vector<int> vals){    
    
        int flat_i = 0;

        for(int i=0; i<vals.size(); i++){

            if(vals[i] >= sizes[i] || vals[i] < 0)
                return -1;

            int multiplier = 1;
            for(int m=(i+1); m<vals.size(); m++){
                multiplier *= sizes[m];
            }

            flat_i += vals[i] * multiplier;
        }

        return flat_i;
    }

    /**
     * Translate index to coordinates.
     * @param index Flat index value
     * @return Coordinates of the point specified by the flat index
     */
    std::vector<int> translate_index_to_coords(int index){
        if(index > index_size){
            std::vector<int> blank_ret;
            return blank_ret;
        }

        std::vector<int> ret;

        for(int i=0; i<sizes.size(); i++){

            int dimsize = 1;
            for(int m=(i+1); m<sizes.size(); m++){
                dimsize *= sizes[m];
            }       
            
            ret.push_back(index / dimsize);
            index = index % dimsize;
        }

        return ret;
    }

    /**
     * Get the value specified by the loc coordinates
     * @param rLoc Coordinates for the value to retrieve
     * @return The value
     */
    int getValue(std::vector<int> rLoc){

        if(rLoc.size() != sizes.size())
            return -1;

        int flat_i = translate_coord_to_index(rLoc);
        return index[flat_i];
    }

    /**
     * Set the value at the given coordinates
     * @param rLoc Location of the value to set
     * @param rVal Value to set
     */
    void setValue(std::vector<int> rLoc, int rVal){
        
        if(rLoc.size() != sizes.size())
            return;

        int flat_i = translate_coord_to_index(rLoc);
        index[flat_i] = rVal;
    }

    /** Flat storage of values. */
    int *index;
    int index_size;
    
    //Dimenional layout
    std::vector<int> sizes;
};



#endif
