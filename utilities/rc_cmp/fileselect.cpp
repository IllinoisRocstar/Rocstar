#include "fileselect.h"

datafile * Fileselect::detectFiletype(ifstream &infile, ofstream &outfile, bool loud, string infile_name, int dim, 
                                      std::vector< std::vector<int> > rFieldMappings, std::vector<index_order> rIndexOrder,
                                      std::vector<adj_map> conv_factor, std::vector<adj_map> norm_val){

    datafile *ret = NULL;

    //Test if file is a tecplot file
    if( testTecplot(infile, string("BLOCK")) ){ 
        ret = new tecplot_data(infile, outfile, loud, infile_name, dim, rFieldMappings, rIndexOrder, conv_factor, norm_val);
    }
    //Filetype not supported, output error message
    else{
        outfile << "ERROR   |  File \"" << infile_name << "\":  ";
        outfile << "Unsupported file format" << "\n";
        
        cout << "ERROR   |  File \"" << infile_name << "\":  ";
        cout << "Unsupported file format" << "\n";

        ret = NULL;
    }

    return ret;

}

bool Fileselect::testTecplot(ifstream &infile, string datapacking){

    //Seek to beginning
    infile.clear();
    infile.seekg(0, std::ios::beg);

    char line[BUFFER_SIZE];
    string fdstr("");

    bool zone_found = false;
    bool zonetype_found = false;
    bool datapacking_found = false;

    bool ordered = false;

    bool ijk[3] = {false, false, false};
    bool nodes_elements[2] = {false, false};

    //Loop through file, check lines for required variable declarations
    while(true){

        infile.getline(line, BUFFER_SIZE-1);
        fdstr = line;



        //
        //Strip beginning whitespace
        //

        int offset = 0;

        while(fdstr.c_str()[offset] == ' ' && offset < fdstr.length())
            offset++;

        if(offset > 0){

            char tmp[fdstr.length()-offset];
            memset(tmp, '\0', fdstr.length()-offset);

            for(int i=0; i<fdstr.length(); i++)
                tmp[i] = fdstr.c_str()[i+offset];

            fdstr = "";
            fdstr = tmp;

        }


        
        //
        //Check header
        //

        if(fdstr.find("ZONE") != string::npos){
            zone_found = true;
        }

        //Try to read file's I value
        if(fdstr.find("I=") != string::npos){
            std::vector<string> i_val = readCommaExpression(fdstr, string("I"));

            if(i_val.size() == 1){
                if( !check_non_numeric_str(string(i_val[0])) ){
                    ijk[0] = true;      
                }
            }   
        }

        //Try to read file's J value
        if(fdstr.find("J=") != string::npos){   
            std::vector<string> j_val = readCommaExpression(fdstr, string("J"));

            if(j_val.size() == 1){
                if( !check_non_numeric_str(string(j_val[0])) ){
                    ijk[1] = true;
                }
            }   
        }

        //Try to read file's K value
        if(fdstr.find("K=") != string::npos){   
            std::vector<string> k_val = readCommaExpression(fdstr, string("K"));

            if(k_val.size() == 1){
                if( !check_non_numeric_str(string(k_val[0])) ){ 
                    ijk[2] = true;
                }
            }
        }

        //Look for a Nodes expression
        if(fdstr.find("Nodes=") != string::npos){
            std::vector<string> node_val = readCommaExpression(fdstr, string("Nodes"));
            
            if(node_val.size() == 1){
                if(atoi(node_val[0].c_str()) >= 0){
                    nodes_elements[0] = true;
                }
            }
        }

        //Look for an Elements expression
        if(fdstr.find("Elements=") != string::npos){
            std::vector<string> elem_val = readCommaExpression(fdstr, string("Elements"));
            
            if(elem_val.size() == 1){
                if(atoi(elem_val[0].c_str()) >= 0){
                    nodes_elements[1] = true;
                }
            }
        }


        //TODO
        //Integrate this with zone class code
        if(fdstr.find("ZONETYPE=Ordered") != string::npos){ 
            zonetype_found = true;
            ordered = true; 
        }
        else if(fdstr.find("ZONETYPE=FEQuadrilateral") != string::npos){
            zonetype_found = true;
            ordered = false;
        }
        //--------
        //

        if(fdstr.find("DATAPACKING=") != string::npos){

            char type[256];
            char space[256];
            sscanf(fdstr.c_str(), "DATAPACKING=%255s", type);

            if(string(type).compare(datapacking) == 0)
                datapacking_found = true;
        }



        //
        //Loop exit conditions
        //

        //Check detected data for valid file
        if(zone_found && zonetype_found && datapacking_found){
            if(ordered && check_bool_arr(ijk, 3))
                return true;
            else if(!ordered && check_bool_arr(nodes_elements, 2))
                return true;
        }

        //If at end of file, and not all required values have been found
        //return false
        if(infile.eof()){   
            return false;
        }

    }

    return false;
}

/*
 * Print list of valid filetypes to stdout
 */
void Fileselect::printFiletypes(){
    cout << "\t" << "Tecplot ASCII Block Datafile - Ordered Data" << endl;
    cout << "\t" << "Tecplot ASCII Block Datafile - FEQuadrilateral" << endl;
}
