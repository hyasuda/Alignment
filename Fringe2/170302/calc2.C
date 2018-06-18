#include <iostream>
#include <fstream>
#include <cmath>
#include <string>


void calc(std::string file = "0um.dat"){

    std::ifstream is(file,std::ios::in);

    if(is){
        is.seekg(0, is.end);
        int length = is.tellg();
        is.seekg(0, is.beg);

        // char * buffer = new char [length];


        std::cout << "Reading " << length << " characters... ";

        //read data as a block;
        is.read (buffer,length);

        if(is)
            std::cout << "all characters read successfully!" << endl;
        else
            std::cout << "error : only " << is.gcount() << " could be read";
        is.close();

x        delete[] buffer;
     }

    else
        std::cout << "File could not be found ... sorry ... " << endl;

    // return 0;
}
