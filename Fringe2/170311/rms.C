#include <fstream>
#include <sstream>
#include <string>



void calc(string file = "0um.dat"){

    std::ifstream infile(file);

    std::string line;
    int lineno = 0;
    double sum = 0;
    double rms = 0;
    double mean = 0;
    double res = 0;
    double nres = 0;
    double sigma = 0;
    std::vector<double> l;

    while (std::getline(infile, line))
    {

        std::istringstream iss(line);
        double length;


        if(!(iss >> length))
        {
            std::cout << "file is empty !!" << endl;
        }

        sum += length;
        l.push_back(length);
        cout << " length " << length << endl;
        lineno++ ;
    }

    for(int i = 0; i < l.size() ; i++){
        cout << "length = " << l[i] << endl;
    }

    mean = sum/lineno;
    cout << " line number =  " << lineno << endl;
    cout << " Mean =  " << mean << endl;

    for(int i = 0;i < l.size() ; i++){
        res = l[i] - mean;
        nres += res*res;
        cout << "Residual = " << res << " || " << "Residua square sum = " << nres << endl;
    }

    rms = sqrt(nres/lineno);
    sigma = sqrt(lineno)/sqrt((lineno-1))*rms;
    cout << " RMS [Measure - Mean] = " << rms << endl;
    cout << " Standard deviation = " << sigma << endl;
    // cout << " Error = " << (double) sqrt(1/(lineno-1))*rms << " = " << sigma/sqrt(lineno) << endl;
    cout << " Mean =  " << mean << endl;
    cout << " Error ( SD of Mean ) = " << sigma/sqrt(lineno) << endl;
    // cout << lineno << " " << sqrt(1/(lineno-1)) << " " << 1/sqrt(lineno-1) << endl;
}
