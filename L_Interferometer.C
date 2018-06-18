#include "fstream"
#include "TH1D.h"
#include "TVirtualfft.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

const int N = 1;
// const int M = 20;
const int M = 1;
const double n_air    = 1.000263926;
// const double f_rep    = 59.45189199998;
// const double f_rep    = 59.4525829999;
// const double f_rep = 1189.02;
const double f_rep = 1200.;


// const double DeltaL   = 94.4246;

const double e_nair   = 0.9e-6; // Fixed
const double e_frep   = 3.0e-12; // Fixed
// const double e_DeltaL = 0.39;

const double c_0 = 299792458; // light speed [m/s]

void L_Interferometer(const double DeltaL = 94.4246, const double e_DeltaL = 0.39){

    double L_0 = 0;
    double   L = 0;
    double e_L = 0;

    L_0 = N*c_0/(2*n_air*M*f_rep);
    L = L_0 + DeltaL;
    e_L = sqrt(pow(e_nair*L,2) + pow(e_frep*L,2) + pow(e_DeltaL,2) );

    cout << setprecision(10);
    cout << "***Setup***" << endl;
    cout << "light speed      :   c_0  = "       <<   c_0    << "[m/s]"   << endl;
    cout << "Integer          :     N  = "       <<     N    << ""        << endl;
    cout << "Integer by etalon:     M  = "       <<     M    << ""        << endl;
    cout << "Refractive index : n_air  = "       << n_air    << ""        << endl;
    cout << "Repetition freq. : f_rep  = "       << f_rep    << "[MHz]"   << endl;
    cout << "Delta L          : DeltaL = "       << DeltaL   << "[um]"   << endl;
    cout << "delta n_air      : delta_n_air = "  << e_nair   << ""       << endl;
    cout << "delta f_rep      : delta_f_rep = "  << e_frep   << "[MHz]"  << endl;
    cout << "delta Delta L    : delta_DeltaL = " << e_DeltaL << "[um]" << endl;
    cout << "***Result***" << endl;
    cout << "Pulse Interval   :      L_0  = "    <<   L_0    << "[um]"  << endl;
    cout << "L_interferomter  :      L    = "    <<     L    << "[um]"  << endl;
    cout << "Uncertainty      :delta_L    = "    <<   e_L    << "[um]"  << endl;
}
