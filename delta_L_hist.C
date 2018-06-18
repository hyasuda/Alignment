// This code is to calculate the uncertainty of L ( Measurement of length in optical comb laser interferometer )
// Valuable is : e_DeltaL , e_feta ( relative )
// e_feta is determined by the calibration measurement with Gauge block.
// Delta_L is basically about 0.3 um ( from experience ) in Agilent oscilloscope.
// There is more one uncertainty for calculation of L uncertainty : nair.
// This uncertainty is limited by the precision of device ( thermometer etc... ).
// So I defined nair uncertainty is as 3.0e-12.
// Pulse interval is defined as the 125 mm with 1.2 GHz etalon.
// by this line, commented by H. Yasuda on 20 Mar. 2018.

// Add the ideal condition uncertainty
// Room temperature uncertainty suppressed by +-1.

#include "fstream"
#include "TH1D.h"
#include "TVirtualfft.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

// const int N = 1;
const int M = 1; //frep equals design value of fiber etalon
// const int M = 20; // Integral multiple of frep

//--- Constant : Physics and Setup ---//
const double c0 = 299792458e+3; // light speed [mm/s]
// const double c0 = 299792458e+6; // light speed [um/s]
// const double feta = 1.200e+9; // etalon repetition frequency [Hz]
// const double feta = 1.199765e+9; // calibrated etalon repetition frequency [Hz]
const double nair = 1.0002647;
const double Le = 700; // end length of the plot [mm]

const double dnair = 0.9e-6; // Uncertainty of Refractive index : relative
// const double efrep = 3.0e-12; // Uncertainty of frep : relative ( not used )
// const double e_DeltaL = 0.34; // Uncertainty of DeltaL


const double dig = 0.1; // The resolution of the the plot

using namespace std;

void delta_L(const double eDeltaL = 0.34, const double dfeta = 1.6e-5, const double feta = 1.199765e+9){

    TCanvas * c1 = new TCanvas("c1","c1");
    TCanvas * c2 = new TCanvas("c2","c2");

    TGraph * g1 = new TGraph();
    TH1D * h1 = new TH1D("h1", "h1", 2, 0, 700);

    double L = 0;
    double eL = 0;

    int N = 0;

    //--- for test ---//
    N = 1;
    double l0 = c0/(2*nair*feta); // [um]
    cout << l0 * N << endl;
    N = 0;
    //----------------//

    for(int i = 0 ; i < (int)Le/dig ; i++ ){
        L = i*dig; // [mm]
        if( L > l0/2 * (1 + 2*N) ) N++;
        eL = sqrt( pow(eDeltaL,2) + pow(l0*N*1000,2) * ( pow(dnair,2) + pow(dfeta,2) ) ); // [um]
        cout << N << " " << i*dig << " " << eL << endl;
        g1->SetPoint(i,i*dig,eL);
    }

    c1->cd();
    g1->Draw("AL");
    g1->GetYaxis()->SetRangeUser(0,12.5);
    g1->GetXaxis()->SetRangeUser(0,700);
    g1->SetLineWidth(3);


    L = 0;
    eL = 0;
    N = 0;

    for(int N = 0 ; L < (int)Le/dig ; N++){
        L = l0 * N;
        eL = sqrt( pow(eDeltaL,2) + pow(l0*N*1000,2) * ( pow(dnair,2) + pow(dfeta,2) ) ); // [um]
        h1->SetBinContent(N, L, eL);
    }

    c2->cd();
    h1->Draw();
}
