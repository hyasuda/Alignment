#include "fstream"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

void NoiseDistribution(){

    TCanvas *c1 = new TCanvas("c1","c1");

    double t,n;
    double RangeMin = -0.04000;
    double RangeMax =  0.04000;

    TH1D *hnoise = new TH1D("hnoise","hnoise",51,-0.4,0.4);

    ifstream data("../../Data/Alignment/170126/noise.dat");
    std::vector<double> vx,vy;

    while(data >> t >> n){

        if(RangeMin < t && t < RangeMax){
            vx.push_back(t);
            vy.push_back(n);
        }

        hnoise->Fill(n);
    }

    hnoise->SetTitle("Scanning Stage");
    hnoise->SetXTitle("Voltage [V]");
    hnoise->SetYTitle("Event");
    hnoise->Fit("gaus");
    gStyle->SetOptFit();

    hnoise->Draw();
}
