#include "fstream"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

void hist( string file = "0um.dat", int nbin = 5, double Ref = 0.){

    std::vector<double> l;
    Double_t length;

    ifstream data(file, std::ios::in);
    // if(!data){
    //     cout<<" ... error, file not found" << ":" << file << endl;
    //     return;
    // }

    TCanvas *c1 = new TCanvas("c1","c1");
    TH1D *hist = new TH1D("hist","Distribution of L(OC) around average",nbin,Ref - 5., Ref + 5.);
    hist->SetXTitle("L(OC) around average [um]");
    hist->GetXaxis()->SetTitleSize(.05);
    hist->SetStats(0);

    while(data >> length){
        hist->Fill( length );
    }



    // for(int ibin = 0; ibin < nbin ; ibin++){
    //     hist->Fill(  l[i] )
    // }
    // gStyle->SetOptFit();
    // hist->Fit("gaus");
    cout << "hist Entries = " << hist->GetEntries() << endl;
    double sum = 0;
    for(int ibin = 0; ibin < hist->GetEntries() ; ibin++){
        sum += hist->GetBinCenter(ibin);
    }
    cout << "Histogram Average =  "  << sum/hist->GetEntries() << endl;
    hist->Draw();

}
