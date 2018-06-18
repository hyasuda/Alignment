#include "fstream"
#include "TH1D.h"
#include "TVirtualfft.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

void Alignment(int fCut = 100){

    TCanvas *c1 = new TCanvas("c1","c1");
    TCanvas *c2 = new TCanvas("c2","c2");
    TCanvas *c3 = new TCanvas("c3","c3");
    TCanvas *c4 = new TCanvas("c4","c4");
    TCanvas *c5 = new TCanvas("c5","c5");
    TCanvas *c6 = new TCanvas("c6","c6");
    TCanvas *c7 = new TCanvas("c7","c7");
    TCanvas *c8 = new TCanvas("c8","c8");
    TCanvas *c9 = new TCanvas("c9","c9");

    int n = 800;
    // int n;
    double t,v,v2,s,m;
    double RangeMin = -0.04000;
    double RangeMax =  0.04000;
    double Fs = 10.000;
    // double TimeInt = 0.0001;
    // double NoiseCut = 0.05;
    double re,im;
    double *re_full = new Double_t[n];
    double *im_full = new Double_t[n];
    double *re_full_cut = new Double_t[n];
    double *im_full_cut = new Double_t[n];
    double diff = 0.;
    double diffsw = 0.;
    // double fCut = 20;
    double HistoMax = 4.;
    double ave;
    double sum;
    // n = (RangeMax - RangeMin) / TimeInt;

    TH1D *hist = new TH1D("RawData","RawData",n,RangeMin,RangeMax);
    TH1D *hist2 = new TH1D("RawSquared","RawSquared",n,RangeMin,RangeMax);
    TGraph *h1 = new TGraph("/Users/YASUDA/Data/Alignment/170124/print_0011.dat");
    // TH1 *hm = 0;
    // TH1 *hmcut = 0;

    TH1 *hm = new TH1D("hm","hm",n,0.,Fs);
    TH1 *hmcut = new TH1D("hmcut","hmcut",n,0.,Fs);

    // TH1 *hm2 = 0;
    // TH1 *hmcut2 = 0;

    TH1 *hbcut = 0;
    // TH1 *hbcut;
    TH1 *hbcut2 = new TH1D("hbcut2","hbcut2",n,RangeMin,RangeMax);

    h_diff   = new TH1D(  "Diff",  "Diff",n,0,800);
    h_diffsw = new TH1D("Diffsw","Diffsw",n,0,800);

    h_diff2   = new TH1D(  "Diff2",  "Diff",n,RangeMin,RangeMax);
    h_diffsw2 = new TH1D("Diffsw2","Diffsw",n,RangeMin,RangeMax);

    // TH1D *hnoise = new TH1D("hnoise","hnoise",1600,0.0,10.0);
    TH1D *hnoise = new TH1D("hnoise","hnoise",400,0.0,10.0);

    ifstream data("/Users/YASUDA/Data/Alignment/170124/print_0011.dat");

    std::vector<double> vx, vy1, vy2;

    while(data >> t >> v >> s >> m){

        if(RangeMin < t && t < RangeMax){
            vx.push_back(t);
            vy1.push_back(m);
            vy2.push_back(s);
        }

        int binn = hist->FindBin(t);
        if(fabs(hist->GetBinLowEdge(binn)-t)>0.00002 && binn!=0 && binn!=n+1){
            binn += 1;
        }
        v2 = v*v; //pow(v,2);
        // if(v2<NoiseCut){
        //     v2 = 0;
        // }
        hist->SetBinContent(binn,v);
        hist2->SetBinContent(binn,v2);

        hnoise->Fill(m);
    }

    // for(int i=0;i<n;i++){
    //     cout << i << " " <<  << endl;
    // }

    TGraph *hms = new TGraph(vx.size(), &vx[0], &vy1[0]);
    TGraph  *hs = new TGraph(vx.size(), &vx[0], &vy2[0]);

    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.1);

    hist->SetXTitle("Second [s]");
    hist->SetYTitle("Voltage [V]");

    hist2->SetXTitle("Second[s]");
    hist2->SetLineColor(2);

    hm = hist2->FFT(hm,"MAG");
    hm->SetName("MAGnoCut");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->SetLineWidth(2);
    // hm->SetXTitle("# of Bin");
    hm->SetXTitle("Frequency [kHz]");
    hm->SetYTitle("Magnitude");

    // for(int i = 0; i < n ; i++){

    // }

    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();

    fft->GetPointComplex(0,re,im);
    printf("1st transform: DC component: %f\n", re);
    fft->GetPointComplex(n/2+1, re, im);
    printf("1st transform: Nyquist harmonic: %f\n", re);

    fft->GetPointsComplex(re_full, im_full);

    TVirtualFFT *fftcut = TVirtualFFT::GetCurrentTransform();

    for(int i=0;i < n ;i++){
        re_full_cut[i] = 0;
        im_full_cut[i] = 0;
    }

    for(int i=0;i < fCut ;i++){
        fftcut->GetPointComplex(  i,re_full_cut[  i], im_full_cut[  i]);
        fftcut->GetPointComplex(n-i,re_full_cut[n-i], im_full_cut[n-i]);
    }

    hmcut = hist2->FFT(hmcut,"MAG");

    for(int i=0;i<n;i++){
        hmcut->SetBinContent(i+1,sqrt(pow(re_full_cut[i],2)+pow(im_full_cut[i],2)));
    }


    TVirtualFFT::SetTransform(0);
    hmcut->SetName("MAGCut");
    hmcut->SetLineColor(2);


    TVirtualFFT *fftBackCut = TVirtualFFT::FFT(1,&n,"C2R M K");
    fftBackCut->SetPointsComplex(re_full_cut,im_full_cut);
    fftBackCut->Transform();

    hbcut = TH1::TransformHisto(fftBackCut,hbcut,"Re");
    hbcut->SetName("LowPathFilter");
    hbcut->SetTitle("Frequency Cut");
    hbcut->SetXTitle("# of Bin");
    // hbcut->SetYTitle("# of Bin");

    // for(int i=0; i < n ; i++){
    //     hbcut2->SetBinContent(i,hbcut->GetBinContent(i));
    // }

    for(int i=0; i < n ; i++){
        hbcut2->SetBinContent(i,hbcut->GetBinContent(i)/n);
    }



    for(int i=0;i<n;i++){
        diff = (hbcut->GetBinContent(i+1)) - ((hbcut->GetBinContent(i)));
        h_diff->SetBinContent(i,diff);

        if(diffsw){
            if(diff < 0.0) diffsw = 0;
        } else {
            if(diff > 30.0) diffsw = HistoMax;
        }

        h_diffsw->SetBinContent(i,diffsw);
        h_diffsw2->SetBinContent(i,diffsw);
    }

    h_diffsw->SetLineColor(2);
    h_diffsw2->SetLineColor(2);

    hnoise->SetTitle("Target Stage");
    hnoise->SetXTitle("Voltage [V]");
    hnoise->SetYTitle("Event");
    hnoise->Fit("gaus");

    hbcut2->SetTitle("Fringe after Low Path Filter");
    hbcut2->SetXTitle("Time[sec]");
    hbcut2->SetYTitle("Rescaled Voltage [V]");
    hs->SetTitle("ScanningStage;Time[t];Voltage[V]");
    // hs->Fit("pol1","","",RangeMin + 0.001,RangeMax - 0.001);
    hs->SetMarkerStyle(20);
    hs->SetMarkerSize(0.1);
    gStyle->SetOptFit();


    // hs->SetTitle("MeasurementStage;Time[t],Voltage[V]");


    c1->cd();
    hist->Draw();
    h1->Draw("PSAME");

    c2->cd();
    hist2->Draw();

    c3->cd();
    hm->Draw();
    hmcut->Draw("SAME");

    c4->cd();
    // hbcut->Draw();
    hbcut2->Draw();
    h_diffsw2->Draw("SAME");
    // h_diffsw->Draw();
    // hbcut->Draw("SAME");


    c5->cd();
    hnoise->Draw();

    c6->cd();
    hs->Draw("AP");
    // h_diffsw2->Draw("SAME");

    c7->cd();
    h_diff->Draw();
    h_diffsw->Draw("SAME");

    c8->cd();
    hbcut->Draw();
    h_diffsw->Draw("SAME");

    c9->cd();
    hbcut2->Draw();
    h_diffsw2->Draw("SAME");
    hs->Draw("P");
}
