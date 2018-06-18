{
    c1 = new TCanvas("c1","(OC - ONOSOKKI) vs. ONOSOKKI",200,10,700,500);
    // c1->SetFillColor(42);
    c1->SetGrid();
    // c1->GetFrame()->SetFillColor(21);
    // c1->GetFrame()->SetBorderSize(12);

    double p0 = 0;

    /* Error -> Total RMS , included data of 40 um */
    // const Int_t n = 6;
    // Double_t x[n]  = {-0.29,  9.42, 19.49, 29.44, 39.35, 49.46};
    // Double_t y[n]  = {95.90, 96.41, 96.69, 96.19, 96.99, 96.40};
    // Double_t ex[n] = {  .03,   .02,   .03,   .05,   .0 ,   .01};
    // Double_t ey[n] = {  .11,   .11,   .15,   .17,   .34,   .20};

    /* Error -> Each RMS , Not included data of 40 um */
    const Int_t n = 5;
    Double_t x[n]  = {-0.29,  9.42, 19.49, 29.44, 49.46};
    Double_t y[n]  = {95.90, 96.41, 96.69, 96.19, 96.40};
    Double_t ex[n] = {  .03,   .02,   .03,   .05,   .01};
    Double_t ey[n] = {  .11,   .15,   .19,   .12,   .19};

    gr = new TGraphErrors(n,x,y,0,ey);
    gr->SetTitle("(OC - ONOSOKKI) vs. ONOSOKKI;ONOSOKKI [um];OC - ONOSOKKI [um]");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    // gStyle->SetOptStat(1);
    // gr->GetFunction("pol1")->SetLineWidth(2);
    gr->Fit("pol1");
    p0 = gr->GetFunction("pol1")->GetParameter(0);
    cout << gr->GetFunction("pol1")->GetParameter(0)<< " " << p0 << endl;
    // p0 = gr->GetFunction("pol0")->GetParameter(0);
    // cout << "pol0 = " << p0 << endl;
    // gr->GetFunction("pol0")->SetOptFit();
    // gr->SetOptFit();
    gStyle->SetOptFit();
    gr->Draw("AP");

    c2 = new TCanvas("c2","(OC - ONOSOKKI) vs. ONOSOKKI",200,10,700,500);
    c2->SetGrid();


    Double_t x_m[n]  = {-0.29-p0,  9.42-p0, 19.49-p0, 29.44-p0, 49.46-p0};
    Double_t y_m[n]  = {95.90-p0, 96.41-p0, 96.69-p0, 96.19-p0, 96.40-p0};

    grm = new TGraphErrors(n,x,y_m,0,ey);
    grm->SetTitle("(OC - ONOSOKKI) vs. ONOSOKKI;ONOSOKKI [um];OC - ONOSOKKI [um]");
    grm->SetMarkerColor(4);
    grm->SetMarkerStyle(21);
    grm->Draw("AP");

    return c1;
}
