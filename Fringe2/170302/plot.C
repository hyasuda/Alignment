{
    c1 = new TCanvas("c1","Displacement of ONOSOKKI");
    c1->SetGrid();

    const Int_t n = 4;
    Double_t x[n] = { 10   , 20   , 30   , 50    };
    Double_t y[n] = {  9.71, 19.78, 29.73, 49.74 };
    Double_t ex[n] = { 0   ,  0   ,  0   ,  0    };
    Double_t ey[n] = { 0.11,  0.12,  0.15,  0.10};

    gr = new TGraphErrors(n,x,y,ex,ey);
    gr->Fit("pol1");
    gStyle->SetOptFit();
    gr->SetTitle("Displacement of ONOSOKKI");
    gr->GetXaxis()->SetTitle("Stage Position [um]");
    gr->GetYaxis()->SetTitle("Displacement [um]");
    gr->SetMarkerColor(4);
    // gr->SetMarkerSize(1);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    // gROOT->SetStyle("plain");
    // gStyle->BuildLegend()->SetFillColor(10);
    // c1->BuildLegend()->SetFillColor(10);

    c2 = new TCanvas("c2","Displacement of Optical Comb");
    c2->SetGrid();
    const Int_t n_oc = 4;
    // Double_t x[n] = { 10   , 20   , 30   , 50    };
    Double_t y_oc[n] = { 10.00, 21.06, 29.31, 49.72 };
    Double_t ex_oc[n] = { 0   ,  0   ,  0   ,  0    };
    Double_t ey_oc[n] = { 2.27,  1.91,  2.18,  1.46 };

    gr_oc = new TGraphErrors(n_oc,x,y_oc,ex_oc,ey_oc);
    gr_oc->Fit("pol1");
    gStyle->SetOptFit();
    gr_oc->SetTitle("Displacement of Optical Comb Interference");
    gr_oc->GetXaxis()->SetTitle("Stage Position [um]");
    gr_oc->GetYaxis()->SetTitle("Displacement [um]");
    gr_oc->SetMarkerColor(4);
    // gr->SetMarkerSize(1);
    gr_oc->SetMarkerStyle(21);
    gr_oc->Draw("AP");

    c3 = new TCanvas("c3","Displacement");
    c3->SetGrid();
    gr->Draw();
    gr_oc->Draw("SAME");
    c3->BuildLegend();
}
