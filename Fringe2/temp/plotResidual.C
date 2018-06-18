void plotResidual() {

   TH1D * h1 = new TH1D("h1","h1",50,-3, 3);

   // fill histogram with gaussian data
   h1->FillRandom("gaus");
   // h1->FillRandom("pol1");

   // fit the histogram
   h1->Fit("gaus");
   // h1->Fit("pol1");

   // make residual plot

   TH1D * h2 = new TH1D("h2","Fit Residuals",50,-3, 3);
   TF1 * fittedFunc = h1->GetFunction("gaus");
   // TF1 * fittedFunc = h1->GetFunction("pol1");

   for (int ibin = 1; ibin <= 50; ++ibin) {
      double res =  (h1->GetBinContent(ibin)- fittedFunc->Eval( h1->GetBinCenter(ibin) ) )/h1->GetBinError(ibin);
      h2->SetBinContent(ibin, res  );
      h2->SetBinError(ibin, 1  );
   }

   TCanvas * c1 = new TCanvas();
   c1->Divide(1,2);
   c1->cd(1);
   h1->Draw();
   c1->cd(2);
   h2->Draw("E");
}
