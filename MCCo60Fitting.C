double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = par[5];
  double sigma2 = par[6]; //sigma1*TMath::Sqrt(E2/E1);
  double A2_erf = A2_gaus*par[3];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  double gauss2 = A2_gaus * TMath::Gaus(x[0],E2,sigma2);
  double erf2 = A2_erf * 0.5 * TMath::Erfc((x[0] - E2) / (TMath::Sqrt(2)*sigma2));
  
  return gauss1 + erf1 + gauss2 + erf2;
}

double fitFunctionMulti(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = 1.1358*E1;
  double sigma2 = sigma1*TMath::Sqrt(E2/E1);
  double A2_erf = A2_gaus*par[3];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  double gauss2 = A2_gaus * TMath::Gaus(x[0],E2,sigma2);
  double erf2 = A2_erf * TMath::Erfc((x[0] - E2) / (TMath::Sqrt(2)*sigma2));
  
  return gauss1 + erf1 + gauss2 + erf2;
}

void MCCo60Fitting()
{
  char *FileName = "../analysis/MC/P2_SourceP4_px_60Co_2.root";

  TFile *f = new TFile(FileName,"READ");
  TTree *t = (TTree*)f->Get("t");

  TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,2000);
  TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,2000);

  t->Draw("ecrec*1.026 - 1.988>>hMean","nsite == 1");
  t->Draw("ecrec>>hMulti","nsite > 1");
  
  /*int nsite;
  double ecrec;
  
  t->SetBranchAddress("ecrec",&ecrec);
  t->SetBranchAddress("nsite",&nsite);
  
  TRandom *rand = new TRandom();
  double r0 = 0.0;
  double r1 = 0.0;
  double r2 = 0.0;
  for (int i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    
    double ecSingle = ecrec*1.012 + 0.321;
    if (nsite == 1) {r0 = 0.01857; r1 = -1.044; r2 = 57.49;}
    else {r0 = 0.02953; r1 = -1.9; r2 = 81.85;}
    //double rr = sqrt((r0*ecrec + r1*sqrt(ecrec) + r2)**2 - 0.017*67*ecrec);
    //double rr = sqrt((r0*ecrec + r1*sqrt(ecrec) + r2)**2 - 15*15*ecrec/2615);
    double rr = TMath::Sqrt((r0*ecSingle + r1*sqrt(ecSingle) + r2)**2 - 14.96*14.96);
    
    ecSingle += rr*rand->Gaus();
    
    if (nsite == 1) {hMean->Fill(ecSingle);}
  }*/
  
  TF1 *fit = new TF1("fit",fitFunction,1000,1400,7);
  fit->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  
  fit->SetParameters(10000,1140,15,0.6,0.8,1310,15);
  
  fit->SetParLimits(1,1100,1180);
  fit->SetParLimits(2,10,200);
  fit->SetParLimits(3,0.01,1);
  fit->SetParLimits(4,0.01,1);
  fit->SetParLimits(5,1280,1380);
  fit->SetParLimits(6,10,200);
  
  hMean->Fit("fit","r");
  
  TF1 *fitMulti = new TF1("fitMulti",fitFunction,1000,1350,7);
  fitMulti->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  
  fitMulti->SetParameters(10000,1173,15,0.2,0.8,1333,15);
  
  fitMulti->SetParLimits(1,1100,1250);
  fitMulti->SetParLimits(2,10,200);
  fitMulti->SetParLimits(3,0.01,0.8);
  fitMulti->SetParLimits(4,0.5,2.0);
  fitMulti->SetParLimits(5,1280,1380);
  fitMulti->SetParLimits(6,10,200);
  
  hMulti->Fit("fitMulti","r");
  
  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  
  TCanvas *c2 = new TCanvas();
  hMulti->Draw();

  return;
}
