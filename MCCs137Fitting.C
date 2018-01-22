double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

double fitFunctionMulti(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

void MCCs137Fitting()
{
  char *FileName = "../analysis/MC/P2_SourceP4_px_137Cs_2.root";

  TFile *f = new TFile(FileName,"READ");
  TTree *t = (TTree*)f->Get("t");

  TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,1500);
  TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,1500);

  t->Draw("ecrec>>hMean","nsite == 1");
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

  TF1 *fit = new TF1("fit",fitFunction,550,700,4);
  fit->SetParNames("A1","E1","#sigma","A2");
  
  fit->SetParameters(500,630,15,0.2);
  
  fit->SetParLimits(1,580,700);
  fit->SetParLimits(2,5,200);
  fit->SetParLimits(3,0.01,1.0);
  
  hMean->Fit("fit","r");
  
  TF1 *fitMulti = new TF1("fitMulti",fitFunctionMulti,550,700,4);
  fitMulti->SetParNames("A","E","#sigma","R");
  
  fitMulti->SetParameters(150,662,15,0.2);
  
  fitMulti->SetParLimits(1,580,700);
  fitMulti->SetParLimits(2,5,200);
  fitMulti->SetParLimits(3,0.0,0.5);
  
  hMulti->Fit("fitMulti","r");
  
  TCanvas *c1 = new TCanvas();
  hMean->Draw();

  return;
}
