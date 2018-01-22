double *PeakPos1;
double *PeakPosErr1;
double *PeakPosMulti;
double *PeakPosErrMulti;

void ProcessRun(int RunID, int ID);

TCanvas *c1 = new TCanvas();

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

void Th228Fitting()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
/*  const int nrRuns = 36;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
*/

// modified list
/*  const int nrRuns = 33;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1702, 1709, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1702, 1709, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
*/
/*  const int nrRuns = 15;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};
*/

  /*const int nrRuns = 13;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};*/
  //int runID[nrRuns] = {2417, 2421, 2422, 2424, 2426};
  //double runIDCopy[nrRuns] = {2417, 2421, 2422, 2424, 2426};
  
  const int nrRuns = 36;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPosMulti = new double[nrRuns];
  PeakPosErrMulti = new double[nrRuns];

  c1->Print("Th228_ss.ps[");
  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }
  c1->Print("Th228_ss.ps]");

  cout << "Mean (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (multi site) = " << TMath::Mean(nrRuns,PeakPosMulti) << " +- " << TMath::RMS(nrRuns,PeakPosMulti) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPosMulti,0,PeakPosErrMulti);

  gr1->SetTitle("Th228 peak (single site)");
  gr2->SetTitle("Th228 peak (multi site)");

  gr1->SetName("Th228_single");
  gr2->SetName("Th228_multi");

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetYaxis()->SetRangeUser(2650,2800);
  gr1->GetXaxis()->SetLimits(1570,3150);
  gr1->SetTitle("Th228 Peak");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetYaxis()->SetRangeUser(2650,2800);
  gr2->GetXaxis()->SetLimits(1570,3150);
  gr2->SetTitle("Th228 Peak (multi site)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");

  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");

  TFile *fOut = new TFile("../analysis/Th228_PeakPosition_070212.root","RECREATE");

  gr1->Write();
  gr2->Write();

  fOut->Close();

  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  sprintf(FileName,"../analysis/V9/%i_noReclustering_Fiducial.root",RunID);
  //sprintf(FileName,"../analysis/V3/%i_9us.root",RunID);

  cout << "Processing run " << RunID << "..." << endl;

  TFile *f = new TFile(FileName,"READ");
/*  TH1F *hMean = (TH1F*)f->Get("hMean");
  TH1F *hN = (TH1F*)f->Get("hN");
  TH1F *hP = (TH1F*)f->Get("hP");
  TH1F *hMulti = (TH1F*)f->Get("hMulti");
*/
  TTree *t = (TTree*)f->Get("t");

  double ecrec;
  int nsite;

  t->SetBranchAddress("ecrec",&ecrec);
  t->SetBranchAddress("nsite",&nsite);

  TH1F *hMean = new TH1F("hMean","Charge Spectrum (single sites)",100,2000,3500);
  TH1F *hMulti = new TH1F("hMulti","Charge Spectrum (multi site)m",100,2000,3500);

  //t->Draw("ecrec_gc>>hMean","nsite == 1");
  //t->Draw("ecrec>>hMean","nsite == 1 && TMath::Abs(zc) < 150");
  //t->Draw("ecrec_gc>>hMulti","nsite == 2");
  
  t->Draw("epcrec>>hMean","nsite == 1");
  t->Draw("epcrec>>hMulti","nsite > 1");

  TF1 *fit = new TF1("fit",fitFunction,2300,3000,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameters(800,2614,100,0.8);

  fit->SetParLimits(1,2500,2900);
  fit->SetParLimits(2,50,200);
  fit->SetParLimits(3,0.2,2.0);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunctionMulti,2300,2850,4);
  //TF1 *fitMulti = new TF1("fitMulti","gaus",2550,2850);
  fitMulti->SetParNames("A","E","#sigma","R");

  fitMulti->SetParameters(800,2614,100,0.2);

  fitMulti->SetParLimits(1,2400,2900);
  fitMulti->SetParLimits(2,50,200);
  fitMulti->SetParLimits(3,0.01,0.8);

  hMulti->Fit("fitMulti","r");

  double par[4];
  double *parErr;
  fit->GetParameters(par);
  parErr = fit->GetParErrors();

  PeakPos1[ID] = par[1];
  PeakPosErr1[ID] = parErr[1];

  double parMulti[4];
  double *parErrMulti;
  fitMulti->GetParameters(parMulti);
  parErrMulti = fitMulti->GetParErrors();

  PeakPosMulti[ID] = parMulti[1];
  PeakPosErrMulti[ID] = parErrMulti[1];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << ")" << "  (" << PeakPosMulti[ID] << " +- " << PeakPosErrMulti[ID] << ")" << endl;

  char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/%i.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  hMean->Write();
  hMulti->Write();
  fHisto->Close();
  
  c1->cd();
  hMean->Draw("EZP");
  c1->Print("Th228_ss.ps");

  hMean->Delete();
  //hN->Delete();
  //hP->Delete();
  hMulti->Delete();
  fit->Delete();
  fitMulti->Delete();

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",2000,3500);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *erf1 = new TF1("erf1",errf,2000,3500,3);
  erf1->SetParameters(par[3]*par[0],par[1],par[2]);
  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  gaus1->Draw("same");
  erf1->Draw("same");

  //TCanvas *c2 = new TCanvas();
  //hN->Draw();

  //TCanvas *c3 = new TCanvas();
  //hP->Draw();

  TCanvas *c4 = new TCanvas();
  hMulti->Draw()*/

  return;
}
