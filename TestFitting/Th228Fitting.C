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
  TChain *t = new TChain("t","tree");

  //t->Add("../../analysis/V2/2417_noReclustering_Fiducial.root");
  //t->Add("../../analysis/V2/2418_noReclustering_Fiducial.root");
  //t->Add("../../analysis/V2/2421_noReclustering_Fiducial.root");
  t->Add("../../analysis/MC/NewFitWindow/Th228_noReclustering_Fiducial.root");

  TH1F *h1 = new TH1F("h1","single site spectrum",200,0,3500);
  TH1F *h2 = new TH1F("h2","multi site spectrum",200,0,3500);
  TH1F *h3 = new TH1F("h3","multi site spectrum (2)",200,0,3500);
  TH1F *h4 = new TH1F("h4","multi site spectrum (3)",200,0,3500);
  TH1F *h5 = new TH1F("h5","multi site spectrum (4)",200,0,3500);

  t->Draw("ecraw>>h1","nsite == 1");
  t->Draw("ecraw>>h2","nsite > 1");
  t->Draw("ecraw>>h3","nsite == 2");
  t->Draw("ecraw>>h4","nsite == 3");
  t->Draw("ecraw>>h5","nsite == 4");

  h2->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h4->SetLineColor(kRed);
  h5->SetLineColor(kRed);

  h2->Scale(0.289);

  TF1 *fit = new TF1("fit",fitFunction,2200,3400,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameters(800,2700,100,0.8);

  fit->SetParLimits(1,2300,2900);
  fit->SetParLimits(2,10,200);
  fit->SetParLimits(3,0.0,2.0);

  h1->Fit("fit","r");

  double par[4];
  fit->GetParameters(par);

  TF1 *fitMulti = new TF1("fitMulti",fitFunctionMulti,2500,2640,4);
  fitMulti->SetParNames("A","E","#sigma","R");

  fitMulti->SetParameters(150,2614,100,0.2);

  fitMulti->SetParLimits(1,2400,2900);
  fitMulti->SetParLimits(2,50,200);
  fitMulti->SetParLimits(3,0.0,0.8);

  h2->Fit("fitMulti","r");

  double parMulti[4];
  fitMulti->GetParameters(parMulti);

  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",2000,3000);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *erf1 = new TF1("erf1",errf,2000,3000,3);
  erf1->SetParameters(par[0]*par[3],par[1],par[2]);

  TF1 *gaus1Multi = new TF1("gaus1Multi","gaus",2000,3000);
  gaus1Multi->SetParameters(parMulti[0],parMulti[1],parMulti[2]);

  TF1 *erf1Multi = new TF1("erf1Multi",errf,2000,3000,3);
  erf1Multi->SetParameters(parMulti[0]*parMulti[3],parMulti[1],parMulti[2]);

  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);

  gaus1Multi->SetLineWidth(1);
  gaus1Multi->SetLineColor(kGreen+1);
  erf1Multi->SetLineWidth(1);
  erf1Multi->SetLineColor(kGreen+1);
  erf1Multi->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  h1->Draw("EZ");
  h2->Draw("same");
  fit->Draw("same");
  fitMulti->Draw("same");
  gaus1->Draw("same");
  gaus1Multi->Draw("same");
  erf1->Draw("same");
  erf1Multi->Draw("same");

  return;
}
