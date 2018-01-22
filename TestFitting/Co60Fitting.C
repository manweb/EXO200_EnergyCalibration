double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = par[5];
  double sigma2 = sigma1*TMath::Sqrt(E2/E1);
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

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

void Co60Fitting()
{
  TChain *t = new TChain("t","tree");

  //t->Add("../../analysis/V2/2412_noReclustering_Fiducial.root");
  //t->Add("../../analysis/V2/2415_noReclustering_Fiducial.root");
  //t->Add("../../analysis/V2/2416_noReclustering_Fiducial.root");
  t->Add("../../analysis/MC/NewFitWindow/Co60_noReclustering_Fiducial.root");

  TH1F *h1 = new TH1F("h1","single site spectrum",200,0,2500);
  TH1F *h2 = new TH1F("h2","multi site spectrum",200,0,2500);
  TH1F *h3 = new TH1F("h3","multi site spectrum (2)",200,0,2500);
  TH1F *h4 = new TH1F("h4","multi site spectrum (3)",200,0,2500);
  TH1F *h5 = new TH1F("h5","multi site spectrum (4)",200,0,2500);

/*  t->Draw("ecrec>>h1","nsite == 1");
  t->Draw("ecrec>>h2","nsite > 1");
  t->Draw("ecrec>>h3","nsite == 2");
  t->Draw("ecrec>>h4","nsite == 3");
  t->Draw("ecrec>>h5","nsite == 4");
*/
  t->Draw("ecraw>>h1","nsite == 1");
  t->Draw("ecraw>>h2","nsite > 1");
  t->Draw("ecraw>>h3","nsite == 2");
  t->Draw("ecraw>>h4","nsite == 3");
  t->Draw("ecraw>>h5","nsite == 4");

  h2->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h4->SetLineColor(kRed);
  h5->SetLineColor(kRed);

  //h2->Scale(0.645);

  TF1 *fit = new TF1("fit",fitFunction,900,1500,6);
  fit->SetParNames("A1","E1","#sigma","R1","R2","E2");

  fit->SetParameters(650,1180,75,0.6,0.8,1370);

  fit->SetParLimits(1,1100,1250);
  fit->SetParLimits(2,10,200);
  fit->SetParLimits(3,0.001,1);
  fit->SetParLimits(4,0.001,1);
  fit->SetParLimits(5,1250,1400);

  h1->Fit("fit","r");

  double par[6];
  fit->GetParameters(par);

  TF1 *fitMulti = new TF1("fitMulti",fitFunctionMulti,800,1400,5);
  fitMulti->SetParNames("A1","E1","#sigma","R1","R2");

  fitMulti->SetParameters(5000,1173,70,0.2,0.8);

  fitMulti->SetParLimits(1,1000,1250);
  fitMulti->SetParLimits(2,10,200);
  fitMulti->SetParLimits(3,0.001,1.0);
  fitMulti->SetParLimits(4,0.001,1.0);

  h2->Fit("fitMulti","r");

  double parMulti[5];
  fitMulti->GetParameters(parMulti);

  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",800,1800);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *gaus2 = new TF1("gaus2","[0]*TMath::Gaus(x,[1],[2])",800,1800);
  gaus2->SetParameters(par[0]*par[4],par[5],par[2]*TMath::Sqrt(par[5]/par[1]));

  TF1 *erf1 = new TF1("erf1",errf,800,1800,3);
  erf1->SetParameters(par[0]*par[3],par[1],par[2]);

  TF1 *erf2 = new TF1("erf2",errf,800,1800,3);
  erf2->SetParameters(par[0]*par[4]*par[3],par[5],par[2]*TMath::Sqrt(par[5]/par[1]));

  TF1 *gaus1Multi = new TF1("gaus1Multi","gaus",600,1500);
  gaus1Multi->SetParameters(parMulti[0],parMulti[1],parMulti[2]);

  TF1 *gaus2Multi = new TF1("gaus2Multi","gaus",600,1500);
  gaus2Multi->SetParameters(parMulti[0]*parMulti[4],parMulti[1]*1.1358,parMulti[2]*TMath::Sqrt(1.1358));

  TF1 *erf1Multi = new TF1("erf1Multi",errf,800,1800,3);
  erf1Multi->SetParameters(parMulti[0]*parMulti[3],parMulti[1],parMulti[2]);

  TF1 *erf2Multi = new TF1("erf2Multi",errf,800,1800,3);
  erf2Multi->SetParameters(parMulti[0]*parMulti[4]*parMulti[3],parMulti[1]*1.1358,par[2]*TMath::Sqrt(1.1358));

  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  gaus2->SetLineWidth(1);
  gaus2->SetLineColor(kBlue);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);
  erf2->SetLineWidth(1);
  erf2->SetLineColor(kBlue);
  erf2->SetLineStyle(2);

  gaus1Multi->SetLineWidth(1);
  gaus1Multi->SetLineColor(kGreen+1);
  gaus2Multi->SetLineWidth(1);
  gaus2Multi->SetLineColor(kYellow+1);
  erf1Multi->SetLineWidth(1);
  erf1Multi->SetLineColor(kGreen+1);
  erf1Multi->SetLineStyle(2);
  erf2Multi->SetLineWidth(1);
  erf2Multi->SetLineColor(kYellow+1);
  erf2Multi->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  h1->Draw("EZ");
  h2->Draw("same");
  fit->Draw("same");
  fitMulti->Draw("same");
  gaus1->Draw("same");
  gaus2->Draw("same");
  gaus1Multi->Draw("same");
  gaus2Multi->Draw("same");
  erf1->Draw("same");
  erf2->Draw("same");
  erf1Multi->Draw("same");
  erf2Multi->Draw("same");

  return;
}
