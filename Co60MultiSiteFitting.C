double fitGauss(double *x, double *par)
{
  double A1 = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A2 = par[0]*par[3];
  double E2 = 1.1358*E1;
  double sigma2 = sigma1*TMath::Sqrt(E2/E1);

  return A1*TMath::Gaus(x[0],E1,sigma1) + A2*TMath::Gaus(x[0],E2,sigma2);
}

double fitFunction(double *x, double *par)
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

void Co60MultiSiteFitting()
{
/*  const int nrRuns = 25;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875};
*/
  const int nrRuns = 8;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690};
  double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690};

  TChain *t = new TChain("t");

  for (int i = 0; i < nrRuns; i++) {
     char fname[100];
     sprintf(fname,"../analysis/Tree/%i_noReclustering_Fiducial.root",runID[i]);

     t->Add(fname);
  }

  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;

  TH1F *h = new TH1F("h","multi site spectrum",200,0,3000);

  t->Draw("ecrec>>h","nsite > 1");

  TF1 *fit = new TF1("fit",fitGauss,1050,1400,4);
  fit->SetParNames("A1","E1","#sigma","R");
  fit->SetParameters(8000,1150,100,0.8);

  fit->SetParLimits(3,0.5,1);

  TF1 *fit2 = new TF1("fit2",fitFunction,900,1500,5);
  fit2->SetParNames("A1","E1","#sigma1","R1","R2");

  fit2->SetParameters(5000,1173,70,0.6,0.8);
  fit2->SetParLimits(1,1100,1250);
  fit2->SetParLimits(2,50,200);
  fit2->SetParLimits(3,0.1,0.5);
  fit2->SetParLimits(4,0.5,1.0);

  h->Fit("fit2","r");

  double par[5];
  double *parErr;
  fit2->GetParameters(par);
  parErr = fit2->GetParErrors();

  TF1 *gaus1 = new TF1("gaus1","gaus",800,1500);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *gaus2 = new TF1("gaus2","gaus",800,1500);
 gaus2->SetParameters(par[0]*par[3],1.1358*par[1],par[2]*TMath::Sqrt(1.1358));

  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);

  gaus2->SetLineWidth(1);
  gaus2->SetLineColor(kBlue);

  cout << "Peak1: " << par[1] << " +- " << parErr[1] << endl;
  cout << "Peak2: " << par[1]*1.1358 << " +- " << parErr[1]*1.1358 << endl;

  TCanvas *c1 = new TCanvas();
  h->Draw();
  gaus1->Draw("same");
  gaus2->Draw("same");

  return;
}
