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

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

void Th228MultiSiteFitting()
{
/*  const int nrRuns = 36;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
*/
  const int nrRuns = 15;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};

  TChain *t = new TChain("t");

  for (int i = 0; i < nrRuns; i++) {
     char fname[100];
     sprintf(fname,"../analysis/Tree/%i_noReclustering_Fiducial.root",runID[i]);

     t->Add(fname);
  }

  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;

  TH1F *h = new TH1F("h","multi site spectrum",200,0,4500);

  t->Draw("ecrec>>h","nsite > 1");

  TF1 *fit = new TF1("fit","gaus",2550,2850);
  //fit->SetParNames("A1","E1","#sigma","R");
  //fit->SetParameters(8000,1150,100,0.8);

  //fit->SetParLimits(3,0.5,1);

  TF1 *fit2 = new TF1("fit2",fitFunction,2000,3000,4);
  fit2->SetParNames("A","E","#sigma","R");

  fit2->SetParameters(150,2600,100,0.2);
  fit2->SetParLimits(1,2400,2900);
  fit2->SetParLimits(2,50,200);
  fit2->SetParLimits(3,0.1,0.5);

  h->Fit("fit2","r");

  double par[4];
  double *parErr;
  fit2->GetParameters(par);
  parErr = fit2->GetParErrors();

  cout << "Peak: " << par[1] << " +- " << parErr[1] << endl;

  TCanvas *c1 = new TCanvas();
  h->Draw();

  return;
}
