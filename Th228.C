double *PeakPos1;
double *PeakPosErr1;
double *PeakPosMulti;
double *PeakPosErrMulti;

void ProcessRun(int RunID, int ID);

double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

void Th228()
{
  const int nrRuns = 41;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940};

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPosMulti = new double[nrRuns];
  PeakPosErrMulti = new double[nrRuns];

  for (int i = 16; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

  cout << "Mean (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (multi site) = " << TMath::Mean(nrRuns,PeakPosMulti) << " +- " << TMath::RMS(nrRuns,PeakPosMulti) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPosMulti,0,PeakPosErrMulti);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetYaxis()->SetRangeUser(2650,2800);
  gr1->SetTitle("Th228 Peak");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetYaxis()->SetRangeUser(2650,2800);
  gr2->SetTitle("Th228 Peak (multi site)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");

  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");

  return;
}

void ProcessRun(int RunID, int ID)
{
  TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,3500);
  TH1F *hN = new TH1F("hN","Charge Energy Spectrum (purity - 2us)",200,0,3500);
  TH1F *hP = new TH1F("hP","Charge Energy Spectrum (purity + 2us)",200,0,3500);
  TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,3500);

  char FileName[100];
  sprintf(FileName,"/EXO200Data/Disk5/processed/%i/recon0000%i-*.root",RunID,RunID);

  cout << "Processing run " << RunID << "..." << endl;

  TChain *t = new TChain("tree");
  t->Add(FileName);

  EXOEventData *ED = 0;
  t->SetBranchAddress("EventBranch",&ED);

  double elife = 0.0;

  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
     t->GetEntry(i);
     if (i%1000 == 0) {cout << i << " events processed (elife = " << elife << ")" << endl;}

     double purFitP0;
     double purFitP1;
     double purFitP2;
     double purFitP3;
     double purFitP4;

     double purTime = double(ED->fEventHeader.fTriggerSeconds - 1304146800.0) / 3600.0 / 24.0;

     if (purTime < 58) {
           purFitP0 = -284.596;
           purFitP1 = 53.6978;
           purFitP2 = -1.88664;
           purFitP3 = 0.0269101;
           purFitP4 = -0.000133772;
        }
        if (purTime >= 58 && purTime < 81.6) {
           purFitP0 = 14068.5;
           purFitP1 = -908.011;
           purFitP2 = 21.8864;
           purFitP3 = -0.230994;
           purFitP4 = 0.00090631;
        }
        if (purTime >= 81.6) {
           purFitP0 = -9011.55;
           purFitP1 = 115.417;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }

     elife = purFitP4*purTime*purTime*purTime*purTime + purFitP3*purTime*purTime*purTime + purFitP2*purTime*purTime + purFitP1*purTime + purFitP0;

     int nsc = ED->GetNumScintillationClusters();
     for (int scID = 0; scID < nsc; scID++) {
        EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(scID);

        if (scint_cluster->fTime > 1928000) {continue;}
        int ncl = scint_cluster->GetNumChargeClusters();

        double fX = 0.0;
        double fY = 0.0;
        double fZ = 0.0;
        double epcl_sum = 0.0;
        for (int clID = 0; clID < ncl; clID++) {
           EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);

           fX = charge_cluster->fX;
           fY = charge_cluster->fY;
           fZ = charge_cluster->fZ;
           if (TMath::Sqrt(fX*fX + fY*fY) > 163) {continue;}
           if (fZ > 172) {continue;}
           if (fZ < -172) {continue;}

           double dtcl = charge_cluster->fDriftTime / 1000.0;
           epcl_sum += charge_cluster->fCorrectedEnergy * TMath::Exp(dtcl/elife);
        }

        if (epcl_sum > 0.0) {hMulti->Fill(epcl_sum);}

        if (ncl != 1) {continue;}
        EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(0);

        fX = charge_cluster->fX;
        fY = charge_cluster->fY;
        fZ = charge_cluster->fZ;
        if (TMath::Sqrt(fX*fX + fY*fY) > 163) {continue;}
        if (fZ > 172) {continue;}
        if (fZ < -172) {continue;}

        double dtcl = scint_cluster->GetChargeClusterAt(0)->fDriftTime / 1000.0;
        double epcl = scint_cluster->GetChargeClusterAt(0)->fCorrectedEnergy * TMath::Exp(dtcl/elife);
        double epclN = scint_cluster->GetChargeClusterAt(0)->fCorrectedEnergy * TMath::Exp(dtcl/(elife - 2.0));
        double epclP = scint_cluster->GetChargeClusterAt(0)->fCorrectedEnergy * TMath::Exp(dtcl/(elife + 2.0));

        hMean->Fill(epcl);
        hN->Fill(epclN);
        hP->Fill(epclP);
     }
  }

  TF1 *fit = new TF1("fit",fitFunction,2000,3400,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameter(1,2700);

  fit->SetParLimits(1,2500,2900);
  fit->SetParLimits(2,50,200);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti","gaus",2500,2870);

  hMulti->Fit("fitMulti","r");

  double par[4];
  double *parErr;
  fit->GetParameters(par);
  parErr = fit->GetParErrors();

  PeakPos1[ID] = par[1];
  PeakPosErr1[ID] = parErr[1];

  double parMulti[3];
  double *parErrMulti;
  fitMulti->GetParameters(parMulti);
  parErrMulti = fitMulti->GetParErrors();

  PeakPosMulti[ID] = parMulti[1];
  PeakPosErrMulti[ID] = parErrMulti[1];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << ")" << "  (" << PeakPosMulti[ID] << " +- " << PeakPosErrMulti[ID] << ")" << endl;

  char oName[100];
  sprintf(oName,"../analysis/%i.root",RunID);

  TFile *oFile = new TFile(oName,"RECREATE");
  hMean->Write();
  hN->Write();
  hP->Write();
  hMulti->Write();

  oFile->Close();

  hMean->Delete();
  hN->Delete();
  hP->Delete();
  hMulti->Delete();
  fit->Delete();

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",2000,3500);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *erf1 = new TF1("erf1",errf,2000,3500,3);
  erf1->SetParameters(par[3],par[1],par[2]);
  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  gaus1->Draw("same");
  erf1->Draw("same");

  TCanvas *c2 = new TCanvas();
  hN->Draw();

  TCanvas *c3 = new TCanvas();
  hP->Draw();*/

  return;
}
