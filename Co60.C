double *PeakPos1;
double *PeakPos2;
double *PeakPosErr1;
double *PeakPosErr2;
double *PeakPos1Multi;
double *PeakPos2Multi;
double *PeakPosErr1Multi;
double *PeakPosErr2Multi;

void ProcessRun(int RunID, int ID);

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

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

void Co60()
{
  const int nrRuns = 37;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};
  double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};

  PeakPos1 = new double[nrRuns];
  PeakPos2 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];
  PeakPosErr2 = new double[nrRuns];

  PeakPos1Multi = new double[nrRuns];
  PeakPos2Multi = new double[nrRuns];
  PeakPosErr1Multi = new double[nrRuns];
  PeakPosErr2Multi = new double[nrRuns];

  for (int i = 25; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

  cout << "Mean1 (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean2 (single site) = " << TMath::Mean(nrRuns,PeakPos2) << " +- " << TMath::RMS(nrRuns,PeakPos2) << endl;
  cout << "Mean1 (multi site) = " << TMath::Mean(nrRuns,PeakPos1Multi) << " +- " << TMath::RMS(nrRuns,PeakPos1Multi) << endl;
  cout << "Mean2 (multi site) = " << TMath::Mean(nrRuns,PeakPos2Multi) << " +- " << TMath::RMS(nrRuns,PeakPos2Multi) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2,0,PeakPosErr2);
  TGraphErrors *gr3 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1Multi,0,PeakPosErr1Multi);
  TGraphErrors *gr4 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2Multi,0,PeakPosErr2Multi);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);
  gr2->SetMarkerColor(kRed);

  gr1->GetYaxis()->SetRangeUser(1150,1400);
  gr1->SetTitle("Co60 Peaks");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.6);

  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.6);
  gr4->SetMarkerColor(kRed);

  gr3->GetYaxis()->SetRangeUser(1150,1400);
  gr3->SetTitle("Co60 Peaks (multi site");
  gr3->GetXaxis()->SetTitle("run number");
  gr3->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");
  gr2->Draw("Psame");

  TCanvas *c2 = new TCanvas();
  gr3->Draw("AP");
  gr4->Draw("Psame");

  return;
}

void ProcessRun(int RunID, int ID)
{
  TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,2000);
  TH1F *hN = new TH1F("hN","Charge Energy Spectrum (purity - 2us)",200,0,2000);
  TH1F *hP = new TH1F("hP","Charge Energy Spectrum (purity + 2us)",200,0,2000);
  TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,2000);

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

  TF1 *fit = new TF1("fit",fitFunction,800,1800,6);
  fit->SetParNames("A1","E1","#sigma","R1","R2","E2");

  fit->SetParameters(650,1180,75,0.6,0.8,1370);

  fit->SetParLimits(1,1100,1250);
  fit->SetParLimits(2,10,200);
  fit->SetParLimits(3,0.2,1);
  fit->SetParLimits(4,0.2,1);
  fit->SetParLimits(5,1300,1500);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunction,800,1500,6);
  fitMulti->SetParNames("A1","E1","#sigma","R1","R2","E2");

  fitMulti->SetParameters(650,1180,75,0.6,0.8,1370);

  fitMulti->SetParLimits(1,1000,1250);
  fitMulti->SetParLimits(2,10,200);
  fitMulti->SetParLimits(3,0.5,1);
  fitMulti->SetParLimits(4,0.5,1);
  fitMulti->SetParLimits(5,1250,1400);

  hMulti->Fit("fitMulti","r");

  double par[6];
  double *parErr;
  fit->GetParameters(par);
  parErr = fit->GetParErrors();

  PeakPos1[ID] = par[1];
  PeakPos2[ID] = par[5];
  PeakPosErr1[ID] = parErr[1];
  PeakPosErr2[ID] = parErr[5];

  double parMulti[6];
  double *parErrMulti;
  fitMulti->GetParameters(parMulti);
  parErrMulti = fitMulti->GetParErrors();

  PeakPos1Multi[ID] = parMulti[1];
  PeakPos2Multi[ID] = parMulti[5];
  PeakPosErr1Multi[ID] = parErrMulti[1];
  PeakPosErr2Multi[ID] = parErrMulti[5];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << "  " << PeakPos2[ID] << " +- " << PeakPosErr2[ID] << ")  (" << PeakPos1Multi[ID] << " +- " << PeakPosErr1Multi[ID] << "  " << PeakPos2Multi[ID] << " +- " << PeakPosErr2Multi[ID] << ")" << endl;

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

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",800,180);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *gaus2 = new TF1("gaus2","[0]*TMath::Gaus(x,[1],[2])",800,1800);
  gaus2->SetParameters(par[0]*par[4],par[5],par[2]*TMath::Sqrt(par[5]/par[1]));

  TF1 *erf1 = new TF1("erf1",errf,800,1800,3);
  erf1->SetParameters(par[0]*par[3],par[1],par[2]);

  TF1 *erf2 = new TF1("erf2",errf,800,1800,3);
  erf2->SetParameters(par[0]*par[4]*par[3],par[5],par[2]*TMath::Sqrt(par[5]/par[1]));

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

  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  gaus1->Draw("same");
  gaus2->Draw("same");
  erf1->Draw("same");
  erf2->Draw("same");

  TCanvas *c2 = new TCanvas();
  hN->Draw();

  TCanvas *c3 = new TCanvas();
  hP->Draw();*/

  return;
}
