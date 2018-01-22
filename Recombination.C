void Recombination()
{
  TChain *t = new TChain("tree","t");

  t->Add("/EXO200Data/Disk5/processed/2424/recon00002424-*.root");

  EXOEventData *ED = 0;
  t->SetBranchAddress("EventBranch",&ED);

  TH1F *hChargeSingle = new TH1F("hChargeSingle","single site charge spectrum",200,0,3500);
  TH1F *hChargeMulti = new TH1F("hChargeMulti","multi site charge spectrum",200,0,3500);
  TH1F *hScintSingle = new TH1F("hScintSingle","single site scintillation spectrum",200,0,20000);
  TH1F *hScintMulti = new TH1F("hScintMulti","multo site scintillation spectrum",200,0,20000);

  int nentries = t->GetEntries();
  cout << "number of events: " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
     t->GetEntry(i);

     if (i%1000 == 0) {cout << i << " entries processed" << endl;}
     int nsc = ED->GetNumScintillationClusters();

     for (int scID = 0; scID < nsc; scID++) {
        EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(scID);

        double ecrec = 0.0;
        int ncl = scint_cluster->GetNumChargeClusters();
        for (int clID = 0; clID < ncl; clID++) {
           EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);

           double cX = charge_cluster->fX;
           double cY = charge_cluster->fY;
           double cZ = charge_cluster->fZ;

           bool fiducial = true;
           if (TMath::Sqrt(cX*cX + cY*cY) > 80) {fiducial = false; continue;}
           if (TMath::Abs(cZ) > 142) {fiducial = false; continue;}
           if (TMath::Abs(cZ) < 50) {fiducial = false; continue;}

           ecrec += charge_cluster->fCorrectedEnergy * TMath::Exp(charge_cluster->fDriftTime / 1000.0 / 5000.0);
        }

        if (!fiducial) {continue;}
        double csc = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;

        if (ncl == 1) {hChargeSingle->Fill(ecrec); hScintSingle->Fill(csc);}
        if (ncl > 1) {hChargeMulti->Fill(ecrec); hScintMulti->Fill(csc);}
     }
  }

  hChargeMulti->SetLineColor(kRed);
  hScintMulti->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas();
  hChargeSingle->Draw();
  hChargeMulti->Draw("same");

  TCanvas *c2 = new TCanvas();
  hScintSingle->Draw();
  hScintMulti->Draw("same");

  return;
}
