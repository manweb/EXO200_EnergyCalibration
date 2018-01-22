{
  gSystem->Load("libEXOUtilities");
  gROOT->SetStyle("Plain");
  
  double E[1000];
  int nsite = 0;
  
  TChain *t = new TChain("tree");
  
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2448/recon00002448-*.root");
  /*t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2719/recon00002719-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2725/recon00002725-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2732/recon00002732-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2737/recon00002737-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2748/recon00002748-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2766/recon00002766-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2771/recon00002771-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2785/recon00002785-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2804/recon00002804-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2811/recon00002811-*.root");
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/2817/recon00002817-*.root");
*/  
  EXOEventData *ED = 0;
  
  t->SetBranchAddress("EventBranch",&ED);
  
  TH1F *h = new TH1F("h","h",100,0,150);
  
  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (int i = 0; i < nentries; i++) {
    if (i%1000 == 0) {cout << i << " events processed" << endl;}
    t->GetEntry(i);
    
    int nsc = ED->GetNumScintillationClusters();
    for (int scID = 0; scID < nsc; scID++) {
      EXOScintillationCluster *sc = ED->GetScintillationCluster(scID);
      
      int ncl = sc->GetNumChargeClusters();
      double Etotal = 0.0;
      bool fiducial = true;
      nsite = 0;
      for (int clID = 0; clID < ncl; clID++) {
        EXOChargeCluster *cc = sc->GetChargeClusterAt(clID);
        
        if (TMath::Sqrt(cc->fX*cc->fX + cc->fY*cc->fY) > 163.0) {fiducial = false; break;}
        if (TMath::Abs(cc->fZ) > 172.0 || TMath::Abs(cc->fZ) < 20.0) {fiducial = false; break;}
        
        E[clID] = cc->fRawEnergy;
        Etotal += cc->fRawEnergy;
        nsite++;
      }

      if (nsite < 2) {continue;}
      if (!fiducial) {continue;}
      if (Etotal < 1000.0) {continue;}
      for (int k = 0; k < nsite; k++) {if (E[k] < 150.0) {h->Fill(E[k]);}}
    }
  }
  
  h->Draw("EZP");
}
