double MC_erf_p0 = 120.8;
double MC_erf_p1 = 33.18;
double Data_erf_p0 = 165.4;
double Data_erf_p1 = 68.65;

void GenerateTree()
{
  TRandom3 *rand = new TRandom3(0);
  
  TF1 *fit = new TF1("fit","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit->SetParameters(MC_erf_p0,MC_erf_p1);
  
  int n = 100000;
  
  double e;
  double vc;
  
  TTree *t = new TTree("t","tree");
  t->Branch("e",&e,"e/D");
  t->Branch("vc",&vc,"vc/D");
  
  for (int i = 0; i < n; i++) {
    e = 1000*rand->Uniform();
    
    double r = rand->Uniform();
    if (r > fit->Eval(e)) {vc = -999;}
    else {vc = 0;}
    
    t->Fill();
  }
  
  TFile *f = new TFile("VThresholdTest.root","RECREATE");
  t->Write();
  f->Close();
  
  return;
}

void PlotEfficiency()
{
  TChain *t = new TChain("t","tree");
  
  t->Add("VThresholdTest.root");
  
  t->Draw("e>>h(200,0,1000)","","goff");
  t->Draw("e>>h2(200,0,1000)","vc != -999","goff");
  
  TH1F *hDiff = h2->Clone("hDiff");
  hDiff->Sumw2();
  hDiff->Divide(h);
  
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(0.8);
  
  TF1 *fit = new TF1("fit","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit->SetParameters(120,30);
  
  fit->SetLineColor(kBlue);
  fit->SetLineWidth(2);
  
  hDiff->Fit("fit","r");
  
  TF1 *fit2 = new TF1("fit2","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit2->SetParameters(Data_erf_p0,Data_erf_p1);
  
  fit2->SetLineColor(kRed);
  fit2->SetLineWidth(2);
  fit2->SetLineStyle(2);
  
  TCanvas *c1 = new TCanvas();
  hDiff->Draw();
  fit2->Draw("same");
  
  return;
}

void Tweak()
{
  TChain *t = new TChain("t","tree");
  
  t->Add("VThresholdTest.root");
  
  TRandom3 *rand = new TRandom3(0);
  
  double e;
  double vc;
  
  t->SetBranchAddress("e",&e);
  t->SetBranchAddress("vc",&vc);
  
  TH1F *h = new TH1F("h","h",200,0,1000);
  TH1F *h2 = new TH1F("h2","h2",200,0,1000);
  
  TF1 *fit = new TF1("fit","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit->SetParameters(MC_erf_p0,MC_erf_p1);
  
  fit->SetLineColor(kBlue);
  fit->SetLineWidth(2);
  fit->SetLineStyle(2);
  
  TF1 *fit2 = new TF1("fit2","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit2->SetParameters(Data_erf_p0,Data_erf_p1);
  
  fit2->SetLineColor(kRed);
  fit2->SetLineWidth(2);
  
  int n = t->GetEntries();
  for (int i = 0; i < n; i++) {
    t->GetEntry(i);
    
    h->Fill(e);
    if (vc == -999) {continue;}
    
    double r = fit->Eval(e)*rand->Uniform();
    if (r > fit2->Eval(e) && r < fit->Eval(e)) {continue;}
    
    h2->Fill(e);
  }
  
  TH1F *hDiff = h2->Clone("hDiff");
  hDiff->Sumw2();
  hDiff->Divide(h);
  
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(0.8);
  
  TF1 *fit3 = new TF1("fit3","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit3->SetParameters(160,60);
  
  fit3->SetLineColor(kBlue);
  fit3->SetLineWidth(2);
  
  hDiff->Fit("fit3","r");
  
  TCanvas *c1 = new TCanvas();
  hDiff->Draw();
  fit->Draw("same");
  fit2->Draw("same");
}

void TweakPlot()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  TChain *t = new TChain("t");
  
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/VThresh/0_noReclustering_Fiducial_erf_new.root");
  
  TF1 *fit = new TF1("fit","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit->SetParameters(MC_erf_p0,MC_erf_p1);
  
  fit->SetLineColor(kBlue);
  fit->SetLineWidth(2);
  fit->SetLineStyle(2);
  
  TF1 *fit2 = new TF1("fit2","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit2->SetParameters(Data_erf_p0,Data_erf_p1);
  
  fit2->SetLineColor(kRed);
  fit2->SetLineWidth(2);
  
  t->Draw("ecraw>>h(200,0,3000)","nsite == 1 && nWires <= 2","goff");
  t->Draw("ecraw>>h2(200,0,3000)","nsite == 1 && nWires <= 2 && vc != -999","goff");
  
  TH1F *hDiff = h2->Clone("hDiff");
  hDiff->Sumw2();
  hDiff->Divide(h);
  
  hDiff->SetMarkerStyle(20);
  hDiff->SetMarkerSize(0.8);
  
  TF1 *fit3 = new TF1("fit3","0.5*(1.0 + TMath::Erf((x - [0]) / (TMath::Sqrt(2) * [1])))",0,1000);
  
  fit3->SetParameters(160,60);
  
  fit3->SetLineColor(kBlue);
  fit3->SetLineWidth(2);
  
  hDiff->Fit("fit3","r");
  
  TCanvas *c1 = new TCanvas();
  hDiff->Draw();
  fit->Draw("same");
  fit2->Draw("same");
}