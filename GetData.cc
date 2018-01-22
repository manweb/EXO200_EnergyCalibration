#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TH1.h>
#include <TGraphErrors.h>

#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOChargeCluster.hh"
#include "EXOUtilities/EXOScintillationCluster.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include "EXOUtilities/EXO3DView.hh"

using namespace ROOT::Math;

double xcl[1000];
double ycl[1000];
double zcl[1000];
double ucl[1000];
double vcl[1000];
double dtcl[1000];
double ercl[1000];
double eccl[1000];
double epcl[1000];
double tcl[1000];
double uWireRatio[1000];
int dhalfcl[1000];
int clusters[1000];
int nUWires[1000];
int IsCulled[1000];
int nremoved;
bool fiducial;

/*double charge_p0_1U = 65.05;
double charge_p1_1U = 0.9266;
double charge_p0_2U = 81.0;
double charge_p1_2U = 0.9347;
double charge_p0_3U = 105.6;
double charge_p1_3U = 0.9378;*/

double charge_p0_1U = 0.0;
double charge_p1_1U = 1.0;
double charge_p0_2U = 0.0;
double charge_p1_2U = 1.0;
double charge_p0_3U = 0.0;
double charge_p1_3U = 1.0;

bool ApplyThreshold = false;
bool ApplyHexCulling = true;
bool ApplyFiducialCut = true;
bool isMC = false;

TRandom3 *rand3 = new TRandom3(0);

//double MC_erf_p0 = 164.9;
//double MC_erf_p1 = 36.36;
//double Data_erf_p0 = 247.8;
//double Data_erf_p1 = 67.12;

// new values
double MC_erf_p0 = 120.8;
double MC_erf_p1 = 33.18;
double Data_erf_p0 = 165.4;
double Data_erf_p1 = 68.65;

// test
//double MC_erf_p0 = 120.8;
//double MC_erf_p1 = 33.18;
//double Data_erf_p0 = 200.0;
//double Data_erf_p1 = 60.0;

//TH1F *h = new TH1F("h","h",200,0,3000);
//TH1F *hCulled = new TH1F("hCulled","hCulled",200,0,3000);
//TH2F *h2Culled = new TH2F("h2Culled","h2Culled",200,-200,200,200,-200,200);

void ProcessRun(int RunID);
bool HexagonalCut(double x, double y, double u, double v, double R);

int main(int argc, char* argv[])
{
  if (argc != 2) {std::cout << "Incorrect number of arguments" << std::endl; return 0;}

  int RunID = atoi(argv[1]);

  ProcessRun(RunID);

  return 1;
}

void ProcessRun(int RunID)
{
  //char *FileName = Form("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/%i/recon0000%i-*.root",RunID,RunID);
  char *FileName = Form("/nfs/slac/g/exo/cgd8d/LowFieldStudy_June/Run%i/repro0000%i-*.root",RunID,RunID);

  //char *FileName = "/nfs/slac/g/exo/MC/P2_SourceP4_px_228Th/P2_SourceP4_px_228Th_*.root";
  //char *FileName = "/nfs/slac/g/exo/MC/P2_SourceP4_px_60Co/P2_SourceP4_px_60Co_*.root";
  //char *FileName = "/nfs/slac/g/exo/MC/P2_SourceP4_px_137Cs/P2_SourceP4_px_137Cs_*.root";
  
  // open purity file
  //TFile *f1 = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/scripts/20120320_purity_tpc1.root","READ");
  //TFile *f2 = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/scripts/20120320_purity_tpc2.root","READ");
  
  //TGraphErrors *gr_pur_tpc1 = (TGraphErrors*)f1->Get("Graph");
  //TGraphErrors *gr_pur_tpc2 = (TGraphErrors*)f2->Get("Graph");
  
  //f1->Close();
  //f2->Close();
  
  //delete f1;
  //delete f2;

  std::cout << "Processing run " << RunID << "..." << std::endl;

  double ecraw;
  double ecrec;
  double epcrec;
  double csc;
  double cscraw;
  int algsc;
  double xc;
  double yc;
  double zc;
  double uc;
  double vc;
  double tc;
  double R;
  int nsite;
  int nWires;

  TTree *tOut = new TTree("t","tree");
  tOut->Branch("ecraw",&ecraw,"ecraw/D");
  tOut->Branch("ecrec",&ecrec,"ecrec/D");
  tOut->Branch("epcrec",&epcrec,"epcrec/D");
  tOut->Branch("csc",&csc,"csc/D");
  tOut->Branch("cscraw",&cscraw,"cscraw/D");
  tOut->Branch("algsc",&algsc,"algsc/I");
  tOut->Branch("xc",&xc,"xc/D");
  tOut->Branch("yc",&yc,"yc/D");
  tOut->Branch("zc",&zc,"zc/D");
  tOut->Branch("uc",&uc,"uc/D");
  tOut->Branch("vc",&vc,"vc/D");
  tOut->Branch("tc",&tc,"tc/D");
  tOut->Branch("R",&R,"R/D");
  tOut->Branch("nsite",&nsite,"nsite/I");
  tOut->Branch("nWires",&nWires,"nWires/I");

  TChain *t = new TChain("tree");
  t->Add(FileName);

  EXOEventData *ED = 0;
  t->SetBranchAddress("EventBranch",&ED);
  
  double elife = 0.0;

  int nentries = t->GetEntries();
  std::cout << "nentries = " << nentries << std::endl;
  int count = 0;
  for (int evtID = 0; evtID < nentries; evtID++) {
     t->GetEntry(evtID);
     if (evtID%1000 == 0) {std::cout << evtID << " events processed (elife = " << elife << ")" << std::endl;}

     int nsc = ED->GetNumScintillationClusters();
     for (int scID = 0; scID < nsc; scID++) {
        EXOScintillationCluster *scint_cluster = ED->GetScintillationCluster(scID);

        double tsc1 = scint_cluster->fTime;

        bool GoodScintCluster = true;
        // make sure no other scintillation cluster is within 110us
        for (int i = 0; i < nsc; i++) {
           if (i == scID) {continue;}

           double tsc2 = ED->GetScintillationCluster(i)->fTime;
           if (TMath::Abs(tsc2 - tsc1) < 110000) {GoodScintCluster = false; break;}
        }

        if (!GoodScintCluster) {continue;}
        if (scint_cluster->fTime > 1928000) {continue;}
        if (scint_cluster->fTime < 1000) {continue;}
        int ncl = scint_cluster->GetNumChargeClusters();

        csc = scint_cluster->fRawEnergy;
        cscraw = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;
        algsc = scint_cluster->fAlgorithmUsed;

        // safe all charge clusters in array
        nsite = 0;
        bool HasMissingPositions = false;
        for (int clID = 0; clID < ncl; clID++) {
           EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);
          
           double E_calib = charge_cluster->fRawEnergy;
          
           double Perf_data = 0.5 * (1.0 + TMath::Erf((E_calib - Data_erf_p0) / (TMath::Sqrt(2) * Data_erf_p1)));
           double Perf_MC = 0.5 * (1.0 + TMath::Erf((E_calib - MC_erf_p0) / (TMath::Sqrt(2) * MC_erf_p1)));
          
           double rn = Perf_MC*rand3->Uniform();

           if (charge_cluster->fX == -999 || charge_cluster->fY == -999) {HasMissingPositions = true; break;}
           //if (isMC && rn > Perf_data) {HasMissingPositions = true; break;}
           if (isMC && rn > Perf_data) {charge_cluster->fX = -999; charge_cluster->fY = -999; charge_cluster->fV = -999;}

           //if (TMath::Sqrt(charge_cluster->fX * charge_cluster->fX + charge_cluster->fY * charge_cluster->fY) > 163.0) {continue;}

           if (HexagonalCut(charge_cluster->fX,charge_cluster->fY,charge_cluster->fU,charge_cluster->fV,163.0) && ApplyHexCulling) {continue;}

           //if (charge_cluster->fDetectorHalf == 0 || charge_cluster->fDetectorHalf == 2) {elife = gr_pur_tpc1->Eval(ED->fEventHeader.fTriggerSeconds);}
           //if (charge_cluster->fDetectorHalf == 1 || charge_cluster->fDetectorHalf == 3) {elife = gr_pur_tpc2->Eval(ED->fEventHeader.fTriggerSeconds);}

           xcl[nsite] = charge_cluster->fX;
           ycl[nsite] = charge_cluster->fY;
           zcl[nsite] = charge_cluster->fZ;
           ucl[nsite] = charge_cluster->fU;
           vcl[nsite] = charge_cluster->fV;
           dtcl[nsite] = charge_cluster->fDriftTime / 1000.0;
           ercl[nsite] = charge_cluster->fRawEnergy;
           eccl[nsite] = charge_cluster->fCorrectedEnergy;
           epcl[nsite] = charge_cluster->fPurityCorrectedEnergy;
           tcl[nsite] = charge_cluster->fCollectionTime;
           dhalfcl[nsite] = charge_cluster->fDetectorHalf;
           clusters[nsite] = 1;
          
           size_t nU = charge_cluster->GetNumUWireSignals();
           std::vector<int> UChannels;
           for (size_t i = 0; i < nU; i++) {
             int UCH = ((EXOUWireSignal*)charge_cluster->GetUWireSignalAt(i))->fChannel;
             if (std::find(UChannels.begin(), UChannels.end(), UCH) == UChannels.end()) {UChannels.push_back(UCH);}
           }
          
           nUWires[nsite] = UChannels.size();
          
           uWireRatio[nsite] = -999.9;
           if (nUWires[nsite] == 2) {
             double uEnergy1 = ((EXOUWireSignal*)charge_cluster->GetUWireSignalAt(0))->fCorrectedEnergy;
             double uEnergy2 = ((EXOUWireSignal*)charge_cluster->GetUWireSignalAt(1))->fCorrectedEnergy;
             if (uEnergy2 == 0.0) {continue;}
             
             uWireRatio[nsite] = uEnergy1 / uEnergy2;
           }

           nsite++;
        }

        if (HasMissingPositions) {continue;}

        // recluster
        nremoved = 0;
/*        for (int i = 0; i < ncl; i++) {
           if (clusters[i] == 0) {continue;}
           for (int k = i+1; k < ncl; k++) {
              if (clusters[k] == 0) {continue;}
              if (TMath::Abs(ucl[i] - ucl[k]) < 10.0 && TMath::Abs(tcl[i] - tcl[k]) < 10000.0 && dhalfcl[i] == dhalfcl[k]) {
                 epcl[i] += epcl[k];
                 epcl_p[i] += epcl_p[k];
                 epcl_n[i] += epcl_n[k];
                 epcl_gc[i] += epcl_gc[k];
                 epcl_p_gc[i] += epcl_p_gc[k];
                 epcl_n_gc[i] += epcl_n_gc[k];
                 nUWires[i] += nUWires[k];
                 xcl[i] = (xcl[i] + xcl[k]) / 2.0;
                 ycl[i] = (ycl[i] + ycl[k]) / 2.0;
                 zcl[i] = (zcl[i] + zcl[k]) / 2.0;
                 clusters[k] = 0;
                 nremoved++;
              }
           }
        }
*/
        //for (int i = 0; i < ncl; i++) {
        //   if (ercl[i] < 100.0) {clusters[i] = 0; nremoved++;}
        //}

        // apply fiducial cut
        fiducial = true;
        ecraw = 0.0;
        ecrec = 0.0;
        epcrec = 0.0;
        nWires = -999;
        xc = -999;
        yc = -999;
        zc = -999;
        uc = -999;
        vc = -999;
        tc = -999;
        for (int i = 0; i < nsite; i++) {
           if (clusters[i] == 0) {continue;}
           ecraw += ercl[i];
           ecrec += eccl[i];
           if (nUWires[i] == 1) {epcrec += epcl[i]*charge_p1_1U + charge_p0_1U;}
           if (nUWires[i] == 2) {epcrec += epcl[i]*charge_p1_2U + charge_p0_2U;}
           if (nUWires[i] > 2) {epcrec += epcl[i]*charge_p1_3U + charge_p0_3U;}
           R = uWireRatio[i];

           //if (IsCulled[i] == 1) {hCulled->Fill(ercl[i]); h2Culled->Fill(xcl[i],ycl[i]);}
           //h->Fill(ercl[i]);
           //count++;
          
           if (!ApplyFiducialCut) {continue;}
          
           //epcrec += epcl[i] * Correction;
           if (TMath::Sqrt(xcl[i]*xcl[i] + ycl[i]*ycl[i]) > 163) {fiducial = false; continue;}
           //if (ucl[i] > 163.0) {fiducial = false; continue;}
           if (zcl[i] > 172 || zcl[i] < -172) {fiducial = false; continue;}
           if (zcl[i] > -20 && zcl[i] < 20) {fiducial = false; continue;}
           //if (epcl[i] < 100.0) {fiducial = false; continue;} // for Cs137 every charge cluster must at least have 100 keV
           //if (nUWires[i] != 1) {fiducial = false; continue;} // take only events where all clusters have one u wire hi 
        }
        if (ncl == 0) {fiducial = false;}

        //nsite = ncl - nremoved;
        if (nsite == 1) {nWires = nUWires[0]; xc = xcl[0]; yc = ycl[0]; zc = zcl[0]; uc = ucl[0]; vc = vcl[0]; tc = dtcl[0];}
        if (!fiducial) {continue;}
        if (nsite == 0) {continue;}
        tOut->Fill();
     }

     //if (evtID  > 400000) {break;}
  }

  //TH1F *hDiff = (TH1F*)hCulled->Clone("hDiff");
  //hDiff->Sumw2();
  //hDiff->Divide(h);
  //h2Culled->Scale(1.0/count);

  //char *oFile = Form("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V15/%i_noReclustering_Fiducial.root",RunID);
  char *oFile = Form("%i_noReclustering_Fiducial_HexCull.root",RunID);
  
  //char *oFile = "/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/MC/200412/P2_SourceP4_px_228Th.root";
  //char *oFile = "/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/MC/200412/P2_SourceP4_px_60Co.root";
  //char *oFile = "/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/MC/200412/P2_SourceP4_px_137Cs.root";

  TFile *f = new TFile(oFile,"RECREATE");
  tOut->Write();
  f->Close();
  
  // clean up
  tOut->Delete();
  t->Delete();
  
  return;
}

bool HexagonalCut(double x, double y, double u, double v, double R)
{
  double R_teflon = 183.0;
  double RHex = R;
  
  // move all coordinates into first quadrant
  /*double x_TMP = TMath::Abs(x);
  double y_TMP = TMath::Abs(y);
  
  double r = TMath::Sqrt(x*x + y*y);
  
  // cull everything outside the teflon radius
  if (r > R_teflon) {return true;}
  
  // apply quadratic cut first
  if (x_TMP > R || y_TMP > R / TMath::Cos(TMath::Pi() / 6.0)) {return true;}
  
  // now cut out the edges
  double slope = -1.0*TMath::Tan(TMath::Pi() / 3.0);
  double offset = (TMath::Tan(TMath::Pi() / 6.0) + TMath::Tan(TMath::Pi() / 3.0))*R;
  if (y_TMP > slope*x_TMP + offset) {return true;}*/

  double r=sqrt(x*x+y*y);

  // cull everything outside teflon or outside choosen hexagon
  if ( r > R_teflon ||
       fabs(u) > RHex ||
       fabs(v) > RHex ||
       fabs(x) > RHex ) {
    return true;
  }
  
  return false;
}
