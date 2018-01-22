double xcl[1000];
double ycl[1000];
double zcl[1000];
double ucl[1000];
double vcl[1000];
double dtcl[1000];
double ercl[1000];
double epcl[1000];
double epcl2[1000];
double epcl_p[1000];
double epcl_n[1000];
double epcl_gc[1000];
double epcl_p_gc[1000];
double epcl_n_gc[1000];
double tcl[1000];
int dhalfcl[1000];
int clusters[1000];
int nUWires[1000];
int nremoved;
bool fiducial;

void ProcessRun(int RunID, int ID);

double HardwareGains[114] = {365.9,367.7,330.8,323,328.1,318.9,324.9,317.9,302.7,333.4,368.6,357.3,321.1,359.3,354.8,385.9,366.5,324.3,333.8,325.9,328.6,328.7,327.4,327.9,326.9,328.8,329.6,328.2,321.8,332,328.7,328.9,321,317,322.1,326.4,331.4,329.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,323.8,318.3,325.5,337.4,340,328.2,328.5,311.3,311.3,309.6,298.9,308.1,312.8,313.7,314.6,321.8,309,311.1,319,334,322.6,325.7,323.8,326.8,331.5,331.1,332.3,336.8,331,333.4,333.4,329.3,325.6,316.7,301.2,324.6,328.6,333.1};

//double GainCorrection[114] = {1.0038,1.0021,1.0131,1.0022,0.99502,1.008,1.0001,1.0001,0.99318,0.99944,0.99956,0.98753,0.99143,0.9838,0.99642, 1.002,0.99631,1.0014,1.0016,0.98875,0.99517,1.0097,1.0105,0.99809,1.0124,1.0147,1.0055,0.9976,1.0085,1.0133,1.0102,1.0109,1.0112, 1.0047,1.0012,1.0077,0.99451,1.0106,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.98721,1.0069, 1.0157,1.0112,1.0126,1.0095,1.0092,0.98987,0.99004,0.99311,0.97877,0.98696,0.98721,0.98486,0.98922,0.98128,0.98474,0.98069,0.97353, 0.99019,0.99332,0.99796,1.0033,1.0064,1.0118,1.0067,1.0027,1.0043,1.0027,1.0042,1.0005,1.0033,0.99836,1.0009,1.0017,1.008,0.99767,1.0111};

double GainCorrection[114] ={1.0063,0.99631,0.95355,0.98478,1.0008,0.97925,0.97718,0.97029,0.95922,0.9997,1.0057,1.003,0.99045,1.0057,1.0056,1.0054,0.94392,1.0269,1.012,1.0023,0.99201,1.0132,1.003,1.0087,1.0061,1.0144,0.98746,0.99365,0.99455,0.9931,1.0019,0.97818,0.99212,0.97474,0.89061,0.99739,0.99313,0.90054,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0075,1.0083,1.0017,1.0243,1.0395,1.0274,1.0304,1.0064,0.99619,1.027,0.99832,1.0113,0.99871,1.0122,1.0006,0.98523,0.99921,1.0162,1.0269,0.99464,0.97529,0.99449,0.99837,0.98738,1.0005,0.98692,1.0051,0.99625,0.98746,1.028,1.0159,1.0279,1.0015,1.0108,1.0035,1.009,0.99671,0.99012};

void GetData()
{
  // Co60 runs
  //const int nrRuns = 37;
  //int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};
  //double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};

  // Th228 runs
  //const int nrRuns = 41;
  //int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940};
  //double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940};

  // Cs137 runs
  //const int nrRuns = 7;
  //int runID[nrRuns] = {2109, 2110, 2111, 2112, 2113, 2117, 2119};
  //double runIDCopy[nrRuns] = {2109, 2110, 2111, 2112, 2113, 2117, 2119};

  //const int nrRuns = 85;
  //int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965, 1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940, 2109, 2110, 2111, 2112, 2113, 2117, 2119};
  //double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965, 1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1881, 1889, 1894, 1915, 1916, 1917, 1918, 1921, 1922, 1923, 1924, 1925, 1926, 1929, 1930, 1931, 1932, 1937, 1940, 2109, 2110, 2111, 2112, 2113, 2117, 2119};

  // runs with new shaping times
  const int nrRuns = 100;
  int runID[nrRuns] = {2410, 2412, 2415, 2416, 2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2449, 2450, 2469, 2473, 2479, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2674, 2683, 2689, 2708, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865, 2866, 2867, 2883, 2884, 2885, 2886, 2895, 2897, 2898, 2904, 2923, 2927, 2933, 2938, 2943, 2947, 2955, 2966, 2971, 2981, 2986, 2991, 2995, 3007, 3018, 3019, 3024, 3028, 3032, 3034, 3048, 3053, 3057, 3074, 3091, 3099, 3109};
  double runIDCopy[nrRuns] = {2410, 2412, 2415, 2416, 2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2449, 2450, 2469, 2473, 2479, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2674, 2683, 2689, 2708, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865, 2866, 2867, 2883, 2884, 2885, 2886, 2895, 2897, 2898, 2904, 2923, 2927, 2933, 2938, 2943, 2947, 2955, 2966, 2971, 2981, 2986, 2991, 2995, 3007, 3018, 3019, 3024, 3028, 3032, 3034, 3048, 3053, 3057, 3074, 3091, 3099, 3109};

/*  const int nrRuns = 12;
  int runID[nrRuns] = {2424,2426,2431,2432,2433,2434,2448,2469,2473,2496,2543,2578};
  double runIDCopy[nrRuns] = {2424,2426,2431,2432,2433,2434,2448,2469,2473,2496,2543,2578};
*/
  for (int i = 0; i < 100; i++) {
     ProcessRun(runID[i],i);
  }

  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  if (RunID < 10000) {sprintf(FileName,"/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/alpha/%i/recon0000%i-*.root",RunID,RunID);}
  else {sprintf(FileName,"/EXO200Data/Disk3/processed/%i/recon0000%i-*.root",RunID,RunID);}

  //char FileName[100];
  //if (RunID < 10000) {sprintf(FileName,"/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/alpha/apdposcorr/run%i_scintposcorr.root",RunID);}
  //else {sprintf(FileName,"/EXO200Data/Disk3/processed/%i/recon0000%i-*.root",RunID,RunID);}

  //char *FileName = "/EXO200Data/Disk1/MC/ReDigitized/NewFitWindow/Cs137/MC*.root";
  //char *FileName = "/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/alpha/2424/recon00002424-*.root";

  cout << "Processing run " << RunID << "..." << endl;

  double ecraw;
  double ecrec;
  double epcrec;
  double ecrec_p;
  double ecrec_n;
  double ecrec_gc;
  double ecrec_p_gc;
  double ecrec_n_gc;
  double csc;
  double cscraw;
  int algsc;
  double xc;
  double yc;
  double zc;
  double tc;
  int nsite;
  int nWires;

  TTree *tOut = new TTree("t","tree");
  tOut->Branch("ecraw",&ecraw,"ecraw/D");
  tOut->Branch("ecrec",&ecrec,"ecrec/D");
  tOut->Branch("epcrec",&epcrec,"epcrec/D");
  tOut->Branch("ecrec_p",&ecrec_p,"ecrec_p/D");
  tOut->Branch("ecrec_n",&ecrec_n,"ecrec_n/D");
  tOut->Branch("ecrec_gc",&ecrec_gc,"ecrec_gc/D");
  tOut->Branch("ecrec_p_gc",&ecrec_p_gc,"ecrec_p_gc/D");
  tOut->Branch("ecrec_n_gc",&ecrec_n_gc,"ecrec_n_gc/D");
  tOut->Branch("csc",&csc,"csc/D");
  tOut->Branch("cscraw",&cscraw,"cscraw/D");
  tOut->Branch("algsc",&algsc,"algsc/I");
  tOut->Branch("xc",&xc,"xc/D");
  tOut->Branch("yc",&yc,"yc/D");
  tOut->Branch("zc",&zc,"zc/D");
  tOut->Branch("tc",&tc,"tc/D");
  tOut->Branch("nsite",&nsite,"nsite/I");
  tOut->Branch("nWires",&nWires,"nWires/I");

  TChain *t = new TChain("tree");
  t->Add(FileName);

  EXOEventData *ED = 0;
  t->SetBranchAddress("EventBranch",&ED);

  double elife = 0.0;

  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (int evtID = 0; evtID < nentries; evtID++) {
     t->GetEntry(evtID);
     if (evtID%1000 == 0) {cout << evtID << " events processed (elife = " << elife << ")" << endl;}

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
        if (purTime >= 81.6 && purTime < 94.0) {
           purFitP0 = -9011.55;
           purFitP1 = 115.417;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 94.0 && purTime < 102.5) {
           purFitP0 = 2000.0;
           purFitP1 = 0.0;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 102.5 && purTime < 113.0) {
           purFitP0 = -1208000.0;
           purFitP1 = 34380.0;
           purFitP2 = -325.9;
           purFitP3 = 1.03;
           purFitP4 = 0.0;
        }
        if (purTime >= 113.0 && purTime < 129.6) {
           purFitP0 = -48740.0;
           purFitP1 = 805.0;
           purFitP2 = -3.259;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 129.6 && purTime < 142.0) {
           purFitP0 = -29510.0;
           purFitP1 = 230.1;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }
        if (purTime >= 142.0) {
           purFitP0 = 3300.0;
           purFitP1 = 0.0;
           purFitP2 = 0.0;
           purFitP3 = 0.0;
           purFitP4 = 0.0;
        }

     elife = purFitP4*purTime*purTime*purTime*purTime + purFitP3*purTime*purTime*purTime + purFitP2*purTime*purTime + purFitP1*purTime + purFitP0;

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
        int ncl = scint_cluster->GetNumChargeClusters();

        csc = scint_cluster->fRawEnergy;
        cscraw = scint_cluster->fCountsSumOnAPDPlaneOne + scint_cluster->fCountsSumOnAPDPlaneTwo;
        algsc = scint_cluster->fAlgorithmUsed;

        // safe all charge clusters in array
        for (int clID = 0; clID < ncl; clID++) {
           EXOChargeCluster *charge_cluster = scint_cluster->GetChargeClusterAt(clID);

           xcl[clID] = charge_cluster->fX;
           ycl[clID] = charge_cluster->fY;
           zcl[clID] = charge_cluster->fZ;
           ucl[clID] = charge_cluster->fU;
           vcl[clID] = charge_cluster->fV;
           dtcl[clID] = charge_cluster->fDriftTime / 1000.0;
           ercl[clID] = charge_cluster->fCorrectedEnergy;
           epcl[clID] = charge_cluster->fCorrectedEnergy * TMath::Exp(dtcl[clID]/elife);
           epcl_p[clID] = charge_cluster->fCorrectedEnergy * TMath::Exp(dtcl[clID]/(elife+2.0));
           epcl_n[clID] = charge_cluster->fCorrectedEnergy * TMath::Exp(dtcl[clID]/(elife-2.0));
           epcl2[clID] = charge_cluster->fPurityCorrectedEnergy;
           tcl[clID] = charge_cluster->fCollectionTime;
           dhalfcl[clID] = charge_cluster->fDetectorHalf;
           clusters[clID] = 1;

           int nU = charge_cluster->GetNumUWireSignals();
           epcl_gc[clID] = 0.0;
           epcl_p_gc[clID] = 0.0;
           epcl_n_gc[clID] = 0.0;
           for (int UID = 0; UID < nU; UID++) {
              EXOUWireSignal *u_signal = charge_cluster->GetUWireSignalAt(UID);

              int UCHID = u_signal->fChannel;
              double U_corr = 1.0;
              double HW_gain = 300.0;
              if (UCHID >= 0 && UCHID < 114) {HW_gain = HardwareGains[UCHID] / 300.0; U_corr = GainCorrection[UCHID];}
              //double eU = u_signal->fRawEnergy * HW_gain * U_corr;
              double eU = u_signal->fCorrectedEnergy * U_corr;
              epcl_gc[clID] += eU * TMath::Exp(dtcl[clID]/elife);
              epcl_p_gc[clID] += eU * TMath::Exp(dtcl[clID]/(elife+2.0));
              epcl_n_gc[clID] += eU * TMath::Exp(dtcl[clID]/(elife-2.0));
           }
           nUWires[clID] = nU;
        }

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
        ecrec_p = 0.0;
        ecrec_n = 0.0;
        ecrec_gc = 0.0;
        ecrec_p_gc = 0.0;
        ecrec_n_gc = 0.0;
        nWires = -999;
        xc = -999;
        yc = -999;
        zc = -999;
        tc = -999;
        for (int i = 0; i < ncl; i++) {
           if (clusters[i] == 0) {continue;}
           ecraw += ercl[i];
           ecrec += epcl[i];
           epcrec += epcl2[i];
           ecrec_p += epcl_p[i];
           ecrec_n += epcl_n[i];
           ecrec_gc += epcl_gc[i];
           ecrec_p_gc += epcl_p_gc[i];
           ecrec_n_gc += epcl_n_gc[i];
           if (TMath::Sqrt(xcl[i]*xcl[i] + ycl[i]*ycl[i]) > 163) {fiducial = false; continue;}
           if (zcl[i] > 152 || zcl[i] < -152) {fiducial = false; continue;}
           if (zcl[i] > -20 && zcl[i] < 20) {fiducial = false; continue;}
           //if (epcl[i] < 100.0) {fiducial = false; continue;} // for Cs137 every charge cluster must at least have 100 keV
           //if (nUWires[i] != 1) {fiducial = false; continue;} // take only events where all clusters have one u wire hit
        }
        if (ncl == 0) {fiducial = false;}

        nsite = ncl - nremoved;
        if (nsite == 1) {nWires = nUWires[0]; xc = xcl[0]; yc = ycl[0]; zc = zcl[0]; tc = dtcl[0];}
        if (!fiducial) {continue;}
        if (nsite == 0) {continue;}
        tOut->Fill();
     }
  }
  //char oFile[100];
  //sprintf(oFile,"../analysis/V5/%i_noReclustering_Fiducial.root",RunID);

  char oFile[100];
  sprintf(oFile,"/nfs/slac/g/exo/maweber/EnergyCalibration/analysis/alpha/V3/%i_noReclustering_Fiducial.root",RunID);

  //char *oFile = "../analysis/V2/2424_noReclustering_Fiducial.root";
  TFile *f = new TFile(oFile,"RECREATE");
  tOut->Write();
  f->Close();

  // clean up
  tOut->Delete();
  t->Delete();

  return;
}
