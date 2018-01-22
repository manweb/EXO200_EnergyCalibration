double LineFitError(double *x, double *par)
{
  double sy1 = par[7];
  double sy2 = par[8];
  double sy3 = par[9];
  double sy4 = par[10];
  double sy5 = par[11];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];
  double x5 = par[6];

  double bracketX = 1.0/5.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4) + x5/(sy5*sy5)); // [x]
  double bracketX2 = 1.0/5.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4) + x5*x5/(sy5*sy5)); // [x2]
  double bracket1 = 1.0/5.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4) + 1.0/(sy5*sy5)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/5.0 / bracket1;
  double sb2 = 1.0/5.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[0] + par[1]*x[0] + par[12]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
}

double LineFitErrorMulti(double *x, double *par)
{
  double sy1 = par[6];
  double sy2 = par[7];
  double sy3 = par[8];
  double sy4 = par[9];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];

  double bracketX = 1.0/4.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4)); // [x]
  double bracketX2 = 1.0/4.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4)); // [x2]
  double bracket1 = 1.0/4.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/4.0 / bracket1;
  double sb2 = 1.0/4.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[0] + par[1]*x[0] + par[10]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
}

double LineFitResidual(double *x, double *par)
{
  double sy1 = par[7];
  double sy2 = par[8];
  double sy3 = par[9];
  double sy4 = par[10];
  double sy5 = par[11];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];
  double x5 = par[6];

  double bracketX = 1.0/5.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4) + x5/(sy5*sy5)); // [x]
  double bracketX2 = 1.0/5.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4) + x5*x5/(sy5*sy5)); // [x2]
  double bracket1 = 1.0/5.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4) + 1.0/(sy5*sy5)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/5.0 / bracket1;
  double sb2 = 1.0/5.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[12]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2) / x[0] * 100.0;
}

double LineFitResidualMulti(double *x, double *par)
{
  double sy1 = par[8];
  double sy2 = par[9];
  double sy3 = par[10];
  double sy4 = par[11];

  // create covariance of the line fit
  double x1 = par[4];
  double x2 = par[5];
  double x3 = par[6];
  double x4 = par[7];

  double bracketX = 1.0/4.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4)); // [x]
  double bracketX2 = 1.0/4.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4)); // [x2]
  double bracket1 = 1.0/4.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/4.0 / bracket1;
  double sb2 = 1.0/4.0* bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return (par[2] + par[3]*x[0] + par[12]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2) - par[0] - par[1]*x[0]) / x[0] * 100.0;
}

void EnergyCalibration_new()
{
/*  double Etrue_single[4] = {511.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {468.03, 1203.23, 1377.2, 2751.33};
  double ErecErr_single[4] = {5.58, 8.74119, 7.51124, 15.9804};

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1702.1};
  double ErecErr_pp[1] = {12.67};
*/
/*  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  //double Erec_single[5] = {485.4, 646.75, 1202.16, 1376.847, 2745.87}; // no reclustering (fiducial volume)
  //double Erec_single[5] = {485.4, 641.76, 1196.62, 1373.14, 2759.15}; // reclustering (fiducial volume)
  //double Erec_single[5] = {463.86, 661.7, 1200.46, 1374.22, 2748.94}; // old value for 511 peak
  double Erec_single[5] = {485.4, 642.25, 1191.34, 1364.67, 2725.28}; // no reclustering (fiducial volume)
  //double ErecErr_single[5] = {4.88, 5.27, 8.47, 7.7, 15.97}; // no reclustering (fiducial volume)
  //double ErecErr_single[5] = {4.88, 5.18, 9.1, 9.06, 19.92}; // reclustering (fiducial volume)
  //double ErecErr_single[5] = {5.46, 5.011, 11.27, 9.21, 15.16}; // old value for 511 peak
  double ErecErr_single[5] = {4.88, 4.94, 3.51, 3.74, 11.05}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {584.34, 1158.36, 1315.66, 2705.29}; // gaus + erf fit
  //double Erec_multi[4] = {570.288, 1153.16, 1309.76, 2677.44}; // reclustering
  double ErecErr_multi[4] = {6.23, 8.01, 9.1, 14.18}; // gaus + erf fit
  //double ErecErr_multi[4] = {6.04123, 6.4513, 7.32739, 11.9454}; // reclustering

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1703.81};
  double ErecErr_pp[1] = {10.59};
*/
// **** New shaping times *************************************************************************

/*  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {459.88, 608.73, 1137.75, 1305.1, 2640.21}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {3.13, 2.35, 1.48, 2.62, 8.73}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {565.47, 1104.27, 1254.23, 2593.04}; // gauss + erf fit
  //double Erec_multi[4] = {567.92, 1112.38, 1263.44, 2611.12}; // gauss + erf fit (only 2 site events)
  //double Erec_multi[4] = {547.23, 1085.2, 1232.56, 2575.35}; // gauss + erf fit (only 3 site events)
  double ErecErr_multi[4] = {2.39, 4.17, 4.74, 7.24}; // gauss + erf fit
  //double ErecErr_multi[4] = {0.76, 4.13, 4.69, 6.97}; // gauss + erf fit (only 2 site events)
  //double ErecErr_multi[4] = {11.58, 2.9, 3.3, 7.88}; // gauss + erf fit (only 2 site events)

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1622.0};
  double ErecErr_pp[1] = {8.53};
*/
  // Redoing calibration with wire gain corrections (12/27/11)
/*  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {460.77, 618.46, 1143.74, 1309.29, 2645.44}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {2.27, 5.2, 3.51, 3.29, 6.42}; // no reclustering (fiducial volume)
  
  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {570.63, 1106.86, 1257.17, 2618.55}; // gauss + erf fit
  double ErecErr_multi[4] = {3.16, 2.27, 2.57, 6.56}; // gauss + erf fit
  
  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1624.26};
  double ErecErr_pp[1] = {4.97};
  */
  // Redoing calibration with wire gain corrections and measured shaping times (01/30/12)
  /*double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {465.05, 619.74, 1156.88, 1324.39, 2684.51}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {3.71, 6.84, 4.25, 4.6, 3.46}; // no reclustering (fiducial volume)
  
  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {573.77, 1122.09, 1274.47, 2636.09}; // gauss + erf fit
  double ErecErr_multi[4] = {2.74, 4.88, 5.55, 4.04}; // gauss + erf fit
  
  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1648.12};
  double ErecErr_pp[1] = {7.23};*/
  
  // Redoing calibration with wire gain corrections and measured shaping times (02/07/12)
  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {471.52, 630.38, 1168.79, 1337.6, 2706.22}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {2.19, 7.4, 2.79, 3.72, 3.44}; // no reclustering (fiducial volume)
  
  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {579.75, 1130.61, 1284.15, 2650.09}; // gauss + erf fit
  double ErecErr_multi[4] = {2.3, 4.05, 4.6, 4.5}; // gauss + erf fit
  
  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1667.23};
  double ErecErr_pp[1] = {5.95};

// ************************************************************************************************

// **** MC new shaping times **********************************************************************

/*  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {503.95, 653.27, 1160.0, 1317.28, 2582.6}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {0.405, 0.09, 0.135, 0.139, 0.161}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {652.55, 1160.91, 1318.56, 2586.53}; // gauss + erf fit
  double ErecErr_multi[4] = {0.104, 0.098, 0.112, 0.249}; // gauss + erf fit

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1583.62};
  double ErecErr_pp[1] = {0.489};
*/
// ************************************************************************************************
  
// **** New MC (01/26/11) files in /nfs/slac/g/exo/exo_data/test/MC *******************************
  
/*  double Etrue_single[5] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {646.455, 1147.14, 1301.49, 2541.13}; // no reclustering (fiducial volume)
  double ErecErr_single[5] = {0.05, 0.11, 0.13, 0.38}; // no reclustering (fiducial volume)
  
  double Etrue_multi[4] = {661.7, 1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {644.961, 1146.63, 1302.55, 2559.44}; // gauss + erf fit
  double ErecErr_multi[4] = {0.08, 0.11, 0.11, 0.18}; // gauss + erf fit
  
  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1583.62};
  double ErecErr_pp[1] = {0.489};
*/
// ************************************************************************************************

  TGraphErrors *gr1 = new TGraphErrors(5,Erec_single,Etrue_single,ErecErr_single,0); // E_true vs E_rec (single site)
  TGraphErrors *gr2 = new TGraphErrors(5,Etrue_single,Erec_single,0,ErecErr_single); // E_rec vs E_true (singel site)
  TGraphErrors *gr3 = new TGraphErrors(1,Erec_pp,Etrue_pp,ErecErr_pp,0); // double escape peak
  TGraphErrors *gr4 = new TGraphErrors(4,Erec_multi,Etrue_multi,ErecErr_multi,0); // E_true vs E_rec (multi site)
  TGraphErrors *gr5 = new TGraphErrors(4,Etrue_multi,Erec_multi,0,ErecErr_multi); // E_rec vs E_true (multi site)

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetXaxis()->SetTitle("rec energy");
  gr1->GetYaxis()->SetTitle("true energy");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetXaxis()->SetTitle("true energy");
  gr2->GetYaxis()->SetTitle("rec energy");

  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.6);
  gr3->SetMarkerColor(kRed);

  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.6);

  gr4->GetXaxis()->SetTitle("rec energy");
  gr4->GetYaxis()->SetTitle("true energy");

  gr5->SetMarkerStyle(20);
  gr5->SetMarkerSize(0.6);

  gr5->GetXaxis()->SetTitle("true energy");
  gr5->GetYaxis()->SetTitle("rec energy");

// **** Single Site ***************************************************************************************
// Fit E_rec vs E_true
  TF1 *fit = new TF1("fit","[0]+[1]*x",0,3000);
  fit->SetLineWidth(1);

  gr2->Fit("fit","r");

  double par[2];
  fit->GetParameters(par);

  double a = 1.0/par[1];
  double b = -1.0*a*par[0];

  TF1 *ErrLine1_1 = new TF1("ErrLine1_1",LineFitError,0,3000,13);

  ErrLine1_1->SetParameter(0,par[0]);
  ErrLine1_1->SetParameter(1,par[1]);
  ErrLine1_1->SetParameter(2,Etrue_single[0]);
  ErrLine1_1->SetParameter(3,Etrue_single[1]);
  ErrLine1_1->SetParameter(4,Etrue_single[2]);
  ErrLine1_1->SetParameter(5,Etrue_single[3]);
  ErrLine1_1->SetParameter(6,Etrue_single[4]);
  ErrLine1_1->SetParameter(7,ErecErr_single[0]);
  ErrLine1_1->SetParameter(8,ErecErr_single[1]);
  ErrLine1_1->SetParameter(9,ErecErr_single[2]);
  ErrLine1_1->SetParameter(10,ErecErr_single[3]);
  ErrLine1_1->SetParameter(11,ErecErr_single[4]);
  ErrLine1_1->SetParameter(12,1);

  ErrLine1_1->SetLineWidth(1);

  TF1 *ErrLine2_1 = new TF1("ErrLine2_1",LineFitError,0,3000,13);

  ErrLine2_1->SetParameter(0,par[0]);
  ErrLine2_1->SetParameter(1,par[1]);
  ErrLine2_1->SetParameter(2,Etrue_single[0]);
  ErrLine2_1->SetParameter(3,Etrue_single[1]);
  ErrLine2_1->SetParameter(4,Etrue_single[2]);
  ErrLine2_1->SetParameter(5,Etrue_single[3]);
  ErrLine2_1->SetParameter(6,Etrue_single[4]);
  ErrLine2_1->SetParameter(7,ErecErr_single[0]);
  ErrLine2_1->SetParameter(8,ErecErr_single[1]);
  ErrLine2_1->SetParameter(9,ErecErr_single[2]);
  ErrLine2_1->SetParameter(10,ErecErr_single[3]);
  ErrLine2_1->SetParameter(11,ErecErr_single[4]);
  ErrLine2_1->SetParameter(12,-1);

  ErrLine2_1->SetLineWidth(1);

// Fit E_true vs E_rec
  TF1 *LineFit = new TF1("LineFit","[0]+[1]*x",0,3000);

  LineFit->SetLineWidth(1);

  LineFit->SetParameter(0,b);
  LineFit->SetParameter(1,a);

  gr1->Fit("LineFit","r");
  
  gMinuit->mnmatu(1);

  double par2[2];
  LineFit->GetParameters(par2);

  TF1 *ErrLine1 = new TF1("ErrLine1",LineFitError,0,3000,13);

  ErrLine1->SetParameter(0,par2[0]);
  ErrLine1->SetParameter(1,par2[1]);
  ErrLine1->SetParameter(2,Etrue_single[0]);
  ErrLine1->SetParameter(3,Etrue_single[1]);
  ErrLine1->SetParameter(4,Etrue_single[2]);
  ErrLine1->SetParameter(5,Etrue_single[3]);
  ErrLine1->SetParameter(6,Etrue_single[4]);
  ErrLine1->SetParameter(7,ErecErr_single[0]*par2[1]);
  ErrLine1->SetParameter(8,ErecErr_single[1]*par2[1]);
  ErrLine1->SetParameter(9,ErecErr_single[2]*par2[1]);
  ErrLine1->SetParameter(10,ErecErr_single[3]*par2[1]);
  ErrLine1->SetParameter(11,ErecErr_single[4]*par2[1]);
  ErrLine1->SetParameter(12,1);

  ErrLine1->SetLineWidth(1);

  TF1 *ErrLine2 = new TF1("ErrLine2",LineFitError,0,3000,13);

  ErrLine2->SetParameter(0,par2[0]);
  ErrLine2->SetParameter(1,par2[1]);
  ErrLine2->SetParameter(2,Etrue_single[0]);
  ErrLine2->SetParameter(3,Etrue_single[1]);
  ErrLine2->SetParameter(4,Etrue_single[2]);
  ErrLine2->SetParameter(5,Etrue_single[3]);
  ErrLine2->SetParameter(6,Etrue_single[4]);
  ErrLine2->SetParameter(7,ErecErr_single[0]*par2[1]);
  ErrLine2->SetParameter(8,ErecErr_single[1]*par2[1]);
  ErrLine2->SetParameter(9,ErecErr_single[2]*par2[1]);
  ErrLine2->SetParameter(10,ErecErr_single[3]*par2[1]);
  ErrLine2->SetParameter(11,ErecErr_single[4]*par2[1]);
  ErrLine2->SetParameter(12,-1);

  ErrLine2->SetLineWidth(1);

// **** Multi site ****************************************************************************************
// Fit E_rec vs E_true
  TF1 *fitMulti = new TF1("fitMulti","[0]+[1]*x",0,3000);
  fitMulti->SetLineWidth(1);

  gr5->Fit("fitMulti","r");

  double parMulti[2];
  fitMulti->GetParameters(parMulti);

  double aMulti = 1.0/parMulti[1];
  double bMulti = -1.0*aMulti*parMulti[0];

  TF1 *ErrLineMulti1_1 = new TF1("ErrLineMulti1_1",LineFitErrorMulti,0,3000,11);

  ErrLineMulti1_1->SetParameter(0,parMulti[0]);
  ErrLineMulti1_1->SetParameter(1,parMulti[1]);
  ErrLineMulti1_1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti1_1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti1_1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti1_1->SetParameter(5,Etrue_multi[3]);
  ErrLineMulti1_1->SetParameter(6,ErecErr_multi[0]);
  ErrLineMulti1_1->SetParameter(7,ErecErr_multi[1]);
  ErrLineMulti1_1->SetParameter(8,ErecErr_multi[2]);
  ErrLineMulti1_1->SetParameter(9,ErecErr_multi[3]);
  ErrLineMulti1_1->SetParameter(10,1);

  ErrLineMulti1_1->SetLineWidth(1);

  TF1 *ErrLineMulti2_1 = new TF1("ErrLineMulti2_1",LineFitErrorMulti,0,3000,11);

  ErrLineMulti2_1->SetParameter(0,parMulti[0]);
  ErrLineMulti2_1->SetParameter(1,parMulti[1]);
  ErrLineMulti2_1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti2_1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti2_1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti2_1->SetParameter(5,Etrue_multi[3]);
  ErrLineMulti2_1->SetParameter(6,ErecErr_multi[0]);
  ErrLineMulti2_1->SetParameter(7,ErecErr_multi[1]);
  ErrLineMulti2_1->SetParameter(8,ErecErr_multi[2]);
  ErrLineMulti2_1->SetParameter(9,ErecErr_multi[3]);
  ErrLineMulti2_1->SetParameter(10,-1);

  ErrLineMulti2_1->SetLineWidth(1);

// Fit E_true vs E_rec
  TF1 *LineFitMulti = new TF1("LineFitMulti","[0]+[1]*x",0,3000);

  LineFitMulti->SetLineWidth(1);

  LineFitMulti->SetParameter(0,bMulti);
  LineFitMulti->SetParameter(1,aMulti);

  gr4->Fit("LineFit","r");
  
  gMinuit->mnmatu(1);

  double parMulti2[2];
  LineFitMulti->GetParameters(parMulti2);

  TF1 *ErrLineMulti1 = new TF1("ErrLineMulti1",LineFitErrorMulti,0,3000,11);

  ErrLineMulti1->SetParameter(0,parMulti2[0]);
  ErrLineMulti1->SetParameter(1,parMulti2[1]);
  ErrLineMulti1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti1->SetParameter(5,Etrue_multi[3]);
  ErrLineMulti1->SetParameter(6,ErecErr_multi[0]*parMulti2[1]);
  ErrLineMulti1->SetParameter(7,ErecErr_multi[1]*parMulti2[1]);
  ErrLineMulti1->SetParameter(8,ErecErr_multi[2]*parMulti2[1]);
  ErrLineMulti1->SetParameter(9,ErecErr_multi[3]*parMulti2[1]);
  ErrLineMulti1->SetParameter(10,1);

  ErrLineMulti1->SetLineWidth(1);

  TF1 *ErrLineMulti2 = new TF1("ErrLineMulti2",LineFitErrorMulti,0,3000,11);

  ErrLineMulti2->SetParameter(0,parMulti2[0]);
  ErrLineMulti2->SetParameter(1,parMulti2[1]);
  ErrLineMulti2->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti2->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti2->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti2->SetParameter(5,Etrue_multi[3]);
  ErrLineMulti2->SetParameter(6,ErecErr_multi[0]*parMulti2[1]);
  ErrLineMulti2->SetParameter(7,ErecErr_multi[1]*parMulti2[1]);
  ErrLineMulti2->SetParameter(8,ErecErr_multi[2]*parMulti2[1]);
  ErrLineMulti2->SetParameter(9,ErecErr_multi[3]*parMulti2[1]);
  ErrLineMulti2->SetParameter(10,-1);

  ErrLineMulti2->SetLineWidth(1);

// **** Residuals ***********************************************************************************************
// Single site
  // calculate residual for each data point
  double r1 = (Erec_single[0] - fit->Eval(Etrue_single[0]))/fit->Eval(Etrue_single[0]) * 100.0;
  double r2 = (Erec_single[1] - fit->Eval(Etrue_single[1]))/fit->Eval(Etrue_single[1]) * 100.0;
  double r3 = (Erec_single[2] - fit->Eval(Etrue_single[2]))/fit->Eval(Etrue_single[2]) * 100.0;
  double r4 = (Erec_single[3] - fit->Eval(Etrue_single[3]))/fit->Eval(Etrue_single[3]) * 100.0;
  double r5 = (Erec_single[4] - fit->Eval(Etrue_single[4]))/fit->Eval(Etrue_single[4]) * 100.0;

  double rErr1 = ErecErr_single[0]/fit->Eval(Etrue_single[0]) * 100.0;
  double rErr2 = ErecErr_single[1]/fit->Eval(Etrue_single[1]) * 100.0;
  double rErr3 = ErecErr_single[2]/fit->Eval(Etrue_single[2]) * 100.0;
  double rErr4 = ErecErr_single[3]/fit->Eval(Etrue_single[3]) * 100.0;
  double rErr5 = ErecErr_single[4]/fit->Eval(Etrue_single[4]) * 100.0;

  double residual[5] = {r1, r2, r3, r4, r5};
  double residualErr[5] = {rErr1, rErr2, rErr3, rErr4, rErr5};
  TGraphErrors *grResidual_single = new TGraphErrors(5,Etrue_single,residual,0,residualErr);

  grResidual_single->SetMarkerStyle(20);
  grResidual_single->SetMarkerSize(0.8);
  grResidual_single->SetLineStyle(21);

  grResidual_single->GetXaxis()->SetTitle("true energy [keV]");
  grResidual_single->GetYaxis()->SetTitle("residual [%]");

  grResidual_single->GetYaxis()->SetRangeUser(-8.0,8.0);

// Multi site
  // calculate residual for each data point
  double rMulti1 = (Erec_multi[0] - fit->Eval(Etrue_multi[0]))/fit->Eval(Etrue_multi[0]) * 100.0;
  double rMulti2 = (Erec_multi[1] - fit->Eval(Etrue_multi[1]))/fit->Eval(Etrue_multi[1]) * 100.0;
  double rMulti3 = (Erec_multi[2] - fit->Eval(Etrue_multi[2]))/fit->Eval(Etrue_multi[2]) * 100.0;
  double rMulti4 = (Erec_multi[3] - fit->Eval(Etrue_multi[3]))/fit->Eval(Etrue_multi[3]) * 100.0;

  double rErrMulti1 = ErecErr_multi[0]/fit->Eval(Etrue_multi[0]) * 100.0;
  double rErrMulti2 = ErecErr_multi[1]/fit->Eval(Etrue_multi[1]) * 100.0;
  double rErrMulti3 = ErecErr_multi[2]/fit->Eval(Etrue_multi[2]) * 100.0;
  double rErrMulti4 = ErecErr_multi[3]/fit->Eval(Etrue_multi[3]) * 100.0;

  double residualMulti[4] = {rMulti1, rMulti2, rMulti3, rMulti4};
  double residualErrMulti[4] = {rErrMulti1, rErrMulti2, rErrMulti3, rErrMulti4};
  TGraphErrors *grResidual_multi = new TGraphErrors(4,Etrue_multi,residualMulti,0,residualErrMulti);

  grResidual_multi->SetMarkerStyle(22);
  grResidual_multi->SetMarkerSize(0.8);
  grResidual_multi->SetLineStyle(21);

// Double escape Peak
  double r_pp = (Erec_pp[0] - fit->Eval(Etrue_pp[0]))/fit->Eval(Etrue_pp[0]) * 100.0;
  double rErr_pp = ErecErr_pp[0]/fit->Eval(Etrue_pp[0]) * 100.0;

  double residual_pp[1] = {r_pp};
  double residualErr_pp[1] = {rErr_pp};
  TGraphErrors *grResidual_pp = new TGraphErrors(1,Etrue_pp,residual_pp,0,residualErr_pp);

  grResidual_pp->SetMarkerStyle(21);

// Error band single site
  TF1 *ResidualLine1 = new TF1("ResidualLine1",LineFitResidual,0,3000,13);

  ResidualLine1->SetParameter(0,par[0]);
  ResidualLine1->SetParameter(1,par[1]);
  ResidualLine1->SetParameter(2,Etrue_single[0]);
  ResidualLine1->SetParameter(3,Etrue_single[1]);
  ResidualLine1->SetParameter(4,Etrue_single[2]);
  ResidualLine1->SetParameter(5,Etrue_single[3]);
  ResidualLine1->SetParameter(6,Etrue_single[4]);
  ResidualLine1->SetParameter(7,ErecErr_single[0]);
  ResidualLine1->SetParameter(8,ErecErr_single[1]);
  ResidualLine1->SetParameter(9,ErecErr_single[2]);
  ResidualLine1->SetParameter(10,ErecErr_single[3]);
  ResidualLine1->SetParameter(11,ErecErr_single[4]);
  ResidualLine1->SetParameter(12,1);

  ResidualLine1->SetLineWidth(1);

  TF1 *ResidualLine2 = new TF1("ResidualLine2",LineFitResidual,0,3000,13);

  ResidualLine2->SetParameter(0,par[0]);
  ResidualLine2->SetParameter(1,par[1]);
  ResidualLine2->SetParameter(2,Etrue_single[0]);
  ResidualLine2->SetParameter(3,Etrue_single[1]);
  ResidualLine2->SetParameter(4,Etrue_single[2]);
  ResidualLine2->SetParameter(5,Etrue_single[3]);
  ResidualLine2->SetParameter(6,Etrue_single[4]);
  ResidualLine2->SetParameter(7,ErecErr_single[0]);
  ResidualLine2->SetParameter(8,ErecErr_single[1]);
  ResidualLine2->SetParameter(9,ErecErr_single[2]);
  ResidualLine2->SetParameter(10,ErecErr_single[3]);
  ResidualLine2->SetParameter(11,ErecErr_single[4]);
  ResidualLine2->SetParameter(12,-1);

  ResidualLine2->SetLineWidth(1);

// Error band multi site
  TF1 *ResidualLine1Multi = new TF1("ResidualLine1Multi",LineFitResidualMulti,0,3000,13);

  ResidualLine1Multi->SetParameter(0,par[0]);
  ResidualLine1Multi->SetParameter(1,par[1]);
  ResidualLine1Multi->SetParameter(2,parMulti[0]);
  ResidualLine1Multi->SetParameter(3,parMulti[1]);
  ResidualLine1Multi->SetParameter(4,Etrue_multi[0]);
  ResidualLine1Multi->SetParameter(5,Etrue_multi[1]);
  ResidualLine1Multi->SetParameter(6,Etrue_multi[2]);
  ResidualLine1Multi->SetParameter(7,Etrue_multi[3]);
  ResidualLine1Multi->SetParameter(8,ErecErr_multi[0]);
  ResidualLine1Multi->SetParameter(9,ErecErr_multi[1]);
  ResidualLine1Multi->SetParameter(10,ErecErr_multi[2]);
  ResidualLine1Multi->SetParameter(11,ErecErr_multi[3]);
  ResidualLine1Multi->SetParameter(12,1);

  ResidualLine1Multi->SetLineWidth(1);
  ResidualLine1Multi->SetLineStyle(2);

  TF1 *ResidualLine2Multi = new TF1("ResidualLine2Multi",LineFitResidualMulti,0,3000,13);

  ResidualLine2Multi->SetParameter(0,par[0]);
  ResidualLine2Multi->SetParameter(1,par[1]);
  ResidualLine2Multi->SetParameter(2,parMulti[0]);
  ResidualLine2Multi->SetParameter(3,parMulti[1]);
  ResidualLine2Multi->SetParameter(4,Etrue_multi[0]);
  ResidualLine2Multi->SetParameter(5,Etrue_multi[1]);
  ResidualLine2Multi->SetParameter(6,Etrue_multi[2]);
  ResidualLine2Multi->SetParameter(7,Etrue_multi[3]);
  ResidualLine2Multi->SetParameter(8,ErecErr_multi[0]);
  ResidualLine2Multi->SetParameter(9,ErecErr_multi[1]);
  ResidualLine2Multi->SetParameter(10,ErecErr_multi[2]);
  ResidualLine2Multi->SetParameter(11,ErecErr_multi[3]);
  ResidualLine2Multi->SetParameter(12,-1);

  ResidualLine2Multi->SetLineWidth(1);
  ResidualLine2Multi->SetLineStyle(2);

  gr4->SetMarkerColor(kBlue);

  TCanvas *c1 = new TCanvas();
  gr2->Draw("AP");
  ErrLine1_1->Draw("same");
  ErrLine2_1->Draw("same");

  TCanvas *c2 = new TCanvas();
  gr1->Draw("AP");
  gr3->Draw("Psame");
  gr4->Draw("Psame");
  //LineFit->Draw("same");
  ErrLine1->Draw("same");
  ErrLine2->Draw("same");

  TCanvas *c3 = new TCanvas();
  grResidual_single->Draw("AP");
  grResidual_multi->Draw("Psame");
  grResidual_pp->Draw("Psame");
  ResidualLine1->Draw("same");
  ResidualLine2->Draw("same");
  ResidualLine1Multi->Draw("same");
  ResidualLine2Multi->Draw("same");

  grResidual_single->GetXaxis()->SetLimits(0,3000);
  c3->Update();

  TCanvas *c4 = new TCanvas();
  gr5->Draw("AP");
  ErrLineMulti1_1->Draw("same");
  ErrLineMulti2_1->Draw("same");

  TCanvas *c5 = new TCanvas();
  gr4->Draw("AP");
  ErrLineMulti1->Draw("same");
  ErrLineMulti2->Draw("same");

  return;
}
