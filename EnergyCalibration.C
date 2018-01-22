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

void EnergyCalibration()
{
  double Etrue_single[5] = {511.0, 661.7, 1173.2, 1332.5, 2614.0};
  double Erec_single[5] = {468.03, 646.67, 1203.23, 1377.2, 2751.33};
  double ErecErr_single[5] = {5.58, 5.011, 8.74119, 7.51124, 15.9804};

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1702.1};
  double ErecErr_pp[1] = {12.67};

  TGraphErrors *gr1 = new TGraphErrors(5,Erec_single,Etrue_single,ErecErr_single,0); // E_true vs E_rec
  TGraphErrors *gr2 = new TGraphErrors(5,Etrue_single,Erec_single,0,ErecErr_single); // E_rec vs E_true
  TGraphErrors *gr3 = new TGraphErrors(1,Erec_pp,Etrue_pp,ErecErr_pp,0);

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

  TF1 *fit = new TF1("fit","[0]+[1]*x",0,3000);
  fit->SetLineWidth(1);

  gr2->Fit("fit","r");

  double par[2];
  fit->GetParameters(par);

  double a = 1.0/par[1];
  double b = -1.0*a*par[0];

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

  double r_pp = (Erec_pp[0] - fit->Eval(Etrue_pp[0]))/fit->Eval(Etrue_pp[0]) * 100.0;
  double rErr_pp = ErecErr_pp[0]/fit->Eval(Etrue_pp[0]) * 100.0;

  double residual_pp[1] = {r_pp};
  double residualErr_pp[1] = {rErr_pp};
  TGraphErrors *grResidual_pp = new TGraphErrors(1,Etrue_pp,residual_pp,0,residualErr_pp);

  grResidual_pp->SetMarkerStyle(21);

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

  TF1 *LineFit = new TF1("LineFit","[0]+[1]*x",0,3000);

  LineFit->SetLineWidth(1);

  LineFit->SetParameter(0,b);
  LineFit->SetParameter(1,a);

  gr1->Fit("LineFit","r");

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

  TCanvas *c1 = new TCanvas();
  gr2->Draw("AP");
  ErrLine1_1->Draw("same");
  ErrLine2_1->Draw("same");

  TCanvas *c2 = new TCanvas();
  gr1->Draw("AP");
  gr3->Draw("Psame");
  //LineFit->Draw("same");
  ErrLine1->Draw("same");
  ErrLine2->Draw("same");

  TCanvas *c3 = new TCanvas();
  grResidual_single->Draw("AP");
  grResidual_pp->Draw("Psame");
  ResidualLine1->Draw("same");
  ResidualLine2->Draw("same");

  grResidual_single->GetXaxis()->SetLimits(0,3000);
  c3->Update();

  return;
}
