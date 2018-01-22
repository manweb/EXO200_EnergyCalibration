 using namespace RooFit;

void RooTest()
{
  RooRealVar x("x","x",0,-10,10);
  RooRealVar xHist("xHist","ecrec",0,3500);
  RooRealVar mean("mean","Gaussian mean",0);
  RooRealVar sigma("sigma","Gaussian sigma",65,0,200);
  RooRealVar mean2("mean2","Gaussian mean",0);
  RooRealVar sigma2("sigma2","Gaussian sigma",65,0,200);

  RooGaussian gauss("gauss","Gaussian pdf",xHist,mean,sigma);
  RooGaussian gauss2("gauss2","Gaussian pdf",xHist,mean2,sigma2);

  // create fit function
/*  RooRealVar mean1("mean1","mean of first peak",1170,1100,1250);
  RooRealVar mean2("mean2","mean of second peak",1315,1250,1400);
  RooRealVar sigma1("sigma1","sigma of first peak",60,20,200);
  RooRealVar sigma2("sigma2","sigma of second peak",60,20,200);
  RooRealVar r1("r1","fraction 1",0.05,0.0,1.0);
  RooRealVar r2("r2","fraction 2",0.05,0.0,1.0);
  RooRealVar r("r","fraction",0.5,0.0,1.0);

  RooGaussian gauss1("gauss1","gauss of first peak",xHist,mean1,sigma1);
  RooGaussian gauss2("gauss2","gauss of second peak",xHist,mean2,sigma2);
  RooGenericPdf erf1("erf1","0.5*(1.0+TMath::Erfc((xHist-mean1)/(sqrt(2)*sigma1)))",RooArgList(xHist,mean1,sigma1));
  RooGenericPdf erf2("erf2","0.5*(1.0+TMath::Erfc((xHist-mean2)/(sqrt(2)*sigma2)))",RooArgList(xHist,mean2,sigma2));

  RooAddPdf peak1("peak1","peak 1",RooArgList(gauss1,erf1),r1);
  RooAddPdf peak2("peak2","peak 2",RooArgList(gauss2,erf2),r2);
  RooAddPdf fit("fit","fit function",RooArgList(peak1,peak2),r);
*/
  // load data from file
  TChain *t = new TChain("t","tree");
  t->Add("../../analysis/MC/NewFitWindow/Th228_noReclustering_Fiducial.root");

  TH1F *h1 = new TH1F("h1","single site spectrum",200,0,3500);
  t->Draw("ecraw>>h1","nsite == 1");

  TH1F *h12 = new TH1F("h12","multi site spectrum",200,0,3500);
  t->Draw("ecraw>>h12","nsite > 1");

  // load data from file
  TChain *t2 = new TChain("t","tree");
  t2->Add("../../analysis/V4/2424_noReclustering_Fiducial.root");

  TH1F *h2 = new TH1F("h2","single site spectrum",200,0,3500);
  t2->Draw("ecrec*0.967+71.64>>h2","nsite == 1");

  TH1F *h22 = new TH1F("h22","multi site spectrum",200,0,3500);
  t2->Draw("ecrec*0.9634+116.5>>h22","nsite > 1");

  RooDataHist RealData("RealData","Real data (single site)",xHist,h2);
  RooDataHist RealData2("RealData2","Real data (multi site)",xHist,h22);

  //RooRealVar xHist("xHist","histogram axis",0,2500);
  RooDataHist data("data","single site spectrum",xHist,h1);
  RooDataHist data2("data2","multi site spectrum",xHist,h12);

  RooHistPdf dataPDF("dataPDF","pdf from data",xHist,data);
  RooHistPdf data2PDF("data2PDF","pdf from data2",xHist,data2);

  RooFFTConvPdf dxg("dxg","data convolved with gaussian",xHist,dataPDF,gauss);
  RooFFTConvPdf dxg2("dxg2","data2 convolved with gaussian",xHist,data2PDF,gauss2);

  dxg.fitTo(RealData, Range(720,3000));
  dxg2.fitTo(RealData2, Range(720,3000));

  //mean.Print();
  //sigma.Print();

  RooPlot *histFrame = xHist.frame();
  RealData.plotOn(histFrame, MarkerSize(0.6));
  dxg.plotOn(histFrame);
  //fit.plotOn(histFrame);
  //fit.plotOn(histFrame, Components("erf1"), LineStyle(kDashed));
  //fit.plotOn(histFrame, Components("erf2"), LineStyle(kDashed));
  //gauss.plotOn(histFrame);

  TCanvas *c1 = new TCanvas();
  histFrame->Draw();

  RooPlot *histFrame2 = xHist.frame();
  RealData2.plotOn(histFrame2, MarkerSize(0.6));
  dxg2.plotOn(histFrame2);

  TCanvas *c2 = new TCanvas();
  histFrame2->Draw();

/*  RooPlot *xframe = x.frame();
  gauss.plotOn(xframe);
  xframe->Draw();
*/
  return;
}
