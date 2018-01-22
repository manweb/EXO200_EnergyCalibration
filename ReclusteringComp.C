void ReclusteringComp()
{
   double Etrue[4] = {661.7, 1173.2, 1332.5, 2614.0};
   double SingleSite_0us[4] = {622.65, 1147.75, 1315.25, 2658.49};
   double SingleSite_1us[4] = {618.89, 1143.74, 1310.89, 2650.57};
   double SingleSite_2us[4] = {616.22, 1140.58, 1307.72, 2645.95};
   double SingleSite_3us[4] = {614.05, 1138.39, 1305.69, 2642.81};
   double SingleSite_5us[4] = {611.73, 1136.0, 1303.46, 2639.59};
   double SingleSite_7us[4] = {610.32, 1134.76, 1302.29, 2638.71};
   double SingleSite_9us[4] = {609.3, 1134.55, 1302.35, 2638.6};

   double MultiSite_0us[4] = {574.5, 1110.96, 1261.83, 2629.11};
   double MultiSite_1us[4] = {570.53, 1105.11, 1255.19, 2622.78};
   double MultiSite_2us[4] = {568.35, 1106.77, 1257.07, 2618.2};
   double MultiSite_3us[4] = {566.28, 1105.19, 1255.27, 2614.82};
   double MultiSite_5us[4] = {562.09, 1102.32, 1252.01, 2610.32};
   double MultiSite_7us[4] = {560.52, 1100.41, 1249.85, 2607.01};
   double MultiSite_9us[4] = {560.37, 1098.93, 1248.16, 2604.51};

   double Residual_0us[4];
   double Residual_1us[4];
   double Residual_2us[4];
   double Residual_3us[4];
   double Residual_5us[4];
   double Residual_7us[4];
   double Residual_9us[4];

   for (int i = 0; i < 5; i++) {
      Residual_0us[i] = (SingleSite_0us[i] - MultiSite_0us[i]) / SingleSite_0us[i] * 100.0;
      Residual_1us[i] = (SingleSite_1us[i] - MultiSite_1us[i]) / SingleSite_1us[i] * 100.0;
      Residual_2us[i] = (SingleSite_2us[i] - MultiSite_2us[i]) / SingleSite_2us[i] * 100.0;
      Residual_3us[i] = (SingleSite_3us[i] - MultiSite_3us[i]) / SingleSite_3us[i] * 100.0;
      Residual_5us[i] = (SingleSite_5us[i] - MultiSite_5us[i]) / SingleSite_5us[i] * 100.0;
      Residual_7us[i] = (SingleSite_7us[i] - MultiSite_7us[i]) / SingleSite_7us[i] * 100.0;
      Residual_9us[i] = (SingleSite_9us[i] - MultiSite_9us[i]) / SingleSite_9us[i] * 100.0;
   }

   TGraphErrors *gr0us = new TGraphErrors(4,Etrue,Residual_0us,0,0);
   TGraphErrors *gr1us = new TGraphErrors(4,Etrue,Residual_1us,0,0);
   TGraphErrors *gr2us = new TGraphErrors(4,Etrue,Residual_2us,0,0);
   TGraphErrors *gr3us = new TGraphErrors(4,Etrue,Residual_3us,0,0);
   TGraphErrors *gr5us = new TGraphErrors(4,Etrue,Residual_5us,0,0);
   TGraphErrors *gr7us = new TGraphErrors(4,Etrue,Residual_7us,0,0);
   TGraphErrors *gr9us = new TGraphErrors(4,Etrue,Residual_9us,0,0);

   gr0us->SetTitle("single - multi site residual");
   gr0us->GetXaxis()->SetTitle("true energy [keV]");
   gr0us->GetYaxis()->SetTitle("residual [%]");

   gr0us->SetMarkerStyle(20);
   gr1us->SetMarkerStyle(21);
   gr2us->SetMarkerStyle(22);
   gr3us->SetMarkerStyle(23);
   gr5us->SetMarkerStyle(24);
   gr7us->SetMarkerStyle(25);
   gr9us->SetMarkerStyle(26);

   gr0us->SetMarkerSize(0.6);
   gr1us->SetMarkerSize(0.6);
   gr2us->SetMarkerSize(0.6);
   gr3us->SetMarkerSize(0.6);
   gr5us->SetMarkerSize(0.6);
   gr7us->SetMarkerSize(0.6);
   gr9us->SetMarkerSize(0.6);
  
   TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
   l->AddEntry(gr0us,"0 #mus","lp");
   l->AddEntry(gr1us,"1 #mus","lp");
   l->AddEntry(gr2us,"2 #mus","lp");
   l->AddEntry(gr3us,"3 #mus","lp");
   l->AddEntry(gr5us,"5 #mus","lp");
   l->AddEntry(gr7us,"7 #mus","lp");
   l->AddEntry(gr9us,"9 #mus","lp");

   TCanvas *c1 = new TCanvas();
   gr0us->Draw("AP");
   gr1us->Draw("Psame");
   gr2us->Draw("Psame");
   gr3us->Draw("Psame");
   gr5us->Draw("Psame");
   gr7us->Draw("Psame");
   gr9us->Draw("Psame");
   l->Draw("same");

   return;
}
