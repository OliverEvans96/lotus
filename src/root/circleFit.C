{
//=========Macro generated from canvas: c1_n106/c1_n106
//=========  (Mon Jan 18 12:25:24 2016) by ROOT version5.34/00

   gStyle->SetCanvasPreferGL(kTRUE);

   TCanvas *c1_n106 = new TCanvas("c1_n106", "c1_n106",0,0,700,500);
   gStyle->SetOptStat(0);
   c1_n106->SetHighLightColor(2);
   c1_n106->Range(0,0,1,1);
   c1_n106->SetFillColor(0);
   c1_n106->SetBorderMode(0);
   c1_n106->SetBorderSize(2);
   c1_n106->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(47);
   graph->SetName("a");
   graph->SetTitle("a");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetPoint(0,119.4395772,28.75);
   graph->SetPoint(1,73.32262862,30.25);
   graph->SetPoint(2,69.00551622,31.75);
   graph->SetPoint(3,62.74957321,33.25);
   graph->SetPoint(4,59.33270821,34.75);
   graph->SetPoint(5,55.24307096,36.25);
   graph->SetPoint(6,50.91900866,37.75);
   graph->SetPoint(7,45.83693601,39.25);
   graph->SetPoint(8,41.53389548,40.75);
   graph->SetPoint(9,33.56399455,42.25);
   graph->SetPoint(10,24.08492339,43.75);
   graph->SetPoint(11,14.57183901,45.25);
   graph->SetPoint(12,6.307831305,46.53561917);
   graph->SetPoint(13,15.22845189,45.3016202);
   graph->SetPoint(14,19.84610489,43.99796855);
   graph->SetPoint(15,23.54114692,43.70871626);
   graph->SetPoint(16,26.7204022,43.28636112);
   graph->SetPoint(17,29.55570767,43.20204742);
   graph->SetPoint(18,32.13992103,42.80140761);
   graph->SetPoint(19,34.53019411,42.15839329);
   graph->SetPoint(20,36.76473508,42.00653134);
   graph->SetPoint(21,38.87060794,41.45921015);
   graph->SetPoint(22,40.8678237,40.98083851);
   graph->SetPoint(23,42.77167829,40.36466825);
   graph->SetPoint(24,44.59417782,39.59276763);
   graph->SetPoint(25,46.3449528,39.13897401);
   graph->SetPoint(26,48.03186919,38.74886353);
   graph->SetPoint(27,49.66145082,38.14511536);
   graph->SetPoint(28,51.23917996,37.68674817);
   graph->SetPoint(29,52.76971648,36.78360812);
   graph->SetPoint(30,54.25706095,36.44518759);
   graph->SetPoint(31,55.70467839,36.1659883);
   graph->SetPoint(32,57.1155936,35.37246283);
   graph->SetPoint(33,58.49246579,34.70244038);
   graph->SetPoint(34,59.83764758,34.28867404);
   graph->SetPoint(35,61.15323238,33.89412524);
   graph->SetPoint(36,62.44109269,33.62665227);
   graph->SetPoint(37,63.70291144,32.93840411);
   graph->SetPoint(38,64.94020783,32.63738554);
   graph->SetPoint(39,66.15435881,32.52041879);
   graph->SetPoint(40,67.34661704,32.11509182);
   graph->SetPoint(41,68.5181261,32.1201284);
   graph->SetPoint(42,69.66993329,30.25);
   graph->SetPoint(43,70.80300067,30.25);
   graph->SetPoint(44,74.09829416,30.25);
   graph->SetPoint(45,76.21602774,30.25);
   graph->SetPoint(46,82.24263506,30.25);
   graph->Draw("al");
   c1_n106->Modified();
   c1_n106->cd();
   c1_n106->SetSelected(c1_n106);
}
