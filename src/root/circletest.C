{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jan 18 15:08:12 2016) by ROOT version5.34/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   c1->SetHighLightColor(2);
   c1->Range(-7.947988,-187.5,77.58692,187.5);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(36);
   graph->SetName("Graph");
   graph->SetTitle("circlepoints.txt");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(0.5);
   graph->SetPoint(0,63.3311,34.75);
   graph->SetPoint(1,59.1755,36.25);
   graph->SetPoint(2,56.2773,37.75);
   graph->SetPoint(3,52.4506,39.25);
   graph->SetPoint(4,47.6442,40.75);
   graph->SetPoint(5,43.2403,42.25);
   graph->SetPoint(6,39.5162,43.75);
   graph->SetPoint(7,32.6287,45.25);
   graph->SetPoint(8,26.0901,46.75);
   graph->SetPoint(9,18.7912,48.25);
   graph->SetPoint(10,11.8577,49.75);
   graph->SetPoint(11,6.30783,50.5381);
   graph->SetPoint(12,15.2285,48.9737);
   graph->SetPoint(13,19.8461,48.0589);
   graph->SetPoint(14,23.5411,47.3785);
   graph->SetPoint(15,26.7204,46.7849);
   graph->SetPoint(16,29.5557,45.8427);
   graph->SetPoint(17,32.1399,45.2828);
   graph->SetPoint(18,34.5302,44.744);
   graph->SetPoint(19,36.7647,44.2241);
   graph->SetPoint(20,38.8706,43.7517);
   graph->SetPoint(21,40.8678,43.174);
   graph->SetPoint(22,42.7717,42.5022);
   graph->SetPoint(23,44.5942,41.7486);
   graph->SetPoint(24,46.345,41.6708);
   graph->SetPoint(25,48.0319,40.6886);
   graph->SetPoint(26,49.6615,40.2257);
   graph->SetPoint(27,51.2392,39.7846);
   graph->SetPoint(28,52.7697,38.8816);
   graph->SetPoint(29,54.2571,38.5548);
   graph->SetPoint(30,55.7047,37.934);
   graph->SetPoint(31,57.1156,37.4675);
   graph->SetPoint(32,58.4925,36.562);
   graph->SetPoint(33,59.8376,35.9542);
   graph->SetPoint(34,61.1532,35.7557);
   graph->SetPoint(35,62.4411,34.589);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","circlepoints.txt",100,0.605503,69.03343);
   Graph_Graph1->SetMinimum(-150);
   Graph_Graph1->SetMaximum(150);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ap");
   
   TEllipse *ellipse = new TEllipse(-15.0544,-140.464,60.38148,60.38148,0,360,0);
   ellipse->SetFillStyle(0);
   ellipse->SetLineWidth(3);
   ellipse->Draw();
   
   TPaveText *pt = new TPaveText(0.3671552,0.9365254,0.6328448,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("circlepoints.txt");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
