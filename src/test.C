{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Mar 31 18:53:43 2016) by ROOT version5.34/00

   gStyle->SetCanvasPreferGL(kTRUE);

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   gStyle->SetOptStat(0);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(25);
   graph->SetName("a");
   graph->SetTitle("a");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetPoint(0,1,6.677567263);
   graph->SetPoint(1,3,6.638419672);
   graph->SetPoint(2,5,6.762065847);
   graph->SetPoint(3,7,6.955666264);
   graph->SetPoint(4,9,6.942639713);
   graph->SetPoint(5,11,6.999792414);
   graph->SetPoint(6,13,7.158171997);
   graph->SetPoint(7,15,7.414726814);
   graph->SetPoint(8,17,7.626986043);
   graph->SetPoint(9,19,8.239227474);
   graph->SetPoint(10,21,7.894423957);
   graph->SetPoint(11,23,8.376321095);
   graph->SetPoint(12,25,8.283345532);
   graph->SetPoint(13,27,8.835146542);
   graph->SetPoint(14,29,8.470757218);
   graph->SetPoint(15,31,7.884155547);
   graph->SetPoint(16,33,8.699373476);
   graph->SetPoint(17,35,7.037448161);
   graph->SetPoint(18,36.36269779,7);
   graph->SetPoint(19,37,6.935467774);
   graph->SetPoint(20,39,6.332075309);
   graph->SetPoint(21,41,5.189228893);
   graph->SetPoint(22,43,5.5327575);
   graph->SetPoint(23,45,5.255041526);
   graph->SetPoint(24,45.77475707,5);
   graph->Draw("apl");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
