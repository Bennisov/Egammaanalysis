void start()
{
    double crossection = 3714 * 1000; //nb
    double Ngen = 99500;
    double Ldat = 1.7; //nb-1

    double w2 = Ldat / (Ngen / crossection);

    TFile *f = new TFile("output-egamma.root");
    TH1F *mca = (TH1F*)f->Get("d_tp_mc_aco");
    TH1F *mc = (TH1F*)f->Get("d_tp_mc_minv");
    TH1F *dataa = (TH1F*)f->Get("d_tp_data_aco");
    TH1F *data = (TH1F*)f->Get("d_tp_data_minv");

    // Check if histograms have entries
    if (data->GetEntries() == 0 || mc->GetEntries() == 0) {
        std::cerr << "Error: Histograms have no entries." << std::endl;
        return;
    }

    mc->Scale(w2);
    data->Sumw2();

    TCanvas *canvas = new TCanvas("c1", "", 800, 800);
    canvas->Divide(1, 2); // Divide the canvas into two pads vertically

    // Upper pad for main histogram plot
    canvas->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0); // Set position and size
    gPad->SetLogy(1);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
    data -> GetXaxis() -> SetRangeUser(0,20);
    data->Draw("E1");
    mc->SetLineColor(kRed);
    mc->SetMarkerStyle(0);
    mc->Draw("SAME L");
    TLegend* legend = new TLegend();
    legend->AddEntry(data, "Data", "lep");
    legend->AddEntry(mc, "MC", "lep");
    legend->Draw();

    // Lower pad for ratio plot
    canvas->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3); // Set position and size
    gPad->SetTopMargin(0); // No margin at the top
    gPad->SetBottomMargin(0.4); // Adjust bottom margin for ratio plot
    TH1F *ratioHist = (TH1F*)data->Clone("ratioHist");
    ratioHist->Divide(mc);
    ratioHist->SetStats(0);
    ratioHist->SetTitle("");
    ratioHist->SetMarkerStyle(20);
    ratioHist->SetMarkerSize(0.8);
    ratioHist->SetMarkerColor(kBlack);
    ratioHist->SetLineColor(kBlack);
    ratioHist->GetYaxis()->SetTitle("Data / MC");
    ratioHist->GetYaxis()->SetRangeUser(0.5, 1.5); // Adjust y-axis range to show ratios around 1
    ratioHist->GetXaxis()->SetTitle("M_{inv} of tag and probe [GeV]");
    ratioHist->GetXaxis()->SetRangeUser(0, 20);
    ratioHist->Draw("E1");

    canvas->Update();
    canvas->SaveAs("tagprobeminv.png");

    // Check if histograms have entries
    if (dataa->GetEntries() == 0 || mca->GetEntries() == 0) {
        std::cerr << "Error: Histograms have no entries." << std::endl;
        return;
    }

    mca->Scale(w2);
    dataa->Sumw2();

    TCanvas *canvas1 = new TCanvas("c2", "", 800, 800);
    canvas1->Divide(1, 2); // Divide the canvas into two pads vertically

    // Upper pad for main histogram plot
    canvas1->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0); // Set position and size
    gPad->SetLogy(1);
    dataa->SetMarkerStyle(20);
    dataa->SetMarkerSize(0.8);
    dataa -> GetXaxis() -> SetRangeUser(0,0.01);
    dataa->Draw("E1");
    mca->SetLineColor(kRed);
    mca->SetMarkerStyle(0);
    mca->Draw("SAME L");
    TLegend* legend1 = new TLegend();
    legend1->AddEntry(dataa, "data", "lep");
    legend1->AddEntry(mca, "MC", "lep");
    legend1->Draw();

    // Lower pad for ratio plot
    canvas1->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3); // Set position and size
    gPad->SetTopMargin(0); // No margin at the top
    gPad->SetBottomMargin(0.4); // Adjust bottom margin for ratio plot
    TH1F *ratioHist1 = (TH1F*)dataa->Clone("ratioHist");
    ratioHist1->Divide(mca);
    ratioHist1->SetStats(0);
    ratioHist1->SetTitle("");
    ratioHist1->SetMarkerStyle(20);
    ratioHist1->SetMarkerSize(0.8);
    ratioHist1->SetMarkerColor(kBlack);
    ratioHist1->SetLineColor(kBlack);
    ratioHist1->GetYaxis()->SetTitle("Data / MC");
    ratioHist1->GetYaxis()->SetRangeUser(0., 2.); // Adjust y-axis range to show ratios around 1
    ratioHist1->GetXaxis()->SetTitle("Aco tagprobe");
    ratioHist->GetXaxis()->SetRangeUser(0,0.01);
    ratioHist1->Draw("E1");

    canvas1->Update();
    canvas1 -> SaveAs("tagprobeaco.png");
}