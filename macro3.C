void macro3()
{
    TFile *f = new TFile("output-egamma.root");
    TH1F *mc_eta = (TH1F*)f->Get("d_rec_eff_mc_eta");
    TH1F *mc_pt = (TH1F*)f->Get("d_rec_eff_mc_pt");
    TH1F *data_eta = (TH1F*)f->Get("d_rec_eff_data_eta");
    TH1F *data_pt = (TH1F*)f->Get("d_rec_eff_data_pt");

    if (data_pt->GetEntries() == 0 || mc_pt->GetEntries() == 0) {
        std::cerr << "Error: Histograms have no entries." << std::endl;
        return;
    }

    TCanvas *canvas = new TCanvas("c1", "", 800, 800);
    canvas->Divide(1, 2);

    canvas->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    data_pt->SetMarkerStyle(20);
    data_pt->SetMarkerSize(0.8);
    data_pt->SetMarkerColor(kBlack);
    data_pt->SetLineColor(kBlack);
    data_pt->GetYaxis()->SetRangeUser(0.,1.);
    //data_pt->GetXaxis()->SetRangeUser(0, 20);
    //data_pt->Draw("E1");
    mc_pt->SetMarkerStyle(20);
    mc_pt->SetMarkerSize(0.8);
    mc_pt->SetMarkerColor(kBlue);
    mc_pt->SetLineColor(kBlue);
    mc_pt->Draw("SAME E1");
    TLegend* legend = new TLegend();
    legend->AddEntry(data_pt, "Data", "lep");
    legend->AddEntry(mc_pt, "MC", "lep");
    legend->Draw();

    canvas->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3); 
    gPad->SetTopMargin(0); 
    gPad->SetBottomMargin(0.4); 
    // TH1F *ratioHist = (TH1F*)data_pt->Clone("ratioHist");
    // ratioHist->Divide(mc_pt);
    // ratioHist->SetStats(0);
    // ratioHist->SetTitle("");
    // ratioHist->SetMarkerStyle(20);
    // ratioHist->SetMarkerSize(0.8);
    // ratioHist->SetMarkerColor(kBlack);
    // ratioHist->SetLineColor(kBlack);
    // ratioHist->GetYaxis()->SetTitle("DATA / MC");
    // ratioHist->GetYaxis()->SetRangeUser(0.5, 1.5);
    // ratioHist->GetXaxis()->SetTitle("pT of pixtrk [GeV]");
    // //ratioHist->GetXaxis()->SetRangeUser(0, 20);
    // ratioHist->Draw("E1");

    canvas->Update();
    canvas->SaveAs("eff_pt.png");

    if (data_eta->GetEntries() == 0 || mc_eta->GetEntries() == 0) {
        std::cerr << "Error: Histograms have no entries." << std::endl;
        return;
    }

    TCanvas *canvas1 = new TCanvas("c2", "", 800, 800);
    canvas1->Divide(1, 2); 

    canvas1->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0); 
    data_eta->SetMarkerStyle(20);
    data_eta->SetMarkerSize(0.8);
    data_eta->SetMarkerColor(kBlack);
    data_eta->SetLineColor(kBlack);
    data_eta->GetYaxis()->SetRangeUser(0.,1.);
    //data_eta->GetXaxis()->SetRangeUser(0, 0.01);
    //data_eta->Draw("E1");
    mc_eta->SetMarkerStyle(20);
    mc_eta->SetMarkerSize(0.8);
    mc_eta->SetMarkerColor(kBlue);
    mc_eta->SetLineColor(kBlue);
    mc_eta->Draw("SAME E1");
    TLegend* legend1 = new TLegend();
    legend1->AddEntry(data_eta, "Data", "lep");
    legend1->AddEntry(mc_eta, "MC", "lep");
    legend1->Draw();


    canvas1->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3); 
    gPad->SetTopMargin(0); 
    gPad->SetBottomMargin(0.4); 
    // TH1F *ratioHist1 = (TH1F*)data_eta->Clone("ratioHist1");
    // ratioHist1->Divide(mc_eta);
    // ratioHist1->SetStats(0);
    // ratioHist1->SetTitle("");
    // ratioHist1->SetMarkerStyle(20);
    // ratioHist1->SetMarkerSize(0.8);
    // ratioHist1->SetMarkerColor(kBlack);
    // ratioHist1->SetLineColor(kBlack);
    // ratioHist1->GetYaxis()->SetTitle("DATA / MC");
    // ratioHist1->GetYaxis()->SetRangeUser(0.5, 1.5); 
    // ratioHist1->GetXaxis()->SetTitle("Aco tagprobe");
    // //ratioHist1->GetXaxis()->SetRangeUser(0, 0.01);
    // ratioHist1->Draw("E1");

    canvas1->Update();
    canvas1->SaveAs("eff_eta.png");
}