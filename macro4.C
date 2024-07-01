void macro4() {
    // Open the ROOT file and retrieve histograms
    TFile *file = new TFile("output-egamma.root");
    TH1F *vl = (TH1F*)file->Get("d_rec_eff_data_pt_vl");
    TH1F *lo = (TH1F*)file->Get("d_rec_eff_data_pt_lo");
    TH1F *me = (TH1F*)file->Get("d_rec_eff_data_pt_me");
    TH1F *ti = (TH1F*)file->Get("d_rec_eff_data_pt_ti");

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1 -> Divide(2,2);
    // Draw histograms on the same canvas
    c1 -> cd(1);
    vl->SetLineColor(kRed);
    vl->Draw("E");  // E option to draw error bars
    c1 -> cd(2);
    lo->SetLineColor(kBlue);
    lo->Draw("E SAME");  // SAME option to overlay on the same canvas
    c1 -> cd(3);
    me->SetLineColor(kGreen);
    me->Draw("E SAME");
    c1 -> cd(4);
    ti->SetLineColor(kMagenta);
    ti->Draw("E SAME");

    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(vl, "Very Loose", "l");
    legend->AddEntry(lo, "Loose", "l");
    legend->AddEntry(me, "Medium", "l");
    legend->AddEntry(ti, "Tight", "l");

    // Draw the legend
    legend->Draw();

    // Update and save the canvas
    c1->Update();
    c1->SaveAs("PktPracy.jpg");
}