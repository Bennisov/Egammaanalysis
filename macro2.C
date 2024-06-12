void macro2() {
    TFile *f = new TFile("output-egamma.root");
    TH1F *mc_eta_p = (TH1F*)f->Get("d_probe_mc_eta");
    TH1F *mc_pt_p = (TH1F*)f->Get("d_probe_mc_pt");
    TH1F *data_eta_p = (TH1F*)f->Get("d_probe_data_eta");
    TH1F *data_pt_p = (TH1F*)f->Get("d_probe_data_pt");

    TH1F *mc_eta_r = (TH1F*)f->Get("d_rec_mc_eta");
    TH1F *mc_pt_r = (TH1F*)f->Get("d_rec_mc_pt");
    TH1F *data_eta_r = (TH1F*)f->Get("d_rec_data_eta");
    TH1F *data_pt_r = (TH1F*)f->Get("d_rec_data_pt");

    if (data_pt_r->GetEntries() == 0 || mc_pt_r->GetEntries() == 0 || data_eta_r->GetEntries() == 0 || mc_eta_r->GetEntries() == 0) {
        std::cerr << "Error: Histograms have no entries." << std::endl;
        return;
    }

    auto calculateEfficiency = [](TH1F* num, TH1F* den) {
        int nBins = num->GetNbinsX();
        std::vector<double> eff(nBins), eff_err(nBins);
        for (int i = 1; i <= nBins; ++i) {
            double n = num->GetBinContent(i);
            double d = den->GetBinContent(i);
            if (d == 0) {
                eff[i-1] = 0;
                eff_err[i-1] = 0;
            } else {
                eff[i-1] = n / d;
                eff_err[i-1] = sqrt(eff[i-1] * (1 - eff[i-1]) / d);
            }
        }
        return std::make_pair(eff, eff_err);
    };

    auto mc_eta_eff = calculateEfficiency(mc_eta_r, mc_eta_p);
    auto mc_pt_eff = calculateEfficiency(mc_pt_r, mc_pt_p);
    auto data_eta_eff = calculateEfficiency(data_eta_r, data_eta_p);
    auto data_pt_eff = calculateEfficiency(data_pt_r, data_pt_p);

    auto createTGraphErrors = [](TH1F* hist, const std::vector<double>& eff, const std::vector<double>& eff_err) {
        int nBins = hist->GetNbinsX();
        TGraphErrors* graph = new TGraphErrors(nBins);
        for (int i = 1; i <= nBins; ++i) {
            double x = hist->GetBinCenter(i);
            double ex = hist->GetBinWidth(i) / 2.0;
            graph->SetPoint(i-1, x, eff[i-1]);
            graph->SetPointError(i-1, ex, eff_err[i-1]);
        }
        return graph;
    };

    TGraphErrors* gr_mc_eta_eff = createTGraphErrors(mc_eta_r, mc_eta_eff.first, mc_eta_eff.second);
    TGraphErrors* gr_mc_pt_eff = createTGraphErrors(mc_pt_r, mc_pt_eff.first, mc_pt_eff.second);
    TGraphErrors* gr_data_eta_eff = createTGraphErrors(data_eta_r, data_eta_eff.first, data_eta_eff.second);
    TGraphErrors* gr_data_pt_eff = createTGraphErrors(data_pt_r, data_pt_eff.first, data_pt_eff.second);

    auto calculateRatio = [](const std::vector<double>& num_eff, const std::vector<double>& num_err, const std::vector<double>& den_eff, const std::vector<double>& den_err) {
        int nBins = num_eff.size();
        std::vector<double> ratio(nBins), ratio_err(nBins);
        for (int i = 0; i < nBins; ++i) {
            if (den_eff[i] == 0) {
                ratio[i] = 0;
                ratio_err[i] = 0;
            } else {
                ratio[i] = num_eff[i] / den_eff[i];
                ratio_err[i] = ratio[i] * sqrt(pow(num_err[i] / num_eff[i], 2) + pow(den_err[i] / den_eff[i], 2));
            }
        }
        return std::make_pair(ratio, ratio_err);
    };

    auto eta_ratio = calculateRatio(data_eta_eff.first, data_eta_eff.second, mc_eta_eff.first, mc_eta_eff.second);
    auto pt_ratio = calculateRatio(data_pt_eff.first, data_pt_eff.second, mc_pt_eff.first, mc_pt_eff.second);

    TGraphErrors* gr_eta_ratio = createTGraphErrors(mc_eta_r, eta_ratio.first, eta_ratio.second);
    TGraphErrors* gr_pt_ratio = createTGraphErrors(mc_pt_r, pt_ratio.first, pt_ratio.second);

    // Create canvases for plotting
    TCanvas *canvas1 = new TCanvas("c1", "", 800, 800);
    canvas1->Divide(1, 2);

    // Plot eta efficiency
    canvas1->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    gr_mc_eta_eff->SetMarkerStyle(20);
    gr_mc_eta_eff->SetMarkerColor(kBlue);
    gr_mc_eta_eff->SetLineColor(kBlue);
    gr_mc_eta_eff->Draw("AP");
    gr_data_eta_eff->SetMarkerStyle(21);
    gr_data_eta_eff->SetMarkerColor(kRed);
    gr_data_eta_eff->SetLineColor(kRed);
    gr_data_eta_eff->Draw("P SAME");

    TLegend *legend_eta = new TLegend();
    legend_eta->AddEntry(gr_mc_eta_eff, "MC", "lep");
    legend_eta->AddEntry(gr_data_eta_eff, "Data", "lep");
    legend_eta->Draw();

    // Plot eta ratio
    canvas1->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.4);

    gr_eta_ratio->SetMarkerStyle(20);
    gr_eta_ratio->SetMarkerColor(kBlack);
    gr_eta_ratio->SetLineColor(kBlack);
    gr_eta_ratio->Draw("AP");
    gr_eta_ratio->GetYaxis()->SetTitle("DATA / MC");
    gr_eta_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    gr_eta_ratio->GetXaxis()->SetTitle("#eta");

    canvas1->Update();
    canvas1->SaveAs("eff_eta.png");

    TCanvas *canvas2 = new TCanvas("c2", "", 800, 800);
    canvas2->Divide(1, 2);

    // Plot pt efficiency
    canvas2->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    gr_mc_pt_eff->SetMarkerStyle(20);
    gr_mc_pt_eff->SetMarkerColor(kBlue);
    gr_mc_pt_eff->SetLineColor(kBlue);
    gr_mc_pt_eff->Draw("AP");
    gr_data_pt_eff->SetMarkerStyle(21);
    gr_data_pt_eff->SetMarkerColor(kRed);
    gr_data_pt_eff->SetLineColor(kRed);
    gr_data_pt_eff->Draw("P SAME");

    TLegend *legend_pt = new TLegend();
    legend_pt->AddEntry(gr_mc_pt_eff, "MC", "lep");
    legend_pt->AddEntry(gr_data_pt_eff, "Data", "lep");
    legend_pt->Draw();

    // Plot pt ratio
    canvas2->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.4);

    gr_pt_ratio->SetMarkerStyle(20);
    gr_pt_ratio->SetMarkerColor(kBlack);
    gr_pt_ratio->SetLineColor(kBlack);
    gr_pt_ratio->Draw("AP");
    gr_pt_ratio->GetYaxis()->SetTitle("DATA / MC");
    gr_pt_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    gr_pt_ratio->GetXaxis()->SetTitle("p_{T} [GeV]");

    canvas2->Update();
    canvas2->SaveAs("eff_pt.png");
}