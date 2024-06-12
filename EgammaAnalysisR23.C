#define EgammaAnalysisR23_cxx
#include "EgammaAnalysisR23.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

bool EgammaAnalysisR23 :: ESel (float el_pt, float el_eta)
{
   float el_abseta =  TMath::Abs( el_eta);

   if( el_abseta  > 2.47 ) return false;

   if( el_abseta  > 1.37 and  el_abseta < 1.52) return false;

   if(el_pt < 2.) return false;

   return true;
}

bool EgammaAnalysisR23 :: SelectTruthElectron (float el_eta, float el_pt){

   //eta requirement
   float el_abseta =  TMath::Abs( el_eta);

   if( el_abseta  > 2.47 ) return false;

   if( el_abseta  > 1.37 and  el_abseta < 1.52) return false;

   //pt photon
   //if( el_pt < 1. ) return false;

   return true;

}
/*********************************/
bool EgammaAnalysisR23 :: SelectRecoElectron (float el_eta, float el_pt){

   //eta requirement
   float el_abseta =  TMath::Abs( el_eta);

   if( el_abseta  > 2.47 ) return false;

   if( el_abseta  > 1.37 and  el_abseta < 1.52) return false;

   //pt photon
   if( el_pt < 2.0 ) return false;

   return true;

}
/*********************************/
double EgammaAnalysisR23 :: DeltaPhi(double phi1, double phi2){

  double smaller = std::min(phi1,phi2);
  double bigger  = std::max(phi1,phi2);
  double diff = bigger - smaller;

  if ( diff < TMath::Pi()) return diff;
  diff = TMath::TwoPi() - diff;
  return diff;
}
/*********************************/
double  EgammaAnalysisR23 :: effErf (double x, double p1, double p2) {
      return (0.5*(TMath::Erf((x - p1) / p2) + 1.));
}

double EgammaAnalysisR23 :: L1TriggerWeight(double et1, double et2){
        //parameters for nominal selection
        double p1 = 4.56279;
        double p2 = 2.83293;
        return effErf(et1+et2,p1,p2);
}

double EgammaAnalysisR23 :: L1TriggerWeight2023(double et1, double et2){
        //parameters for nominal selection
        double p1 =  5.8805;
        double p2 =  2.34756;
        double eff = 0.5*(TMath::Erf((et1+et2-p1)/p2)+1); 
        return eff;
}


bool EgammaAnalysisR23 :: FitToTruth (float truth_pt, float truth_eta, float truth_phi, vector<float> *fit_pt, vector<float> *fit_eta, vector<float> *fit_phi, bool probe, bool eg)
{
   double dR_cut = 0.2;
   if (truth_pt>5) dR_cut = 0.04;
   if (eg) dR_cut =0.1;
   if (probe) dR_cut = 0.02;
   int size_fit = fit_pt->size();
   int size_fit1 = fit_eta->size();
   int size_fit2 = fit_phi->size();
   if (size_fit != size_fit1 or size_fit != size_fit2 or size_fit1 != size_fit2) return 0;
   double dR_min = 999.;
   bool var = false;
      for(int j=0; j<size_fit;j++)
      {
         double eta = 0;
         double phi = 0;
         //try
         //{

            eta=fit_eta->at(j);
            phi=fit_phi->at(j);
            
         //}
         //catch(const std::exception& e)
         //{
            //std::cout<<"ERRinF"<< '\n';
         //}
         double del_eta = truth_eta-eta;
         double del_phi = truth_phi-phi;
         double dR = TMath::Sqrt(del_eta*del_eta+del_phi*del_phi);
         
         if(dR<dR_min)
         {
            dR_min = dR;
         }
         

      }
      if(dR_min<dR_cut) var = true;
      return var;
}


void EgammaAnalysisR23::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EgammaAnalysisR23.C
//      root> EgammaAnalysisR23 t
//      root> t.GetEntry(12); // Fill t S members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

//output file
TFile* outputFile = TFile::Open("output-egamma.root", "RECREATE" );

//vectors
Double_t bins_pt[12] = {0., 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 50. };
Double_t bins_mass[14] = {0., 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 50., 100. ,200.};
Double_t bins_eta[21] = {-2.6,-2.32,-2.0,-1.77,-1.52,-1.37,-0.96,-0.68,-0.41,-0.14,0.14,0.41,0.68,0.96,1.37,1.52,1.77,2.0,2.33,2.6};

//histograms
TH1F *h_truth_e_pt = new TH1F("h_truth_e_pt", ";electron p_{T}^{truth} [GeV]; Events", 11, bins_pt);
TH1F *h_match_e_pt = new TH1F("h_match_e_pt", ";electron p_{T}^{truth} [GeV]; Events", 11, bins_pt);
TH1F *h_match_pixtrk_pt = new TH1F("h_match_pixtrk_pt", ";electron p_{T}^{truth} [GeV]; Events", 11, bins_pt);
TH1F *h_lhveryloose_e_pt = new TH1F("h_lhveryloose_e_pt", ";electron p_{T}^{truth} [GeV]; Events", 11, bins_pt);
TH1F *h_truth_e_eta = new TH1F("h_truth_e_eta", ";electron #eta; Events", 100, -2.6, 2.6);
TH1F *h_match_e_eta = new TH1F("h_match_e_eta", ";electron #eta; Events", 100, -2.6, 2.6);
TH1F *h_truth_e_eff_pt = new TH1F("h_truth_e_eff_pt", ";electron p_{T}^{truth} [GeV]; Efficiency", 11, bins_pt);
TH1F *h_lhveryloose_e_eff_pt = new TH1F("h_lhveryloose_e_eff_pt", ";electron p_{T}^{truth} [GeV]; Efficiency", 11, bins_pt);
TH1F *h_truth_e_eff_eta = new TH1F("h_truth_e_eff_eta", ";electron #eta; Efficiency", 100, -2.6, 2.6);
TH1F *h_truth_pixtrk_eff_pt = new TH1F("h_truth_pixtrk_eff_pt", ";electron p_{T}^{truth} [GeV]; Efficiency", 11, bins_pt);

//general
TH2F *h_hlt_fcal_a_c = new TH2F("h_hlt_fcal_a_c", ";FCal E_{T} on side A; FCal E_{T} on side C; Events", 100, -10., 10., 100, -10., 10.);
TH1F *h_trk_eta = new TH1F("h_trk_eta", ";#eta^{track}; Events", 50, -2.5, 2.5);
TH1F *h_trk_phi = new TH1F("h_trk_phi", ";#eta^{track}; Events", 50, -3.2, 3.2);
TH1F *h_trk_n = new TH1F("h_trk_n", ";offline track multiplicity; Events", 20, 0, 20);
TH1F *h_trk_n_sel = new TH1F("h_trk_n_sel", ";offline track multiplicity with p_{T} > 1 GeV; Events", 20, 0, 20);
TH1F *h_trk_n_trig = new TH1F("h_trk_n_trig", ";offline track multiplicity", 20, 0, 20);
TH1F *h_trk_n_sel_trig = new TH1F("h_trk_n_sel_trig", ";offline track multiplicity with p_{T} > 1 GeV; Events", 20, 0, 20);

//yy2ee histos
TH1F *h_aco = new TH1F("h_aco", ";Acoplanarity; Events", 100, 0, 0.1);
TH1F *h_aco_ee = new TH1F("h_aco_ee", ";Acoplanarity; Events", 100, 0, 0.1);
TH1F *h_sumEt_egCluster_ee = new TH1F("h_sumEt_egCluster_ee", ";E_{T,1}^{egCluster}+E_{T,2}^{egCluster} [GeV]; Events", 11, bins_pt);
TH1F *h_sumEt_egCluster_ee_trig = new TH1F("h_sumEt_egCluster_ee_trig", ";E_{T,1}^{egCluster}+E_{T,2}^{egCluster} [GeV]; Events", 11, bins_pt);

//yy2ee histos with electrons
TH1F *h_cutflow_sel_ee = new TH1F("h_cutflow_sel_ee", ";; Events", 10, 0, 10);
TH1F *h_aco_sel_ee = new TH1F("h_aco_sel_ee", ";Acoplanarity; Events", 100, 0, 0.1);
TH1F *h_reco_pt_sel_ee = new TH1F("h_reco_pt_sel_ee", ";electron <p_{T}^{reco}> [GeV]; Events", 11, bins_pt);
TH1F *h_reco_mass_sel_ee = new TH1F("h_reco_mass_sel_ee", "; M_{ee} [GeV]; Events", 13, bins_mass);
TH1F *h_reco_ptee_sel_ee = new TH1F("h_reco_ptee_sel_ee", "; p_{T,ee} [GeV]; Events", 30, 0.,3.);
TH1F *h_reco_costheta_sel_ee = new TH1F("h_reco_costheta_sel_ee", "; |cos(#theta^*)|; Events", 20, 0.,1.);
TH1F *h_reco_yee_sel_ee = new TH1F("h_reco_yee_sel_ee", "; y_{ee}; Events", 100, -3., 3.);
//2D histos for e detection
TH2F *h_truth_e_pt_eta = new TH2F("d_truth_e_pt_eta",";electron p_{T}^{truth} [GeV];electron #eta;Events [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_match_e_pt_eta = new TH2F("d_match_e_pt_eta",";electron p_{T}^{match} [GeV];electron #eta;Events [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_pixtrk_e_pt_eta = new TH2F("d_pixtrk_e_pt_eta",";electron p_{T}^{pixtrack} [GeV];electron #eta;Events [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_eff_pixtrk = new TH2F("d_eff_pixtrk",";electron p_{T}^{match} [GeV];electron #eta;Eff [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_eff_electron = new TH2F("d_eff_electron",";electron p_{T}^{pixtrack} [GeV];electron #eta;Eff [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_trk_e_pt_eta = new TH2F("d_trk_e_pt_eta",";electron p_{T}^{track} [GeV];electron #eta;Events [-]",11,bins_pt,20, -2.6, 2.6);
TH2F *h_eff_trk = new TH2F("d_eff_trk",";electron p_{T}^{track} [GeV];electron #eta;Eff [-]",11,bins_pt,20, -2.6, 2.6);
//1D histos to compare to publication
TH1F *h_egamma_pt = new TH1F("d_egamma_pt", ";electron p_{T}^{egamma} [GeV]; Events", 11, bins_pt);
TH1F *h_egamma_eta = new TH1F("h_egamma_eta", ";electron #eta; Events", 25,  -2.6, 2.6);
TH1F *h_eff_egamma_pt = new TH1F("d_eff_egamma_pt", ";electron p_{T}^{egamma} [GeV]; Eff", 11, bins_pt);
TH1F *h_eff_egamma_eta = new TH1F("d_eff_egamma_etaa", ";electron #eta; Eff", 25,  -2.6, 2.6);
TH1F *h_pixtrk_pt = new TH1F("d_pixtrk_pt", ";electron p_{T}^{pixtrk} [GeV]; Events", 11, bins_pt);
TH1F *h_pixtrk_eta = new TH1F("d_pixtrk_eta", ";electron #eta; Events", 25,  -2.6, 2.6);
TH1F *h_eff_pixtrk_pt = new TH1F("d_eff_pixtrk_pt", ";electron p_{T}^{truth} [GeV]; Eff", 11, bins_pt);
TH1F *h_eff_pixtrk_eta = new TH1F("d_eff_pixtrk_eta", ";electron #eta; Eff", 25,  -2.6, 2.6);
TH1F *h_e_pt  = new TH1F("d_e_pt", ";electron p_{T}^{match} [GeV]; Events", 11, bins_pt);
TH1F *h_e_eta = new TH1F("d_e_eta", ";electron #eta; Events", 25, -2.6,2.6);
TH1F *h_eff_e_pt = new TH1F("d_eff_e_pt", ";electron p_{T}^{match} [GeV]; Eff", 11, bins_pt);
TH1F *h_eff_e_eta = new TH1F("d_eff_e_eta", ";electron #eta; Eff", 25, -2.6,2.6);
TH1F *h_truth_pt = new TH1F("d_truth_e_pt", ";electron p_{T}^{truth} [GeV]; Events", 11, bins_pt);
TH1F *h_truth_eta = new TH1F("d_truth_eta", ";electron #eta; Events", 25, -2.6,2.6);
//histos tagprobe
TH1F *h_tp_mc_aco = new TH1F("d_tp_mc_aco", ";tagprobe Aco; Events", 100, 0., 0.1);
TH1F *h_tp_mc_minv = new TH1F("d_tp_mc_minv", ";tagprobe inv mass; Events", 201,0.,200.);

TH1F *h_tp_data_aco = new TH1F("d_tp_data_aco", ";tagprobe Aco; Events", 100, 0., 0.1);
TH1F *h_tp_data_minv = new TH1F("d_tp_data_minv", ";tagprobe inv mass; Events", 201, 0., 200.);

auto *h_dr_data_probe = new TH1F("d_dr_data_probe", ";probe and e dr; Events", 101, 0., 0.02);

auto *h_dr_mc_probe = new TH1F("d_dr_mc_probe", ";probe and e dr; Events", 101, 0., 4.);

auto *h_rec_data_pt = new TH1F("d_rec_data_pt", ";PixTrk P_T;Events", 11, bins_pt);
auto *h_rec_data_eta = new TH1F("d_rec_data_eta", ";PixTrk Eta;Events", 26, -2.6, 2.6);

auto *h_rec_eff_data_pt = new TH1F("d_rec_eff_data_pt",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_data_eta = new TH1F("d_rec_eff_data_eta",";PixTrk Eta;Eff", 26, -2.6, 2.6);

auto *h_rec_mc_pt = new TH1F("d_rec_mc_pt", ";PixTrk P_T;Events", 11, bins_pt);
auto *h_rec_mc_eta = new TH1F("d_rec_mc_eta", ";PixTrk Eta;Events", 26, -2.6, 2.6);

auto *h_rec_eff_mc_pt = new TH1F("d_rec_eff_mc_pt ",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_eta = new TH1F("d_rec_eff_mc_eta",";PixTrk Eta;Eff", 26, -2.6, 2.6);

auto *h_probe_mc_pt = new TH1F("d_probe_mc_pt ",";PixTrk P_T;Events", 11, bins_pt);
auto *h_probe_mc_eta = new TH1F("d_probe_mc_eta",";PixTrk Eta;Events", 26, -2.6, 2.6);

auto *h_probe_data_pt = new TH1F("d_probe_data_pt",";PixTrk P_T;Events", 11, bins_pt);
auto *h_probe_data_eta = new TH1F("d_probe_data_eta",";PixTrk Eta;Events", 26, -2.6, 2.6);




h_cutflow_sel_ee->GetXaxis()->SetBinLabel(1,"All");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(2,"Trigger");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(3,"2 electrons");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(4,"OS electrons");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(5,"Acoplanarity");



double weight_mc = 1.;


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //check if MC
      bool isMC = false;
      if(mc_channel_number) isMC = true;

   h_cutflow_sel_ee->Fill(0);
   //trigger
   //if(not passed_HLT_mb_sp_vpix30_hi_FgapAC5_L1TAU1_TE4_VTE200) continue;
   //if(not passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L1TAU1_TE4_VTE200) continue;
   //for 1trk5_pt1 studies
   //if(not passed_HLT_mb_sptrk_hi_FgapAC5_L12TAU1_VTE200) continue;
   //if(not passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L1TAU1_TE4_VTE200_EMPTY) continue;
   //for LbyL in R21
   //if(not (passed_HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200  
   //   or passed_HLT_hi_upc_FgapAC3_hi_gg_upc_L12TAU1_VTE50)) continue;
   //for LbyL in R23
   
   if(not isMC and not (passed_HLT_mb_sp_vpix30_hi_FgapAC5_L1TAU1_TE4_VTE200 
      or passed_HLT_mb_sp_vpix30_hi_FgapAC5_L12TAU1_VTE200
      or passed_HLT_mb_sp_vpix30_hi_FgapAC5_2g0_etcut_L12TAU1_VTE200)) continue;
   //if(not  passed_HLT_mb_excl_1trk5_pt2_L1eEM1_VTE200) continue;
   //std::cout << "trigger passed "  << std::endl;

   h_cutflow_sel_ee->Fill(1);
   h_hlt_fcal_a_c->Fill(HLT_FCalAET,HLT_FCalCET);

   //attempt to 1trk5_pt1 efficiency
   h_trk_n->Fill(track_n);

   Int_t track_n_sel_pt1 = 0;
   Int_t track_n_sel_pt0p2 = 0;
   for(Int_t i=0; i<track_n; i++){
      if(track_pt->at(i)>=1) track_n_sel_pt1++;
      if(track_pt->at(i)>=0.2) track_n_sel_pt0p2++;
      h_trk_eta->Fill(track_eta->at(i));
      h_trk_phi->Fill(track_phi->at(i));
   }

   if(track_n_sel_pt0p2<=15) h_trk_n_sel->Fill(track_n_sel_pt1);
   
   if(passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L12TAU1_VTE200){
      h_trk_n_trig->Fill(track_n);
      if(track_n_sel_pt0p2<=15) h_trk_n_sel_trig->Fill(track_n_sel_pt1);
   }
   

   // if(isMC)
   // {
   //    Int_t size_truth_electron = truth_electron_pt->size();
   //    Int_t size_reco_electron = electron_pt->size();
   //    Int_t size_reco_pixtrk = pix_track_pt->size();

   //    for (Int_t i=0; i < size_truth_electron; i++)
   //    {
   //       // if(not SelectTruthElectron(truth_electron_eta->at(i), truth_electron_pt->at(i))) continue;

   //       // h_truth_e_pt->Fill(truth_electron_pt->at(i));


   //       if(ESel(truth_electron_pt->at(i),truth_electron_eta->at(i))) 
   //       {

   //          h_truth_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //          h_truth_pt->Fill(truth_electron_pt->at(i));
   //          h_truth_eta->Fill(truth_electron_eta->at(i));


            
   //          if (FitToTruth(0 ,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),electron_pt,electron_eta,electron_phi)) 
   //          {
   //             h_match_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //             h_e_pt->Fill(truth_electron_pt->at(i));
   //             h_e_eta->Fill(truth_electron_eta->at(i));
   //          }
            
   //          if (FitToTruth(0,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),pix_track_pt,pix_track_eta,pix_track_phi)) 
   //          {
   //             h_pixtrk_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //             h_pixtrk_pt->Fill(truth_electron_pt->at(i));
   //             h_pixtrk_eta->Fill(truth_electron_eta->at(i));

   //          }
   //          if (FitToTruth(0,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),track_pt,track_eta,track_phi)) 
   //          {
   //             h_trk_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //          }
   //          if (FitToTruth(1,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),eg_cluster_et,eg_cluster_eta,eg_cluster_phi))
   //          {
   //             h_egamma_pt->Fill(truth_electron_pt->at(i));
   //             h_egamma_eta->Fill(truth_electron_eta->at(i));
   //          }


   //       }
         
         
         
         // if(truth_electron_pt->at(i) > 2.5) h_truth_e_eta->Fill(truth_electron_eta->at(i));







         // Double_t dR_el_cut = 0.03; //Fig. 106 of https://cds.cern.ch/record/2703499/files/ATL-COM-PHYS-2019-1446.pdf
         // Int_t index_e_match = -1;
         // Int_t index_e_reco_match = -1;
         // Double_t dR_el_min = 999.;
         // for (Int_t j=0; j < size_reco_electron; j++){
         //    Double_t delta_e_eta = truth_electron_eta->at(i)-electron_eta->at(j);
         //    Double_t delta_e_phi = truth_electron_phi->at(i)-electron_phi->at(j);
         //    Double_t dR_el = TMath::Sqrt( delta_e_eta*delta_e_eta + delta_e_phi*delta_e_phi);
            
         //    if(dR_el<dR_el_min) {
         //       dR_el_min = dR_el;
         //        index_e_match = i; //revisit 
         //        index_e_reco_match = j; //revisit 
         //    }
         //}//reco ele

         
         

         // Int_t index_pixtrk_match = -1;
         // Int_t index_pixtrk_reco_match = -1;
         // Double_t dR_pixtrk_min = 999.;
         // for (Int_t j=0; j < size_reco_pixtrk; j++){
         //    Double_t delta_pixtrk_eta = truth_electron_eta->at(i)-pix_track_eta->at(j);
         //    Double_t delta_pixtrk_phi = truth_electron_phi->at(i)-pix_track_phi->at(j);
         //    Double_t dR_pixtrk = TMath::Sqrt( delta_pixtrk_eta*delta_pixtrk_eta + delta_pixtrk_phi*delta_pixtrk_phi);
            
         //    if(dR_pixtrk<dR_pixtrk_min) {
         //       dR_pixtrk_min = dR_pixtrk;
         //        index_pixtrk_match = i; //revisit 
         //        index_pixtrk_reco_match = j; //revisit 
         //    }
         //}//reco pixel tracks

         
      
      // if(dR_el_min<dR_el_cut){
      //    h_match_e_pt->Fill(truth_electron_pt->at(index_e_match));
      //    if(electron_is_LHVeryLoose->at(index_e_reco_match)) h_lhveryloose_e_pt->Fill(truth_electron_pt->at(index_e_match));
      //    if(truth_electron_pt->at(i) > 2.5) 
      //    {
      //       h_match_e_eta->Fill(truth_electron_eta->at(index_e_match));       
      //    }
      // }  

      // if(dR_pixtrk_min<0.2){
      //    h_match_pixtrk_pt->Fill(truth_electron_pt->at(index_pixtrk_match));
      // }  
      
   //}//truth
   
   // h_truth_e_eff_pt->Divide(h_match_e_pt,h_truth_e_pt, 1., 1., "b" );
   // h_lhveryloose_e_eff_pt->Divide(h_lhveryloose_e_pt,h_truth_e_pt, 1., 1., "b" );
   // h_truth_e_eff_eta->Divide(h_match_e_eta,h_truth_e_eta, 1., 1., "b" );
   // h_truth_pixtrk_eff_pt->Divide(h_match_pixtrk_pt,h_truth_e_pt, 1., 1., "b" );
   // h_eff_pixtrk->Divide(h_pixtrk_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   // h_eff_electron->Divide(h_match_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   // h_eff_trk->Divide(h_trk_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   // h_eff_egamma_pt->Divide(h_egamma_pt,h_truth_pt,1.,1.,"b");
   // h_eff_egamma_eta->Divide(h_egamma_eta,h_truth_eta,1.,1.,"b");
   // h_eff_pixtrk_pt->Divide(h_pixtrk_pt,h_truth_pt,1.,1.,"b");
   // h_eff_pixtrk_eta->Divide(h_pixtrk_eta,h_truth_eta,1.,1.,"b");
   // h_eff_e_pt->Divide(h_e_pt,h_truth_pt,1.,1.,"b");
   // h_eff_e_eta->Divide(h_e_eta,h_truth_eta,1.,1.,"b");
   //}//isMC

   int pix_track_n = pix_track_charge -> size();
   if(eg_cluster_n > 0 and eg_cluster_n < 3 and track_n > 0 and track_n < 3 and pix_track_n > 0 and pix_track_n <3)
      {
         for(int i = 0; i<electron_charge->size(); i++)
         {
            if(electron_is_LHMedium->at(i) and electron_pt->at(i) > 2.0)
            {
               for (int j = 0; j<pix_track_charge->size(); j++)
               {
                  if(pix_track_hits->at(j) > 3 and electron_charge->at(i)*pix_track_charge->at(j) < 0 and TMath::Abs(pix_track_eta->at(j)) < 2.5 )
                  {
                     TLorentzVector tag,probe;
                     tag.SetPtEtaPhiM(electron_pt -> at(i), electron_eta -> at(i), electron_phi->at(i),0.511);
                     probe.SetPtEtaPhiM(pix_track_pt -> at(j), pix_track_eta -> at(j), pix_track_phi->at(j), 0.511);
                     Double_t Aco = 1 - (TMath::Abs(tag.Phi() - probe.Phi())/TMath::Pi());
                     double Minv = (tag+probe).M();
                     //Double_t weight = 1.;
                     if(Aco < 0.1 and Minv > 4.0)
                     {
                        Double_t weight = L1TriggerWeight2023(tag.E(),probe.E()); 
                        
                        if(isMC)
                        {
                           h_probe_mc_pt -> Fill(probe.Pt(),weight);
                           h_probe_mc_eta -> Fill(probe.Eta(),weight);
                        }
                        else
                        {
                           h_probe_data_pt -> Fill(probe.Pt());
                           h_probe_data_eta -> Fill(probe.Eta());
                        }
                        if(isMC)
                        {
                           //Double_t weight = L1TriggerWeight2023(tag.E(),probe.E());
                           h_tp_mc_aco -> Fill(Aco, weight);
                           h_tp_mc_minv -> Fill(Minv, weight);
                        }
                        else
                        {
                        h_tp_data_aco -> Fill(Aco);
                        h_tp_data_minv -> Fill(Minv);
                        }
                        double dR_min = 999.;
                        int ind = -1;
                        for(int k = 0; k<electron_charge->size();k++)
                        {
                           TLorentzVector rec;
                           rec.SetPtEtaPhiM(electron_pt->at(k), electron_eta->at(k), electron_phi->at(k), 0.511);
                           double dR = sqrt((probe.Phi()-rec.Phi())*(probe.Phi()-rec.Phi())+(probe.Eta()-rec.Eta())*(probe.Eta()-rec.Eta()));
                           if(isMC) h_dr_mc_probe -> Fill(dR, weight);
                           else h_dr_data_probe -> Fill(dR);
                           if(dR_min > dR) 
                           {
                              dR_min = dR;
                              ind = k;
                           }
                        }
                        if(dR_min < 0.002)
                        {
                           if(isMC)
                           {
                              TLorentzVector rec;
                              rec.SetPtEtaPhiM(electron_pt->at(ind), electron_eta->at(ind), electron_phi->at(ind), 0.511);
                              Double_t weight1 = L1TriggerWeight2023(probe.E(), rec.E());
                              h_rec_mc_pt -> Fill(probe.Pt(), weight1);
                              h_rec_mc_eta -> Fill(probe.Eta(), weight1);
                           }
                           else
                           {
                              h_rec_data_pt -> Fill(probe.Pt());
                              h_rec_data_eta -> Fill(probe.Eta());
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      
   // //this part does not work in R21
   // //yy2ee selection with pixel tracks
   // Int_t pix_trk_size = pix_track_pt->size();
   // Int_t egCluster_size = eg_cluster_et->size();

   // if(pix_trk_size==2 and egCluster_size>=2){
      
   //    if(pix_track_charge->at(0)*pix_track_charge->at(1)>0) continue;
   //    if(pix_track_pt->at(0)< 0.5 or pix_track_pt->at(1)<0.5) continue;

   //    double Aco = 1 - TMath::Abs(DeltaPhi(pix_track_phi->at(0), pix_track_phi->at(1)))/TMath::Pi();
   //    h_aco->Fill(Aco);
   //    if(Aco>0.04) continue;

   //    //matching pix track - egCluster
   //    Int_t index_pixeg_match = -1;
   //    Int_t index_pixeg_reco_match = -1;
   //    Double_t dR_pixeg_min = 999.;
   //    Int_t index_pixeg_match_1 = -1;
   //    Int_t index_pixeg_match_2 = -1;
   //    for(Int_t i = 0; i < pix_trk_size; i++){ 
         
   //       dR_pixeg_min = 999; //needs to be reset

   //       for(Int_t j = 0; j < egCluster_size; j++){
   //          Double_t delta_pixeg_eta = eg_cluster_eta->at(j)-pix_track_eta->at(i);
   //          Double_t delta_pixeg_phi = DeltaPhi(eg_cluster_phi->at(j),pix_track_phi->at(i));
   //          Double_t dR_pixeg = TMath::Sqrt( delta_pixeg_eta*delta_pixeg_eta + delta_pixeg_phi*delta_pixeg_phi);
            
   //          if(dR_pixeg<dR_pixeg_min) {
   //             dR_pixeg_min = dR_pixeg;
   //              index_pixeg_match = j; //revisit 
   //              index_pixeg_reco_match = i; //revisit 
   //          }
   //       }//eg cluster
        
   //       if(dR_pixeg_min<0.8){
   //          if(i == 0) index_pixeg_match_1 = index_pixeg_match;
   //          if(i == 1) index_pixeg_match_2 = index_pixeg_match;
   //       } 
   //    } //pix tracks

   // //std::cout << "index_pixeg_match_1 = " << index_pixeg_match_1 << std::endl;
   // //std::cout << "index_pixeg_match_2 = " << index_pixeg_match_2 << std::endl;

   // if(index_pixeg_match_1 < 0 or index_pixeg_match_2 < 0) continue;

   // double sum_et_eg = eg_cluster_et->at(index_pixeg_match_1)+eg_cluster_et->at(index_pixeg_match_2);
   // //if(isMC) weight_mc=L1TriggerWeight(eg_cluster_et->at(index_pixeg_match_1),eg_cluster_et->at(index_pixeg_match_2));
   // if(isMC) weight_mc=L1TriggerWeight2023(eg_cluster_et->at(index_pixeg_match_1),eg_cluster_et->at(index_pixeg_match_2));
   // //weight_mc = 1.;
   // //std::cout << "weight = " << weight_mc << std::endl;

   // h_sumEt_egCluster_ee->Fill(sum_et_eg,weight_mc);
   // h_aco_ee->Fill(Aco,weight_mc);

   // if(passed_L1_eEM1_VTE200_TBP) h_sumEt_egCluster_ee_trig->Fill(sum_et_eg,weight_mc);
   
   // }//yy2ee
   

   // //yy2ee selection with electrons
   // Int_t electron_size=electron_pt->size();
   // Int_t good_electrons = 0;
   // std::vector<Int_t>  good_electron_index;
   // for(Int_t i=0; i<electron_size; i++){
   //    if(SelectRecoElectron(electron_eta->at(i), electron_pt->at(i))){
   //      good_electrons++;
   //      good_electron_index.push_back(i);
   //    } 
   // }//yy2ee with electrons

   // if(good_electrons!=2) continue;

   // h_cutflow_sel_ee->Fill(2);

   // Int_t index_e1 = good_electron_index.at(0);
   // Int_t index_e2 = good_electron_index.at(1);
   // double dphi_ee = DeltaPhi(electron_phi->at(index_e1), electron_phi->at(index_e2));
   // double aco_ee=1 - TMath::Abs(dphi_ee)/TMath::Pi(); 
   // TLorentzVector v_e1;
   // TLorentzVector v_e2;
   // v_e1.SetPtEtaPhiM(electron_pt->at(index_e1),electron_eta->at(index_e1),electron_phi->at(index_e1),0.511);
   // v_e2.SetPtEtaPhiM(electron_pt->at(index_e2),electron_eta->at(index_e2),electron_phi->at(index_e2),0.511);
   // double mass_ee = (v_e1+v_e2).M();
   // double rapidity_ee = (v_e1+v_e2).Rapidity();
   // double pt_ee = (v_e1+v_e2).Pt();
   // double costheta_ee = TMath::Abs(tanh((electron_eta->at(index_e1)-electron_eta->at(index_e2))/2.));

   // if(electron_charge->at(index_e1)*electron_charge->at(index_e2)>0) continue;
   // h_cutflow_sel_ee->Fill(3);

   // if(aco_ee>0.02) continue;
   // h_cutflow_sel_ee->Fill(4);
   
   // h_aco_sel_ee->Fill(aco_ee,weight_mc);

   // double av_pt = (electron_pt->at(index_e1)+electron_pt->at(index_e2))/2;
   // h_reco_pt_sel_ee->Fill(av_pt,weight_mc);
   // h_reco_mass_sel_ee->Fill(mass_ee,weight_mc);
   // h_reco_yee_sel_ee->Fill(rapidity_ee,weight_mc);
   // h_reco_ptee_sel_ee->Fill(pt_ee,weight_mc);
   // h_reco_costheta_sel_ee->Fill(costheta_ee,weight_mc);

   }//event loop
   // Save the output
   h_rec_eff_data_pt -> Divide(h_rec_data_pt, h_probe_data_pt, 1. ,1., "b");
   h_rec_eff_data_eta -> Divide(h_rec_data_eta, h_probe_data_eta, 1. ,1., "b");
   h_rec_eff_mc_pt -> Divide(h_rec_mc_pt, h_probe_mc_pt, 1. ,1., "b");
   h_rec_eff_mc_eta -> Divide(h_rec_mc_eta, h_probe_mc_eta, 1. ,1., "b");
   outputFile->Write();
   outputFile->Close();
}

//TEffitiency