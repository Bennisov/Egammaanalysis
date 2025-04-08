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
double EgammaAnalysisR23 :: effErf(double x, double p1, double p2, double p3) {
    double arg = (x - p1) / (p2 + p3 * x);
    return 0.5 * (TMath::Erf(arg) + 1.0);
}
double EgammaAnalysisR23 :: L1TriggerWeight(double et1, double et2){
	//parameters for nominal selection
	double p1 = 5.63188; 
	double p2 = 1.63362;
	double p3 = 0.223143;
	double sumEt = et1+et2;

	return effErf(sumEt, p1, p2, p3);
}
double EgammaAnalysisR23 :: L1TriggerWeightUp(double et1, double et2){
	//parameters and their errors for nominal selection
	double p1 = 5.48885;
	double p2 = 1.67934;
	double p3 = 0.193579;
	double sumEt = et1+et2;
	
	double uncertUp = effErf(sumEt, p1, p2, p3);
	return uncertUp;
	
}
double EgammaAnalysisR23 :: L1TriggerWeightDown(double et1, double et2){
	//parameters and their errors for nominal selection
	double p1 = 5.80343;
	double p2 = 1.70899;
	double p3 = 0.219453;
	double sumEt = et1+et2;

	double uncertDown = effErf(sumEt, p1, p2, p3);
	return uncertDown;
}
bool EgammaAnalysisR23 :: FitToTruth(bool eg, float truth_pt, float truth_eta, float truth_phi, vector<float> *fit_pt, vector<float> *fit_eta, vector<float> *fit_phi)
{
   double dR_cut = 0.1;
   if (truth_pt>5) dR_cut = 0.05;
   if (eg) dR_cut = 0.4;
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
TFile* outputFile = TFile::Open("output-basev2.root", "RECREATE" );

//vectors
    Double_t bins_pt[12] = {0., 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 50. };
    Double_t bins_mass[14] = {0., 1., 2., 3., 4., 5., 7., 10., 15., 20., 25., 50., 100. ,200.};
    Double_t bins_eta[11] = {-2.47, -1.52, -1.37 , -1.,  -0.5, 0., 0.5, 1., 1.37, 1.52, 2.47}; //bin_1D
    Double_t bins_pt2[7] = {0., 2., 3., 5., 7., 10. ,50.}; //bin_2D


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

auto *h_dr_data_probe = new TH1F("d_dr_data_probe", ";probe and e dr; Events", 101, 0., 0.1);

auto *h_dr_mc_probe = new TH1F("d_dr_mc_probe", ";probe and e dr; Events", 101, 0., 0.1);

auto *h_rec_data_pt = new TH1F("d_rec_data_pt", ";PixTrk P_T;Events", 11, bins_pt);
auto *h_rec_data_eta = new TH1F("d_rec_data_eta", ";PixTrk Eta;Events", 10, bins_eta);

auto *h_rec_eff_data_pt = new TH1F("d_rec_eff_data_pt",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_data_pt_vl = new TH1F("d_rec_eff_data_pt_vl",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_data_pt_lo = new TH1F("d_rec_eff_data_pt_lo",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_data_pt_me = new TH1F("d_rec_eff_data_pt_me",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_data_pt_ti = new TH1F("d_rec_eff_data_pt_ti",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_pt_vl = new TH1F("d_rec_eff_mc_pt_vl",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_pt_lo = new TH1F("d_rec_eff_mc_pt_lo",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_pt_me = new TH1F("d_rec_eff_mc_pt_me",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_pt_ti = new TH1F("d_rec_eff_mc_pt_ti",";PixTrk P_T;Eff", 11, bins_pt);

auto *h_rec_eff_data_eta_vl = new TH1F("d_rec_eff_data_eta_vl",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_data_eta_lo = new TH1F("d_rec_eff_data_eta_lo",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_data_eta_me = new TH1F("d_rec_eff_data_eta_me",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_data_eta_ti = new TH1F("d_rec_eff_data_eta_ti",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_mc_eta_vl = new TH1F("d_rec_eff_mc_eta_vl",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_mc_eta_lo = new TH1F("d_rec_eff_mc_eta_lo",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_mc_eta_me = new TH1F("d_rec_eff_mc_eta_me",";PixTrk P_T;Eff", 10, bins_eta);
auto *h_rec_eff_mc_eta_ti = new TH1F("d_rec_eff_mc_eta_ti",";PixTrk P_T;Eff", 10, bins_eta);

auto *h_rec_op_data_pt_vl = new TH1F("d_rec_op_data_pt_vl",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_data_pt_lo = new TH1F("d_rec_op_data_pt_lo",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_data_pt_me = new TH1F("d_rec_op_data_pt_me",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_data_pt_ti = new TH1F("d_rec_op_data_pt_ti",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_mc_pt_vl = new TH1F("d_rec_op_mc_pt_vl",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_mc_pt_lo = new TH1F("d_rec_op_mc_pt_lo",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_mc_pt_me = new TH1F("d_rec_op_mc_pt_me",";PixTrk P_T;op", 11, bins_pt);
auto *h_rec_op_mc_pt_ti = new TH1F("d_rec_op_mc_pt_ti",";PixTrk P_T;op", 11, bins_pt);

auto *h_rec_op_data_eta_vl = new TH1F("d_rec_op_data_eta_vl",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_data_eta_lo = new TH1F("d_rec_op_data_eta_lo",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_data_eta_me = new TH1F("d_rec_op_data_eta_me",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_data_eta_ti = new TH1F("d_rec_op_data_eta_ti",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_mc_eta_vl = new TH1F("d_rec_op_mc_eta_vl",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_mc_eta_lo = new TH1F("d_rec_op_mc_eta_lo",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_mc_eta_me = new TH1F("d_rec_op_mc_eta_me",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_op_mc_eta_ti = new TH1F("d_rec_op_mc_eta_ti",";PixTrk P_T;op", 10, bins_eta);
auto *h_rec_eff_data_eta = new TH1F("d_rec_op_data_eta",";PixTrk Eta;op", 10, bins_eta);

auto *h_rec_mc_pt = new TH1F("d_rec_mc_pt", ";PixTrk P_T;Events", 11, bins_pt);
auto *h_rec_mc_eta = new TH1F("d_rec_mc_eta", ";PixTrk Eta;Events", 10, bins_eta);

auto *h_rec_eff_mc_pt = new TH1F("d_rec_eff_mc_pt",";PixTrk P_T;Eff", 11, bins_pt);
auto *h_rec_eff_mc_eta = new TH1F("d_rec_eff_mc_eta",";PixTrk Eta;Eff", 10, bins_eta);

auto *h_probe_mc_pt = new TH1F("d_probe_mc_pt ",";PixTrk P_T;Events", 11, bins_pt);
auto *h_probe_mc_eta = new TH1F("d_probe_mc_eta",";PixTrk Eta;Events", 10, bins_eta);

auto *h_probe_data_pt = new TH1F("d_probe_data_pt",";PixTrk P_T;Events", 11, bins_pt);
auto *h_probe_data_eta = new TH1F("d_probe_data_eta",";PixTrk Eta;Events", 10, bins_eta);

auto *h_ratio_pt = new TH1F("d_ratio_pt", "", 11, bins_pt);
auto *h_ratio_eta = new TH1F("d_ratio_eta", "", 10, bins_eta);

auto *h_cut_data_aco = new TH1F("d_cut_data_aco", "", 100, 0., 0.1);
auto *h_tclu_mult = new TH1F("d_tclu_mult", "", 10, 0, 10);
auto *h_tclu_mult_tp_data = new TH1F("d_tclu_mult_tp_data", "", 10, 0, 10);
auto *h_tclu_mult_tp_mc = new TH1F("d_tclu_mult_tp_mc", "", 10, 0, 10);
auto *h_tclu_pt_data = new TH1F("d_tclu_pt_data", "", 100, 0, 10);
auto *h_tclu_pt_mc = new TH1F("d_tclu_pt_mc", "", 100, 0, 10);

auto *h_rec_data_2d = new TH2F("d_rec_data_2d", "", 6, bins_pt2, 10, bins_eta);
auto *h_rec_mc_2d = new TH2F("d_rec_mc_2d", "", 6, bins_pt2, 10, bins_eta);

auto *h_probe_data_2d = new TH2F("d_probe_data_2d", "", 6, bins_pt2, 10, bins_eta);
auto *h_probe_mc_2d = new TH2F("d_probe_mc_2d", "", 6, bins_pt2, 10, bins_eta);

auto *h_rec_eff_mc_2d = new TH2F("d_rec_eff_mc_2d", "", 6, bins_pt2, 10, bins_eta);
auto *h_rec_eff_data_2d = new TH2F("d_rec_eff_data_2d", "", 6, bins_pt2, 10, bins_eta);
auto *h_ratio_2d = new TH2F("d_ratio_2d", "", 6, bins_pt2, 10, bins_eta);

auto *h_dr_l5_pixtrk = new TH1F("d_dr_l5_pixtrk", "", 100, 0., 0.5);
auto *h_dr_g5_pixtrk = new TH1F("d_dr_g5_pixtrk", "", 100, 0., 0.5);
auto *h_dr_egamma = new TH1F("d_dr_egamma", "", 100, 0., 0.5);
auto *h_dr_l5_e = new TH1F("d_dr_l5_e", "", 100, 0., 0.5);
auto *h_dr_g5_e = new TH1F("d_dr_g5_e", "", 100, 0., 0.5);
auto *h_dr_e_pixtrk_l3 = new TH1F("d_dr_e_pixtrk_l3", "", 100, 0., 0.2);
auto *h_dr_e_pixtrk_g3 = new TH1F("d_dr_e_pixtrk_g3", "", 100, 0., 0.2);
auto *h_tp_mc_pt = new TH1F("d_tp_mc_pt", "", 100, 0, 10);
auto *h_tp_data_pt = new TH1F("d_tp_data_pt", "", 100, 0, 10);

auto zdc = new TH2F("d_zdc", "", 400, 0, 10, 400, 0, 10);

auto probe_d0_mc = new TH1F("d_probe_d0_mc", "", 100, -20, 20);
auto probe_z0_mc = new TH1F("d_probe_z0_mc", "", 100, -200, 1200);
auto probe_d0_data = new TH1F("d_probe_d0_data", "", 100, -20, 20);
auto probe_z0_data = new TH1F("d_probe_z0_data", "", 100, -200, 1200);



h_cutflow_sel_ee->GetXaxis()->SetBinLabel(1,"All");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(2,"Trigger");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(3,"2 electrons");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(4,"OS electrons");
h_cutflow_sel_ee->GetXaxis()->SetBinLabel(5,"Acoplanarity");



double weight_mc = 1.;
int tag_probe_counter = 0;

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
   
   for (size_t i = 0; i < pix_track_pt->size(); i++)
   {
      for (size_t j = 0; j < electron_pt->size(); j++)
      {
         double dphi = DeltaPhi(pix_track_phi->at(i), electron_phi->at(j));
         double dR = TMath::Sqrt(dphi*dphi+(pix_track_eta->at(i)-electron_eta->at(j))*(pix_track_eta->at(i)-electron_eta->at(j)));
         if (pix_track_pt->at(i) < 3.) 
         {
            h_dr_e_pixtrk_l3 -> Fill(dR);
         }
         else
         {
            h_dr_e_pixtrk_g3 -> Fill(dR);
         }
      }
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

   //          // h_truth_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //          h_truth_pt->Fill(truth_electron_pt->at(i));
   //          h_truth_eta->Fill(truth_electron_eta->at(i));
   //          for (size_t it = 0; it < size_reco_electron; it++)
   //          {
   //             double dr = TMath::Sqrt((electron_phi->at(it)-truth_electron_phi->at(i))*(electron_phi->at(it)-truth_electron_phi->at(i))+(electron_eta->at(it)-truth_electron_eta->at(i))*(electron_eta->at(it)-truth_electron_eta->at(i)));
   //             if (truth_electron_pt->at(i)>5.) h_dr_g5_e->Fill(dr);
   //             else h_dr_l5_e->Fill(dr);             
   //          }
   //          for (size_t it = 0; it < size_reco_pixtrk; it++)
   //          {
   //             double dr = TMath::Sqrt((pix_track_phi->at(it)-truth_electron_phi->at(i))*(pix_track_phi->at(it)-truth_electron_phi->at(i))+(pix_track_eta->at(it)-truth_electron_eta->at(i))*(pix_track_eta->at(it)-truth_electron_eta->at(i)));
   //             if (truth_electron_pt->at(i)>5.) h_dr_g5_pixtrk->Fill(dr);
   //             else h_dr_l5_pixtrk->Fill(dr);    
   //          }
   //          for (size_t it = 0; it < eg_cluster_phi->size(); it++)
   //          {
   //             double dr = TMath::Sqrt((eg_cluster_phi->at(it)-truth_electron_phi->at(i))*(eg_cluster_phi->at(it)-truth_electron_phi->at(i))+(eg_cluster_eta->at(it)-truth_electron_eta->at(i))*(eg_cluster_eta->at(it)-truth_electron_eta->at(i)));
   //             h_dr_egamma->Fill(dr);            
   //          }
            

            
   //          if (FitToTruth(0 ,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),electron_pt,electron_eta,electron_phi)) 
   //          {
   //             // h_match_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //             h_e_pt->Fill(truth_electron_pt->at(i));
   //             h_e_eta->Fill(truth_electron_eta->at(i));
   //          }
            
   //          if (FitToTruth(0,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),pix_track_pt,pix_track_eta,pix_track_phi)) 
   //          {
   //             // h_pixtrk_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //             h_pixtrk_pt->Fill(truth_electron_pt->at(i));
   //             h_pixtrk_eta->Fill(truth_electron_eta->at(i));

   //          }
   //          // if (FitToTruth(0,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),track_pt,track_eta,track_phi)) 
   //          // {
   //          //    h_trk_e_pt_eta->Fill(truth_electron_pt->at(i),truth_electron_eta->at(i));
   //          // }
   //          if (FitToTruth(1,truth_electron_pt->at(i),truth_electron_eta->at(i),truth_electron_phi->at(i),eg_cluster_et,eg_cluster_eta,eg_cluster_phi))
   //          {
   //             h_egamma_pt->Fill(truth_electron_pt->at(i));
   //             h_egamma_eta->Fill(truth_electron_eta->at(i));
   //          }


   //       }
         
         
         
   //       // if(truth_electron_pt->at(i) > 2.5) h_truth_e_eta->Fill(truth_electron_eta->at(i));







   // //       // Double_t dR_el_cut = 0.03; //Fig. 106 of https://cds.cern.ch/record/2703499/files/ATL-COM-PHYS-2019-1446.pdf
   // //       // Int_t index_e_match = -1;
   // //       // Int_t index_e_reco_match = -1;
   // //       // Double_t dR_el_min = 999.;
   // //       // for (Int_t j=0; j < size_reco_electron; j++){
   // //       //    Double_t delta_e_eta = truth_electron_eta->at(i)-electron_eta->at(j);
   // //       //    Double_t delta_e_phi = truth_electron_phi->at(i)-electron_phi->at(j);
   // //       //    Double_t dR_el = TMath::Sqrt( delta_e_eta*delta_e_eta + delta_e_phi*delta_e_phi);
            
   // //       //    if(dR_el<dR_el_min) {
   // //       //       dR_el_min = dR_el;
   // //       //        index_e_match = i; //revisit 
   // //       //        index_e_reco_match = j; //revisit 
   // //       //    }
   // //       // }//reco ele

         
         

   // //       // Int_t index_pixtrk_match = -1;
   // //       // Int_t index_pixtrk_reco_match = -1;
   // //       // Double_t dR_pixtrk_min = 999.;
   // //       // for (Int_t j=0; j < size_reco_pixtrk; j++){
   // //       //    Double_t delta_pixtrk_eta = truth_electron_eta->at(i)-pix_track_eta->at(j);
   // //       //    Double_t delta_pixtrk_phi = truth_electron_phi->at(i)-pix_track_phi->at(j);
   // //       //    Double_t dR_pixtrk = TMath::Sqrt( delta_pixtrk_eta*delta_pixtrk_eta + delta_pixtrk_phi*delta_pixtrk_phi);
            
   // //       //    if(dR_pixtrk<dR_pixtrk_min) {
   // //       //       dR_pixtrk_min = dR_pixtrk;
   // //       //        index_pixtrk_match = i; //revisit 
   // //       //        index_pixtrk_reco_match = j; //revisit 
   // //       //    }
   // //       // }//reco pixel tracks

         
      
   // //    // if(dR_el_min<dR_el_cut){
   // //    //    h_match_e_pt->Fill(truth_electron_pt->at(index_e_match));
   // //    //    if(electron_is_LHVeryLoose->at(index_e_reco_match)) h_lhveryloose_e_pt->Fill(truth_electron_pt->at(index_e_match));
   // //    //    if(truth_electron_pt->at(i) > 2.5) 
   // //    //    {
   // //    //       h_match_e_eta->Fill(truth_electron_eta->at(index_e_match));       
   // //    //    }
   // //    // }  

   // //    // if(dR_pixtrk_min<0.2){
   // //    //    h_match_pixtrk_pt->Fill(truth_electron_pt->at(index_pixtrk_match));
   // //    // }  
      
   // }//truth
   

   // }//isMC
   int tclu_mult = 0;
   int tclu_mult_tp_mc = 0;
   int tclu_mult_tp_data = 0;
   int pix_track_n = pix_track_charge -> size();
   // float beam_e = 2.68 * 1000;
   // zdc -> Fill(zdc_ene_a / beam_e, zdc_ene_c / beam_e);

   // int tclu_mult_nvl = 0;
   // //if((zdc_ene_a > 1000) and (zdc_ene_c > 1000)) continue;
   for (size_t i = 0; i < topo_cluster_pass_sig_cut->size(); i++)
      {
      if (topo_cluster_pass_sig_cut->at(i))
         {
            tclu_mult++;
         }
      }
   if (tclu_mult > 9) continue;
   int electron_n = electron_charge->size();
   if(electron_n > 0 and electron_n < 3)
      {
         for(int i = 0; i<electron_charge->size(); i++)
         {
            if(electron_is_LHTight->at(i)/* and electron_pt->at(i) > 2.0*/)
            {
               for (int j = 0; j<pix_track_charge->size(); j++)
               {
                  if((pix_track_hits->at(j) > 2) and ((electron_charge->at(i)*pix_track_charge->at(j)) < 0) and (TMath::Abs(pix_track_eta->at(j)) < 2.5 ))
                  {
                     TLorentzVector tag,probe;
                     tag.SetPtEtaPhiM(electron_pt -> at(i), electron_eta -> at(i), electron_phi->at(i),0.000511);
                     probe.SetPtEtaPhiM(pix_track_pt -> at(j), pix_track_eta -> at(j), pix_track_phi->at(j), 0.000511);
                     Double_t Aco = 1 - (TMath::Abs(DeltaPhi(tag.Phi(),probe.Phi()))/TMath::Pi());
                     double Minv = (tag+probe).M();
                     //Double_t weight = 1.;
                     double pt = (tag+probe).Pt();
                     if(Aco < 0.1 and Minv > 4.0 and pt < 2.)
                     {
                        if(isMC)
                        {
                        probe_d0_mc->Fill(pix_track_d0->at(j));
                        probe_z0_mc->Fill(pix_track_z0->at(j));
                        }
                        else
                        {
                        probe_d0_data->Fill(pix_track_d0->at(j));
                        probe_z0_data->Fill(pix_track_z0->at(j));
                        }
                        tag_probe_counter ++;
                        Double_t weight = L1TriggerWeightDown(tag.E(),probe.E()); 
                        if(isMC)
                        {
                           h_probe_mc_pt -> Fill(probe.Pt(),weight);
                           h_probe_mc_eta -> Fill(probe.Eta(),weight);
                           h_probe_mc_2d -> Fill(probe.Pt(), probe.Eta(), weight);
                        }
                        else
                        {
                           h_probe_data_pt -> Fill(probe.Pt());
                           h_probe_data_eta -> Fill(probe.Eta());
                           h_probe_data_2d -> Fill(probe.Pt(), probe.Eta());
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
                        for (size_t ite = 0; ite < topo_cluster_pass_sig_cut->size(); ite++)
                        {
                        if (topo_cluster_pass_sig_cut->at(ite))
                           {
                              if(isMC) tclu_mult_tp_mc++;
                              else tclu_mult_tp_data++;
                              if(isMC) h_tclu_pt_mc -> Fill(topo_cluster_pt->at(ite));
                              else h_tclu_pt_data -> Fill(topo_cluster_pt->at(ite));
                           }
                        }

                        double dR_min = 999.;
                        int ind = -1;
                        double dR_cut = 0.2;
                        if(probe.E() < 2.) dR_cut = 0.5;
                        for(int k = 0; k<electron_charge->size();k++)
                        {
                           TLorentzVector rec;
                           rec.SetPtEtaPhiM(electron_pt->at(k), electron_eta->at(k), electron_phi->at(k), 0.000511);
                           double dphi = DeltaPhi(probe.Phi(),rec.Phi());
                           double dR = sqrt(dphi*dphi+(probe.Eta()-rec.Eta())*(probe.Eta()-rec.Eta()));
                           if(isMC) h_dr_mc_probe -> Fill(dR, weight);
                           else h_dr_data_probe -> Fill(dR);
                           if(dR_min > dR) 
                           {
                              dR_min = dR;
                              ind = k;
                           }
                        }
                        if(dR_min < dR_cut)
                        {
                           if(isMC)
                           {
                              h_rec_mc_pt -> Fill(probe.Pt(), weight);
                              h_rec_mc_eta -> Fill(probe.Eta(), weight);
                              h_rec_mc_2d -> Fill(probe.Pt(), probe.Eta(), weight);
                              if(electron_is_LHVeryLoose->at(ind)){
                                 h_rec_op_mc_pt_vl -> Fill(probe.Pt(), weight);
                                 h_rec_op_mc_eta_vl -> Fill(probe.Eta(), weight);
                              }
                              if(electron_is_LHLoose->at(ind)){
                                 h_rec_op_mc_pt_lo -> Fill(probe.Pt(), weight);
                                 h_rec_op_mc_eta_lo -> Fill(probe.Eta(), weight);
                              }
                              if(electron_is_LHMedium->at(ind)) {
                                 h_rec_op_mc_pt_me -> Fill(probe.Pt(), weight);
                                 h_rec_op_mc_eta_me -> Fill(probe.Eta(), weight);
                              }
                              if(electron_is_LHTight->at(ind)) 
                              {
                                 h_rec_op_mc_pt_ti -> Fill(probe.Pt(), weight);
                                 h_rec_op_mc_eta_ti -> Fill(probe.Eta(), weight);
                              }
                           }
                           else
                           {
                              if(electron_is_LHVeryLoose->at(ind)) 
                              {
                              h_rec_op_data_pt_vl -> Fill(probe.Pt());
                              h_rec_op_data_eta_vl -> Fill(probe.Eta());
                              }
                              else
                              {
                                 h_cut_data_aco -> Fill(Aco);
                              }
                              if(electron_is_LHLoose->at(ind)){
                                 h_rec_op_data_pt_lo -> Fill(probe.Pt());
                                 h_rec_op_data_eta_lo -> Fill(probe.Eta());
                              }
                              if(electron_is_LHMedium->at(ind)) {
                                 h_rec_op_data_pt_me -> Fill(probe.Pt());
                                 h_rec_op_data_eta_me -> Fill(probe.Eta());
                              }
                              if(electron_is_LHTight->at(ind)) 
                              {
                                 h_rec_op_data_pt_ti -> Fill(probe.Pt());
                                 h_rec_op_data_eta_ti -> Fill(probe.Eta());
                              }
                              h_rec_data_pt -> Fill(probe.Pt());
                              h_rec_data_eta -> Fill(probe.Eta());
                              h_rec_data_2d -> Fill(probe.Pt(), probe.Eta());
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   //    h_tclu_mult -> Fill(tclu_mult);
      h_tclu_mult_tp_data -> Fill(tclu_mult_tp_data);
      h_tclu_mult_tp_mc -> Fill(tclu_mult_tp_mc);
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

   }

   
   h_rec_eff_data_pt -> Divide(h_rec_data_pt, h_probe_data_pt, 1. ,1., "b");
   h_rec_eff_data_eta -> Divide(h_rec_data_eta, h_probe_data_eta, 1. ,1., "b");
   h_rec_eff_mc_pt -> Divide(h_rec_mc_pt, h_probe_mc_pt, 1. ,1., "b");
   h_rec_eff_mc_eta -> Divide(h_rec_mc_eta, h_probe_mc_eta, 1. ,1., "b");
   h_ratio_pt -> Divide(h_rec_eff_data_pt, h_rec_eff_mc_pt, 1., 1., "b");
   h_ratio_eta -> Divide(h_rec_eff_data_eta, h_rec_eff_mc_eta, 1., 1., "b");
   h_rec_eff_data_pt_vl -> Divide(h_rec_op_data_pt_vl, h_rec_data_pt, 1., 1., "b");
   h_rec_eff_data_pt_lo -> Divide(h_rec_op_data_pt_lo, h_rec_data_pt, 1., 1., "b");
   h_rec_eff_data_pt_me -> Divide(h_rec_op_data_pt_me, h_rec_data_pt, 1., 1., "b");
   h_rec_eff_data_pt_ti -> Divide(h_rec_op_data_pt_ti, h_rec_data_pt, 1., 1., "b");
   h_rec_eff_mc_pt_vl -> Divide(h_rec_op_mc_pt_vl, h_rec_mc_pt, 1., 1., "b");
   h_rec_eff_mc_pt_lo -> Divide(h_rec_op_mc_pt_lo, h_rec_mc_pt, 1., 1., "b");
   h_rec_eff_mc_pt_me -> Divide(h_rec_op_mc_pt_me, h_rec_mc_pt, 1., 1., "b");
   h_rec_eff_mc_pt_ti -> Divide(h_rec_op_mc_pt_ti, h_rec_mc_pt, 1., 1., "b");
   h_rec_eff_data_eta_vl -> Divide(h_rec_op_data_eta_vl, h_rec_data_eta, 1., 1., "b");
   h_rec_eff_data_eta_lo -> Divide(h_rec_op_data_eta_lo, h_rec_data_eta, 1., 1., "b");
   h_rec_eff_data_eta_me -> Divide(h_rec_op_data_eta_me, h_rec_data_eta, 1., 1., "b");
   h_rec_eff_data_eta_ti -> Divide(h_rec_op_data_eta_ti, h_rec_data_eta, 1., 1., "b");
   h_rec_eff_mc_eta_vl -> Divide(h_rec_op_mc_eta_vl, h_rec_mc_eta, 1., 1., "b");
   h_rec_eff_mc_eta_lo -> Divide(h_rec_op_mc_eta_lo, h_rec_mc_eta, 1., 1., "b");
   h_rec_eff_mc_eta_me -> Divide(h_rec_op_mc_eta_me, h_rec_mc_eta, 1., 1., "b");
   h_rec_eff_mc_eta_ti -> Divide(h_rec_op_mc_eta_ti, h_rec_mc_eta, 1., 1., "b");
   h_rec_eff_mc_2d -> Divide(h_rec_mc_2d, h_probe_mc_2d, 1., 1., "b");
   h_rec_eff_data_2d -> Divide(h_rec_data_2d, h_probe_data_2d, 1., 1., "b");
   h_ratio_2d -> Divide(h_rec_eff_data_2d, h_rec_eff_mc_2d, 1., 1., "b");
   h_truth_e_eff_pt->Divide(h_match_e_pt,h_truth_e_pt, 1., 1., "b" );
   h_lhveryloose_e_eff_pt->Divide(h_lhveryloose_e_pt,h_truth_e_pt, 1., 1., "b" );
   h_truth_e_eff_eta->Divide(h_match_e_eta,h_truth_e_eta, 1., 1., "b" );
   h_truth_pixtrk_eff_pt->Divide(h_match_pixtrk_pt,h_truth_e_pt, 1., 1., "b" );
   h_eff_pixtrk->Divide(h_pixtrk_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   h_eff_electron->Divide(h_match_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   h_eff_trk->Divide(h_trk_e_pt_eta,h_truth_e_pt_eta, 1., 1., "b");
   h_eff_egamma_pt->Divide(h_egamma_pt,h_truth_pt,1.,1.,"b");
   h_eff_egamma_eta->Divide(h_egamma_eta,h_truth_eta,1.,1.,"b");
   h_eff_pixtrk_pt->Divide(h_pixtrk_pt,h_truth_pt,1.,1.,"b");
   h_eff_pixtrk_eta->Divide(h_pixtrk_eta,h_truth_eta,1.,1.,"b");
   h_eff_e_pt->Divide(h_e_pt,h_truth_pt,1.,1.,"b");
   h_eff_e_eta->Divide(h_e_eta,h_truth_eta,1.,1.,"b");
   outputFile->Write();
   outputFile->Close();
   //std::cout<<tag_probe_counter;
}

//TEffitiency