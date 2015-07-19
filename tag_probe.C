#include <iostream>
#include <vector>

#include<TFile.h>
#include <TChain.h>
#include <TMath.h>
#include<TH1F.h>
#include<TCanvas.h>
#include<TROOT.h>
#include<TGaxis.h>
#include<TLorentzVector.h>
#include<TLegend.h>
#include<TLatex.h>
void tag_probe(){

//  SetStyle();
//  SetPrelimStyle();
//  SetMeetingStyle();
//  gStyle->SetOptStat(0) ;
//  gStyle->SetPadTickX(1) ;
//  gStyle->SetPadTickY(1) ;
//  gStyle->SetFillStyle(ROOT.kWhite) ;
  TChain *chain = new TChain("IIHEAnalysis") ;
  chain->Add("/user/aidan/public/HEEP/samples/RunIISpring15DR74/DYToEE_50_50ns/outfile_*.root") ;
  
  int gsf_n=0;
  float ev_fixedGridRhoAll=0;
  std::vector<int>* gsf_charge = 0 ;
  std::vector<float>* gsf_px=0;
  std::vector<float>* gsf_py=0;
  std::vector<float>* gsf_pz=0;
  std::vector<float>* gsf_energy=0;
  std::vector<float>* gsf_pt=0;
  std::vector<float>* gsf_eta=0;
  std::vector<float>* gsf_phi=0;
//acceptance
  std::vector<float>* HEEP_cutflow60_Et_value=0;
  std::vector<float>* HEEP_cutflow60_eta_value=0;
  std::vector<bool> * HEEP_cutflow60_acceptance=0;                  
//ID
  std::vector<float>* HEEP_cutflow60_EcalDriven_value=0;
  std::vector<float>* HEEP_cutflow60_dEtaIn_value=0;
  std::vector<float>* HEEP_cutflow60_dPhiIn_value=0;
  std::vector<float>* HEEP_cutflow60_HOverE_value=0;
  std::vector<float>* HEEP_cutflow60_SigmaIetaIeta_value=0;
  std::vector<float>* HEEP_cutflow60_E1x5OverE5x5_value=0;
  std::vector<float>* HEEP_cutflow60_E2x5OverE5x5_value=0;
  std::vector<float>* HEEP_cutflow60_missingHits_value=0;
  std::vector<float>* HEEP_cutflow60_dxyFirstPV_value=0;
  std::vector<bool>* HEEP_cutflow60_ID=0;
//Iso
  std::vector<float>* HEEP_cutflow60_isolEMHadDepth1_value=0;
  std::vector<float>* HEEP_cutflow60_IsolPtTrks_value=0;
  std::vector<bool>* HEEP_cutflow60_isolation=0;

  std::vector<bool>* HEEP_cutflow60_total=0;


  std::vector<int>*  mc_pdgId=0;
  std::vector<int>*  mc_numberOfDaughters=0;
  std::vector<int>*  mc_mother_pdgId=0;
  std::vector<float>* mc_pt=0;
  std::vector<float>* mc_eta=0;  
  std::vector<float>* mc_phi=0;
  std::vector<float>* mc_energy=0;
                      

  TLorentzVector gsf_1;
  TLorentzVector gsf_2;
  TLorentzVector Z;
    
  TLorentzVector mc_e;
  TLorentzVector mc_E;
  TLorentzVector mc_Z;

  TLorentzVector HEEP_tag;
  TLorentzVector HEEP_probe;
  TLorentzVector HEEP_Z;

  float HEEP_theta_tag=0;
  float HEEP_E_tag=0;
  float HEEP_theta_probe=0;
  float HEEP_E_probe=0;

  bool HEEP_total_op=0;
  bool HEEP_accept_op=0;

  float max=0;
  int finde=0;
  int findE=0;
 
  int ID_j=0;
  int ID_k=0;
  int Iso_j=0;
  int Iso_k=0;
  int N_total=0;
  
  float moreone_Et=0;
  float moreone_M=0;
  float moreone_Eta=0;
  float two_Et=0;
  float two_M=0;
  float two_Eta=0;

  TH2F *H_M_Et_moreone = new TH2F("M_Et_moreone", "", 60, 60, 120, 5000, 0, 5000) ; 
  TH2F *H_M_Et_two = (TH2F*) H_M_Et_moreone->Clone("M_Et_two");
  TH2F *H_M_Eta_moreone = new TH2F("M_Eta_moreone", "", 60, 60, 120, 100, -2.5, 2.5) ;  
  TH2F *H_M_Eta_two = (TH2F*) H_M_Eta_moreone->Clone("M_Eta_two");  
  
  TH1F *H_M_moreone=new TH1F("M_moreone","",5001,-1,5000);
  TH1F *H_M_two=new TH1F("M_two","",5001,-1,5000);
  TH1F *H_Pt_gsf=new TH1F("Pt_gsf","",5001,-1,5000);
  TH1F *H_Eta_gsf=new TH1F("Eta_gsf","",100,-2.5,2.5);
  TH1F *H_Phi_gsf=new TH1F("Phi_gsf","",100,-5,5);

  TH1F *H_Et_HEEP = (TH1F*) H_Pt_gsf->Clone("Et_HEEP");
  TH1F *H_Eta_HEEP = (TH1F*) H_Eta_gsf->Clone("Eta_HEEP");
  TH1F *H_Phi_HEEP = (TH1F*) H_Phi_gsf->Clone("Phi_HEEP");

  TTree *T_moreone = new TTree("t_moreone","t_moreone");
  TTree *T_two = new TTree("t_two","t_two");
  
  T_moreone->Branch("moreone_Et",&moreone_Et,"moreone_Et/F");
  T_moreone->Branch("moreone_Eta",&moreone_Eta,"moreone_Eta/F");  
  T_moreone->Branch("moreone_M",&moreone_M,"moreone_M/F");
  T_two->Branch("two_Et",&two_Et,"two_Et/F");
  T_two->Branch("two_Eta",&two_Eta,"two_Eta/F");
  T_two->Branch("two_M",&two_M,"two_M/F"); 
//  chain->SetBranchAddress("gsf_charge", &gsf_charge) ;
  chain->SetBranchAddress("gsf_px", &gsf_px) ;
  chain->SetBranchAddress("gsf_py", &gsf_py) ;
  chain->SetBranchAddress("gsf_pz", &gsf_pz) ;
  chain->SetBranchAddress("gsf_energy", &gsf_energy) ;
  chain->SetBranchAddress("gsf_pt", &gsf_pt) ;
  chain->SetBranchAddress("gsf_eta", &gsf_eta) ;
  chain->SetBranchAddress("gsf_phi", &gsf_phi) ;
  chain->SetBranchAddress("gsf_n", &gsf_n) ;

  chain->SetBranchAddress("HEEP_cutflow60_Et_value", &HEEP_cutflow60_Et_value);
  chain->SetBranchAddress("HEEP_cutflow60_eta_value", &HEEP_cutflow60_eta_value);
  chain->SetBranchAddress("HEEP_cutflow60_acceptance", &HEEP_cutflow60_acceptance);
/*
  chain->SetBranchAddress("HEEP_cutflow60_EcalDriven_value", &HEEP_cutflow60_EcalDriven_value);
  chain->SetBranchAddress("HEEP_cutflow60_dEtaIn_value",&HEEP_cutflow60_dEtaIn_value);
  chain->SetBranchAddress("HEEP_cutflow60_dPhiIn_value", &HEEP_cutflow60_dPhiIn_value);
  chain->SetBranchAddress("HEEP_cutflow60_HOverE_value", &HEEP_cutflow60_HOverE_value);
  chain->SetBranchAddress("HEEP_cutflow60_SigmaIetaIeta_value", &HEEP_cutflow60_SigmaIetaIeta_value);
  chain->SetBranchAddress("HEEP_cutflow60_E1x5OverE5x5_value", &HEEP_cutflow60_E1x5OverE5x5_value);
  chain->SetBranchAddress("HEEP_cutflow60_E2x5OverE5x5_value", &HEEP_cutflow60_E2x5OverE5x5_value);
  chain->SetBranchAddress("HEEP_cutflow60_missingHits_value", &HEEP_cutflow60_missingHits_value);
  chain->SetBranchAddress("HEEP_cutflow60_dxyFirstPV_value", &HEEP_cutflow60_dxyFirstPV_value);
  chain->SetBranchAddress("HEEP_cutflow60_ID", &HEEP_cutflow60_ID);

  chain->SetBranchAddress("HEEP_cutflow60_isolEMHadDepth1_value", &HEEP_cutflow60_isolEMHadDepth1_value);
  chain->SetBranchAddress("HEEP_cutflow60_IsolPtTrks_value", &HEEP_cutflow60_IsolPtTrks_value);
  chain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll) ;
  chain->SetBranchAddress("HEEP_cutflow60_isolation", &HEEP_cutflow60_isolation);
*/
  chain->SetBranchAddress("HEEP_cutflow60_total", &HEEP_cutflow60_total);
/*
  chain->SetBranchAddress("mc_pdgId", &mc_pdgId) ;
  chain->SetBranchAddress("mc_numberOfDaughters", &mc_numberOfDaughters) ;
  chain->SetBranchAddress("mc_mother_pdgId", &mc_mother_pdgId) ;
  chain->SetBranchAddress("mc_energy", &mc_energy) ;
  chain->SetBranchAddress("mc_pt", &mc_pt) ;
  chain->SetBranchAddress("mc_eta", &mc_eta) ;
  chain->SetBranchAddress("mc_phi", &mc_phi) ;
*/



    long nEntries = chain->GetEntries() ;
    std::cout<<"Entries="<<nEntries<<std::endl;
    for(long long i=0 ; i<nEntries ; i++){
    chain->GetEntry(i) ;
    if(i%10000==0){std::cout<<"N="<<i<<std::endl;}


    for(int j=0; j<gsf_n;j++){
    
    H_Pt_gsf->Fill(gsf_pt->at(j));
    H_Eta_gsf->Fill(gsf_eta->at(j));
    H_Phi_gsf->Fill(gsf_phi->at(j));

    if(HEEP_cutflow60_total->at(j)){ 
    
    H_Et_HEEP->Fill(HEEP_cutflow60_Et_value->at(j));
    H_Eta_HEEP->Fill(HEEP_cutflow60_eta_value->at(j));
    H_Phi_HEEP->Fill(gsf_phi->at(j));
  
    HEEP_theta_tag = 2*atan(exp(-HEEP_cutflow60_eta_value->at(j)));   
    
    HEEP_E_tag = HEEP_cutflow60_Et_value->at(j)/(sin(HEEP_theta_tag));

    HEEP_tag.SetPtEtaPhiE(gsf_pt->at(j),HEEP_cutflow60_eta_value->at(j),gsf_phi->at(j),HEEP_E_tag);

    for(int k=0;k<gsf_n;k++){

    if(k==j) continue;

    if(HEEP_cutflow60_acceptance->at(k)){ 

    HEEP_theta_probe = 2*atan(exp(-HEEP_cutflow60_eta_value->at(k)));

    HEEP_E_probe = HEEP_cutflow60_Et_value->at(k)/(sin(HEEP_theta_probe));

    HEEP_probe.SetPtEtaPhiE(gsf_pt->at(k),HEEP_cutflow60_eta_value->at(k),gsf_phi->at(k),HEEP_E_probe);

    HEEP_Z = HEEP_tag + HEEP_probe;
   
//    if( (60 > HEEP_Z.M()) || (120 < HEEP_Z.M()) ) continue;
    
    H_M_moreone->Fill(HEEP_Z.M());  
    moreone_M=HEEP_Z.M(); 
    moreone_Et=HEEP_cutflow60_Et_value->at(k);
    moreone_Eta=HEEP_cutflow60_eta_value->at(k);
    T_moreone->Fill();
    H_M_Eta_moreone->Fill(HEEP_Z.M(),HEEP_probe.Eta());
    H_M_Et_moreone->Fill(HEEP_Z.M(),HEEP_probe.Et());
    if(HEEP_cutflow60_total->at(k)){ 
                                     H_M_two->Fill(HEEP_Z.M());
                                     H_M_Eta_two->Fill(HEEP_Z.M(),HEEP_probe.Eta());
                                     H_M_Et_two->Fill(HEEP_Z.M(),HEEP_probe.Et());
                                     
                                     two_M=HEEP_Z.M();
                                     two_Et=HEEP_cutflow60_Et_value->at(k);
                                     two_Eta=HEEP_cutflow60_eta_value->at(k);
                                     T_two->Fill();
                                    }
    
                                 }//accept
                                 }//for k
                                }//total
                             }//for j
 } //for entries
   TFile *f1=new TFile("Tag_probe_DYToEE_50ns.root","UPDATE");

   H_M_Et_moreone->Write();
   H_M_Eta_moreone->Write();
   H_M_Et_two->Write();
   H_M_Eta_two->Write();

   H_M_moreone->Write();
   H_M_two->Write();
   H_Pt_gsf->Write();
   H_Eta_gsf->Write();
   H_Phi_gsf->Write();
   H_Et_HEEP->Write();
   H_Eta_HEEP->Write();
   H_Phi_HEEP->Write();

   T_moreone->Write();
   T_two->Write();
   f1->Write();
}     
