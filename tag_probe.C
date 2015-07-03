#include"/home/geant4/plotstyle_v4/bes3plotstyle_new.h"
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
void tag_probe_new(){

//  SetStyle();
//  SetPrelimStyle();
//  SetMeetingStyle();
//  gStyle->SetOptStat(0) ;
//  gStyle->SetPadTickX(1) ;
//  gStyle->SetPadTickY(1) ;
//  gStyle->SetFillStyle(ROOT.kWhite) ;
  TChain *chain = new TChain("IIHEAnalysis") ;
  chain->Add("DYtoEE_20bx25_0630.root") ;
  
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
  std::vector<float>* HEEP_cutflow51_Et_value=0;
  std::vector<float>* HEEP_cutflow51_eta_value=0;
  std::vector<bool> * HEEP_cutflow51_acceptance=0;                  
//ID
  std::vector<float>* HEEP_cutflow51_EcalDriven_value=0;
  std::vector<float>* HEEP_cutflow51_dEtaIn_value=0;
  std::vector<float>* HEEP_cutflow51_dPhiIn_value=0;
  std::vector<float>* HEEP_cutflow51_HOverE_value=0;
  std::vector<float>* HEEP_cutflow51_SigmaIetaIeta_value=0;
  std::vector<float>* HEEP_cutflow51_E1x5OverE5x5_value=0;
  std::vector<float>* HEEP_cutflow51_E2x5OverE5x5_value=0;
  std::vector<float>* HEEP_cutflow51_missingHits_value=0;
  std::vector<float>* HEEP_cutflow51_dxyFirstPV_value=0;
  std::vector<bool>* HEEP_cutflow51_ID=0;
//Iso
  std::vector<float>* HEEP_cutflow51_isolEMHadDepth1_value=0;
  std::vector<float>* HEEP_cutflow51_IsolPtTrks_value=0;
  std::vector<bool>* HEEP_cutflow51_isolation=0;

  std::vector<bool>* HEEP_cutflow51_total=0;


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
//  TH1F *H_gsf_E=new TH1F("gsf_E","energy of gsf_ Ee vs mc_Ee",750,0,1500);
//  TH1F *H_gsf_pt=new TH1F("gsf_pt_stotal","",100,0,100);
//  TH1F *H_gsf_eta=new TH1F("gsf_eta_stotal","",50,-2.5,2.5);

//  TH1F *H_gsf_phi=new TH1F("gsf_phi","phi of gsf_Ee vs mc_Ee",50,-5,5);
//  TH1F *hhlt2E1 = (TH1F*) hhlt_base->Clone("hhlt2E1");
//  TH1F *H_Z_m=new TH1F("Z_m_moreone","",60,60,120);
//  TH1F *H_Z_m_two=new TH1F("Z_m_two","",60,60,120);

//  TH1F *H_mc_Z=new TH1F("Z_mc","",120,0,120);
//  TH1F *H_HEEP_theta=new TH1F("HEEP_theta","",201,-1,200);

  TH1F *H_moreone_35_50=new TH1F("moreone_35_50","",60,60,120);
  TH1F *H_moreone_50_100 = (TH1F*) H_moreone_35_50->Clone("moreone_50_100");
  TH1F *H_moreone_100_5000 = (TH1F*) H_moreone_35_50->Clone("moreone_100_5000");
  TH1F *H_two_35_50 = (TH1F*) H_moreone_35_50->Clone("two_35_50");
  TH1F *H_two_50_100 = (TH1F*) H_moreone_35_50->Clone("two_50_100");
  TH1F *H_two_100_5000 = (TH1F*) H_moreone_35_50->Clone("two_100_5000");

  TH1F *H_moreone_barrel = (TH1F*) H_moreone_35_50->Clone("moreone_barrel");
  TH1F *H_moreone_endcap = (TH1F*) H_moreone_35_50->Clone("moreone_endcap");
  TH1F *H_two_barrel = (TH1F*) H_moreone_35_50->Clone("two_barrel");
  TH1F *H_two_endcap = (TH1F*) H_moreone_35_50->Clone("two_endcap");

  chain->SetBranchAddress("gsf_charge", &gsf_charge) ;
  chain->SetBranchAddress("gsf_px", &gsf_px) ;
  chain->SetBranchAddress("gsf_py", &gsf_py) ;
  chain->SetBranchAddress("gsf_pz", &gsf_pz) ;
  chain->SetBranchAddress("gsf_energy", &gsf_energy) ;
  chain->SetBranchAddress("gsf_pt", &gsf_pt) ;
  chain->SetBranchAddress("gsf_eta", &gsf_eta) ;
  chain->SetBranchAddress("gsf_phi", &gsf_phi) ;
  chain->SetBranchAddress("gsf_n", &gsf_n) ;

  chain->SetBranchAddress("HEEP_cutflow51_Et_value", &HEEP_cutflow51_Et_value);
  chain->SetBranchAddress("HEEP_cutflow51_eta_value", &HEEP_cutflow51_eta_value);
  chain->SetBranchAddress("HEEP_cutflow51_acceptance", &HEEP_cutflow51_acceptance);

  chain->SetBranchAddress("HEEP_cutflow51_EcalDriven_value", &HEEP_cutflow51_EcalDriven_value);
  chain->SetBranchAddress("HEEP_cutflow51_dEtaIn_value",&HEEP_cutflow51_dEtaIn_value);
  chain->SetBranchAddress("HEEP_cutflow51_dPhiIn_value", &HEEP_cutflow51_dPhiIn_value);
  chain->SetBranchAddress("HEEP_cutflow51_HOverE_value", &HEEP_cutflow51_HOverE_value);
  chain->SetBranchAddress("HEEP_cutflow51_SigmaIetaIeta_value", &HEEP_cutflow51_SigmaIetaIeta_value);
  chain->SetBranchAddress("HEEP_cutflow51_E1x5OverE5x5_value", &HEEP_cutflow51_E1x5OverE5x5_value);
  chain->SetBranchAddress("HEEP_cutflow51_E2x5OverE5x5_value", &HEEP_cutflow51_E2x5OverE5x5_value);
  chain->SetBranchAddress("HEEP_cutflow51_missingHits_value", &HEEP_cutflow51_missingHits_value);
  chain->SetBranchAddress("HEEP_cutflow51_dxyFirstPV_value", &HEEP_cutflow51_dxyFirstPV_value);
  chain->SetBranchAddress("HEEP_cutflow51_ID", &HEEP_cutflow51_ID);

  chain->SetBranchAddress("HEEP_cutflow51_isolEMHadDepth1_value", &HEEP_cutflow51_isolEMHadDepth1_value);
  chain->SetBranchAddress("HEEP_cutflow51_IsolPtTrks_value", &HEEP_cutflow51_IsolPtTrks_value);
  chain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll) ;
  chain->SetBranchAddress("HEEP_cutflow51_isolation", &HEEP_cutflow51_isolation);

  chain->SetBranchAddress("HEEP_cutflow51_total", &HEEP_cutflow51_total);

  chain->SetBranchAddress("mc_pdgId", &mc_pdgId) ;
  chain->SetBranchAddress("mc_numberOfDaughters", &mc_numberOfDaughters) ;
  chain->SetBranchAddress("mc_mother_pdgId", &mc_mother_pdgId) ;
  chain->SetBranchAddress("mc_energy", &mc_energy) ;
  chain->SetBranchAddress("mc_pt", &mc_pt) ;
  chain->SetBranchAddress("mc_eta", &mc_eta) ;
  chain->SetBranchAddress("mc_phi", &mc_phi) ;




    long nEntries = chain->GetEntries() ;
    std::cout<<"Entries="<<nEntries<<std::endl;
    for(long long i=0 ; i<nEntries ; i++){
    chain->GetEntry(i) ;
    if(i%10000==0){std::cout<<"N="<<i<<std::endl;}


    for(int j=0; j<gsf_n;j++){


    if(HEEP_cutflow51_total->at(j)){ 

    HEEP_theta_tag = 2*atan(exp(-HEEP_cutflow51_eta_value->at(j)));   
    
    HEEP_E_tag = HEEP_cutflow51_Et_value->at(j)/(sin(HEEP_theta_tag));

    HEEP_tag.SetPtEtaPhiE(gsf_pt->at(j),HEEP_cutflow51_eta_value->at(j),gsf_phi->at(j),HEEP_E_tag);

    for(int k=0;k<gsf_n;k++){

    if(k==j) continue;

    if(HEEP_cutflow51_acceptance->at(k)){ 

    HEEP_theta_probe = 2*atan(exp(-HEEP_cutflow51_eta_value->at(k)));

    HEEP_E_probe = HEEP_cutflow51_Et_value->at(k)/(sin(HEEP_theta_probe));

    HEEP_probe.SetPtEtaPhiE(gsf_pt->at(k),HEEP_cutflow51_eta_value->at(k),gsf_phi->at(k),HEEP_E_probe);

    HEEP_Z = HEEP_tag + HEEP_probe;
   
    if( (60 > HEEP_Z.M()) || (120 < HEEP_Z.M()) ) continue;

    if( fabs(HEEP_probe.Eta()) <1.5){H_moreone_barrel->Fill(HEEP_Z.M());
	                              if(HEEP_cutflow51_total->at(k)) H_two_barrel->Fill(HEEP_Z.M());
                                    }
    
    if( fabs(HEEP_probe.Eta())> 1.5 && fabs(HEEP_probe.Eta())<2.5 ){H_moreone_endcap->Fill(HEEP_Z.M()); 
	                                                             if(HEEP_cutflow51_total->at(k)) H_two_endcap->Fill(HEEP_Z.M());
							                          }



    if(HEEP_probe.Et()>35 && HEEP_probe.Et()<50){ 

    H_moreone_35_50->Fill(HEEP_Z.M());

    if(HEEP_cutflow51_total->at(k)) H_two_35_50->Fill(HEEP_Z.M()); 

    continue;
                                                }
     
    if(HEEP_probe.Et()>50 && HEEP_probe.Et()<100){

    H_moreone_50_100->Fill(HEEP_Z.M());

    if(HEEP_cutflow51_total->at(k)) H_two_50_100->Fill(HEEP_Z.M());

    continue;                                           
                                                 }


    if(HEEP_probe.Et()>100 && HEEP_probe.Et()<5000){

    H_moreone_100_5000->Fill(HEEP_Z.M());

    if(HEEP_cutflow51_total->at(k)) H_two_100_5000->Fill(HEEP_Z.M());

    continue;
                                                   }
                                 }//accept
                                 }//for k
                                }//total
                             }//for j
/*
  for(int ii=0; ii<mc_pdgId->size();ii++){
                                          if( mc_pdgId->at(ii)==11 && mc_numberOfDaughters->at(ii)==0 ){ 
					   mc_e.SetPtEtaPhiE(mc_pt->at(ii),mc_eta->at(ii),mc_phi->at(ii),mc_energy->at(ii));
					   finde=1;
							                                                             }
							if( mc_pdgId->at(ii)==-11 && mc_numberOfDaughters->at(ii)==0 ){
                                 mc_E.SetPtEtaPhiE(mc_pt->at(ii),mc_eta->at(ii),mc_phi->at(ii),mc_energy->at(ii));
					   findE=1;                                        
						      	                                                         }
							if(mc_pdgId->at(ii)==23) {
       				   mc_Z.SetPtEtaPhiE(mc_pt->at(ii),mc_eta->at(ii),mc_phi->at(ii),mc_energy->at(ii));
					   H_mc_Z->Fill(mc_Z.M());
							                         }
                                         }
*/  
/*  
if(finde && findE) {  mc_Z=mc_e+mc_E;
                      H_mc_Z->Fill(mc_Z.M());
                   }
finde=0;
findE=0;
*/
 } //for entries
/*
  H_gsf_pt->SetMarkerStyle(21) ;
  H_gsf_pt->SetMarkerColor(kRed) ;
  H_gsf_pt->SetLineColor(kRed) ;
  H_gsf_pt->Sumw2() ;

  H_gsf_eta->SetMarkerStyle(21) ;
  H_gsf_eta->SetMarkerColor(kRed) ;
  H_gsf_eta->SetLineColor(kRed) ;
  H_gsf_eta->Sumw2() ;
*/
/*
  H_Z_m->SetMarkerStyle(21) ;
  H_Z_m->SetMarkerColor(kRed) ;
  H_Z_m->SetLineColor(kRed) ;
  H_Z_m->Sumw2() ;
*/
   TFile *f1=new TFile("Tag_probe_0703.root","UPDATE");

/*
   TCanvas* c = new TCanvas("c", "c");
   NameAxes(H_gsf_pt, "Pt^{gsf} (GeV/c)", "Number / 1 GeV/c");
   H_gsf_pt->Draw();
   H_gsf_pt->Write();

   TCanvas* c1 = new TCanvas("c1", "c1");
   NameAxes(H_gsf_eta, "#eta^{gsf}", "Number / 0.1");
   H_gsf_eta->Draw();
   H_gsf_eta->Write();
*/
   
//    TCanvas* c = new TCanvas("c", "c");
//    NameAxes(H_Z_m, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
//    H_Z_m->Draw();
//    H_Z_m->Write();
/*
    TLegend legend(0.6, 0.65, 0.85, 0.5) ;
    legend.SetFillColor(kWhite) ;
    legend.SetBorderSize(0) ;
    legend.AddEntry(H_Z_m,  "Invariant mass of two opposite sign HEEP5.1(Accept+ID+Iso) electrons", "pe");
//    legend.SetTextSize(0.06);
    legend.SetTextFont(62);
    legend.Draw();

    TLatex* label = new TLatex( 0.1, 0.955,"CMS Interal #sqrt{s} = 13 TeV, L = 1 fb^{-1}") ;
    label->SetNDC() ;
    label->Draw();
*/
     TCanvas* c1 = new TCanvas("c1", "c1");
     NameAxes(H_moreone_35_50, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_moreone_35_50->Draw();
     H_moreone_35_50->Write();

     TCanvas* c2 = new TCanvas("c2", "c2");
     NameAxes(H_moreone_50_100, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_moreone_50_100->Draw();
     H_moreone_50_100->Write();

     TCanvas* c3 = new TCanvas("c3", "c3");
     NameAxes(H_moreone_100_5000, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_moreone_100_5000->Draw();
     H_moreone_100_5000->Write();

     TCanvas* c4 = new TCanvas("c4", "c4");
     NameAxes(H_two_35_50, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_two_35_50->Draw();
     H_two_35_50->Write();

     TCanvas* c5 = new TCanvas("c5", "c5");
     NameAxes(H_two_50_100, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_two_50_100->Draw();
     H_two_50_100->Write();

     TCanvas* c6 = new TCanvas("c6", "c6");
     NameAxes(H_two_100_5000, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_two_100_5000->Draw();
     H_two_100_5000->Write();


     TCanvas* c7 = new TCanvas("c7", "c7");
     NameAxes(H_moreone_barrel, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_moreone_barrel->Draw();
     H_moreone_barrel->Write();

     TCanvas* c8 = new TCanvas("c8", "c8");
     NameAxes(H_moreone_endcap, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_moreone_endcap->Draw();
     H_moreone_endcap->Write();

     TCanvas* c9 = new TCanvas("c9", "c9");
     NameAxes(H_two_barrel, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_two_barrel->Draw();
     H_two_barrel->Write();

     TCanvas* c10 = new TCanvas("c10", "c10");
     NameAxes(H_two_endcap, "M^{Z} (GeV/c^{2})", "Number / 1 GeV/c^{2}");
     H_two_endcap->Draw();
     H_two_endcap->Write();

    
   f1->Write();
}     
