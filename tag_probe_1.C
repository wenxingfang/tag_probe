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
void tag_probe_1(){

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


  float max=0;
  int finde=0;
  int findE=0;
 
  int ID_j=0;
  int ID_k=0;
  int Iso_j=0;
  int Iso_k=0;

  TH1F *H_gsf_E=new TH1F("gsf_E","energy of gsf_ Ee vs mc_Ee",750,0,1500);
  TH1F *H_gsf_pt=new TH1F("gsf_pt_accept","",100,0,100);
  TH1F *H_gsf_eta=new TH1F("gsf_eta_accept","",50,-2.5,2.5);
  TH1F *H_gsf_phi=new TH1F("gsf_phi","phi of gsf_Ee vs mc_Ee",50,-5,5);
//  TH1F *hhlt2E1 = (TH1F*) hhlt_base->Clone("hhlt2E1");
  TH1F *H_Z_m=new TH1F("Z_m_moreone","",120,0,120);

  TH1F *H_mc_Z=new TH1F("Z_mc","",120,0,120);

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
  chain->SetBranchAddress("HEEP_cutflow51_isolation", &HEEP_cutflow51_isolation);
  chain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll) ;

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
       if(HEEP_cutflow51_Et_value->at(j) < 35 || abs(HEEP_cutflow51_eta_value->at(j))>2.5 || ( abs(HEEP_cutflow51_eta_value->at(j))>1.442 && abs(HEEP_cutflow51_eta_value->at(j)) <1.56) ) continue;

	 if( abs(HEEP_cutflow51_eta_value->at(j)) < 1.442){
	    if( HEEP_cutflow51_EcalDriven_value->at(j)==1 && HEEP_cutflow51_dEtaIn_value->at(j)<0.004 && HEEP_cutflow51_dPhiIn_value->at(j)< 0.06 && HEEP_cutflow51_HOverE_value->at(j)< (0.5+2.5/gsf_energy->at(j)) && ( HEEP_cutflow51_E1x5OverE5x5_value->at(j)> 0.83 || HEEP_cutflow51_E2x5OverE5x5_value->at(j)> 0.94 ) && HEEP_cutflow51_missingHits_value->at(j)< 2 && HEEP_cutflow51_dxyFirstPV_value->at(j)<0.02 ) ID_j=1;
	    if(HEEP_cutflow51_IsolPtTrks_value->at(j)< 5){if( HEEP_cutflow51_isolEMHadDepth1_value->at(j) < (2+0.03*HEEP_cutflow51_Et_value->at(j)+0.28*ev_fixedGridRhoAll) ) Iso_j=1;}
	                                                 }
    
	 if( abs(HEEP_cutflow51_eta_value->at(j))>1.56 && abs(HEEP_cutflow51_eta_value->at(j))<2.5 ){
	    if( HEEP_cutflow51_EcalDriven_value->at(j)==1 && HEEP_cutflow51_dEtaIn_value->at(j)<0.006 && HEEP_cutflow51_dPhiIn_value->at(j)< 0.06 && HEEP_cutflow51_HOverE_value->at(j)< (0.5+12.5/gsf_energy->at(j)) && HEEP_cutflow51_SigmaIetaIeta_value->at(j)< 0.03 && HEEP_cutflow51_missingHits_value->at(j)< 2 && HEEP_cutflow51_dxyFirstPV_value->at(j)<0.05 ) ID_j=1;
	             if(HEEP_cutflow51_IsolPtTrks_value->at(j)< 5) {if( ( HEEP_cutflow51_Et_value->at(j)<50 && HEEP_cutflow51_isolEMHadDepth1_value->at(j)<(2.5+0.28*ev_fixedGridRhoAll) ) || (HEEP_cutflow51_Et_value->at(j)>50 && HEEP_cutflow51_isolEMHadDepth1_value->at(j)<(2.5+0.28*ev_fixedGridRhoAll+0.03*(HEEP_cutflow51_isolEMHadDepth1_value->at(j)-50)) ) ) Iso_j=1;}
			                                                                               }




//	                         if( (ID_j+Iso_j)==2 ){
					                       H_gsf_pt->Fill(HEEP_cutflow51_Et_value->at(j));
					                       H_gsf_eta->Fill(HEEP_cutflow51_eta_value->at(j));
//					               }

	                       ID_j=0;
				     Iso_j=0;
                             }
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

   TFile *f1=new TFile("Tag_probe_gsf_0630.root","UPDATE");


   TCanvas* c = new TCanvas("c", "c");
   NameAxes(H_gsf_pt, "Pt^{gsf} (GeV/c)", "Number / 1 GeV/c");
   H_gsf_pt->Draw();
   H_gsf_pt->Write();

   TCanvas* c1 = new TCanvas("c1", "c1");
   NameAxes(H_gsf_eta, "#eta^{gsf}", "Number / 0.1");
   H_gsf_eta->Draw();
   H_gsf_eta->Write();

   
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

//    TCanvas* c = new TCanvas("c", "c");
//    H_mc_Z->Draw();

   f1->Write();
}     
