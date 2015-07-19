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
void produce_hist(){

//  SetStyle();
//  SetPrelimStyle();
//  SetMeetingStyle();
//  gStyle->SetOptStat(0) ;
//  gStyle->SetPadTickX(1) ;
//  gStyle->SetPadTickY(1) ;
//  gStyle->SetFillStyle(ROOT.kWhite) ;
  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_probe_EgammaB2015.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_probe_EgammaB2015.root") ;

  

  float moreone_Et=0;
  float moreone_Eta=0;
  float moreone_M=0;
  float two_Et=0;
  float two_Eta=0;
  float two_M=0;

  TH1F *H_moreone_35_50_B=new TH1F("moreone_35_50_B","",60,60,120);
  TH1F *H_moreone_50_80_B = (TH1F*) H_moreone_35_50_B->Clone("moreone_50_80_B");
  TH1F *H_moreone_80_120_B = (TH1F*) H_moreone_35_50_B->Clone("moreone_80_120_B");
  TH1F *H_moreone_120_200_B = (TH1F*) H_moreone_35_50_B->Clone("moreone_120_200_B");
  TH1F *H_moreone_200_plus_B = (TH1F*) H_moreone_35_50_B->Clone("moreone_200_plus_B");
  TH1F *H_moreone_35_50_E=new TH1F("moreone_35_50_E","",60,60,120);
  TH1F *H_moreone_50_80_E = (TH1F*) H_moreone_35_50_E->Clone("moreone_50_80_E");
  TH1F *H_moreone_80_120_E = (TH1F*) H_moreone_35_50_E->Clone("moreone_80_120_E");
  TH1F *H_moreone_120_200_E = (TH1F*) H_moreone_35_50_E->Clone("moreone_120_200_E");
  TH1F *H_moreone_200_plus_E = (TH1F*) H_moreone_35_50_E->Clone("moreone_200_plus_E");


  TH1F *H_two_35_50_B=new TH1F("two_35_50_B","",60,60,120);
  TH1F *H_two_50_80_B = (TH1F*) H_two_35_50_B->Clone("two_50_80_B");
  TH1F *H_two_80_120_B = (TH1F*) H_two_35_50_B->Clone("two_80_120_B");
  TH1F *H_two_120_200_B = (TH1F*) H_two_35_50_B->Clone("two_120_200_B");
  TH1F *H_two_200_plus_B = (TH1F*) H_two_35_50_B->Clone("two_200_plus_B");
  TH1F *H_two_35_50_E=new TH1F("two_35_50_E","",60,60,120);
  TH1F *H_two_50_80_E = (TH1F*) H_two_35_50_E->Clone("two_50_80_E");
  TH1F *H_two_80_120_E = (TH1F*) H_two_35_50_E->Clone("two_80_120_E");
  TH1F *H_two_120_200_E = (TH1F*) H_two_35_50_E->Clone("two_120_200_E");
  TH1F *H_two_200_plus_E = (TH1F*) H_two_35_50_E->Clone("two_200_plus_E");


  chain_moreone->SetBranchAddress("moreone_Et", &moreone_Et) ;
  chain_moreone->SetBranchAddress("moreone_Eta", &moreone_Eta) ;
  chain_moreone->SetBranchAddress("moreone_M", &moreone_M) ;

  chain_two->SetBranchAddress("two_Et", &two_Et) ;
  chain_two->SetBranchAddress("two_Eta", &two_Eta) ;
  chain_two->SetBranchAddress("two_M", &two_M) ;


    long moreone_nEntries = chain_moreone->GetEntries() ;
    std::cout<<"moreone_Entries="<<moreone_nEntries<<std::endl;
    for(long long i=0 ; i<moreone_nEntries ; i++){
    chain_moreone->GetEntry(i) ;
    if(i%10000==0) std::cout<<"N="<<i<<std::endl;
     
    if(abs(moreone_Eta)<1.442){if( moreone_Et>35 && moreone_Et<50 ) H_moreone_35_50_B->Fill(moreone_M);
	                         if( moreone_Et>50 && moreone_Et<80 ) H_moreone_50_80_B->Fill(moreone_M);
					 if( moreone_Et>80 && moreone_Et<120 ) H_moreone_80_120_B->Fill(moreone_M);
					 if( moreone_Et>120 && moreone_Et<200 ) H_moreone_120_200_B->Fill(moreone_M); 
					 if( moreone_Et>200 ) H_moreone_200_plus_B->Fill(moreone_M);
                              }

    if( (abs(moreone_Eta)>1.56) && (abs(moreone_Eta)<2.5) ){
                               if( moreone_Et>35 && moreone_Et<50 ) H_moreone_35_50_E->Fill(moreone_M);
                               if( moreone_Et>50 && moreone_Et<80 ) H_moreone_50_80_E->Fill(moreone_M);
                               if( moreone_Et>80 && moreone_Et<120 ) H_moreone_80_120_E->Fill(moreone_M);
                               if( moreone_Et>120 && moreone_Et<200 ) H_moreone_120_200_E->Fill(moreone_M);
                               if( moreone_Et>200 ) H_moreone_200_plus_E->Fill(moreone_M);
                                                           }
                                         } //for moreone_entries


   long two_nEntries = chain_two->GetEntries() ;
   std::cout<<"two_Entries="<<two_nEntries<<std::endl;
   for(long long j=0 ; j<two_nEntries ; j++){
   chain_two->GetEntry(j) ;
   if(j%10000==0) std::cout<<"N="<<j<<std::endl;

   if(abs(two_Eta)<1.442){
	                        if( two_Et>35 && two_Et<50 ) H_two_35_50_B->Fill(two_M);
                              if( two_Et>50 && two_Et<80 ) H_two_50_80_B->Fill(two_M);
                              if( two_Et>80 && two_Et<120 ) H_two_80_120_B->Fill(two_M);
                              if( two_Et>120 && two_Et<200 ) H_two_120_200_B->Fill(two_M);
                              if( two_Et>200 ) H_two_200_plus_B->Fill(two_M);
                             }

   if( (abs(two_Eta)>1.56) && (abs(two_Eta)<2.5) ){
                              if( two_Et>35 && two_Et<50 ) H_two_35_50_E->Fill(two_M);
                              if( two_Et>50 && two_Et<80 ) H_two_50_80_E->Fill(two_M);
                              if( two_Et>80 && two_Et<120 ) H_two_80_120_E->Fill(two_M);
                              if( two_Et>120 && two_Et<200 ) H_two_120_200_E->Fill(two_M);
                              if( two_Et>200 ) H_two_200_plus_E->Fill(two_M);
                                                           }
					                     } //for two_entries



   TFile *f1=new TFile("Tag_probe_Data2015B_hist.root","UPDATE");
   
   H_moreone_35_50_B->Write();
   H_moreone_50_80_B->Write();
   H_moreone_80_120_B->Write();
   H_moreone_120_200_B->Write();
   H_moreone_200_plus_B->Write();
   H_moreone_35_50_E->Write();
   H_moreone_50_80_E->Write();
   H_moreone_80_120_E->Write();
   H_moreone_120_200_E->Write();
   H_moreone_200_plus_E->Write();

   H_two_35_50_B->Write();
   H_two_50_80_B->Write();
   H_two_80_120_B->Write();
   H_two_120_200_B->Write();
   H_two_200_plus_B->Write();
   H_two_35_50_E->Write();
   H_two_50_80_E->Write();
   H_two_80_120_E->Write();
   H_two_120_200_E->Write();
   H_two_200_plus_E->Write();
   f1->Write();
}     
