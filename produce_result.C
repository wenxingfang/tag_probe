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
void produce_result(){

//  SetStyle();
//  SetPrelimStyle();
//  SetMeetingStyle();
//  gStyle->SetOptStat(0) ;
//  gStyle->SetPadTickX(1) ;
//  gStyle->SetPadTickY(1) ;
//  gStyle->SetFillStyle(ROOT.kWhite) ;

  float moreone_Et=0;
  float moreone_Eta=0;
  float moreone_M=0;
  float two_Et=0;
  float two_Eta=0;
  float two_M=0;
  
  float B_bin_low=0;
  float B_bin_high=0;
  float E_bin_low=0;
  float E_bin_high=0;

  const int B_N_bin=10;
  const int E_N_bin=4;
  float moreone_N_data_B[B_N_bin],moreone_N_data_B_err[B_N_bin],two_N_data_B[B_N_bin],two_N_data_B_err[B_N_bin];
  float moreone_N_data_E[E_N_bin],moreone_N_data_E_err[E_N_bin],two_N_data_E[E_N_bin],two_N_data_E_err[E_N_bin];
  int moreone_N_DY_B[B_N_bin],moreone_N_DY_B_err[B_N_bin],two_N_DY_B[B_N_bin],two_N_DY_B_err[B_N_bin];
  int moreone_N_DY_E[E_N_bin],moreone_N_DY_E_err[E_N_bin],two_N_DY_E[E_N_bin],two_N_DY_E_err[E_N_bin];
  int moreone_N_TT_1_B[B_N_bin],moreone_N_TT_1_B_err[B_N_bin],two_N_TT_1_B[B_N_bin],two_N_TT_1_B_err[B_N_bin];
  int moreone_N_TT_1_E[E_N_bin],moreone_N_TT_1_E_err[E_N_bin],two_N_TT_1_E[E_N_bin],two_N_TT_1_E_err[E_N_bin];
  int moreone_N_TT_2_B[B_N_bin],moreone_N_TT_2_B_err[B_N_bin],two_N_TT_2_B[B_N_bin],two_N_TT_2_B_err[B_N_bin];
  int moreone_N_TT_2_E[E_N_bin],moreone_N_TT_2_E_err[E_N_bin],two_N_TT_2_E[E_N_bin],two_N_TT_2_E_err[E_N_bin];
  int moreone_N_WJet_1_B[B_N_bin],moreone_N_WJet_1_B_err[B_N_bin],two_N_WJet_1_B[B_N_bin],two_N_WJet_1_B_err[B_N_bin];
  int moreone_N_WJet_1_E[E_N_bin],moreone_N_WJet_1_E_err[E_N_bin],two_N_WJet_1_E[E_N_bin],two_N_WJet_1_E_err[E_N_bin];
  int moreone_N_WJet_2_B[B_N_bin],moreone_N_WJet_2_B_err[B_N_bin],two_N_WJet_2_B[B_N_bin],two_N_WJet_2_B_err[B_N_bin];
  int moreone_N_WJet_2_E[E_N_bin],moreone_N_WJet_2_E_err[E_N_bin],two_N_WJet_2_E[E_N_bin],two_N_WJet_2_E_err[E_N_bin];


  float moreone_data_ex_B[B_N_bin],moreone_data_ex_B_err[B_N_bin],two_data_ex_B[B_N_bin],two_data_ex_B_err[B_N_bin];
  float moreone_data_ex_E[E_N_bin],moreone_data_ex_E_err[E_N_bin],two_data_ex_E[E_N_bin],two_data_ex_E_err[E_N_bin];
  float moreone_MC_ex_B[B_N_bin],moreone_MC_ex_B_err[B_N_bin],two_MC_ex_B[B_N_bin],two_MC_ex_B_err[B_N_bin];
  float moreone_MC_ex_E[E_N_bin],moreone_MC_ex_E_err[E_N_bin],two_MC_ex_E[E_N_bin],two_MC_ex_E_err[E_N_bin];
  float moreone_MC_in_B[B_N_bin],moreone_MC_in_B_err[B_N_bin],two_MC_in_B[B_N_bin],two_MC_in_B_err[B_N_bin];
  float moreone_MC_in_E[E_N_bin],moreone_MC_in_E_err[E_N_bin],two_MC_in_E[E_N_bin],two_MC_in_E_err[E_N_bin];

  float eff_data_in_B[B_N_bin],eff_data_in_B_err[B_N_bin],eff_data_ex_B[B_N_bin],eff_data_ex_B_err[B_N_bin];
  float eff_data_in_E[E_N_bin],eff_data_in_E_err[E_N_bin],eff_data_ex_E[E_N_bin],eff_data_ex_E_err[E_N_bin];
  float eff_MC_in_B[B_N_bin],eff_MC_in_B_err[B_N_bin],eff_MC_ex_B[B_N_bin],eff_MC_ex_B_err[B_N_bin];
  float eff_MC_in_E[E_N_bin],eff_MC_in_E_err[E_N_bin],eff_MC_ex_E[E_N_bin],eff_MC_ex_E_err[E_N_bin];

  float ratio_B_in[B_N_bin],ratio_B_in_err[B_N_bin],ratio_B_ex[B_N_bin],ratio_B_ex_err[B_N_bin];
  float ratio_E_in[E_N_bin],ratio_E_in_err[E_N_bin],ratio_E_ex[E_N_bin],ratio_E_ex_err[E_N_bin];

  float x_B[B_N_bin]={36.5,39.5,42.5,45.5,48.5,55,67.5,87.5,150,350};
  float x_B_err[B_N_bin]={1.5,1.5,1.5,1.5,1.5,5,7.5,12.5,50,150};

  float x_E[E_N_bin]={37.5,45,60,285};
  float x_E_err[E_N_bin]={2.5,5,10,215};


  for(int N_f=0;N_f<6;N_f++){

  if(N_f==0){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_Data2015B_DoubleEG_PromptReco_50ns.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_Data2015B_DoubleEG_PromptReco_50ns.root") ;
            }

  if(N_f==1){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_DYJetToLL_50ns.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_DYJetToLL_50ns.root") ;
            }

  if(N_f==2){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_TTbar_50ns_half1.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_TTbar_50ns_half1.root") ;
            }

  if(N_f==3){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_TTbar_50ns_half2.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_TTbar_50ns_half2.root") ;
            }

  if(N_f==4){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_WJetLNv_50ns_half1.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_WJetLNv_50ns_half1.root") ;
            }

  if(N_f==5){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_WJetLNv_50ns_half2.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_WJetLNv_50ns_half2.root") ;
            }

  chain_moreone->SetBranchAddress("moreone_Et", &moreone_Et) ;
  chain_moreone->SetBranchAddress("moreone_Eta", &moreone_Eta) ;
  chain_moreone->SetBranchAddress("moreone_M", &moreone_M) ;

  chain_two->SetBranchAddress("two_Et", &two_Et) ;
  chain_two->SetBranchAddress("two_Eta", &two_Eta) ;
  chain_two->SetBranchAddress("two_M", &two_M) ;


  for(int N_bin=0;N_bin<10;N_bin++){
   switch(N_bin){
                case 0:
			 B_bin_low=35;
			 B_bin_high=38;
			 E_bin_low=35;
			 E_bin_high=40;
			 break;
                case 1:
			 B_bin_low=38;
			 B_bin_high=41;
			 E_bin_low=40;
			 E_bin_high=50;
			 break;
                case 2:
			 B_bin_low=41;
			 B_bin_high=44;
			 E_bin_low=50;
			 E_bin_high=70;
			 break;
                case 3:
			 B_bin_low=44;
			 B_bin_high=47;
			 E_bin_low=70;
			 E_bin_high=5000;
			 break;
                case 4:
			 B_bin_low=47;
			 B_bin_high=50;
			 break;
                case 5:
			 B_bin_low=50;
			 B_bin_high=60;
			 break;
                case 6:
			 B_bin_low=60;
			 B_bin_high=75;
			 break;
                case 7:
			 B_bin_low=75;
			 B_bin_high=100;
			 break;
                case 8:
			 B_bin_low=100;
			 B_bin_high=200;
			 break;
                case 9:
			 B_bin_low=200;
			 B_bin_high=5000;
			 break;
                }

  TH1F *H_moreone_B=new TH1F("moreone_B","",60,60,120);
  TH1F *H_moreone_E=new TH1F("moreone_E","",60,60,120);

  TH1F *H_two_B=new TH1F("two_B","",60,60,120);
  TH1F *H_two_E=new TH1F("two_E","",60,60,120);

    long moreone_nEntries = chain_moreone->GetEntries() ;
    std::cout<<"moreone_Entries="<<moreone_nEntries<<std::endl;
    for(long long i=0 ; i<moreone_nEntries ; i++){
    chain_moreone->GetEntry(i) ;
//    if(i%10000==0) std::cout<<"N="<<i<<std::endl;
     
    if(abs(moreone_Eta)<1.442){if( moreone_Et>B_bin_low && moreone_Et<B_bin_high ) H_moreone_B->Fill(moreone_M);
                              }

    if( N_bin<4 && (abs(moreone_Eta)>1.56) && (abs(moreone_Eta)<2.5) ){
                               if( moreone_Et>E_bin_low && moreone_Et<E_bin_high ) H_moreone_E->Fill(moreone_M);
                                                           }
                                         } //for moreone_entries


   long two_nEntries = chain_two->GetEntries() ;
   std::cout<<"two_Entries="<<two_nEntries<<std::endl;
   for(long long j=0 ; j<two_nEntries ; j++){
   chain_two->GetEntry(j) ;
//   if(j%10000==0) std::cout<<"N="<<j<<std::endl;

    if(abs(two_Eta)<1.442){if( two_Et>B_bin_low && two_Et<B_bin_high ) H_two_B->Fill(two_M);
                              }

    if( N_bin<4 && (abs(two_Eta)>1.56) && (abs(two_Eta)<2.5) ){
                               if( two_Et>E_bin_low && two_Et<E_bin_high ) H_two_E->Fill(two_M);
                                                           }

  					                     } //for two_entries
   if(N_bin<10){
	         switch(N_f){
	                       case 0:
					  moreone_N_data_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_data_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_data_B[N_bin]=H_two_B->Integral();
					  two_N_data_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
	                       case 1:
					  moreone_N_DY_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_DY_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_DY_B[N_bin]=H_two_B->Integral();
					  two_N_DY_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
	                       case 2:
					  moreone_N_TT_1_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_TT_1_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_TT_1_B[N_bin]=H_two_B->Integral();
					  two_N_TT_1_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
	                       case 3:
					  moreone_N_TT_2_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_TT_2_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_TT_2_B[N_bin]=H_two_B->Integral();
					  two_N_TT_2_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
	                       case 4:
					  moreone_N_WJet_1_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_WJet_1_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_WJet_1_B[N_bin]=H_two_B->Integral();
					  two_N_WJet_1_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
	                       case 5:
					  moreone_N_WJet_2_B[N_bin]=H_moreone_B->Integral();
					  moreone_N_WJet_2_B_err[N_bin]=sqrt(H_moreone_B->Integral());
					  two_N_WJet_2_B[N_bin]=H_two_B->Integral();
					  two_N_WJet_2_B_err[N_bin]=sqrt(H_two_B->Integral());
					  break;
                           }
               } 
   if(N_bin<4){
	         switch(N_f){
	                       case 0:
					  moreone_N_data_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_data_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_data_E[N_bin]=H_two_E->Integral();
					  two_N_data_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
	                       case 1:
					  moreone_N_DY_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_DY_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_DY_E[N_bin]=H_two_E->Integral();
					  two_N_DY_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
	                       case 2:
					  moreone_N_TT_1_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_TT_1_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_TT_1_E[N_bin]=H_two_E->Integral();
					  two_N_TT_1_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
	                       case 3:
					  moreone_N_TT_2_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_TT_2_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_TT_2_E[N_bin]=H_two_E->Integral();
					  two_N_TT_2_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
	                       case 4:
					  moreone_N_WJet_1_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_WJet_1_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_WJet_1_E[N_bin]=H_two_E->Integral();
					  two_N_WJet_1_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
	                       case 5:
					  moreone_N_WJet_2_E[N_bin]=H_moreone_E->Integral();
					  moreone_N_WJet_2_E_err[N_bin]=sqrt(H_moreone_E->Integral());
					  two_N_WJet_2_E[N_bin]=H_two_E->Integral();
					  two_N_WJet_2_E_err[N_bin]=sqrt(H_two_E->Integral());
					  break;
                           }
              }
  H_moreone_B->Delete();
  H_moreone_E->Delete();
  H_two_B->Delete();
  H_two_E->Delete();
  }//for N_bin
  
  
  }

 for(int L=0;L<10;L++){
                      moreone_data_ex_B[L]=moreone_N_data_B[L]-0.008366*(moreone_N_TT_1_B[L]+moreone_N_TT_2_B[L])-0.127428*(moreone_N_WJet_1_B[L]+moreone_N_WJet_2_B[L]);
			    moreone_data_ex_B_err[L]=sqrt(moreone_data_ex_B[L]);
                      two_data_ex_B[L]=two_N_data_B[L]-0.008366*(two_N_TT_1_B[L]+two_N_TT_2_B[L])-0.127428*(two_N_WJet_1_B[L]+two_N_WJet_2_B[L]);
			    two_data_ex_B_err[L]=sqrt(two_data_ex_B[L]);
                      moreone_MC_in_B[L]=0.02085*moreone_N_DY_B[L]+0.008366*(moreone_N_TT_1_B[L]+moreone_N_TT_2_B[L])+0.127428*(moreone_N_WJet_1_B[L]+moreone_N_WJet_2_B[L]);
			    moreone_MC_in_B_err[L]=sqrt(moreone_MC_in_B[L]);
                      two_MC_in_B[L]=0.02085*two_N_DY_B[L]+0.008366*(two_N_TT_1_B[L]+two_N_TT_2_B[L])+0.127428*(two_N_WJet_1_B[L]+two_N_WJet_2_B[L]);
			    two_MC_in_B_err[L]=sqrt(two_MC_in_B[L]);
                      moreone_MC_ex_B[L]=0.02085*moreone_N_DY_B[L]+0.016732*(moreone_N_TT_1_B[L]-moreone_N_TT_2_B[L])+0.254856*(moreone_N_WJet_1_B[L]-moreone_N_WJet_2_B[L]);
			    moreone_MC_ex_B_err[L]=sqrt(moreone_MC_ex_B[L]);
                      two_MC_ex_B[L]=0.02085*two_N_DY_B[L]+0.016732*(two_N_TT_1_B[L]-two_N_TT_2_B[L])+0.254856*(two_N_WJet_1_B[L]-two_N_WJet_2_B[L]);
			    two_MC_ex_B_err[L]=sqrt(two_MC_ex_B[L]);
                      }

 for(int M=0;M<4;M++){
                      moreone_data_ex_E[M]=moreone_N_data_E[M]-0.008366*(moreone_N_TT_1_E[M]+moreone_N_TT_2_E[M])-0.127428*(moreone_N_WJet_1_E[M]+moreone_N_WJet_2_E[M]);
			    moreone_data_ex_E_err[M]=sqrt(moreone_data_ex_E[M]);
                      two_data_ex_E[M]=two_N_data_E[M]-0.008366*(two_N_TT_1_E[M]+two_N_TT_2_E[M])-0.127428*(two_N_WJet_1_E[M]+two_N_WJet_2_E[M]);
			    two_data_ex_E_err[M]=sqrt(two_data_ex_E[M]);
                      moreone_MC_in_E[M]=0.02085*moreone_N_DY_E[M]+0.008366*(moreone_N_TT_1_E[M]+moreone_N_TT_2_E[M])+0.127428*(moreone_N_WJet_1_E[M]+moreone_N_WJet_2_E[M]);
			    moreone_MC_in_E_err[M]=sqrt(moreone_MC_in_E[M]);
                      two_MC_in_E[M]=0.02085*two_N_DY_E[M]+0.008366*(two_N_TT_1_E[M]+two_N_TT_2_E[M])+0.127428*(two_N_WJet_1_E[M]+two_N_WJet_2_E[M]);
			    two_MC_in_E_err[M]=sqrt(two_MC_in_E[M]);
                      moreone_MC_ex_E[M]=0.02085*moreone_N_DY_E[M]+0.016732*(moreone_N_TT_1_E[M]-moreone_N_TT_2_E[M])+0.254856*(moreone_N_WJet_1_E[M]-moreone_N_WJet_2_E[M]);
			    moreone_MC_ex_E_err[M]=sqrt(moreone_MC_ex_E[M]);
                      two_MC_ex_E[M]=0.02085*two_N_DY_E[M]+0.016732*(two_N_TT_1_E[M]-two_N_TT_2_E[M])+0.254856*(two_N_WJet_1_E[M]-two_N_WJet_2_E[M]);
			    two_MC_ex_E_err[M]=sqrt(two_MC_ex_E[M]);
                      }
//////////////////////////////////////efficiency///////////////////////////////////////////////////////
 for(int N=0;N<10;N++){
                    if(moreone_N_data_B[N]!=0){ 
			    eff_data_in_B[N]=two_N_data_B[N]/moreone_N_data_B[N];
			    eff_data_in_B_err[N]=sqrt(eff_data_in_B[N]*(1-eff_data_in_B[N])/moreone_N_data_B[N]);
			                            }
                    else{eff_data_in_B[N]=0;
			       eff_data_in_B_err[N]=0;
			      }

                    if(moreone_data_ex_B[N]!=0){ 
			    eff_data_ex_B[N]=two_data_ex_B[N]/moreone_data_ex_B[N];
			    eff_data_ex_B_err[N]=sqrt(eff_data_ex_B[N]*(1-eff_data_ex_B[N])/moreone_data_ex_B[N]);
			                            }
                    else{eff_data_ex_B[N]=0;
			       eff_data_ex_B_err[N]=0;
			      }
                    
			  if(moreone_MC_ex_B[N]!=0){ 
			    eff_MC_ex_B[N]=two_MC_ex_B[N]/moreone_MC_ex_B[N];
			    eff_MC_ex_B_err[N]=sqrt(eff_MC_ex_B[N]*(1-eff_MC_ex_B[N])/moreone_MC_ex_B[N]);
			                            }
                    else{eff_MC_ex_B[N]=0;
			       eff_MC_ex_B_err[N]=0;
			      }
			 
			  if(moreone_MC_in_B[N]!=0){ 
			    eff_MC_in_B[N]=two_MC_in_B[N]/moreone_MC_in_B[N];
			    eff_MC_in_B_err[N]=sqrt(eff_MC_in_B[N]*(1-eff_MC_in_B[N])/moreone_MC_in_B[N]);
			                            }
                    else{eff_MC_in_B[N]=0;
			       eff_MC_in_B_err[N]=0;
			      }

                      }


 for(int P=0;P<4;P++){
                    if(moreone_N_data_E[P]!=0){ 
			    eff_data_in_E[P]=two_N_data_E[P]/moreone_N_data_E[P];
			    eff_data_in_E_err[P]=sqrt(eff_data_in_E[P]*(1-eff_data_in_E[P])/moreone_N_data_E[P]);
			                            }
                    else{eff_data_in_E[P]=0;
			       eff_data_in_E_err[P]=0;
			      }

                    if(moreone_data_ex_E[P]!=0){ 
			    eff_data_ex_E[P]=two_data_ex_E[P]/moreone_data_ex_E[P];
			    eff_data_ex_E_err[P]=sqrt(eff_data_ex_E[P]*(1-eff_data_ex_E[P])/moreone_data_ex_E[P]);
			                            }
                    else{eff_data_ex_E[P]=0;
			       eff_data_ex_E_err[P]=0;
			      }
                    
			  if(moreone_MC_ex_E[P]!=0){ 
			    eff_MC_ex_E[P]=two_MC_ex_E[P]/moreone_MC_ex_E[P];
			    eff_MC_ex_E_err[P]=sqrt(eff_MC_ex_E[P]*(1-eff_MC_ex_E[P])/moreone_MC_ex_E[P]);
			                            }
                    else{eff_MC_ex_E[P]=0;
			       eff_MC_ex_E_err[P]=0;
			      }
			 
			  if(moreone_MC_in_E[P]!=0){ 
			    eff_MC_in_E[P]=two_MC_in_E[P]/moreone_MC_in_E[P];
			    eff_MC_in_E_err[P]=sqrt(eff_MC_in_E[P]*(1-eff_MC_in_E[P])/moreone_MC_in_E[P]);
			                            }
                    else{eff_MC_in_E[P]=0;
			       eff_MC_in_E_err[P]=0;
			      }

                      }
///////////////////////////////////////////////Data/MC ratio////////////////////////////////////////////////
 for(int RB=0;RB<10;RB++){
                         if(eff_MC_in_B[RB]!=0){
				    ratio_B_in[RB]=eff_data_in_B[RB]/eff_MC_in_B[RB];
                            ratio_B_in_err[RB]=sqrt(pow(eff_data_in_B_err[RB]/eff_MC_in_B[RB],2)+pow(eff_data_in_B[RB]*eff_MC_in_B_err[RB],2)/pow(eff_MC_in_B[RB],4));
				                        }                  
                         else{ratio_B_in[RB]=0;
				      ratio_B_in_err[RB]=0;
				     }

                         if(eff_MC_ex_B[RB]!=0){
				    ratio_B_ex[RB]=eff_data_ex_B[RB]/eff_MC_ex_B[RB];
                            ratio_B_ex_err[RB]=sqrt(pow(eff_data_ex_B_err[RB]/eff_MC_ex_B[RB],2)+pow(eff_data_ex_B[RB]*eff_MC_ex_B_err[RB],2)/pow(eff_MC_ex_B[RB],4));
				                        }                  
                         else{ratio_B_ex[RB]=0;
				      ratio_B_ex_err[RB]=0;
				     }
                         }

 for(int RE=0;RE<4;RE++){
                         if(eff_MC_in_E[RE]!=0){
				    ratio_E_in[RE]=eff_data_in_E[RE]/eff_MC_in_E[RE];
                            ratio_E_in_err[RE]=sqrt(pow(eff_data_in_E_err[RE]/eff_MC_in_E[RE],2)+pow(eff_data_in_E[RE]*eff_MC_in_E_err[RE],2)/pow(eff_MC_in_E[RE],4));
				                        }                  
                         else{ratio_E_in[RE]=0;
				      ratio_E_in_err[RE]=0;
				     }

                         if(eff_MC_ex_E[RE]!=0){
				    ratio_E_ex[RE]=eff_data_ex_E[RE]/eff_MC_ex_E[RE];
                            ratio_E_ex_err[RE]=sqrt(pow(eff_data_ex_E_err[RE]/eff_MC_ex_E[RE],2)+pow(eff_data_ex_E[RE]*eff_MC_ex_E_err[RE],2)/pow(eff_MC_ex_E[RE],4));
				                        }                  
                         else{ratio_E_ex[RE]=0;
				      ratio_E_ex_err[RE]=0;
				     }
                         }



///////////////////////////////////////////////////////


//////////////////////////////////////////////////Draw/////////////////////////////////////////////////////////////

float size = 0.25 ;

TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
c1->SetFillColor(42);
c1->SetGrid();

c1->cd() ;
TPad *pad1 = new TPad("pad1", "", 0.0, size, 1.0, 1.0, 0) ;
TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, size, 0) ;
      
pad1->Draw() ;
pad2->Draw() ;

pad1->cd() ;
pad1->SetGridy() ;
pad1->SetGridx() ;
pad1->SetLogx();

TGraphErrors *gr1 = new TGraphErrors(E_N_bin,x_E,eff_MC_in_E,x_E_err,eff_MC_in_E_err);
gr1->SetName("gr1");
gr1->SetTitle("");
gr1->SetMarkerColor(1);
gr1->SetLineColor(2);
gr1->SetMarkerStyle(26);
gr1->SetLineWidth(2);
gr1->SetMarkerSize(1.3);
gr1->SetMaximum(1);
gr1->SetMinimum(0.5);
gr1->GetXaxis()->SetMoreLogLabels();

gr1->Draw("AP");

TGraphErrors *gr2 = new TGraphErrors(E_N_bin,x_E,eff_MC_ex_E,x_E_err,eff_MC_ex_E_err);
gr2->SetMarkerColor(1);
gr2->SetLineColor(4);
gr2->SetName("gr2");
gr2->SetMarkerStyle(32);
gr2->SetLineWidth(2);
gr2->SetMarkerSize(1.3);
gr2->Draw("P");

TGraphErrors *gr3 = new TGraphErrors(E_N_bin,x_E,eff_data_in_E,x_E_err,eff_data_in_E_err);
gr3->SetMarkerColor(1);
gr3->SetLineColor(1);
gr3->SetName("gr3");
gr3->SetMarkerStyle(24);
gr3->SetLineWidth(2);
gr3->SetMarkerSize(1.3);
gr3->Draw("P");

TGraphErrors *gr4 = new TGraphErrors(E_N_bin,x_E,eff_data_ex_E,x_E_err,eff_data_ex_E_err);
gr4->SetMarkerColor(1);
gr4->SetLineColor(6);
gr4->SetName("gr4");
gr4->SetMarkerStyle(25);
gr4->SetLineWidth(2);
gr4->SetMarkerSize(1.3);
gr4->Draw("P");

leg = new TLegend(0.1,0.7,0.3,0.9);

leg->AddEntry("gr1","MC inclusive","lep");
leg->AddEntry("gr2","MC exclusive","lep");
leg->AddEntry("gr3","data","lep");
leg->AddEntry("gr4","data exclusive","lep");
leg->Draw();

TLatex* Label = new TLatex(0.1 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV, 50 ns, Run 2015B DoubleEG PromptReco 49.98 pb^{-1}");
Label->SetTextSize(0.05);
Label->SetNDC();
Label->Draw();

TLatex* Label1 = new TLatex(0.2 , 0.85 , "Endcap E_{T}>35GeV");
Label1->SetTextSize(0.05);
Label1->SetNDC();
Label1->Draw();

pad2->cd() ;
pad2->SetGridy() ;
pad2->SetGridx() ;
pad2->SetLogx();

TGraphErrors *gr5 = new TGraphErrors(E_N_bin,x_E,ratio_E_in,x_E_err,ratio_E_in_err);
gr5->SetTitle("");
gr5->SetMarkerColor(1);
gr5->SetLineColor(2);
gr5->SetName("gr5");
gr5->SetMarkerStyle(26);
gr5->SetLineWidth(2);
gr5->SetMarkerSize(1.3);
gr5->SetMinimum(0.8);
gr5->Draw("AP");

TGraphErrors *gr6 = new TGraphErrors(E_N_bin,x_E,ratio_E_ex,x_E_err,ratio_E_ex_err);
gr6->SetMarkerColor(1);
gr6->SetLineColor(4);
gr6->SetName("gr6");
gr6->SetMarkerStyle(32);
gr6->SetLineWidth(2);
gr6->SetMarkerSize(1.3);
gr6->Draw("P");

leg1 = new TLegend(0.1,0.5,0.4,0.9);
leg1->AddEntry("gr5","Data / MC inclusive","lep");
leg1->AddEntry("gr6","Data exclusive / MC exclusive","lep");
leg1->Draw();
 
 
 // TFile *f1=new TFile("Tag_B_probe_Data2015B_DoubleEGPromptReco_hist.root","UPDATE");
  // f1->Write();
}     
