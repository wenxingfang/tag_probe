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
void produce_result_eta(){

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
  
  float eta_low=0;
  float eta_high=0;

  const int N_bin_eta=30;
  float moreone_N_data[N_bin_eta],moreone_N_data_err[N_bin_eta],two_N_data[N_bin_eta],two_N_data_err[N_bin_eta];
  int moreone_N_DY[N_bin_eta],moreone_N_DY_err[N_bin_eta],two_N_DY[N_bin_eta],two_N_DY_err[N_bin_eta];
  int moreone_N_TT_1[N_bin_eta],moreone_N_TT_1_err[N_bin_eta],two_N_TT_1[N_bin_eta],two_N_TT_1_err[N_bin_eta];
  int moreone_N_TT_2[N_bin_eta],moreone_N_TT_2_err[N_bin_eta],two_N_TT_2[N_bin_eta],two_N_TT_2_err[N_bin_eta];
  int moreone_N_WJet_1[N_bin_eta],moreone_N_WJet_1_err[N_bin_eta],two_N_WJet_1[N_bin_eta],two_N_WJet_1_err[N_bin_eta];
  int moreone_N_WJet_2[N_bin_eta],moreone_N_WJet_2_err[N_bin_eta],two_N_WJet_2[N_bin_eta],two_N_WJet_2_err[N_bin_eta];


  float moreone_data_ex[N_bin_eta],moreone_data_ex_err[N_bin_eta],two_data_ex[N_bin_eta],two_data_ex_err[N_bin_eta];
  float moreone_MC_ex[N_bin_eta],moreone_MC_ex_err[N_bin_eta],two_MC_ex[N_bin_eta],two_MC_ex_err[N_bin_eta];
  float moreone_MC_in[N_bin_eta],moreone_MC_in_err[N_bin_eta],two_MC_in[N_bin_eta],two_MC_in_err[N_bin_eta];

  float eff_data_in[N_bin_eta],eff_data_in_err[N_bin_eta],eff_data_ex[N_bin_eta],eff_data_ex_err[N_bin_eta];
  float eff_MC_in[N_bin_eta],eff_MC_in_err[N_bin_eta],eff_MC_ex[N_bin_eta],eff_MC_ex_err[N_bin_eta];

  float ratio_in[N_bin_eta],ratio_in_err[N_bin_eta],ratio_ex[N_bin_eta],ratio_ex_err[N_bin_eta];

  float x[N_bin_eta];
  float x_err[N_bin_eta];

  for(int ii=0;ii<N_bin_eta;ii++){x[ii]=-2.9+0.2*ii;
                                  x_err[ii]=0.1;
                                 }


  for(int N_f=0;N_f<6;N_f++){

  if(N_f==0){

  TChain *chain_moreone = new TChain("t_moreone") ;
  chain_moreone->Add("Tag_B_probe_Data2015B_DoubleEG_goldenLumi_PromptReco_50ns.root") ;

  TChain *chain_two = new TChain("t_two") ;
  chain_two->Add("Tag_B_probe_Data2015B_DoubleEG_goldenLumi_PromptReco_50ns.root") ;
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


  for(int N_bin=0;N_bin<N_bin_eta;N_bin++){
  eta_low=-3+0.2*N_bin;
  eta_high=eta_low+0.2;

  TH1F *H_moreone=new TH1F("moreone","",60,60,120);

  TH1F *H_two=new TH1F("two","",60,60,120);

    long moreone_nEntries = chain_moreone->GetEntries() ;
//    std::cout<<"moreone_Entries="<<moreone_nEntries<<std::endl;
    for(long long i=0 ; i<moreone_nEntries ; i++){
    chain_moreone->GetEntry(i) ;
//    if(i%10000==0) std::cout<<"N="<<i<<std::endl;
     

    if( moreone_Eta>eta_low && moreone_Eta<eta_high ) H_moreone->Fill(moreone_M);
                                                           
                                         } //for moreone_entries


   long two_nEntries = chain_two->GetEntries() ;
//   std::cout<<"two_Entries="<<two_nEntries<<std::endl;
   for(long long j=0 ; j<two_nEntries ; j++){
   chain_two->GetEntry(j) ;
//   if(j%10000==0) std::cout<<"N="<<j<<std::endl;

    if( two_Eta>eta_low && two_Eta<eta_high ) H_two->Fill(two_M);


  					                     } //for two_entries
	         switch(N_f){
	                       case 0:
					  moreone_N_data[N_bin]=H_moreone->Integral();
					  moreone_N_data_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_data[N_bin]=H_two->Integral();
					  two_N_data_err[N_bin]=sqrt(H_two->Integral());
					  break;
	                       case 1:
					  moreone_N_DY[N_bin]=H_moreone->Integral();
					  moreone_N_DY_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_DY[N_bin]=H_two->Integral();
					  two_N_DY_err[N_bin]=sqrt(H_two->Integral());
					  break;
	                       case 2:
					  moreone_N_TT_1[N_bin]=H_moreone->Integral();
					  moreone_N_TT_1_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_TT_1[N_bin]=H_two->Integral();
					  two_N_TT_1_err[N_bin]=sqrt(H_two->Integral());
					  break;
	                       case 3:
					  moreone_N_TT_2[N_bin]=H_moreone->Integral();
					  moreone_N_TT_2_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_TT_2[N_bin]=H_two->Integral();
					  two_N_TT_2_err[N_bin]=sqrt(H_two->Integral());
					  break;
	                       case 4:
					  moreone_N_WJet_1[N_bin]=H_moreone->Integral();
					  moreone_N_WJet_1_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_WJet_1[N_bin]=H_two->Integral();
					  two_N_WJet_1_err[N_bin]=sqrt(H_two->Integral());
					  break;
	                       case 5:
					  moreone_N_WJet_2[N_bin]=H_moreone->Integral();
					  moreone_N_WJet_2_err[N_bin]=sqrt(H_moreone->Integral());
					  two_N_WJet_2[N_bin]=H_two->Integral();
					  two_N_WJet_2_err[N_bin]=sqrt(H_two->Integral());
					  break;
                           }
  H_moreone->Delete();
  H_two->Delete();
  }//for N_bin
  
  std::cout<<"file="<<N_f<<std::endl;  
  }

 for(int L=0;L<N_bin_eta;L++){
                      moreone_data_ex[L]=moreone_N_data[L]-0.006735*(moreone_N_TT_1[L]+moreone_N_TT_2[L])-0.10259*(moreone_N_WJet_1[L]+moreone_N_WJet_2[L]);
			    moreone_data_ex_err[L]=sqrt(moreone_data_ex[L]);
                      two_data_ex[L]=two_N_data[L]-0.006735*(two_N_TT_1[L]+two_N_TT_2[L])-0.10259*(two_N_WJet_1[L]+two_N_WJet_2[L]);
			    two_data_ex_err[L]=sqrt(two_data_ex[L]);
                      moreone_MC_in[L]=0.01679*moreone_N_DY[L]+0.006735*(moreone_N_TT_1[L]+moreone_N_TT_2[L])+0.10259*(moreone_N_WJet_1[L]+moreone_N_WJet_2[L]);
			    moreone_MC_in_err[L]=sqrt(moreone_MC_in[L]);
                      two_MC_in[L]=0.01679*two_N_DY[L]+0.006735*(two_N_TT_1[L]+two_N_TT_2[L])+0.10259*(two_N_WJet_1[L]+two_N_WJet_2[L]);
			    two_MC_in_err[L]=sqrt(two_MC_in[L]);
                      moreone_MC_ex[L]=0.01679*moreone_N_DY[L]+0.01347*(moreone_N_TT_1[L]-moreone_N_TT_2[L])+0.20518*(moreone_N_WJet_1[L]-moreone_N_WJet_2[L]);
			    moreone_MC_ex_err[L]=sqrt(moreone_MC_ex[L]);
                      two_MC_ex[L]=0.01679*two_N_DY[L]+0.01347*(two_N_TT_1[L]-two_N_TT_2[L])+0.20518*(two_N_WJet_1[L]-two_N_WJet_2[L]);
			    two_MC_ex_err[L]=sqrt(two_MC_ex[L]);
                      }

//////////////////////////////////////efficiency///////////////////////////////////////////////////////
 for(int N=0;N<N_bin_eta;N++){
                    if(moreone_N_data[N]!=0){ 
			    eff_data_in[N]=two_N_data[N]/moreone_N_data[N];
			    eff_data_in_err[N]=sqrt(eff_data_in[N]*(1-eff_data_in[N])/moreone_N_data[N]);
			                            }
                    else{eff_data_in[N]=0;
			       eff_data_in_err[N]=0;
			      }

                    if(moreone_data_ex[N]!=0){ 
			    eff_data_ex[N]=two_data_ex[N]/moreone_data_ex[N];
			    eff_data_ex_err[N]=sqrt(eff_data_ex[N]*(1-eff_data_ex[N])/moreone_data_ex[N]);
			                            }
                    else{eff_data_ex[N]=0;
			       eff_data_ex_err[N]=0;
			      }
                    
			  if(moreone_MC_ex[N]!=0){ 
			    eff_MC_ex[N]=two_MC_ex[N]/moreone_MC_ex[N];
			    eff_MC_ex_err[N]=sqrt(eff_MC_ex[N]*(1-eff_MC_ex[N])/moreone_MC_ex[N]);
			                            }
                    else{eff_MC_ex[N]=0;
			       eff_MC_ex_err[N]=0;
			      }
			 
			  if(moreone_MC_in[N]!=0){ 
			    eff_MC_in[N]=two_MC_in[N]/moreone_MC_in[N];
			    eff_MC_in_err[N]=sqrt(eff_MC_in[N]*(1-eff_MC_in[N])/moreone_MC_in[N]);
			                            }
                    else{eff_MC_in[N]=0;
			       eff_MC_in_err[N]=0;
			      }

                      }


///////////////////////////////////////////////Data/MC ratio////////////////////////////////////////////////
 for(int R=0;R<N_bin_eta;R++){
                         if(eff_MC_in[R]!=0){
				    ratio_in[R]=eff_data_in[R]/eff_MC_in[R];
                            ratio_in_err[R]=sqrt(pow(eff_data_in_err[R]/eff_MC_in[R],2)+pow(eff_data_in[R]*eff_MC_in_err[R],2)/pow(eff_MC_in[R],4));
				                        }                  
                         else{ratio_in[R]=0;
				      ratio_in_err[R]=0;
				     }

                         if(eff_MC_ex[R]!=0){
				    ratio_ex[R]=eff_data_ex[R]/eff_MC_ex[R];
                            ratio_ex_err[R]=sqrt(pow(eff_data_ex_err[R]/eff_MC_ex[R],2)+pow(eff_data_ex[R]*eff_MC_ex_err[R],2)/pow(eff_MC_ex[R],4));
				                        }                  
                         else{ratio_ex[R]=0;
				      ratio_ex_err[R]=0;
				     }
                         }




///////////////////////////////////////////////////////


//////////////////////////////////////////////////Draw/////////////////////////////////////////////////////////////

float size = 0.25 ;

TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
//c1->SetFillColor(42);
c1->SetGrid();

c1->cd() ;
TPad *pad1 = new TPad("pad1", "", 0.0, size, 1.0, 1.0, 0) ;
TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, size, 0) ;
      
pad1->Draw() ;
pad2->Draw() ;

pad1->cd() ;
pad1->SetGridy() ;
pad1->SetGridx() ;
//pad1->SetLogx();

TGraphErrors *gr1 = new TGraphErrors(N_bin_eta,x,eff_MC_in,x_err,eff_MC_in_err);
gr1->SetName("gr1");
gr1->SetTitle("");
gr1->SetMarkerColor(1);
gr1->SetLineColor(2);
gr1->SetMarkerStyle(26);
gr1->SetLineWidth(2);
gr1->SetMarkerSize(1.3);
gr1->SetMaximum(1);
gr1->SetMinimum(0.5.5.5.5.5);
gr1->GetXaxis()->SetMoreLogLabels();

gr1->Draw("AP");

TGraphErrors *gr2 = new TGraphErrors(N_bin_eta,x,eff_MC_ex,x_err,eff_MC_ex_err);
gr2->SetMarkerColor(1);
gr2->SetLineColor(4);
gr2->SetName("gr2");
gr2->SetMarkerStyle(32);
gr2->SetLineWidth(2);
gr2->SetMarkerSize(1.3);
gr2->Draw("P");

TGraphErrors *gr3 = new TGraphErrors(N_bin_eta,x,eff_data_in,x_err,eff_data_in_err);
gr3->SetMarkerColor(1);
gr3->SetLineColor(1);
gr3->SetName("gr3");
gr3->SetMarkerStyle(24);
gr3->SetLineWidth(2);
gr3->SetMarkerSize(1.3);
gr3->Draw("P");

TGraphErrors *gr4 = new TGraphErrors(N_bin_eta,x,eff_data_ex,x_err,eff_data_ex_err);
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

TLatex* Label = new TLatex(0.1 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV, 50 ns, Run 2015B DoubleEG PromptReco 40.24 pb^{-1}");
Label->SetTextSize(0.05);
Label->SetNDC();
Label->Draw();

TLatex* Label1 = new TLatex(0.2 , 0.85 , "E_{T}>35GeV");
Label1->SetTextSize(0.05);
Label1->SetNDC();
Label1->Draw();

pad2->cd() ;
pad2->SetGridy() ;
pad2->SetGridx() ;
//pad2->SetLogx();

TGraphErrors *gr5 = new TGraphErrors(N_bin_eta,x,ratio_in,x_err,ratio_in_err);
gr5->SetTitle("");
gr5->SetMarkerColor(1);
gr5->SetLineColor(2);
gr5->SetName("gr5");
gr5->SetMarkerStyle(26);
gr5->SetLineWidth(2);
gr5->SetMarkerSize(1.3);
gr5->SetMinimum(0.8);
gr5->Draw("AP");

TGraphErrors *gr6 = new TGraphErrors(N_bin_eta,x,ratio_ex,x_err,ratio_ex_err);
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
