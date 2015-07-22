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

void Draw_graph(){

float MC_in_moreone_nsig[10]={7056,1609,169,58.8,12.1,903,196,18.1,4.4,0.6};
float MC_in_two_nsig[10]={5786,1383,143,50,11.3,610,142,12.1,2.8,0.37};
float MC_ex_moreone_nsig[10]={6731,1500,149,54.9,12.3,864,197,17.2,4.2,0.6};
float MC_ex_two_nsig[10]={5749,1382,139,49.9,11.3,608,141,11.3,2.7,0.33};

float MC_in_moreone_nsig_err[10]={95,44,15,9.9,4,30,19,4.6,3.2,1};
float MC_in_two_nsig_err[10]={76,42,13,10,3.7,25,12,4.5,2,0.72};
float MC_ex_moreone_nsig_err[10]={82,58,14,9.6,3.5,29,14,4.4,4.4,2.1};
float MC_ex_two_nsig_err[10]={76,42,18,9.6,3.8,25,12,5.1,1.8,0.67};

float data_in_moreone_nsig[10]={11447,1210,105,41.9,11,699,53,5,0,0};
float data_in_two_nsig[10]={10026,1067,99,34,7.5,515,37,0,0,0};
float data_ex_moreone_nsig[10]={11164,1155,109,38.7,10.6,658,28.1,3,0,0};
float data_ex_two_nsig[10]={9989,1066,108,32.8,5.9,513,25.4,0,0,0};

float data_in_moreone_nsig_err[10]={107,43,11,6.6,3.9,26,7.3,2.2,1.5,0.5};
float data_in_two_nsig_err[10]={100,47,10,5.8,3.4,23,6.1,0.5,1.2,0.5};
float data_ex_moreone_nsig_err[10]={106,41,10,6.2,3.7,26,8,1.7,1,0};
float data_ex_two_nsig_err[10]={100,47,11,5.7,2.8,23,6.7,0,0.93,0.5};

float MC_in_eff[10]=0;
float MC_ex_eff[10]=0;
float data_in_eff[10]=0;
float data_ex_eff[10]=0;
float MC_in_eff_err[10]=0;
float MC_ex_eff_err[10]=0;
float data_in_eff_err[10]=0;
float data_ex_eff_err[10]=0;

float x[5]={42.5,65,100,160,350};
float x_err[5]={7.5,15,20,40,150};

float B_MC_in_eff[5]=0;
float B_MC_ex_eff[5]=0;
float B_data_in_eff[5]=0;
float B_data_ex_eff[5]=0;
float B_MC_in_eff_err[5]=0;
float B_MC_ex_eff_err[5]=0;
float B_data_in_eff_err[5]=0;
float B_data_ex_eff_err[5]=0;

float E_MC_in_eff[5]=0;
float E_MC_ex_eff[5]=0;
float E_data_in_eff[5]=0;
float E_data_ex_eff[5]=0;
float E_MC_in_eff_err[5]=0;
float E_MC_ex_eff_err[5]=0;
float E_data_in_eff_err[5]=0;
float E_data_ex_eff_err[5]=0;

for(int i=0;i<10;i++){ 
                       if(MC_in_moreone_nsig[i]!=0) {
				                             MC_in_eff[i]=MC_in_two_nsig[i]/MC_in_moreone_nsig[i];
								     MC_in_eff_err[i]=sqrt( pow(MC_in_two_nsig_err[i]/MC_in_moreone_nsig[i],2)+pow(MC_in_two_nsig[i]*MC_in_moreone_nsig_err[i]/(pow(MC_in_moreone_nsig[i],2)),2)+pow(MC_in_two_nsig_err[i],2)*(-2*MC_in_two_nsig[i]/pow(MC_in_moreone_nsig[i],3)) );
								 //    std::cout<<MC_in_eff_err[i]<<std::endl;
			                                  }
                       else {MC_in_eff[i]=0;
				     MC_in_eff_err[i]=0; 
			          }
                       if(MC_ex_moreone_nsig[i]!=0) {
				                            MC_ex_eff[i]=MC_ex_two_nsig[i]/MC_ex_moreone_nsig[i];
								    MC_ex_eff_err[i]=sqrt( pow(MC_ex_two_nsig_err[i]/MC_ex_moreone_nsig[i],2)+pow(MC_ex_two_nsig[i]*MC_ex_moreone_nsig_err[i]/(pow(MC_ex_moreone_nsig[i],2)),2)+pow(MC_ex_two_nsig_err[i],2)*(-2*MC_ex_two_nsig[i]/pow(MC_ex_moreone_nsig[i],3)) );

			                                  } 
                       else {
				     MC_ex_eff[i]=0;
				     MC_ex_eff_err[i]=0;
			          }
                       if(data_in_moreone_nsig[i]!=0) {
				                            data_in_eff[i]=data_in_two_nsig[i]/data_in_moreone_nsig[i];
								    data_in_eff_err[i]=sqrt( pow(data_in_two_nsig_err[i]/data_in_moreone_nsig[i],2)+pow(data_in_two_nsig[i]*data_in_moreone_nsig_err[i]/(pow(data_in_moreone_nsig[i],2)),2)+pow(data_in_two_nsig_err[i],2)*(-2*data_in_two_nsig[i]/pow(data_in_moreone_nsig[i],3)) );
								      }
                       else {
				     data_in_eff[i]=0;
				     data_in_eff_err[i]=0;
			          }
                       if(data_ex_moreone_nsig[i]!=0) {
				                            data_ex_eff[i]=data_ex_two_nsig[i]/data_ex_moreone_nsig[i];
								    data_ex_eff_err[i]=sqrt( pow(data_ex_two_nsig_err[i]/data_ex_moreone_nsig[i],2)+pow(data_ex_two_nsig[i]*data_ex_moreone_nsig_err[i]/(pow(data_ex_moreone_nsig[i],2)),2)+pow(data_ex_two_nsig_err[i],2)*(-2*data_ex_two_nsig[i]/pow(data_ex_moreone_nsig[i],3)) );

								      }
                       else {
				     data_ex_eff[i]=0;
				     data_ex_eff_err[i]=0;
			          }
                     }


for(int j=0;j<5;j++){

 B_MC_in_eff[j]=MC_in_eff[j];
 B_MC_ex_eff[j]=MC_ex_eff[j];
 B_data_in_eff[j]=data_in_eff[j];
 B_data_ex_eff[j]=data_ex_eff[j];
 B_MC_in_eff_err[j]=MC_in_eff_err[j];
 B_MC_ex_eff_err[j]=MC_ex_eff_err[j];
 B_data_in_eff_err[j]=data_in_eff_err[j];
 B_data_ex_eff_err[j]=data_ex_eff_err[j];
}

for(int k=5;k<10;k++){

 E_MC_in_eff[k-5]=MC_in_eff[k];
 E_MC_ex_eff[k-5]=MC_ex_eff[k];
 E_data_in_eff[k-5]=data_in_eff[k];
 E_data_ex_eff[k-5]=data_ex_eff[k];
 E_MC_in_eff_err[k-5]=MC_in_eff_err[k];
 E_MC_ex_eff_err[k-5]=MC_ex_eff_err[k];
 E_data_in_eff_err[k-5]=data_in_eff_err[k];
 E_data_ex_eff_err[k-5]=data_ex_eff_err[k];
}

int n=5;

TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
c1->SetFillColor(42);
c1->SetGrid();


TGraphErrors *gr1 = new TGraphErrors(n,x,E_MC_in_eff,x_err,E_MC_in_eff_err);
gr1->SetName("gr1");
gr1->SetMarkerColor(2);
gr1->SetMarkerStyle(21);
gr1->SetLineWidth(1);
gr1->SetMarkerSize(1.3);
gr1->SetMaximum(1);
gr1->SetMinimum(0.);
gr1->Draw("AP");

TGraphErrors *gr2 = new TGraphErrors(n,x,E_MC_ex_eff,x_err,E_MC_ex_eff_err);
gr2->SetMarkerColor(3);
gr2->SetName("gr2");
gr2->SetMarkerStyle(22);
gr2->SetLineWidth(2);
gr2->SetMarkerSize(1.3);
gr2->Draw("P");

TGraphErrors *gr3 = new TGraphErrors(n,x,E_data_in_eff,x_err,E_data_in_eff_err);
gr3->SetMarkerColor(4);
gr3->SetName("gr3");
gr3->SetMarkerStyle(23);
gr3->SetLineWidth(3);
gr3->SetMarkerSize(1.3);
gr3->Draw("P");

TGraphErrors *gr4 = new TGraphErrors(n,x,E_data_ex_eff,x_err,E_data_ex_eff_err);
gr4->SetMarkerColor(5);
gr4->SetName("gr4");
gr4->SetMarkerStyle(20);
gr4->SetLineWidth(1);
gr4->SetMarkerSize(1.3);
gr4->Draw("P");

leg = new TLegend(0.1,0.7,0.48,0.9);

leg->AddEntry("gr1","MC inclusive","lep");
leg->AddEntry("gr2","MC exclusive","lep");
leg->AddEntry("gr3","data","lep");
leg->AddEntry("gr4","data exclusive","lep");
leg->Draw();

TLatex* Label = new TLatex(0.2 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV, 50 ns, 49.98 pb^{-1}");
Label->SetTextSize(0.05);
Label->SetNDC();
Label->Draw();

TLatex* Label1 = new TLatex(0.8 , 0.85 , "Endcap");
Label1->SetTextSize(0.05);
Label1->SetNDC();
Label1->Draw();
}
