//Fitting with the RooFit package
//Author: Chengping Shen 
void Fit()
{
  gSystem->Load("libRooFit") ;
  using namespace RooFit ;

  // --- Build Gaussian PDFs ---
  RooRealVar x("x","M_{Z} (GeV/c^{2})", 60, 120) ;

  RooRealVar mean("mean","",91,85,95) ;

  RooRealVar mean5("mean5","",89,85,95) ;

  RooRealVar width1("width1","", 4,0.3,8) ;
  RooRealVar width2("width2","", 4,0.3,8) ;
  RooRealVar width3("width3","", 4,0.3,8) ;
  RooRealVar width4("width4","", 4,0.3,8) ;
  RooRealVar width5("width5","", 5,0.2,8) ;

  RooBreitWigner bw1("bw1","Breit Wigner", x, mean, width1 );
  RooBreitWigner bw2("bw2","Breit Wigner", x, mean, width2 );
  RooBreitWigner bw3("bw3","Breit Wigner", x, mean, width3 );
  RooBreitWigner bw4("bw4","Breit Wigner", x, mean, width4 );
  RooBreitWigner bw5("bw5","Breit Wigner", x, mean5, width5 );

  // Construct px = 1 (flat in x)
  RooRealVar c1("c1","c1",0.01.,-1,1);
//  RooRealVar c2("c2","c2",1.,-5,5);
//  RooRealVar c3("c3","c3",0.,-50,50);
  RooPolynomial px("px","px",x,RooArgSet(c1)) ;
//  RooPolynomial px("px","px",x,RooArgSet(c1,c2)) ;
//  RooPolynomial px("px","px",x) ;

  RooRealVar c_e_1("c_e_1","",1.,-5,5);
  RooRealVar c_e_2("c_e_2","",1.,-5,5);

  RooPolynomial px8("px8","px",x,RooArgSet(c_e_1,c_e_2)) ;
//*** single gaussion
 RooRealVar sigmean("sigmean","x",0) ;
 RooRealVar sigwidth1("sigwidth1","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth2("sigwidth2","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth3("sigwidth3","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth4("sigwidth4","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth5("sigwidth5","x resolution", 4,0.15,10) ;

 RooGaussian gauss1("gauss1","gaussian PDF",x,sigmean,sigwidth1) ;
 RooGaussian gauss2("gauss2","gaussian PDF",x,sigmean,sigwidth2) ;
 RooGaussian gauss3("gauss3","gaussian PDF",x,sigmean,sigwidth3) ;
 RooGaussian gauss4("gauss4","gaussian PDF",x,sigmean,sigwidth4) ;
 RooGaussian gauss5("gauss5","gaussian PDF",x,sigmean,sigwidth5) ;
//*** eff curve ****
// RooGenericPdf eff("eff", "eff curve", "-0.8262+0.4445*x-0.0535*x**2", RooArgSet(x));


//Fitting pdf=Breit Wiger convolve Gaussion
// RooNumConvPdf bwxgaus("bw","bwxgaus",x,bw,gauss);
//  RooVoigtian bwxxx("bw", "bw (X) gauss", x, mean, width, sigwidth);
//x.setBins("fft",10000);
RooFFTConvPdf bxg1("bxg1","",x,bw1,gauss1); 
RooFFTConvPdf bxg2("bxg2","",x,bw2,gauss2);
RooFFTConvPdf bxg3("bxg3","",x,bw3,gauss3);
RooFFTConvPdf bxg4("bxg4","",x,bw4,gauss4);
RooFFTConvPdf bxg5("bxg5","",x,bw5,gauss5);
//*** BW * eff curve *****
// RooProdPdf bwtmp("bwtmp", "bwtmp", RooArgSet(bwxxx, eff));

// RooRealVar nsig("nsig","number of signal events",Nsig) ;
 RooRealVar nsig1("nsig1","number of signal events",4000,1,7000) ;
 RooRealVar nsig2("nsig2","number of signal events",500,1,7000) ;
 RooRealVar nsig3("nsig3","number of signal events",10,1,7000) ;
 RooRealVar nsig4("nsig4","number of signal events",100,1,7000) ;
 RooRealVar nsig5("nsig5","number of signal events",100,1,7000) ;

// RooRealVar nsig("nsig","number of signal events",0) ;
 RooRealVar nbkg1("nbkg1","number of background events",1,1,100) ;
 RooRealVar nbkg2("nbkg2","number of background events",1,1,1000) ;
 RooRealVar nbkg3("nbkg3","number of background events",1,1,100) ;
 RooRealVar nbkg4("nbkg4","number of background events",1,1,100) ;
 RooRealVar nbkg5("nbkg5","number of background events",1,1,100) ;

 RooAddPdf model1("sum1","bkg+bw",RooArgList(bxg1,px),RooArgList(nsig1,nbkg1)) ;
 RooAddPdf model2("sum2","bkg+bw",RooArgList(bxg2,px),RooArgList(nsig2,nbkg2)) ;
 RooAddPdf model3("sum3","bkg+bw",RooArgList(bxg3,px),RooArgList(nsig3,nbkg3)) ;
 RooAddPdf model4("sum4","bkg+bw",RooArgList(bxg4,px),RooArgList(nsig4,nbkg4)) ;
 RooAddPdf model5("sum5","bkg+bw",RooArgList(bxg5,px),RooArgList(nsig5,nbkg5)) ;
 



// RooAddPdf model("sum","bw",RooArgList(bwxxx),RooArgList(nsig)) ;
//  RooRealVar nsig("f","f",0.1,0.,1.) ;
//  RooAddPdf model("model","model",RooArgList(bwtmp,px),f) ;//
//




//
  TFile *f1 = new TFile("./Tag_probe_0703.root");
  TH1F *h_two_35_50=(TH1F*)f1->Get("two_35_50");
  TH1F *h_two_50_100=(TH1F*)f1->Get("two_50_100");
  TH1F *h_two_100_5000=(TH1F*)f1->Get("two_100_5000");
  TH1F *h_two_barrel=(TH1F*)f1->Get("two_barrel");
  TH1F *h_two_endcap=(TH1F*)f1->Get("two_endcap");

  RooDataHist H_two_35_50("H_two_35_50","",x,h_two_35_50);
  RooDataHist H_two_50_100("H_two_50_100","",x,h_two_50_100);
  RooDataHist H_two_100_5000("H_two_100_5000","",x,h_two_100_5000); 
  RooDataHist H_two_barrel("H_two_barrel","",x,h_two_barrel);
  RooDataHist H_two_endcap("H_two_endcap","",x,h_two_endcap);



  //// Fix some parameters 
   sigmean.setConstant(kTRUE);
//   sigwidth.setConstant(kTRUE);
//   mean.setConstant(kTRUE);
//   width.setConstant(kTRUE);  
//
   
//   nsig.setConstant(kTRUE);
//   rr = (RooFitResult*) model.fitTo(hmc,Extended(),Minos(kFALSE),Save(kTRUE));
    RooFitResult* m_35_50 = model1.fitTo(H_two_35_50,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooFitResult* m_50_100 = model2.fitTo(H_two_50_100,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooFitResult* m_100_5000 = model3.fitTo(H_two_100_5000,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooFitResult* m_barrel = model4.fitTo(H_two_barrel,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooFitResult* m_endcap = model5.fitTo(H_two_endcap,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    //  sum1.fitTo(hdata,Extended()) ;
      // --- Plot toy data and composite PDF overlaid ---
   


//    RooPlot* xframe = x.frame();



	   // set Margins for labels etc.
/*    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadBottomMargin(0.17);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
*/
  
    RooPlot* xframe_m_35_50 = x.frame();
    H_two_35_50.plotOn(xframe_m_35_50) ;
    model1.plotOn(xframe_m_35_50) ;

    RooPlot* xframe_m_50_100 = x.frame();
    H_two_50_100.plotOn(xframe_m_50_100) ;
    model2.plotOn(xframe_m_50_100) ;

      
    RooPlot* xframe_m_100_5000 = x.frame();
    H_two_100_5000.plotOn(xframe_m_100_5000) ;
    model3.plotOn(xframe_m_100_5000) ;



    RooPlot* xframe_m_barrel = x.frame();
    H_two_barrel.plotOn(xframe_m_barrel) ;
    model4.plotOn(xframe_m_barrel) ;

    RooPlot* xframe_m_endcap = x.frame();
    H_two_endcap.plotOn(xframe_m_endcap) ;
    model5.plotOn(xframe_m_endcap) ;
//    px8.plotOn(xframe_m_endcap) ;





  // --- Plot frame on canvas ---
/*    xframe->SetTitle("");
    xframe->SetYTitle("Number / (GeV/c^{2})");
    xframe->GetXaxis()->CenterTitle(true);
    xframe->GetXaxis()->SetTitleSize(0.08);
    xframe->GetXaxis()->SetLabelSize(0.06);
    xframe->SetXTitle("M^{Z}(GeV/c^{2})");
    xframe->GetYaxis()->CenterTitle(true);
    xframe->GetYaxis()->SetTitleSize(0.08);
    xframe->GetYaxis()->SetLabelSize(0.06);
    xframe->SetNdivisions(510,"X");
*/
/*
    TCanvas* c = new TCanvas("b0-jp-mc","b0-jp-mc",600,400) ;
  
    TLatex* N_lab = new TLatex(0.2 , 0.85 , "N_{signal}=5600.04#pm74.82");

    N_lab->SetTextSize(0.06);
    N_lab->SetNDC();

   TLatex* Label = new TLatex(0.2 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV 25 ns, L = 14.37 pb^{-1}");

    Label->SetTextSize(0.05);
    Label->SetNDC();


    xframe->addObject(N_lab);
    xframe->addObject(Label);
    xframe->Draw() ;
*/

    TLatex* Label = new TLatex(0.2 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV 25 ns, L = 14.37 pb^{-1}");
    Label->SetTextSize(0.05);
    Label->SetNDC();

   TCanvas* ca1 = new TCanvas("ca1","ca1",600,400) ;
//   TLatex* N_1 = new TLatex(0.2 , 0.85 , "N_{signal}=%.2f#pm%.2f",&nsig1.getVal(),&nsig1.getError());
   TLatex* N_1 = new TLatex(0.15 , 0.8 , "E_{T}^{probe}(35-50GeV)");
   N_1->SetTextSize(0.05);
   N_1->SetNDC();
   TLatex* L_1 = new TLatex(0.15 , 0.7 , "N=4671 #pm 70");
   L_1->SetTextSize(0.05);
   L_1->SetNDC();
   xframe_m_35_50->addObject(L_1);
   xframe_m_35_50->addObject(N_1);
   xframe_m_35_50->addObject(Label);
   xframe_m_35_50->Draw();


   TCanvas* ca2 = new TCanvas("ca2","ca2",600,400) ;
   TLatex* N_2 = new TLatex(0.15 , 0.8 , "E_{T}^{probe}(50-100GeV)");
   N_2->SetTextSize(0.06);
   N_2->SetNDC();
   TLatex* L_2 = new TLatex(0.15 , 0.7 , " N=892 #pm 58");
   L_2->SetTextSize(0.06);
   L_2->SetNDC();
   xframe_m_50_100->addObject(L_2);
   xframe_m_50_100->addObject(N_2);
   xframe_m_50_100->addObject(Label);
   xframe_m_50_100->Draw();

   TCanvas* ca3 = new TCanvas("ca3","ca3",600,400) ;
   TLatex* N_3 = new TLatex(0.15 , 0.8 , "E_{T}^{probe}(100-5000GeV)");
   N_3->SetTextSize(0.06);
   N_3->SetNDC();
   TLatex* L_3 = new TLatex(0.15 , 0.7 , "N=37 #pm 11");
   L_3->SetTextSize(0.06);
   L_3->SetNDC();
   xframe_m_100_5000->addObject(L_3);
   xframe_m_100_5000->addObject(N_3);
   xframe_m_100_5000->addObject(Label);
   xframe_m_100_5000->Draw();
  
   TCanvas* ca4 = new TCanvas("ca4","ca4",600,400) ;
   TLatex* N_4 = new TLatex(0.15 , 0.8 , "probe in barrel");
   N_4->SetTextSize(0.06);
   N_4->SetNDC();
   TLatex* L_4 = new TLatex(0.15 , 0.7 , "N=4104 #pm 72");
   L_4->SetTextSize(0.06);
   L_4->SetNDC();
   xframe_m_barrel->addObject(L_4);
   xframe_m_barrel->addObject(N_4);
   xframe_m_barrel->addObject(Label);
   xframe_m_barrel->Draw();

   TCanvas* ca5 = new TCanvas("ca5","ca5",600,400) ;
   TLatex* N_5 = new TLatex(0.15 , 0.8 , "probe in endcap");
   N_5->SetTextSize(0.06);
   N_5->SetNDC();
   TLatex* L_5 = new TLatex(0.15 , 0.7 , "N=1495 #pm 49");
   L_5->SetTextSize(0.06);
   L_5->SetNDC();
   xframe_m_endcap->addObject(L_5);
   xframe_m_endcap->addObject(N_5);
   xframe_m_endcap->addObject(Label);
   xframe_m_endcap->Draw();



   std::cout<<"N_m_35_50="<<nsig1.getVal()<<" Error="<<nsig1.getError()<<"bkg="<<nbkg1.getVal()<<"error="<<nbkg1.getError()<<std::endl;
   std::cout<<"N_m_50_100="<<nsig2.getVal()<<" Error="<<nsig2.getError()<<"bkg="<<nbkg2.getVal()<<"error="<<nbkg2.getError()<<std::endl;
   std::cout<<"N_m_100_5000="<<nsig3.getVal()<<" Error="<<nsig3.getError()<<"bkg="<<nbkg3.getVal()<<"error="<<nbkg3.getError()<<std::endl;
   std::cout<<"N_m_barrel="<<nsig4.getVal()<<" Error="<<nsig4.getError()<<"bkg="<<nbkg4.getVal()<<"error="<<nbkg4.getError()<<std::endl;
   std::cout<<"N_m_endcap="<<nsig5.getVal()<<" Error="<<nsig5.getError()<<"bkg="<<nbkg5.getVal()<<"error="<<nbkg5.getError()<<std::endl;

   /*
   
N_m_35_50=5478.21 Error=74.0124bkg=1error=2.50334
N_m_50_100=894.343 Error=37.0291bkg=112.65error=24.2647
N_m_100_5000=41.9847 Error=6.56334bkg=1.0076error=10.368
N_m_barrel=4713.69 Error=68.6724bkg=1.09938error=6.9623
N_m_endcap=1814.02 Error=42.6083bkg=1.00055error=17.8272


N_m_35_50=4670.23 Error=68.34bkg=1error=2.13229
N_m_50_100=783.263 Error=34.7269bkg=108.875error=23.0412
N_m_100_5000=35.1818 Error=6.68252bkg=1.78685error=3.50993
N_m_barrel=4103.02 Error=64.0624bkg=1error=8.27622
N_m_endcap=1494.96 Error=38.6774bkg=1.00013error=10.4842
   */
}


