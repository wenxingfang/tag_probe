//Fitting with the RooFit package
//Author: Chengping Shen 
void Fit()
{
  gSystem->Load("libRooFit") ;
  using namespace RooFit ;

  RooRealVar x("x","M_{ee} (GeV/c^{2})", 60, 120) ;

  //BreitWinger
  
  RooRealVar mean("mean","",91,80,95) ;

  RooRealVar mean2("mean2","",91.1,80.1,95) ;

  RooRealVar mean5("mean5","",91,80,95) ;

  RooRealVar width("width","", 4,0.3,8) ;
  RooRealVar width1("width1","", 4,0.3,8) ;
  RooRealVar width2("width2","", 2.5,0.3,8) ;
  RooRealVar width3("width3","", 4,0.3,8) ;
  RooRealVar width4("width4","", 4,0.1,8) ;
  RooRealVar width5("width5","", 6,0.2,8) ;

  RooBreitWigner bw1("bw1","Breit Wigner", x, mean, width1 );
  RooBreitWigner bw2("bw2","Breit Wigner", x, mean, width );
  RooBreitWigner bw3("bw3","Breit Wigner", x, mean, width3 );
  RooBreitWigner bw4("bw4","Breit Wigner", x, mean, width4 );
  RooBreitWigner bw5("bw5","Breit Wigner", x, mean5, width5 );

  //Novosibirsk

  RooRealVar meanNb("meanNb","meanNb", 0,-1,1) ;
  RooRealVar widthNb("widthNb","widthNb", 0.1, 0.001, 5) ;
  RooRealVar tailNb("tailNb","tailNb", 0.3, -5, 4) ;
  RooNovosibirsk Nb("Nb","Nb",x,meanNb,widthNb,tailNb);

  //Landau
  RooRealVar meanLa("meanLa","meanLa", 0) ;
  RooRealVar sigmaLa2("sigmaLa","x resolution", 1.0,0.15,10) ;
  RooLandau La2("La2","La2",x,meanLa,sigmaLa2);

  RooRealVar sigmaLa5("sigmaLa5","x resolution", 2,0.01,10) ;
  RooLandau La5("La5","La5",x,meanLa,sigmaLa5);

  //Polynomial

RooRealVar c1("c1","c1",0.001,-0.5,0.01);
RooRealVar c2("c2","c2",0.001,-1.5,0.01);
//  RooRealVar c3("c3","c3",0.,-50,50);
RooPolynomial px("px","px",x,RooArgSet(c1,c2)) ;
  //Chebychev
RooChebychev Cheb("Cheb","Cheb",x,RooArgSet(c1,c2));
//RooCBShape

  RooRealVar CB_m0("CB_m0","x",0) ;
  RooRealVar CB_sigma("CB_sigma","x resolution", 2,0.15,10) ;
  RooRealVar CB_alpha("CB_alpha","x resolution", 2,-10,10) ;
  RooRealVar CB_n("CB_n","x resolution", 2,1,130) ;
  
  RooCBShape CB("CB","CB",x,CB_m0,CB_sigma,CB_alpha,CB_n);
   //gaussion

 RooRealVar sigmean("sigmean","x",0) ;
 RooRealVar sigwidth1("sigwidth1","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth2("sigwidth2","x resolution", 3,0.15,10) ;
 RooRealVar sigwidth3("sigwidth3","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth4("sigwidth4","x resolution", 2,0.15,10) ;
 RooRealVar sigwidth5("sigwidth5","x resolution", 4,0.15,10) ;

 RooGaussian gauss1("gauss1","gaussian PDF",x,sigmean,sigwidth1) ;
 RooGaussian gauss2("gauss2","gaussian PDF",x,sigmean,sigwidth2) ;
 RooGaussian gauss3("gauss3","gaussian PDF",x,sigmean,sigwidth3) ;
 RooGaussian gauss4("gauss4","gaussian PDF",x,sigmean,sigwidth4) ;
 RooGaussian gauss5("gauss5","gaussian PDF",x,sigmean,sigwidth5) ;

 //Double Gaussin

 RooRealVar mean1("mean1","#eta_{c} mass", 0.0) ;
 RooRealVar sigma1("sigma1", "eta_{c} resolution1",0.07 , 0.001,1.);
 RooRealVar ratio2("ratio2", "ratio of the second gaussion", 7.e-01,0.01,1.);
 RooRealVar diff("diff", "mean difference", 0.02 , -10.0, 10.0) ;
 RooRealVar sigratio2("sigratio2", "the second gaussion sigma", 4.0e-01, 0.0, 10.0);

 RooGenericPdf dgaus("dgaus", "double gaussion", "(1.-ratio2)*(1./sqrt(2.*3.1415926))*exp(-0.5*((x-mean1)/sigma1)**2)/sigma1 + ratio2*(1./sqrt(2.*3.1415926))*exp(-0.5*((x-mean1-diff)/(sigma1*sigratio2))**2)/(sigma1*sigratio2)", RooArgSet(x,mean1,sigma1,ratio2, diff, sigratio2));


//*** eff curve ****
// RooGenericPdf eff("eff", "eff curve", "-0.8262+0.4445*x-0.0535*x**2", RooArgSet(x));


//Fitting pdf=Breit Wiger convolve Gaussion
// RooNumConvPdf bwxgaus("bw","bwxgaus",x,bw,gauss);
//  RooVoigtian bwxxx("bw", "bw (X) gauss", x, mean, width, sigwidth);
//x.setBins("fft",10000);

 //bwXgauss
RooFFTConvPdf bxg1("bxg1","",x,bw1,gauss1); 
RooFFTConvPdf bxg2("bxg2","",x,bw2,gauss2);
RooFFTConvPdf bxg3("bxg3","",x,bw3,gauss3);
RooFFTConvPdf bxg4("bxg4","",x,bw4,gauss4);
RooFFTConvPdf bxg5("bxg5","",x,bw5,gauss5);

 //bwXNovosibirsk

RooFFTConvPdf bxN2("bxN2","",x,bw2,Nb);

 //bwXLandau

RooFFTConvPdf bxLa2("bxLa2","",x,bw2,La2);
RooFFTConvPdf bxLa5("bxLa5","",x,bw2,La5);


 //bwXCB
RooFFTConvPdf bxCB("bxCB","",x,bw1,CB);


//*** BW * eff curve *****
// RooProdPdf bwtmp("bwtmp", "bwtmp", RooArgSet(bwxxx, eff));

// RooRealVar nsig("nsig","number of signal events",Nsig) ;
 RooRealVar nsig("nsig","number of signal",10000,0,200000) ;
 RooRealVar nsig1("nsig1","number of signal",10000,1,200000) ;
 RooRealVar nsig2("nsig2","number of signal",10000,1,200000) ;
 RooRealVar nsig3("nsig3","number of signal",10000,1,170000) ;
 RooRealVar nsig4("nsig4","number of signal",10000,1,200000) ;
 RooRealVar nsig5("nsig5","number of signal",10000,1,70000) ;
 RooRealVar nsig6("nsig6","number of signal",10000,1,200000) ;
 RooRealVar nsig7("nsig7","number of signal",10000,1,200000) ;
 RooRealVar nsig8("nsig8","number of signal",10000,1,170000) ;
 RooRealVar nsig9("nsig9","number of signal",10000,1,200000) ;
 RooRealVar nsig10("nsig10","number of signal",10000,1,70000) ;

// RooRealVar nsig("nsig","number of signal events",0) ;
 RooRealVar nbkg("nbkg","number of background events",1,0,10000) ;
 RooRealVar nbkg1("nbkg1","number of background events",1,0,100) ;
 RooRealVar nbkg2("nbkg2","number of background events",1,1,10000) ;
 RooRealVar nbkg3("nbkg3","number of background events",1,0,1000) ;
 RooRealVar nbkg4("nbkg4","number of background events",1,1,500) ;
 RooRealVar nbkg5("nbkg5","number of background events",1,0,100) ;

 RooAddPdf model1("sum1","bkg+bw",RooArgList(bxCB,Cheb),RooArgList(nsig,nbkg)) ;
 RooAddPdf model2("sum2","bkg+bw",RooArgList(bxLa2,Cheb),RooArgList(nsig,nbkg)) ;
 RooAddPdf model3("sum3","bkg+bw",RooArgList(bxg3,Cheb),RooArgList(nsig,nbkg)) ;
 RooAddPdf model4("sum4","bkg+bw",RooArgList(bxg4,Cheb),RooArgList(nsig,nbkg)) ;
 RooAddPdf model5("sum5","bkg+bw",RooArgList(bxLa5,Cheb),RooArgList(nsig,nbkg)) ;

/* 
 RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.9,0.,1.) ;
 RooAddPdf bxCB_Cheb("bxCB+cheb","bxCB+cheb",RooArgList(bxCB,Cheb),sig1frac);
 RooAddPdf bxLa2_Cheb("bxLa2+cheb","bxLa2+cheb",RooArgList(bxLa2,Cheb),sig1frac);
 RooAddPdf bxg3_Cheb("bxg3+cheb","bxg3+cheb",RooArgList(bxg3,Cheb),sig1frac);
 RooAddPdf bxg4_Cheb("bxg4+cheb","bxg4+cheb",RooArgList(bxg4,Cheb),sig1frac);
 RooAddPdf bxLa5_Cheb("bxLa5+cheb","bxLa5+cheb",RooArgList(bxLa5,Cheb),sig1frac);

 RooExtendPdf model1("esig","extended signal p.d.f",bxCB_Cheb,n1) ;
 RooExtendPdf model2("esig","extended signal p.d.f",bxLa2_Cheb,n2) ;
 RooExtendPdf model3("esig","extended signal p.d.f",bxg3_Cheb,n3) ;
 RooExtendPdf model4("esig","extended signal p.d.f",bxg4_Cheb,n4) ;
 RooExtendPdf model5("esig","extended signal p.d.f",bxLa5_Cheb,n5) ;

 RooExtendPdf model6("esig","extended signal p.d.f",bxLa2_Cheb,n6) ;
 */
// RooAddPdf model("sum","bw",RooArgList(bwxxx),RooArgList(nsig)) ;
//  RooRealVar nsig("f","f",0.1,0.,1.) ;
//  RooAddPdf model("model","model",RooArgList(bwtmp,px),f) ;//
//




//
  TFile *f1 = new TFile("./Tag_probe_Data2015B_hist.root");
  TH1F *h_moreone_35_50_B=(TH1F*)f1->Get("moreone_35_50_B");
  TH1F *h_moreone_50_80_B=(TH1F*)f1->Get("moreone_50_80_B");
  TH1F *h_moreone_80_120_B=(TH1F*)f1->Get("moreone_80_120_B");
  TH1F *h_moreone_120_200_B=(TH1F*)f1->Get("moreone_120_200_B");
  TH1F *h_moreone_200_plus_B=(TH1F*)f1->Get("moreone_200_plus_B");
  TH1F *h_moreone_35_50_E=(TH1F*)f1->Get("moreone_35_50_E");
  TH1F *h_moreone_50_80_E=(TH1F*)f1->Get("moreone_50_80_E");
  TH1F *h_moreone_80_120_E=(TH1F*)f1->Get("moreone_80_120_E");
  TH1F *h_moreone_120_200_E=(TH1F*)f1->Get("moreone_120_200_E");
  TH1F *h_moreone_200_plus_E=(TH1F*)f1->Get("moreone_200_plus_E");


  RooDataHist H_moreone_35_50_B("H_moreone_35_50_B","",x,h_moreone_35_50_B);
  RooDataHist H_moreone_50_80_B("H_moreone_50_80_B","",x,h_moreone_50_80_B);
  RooDataHist H_moreone_80_120_B("H_moreone_80_120_B","",x,h_moreone_80_120_B); 
  RooDataHist H_moreone_120_200_B("H_moreone_120_200_B","",x,h_moreone_120_200_B);
  RooDataHist H_moreone_200_plus_B("H_moreone_200_plus_B","",x,h_moreone_200_plus_B);
  RooDataHist H_moreone_35_50_E("H_moreone_35_50_E","",x,h_moreone_35_50_E);
  RooDataHist H_moreone_50_80_E("H_moreone_50_80_E","",x,h_moreone_50_80_E);
  RooDataHist H_moreone_80_120_E("H_moreone_80_120_E","",x,h_moreone_80_120_E); 
  RooDataHist H_moreone_120_200_E("H_moreone_120_200_E","",x,h_moreone_120_200_E);
  RooDataHist H_moreone_200_plus_E("H_moreone_200_plus_E","",x,h_moreone_200_plus_E);


  //// Fix some parameters 
   sigmean.setConstant(kTRUE);
//   sigwidth.setConstant(kTRUE);
   
//   rr = (RooFitResult*) model.fitTo(hmc,Extended(),Minos(kFALSE),Save(kTRUE));


    //  sum1.fitTo(hdata,Extended()) ;
      // --- Plot toy data and composite PDF overlaid ---
	   // set Margins for labels etc.
/*    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadBottomMargin(0.17);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
*/

    TLatex* Label = new TLatex(0.2 , 0.95 , "CMS Interal, #sqrt{s} = 13 TeV, 50 ns, 49.98 pb^{-1}");
    Label->SetTextSize(0.05);
    Label->SetNDC();

    TLatex* Label_BWXLa = new TLatex(0.15 , 0.7 , "BW#otimesLandau + Chebychev");
    Label_BWXLa->SetTextSize(0.05);
    Label_BWXLa->SetNDC();
    
    TLatex* Label_BWXCB = new TLatex(0.15 , 0.7 , "BW#otimesCB + Chebychev");
    Label_BWXCB->SetTextSize(0.05);
    Label_BWXCB->SetNDC();

    TLatex* Label_BWXGauss = new TLatex(0.15 , 0.7 , "BW#otimesGauss + Chebychev");
    Label_BWXGauss->SetTextSize(0.05);
    Label_BWXGauss->SetNDC();

/*
    RooFitResult* m_35_50_B = model1.fitTo(H_moreone_35_50_B,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_35_50_B = x.frame();
    H_moreone_35_50_B.plotOn(xframe_m_35_50_B,DataError(RooAbsData::SumW2)) ;
    model1.plotOn(xframe_m_35_50_B) ;
    model1.plotOn(xframe_m_35_50_B,Components(Cheb),LineStyle(kDashed),LineColor(kRed));
    model1.paramOn(xframe_m_35_50_B,Layout(0.55)) ;

   TCanvas* ca1 = new TCanvas("ca1","ca1",600,400) ;
   TLatex* N_1 = new TLatex(0.15 , 0.8 , "Barrel, E_{T}^{probe}(35-50GeV)");
   N_1->SetTextSize(0.05);
   N_1->SetNDC();
   xframe_m_35_50_B->addObject(N_1);
   xframe_m_35_50_B->addObject(Label);
   xframe_m_35_50_B->addObject(Label_BWXCB);
   xframe_m_35_50_B->SetTitle("");
   xframe_m_35_50_B->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_35_50_B->GetXaxis()->CenterTitle(true);
   xframe_m_35_50_B->GetXaxis()->SetTitleSize(0.08);
   xframe_m_35_50_B->GetXaxis()->SetLabelSize(0.06);
   xframe_m_35_50_B->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_35_50_B->GetYaxis()->CenterTitle(true);
   xframe_m_35_50_B->GetYaxis()->SetTitleSize(0.08);
   xframe_m_35_50_B->GetYaxis()->SetLabelSize(0.06);
   xframe_m_35_50_B->SetTitleOffset(0.7, "Y");
   xframe_m_35_50_B->SetTitleOffset(0.8, "X");
   xframe_m_35_50_B->Draw();

    RooFitResult* m_50_80_B = model2.fitTo(H_moreone_50_80_B,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_50_80_B = x.frame();
    H_moreone_50_80_B.plotOn(xframe_m_50_80_B,DataError(RooAbsData::SumW2)) ;
    model2.plotOn(xframe_m_50_80_B) ;
    model2.plotOn(xframe_m_50_80_B,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model2.paramOn(xframe_m_50_80_B,Layout(0.55)) ;

   TCanvas* ca2 = new TCanvas("ca2","ca2",600,400) ;
   TLatex* N_2 = new TLatex(0.15 , 0.8 , "Barrel, E_{T}^{probe}(50-80GeV)");
   N_2->SetTextSize(0.06);
   N_2->SetNDC();
   xframe_m_50_80_B->addObject(N_2);
   xframe_m_50_80_B->addObject(Label);
   xframe_m_50_80_B->addObject(Label_BWXLa);
   xframe_m_50_80_B->SetTitle("");
   xframe_m_50_80_B->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_50_80_B->GetXaxis()->CenterTitle(true);
   xframe_m_50_80_B->GetXaxis()->SetTitleSize(0.08);
   xframe_m_50_80_B->GetXaxis()->SetLabelSize(0.06);
   xframe_m_50_80_B->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_50_80_B->GetYaxis()->CenterTitle(true);
   xframe_m_50_80_B->GetYaxis()->SetTitleSize(0.08);
   xframe_m_50_80_B->GetYaxis()->SetLabelSize(0.06);
   xframe_m_50_80_B->SetTitleOffset(0.7, "Y");
   xframe_m_50_80_B->SetTitleOffset(0.8, "X");
   xframe_m_50_80_B->Draw();
      
    RooFitResult* m_80_120_B = model3.fitTo(H_moreone_80_120_B,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_80_120_B = x.frame();
    H_moreone_80_120_B.plotOn(xframe_m_80_120_B,DataError(RooAbsData::SumW2)) ;
    model3.plotOn(xframe_m_80_120_B) ;
    model3.plotOn(xframe_m_80_120_B,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model3.paramOn(xframe_m_80_120_B,Layout(0.55)) ;

   TCanvas* ca3 = new TCanvas("ca3","ca3",600,400) ;
   TLatex* N_3 = new TLatex(0.15 , 0.8 , "Barrel, E_{T}^{probe}(80-120GeV)");
   N_3->SetTextSize(0.06);
   N_3->SetNDC();
   xframe_m_80_120_B->addObject(N_3);
   xframe_m_80_120_B->addObject(Label);
   xframe_m_80_120_B->addObject(Label_BWXGauss);
   xframe_m_80_120_B->SetTitle("");
   xframe_m_80_120_B->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_80_120_B->GetXaxis()->CenterTitle(true);
   xframe_m_80_120_B->GetXaxis()->SetTitleSize(0.08);
   xframe_m_80_120_B->GetXaxis()->SetLabelSize(0.06);
   xframe_m_80_120_B->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_80_120_B->GetYaxis()->CenterTitle(true);
   xframe_m_80_120_B->GetYaxis()->SetTitleSize(0.08);
   xframe_m_80_120_B->GetYaxis()->SetLabelSize(0.06);
   xframe_m_80_120_B->SetTitleOffset(0.7, "Y");
   xframe_m_80_120_B->SetTitleOffset(0.8, "X");
   xframe_m_80_120_B->Draw();

    RooFitResult* m_120_200_B = model4.fitTo(H_moreone_120_200_B,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_120_200_B = x.frame();
    H_moreone_120_200_B.plotOn(xframe_m_120_200_B,DataError(RooAbsData::SumW2)) ;
    model4.plotOn(xframe_m_120_200_B) ;
    model4.plotOn(xframe_m_120_200_B,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model4.paramOn(xframe_m_120_200_B,Layout(0.55)) ;
    
   TCanvas* ca4 = new TCanvas("ca4","ca4",600,400) ;
   TLatex* N_4 = new TLatex(0.15 , 0.8 , "Barrel, E_{T}^{probe}(120-200GeV)");
   N_4->SetTextSize(0.06);
   N_4->SetNDC();
   xframe_m_120_200_B->addObject(N_4);
   xframe_m_120_200_B->addObject(Label);
   xframe_m_120_200_B->addObject(Label_BWXGauss);
   xframe_m_120_200_B->SetTitle("");
   xframe_m_120_200_B->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_120_200_B->GetXaxis()->CenterTitle(true);
   xframe_m_120_200_B->GetXaxis()->SetTitleSize(0.08);
   xframe_m_120_200_B->GetXaxis()->SetLabelSize(0.06);
   xframe_m_120_200_B->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_120_200_B->GetYaxis()->CenterTitle(true);
   xframe_m_120_200_B->GetYaxis()->SetTitleSize(0.08);
   xframe_m_120_200_B->GetYaxis()->SetLabelSize(0.06);
   xframe_m_120_200_B->SetTitleOffset(0.7, "Y");
   xframe_m_120_200_B->SetTitleOffset(0.8, "X");
   xframe_m_120_200_B->Draw();

    RooFitResult* m_200_plus_B = model5.fitTo(H_moreone_200_plus_B,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_200_plus_B = x.frame();
    H_moreone_200_plus_B.plotOn(xframe_m_200_plus_B,DataError(RooAbsData::SumW2)) ;
    model5.plotOn(xframe_m_200_plus_B) ;
    model5.plotOn(xframe_m_200_plus_B,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model5.paramOn(xframe_m_200_plus_B,Layout(0.55)) ;

   TCanvas* ca5 = new TCanvas("ca5","ca5",600,400) ;
   TLatex* N_5 = new TLatex(0.15 , 0.8 , "Barrel, E_{T}^{probe}(>200GeV)");
   N_5->SetTextSize(0.06);
   N_5->SetNDC();
   xframe_m_200_plus_B->addObject(N_5);
   xframe_m_200_plus_B->addObject(Label);
   xframe_m_200_plus_B->addObject(Label_BWXLa);
   xframe_m_200_plus_B->SetTitle("");
   xframe_m_200_plus_B->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_200_plus_B->GetXaxis()->CenterTitle(true);
   xframe_m_200_plus_B->GetXaxis()->SetTitleSize(0.08);
   xframe_m_200_plus_B->GetXaxis()->SetLabelSize(0.06);
   xframe_m_200_plus_B->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_200_plus_B->GetYaxis()->CenterTitle(true);
   xframe_m_200_plus_B->GetYaxis()->SetTitleSize(0.08);
   xframe_m_200_plus_B->GetYaxis()->SetLabelSize(0.06);
   xframe_m_200_plus_B->SetTitleOffset(0.7, "Y");
   xframe_m_200_plus_B->SetTitleOffset(0.8, "X");
   xframe_m_200_plus_B->Draw();

    RooFitResult* m_35_50_E = model2.fitTo(H_moreone_35_50_E,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_35_50_E = x.frame();
    H_moreone_35_50_E.plotOn(xframe_m_35_50_E,DataError(RooAbsData::SumW2)) ;
    model2.plotOn(xframe_m_35_50_E) ;
    model2.plotOn(xframe_m_35_50_E,Components(Cheb),LineStyle(kDashed),LineColor(kRed));
    model2.paramOn(xframe_m_35_50_E,Layout(0.55)) ;

   TCanvas* ca6 = new TCanvas("ca6","ca6",600,400) ;
   TLatex* N_6 = new TLatex(0.15 , 0.8 , "Earrel, E_{T}^{probe}(35-50GeV)");
   N_6->SetTextSize(0.05);
   N_6->SetNDC();
   xframe_m_35_50_E->addObject(N_6);
   xframe_m_35_50_E->addObject(Label);
   xframe_m_35_50_E->addObject(Label_BWXLa);
   xframe_m_35_50_E->SetTitle("");
   xframe_m_35_50_E->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_35_50_E->GetXaxis()->CenterTitle(true);
   xframe_m_35_50_E->GetXaxis()->SetTitleSize(0.08);
   xframe_m_35_50_E->GetXaxis()->SetLabelSize(0.06);
   xframe_m_35_50_E->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_35_50_E->GetYaxis()->CenterTitle(true);
   xframe_m_35_50_E->GetYaxis()->SetTitleSize(0.08);
   xframe_m_35_50_E->GetYaxis()->SetLabelSize(0.06);
   xframe_m_35_50_E->SetTitleOffset(0.7, "Y");
   xframe_m_35_50_E->SetTitleOffset(0.8, "X");
   xframe_m_35_50_E->Draw();

    RooFitResult* m_50_80_E = model2.fitTo(H_moreone_50_80_E,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_50_80_E = x.frame();
    H_moreone_50_80_E.plotOn(xframe_m_50_80_E,DataError(RooAbsData::SumW2)) ;
    model2.plotOn(xframe_m_50_80_E) ;
    model2.plotOn(xframe_m_50_80_E,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model2.paramOn(xframe_m_50_80_E,Layout(0.55)) ;
      
   TCanvas* ca7 = new TCanvas("ca7","ca7",600,400) ;
   TLatex* N_7 = new TLatex(0.15 , 0.8 , "Earrel, E_{T}^{probe}(50-80GeV)");
   N_7->SetTextSize(0.06);
   N_7->SetNDC();
   xframe_m_50_80_E->addObject(N_7);
   xframe_m_50_80_E->addObject(Label);
   xframe_m_50_80_E->addObject(Label_BWXLa);
   xframe_m_50_80_E->SetTitle("");
   xframe_m_50_80_E->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_50_80_E->GetXaxis()->CenterTitle(true);
   xframe_m_50_80_E->GetXaxis()->SetTitleSize(0.08);
   xframe_m_50_80_E->GetXaxis()->SetLabelSize(0.06);
   xframe_m_50_80_E->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_50_80_E->GetYaxis()->CenterTitle(true);
   xframe_m_50_80_E->GetYaxis()->SetTitleSize(0.08);
   xframe_m_50_80_E->GetYaxis()->SetLabelSize(0.06);
   xframe_m_50_80_E->SetTitleOffset(0.7, "Y");
   xframe_m_50_80_E->SetTitleOffset(0.8, "X");
   xframe_m_50_80_E->Draw();

    RooFitResult* m_80_120_E = model1.fitTo(H_moreone_80_120_E,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_80_120_E = x.frame();
    H_moreone_80_120_E.plotOn(xframe_m_80_120_E,DataError(RooAbsData::SumW2)) ;
    model1.plotOn(xframe_m_80_120_E) ;
    model1.plotOn(xframe_m_80_120_E,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model1.paramOn(xframe_m_80_120_E,Layout(0.55)) ;

   TCanvas* ca8 = new TCanvas("ca8","ca8",600,400) ;
   TLatex* N_8 = new TLatex(0.15 , 0.8 , "Earrel, E_{T}^{probe}(80-120GeV)");
   N_8->SetTextSize(0.06);
   N_8->SetNDC();
   xframe_m_80_120_E->addObject(N_8);
   xframe_m_80_120_E->addObject(Label);
   xframe_m_80_120_E->addObject(Label_BWXCB);
   xframe_m_80_120_E->SetTitle("");
   xframe_m_80_120_E->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_80_120_E->GetXaxis()->CenterTitle(true);
   xframe_m_80_120_E->GetXaxis()->SetTitleSize(0.08);
   xframe_m_80_120_E->GetXaxis()->SetLabelSize(0.06);
   xframe_m_80_120_E->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_80_120_E->GetYaxis()->CenterTitle(true);
   xframe_m_80_120_E->GetYaxis()->SetTitleSize(0.08);
   xframe_m_80_120_E->GetYaxis()->SetLabelSize(0.06);
   xframe_m_80_120_E->SetTitleOffset(0.7, "Y");
   xframe_m_80_120_E->SetTitleOffset(0.8, "X");
   xframe_m_80_120_E->Draw();

    RooFitResult* m_120_200_E = model4.fitTo(H_moreone_120_200_E,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_120_200_E = x.frame();
    H_moreone_120_200_E.plotOn(xframe_m_120_200_E,DataError(RooAbsData::SumW2)) ;
    model4.plotOn(xframe_m_120_200_E) ;
    model4.plotOn(xframe_m_120_200_E,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model4.paramOn(xframe_m_120_200_E,Layout(0.55)) ;
    
   TCanvas* ca9 = new TCanvas("ca9","ca9",600,400) ;
   TLatex* N_9 = new TLatex(0.15 , 0.8 , "Earrel, E_{T}^{probe}(120-200GeV)");
   N_9->SetTextSize(0.06);
   N_9->SetNDC();
   xframe_m_120_200_E->addObject(N_9);
   xframe_m_120_200_E->addObject(Label);
   xframe_m_120_200_E->addObject(Label_BWXGauss);
   xframe_m_120_200_E->SetTitle("");
   xframe_m_120_200_E->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_120_200_E->GetXaxis()->CenterTitle(true);
   xframe_m_120_200_E->GetXaxis()->SetTitleSize(0.08);
   xframe_m_120_200_E->GetXaxis()->SetLabelSize(0.06);
   xframe_m_120_200_E->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_120_200_E->GetYaxis()->CenterTitle(true);
   xframe_m_120_200_E->GetYaxis()->SetTitleSize(0.08);
   xframe_m_120_200_E->GetYaxis()->SetLabelSize(0.06);
   xframe_m_120_200_E->SetTitleOffset(0.7, "Y");
   xframe_m_120_200_E->SetTitleOffset(0.8, "X");
   xframe_m_120_200_E->Draw();
*/
    RooFitResult* m_200_plus_E = model5.fitTo(H_moreone_200_plus_E,Minos(kFALSE),Range(60,120),Save(kTRUE)) ;
    RooPlot* xframe_m_200_plus_E = x.frame();
    H_moreone_200_plus_E.plotOn(xframe_m_200_plus_E,DataError(RooAbsData::SumW2)) ;
    model5.plotOn(xframe_m_200_plus_E) ;
    model5.plotOn(xframe_m_200_plus_E,Components(Cheb),LineStyle(kDashed),LineColor(kRed)) ;
    model5.paramOn(xframe_m_200_plus_E,Layout(0.55)) ;
   
    TCanvas* ca10 = new TCanvas("ca10","ca10",600,400) ;
   TLatex* N_10 = new TLatex(0.15 , 0.8 , "Earrel, E_{T}^{probe}(>200GeV)");
   N_10->SetTextSize(0.06);
   N_10->SetNDC();
   xframe_m_200_plus_E->addObject(N_10);
   xframe_m_200_plus_E->addObject(Label);
   xframe_m_200_plus_E->addObject(Label_BWXLa);
   xframe_m_200_plus_E->SetTitle("");
   xframe_m_200_plus_E->SetYTitle("Number / 1 GeV/c^{2}");
   xframe_m_200_plus_E->GetXaxis()->CenterTitle(true);
   xframe_m_200_plus_E->GetXaxis()->SetTitleSize(0.08);
   xframe_m_200_plus_E->GetXaxis()->SetLabelSize(0.06);
   xframe_m_200_plus_E->SetXTitle("M_{ee}(GeV/c^{2})");
   xframe_m_200_plus_E->GetYaxis()->CenterTitle(true);
   xframe_m_200_plus_E->GetYaxis()->SetTitleSize(0.08);
   xframe_m_200_plus_E->GetYaxis()->SetLabelSize(0.06);
   xframe_m_200_plus_E->SetTitleOffset(0.7, "Y");
   xframe_m_200_plus_E->SetTitleOffset(0.8, "X");
   xframe_m_200_plus_E->Draw();



  // --- Plot frame on canvas ---
/* 
    xframe->SetTitle("");
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

   

    

//   std::cout<<"N_m_35_50_B="<<nsig1.getVal()<<" Error="<<nsig1.getError()<<" bgk="<<nbgk1.getVal()<<"error="<<<<std::endl;
//   std::cout<<"N_m_50_80_B="<<nsig2.getVal()<<" Error="<<nsig2.getError()<<std::endl;
//   std::cout<<"N_m_80_120_B="<<nsig3.getVal()<<" Error="<<nsig3.getError()<<std::endl;
//   std::cout<<"N_m_120_200_B="<<nsig4.getVal()<<" Error="<<nsig4.getError()<<std::endl;
//   std::cout<<"N_m_200_plus_B="<<nsig5.getVal()<<" Error="<<nsig5.getError()<<std::endl;
//   std::cout<<"N_m_35_50_E="<<nsig6.getVal()<<" Error="<<nsig6.getError()<<std::endl;
//   std::cout<<"N_m_50_80_E="<<nsig7.getVal()<<" Error="<<nsig7.getError()<<std::endl;
//   std::cout<<"N_m_80_120_E="<<nsig8.getVal()<<" Error="<<nsig8.getError()<<std::endl;
//   std::cout<<"N_m_120_200_E="<<nsig9.getVal()<<" Error="<<nsig9.getError()<<std::endl;
//   std::cout<<"N_m_200_plus_E="<<nsig10.getVal()<<" Error="<<nsig10.getError()<<std::endl;

}


