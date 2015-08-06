import array, math

debug = True
lumi = 40.24

##########################################################################################
#                             Import ROOT and apply settings                             #
##########################################################################################
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
#ROOT.gStyle.SetFillStyle(ROOT.kWhite)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasBorderMode(ROOT.kWhite)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadBorderMode(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
#ROOT.gStyle.SetErrorX(0)

##########################################################################################
#                                     Create canvas                                      #
##########################################################################################
canvas = ROOT.TCanvas('canvas', '', 100, 100, 1200, 800)
canvas.SetGridx()
canvas.SetGridy()

##########################################################################################
#                                   General settings                                     #
##########################################################################################
nBins_mee = 100
lower_mee =  50
upper_mee = 150

hBase_mee = ROOT.TH1F('hBase_mee', '', nBins_mee, lower_mee, upper_mee)
hBase_mee.Sumw2()
hBase_mee.GetXaxis().SetTitle('m(ee) [GeV]')
hBase_mee.GetYaxis().SetTitle('#varepsilon(HEEP)')
hBase_mee.SetMarkerStyle(20)
hBase_mee.SetMarkerColor(ROOT.kBlack)
hBase_mee.SetLineColor(ROOT.kBlack)

e1 = 'Zee_i1[Zee_highestMass]'
e2 = 'Zee_i2[Zee_highestMass]'

chargeNames   = ['ep','em','ea']
strategyNames = ['data_inc','data_exc','MC_inc','MC_exc']
variableNames = ['Et','eta','phi']
regionNames   = ['Barrel','Endcap']

#chargeNames   = ['ea']
#variableNames = ['Et']
#regionNames   = ['Barrel']

##########################################################################################
#                                        Regions                                         #
##########################################################################################
class region_object:
    def __init__(self, name, eta_lower, eta_upper):
        self.name = name
        self.eta_lower = eta_lower
        self.eta_upper = eta_upper
    def cut(self, ea):
        return 'abs(HEEP_cutflow60_eta_value[%s])>=%.6f && abs(HEEP_cutflow60_eta_value[%s])<%.6f'%(ea, self.eta_lower, ea, self.eta_upper)

regions = {}
regions['Barrel'] = region_object('Barrel', 0,  1.4442)
regions['Endcap'] = region_object('Endcap', 1.556, 2.5)

##########################################################################################
#                                      Create fit                                        #
##########################################################################################
class fitObject:
    def __init__(self):
        self.rrv_mee = ROOT.RooRealVar('m','m(ee) [GeV]', lower_mee, upper_mee)
        
        self.bw_mean  = ROOT.RooRealVar('bw_mean' , '', 91.2,  80, 100)
        self.bw_width = ROOT.RooRealVar('bw_width', '',  7.5, 0.3,  20)
        self.bw = ROOT.RooBreitWigner('bw', 'Breit Wigner', self.rrv_mee, self.bw_mean, self.bw_width)
        
        self.Landau_mean  = ROOT.RooRealVar('Landau_mean' , 'Landau mean' , 0, 0, 100)
        self.Landau_sigma = ROOT.RooRealVar('Landau_sigma', 'Landau sigma', 0.3, 0.01, 10)
        self.Landau       = ROOT.RooLandau ('Landau'      , 'Landau', self.rrv_mee, self.Landau_mean, self.Landau_sigma)
        
        self.Cheb1 = ROOT.RooRealVar('c1', 'c1',  0.71, -1.0, 1.0)
        self.Cheb2 = ROOT.RooRealVar('c2', 'c2', -0.16, -1.0, 1.0)
        self.Cheb3 = ROOT.RooRealVar('c3', 'c3', -0.10, -1.0, 1.0)
        #self.Chebychev = ROOT.RooChebychev('Cheb', 'Cheb', self.rrv_mee, ROOT.RooArgList(self.Cheb1,self.Cheb2,self.Cheb3))
        self.Chebychev = ROOT.RooChebychev('Cheb', 'Cheb', self.rrv_mee, ROOT.RooArgList(self.Cheb1,self.Cheb2))
      
        self.nReal = ROOT.RooRealVar('nReal', 'Number of real electrons', 10, 0, 100000)
        self.nFake = ROOT.RooRealVar('nFake', 'Number of fake electrons', 10, 0, 100000)
        self.bwLandau = ROOT.RooFFTConvPdf('bwxLa', '', self.rrv_mee, self.bw, self.Landau)
      
        self.model = ROOT.RooAddPdf('model', 'bkg+bw', ROOT.RooArgList(self.bwLandau, self.Chebychev), ROOT.RooArgList(self.nReal,self.nFake))
        
    def plot(self, RDH, canvas, printName, nHEEP, vartext, strname, cname, rname):
        mee_frame = self.rrv_mee.frame()
        RDH.plotOn(mee_frame, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        self.model.plotOn(mee_frame)
        self.model.plotOn(mee_frame, ROOT.RooFit.Components('Cheb'), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
        self.model.paramOn(mee_frame, ROOT.RooFit.Layout(0.55))
      
        chi2OverNDof = mee_frame.chiSquare()
        
        HEEP_text = 'HEEP+gsf' if nHEEP==1 else 'HEEP+HEEP'
        HEEP_label = ROOT.TLatex(0.15 , 0.75 , HEEP_text)
        chi_label = ROOT.TLatex(0.15 , 0.7 , '#chi^{2}/ndof = %.2f'%chi2OverNDof)
        var_label  = ROOT.TLatex(0.15 , 0.85 , vartext)
        HEEP_label.SetTextSize(0.04)
        chi_label.SetTextSize(0.04)
        var_label.SetTextSize(0.04)
        HEEP_label.SetNDC()
        chi_label.SetNDC()
        var_label.SetNDC()
        
        CMS_label = ROOT.TLatex(0.7, 0.97, 'CMS internal')
        CMS_label.SetNDC()
        
        x = 0.15
        y = 0.55
        lumi_label = ROOT.TLatex(x, y, '#int L dt = %.1f pb^{-1}'%lumi)
        lumi_label.SetNDC()
        
        beam_label = ROOT.TLatex(x, y-0.1, '#sqrt{s}=13 TeV 50 ns')
        beam_label.SetNDC()
        
        run_label = 0
        if 'data' in strname:
            run_label = ROOT.TLatex(0.05, 0.93, 'DoubleEG PromptReco Run 2015B') ;
        else:
            run_label = ROOT.TLatex(0.05, 0.93, 'RunIISpringDR74 MC') ;
        run_label.SetNDC()
        
        region_label = 0
        region_label = ROOT.TLatex(0.15, y-0.25, rname) ;
        region_label.SetNDC()
        
        e_label = 0
        if cname=='ep':
            e_label = ROOT.TLatex(0.15, y-0.35, 'e^{+}')
        elif cname=='em':
            e_label = ROOT.TLatex(0.15, y-0.35, 'e^{-}')
        else:
            e_label = ROOT.TLatex(0.15, y-0.35, 'e^{+},e^{-}')
        e_label.SetNDC()
        
        mee_frame.addObject(HEEP_label )
        #mee_frame.addObject(chi_label )
        #mee_frame.addObject(Et_label  )
        #mee_frame.addObject(eta_label )
        #mee_frame.addObject(phi_label )
        mee_frame.addObject(var_label )
        mee_frame.addObject(CMS_label )
        mee_frame.addObject(beam_label)
        mee_frame.addObject(lumi_label)
        mee_frame.addObject(run_label )
        mee_frame.addObject(e_label )
        mee_frame.addObject(region_label )
        mee_frame.SetTitle('')
        
        mee_frame.SetXTitle("m(ee) [GeV]")
        mee_frame.SetYTitle("entries per 1 GeV")
        
        mee_frame.GetXaxis().SetTitleSize(0.06)
        mee_frame.GetXaxis().SetLabelSize(0.04)
        
        mee_frame.GetYaxis().SetTitleSize(0.06)
        mee_frame.GetYaxis().SetLabelSize(0.04)
        
        mee_frame.SetTitleOffset(0.7, "Y")
        mee_frame.SetTitleOffset(0.8, "X")
        mee_frame.Draw()
        
        canvas.Print('plots/TagAndProbe_py/%s.png'%printName)

##########################################################################################
#                         Organise slices to optimise CPU time                           #
##########################################################################################
class slice_collection:
    def __init__(self, var, cut, nHEEP):
        NRegions = 2
        NChargeRegionBins = 2*NRegions+1
        
        dmee = upper_mee-lower_mee
        mee_list = []
        for i in range(0,nBins_mee+1):
            mee = lower_mee + dmee*i/nBins_mee
            mee_list.append(mee)
        meeBins = array.array('d', mee_list)
        
        rc_range = range(-NRegions,NRegions+2)
        rc_list = []
        for rc in rc_range:
            rc_list.append(rc-0.5)
        rcBins = array.array('d', rc_list)
        
        self.hBase_3D = ROOT.TH3F('hBase_3D', '', var.nBins, var.binsArray, nBins_mee, meeBins, NChargeRegionBins, rcBins)
        self.hBase_2D = ROOT.TH2F('hBase_2D', '', var.nBins, var.binsArray, nBins_mee, meeBins)
        
        self.hBase_3D.Sumw2()
        self.hBase_2D.Sumw2()
        
        h_3D_raw = {}
        h_3D_tmp = {}
        
        distribution_names = ['data','MCDY','MCOtherA','MCOtherB']
        for dname in distribution_names:
            h_3D_raw[dname] = self.hBase_3D.Clone('h_3D_%s_raw'%dname)
            h_3D_tmp[dname] = self.hBase_3D.Clone('h_3D_%s_tmp'%dname)
        
        cut_OS = 'gsf_charge[%s]*gsf_charge[%s]==-1'%(e1,e2)
        chargeRegionString = 'gsf_charge[%%s]*( 1*(%s) + 2*(%s) )'%(regions['Barrel'].cut('%s'), regions['Endcap'].cut('%s'))
        
        for sname in samples:
            s = samples[sname]
            dnames = []
            for d in distribution_names:
                if s.type=='data':
                    dnames = ['data']
                elif s.type=='MCDY':
                    dnames = ['MCDY']
                elif s.type=='MCOther':
                    dnames = ['MCOtherA','MCOtherB']
            
            for dname in dnames:
                for es in [[e1,e2],[e2,e1]]:
                    # e1 is the tag, e2 is the probe
                    drawX = '%s[%s]'%(var.branch, e2)
                    drawY = 'Zee_mass_HEEP[Zee_highestMass]'
                    drawZ = chargeRegionString%(e2,e2,e2,e2,e2)
                    drawString = '%s:%s:%s>>h_3D_%s_tmp'%(drawZ,drawY,drawX,dname)
                    total_cut = '(%s) && (%s) && (HEEP_cutflow60_Et[%s]) && (HEEP_cutflow60_total[%s]) && (%s)'%(cut, regions['Barrel'].cut(e1), e2, e1, cut_OS)
                    if nHEEP==2:
                        total_cut = '%s && (HEEP_cutflow60_total[%s])'%(total_cut, e2)
                    if dname=='MCOtherA':
                        total_cut = '(%s) && (ev_event%%2==0)'%total_cut
                    elif dname=='MCOtherB':
                        total_cut = '(%s) && (ev_event%%2==1)'%total_cut
                
                    # Fill the histogram and add it.
                    print dname , sname , s.chain.Draw(drawString, total_cut)
                    print h_3D_tmp[dname].GetSumOfWeights()
                    h_3D_raw[dname].Add(h_3D_tmp[dname])
                    
                if s.type!='data':
                    h_3D_raw[dname].Scale(lumi/s.effectiveLumi)
        
        # Now assemble the four different scenarios.
        self.h_3D_strategies = {}
        self.h_3D_strategies['data_inc'] = h_3D_raw['data'].Clone('h_3D_data_inc')
        
        self.h_3D_strategies['data_exc'] = h_3D_raw['data'].Clone('h_3D_data_exc')
        self.h_3D_strategies['data_exc'].Add(h_3D_raw['MCOtherA'], -1)
        self.h_3D_strategies['data_exc'].Add(h_3D_raw['MCOtherB'], -1)
        
        self.h_3D_strategies['MC_inc'  ] = h_3D_raw['MCDY'].Clone('h_3D_MC_inc')
        self.h_3D_strategies['MC_inc'  ].Add(h_3D_raw['MCOtherA'])
        self.h_3D_strategies['MC_inc'  ].Add(h_3D_raw['MCOtherB'])
        
        self.h_3D_strategies['MC_exc'  ] = h_3D_raw['MCDY'].Clone('h_3D_MC_exc')
        self.h_3D_strategies['MC_exc'  ].Add(h_3D_raw['MCOtherA'],  2)
        self.h_3D_strategies['MC_exc'  ].Add(h_3D_raw['MCOtherB'], -2)
        
        centreBin = (NChargeRegionBins+1)/2
        self.slice_lists = {}
        for strname in strategyNames:
            self.slice_lists[strname] = {}
            for rname in regionNames:
                self.slice_lists[strname][rname] = {}
                for cname in ['ep','em']:
                    h = self.hBase_2D.Clone('h_2D_mee_%s_%s_%s_%d'%(var.name, rname, cname, nHEEP))
                    charge = 0
                    if cname=='ep':
                        charge =  1
                    elif cname=='em':
                        charge = -1
                    rBin = 1 if rname=='Barrel' else 2
                    rcBin = centreBin + charge*rBin
                
                    for mBin in range(1, h.GetNbinsY()+1):
                        for vBin in range(1, h.GetNbinsX()+1):
                            h.SetBinContent(vBin, mBin, self.h_3D_strategies[strname].GetBinContent(vBin, mBin, rcBin))
                            h.SetBinError  (vBin, mBin, self.h_3D_strategies[strname].GetBinError  (vBin, mBin, rcBin))
                
                    self.slice_lists[strname][rname][cname] = h
            
                # Combine to get both charges
                h = self.slice_lists[strname][rname]['ep'].Clone('h_2D_mee_%s_%s_ea_%d'%(var.name, rname, nHEEP))
                h.Add(self.slice_lists[strname][rname]['em'])
                self.slice_lists[strname][rname]['ea'] = h
    
    def get_histogram(self, strname, rname, cname):
        return self.slice_lists[strname][rname][cname]
        

##########################################################################################
#                                   Create variables                                     #
##########################################################################################
class variableObject:
    def __init__(self, name, latex, unit, branch, binBoundaries):
        self.name = name
        self.latex = latex
        self.unit = unit
        self.branch = branch
        self.binBoundaries = binBoundaries
    def make_hBases(self):
        self.binBoundaries = sorted(self.binBoundaries)
        self.nBins = len(self.binBoundaries)-1
        self.lower = self.binBoundaries[ 0]
        self.upper = self.binBoundaries[-1]
        self.binsArray = array.array('d', self.binBoundaries)
        self.hBase_1D = ROOT.TH1F('hBase_%s'    %self.name, '', self.nBins, self.binsArray)
        self.hBase_1D.GetXaxis().SetTitle('%s'%self.name)
        self.hBase_1D.GetYaxis().SetTitle('#varepsilon(HEEP)')
        self.hBase_1D.SetLineColor(ROOT.kBlack)
        self.hBase_1D.SetMarkerColor(ROOT.kBlack)
        self.hBase_1D.SetMarkerStyle(20)
        
        self.hBase_2D = ROOT.TH2F('hBase_%s_mee'%self.name, '', self.nBins, self.binsArray, nBins_mee, lower_mee, upper_mee)
        self.hBase_2D.GetXaxis().SetTitle('%s'%self.name)
        self.hBase_2D.GetYaxis().SetTitle('m(ee) [GeV]')

pi = 3.14156

phi_boundaries = []
nPhi = 10
for i in range(-nPhi,nPhi+1):
    phi_boundaries.append(i*pi/nPhi)

Et_boundaries  = [35,40,45,50,60,70,80,90,100,125,150,250]
eta_boundaries = [-2.5,-2.25,-2.0,-1.75,-1.566,-1.4442,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.4442,1.566,1.75,2.0,2.25,2.5]

variables = {}
variables['Et' ] = variableObject('Et' , 'E_{T}(probe)', '[GeV]', 'HEEP_cutflow60_Et_value' , Et_boundaries)
variables['eta'] = variableObject('eta', '#eta(probe)' , ''     , 'HEEP_cutflow60_eta_value', eta_boundaries)
variables['phi'] = variableObject('phi', '#phi(probe)' , '', 'gsf_phi', phi_boundaries)

class sliceObject():
    def __init__(self, name, lower, upper, histogram):
        self.name = name
        self.lower = lower
        self.upper = upper
        self.histogram = histogram
    def combine(self, other):
        self.histogram.Add(other.histogram)

##########################################################################################
#                                     Create samples                                     #
##########################################################################################
class sampleObject:
    def __init__(self, name, filename, crossSection, nEvents, type, color):
        self.name = name
        self.filename = filename
        self.crossSection = crossSection
        self.nEvents = nEvents
        self.effectiveLumi = self.nEvents/self.crossSection
        self.type = type # 'Data', 'MCDY', 'MCOther'
        self.color = color
        self.chain = ROOT.TChain('IIHEAnalysis')
        self.chain.Add('%s*'%self.filename)
    def add_bin(self, bin):
        self.binBoundaries.append(bin)
    

samples = {}
samples['DoubleEG'] = sampleObject('DoubleEG', '../data2015/DoubleEG_Run2015B_GoldenLumimask.root', -1,       -1, 'data'   , ROOT.kBlack )
samples['DY'      ] = sampleObject('DY_M50'  , 'TagAndProbe_DY_M50'  ,  3*2008.4, 14439958, 'MCDY'   , ROOT.kYellow)
samples['ttbar'   ] = sampleObject('ttbar'   , 'TagAndProbe_ttbar'   ,    831.76, 19396788, 'MCOther', ROOT.kBlue  )
samples['WJets'   ] = sampleObject('WJets'   , 'TagAndProbe_WJets'   , 2*20508.9, 24132876, 'MCOther', ROOT.kGreen )

for vname in variables:
    v = variables[vname]
    v.make_hBases()

def cut_probe_HEEP(el):
    return 'HEEP_cutflow60_ID[%s] && HEEP_cutflow60_isolation[%s]'%(el,el)

chargeCuts = {'ep':'gsf_chrage[%s]==1', 'em':'gsf_charge[%s]==-1', 'ea':'gsf_charge[%s]!=0'}
cut_OS = 'gsf_charge[%s]*gsf_charge[%s]<0'%(e1,e2)

eff_histos = {}
#for vname in variables:
for vname in variableNames:
    eff_histos[vname] = {}
    v = variables[vname]
    
    slices_HEEP = {}
    slices_HEEP[1] = slice_collection(v, '1==1', 1)
    slices_HEEP[2] = slice_collection(v, '1==1', 2)
    
    for rname in regionNames:
        eff_histos[vname][rname] = {}
        for cname in chargeNames:
            eff_histos[vname][rname][cname] = {}
            for strname in strategyNames:
                nReal_values = {}
                nReal_errors = {}
                for iHEEP in range(1,3):
                    h_2D = slices_HEEP[iHEEP].get_histogram(strname, rname, cname)
                    nReal_values [iHEEP] = {}
                    nReal_errors [iHEEP] = {}
                    hHEEP = hBase_mee.Clone('h_mee_%s_%s_%d'%(strname, rname, iHEEP))
                    for vBin in range(1, h_2D.GetNbinsX()+1):
                        h_1D = hBase_mee.Clone('h_mee')
                        for mBin in range(1, h_2D.GetNbinsY()+1):
                            h_1D.SetBinContent(mBin, h_2D.GetBinContent(vBin, mBin))
                            h_1D.SetBinError  (mBin, h_2D.GetBinError  (vBin, mBin))
                        
                        # Make the fit object here so that all slices start with exactly the same
                        # starting parameters.
                        f = fitObject()
                        RDH = ROOT.RooDataHist('RDH', '', ROOT.RooArgList(f.rrv_mee), h_1D)
                        result = f.model.fitTo(RDH, ROOT.RooFit.Minos(False), ROOT.RooFit.Range(lower_mee, upper_mee), ROOT.RooFit.Save(True))
                        
                        lower =         h_2D.GetXaxis().GetBinLowEdge(vBin)
                        upper = lower + h_2D.GetXaxis().GetBinWidth  (vBin)
                        varText = '%.2f < %s < %.2f %s'%(lower, v.latex, upper, v.unit)
                        pname = 'fit_%s_%s_%s_%s_%.0f_%.0f__%dHEEP'%(rname, strname, cname, v.name, lower*1000, upper*1000, iHEEP)
                        f.plot(RDH, canvas, pname, iHEEP, varText, strname, cname, rname)
                        
                        nReal_values[iHEEP][vBin] = f.nReal.getVal()
                        nReal_errors[iHEEP][vBin] = f.nReal.getError()
                
                h = v.hBase_1D.Clone('h_tmp_eff_%s_%s_%s_%s_%d'%(vname, strname, rname, cname, iHEEP))
                print '%s  %s  %s  %s'%(rname,strname,cname,vname)
                for vBin in range(1, h.GetNbinsX()+1):
                    value_1HEEP = nReal_values[1][vBin]
                    error_1HEEP = nReal_errors[1][vBin]
                    value_2HEEP = nReal_values[2][vBin]
                    error_2HEEP = nReal_errors[2][vBin]
                
                    r = 0.0 if value_1HEEP<1e-6 else value_2HEEP/value_1HEEP
                    eff = r/(2-r)
                    if eff>1 and value_2HEEP-value_1HEEP < 2*(error_1HEEP+error_2HEEP):
                        eff = 1.0
                    err = 0.0 if value_1HEEP<1e-6 or eff<0 or eff>1 else math.sqrt(eff*(1-eff)/value_1HEEP)
                
                    print '%8.2f +- %8.2f  ,  %8.2f +- %8.2f  , %8.4f +- %8.4f'%(value_1HEEP, error_1HEEP, value_2HEEP, error_2HEEP, eff, err)
        
                    h.SetBinContent(vBin, eff)
                    h.SetBinError  (vBin, err)
            
                h.SetMinimum(0.0)
                h.SetMaximum(1.0)
                h.Draw('pe')
                canvas.Print('plots/TagAndProbe_py/h_eff_%s_%s_%s_%s.png'%(rname,strname,cname,v.name))
            
                eff_histos[vname][rname][cname][strname] = h.Clone('h_eff_%s_%s_%s_%s_%d'%(vname, strname, rname, cname, iHEEP))
                print '$$ ' , strname , eff_histos[vname][rname][cname][strname]

for vname in eff_histos:
    for rname in eff_histos[vname]:
        for cname in eff_histos[vname][rname]:
            print vname , rname , cname
            for strname in eff_histos[vname][rname][cname]:
                print '%% ' , strname , eff_histos[vname][rname][cname][strname]
            h_data_inc = eff_histos[vname][rname][cname]['data_inc']
            h_data_exc = eff_histos[vname][rname][cname]['data_exc']
            h_MC_inc   = eff_histos[vname][rname][cname]['MC_inc'  ]
            h_MC_exc   = eff_histos[vname][rname][cname]['MC_exc'  ]
            
            h_data_inc.SetLineColor(ROOT.kBlack)
            h_data_inc.SetMarkerColor(ROOT.kBlack)
            h_data_inc.SetMarkerStyle(20)
            
            h_data_exc.SetLineColor(ROOT.kRed)
            h_data_exc.SetMarkerColor(ROOT.kRed)
            h_data_exc.SetMarkerStyle(24)
            
            h_MC_inc.SetLineColor(ROOT.kBlue)
            h_MC_inc.SetMarkerColor(ROOT.kBlue)
            h_MC_inc.SetMarkerStyle(21)
            
            h_MC_exc.SetLineColor(ROOT.kMagenta)
            h_MC_exc.SetMarkerColor(ROOT.kMagenta)
            h_MC_exc.SetMarkerStyle(25)
            
            h_data_inc.Draw('pe')
            h_data_exc.Draw('pe:sames')
            h_MC_inc  .Draw('pe:sames')
            h_MC_exc  .Draw('pe:sames')
            h_data_inc.Draw('axis:sames')
            
            legend = ROOT.TLegend(0.5, 0.5, 0.9, 0.1)
            legend.SetFillStyle(0)
            legend.SetShadowColor(0)
            legend.SetBorderSize(0)
            legend.AddEntry(h_data_inc, 'Data (all)'              , 'pl')
            legend.AddEntry(h_data_exc, 'Data (non-DY subtracted)', 'pl')
            legend.AddEntry(h_MC_inc  , 'MC (all)'                , 'pl')
            legend.AddEntry(h_MC_exc  , 'MC (non-DY subtracted)'  , 'pl')
            
            CMS_label = ROOT.TLatex(0.7, 0.93, 'CMS internal')
            CMS_label.SetNDC()
            
            x = 0.60
            y = 0.38
            lumi_label = ROOT.TLatex(x, y, '#int L dt = %.1f pb^{-1}'%lumi)
            lumi_label.SetNDC()
            
            beam_label = ROOT.TLatex(x, y-0.1, '#sqrt{s}=13 TeV 50 ns')
            beam_label.SetNDC()
            
            run_label = ROOT.TLatex(0.05, 0.93, 'DoubleEG PromptReco Run 2015B') ;
            run_label.SetNDC()
            
            CMS_label.Draw()
            lumi_label.Draw()
            beam_label.Draw()
            run_label.Draw()
            
            legend.Draw()
            
            canvas.Print('plots/TagAndProbe_py/h_compare_%s_%s_%s.png'%(vname,rname,cname))

fOut = ROOT.TFile('TagAndProbe_histos.root','RECREATE')
for vname in variableNames:
    for rname in regionNames:
        for cname in chargeNames:
            h_data_inc = eff_histos[vname][rname][cname]['data_inc']
            h_data_exc = eff_histos[vname][rname][cname]['data_exc']
            h_MC_inc   = eff_histos[vname][rname][cname]['MC_inc'  ]
            h_MC_exc   = eff_histos[vname][rname][cname]['MC_exc'  ]
        
            h_scale_inc = h_data_inc.Clone('h_SF_%s_%s_%s_inc'%(vname,rname,cname))
            h_scale_inc.Divide(h_MC_inc)
            h_scale_inc.SetMaximum(2.0)
            
            h_scale_inc.Draw('pe')
            canvas.Print('plots/TagAndProbe_py/h_SF_%s_%s_%s_inc.png'%(vname,rname,cname))
        
            h_scale_exc = h_data_inc.Clone('h_SF_%s_%s_%s_exc'%(vname,rname,cname))
            h_scale_exc.Divide(h_MC_exc)
            h_scale_exc.SetMaximum(2.0)
            
            h_scale_exc.Draw('pe')
            canvas.Print('plots/TagAndProbe_py/h_SF_%s_%s_%s_exc.png'%(vname,rname,cname))
            
            for bin in range(1, h_scale_exc.GetNbinsX()+1):
                print bin , h_scale_exc.GetBinContent(bin) , h_scale_exc.GetBinError(bin) , h_scale_inc.GetBinContent(bin) , h_scale_inc.GetBinError(bin)
            
            h_data_inc.Write()
            h_data_exc.Write()
            h_MC_inc  .Write()
            h_MC_exc  .Write()
            h_scale_inc.Write()
            h_scale_exc.Write()
                

fOut.Close()
            