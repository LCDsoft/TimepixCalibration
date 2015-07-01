# file to plot superimposed spectra of the calibration data

import ROOT as R
from optparse import OptionParser
from os import environ
import Constants as C

parser = OptionParser()
parser.add_option("-b", "--assembly",
                  help="Assembly name", dest="ASSEMBLY")

(options, args) = parser.parse_args()

if(options.ASSEMBLY):
    assembly=options.ASSEMBLY
else :
    print "Please specify assembly"
    print "choose from", C.known_assemblies
    parser.print_help()
    exit()

if assembly not in C.known_assemblies:
    print "Assembly not recognised"
    print "choose from", C.known_assemblies
    exit()

# import CLICdp ROOT style
R.gROOT.ProcessLine(".L CLICdpStyle.C")
R.CLICdpStyle()

# extra style


def plotSpectra(assembly,topTOT):
    c = R.TCanvas("c","",0,0,1000,400)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.1)

    # get data
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    rootfileFe = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    rootfileAm = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    rootfileCd = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    if assembly == "B06-W0125":
        rootfileCuIn = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
    else:
        rootfileCuIn = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    if assembly == "B06-W0125":
        rootfileCo = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
    else:
        rootfileCo = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))

    # trees
    treeFe = rootfileFe.Get("pixels")
    treeAm = rootfileAm.Get("pixels")
    treeCd = rootfileCd.Get("pixels")
    treeCuIn = rootfileCuIn.Get("pixels")
    treeCo = rootfileCo.Get("pixels")

    # hists
    histFe_allp = R.TH1F("histFe_allp","",100,0,topTOT)
    histAm_allp = R.TH1F("histAm_allp","",100,0,topTOT)
    histCd_allp = R.TH1F("histCd_allp","",100,0,topTOT)
    histCuIn_allp = R.TH1F("histCuIn_allp","",100,0,topTOT)
    histCo_allp = R.TH1F("histCo_allp","",100,0,topTOT)

    histFe_onep = R.TH1F("histFe_onep","",20,0,topTOT)
    histAm_onep = R.TH1F("histAm_onep","",20,0,topTOT)
    histCd_onep = R.TH1F("histCd_onep","",20,0,topTOT)
    histCuIn_onep = R.TH1F("histCuIn_onep","",20,0,topTOT)
    histCo_onep = R.TH1F("histCo_onep","",20,0,topTOT)

    treeFe.Draw("tot>>histFe_allp","","goff") 
    treeAm.Draw("tot>>histAm_allp","","goff")
    treeCd.Draw("tot>>histCd_allp","","goff")
    treeCuIn.Draw("tot>>histCuIn_allp","","goff")
    treeCo.Draw("tot>>histCo_allp","","goff")

    treeFe.Draw("tot>>histFe_onep","col==128&&row==128","goff") 
    treeAm.Draw("tot>>histAm_onep","col==128&&row==128","goff")
    treeCd.Draw("tot>>histCd_onep","col==128&&row==128","goff")
    treeCuIn.Draw("tot>>histCuIn_onep","col==128&&row==128","goff")
    treeCo.Draw("tot>>histCo_onep","col==128&&row==128","goff")

    for hist in [histFe_allp,histAm_allp,histCd_allp,histCuIn_allp,histCo_allp,
                 histFe_onep,histAm_onep,histCd_onep,histCuIn_onep,histCo_onep]:
            hist.Scale(1./hist.GetMaximum())

    for hist in [histFe_allp,histFe_onep]:
        hist.GetXaxis().SetTitle("TOT [1/96 MHz]")
        hist.GetYaxis().SetTitle("Pixels / Peak # pixels")
        hist.GetYaxis().SetTitleOffset(0.5)
    
    histFe_allp.SetLineColor(R.kBlue)
    histAm_allp.SetLineColor(R.kGreen+2)
    histCd_allp.SetLineColor(R.kRed+2)
    histCuIn_allp.SetLineColor(R.kOrange)
    histCo_allp.SetLineColor(R.kCyan+1)

    histFe_onep.SetLineColor(R.kBlue)
    histAm_onep.SetLineColor(R.kGreen+2)
    histCd_onep.SetLineColor(R.kRed+2)
    histCuIn_onep.SetLineColor(R.kOrange)
    histCo_onep.SetLineColor(R.kCyan+1)

    leg = R.TLegend(0.6,0.45,0.79,0.89)
    leg.AddEntry(histAm_allp,"^{241}Am","l")
    leg.AddEntry(histCd_allp,"^{109}Cd","l")
    leg.AddEntry(histCo_allp,"^{57}Co","l")
    leg.AddEntry(histCuIn_allp,"CuInXRF","l")
    leg.AddEntry(histFe_allp,"^{55}Fe","l")
    
    histFe_allp.Draw()
    histAm_allp.Draw("same")
    histCd_allp.Draw("same")
    histCuIn_allp.Draw("same")
    histCo_allp.Draw("same")
    leg.Draw()
    c.Update()
    c.SaveAs("plots/Spectra/%s_Spectra_AllPixels.pdf" %assembly)

    histFe_onep.Draw()
    histAm_onep.Draw("same")
    histCd_onep.Draw("same")
    histCuIn_onep.Draw("same")
    histCo_onep.Draw("same")
    leg.Draw()
    c.Update()
    c.SaveAs("plots/Spectra/%s_Spectra_OnePixel.pdf" %assembly)

    # if A06-W0110, plot additional spectra of LNLS data
    if (assembly == "A06-W0110"):
        print "Assembly A06-W0110 has additional data from LNLS, plotting that now"

        # get data
        base = "%s/eos/clicdp/data/VertexCalibration/LNLS_Analysis/SinglePixelAnalysis/root_files" %home
        topTOT = 500

        rootfileCoXRF = R.TFile("%s/A06-W0110-25V_CoXRF_CalibTree.root" %base)
        rootfileCrXRF = R.TFile("%s/A06-W0110-25V_CrXRF_CalibTree.root" %base)
        rootfileCuXRF = R.TFile("%s/A06-W0110-25V_CuXRF_CalibTree.root" %base)
        rootfileFeXRF = R.TFile("%s/A06-W0110-25V_FeXRF_CalibTree.root" %base)
        rootfileMnXRF = R.TFile("%s/A06-W0110-25V_MnXRF_CalibTree.root" %base)
        rootfileNiXRF = R.TFile("%s/A06-W0110-25V_NiXRF_CalibTree.root" %base)
        rootfileTiXRF = R.TFile("%s/A06-W0110-25V_TiXRF_CalibTree.root" %base)
        rootfileVXRF = R.TFile("%s/A06-W0110-25V_VXRF_CalibTree.root" %base)

        # trees
        treeCoXRF = rootfileCoXRF.Get("pixels")
        treeCrXRF = rootfileCrXRF.Get("pixels")
        treeCuXRF = rootfileCuXRF.Get("pixels")
        treeFeXRF = rootfileFeXRF.Get("pixels")
        treeMnXRF = rootfileMnXRF.Get("pixels")
        treeNiXRF = rootfileNiXRF.Get("pixels")
        treeTiXRF = rootfileTiXRF.Get("pixels")
        treeVXRF = rootfileVXRF.Get("pixels")
        print "got trees"

        # hists
        histCoXRF_allp = R.TH1F("histCoXRF_allp","",100,0,topTOT)
        histCrXRF_allp = R.TH1F("histCrXRF_allp","",100,0,topTOT)
        histCuXRF_allp = R.TH1F("histCuXRF_allp","",100,0,topTOT)
        histFeXRF_allp = R.TH1F("histFeXRF_allp","",100,0,topTOT)
        histMnXRF_allp = R.TH1F("histMnXRF_allp","",100,0,topTOT)
        histNiXRF_allp = R.TH1F("histNiXRF_allp","",100,0,topTOT)
        histTiXRF_allp = R.TH1F("histTiXRF_allp","",100,0,topTOT)
        histVXRF_allp = R.TH1F("histVXRF_allp","",100,0,topTOT)
        
        histCoXRF_onep = R.TH1F("histCoXRF_onep","",20,0,topTOT)
        histCrXRF_onep = R.TH1F("histCrXRF_onep","",20,0,topTOT)
        histCuXRF_onep = R.TH1F("histCuXRF_onep","",20,0,topTOT)
        histFeXRF_onep = R.TH1F("histFeXRF_onep","",20,0,topTOT)
        histMnXRF_onep = R.TH1F("histMnXRF_onep","",20,0,topTOT)
        histNiXRF_onep = R.TH1F("histNiXRF_onep","",20,0,topTOT)
        histTiXRF_onep = R.TH1F("histTiXRF_onep","",20,0,topTOT)
        histVXRF_onep = R.TH1F("histVXRF_onep","",20,0,topTOT)

        treeCoXRF.Draw("tot>>histCoXRF_allp","","goff") 
        treeCrXRF.Draw("tot>>histCrXRF_allp","","goff") 
        treeCuXRF.Draw("tot>>histCuXRF_allp","","goff") 
        treeFeXRF.Draw("tot>>histFeXRF_allp","","goff") 
        treeMnXRF.Draw("tot>>histMnXRF_allp","","goff") 
        treeNiXRF.Draw("tot>>histNiXRF_allp","","goff") 
        treeTiXRF.Draw("tot>>histTiXRF_allp","","goff") 
        treeVXRF.Draw("tot>>histVXRF_allp","","goff") 

        treeCoXRF.Draw("tot>>histCoXRF_onep","col==128&&row==128","goff") 
        treeCrXRF.Draw("tot>>histCrXRF_onep","col==128&&row==128","goff") 
        treeCuXRF.Draw("tot>>histCuXRF_onep","col==128&&row==128","goff") 
        treeFeXRF.Draw("tot>>histFeXRF_onep","col==128&&row==128","goff") 
        treeMnXRF.Draw("tot>>histMnXRF_onep","col==128&&row==128","goff") 
        treeNiXRF.Draw("tot>>histNiXRF_onep","col==128&&row==128","goff") 
        treeTiXRF.Draw("tot>>histTiXRF_onep","col==128&&row==128","goff") 
        treeVXRF.Draw("tot>>histVXRF_onep","col==128&&row==128","goff") 

        for hist in [histCoXRF_allp,histCrXRF_allp,histCuXRF_allp,histFeXRF_allp,histMnXRF_allp,
                     histNiXRF_allp,histTiXRF_allp,histVXRF_allp,histCoXRF_onep,histCrXRF_onep,
                     histCuXRF_onep,histFeXRF_onep,histMnXRF_onep,histNiXRF_onep,histTiXRF_onep,histVXRF_onep]:
            hist.Scale(1./hist.GetMaximum())

        for hist in [histCuXRF_allp,histCuXRF_onep]:
            hist.GetXaxis().SetTitle("TOT [1/96 MHz]")
            hist.GetYaxis().SetTitle("Pixels / Peak # pixels")
            hist.GetYaxis().SetTitleOffset(0.5)

        histCoXRF_allp.SetLineColor(R.kBlue)
        histCrXRF_allp.SetLineColor(R.kGreen+2)
        histCuXRF_allp.SetLineColor(R.kRed+2)
        histFeXRF_allp.SetLineColor(R.kOrange)
        histMnXRF_allp.SetLineColor(R.kCyan+1)
        histNiXRF_allp.SetLineColor(R.kViolet)
        histTiXRF_allp.SetLineColor(R.kGray)
        histVXRF_allp.SetLineColor(R.kBlack)

        histCoXRF_onep.SetLineColor(R.kBlue)
        histCrXRF_onep.SetLineColor(R.kGreen+2)
        histCuXRF_onep.SetLineColor(R.kRed+2)
        histFeXRF_onep.SetLineColor(R.kOrange)
        histMnXRF_onep.SetLineColor(R.kCyan+1)
        histNiXRF_onep.SetLineColor(R.kViolet)
        histTiXRF_onep.SetLineColor(R.kGray)
        histVXRF_onep.SetLineColor(R.kBlack)

        leg = R.TLegend(0.8,0.2,0.94,0.89)
        leg.AddEntry(histCoXRF_allp,"CoXRF","l")
        leg.AddEntry(histCrXRF_allp,"CrXRF","l")
        leg.AddEntry(histCuXRF_allp,"CuXRF","l")
        leg.AddEntry(histFeXRF_allp,"FeXRF","l")
        leg.AddEntry(histMnXRF_allp,"MnXRF","l")
        leg.AddEntry(histNiXRF_allp,"NiXRF","l")
        leg.AddEntry(histTiXRF_allp,"TiXRF","l")
        leg.AddEntry(histVXRF_allp,"VXRF","l")

        histCuXRF_allp.Draw()
        histCrXRF_allp.Draw("same")
        histVXRF_allp.Draw("same")
        histFeXRF_allp.Draw("same")
        histMnXRF_allp.Draw("same")
        histNiXRF_allp.Draw("same")
        histTiXRF_allp.Draw("same")
        histCoXRF_allp.Draw("same")
        leg.Draw()
        c.Update()
        c.SaveAs("plots/Spectra/A06-W0110_LNLSSpectra_AllPixels.pdf")

        histCuXRF_onep.Draw()
        histCrXRF_onep.Draw("same")
        histVXRF_onep.Draw("same")
        histFeXRF_onep.Draw("same")
        histMnXRF_onep.Draw("same")
        histNiXRF_onep.Draw("same")
        histTiXRF_onep.Draw("same")
        histCoXRF_onep.Draw("same")
        leg.Draw()
        c.Update()
        c.SaveAs("plots/Spectra/A06-W0110_LNLSSpectra_OnePixel.pdf")

def getTopTOT(assembly):
    
    if assembly == "A06-W0110": topTOT = 1300
    if assembly == "B06-W0125": topTOT = 2700
    if assembly == "B07-W0125": topTOT = 1400
    if assembly == "C04-W0110": topTOT = 1300
    if assembly == "D09-W0126": topTOT = 1700
    if assembly == "L04-W0125": topTOT = 1500

    return topTOT


topTOT = getTopTOT(assembly)
plotSpectra(assembly,topTOT)
