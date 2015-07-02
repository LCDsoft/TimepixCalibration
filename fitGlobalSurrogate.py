# will read the points and uncertainties from text files,
# fit the surrogate function and save the result
# for options run python surrogateFitterGlobal.py -h

import ROOT as R
from array import array
import ast
from optparse import OptionParser
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
R.gStyle.SetOptFit(1)
R.gROOT.SetBatch(True)


def fitGlobalSurrogate(assembly):
    lines = []
    fixed_vals = []
    Thline = 0

    # thresholds
    if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
        fTh = open('results/threshold/%s_GlobalThreshold.txt' %assembly,'r')
        Thline = fTh.readline().split()
        lines.append(Thline)
        fixed_vals.append(0.)

    # CERN data
    fFe = open('results/kde/%s_Fe_GlobalResults.txt' %assembly,'r')
    fAm = open('results/kde/%s_Am_GlobalResults.txt' %assembly,'r')
    fCuInXRF = open('results/kde/%s_CuInXRF_GlobalResults.txt' %assembly,'r')
    fCo = open('results/kde/%s_Co_GlobalResults.txt' %assembly,'r')
    fCd = open('results/kde/%s_Cd_GlobalResults.txt' %assembly,'r')
    Feline = fFe.readline().split()
    Amline = fAm.readline().split()
    CuInXRFline = fCuInXRF.readline().split()
    Coline = fCo.readline().split()
    Cdline = fCd.readline().split()
    lines.append(Feline)
    fixed_vals.append(C.FePeakE)
    lines.append([Amline[3],Amline[4],Amline[5],Amline[9]])
    fixed_vals.append(C.Am2PeakE)
    lines.append([Amline[6],Amline[7],Amline[8],Amline[9]])
    fixed_vals.append(C.Am3PeakE)
    lines.append([CuInXRFline[0],CuInXRFline[1],CuInXRFline[2],CuInXRFline[6]])
    fixed_vals.append(C.CuXRFPeakE)
    lines.append([CuInXRFline[3],CuInXRFline[4],CuInXRFline[5],CuInXRFline[6]])
    fixed_vals.append(C.InXRFPeakE)
    lines.append([Coline[0],Coline[1],Coline[2],Coline[6]])
    fixed_vals.append(C.Co1PeakE)
    lines.append([Coline[3],Coline[4],Coline[5],Coline[6]])
    fixed_vals.append(C.Co2PeakE)
    lines.append(Cdline)
    fixed_vals.append(C.CdPeakE)

    # LNLS data
    if assembly == "A06-W0110":
        fCoXRF = open('results/kde/%s_CoXRF_GlobalResults.txt' %assembly,'r')
        fCrXRF = open('results/kde/%s_CrXRF_GlobalResults.txt' %assembly,'r')
        fCuXRF = open('results/kde/%s_CuXRF_GlobalResults.txt' %assembly,'r')
        fFeXRF = open('results/kde/%s_FeXRF_GlobalResults.txt' %assembly,'r')
        fMnXRF = open('results/kde/%s_MnXRF_GlobalResults.txt' %assembly,'r')
        fNiXRF = open('results/kde/%s_NiXRF_GlobalResults.txt' %assembly,'r')
        fTiXRF = open('results/kde/%s_TiXRF_GlobalResults.txt' %assembly,'r')
        fVXRF = open('results/kde/%s_VXRF_GlobalResults.txt' %assembly,'r')
        CoXRFline = fCoXRF.readline().split()
        CrXRFline = fCrXRF.readline().split()
        CuXRFline = fCuXRF.readline().split()
        FeXRFline = fFeXRF.readline().split()
        MnXRFline = fMnXRF.readline().split()
        NiXRFline = fNiXRF.readline().split()
        TiXRFline = fTiXRF.readline().split()
        VXRFline = fVXRF.readline().split()
        lines.append(CoXRFline)
        fixed_vals.append(C.CoXRFPeakE)
        lines.append(CrXRFline)
        fixed_vals.append(C.CrXRFPeakE)
        lines.append(CuXRFline)
        fixed_vals.append(C.CuXRFPeakE)
        lines.append(FeXRFline)
        fixed_vals.append(C.FeXRFPeakE)
        lines.append(MnXRFline)
        fixed_vals.append(C.MnXRFPeakE)
        lines.append(NiXRFline)
        fixed_vals.append(C.NiXRFPeakE)
        lines.append(TiXRFline)
        fixed_vals.append(C.TiXRFPeakE)
        lines.append(VXRFline)
        fixed_vals.append(C.VXRFPeakE)

    for line in lines:
        for i in xrange(len(line)):
            line[i] = ast.literal_eval(line[i])

    tots = array('f',[])
    tot_lerrs = array('f',[])
    tot_uerrs = array('f',[])

    energies = array('f',[])
    energy_lerrs = array('f',[])
    energy_uerrs = array('f',[])

    for line, fixed_val in zip(lines,fixed_vals):
        if line == Thline:
            tots.append(fixed_val)
            tot_lerrs.append(R.sqrt(10.0**2 + 5.0**2))
            tot_uerrs.append(R.sqrt(10.0**2 + 5.0**2))
            
            energies.append(line[0])
            energy_lerrs.append(line[1])
            energy_uerrs.append(line[2])

        else:
            tots.append(line[0])
            tot_lerrs.append(R.sqrt(line[1]**2 + line[3]**2) + R.sqrt(10.0**2 + 5.0**2))
            tot_uerrs.append(R.sqrt(line[2]**2 + line[3]**2) + R.sqrt(10.0**2 + 5.0**2))

            energies.append(fixed_val)
            energy_lerrs.append(0.)
            energy_uerrs.append(0.)


    canv = R.TCanvas()
    canv.SetRightMargin(0.01)
    canv.SetLeftMargin(0.19)
    gr = R.TGraphAsymmErrors(len(energies),energies,tots,energy_lerrs,energy_uerrs,tot_lerrs,tot_uerrs)
    gr.SetMarkerStyle(20)
    gr.GetXaxis().SetTitle("Energy [keV]")
    gr.GetYaxis().SetTitle("TOT [1/96 MHz]")
    gr.GetYaxis().SetTitleOffset(1.4)
    gr.Draw('AP')

    surrogate = R.TF1("surrogate","[0]*x+[1]-([2]/(x-[3]))",2,60)
    surrogate.SetParName(0,'a')
    surrogate.SetParName(1,'b')
    surrogate.SetParName(2,'c')
    surrogate.SetParName(3,'t')
    surrogate.SetParameters(13.9,333.2,1220,-0.1)
    gr.Fit("surrogate","RQ")
    canv.Update()

    stats = gr.GetListOfFunctions().FindObject("stats")
    stats.SetX1NDC(0.57)
    stats.SetX2NDC(0.96)
    stats.SetY1NDC(0.21)
    stats.SetY2NDC(0.46)
    stats.SetBorderSize(0)
    gr.SetMinimum(0.)
    gr.SetMaximum(surrogate.Eval(65.0))
    canv.Update()
    canv.SaveAs("plots/KDESurrogateFits/Global/%s_GlobalSurrogateFit.pdf" %assembly)

    if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
        fTh.close()
    fFe.close()
    fAm.close()
    fCuInXRF.close()
    fCo.close()
    fCd.close()
    if assembly == "A06-W0110":
        fCoXRF.close()
        fCrXRF.close()
        fCuXRF.close()
        fFeXRF.close()
        fMnXRF.close()
        fNiXRF.close()
        fTiXRF.close()
        fVXRF.close()


fitGlobalSurrogate(assembly)

