# file to read the points and uncertainties from text files,
# fit the surrogate function
# and write the constants and uncertainties to another text file

import ROOT as R
from array import array
import ast
import numpy as np
import matplotlib.pyplot as plt
import itertools
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


def fitPixelSurrogates(assembly,parEsts,parLLims,parULims):

    # thresholds
    Thline = 0
    if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
        fTh = open('results/threshold/%s_PixelThreshold.txt' %assembly,'r')

    # CERN data
    fFe = open('results/kde/%s_Fe_PixelResults.txt' %assembly,'r')
    fAm2 = open('results/kde/%s_Am2_PixelResults.txt' %assembly,'r')
    fAm3 = open('results/kde/%s_Am3_PixelResults.txt' %assembly,'r')
    fCuXRF_CERN = open('results/kde/%s_Cu_PixelResults.txt' %assembly,'r')
    fInXRF = open('results/kde/%s_In_PixelResults.txt' %assembly,'r')
    fCo1 = open('results/kde/%s_Co1_PixelResults.txt' %assembly,'r')
    fCo2 = open('results/kde/%s_Co2_PixelResults.txt' %assembly,'r')
    fCd = open('results/kde/%s_Cd_PixelResults.txt' %assembly,'r')

    # LNLS data
    if assembly == "A06-W0110":
        fCoXRF = open('results/kde/%s_CoXRF_PixelResults.txt' %assembly,'r')
        fCrXRF = open('results/kde/%s_CrXRF_PixelResults.txt' %assembly,'r')
        fCuXRF_LNLS = open('results/kde/%s_CuXRF_LNLS_PixelResults.txt' %assembly,'r')
        fFeXRF = open('results/kde/%s_FeXRF_PixelResults.txt' %assembly,'r')
        fMnXRF = open('results/kde/%s_MnXRF_PixelResults.txt' %assembly,'r')
        fNiXRF = open('results/kde/%s_NiXRF_PixelResults.txt' %assembly,'r')
        fTiXRF = open('results/kde/%s_TiXRF_PixelResults.txt' %assembly,'r')
        fVXRF = open('results/kde/%s_VXRF_PixelResults.txt' %assembly,'r')

    chi2ndf2d = []
    a2d = []
    b2d = []
    c2d = []
    t2d = []

    rootf = R.TFile("results/%s_KDECalibration_Pixels_new.root" %assembly, "RECREATE")
    roott = R.TTree("fitPara","")
    pixx = np.zeros(1, dtype=float)
    pixy = np.zeros(1, dtype=float)
    a = np.zeros(1, dtype=float)
    b = np.zeros(1, dtype=float)
    c = np.zeros(1, dtype=float)
    d = np.zeros(1, dtype=float)
    roott.Branch('pixx', pixx, 'pixx/D')
    roott.Branch('pixy', pixy, 'pixy/D')
    roott.Branch('a', a, 'a/D')
    roott.Branch('b', b, 'b/D')
    roott.Branch('c', c, 'c/D')
    roott.Branch('d', d, 'd/D')

    for j in xrange(C.npixX*C.npixY):
        lines = []
        fixed_vals = []

        if j%C.npixX == 0:
            print "==================================", float(j%C.npixX), float(int(j/float(C.npixY)))
            chi2ndf2d.append([])
            a2d.append([])
            b2d.append([])
            c2d.append([])
            t2d.append([])

        if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
            Thline = fTh.readline().split()
            lines.append(Thline)
            fixed_vals.append(0.)

        Feline = fFe.readline().split()
        Am2line = fAm2.readline().split()
        Am3line = fAm3.readline().split()
        CuXRF_CERNline = fCuXRF_CERN.readline().split()
        InXRFline = fInXRF.readline().split()
        Co1line = fCo1.readline().split()
        Co2line = fCo2.readline().split()
        Cdline = fCd.readline().split()
        lines.append(Feline)
        fixed_vals.append(C.FePeakE)
        lines.append(Am2line)
        fixed_vals.append(C.Am2PeakE)
        lines.append(Am3line)
        fixed_vals.append(C.Am3PeakE)
        lines.append(CuXRF_CERNline)
        fixed_vals.append(C.CuXRFPeakE)
        lines.append(InXRFline)
        fixed_vals.append(C.InXRFPeakE)
        lines.append(Co1line)
        fixed_vals.append(C.Co1PeakE)
        lines.append(Co2line)
        fixed_vals.append(C.Co2PeakE)
        lines.append(Cdline)
        fixed_vals.append(C.CdPeakE)

        if assembly == "A06-W0110":
            CoXRFline = fCoXRF.readline().split()
            CrXRFline = fCrXRF.readline().split()
            CuXRF_LNLSline = fCuXRF_LNLS.readline().split()
            FeXRFline = fFeXRF.readline().split()
            MnXRFline = fMnXRF.readline().split()
            NiXRFline = fNiXRF.readline().split()
            TiXRFline = fTiXRF.readline().split()
            VXRFline = fVXRF.readline().split()
            lines.append(CoXRFline)
            fixed_vals.append(C.CoXRFPeakE)
            lines.append(CrXRFline)
            fixed_vals.append(C.CrXRFPeakE)
            lines.append(CuXRF_LNLSline)
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
            #print line
            for i in xrange(len(line)):
                line[i] = ast.literal_eval(line[i])

        tots = array('f',[])
        tot_lerrs = array('f',[])
        tot_uerrs = array('f',[])

        energies = array('f',[])
        energy_lerrs = array('f',[])
        energy_uerrs = array('f',[])

        for line, fixed_val in zip(lines,fixed_vals):
            #print line
            if line[0] == float(j%C.npixX) and line[1] == float(int(j/float(C.npixY))):
                if line == Thline:
                    tots.append(fixed_val)
                    tot_lerrs.append(0.)
                    tot_uerrs.append(0.)
                    
                    energies.append(line[2])
                    energy_lerrs.append(line[3])
                    energy_uerrs.append(line[4])
                else:
                    if line[2] != 0.0:
                        tots.append(line[2])
                        tot_lerrs.append(R.sqrt(line[3]**2 + line[5]**2))
                        tot_uerrs.append(R.sqrt(line[4]**2 + line[5]**2))

                        energies.append(fixed_val)
                        energy_lerrs.append(0.)
                        energy_uerrs.append(0.)
            else:
                print "lines not going as expected"

        if len(tots) < 5:
            chi2ndf2d[-1].append(0.)
            a2d[-1].append(0.)
            b2d[-1].append(0.)
            c2d[-1].append(0.)
            t2d[-1].append(0.)

            pixx[0] = j%C.npixX
            pixy[0] = int(j/float(C.npixY))
            a[0] = 0.
            b[0] = 0.
            c[0] = 0.
            d[0] = 0.

        else:
            canv = R.TCanvas()
            gr = R.TGraphAsymmErrors(len(energies),energies,tots,energy_lerrs,energy_uerrs,tot_lerrs,tot_uerrs)
            gr.SetMarkerStyle(20)
            gr.Draw('AP')

            surrogate = R.TF1("surrogate","[0]*x+[1]-([2]/(x-[3]))",0,60)
            surrogate.SetParName(0,'a')
            surrogate.SetParName(1,'b')
            surrogate.SetParName(2,'c')
            surrogate.SetParName(3,'t')
            # set starting parameters
            surrogate.SetParameters(parEsts[0],parEsts[1],parEsts[2],parEsts[3])
            # limits
            surrogate.SetParLimits(0,parLLims[0],parULims[0])
            surrogate.SetParLimits(1,parLLims[1],parULims[1])
            surrogate.SetParLimits(2,parLLims[2],parULims[2])
            surrogate.SetParLimits(3,parLLims[3],parULims[3])
            gr.Fit("surrogate","RQB") #range, quiet, bounds
            canv.Update()

            stats = gr.GetListOfFunctions().FindObject("stats")
            stats.SetX1NDC(0.5)
            stats.SetX2NDC(0.89)
            stats.SetY1NDC(0.21)
            stats.SetY2NDC(0.46)
            gr.SetMinimum(0.)
            gr.SetMaximum(surrogate.Eval(65.0))
            canv.Update()

            if j%C.npixX == 0 and int(j/float(C.npixY)) < 10:
                canv.SaveAs("plots/KDESurrogateFits/%s_examplefit_%i_new.pdf" %(assembly,int(j/float(C.npixY))))

            chi2ndf2d[-1].append(surrogate.GetChisquare() / surrogate.GetNDF())
            a2d[-1].append(surrogate.GetParameter(0))
            b2d[-1].append(surrogate.GetParameter(1))
            c2d[-1].append(surrogate.GetParameter(2))
            t2d[-1].append(surrogate.GetParameter(3))

            pixx[0] = j%C.npixX
            pixy[0] = int(j/float(C.npixY))
            a[0] = surrogate.GetParameter(0)
            b[0] = surrogate.GetParameter(1)
            c[0] = surrogate.GetParameter(2)
            d[0] = surrogate.GetParameter(3)

        roott.Fill()

    if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
        fTh.close()
    fFe.close()
    fAm2.close()
    fAm3.close()
    fCuXRF_CERN.close()
    fInXRF.close()
    fCo1.close()
    fCo2.close()
    fCd.close()
    if assembly == "A06-W0110":
        fCoXRF.close()
        fCrXRF.close()
        fCuXRF_LNLS.close()
        fFeXRF.close()
        fMnXRF.close()
        fNiXRF.close()
        fTiXRF.close()
        fVXRF.close()

    print "writing tree"
    rootf.Write()
    rootf.Close()
    print "tree done, making 1D plots"

    # 1D plots
    a1d = list(itertools.chain(*a2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted a')
    ax.hist(a1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.savefig("plots/KDESurrogateFits/%s_a_hist_new.pdf" %assembly)

    b1d = list(itertools.chain(*b2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted b')
    ax.hist(b1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.savefig("plots/KDESurrogateFits/%s_b_hist_new.pdf" %assembly)

    c1d = list(itertools.chain(*c2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted c')
    ax.hist(c1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/4.))
    plt.savefig("plots/KDESurrogateFits/%s_c_hist_new.pdf" %assembly)

    t1d = list(itertools.chain(*t2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted t')
    ax.hist(t1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.savefig("plots/KDESurrogateFits/%s_t_hist_new.pdf" %assembly)

    chi2ndf1d = list(itertools.chain(*chi2ndf2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted chi2')
    ax.hist(chi2ndf1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.savefig("plots/KDESurrogateFits/%s_chi2ndf_hist_new.pdf" %assembly)

    print "making 2D plots"
    #2D
    dx, dy = 1.0, 1.0
    y, x = np.mgrid[slice(0, C.npixY + dy, dy),slice(0, C.npixX + dx, dx)]

    chi2ndf2d = np.ma.masked_equal(chi2ndf2d,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, chi2ndf2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    plt.savefig("plots/KDESurrogateFits/%s_chi2ndf_map_nz_new.pdf" %assembly)

    a2d = np.ma.masked_equal(a2d,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, a2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    plt.savefig("plots/KDESurrogateFits/%s_a_map_nz_new.pdf" %assembly)

    b2d = np.ma.masked_equal(b2d,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, b2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    plt.savefig("plots/KDESurrogateFits/%s_b_map_nz_new.pdf" %assembly)

    c2d = np.ma.masked_equal(c2d,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, c2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    plt.savefig("plots/KDESurrogateFits/%s_c_map_nz_new.pdf" %assembly)

    t2d = np.ma.masked_equal(t2d,0)
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X')
    ax.set_ylabel('Pixel Y')
    myplot = plt.pcolor(x, y, t2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    plt.savefig("plots/KDESurrogateFits/%s_t_map_nz_new.pdf" %assembly)


def getParEstLim(assembly):

    if assembly == "A06-W0110":
        parEsts = [14.0, 350.0, 2000.0, -6.0]
        parLLims = [0.0, 0.0, 0.0, -200.0]
        parULims = [30.0, 3000.0, 200000, 50.0]

    if assembly == "B06-W0125":
        parEsts = [30.0, 600.0, 4000.0, -4.0]
        parLLims = [0.0, 0.0, 0.0, -100.0]
        parULims = [50.0, 5000.0, 200000.0, 20.0]

    if assembly == "B07-W0125":
        parEsts = [14.0, 400.0, 500.0, 4.0]
        parLLims = [0.0, 0.0, 0.0, -500.0]
        parULims = [50.0, 5000.0, 400000.0, 50.0]

    if assembly == "C04-W0110":
        parEsts = [13.0, 400.0, 2000.0, 0.0]
        parLLims = [4.0, 0.0, 0.0, -10.0]
        parULims = [20.0, 3000.0, 20000, 10.0]

    if assembly == "D09-W0126":
        parEsts = [18.0, 400.0, 1000.0, 0.0]
        parLLims = [0.0, 0.0, 0.0, -30.0]
        parULims = [30.0, 1000.0, 25000.0, 10.0]

    if assembly == "L04-W0125":
        parEsts = [14.0, 500.0, 5000.0, -6.0]
        parLLims = [0.0, 0.0, 0.0, -70.0]
        parULims = [30.0, 2000.0, 100000.0, 10.0]
    
    return parEsts, parLLims, parULims


parEsts, parLLims, parULims = getParEstLim(assembly)
fitPixelSurrogates(assembly,parEsts,parLLims,parULims)

