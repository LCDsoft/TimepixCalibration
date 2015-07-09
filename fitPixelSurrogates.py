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


def fitPixelSurrogates(assembly,parEsts,parLLims,parULims,globalPars):

    # thresholds
    Thline = 0
    if assembly == "A06-W0110" or assembly == "B06-W0125" or assembly == "C04-W0110" or assembly == "L04-W0125":
        fTh = open('results/threshold/%s_PixelThreshold.txt' %assembly,'r')

    # CERN data
    fFe = open('results/kde/%s_Fe_PixelResults.txt' %assembly,'r')
    fAm = open('results/kde/%s_Am_PixelResults.txt' %assembly,'r')
    fCuInXRF = open('results/kde/%s_CuInXRF_PixelResults.txt' %assembly,'r')
    fCo = open('results/kde/%s_Co_PixelResults.txt' %assembly,'r')
    fCd = open('results/kde/%s_Cd_PixelResults.txt' %assembly,'r')

    # LNLS data
    if assembly == "A06-W0110":
        fCoXRF = open('results/kde/%s_CoXRF_PixelResults.txt' %assembly,'r')
        fCrXRF = open('results/kde/%s_CrXRF_PixelResults.txt' %assembly,'r')
        fCuXRF = open('results/kde/%s_CuXRF_PixelResults.txt' %assembly,'r')
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
    n_noFitPixels = 0

    rootf = R.TFile("results/%s_KDECalibration_Pixels.root" %assembly, "RECREATE")
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
        Amline = fAm.readline().split()
        CuInXRFline = fCuInXRF.readline().split()
        Coline = fCo.readline().split()
        Cdline = fCd.readline().split()
        lines.append(Feline)
        fixed_vals.append(C.FePeakE)
        lines.append([Amline[0],Amline[1],Amline[2],Amline[3],Amline[4],Amline[8]])
        fixed_vals.append(C.Am2PeakE)
        lines.append([Amline[0],Amline[1],Amline[5],Amline[6],Amline[7],Amline[8]])
        fixed_vals.append(C.Am3PeakE)
        lines.append([CuInXRFline[0],CuInXRFline[1],CuInXRFline[2],CuInXRFline[3],CuInXRFline[4],CuInXRFline[8]])
        fixed_vals.append(C.CuXRFPeakE)
        lines.append([CuInXRFline[0],CuInXRFline[1],CuInXRFline[5],CuInXRFline[6],CuInXRFline[7],CuInXRFline[8]])
        fixed_vals.append(C.InXRFPeakE)
        lines.append([Coline[0],Coline[1],Coline[2],Coline[3],Coline[4],Coline[8]])
        fixed_vals.append(C.Co1PeakE)
        lines.append([Coline[0],Coline[1],Coline[5],Coline[6],Coline[7],Coline[8]])
        fixed_vals.append(C.Co2PeakE)
        lines.append(Cdline)
        fixed_vals.append(C.CdPeakE)

        if assembly == "A06-W0110":
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
            if line[0] == float(j%C.npixX) and line[1] == float(int(j/float(C.npixY))):
                if line == Thline:
                    tots.append(fixed_val)
                    tot_lerrs.append(R.sqrt(5.0**2))
                    tot_uerrs.append(R.sqrt(5.0**2))
                    
                    energies.append(line[2])
                    energy_lerrs.append(line[3] + line[2]*0.03)
                    energy_uerrs.append(line[4] + line[2]*0.03)
                else:
                    if line[2] != 0.0:
                        tots.append(line[2])
                        tot_lerrs.append(R.sqrt(line[3]**2 + line[5]**2) + R.sqrt((line[2]*0.03)**2 + 5.0**2))
                        tot_uerrs.append(R.sqrt(line[4]**2 + line[5]**2) + R.sqrt((line[2]*0.03)**2 + 5.0**2))

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

            n_noFitPixels = n_noFitPixels+1

        else:
            canv = R.TCanvas()
            canv.SetRightMargin(0.01)
            canv.SetLeftMargin(0.19)
            gr = R.TGraphAsymmErrors(len(energies),energies,tots,energy_lerrs,energy_uerrs,tot_lerrs,tot_uerrs)
            gr.SetMarkerStyle(20)
            gr.GetXaxis().SetTitle("Energy [keV]")
            gr.GetYaxis().SetTitle("TOT [1/96 MHz]")
            gr.GetYaxis().SetTitleOffset(1.4)
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
            stats.SetX1NDC(0.57)
            stats.SetX2NDC(0.96)
            stats.SetY1NDC(0.21)
            stats.SetY2NDC(0.46)
            stats.SetBorderSize(0)
            gr.SetMinimum(0.)
            gr.SetMaximum(surrogate.Eval(65.0))
            canv.Update()
            if j%C.npixX == 0 and int(j/float(C.npixY)) < 10:
                canv.SaveAs("plots/KDESurrogateFits/%s_examplefit_%i.pdf" %(assembly,int(j/float(C.npixY))))

            chi2ndf2d[-1].append(surrogate.GetChisquare() / surrogate.GetNDF())
            a2d[-1].append(surrogate.GetParameter(0))
            b2d[-1].append(surrogate.GetParameter(1))
            c2d[-1].append(surrogate.GetParameter(2))
            t2d[-1].append(surrogate.GetParameter(3))


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

    # make chi2 1d plot to find chi2/ndf cut
    chi2ndf1d = list(itertools.chain(*chi2ndf2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted chi2', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    n, bins, patches = ax.hist(chi2ndf1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    # centre of highest bin
    peak_chi2ndf = (bins[np.argmax(n)]+bins[np.argmax(n)+1])/2.
    std_chi2ndf = np.std(chi2ndf1d)
    chi2ndf_cut = peak_chi2ndf + 3*std_chi2ndf

    # edit 2d arrays to remove bad fits
    over_chi2threshold = 0
    for i in xrange(256):
        for j in xrange(256):
            if chi2ndf2d[i][j] > chi2ndf_cut:
                over_chi2threshold = over_chi2threshold+1
                chi2ndf2d[i][j] = 0.
                a2d[i][j] = 0.
                b2d[i][j] = 0.
                c2d[i][j] = 0.
                t2d[i][j] = 0.

    # replace zeros with the constants of the neighbour (if they exist)
    # use the chi2 to check a replaced value isn't propagated
    n_replaced = 0
    for i in xrange(256):
        for j in xrange(256):
            if a2d[i][j] == 0. and b2d[i][j] == 0. and c2d[i][j] == 0. and t2d[i][j] == 0. :
                n_replaced = n_replaced + 1
                if i > 0 and a2d[i-1][j] != 0. and chi2ndf2d[i-1][j] != 0. :
                    a2d[i][j] = a2d[i-1][j]
                    b2d[i][j] = b2d[i-1][j]
                    c2d[i][j] = c2d[i-1][j]
                    t2d[i][j] = t2d[i-1][j]
                elif i < 255 and a2d[i+1][j] != 0. and chi2ndf2d[i+1][j] != 0. :
                    a2d[i][j] = a2d[i+1][j]
                    b2d[i][j] = b2d[i+1][j]
                    c2d[i][j] = c2d[i+1][j]
                    t2d[i][j] = t2d[i+1][j]
                elif j > 0 and a2d[i][j-1] != 0. and chi2ndf2d[i][j-1] != 0. :
                    a2d[i][j] = a2d[i][j-1]
                    b2d[i][j] = b2d[i][j-1]
                    c2d[i][j] = c2d[i][j-1]
                    t2d[i][j] = t2d[i][j-1]
                elif j < 255 and a2d[i][j+1] != 0. and chi2ndf2d[i][j+1] != 0. :
                    a2d[i][j] = a2d[i][j+1]
                    b2d[i][j] = b2d[i][j+1]
                    c2d[i][j] = c2d[i][j+1]
                    t2d[i][j] = t2d[i][j+1]
                else:
                    a2d[i][j] = globalPars[0]
                    b2d[i][j] = globalPars[1]
                    c2d[i][j] = globalPars[2]
                    t2d[i][j] = globalPars[3]

    # make the root tree
    print "making the root tree"
    for i in xrange(256):
        for j in xrange(256):
            pixx[0] = j
            pixy[0] = i
            a[0] = a2d[i][j]
            b[0] = b2d[i][j]
            c[0] = c2d[i][j]
            d[0] = t2d[i][j]
            roott.Fill()
    rootf.Write()
    rootf.Close()

    # save plots
    # 1D plots
    print "making 1D plots"
    chi2ndf1d = list(itertools.chain(*chi2ndf2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted $\chi^{2}$ / NDF', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    n, bins, patches = ax.hist(chi2ndf1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.plot((chi2ndf_cut, chi2ndf_cut), (limits[2], limits[3]), 'r-')
    ax.axis(limits)
    fig.tight_layout()
    fig.savefig("plots/KDESurrogateFits/%s_chi2ndf_hist.pdf" %assembly)

    a1d = list(itertools.chain(*a2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted a', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    ax.hist(a1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.plot((globalPars[0], globalPars[0]), (limits[2], limits[3]), 'b-')
    ax.axis(limits)
    fig.tight_layout()
    fig.savefig("plots/KDESurrogateFits/%s_a_hist.pdf" %assembly)

    b1d = list(itertools.chain(*b2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted b', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    ax.hist(b1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.plot((globalPars[1], globalPars[1]), (limits[2], limits[3]), 'b-')
    ax.axis(limits)
    fig.tight_layout()
    fig.savefig("plots/KDESurrogateFits/%s_b_hist.pdf" %assembly)

    c1d = list(itertools.chain(*c2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted c', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    ax.hist(c1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/4.))
    plt.plot((globalPars[2], globalPars[2]), (limits[2], limits[3]), 'b-')
    ax.axis(limits)
    fig.tight_layout()
    fig.savefig("plots/KDESurrogateFits/%s_c_hist.pdf" %assembly)

    t1d = list(itertools.chain(*t2d))
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    ax.set_xlabel('Fitted t', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    ax.hist(t1d, 100, fc='gray', alpha=0.3,normed=True)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    limits = ax.axis()
    plt.xticks(np.arange(limits[0],limits[1]+1,(limits[1]-limits[0])/5.))
    plt.plot((globalPars[3], globalPars[3]), (limits[2], limits[3]), 'b-')
    ax.axis(limits)
    fig.tight_layout()
    fig.savefig("plots/KDESurrogateFits/%s_t_hist.pdf" %assembly)

    print "making 2D plots"
    #2D
    dx, dy = 1.0, 1.0
    y, x = np.mgrid[slice(0, C.npixY + dy, dy),slice(0, C.npixX + dx, dx)]

    chi2ndf2d = np.ma.masked_equal(chi2ndf2d,0) # mask zeros
    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixel Y', x=1, y=1, verticalalignment='top')
    myplot = plt.pcolor(x, y, chi2ndf2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDESurrogateFits/%s_chi2ndf_map_nz.pdf" %assembly)

    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixel Y', x=1, y=1, verticalalignment='top')
    myplot = plt.pcolor(x, y, a2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDESurrogateFits/%s_a_map.pdf" %assembly)

    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixel Y', x=1, y=1, verticalalignment='top')
    myplot = plt.pcolor(x, y, b2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDESurrogateFits/%s_b_map.pdf" %assembly)

    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixel Y', x=1, y=1, verticalalignment='top')
    myplot = plt.pcolor(x, y, c2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDESurrogateFits/%s_c_map.pdf" %assembly)

    fig, ax = plt.subplots(1,1,figsize=(12, 10))
    ax.set_xlabel('Pixel X', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixel Y', x=1, y=1, verticalalignment='top')
    myplot = plt.pcolor(x, y, t2d, cmap='jet')
    cbar = fig.colorbar(myplot)
    cbar.ax.tick_params(labelsize=40) 
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.savefig("plots/KDESurrogateFits/%s_t_map.pdf" %assembly)

    print "Finished assembly %s" %assembly
    print "Peak chi2ndf", peak_chi2ndf
    print "Standard deviation of chi2ndf", std_chi2ndf
    print "So chi2ndf cut at", chi2ndf_cut
    print "Pixels not fit due to no/not enough data points:", n_noFitPixels
    print "Pixels removed due to bad chi2ndf:", over_chi2threshold
    print "Pixels replaced:", n_replaced

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
        parULims = [30.0, 5000.0, 400000.0, 50.0]

    if assembly == "C04-W0110":
        parEsts = [13.0, 400.0, 2000.0, 0.0]
        parLLims = [4.0, 0.0, 0.0, -10.0]
        parULims = [20.0, 3000.0, 20000, 10.0]

    if assembly == "D09-W0126":
        parEsts = [18.0, 500.0, 1000.0, 0.0]
        parLLims = [0.0, 0.0, 0.0, -30.0]
        parULims = [30.0, 5000.0, 25000.0, 10.0]

    if assembly == "L04-W0125":
        parEsts = [14.0, 500.0, 5000.0, -6.0]
        parLLims = [0.0, 0.0, 0.0, -70.0]
        parULims = [30.0, 2000.0, 100000.0, 10.0]
    
    return parEsts, parLLims, parULims

def getGlobalPar(assembly):

    if assembly == "A06-W0110":
        globalPars = [12.76,399.3,2104.0,-1.663]

    if assembly == "B06-W0125":
        globalPars = [30.8,484.3,1301.0,1.65]

    if assembly == "B07-W0125":
        globalPars = [14.1,405.7,2148.0,-1.408]

    if assembly == "C04-W0110":
        globalPars = [14.56,289.4,869.3,0.4814]

    if assembly == "D09-W0126":
        globalPars = [17.49,449.6,1132.0,1.727]

    if assembly == "L04-W0125":
        globalPars = [15.36,414.0,2026.0,-1.043]

    return globalPars

parEsts, parLLims, parULims = getParEstLim(assembly)
globalPars = getGlobalPar(assembly)
fitPixelSurrogates(assembly,parEsts,parLLims,parULims,globalPars)

