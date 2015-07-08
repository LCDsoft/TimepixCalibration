# will use kde method to finding most likely TOT
# for the global calibration

from optparse import OptionParser
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import ROOT as R
from os import environ
import getpass
import Constants as C

parser = OptionParser()
parser.add_option("-b", "--assembly",
                  help="Assembly name", dest="ASSEMBLY")

parser.add_option("-s", "--source",
                  help="Source name", dest="SOURCE")

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

if(options.SOURCE):
    source=options.SOURCE
else :
    print "Please specify source"
    print "choose from", C.known_sources
    parser.print_help()
    exit()

if source not in C.known_sources:
    print "Source not recognised"
    print "choose from", C.known_sources
    exit()

if source in C.LNLS_sources and assembly != "A06-W0110":
    print "Source only available for assembly A06-W0110"
    print "please reconsider input"
    exit()

# wanted to use kde_scipy function, but wouldn't work wth my scipy version
# this is a workaround which does the same


def findMostLikelyTOT(assembly,source,llim,ulim):

    # Load data    
    home = environ['HOME']
    base = "%s/eos/clicdp/data/VertexCalibration" %home
    assembly_start = assembly.split("-")[0]

    if source == "Fe":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Fe55_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Am":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Am241_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "Cd":
        rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cd109_%s_spc.root"%(base,assembly,assembly_start,assembly))
    elif source == "CuInXRF":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Cu_In_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/CuIn_%s.root" %(base,assembly,assembly))
    elif source == "Co":
        if assembly == "B06-W0125":
            rootfile = R.TFile("%s/%s/%s_SinglePixelCalibration/Co57_%s_spc.root"%(base,assembly,assembly_start,assembly))
        else:
            rootfile = R.TFile("%s/%s/Co57_%s.root" %(base,assembly,assembly))
    elif source == "CoXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CoXRF_CalibTree.root" %(base,assembly))
    elif source == "CrXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CrXRF_CalibTree.root" %(base,assembly))
    elif source == "CuXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_CuXRF_CalibTree.root" %(base,assembly))
    elif source == "FeXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_FeXRF_CalibTree.root" %(base,assembly))
    elif source == "MnXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_MnXRF_CalibTree.root" %(base,assembly))
    elif source == "NiXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_NiXRF_CalibTree.root" %(base,assembly))
    elif source == "TiXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_TiXRF_CalibTree.root" %(base,assembly))
    elif source == "VXRF":
        rootfile = R.TFile("%s/LNLS_Analysis/SinglePixelAnalysis/root_files/%s-25V_VXRF_CalibTree.root" %(base,assembly))

    tree = rootfile.Get("pixels")
    print "got tree"

    # Set up junk file to appease ROOT
    username = getpass.getuser()
    junkfile = R.TFile("/tmp/%s/junkfile_%s_%s.root" %(username,assembly,source),"RECREATE")

    # Just keep the events we need
    t2 = tree.CopyTree("tot < %i && tot > %i" %(ulim,llim))
    print "copied tree"

    peak_tots = []
    peak_amps = []
    peak_is = []
    peak_loweris = []
    peak_upperis = []
    peak_ents = []
    peak_lowersigmas = []
    peak_uppersigmas = []

    nsteps = 100
    x_grid = np.linspace(llim, ulim, nsteps)
    step_size = (ulim-llim)/float(nsteps-1)
    ent = t2.GetEntries()
    print "Will load", ent, "entries"
    tot = []
    if ent > 1:
        tot = np.zeros(ent, dtype=float)
        t2.Branch('tot', tot, 'tot/F')
        for i in xrange(ent):
            if i%1000000==0:
                print ".....loading", i
            t2.GetEvent(i)
            tot[i] = t2.tot

        try:
            print "trying"
            # The KDE calculation
            density = gaussian_kde(tot)
            if density.silverman_factor() > 0.1:
                density.covariance_factor = density.silverman_factor
                print "Bandwidth determined from Silverman factor:", density.silverman_factor()
            else:
                density.covariance_factor = lambda: 0.1
                print "Bandwidth set at 0.1 (Silverman factor too small: %f)" %density.silverman_factor()
            density._compute_covariance()
            workaround = density(x_grid)

            # Find peaks by finding where gradient passes through 0
            grad = np.gradient(workaround)
            last_grad = grad[0]
            for i in xrange(1,len(grad)):
                this_grad = grad[i]
                if (last_grad > 0 and this_grad < 0):
                    peak_tots.append((x_grid[i-1] + x_grid[i]) / 2.)
                    peak_amps.append((workaround[i-1] + workaround[i]) / 2.)
                    peak_is.append(i)
                last_grad = this_grad
                
            for i in xrange(len(peak_tots)):
                print "Peak found at TOT", peak_tots[i], "with amplitude", peak_amps[i], "at xgrid", peak_is[i]

            # Filter peaks to take only highest
            print "pre filter", peak_amps,peak_tots,peak_is
            sorted_peaks = sorted(zip(peak_amps,peak_tots,peak_is),reverse=1)
            if source in ["Fe","Cd","CoXRF","CrXRF","CuXRF","FeXRF","MnXRF","NiXRF","TiXRF","VXRF"]:
                peak_amps = [sorted_peaks[0][0]]
                peak_tots = [sorted_peaks[0][1]]
                peak_is = [sorted_peaks[0][2]]
            if source in ["Co","CuInXRF"]:
                peak_amps = [sorted_peaks[0][0],sorted_peaks[1][0]]
                peak_tots = [sorted_peaks[0][1],sorted_peaks[1][1]]
                peak_is = [sorted_peaks[0][2],sorted_peaks[1][2]]
            if source == "Am":
                peak_amps = [sorted_peaks[0][0],sorted_peaks[1][0],sorted_peaks[2][0]]
                peak_tots = [sorted_peaks[0][1],sorted_peaks[1][1],sorted_peaks[2][1]]
                peak_is = [sorted_peaks[0][2],sorted_peaks[1][2],sorted_peaks[2][2]]
            print "post filter", peak_amps,peak_tots,peak_is

            # Put peaks back on TOT order
            ordered_peaks = sorted(zip(peak_tots,peak_amps,peak_is))
            peak_tots = [peaki[0] for peaki in ordered_peaks]
            peak_amps = [peaki[1] for peaki in ordered_peaks]
            peak_is = [peaki[2] for peaki in ordered_peaks]
            print "post ordering", peak_amps,peak_tots,peak_is

            # Calculate the upper and lower i for each peak
            for i in xrange(len(peak_tots)):
                if i == 0 and i == len(peak_tots)-1:
                    peak_loweris.append(0)
                    peak_upperis.append(nsteps)
                elif i == 0:
                    peak_loweris.append(0)
                    peak_upperis.append((peak_is[i] + peak_is[i+1])/2.)
                elif i == len(peak_tots)-1:
                    peak_loweris.append((peak_is[i-1] + peak_is[i])/2.)
                    peak_upperis.append(nsteps)
                else:
                    peak_loweris.append((peak_is[i-1] + peak_is[i])/2.)
                    peak_upperis.append((peak_is[i] + peak_is[i+1])/2.)
            print "peak_is", peak_is
            print "peak_loweris", peak_loweris
            print "peak_upperis", peak_upperis

            # Calculate the entries in each peak
            for i in xrange(len(peak_tots)):
                if i == 0 and i == len(peak_tots)-1:
                    peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i" %(llim,ulim))))
                elif i == 0:
                    peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i" %(llim,(peak_tots[i] + peak_tots[i+1])/2.))))
                elif i == len(peak_tots)-1:
                    peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i" %((peak_tots[i-1] + peak_tots[i])/2.,ulim))))
                else:
                    peak_ents.append(int(t2.GetEntries("tot > %i && tot < %i" %((peak_tots[i-1] + peak_tots[i])/2.,(peak_tots[i] + peak_tots[i+1])/2.))))
            print "peak_ents =", peak_ents

            # Calculate the uncertainty on each peak position
            for i in xrange(len(peak_tots)):
                maxindex = peak_is[i]
                lowerindex = peak_is[i]-1
                while (np.trapz(workaround[lowerindex:maxindex],x=x_grid[lowerindex:maxindex]) < 0.341*2*np.trapz(workaround[peak_loweris[i]:maxindex],x=x_grid[peak_loweris[i]:maxindex])) and (lowerindex>0):
                    lowerindex=lowerindex-1
                lowersigma = (x_grid[maxindex]-x_grid[lowerindex]) / R.sqrt(peak_ents[i])
                peak_lowersigmas.append(lowersigma)

                maxindex = peak_is[i]+1
                upperindex = peak_is[i]+2
                while (np.trapz(workaround[maxindex:upperindex],x=x_grid[maxindex:upperindex]) < 0.341*2*np.trapz(workaround[maxindex:peak_upperis[i]],x=x_grid[maxindex:peak_upperis[i]])) and (upperindex<(nsteps-1)):
                    upperindex=upperindex+1
                uppersigma = (x_grid[upperindex]-x_grid[maxindex]) / R.sqrt(peak_ents[i])
                peak_uppersigmas.append(uppersigma)

        except np.linalg.linalg.LinAlgError as err:
            peak_tots.append(0.)
            peak_lowersigmas.append(0.)
            peak_uppersigmas.append(0.)

    else:
        peak_tots.append(0.)
        peak_lowersigmas.append(0.)
        peak_uppersigmas.append(0.)

    # Make plots
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    plt.xticks(np.arange(llim,ulim+1,(ulim-llim)/5.))
    ax.set_xlabel('TOT [1/96 MHz]', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('Pixels / Total number of pixels', x=1, y=1, verticalalignment='top')
    for i in xrange(len(peak_tots)):
        ax.text(0.01, 0.99 - (i*0.1), r'Peak: $%i \pm ^{%0.2f} _{%0.2f} \pm %0.2f$' %(peak_tots[i],peak_uppersigmas[i],peak_lowersigmas[i],step_size/2.),
                verticalalignment='top', horizontalalignment='left',
                transform=ax.transAxes,
                fontsize=35)
    ax.hist(tot, bins=100,fc='gray',alpha=0.3,normed=True)
    ax.plot(x_grid, workaround, color='blue', alpha=0.5, lw=3)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fontP = FontProperties()
    fontP.set_size(30)
    blue_line = mlines.Line2D([], [], color='blue', lw=3)
    grey_box = mpatches.Patch(color='gray', alpha=0.3)
    fig.legend(handles=[blue_line, grey_box], labels=['KDE PDF', 'Histogram'], loc=(0.63, 0.85), prop=fontP, frameon=False)
    fig.tight_layout()
    fig.savefig("plots/KDEPeaks/Global/%s_%s_GlobalSpectrum.pdf" %(assembly,source))

    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.tick_params(axis='x', pad=20)
    plt.xticks(np.arange(llim,ulim+1,(ulim-llim)/5.))
    ax.set_xlabel('TOT [1/96 MHz]', x=1, y=1, horizontalalignment='right')
    ax.set_ylabel('$\delta$ KDE PDF / $\delta$ TOT', x=1, y=1, verticalalignment='top')
    ax.plot(x_grid, grad, color='red', alpha=0.5, lw=3)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(40)
    fig.tight_layout()
    fig.savefig("plots/KDEPeaks/Global/%s_%s_GlobalSpectrumDeriv.pdf" %(assembly,source))

    # Write results to txt file
    f = open('results/kde/%s_%s_GlobalResults.txt' %(assembly,source), 'w')
    if len(peak_tots) == 1:
        f.write('%f \t %f \t %f \t %f \n' %(peak_tots[0],peak_lowersigmas[0],peak_uppersigmas[0],step_size/2.))
    elif len(peak_tots) == 2:
        f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \n' %(peak_tots[0],peak_lowersigmas[0],peak_uppersigmas[0],peak_tots[1],peak_lowersigmas[1],peak_uppersigmas[1],step_size/2.))
    else:
        f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' %(peak_tots[0],peak_lowersigmas[0],peak_uppersigmas[0],peak_tots[1],peak_lowersigmas[1],peak_uppersigmas[1],peak_tots[2],peak_lowersigmas[2],peak_uppersigmas[2],step_size/2.))
    f.close()
    print "finished", assembly, source


def getLimits(assembly,source):

    if assembly == "A06-W0110":
        if source == "Fe": limits = [0,400]
        if source == "Am": limits = [0,1500]
        if source == "Cd": limits = [0,800]
        if source == "CuInXRF": limits = [0,1000]
        if source == "Co": limits = [0,700]
        if source == "CoXRF": limits = [0,400]
        if source == "CrXRF": limits = [0,400]
        if source == "CuXRF": limits = [0,400]
        if source == "FeXRF": limits = [0,400]
        if source == "MnXRF": limits = [0,400]
        if source == "NiXRF": limits = [0,400]
        if source == "TiXRF": limits = [0,400]
        if source == "VXRF": limits = [0,400]

    if assembly == "B06-W0125":
        if source == "Fe": limits = [0,700]
        if source == "Am": limits = [0,2800]
        if source == "Cd": limits = [0,1700]
        if source == "CuInXRF": limits = [0,1700]
        if source == "Co": limits = [0,1100]

    if assembly == "B07-W0125":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,1500]
        if source == "Cd": limits = [0,900]
        if source == "CuInXRF": limits = [0,1000]
        if source == "Co": limits = [0,600]

    if assembly == "C04-W0110":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,1400]
        if source == "Cd": limits = [0,800]
        if source == "CuInXRF": limits = [0,700]
        if source == "Co": limits = [0,600]

    if assembly == "D09-W0126":
        if source == "Fe": limits = [0,600]
        if source == "Am": limits = [0,1800]
        if source == "Cd": limits = [0,1200]
        if source == "CuInXRF": limits = [0,1100]
        if source == "Co": limits = [0,800]

    if assembly == "L04-W0125":
        if source == "Fe": limits = [0,500]
        if source == "Am": limits = [0,2000]
        if source == "Cd": limits = [0,1100]
        if source == "CuInXRF": limits = [0,1200]
        if source == "Co": limits = [0,700]

    return limits


limits = getLimits(assembly,source)
findMostLikelyTOT(assembly,source,limits[0],limits[1])

