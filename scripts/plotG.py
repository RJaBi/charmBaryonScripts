"""
Plots multiple correaltors on same plot
Uses subplots
"""


import numpy as np  # type: ignore
import os
import sys
# import pickle
import toml
import gvar as gv  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.ticker import MaxNLocator  # type: ignore
from matplotlib.container import ErrorbarContainer  # type: ignore
import matplotlib as mpl  # type: ignore
# Grabbing some modules from elsewhere
from pathlib import Path
insertPath = os.path.join(Path(__file__).resolve().parent, '..', 'lib')
sys.path.insert(0, insertPath)
import myModules as mo  # type: ignore # noqa: E402
import myIO as io  # type: ignore # noqa: E402
import myGVPlotter as GVP  # type: ignore # noqa: E402


def main(args: list):
    """
    Reads some parameters in from the command line
    Puts that toml file in
    then reads it. Does a bunch of analysis
    """

    params = mo.GetArgs(args)
    # Printing the input toml back out - easier to read than printing as a dictionary
    toml.dump(params, f=sys.stdout)
    # Loading correlators
    cfDict, cfLabelsList = io.initCorrelators(params)

    # Set up where to save
    anaDir = os.path.join(params['analysis']['anaDir'])
    print('Analysis output directory is ', anaDir)
    if not os.path.exists(anaDir):
        os.makedirs(anaDir)

    # Get the subplots set up
    numA = params['analysis']['plotLayout']
    # Sharing x-axis or y-axis
    sharex = True
    if numA[1] > 2:
        sharey = True
    else:
        sharey = False
    print(f'sharex = {sharex}, sharey = {sharey}')
    fig, ax = plt.subplots(numA[0], numA[1], figsize=(16.6, 11.6), sharey=sharey, sharex=sharey, gridspec_kw={'hspace': 0, 'wspace': 0})  # noqa: E501
    setLegend = False
    # Counter for the analysis
    aa = 0
    print('LOOPING OVER THE SUBPLOTS')
    # Now iterate over the different sets
    for hh in range(0, numA[0]):
        # hh = vertical
        # vv = horizontal
        for vv in range(0, numA[1]):
            print(f'hh = {hh}, vv = {vv}, aa = {aa}')
            # Handling 1-dim and 2-dim
            if numA[0] > 1 and numA[1] > 1:
                thisAX = ax[hh][vv]
            elif numA[0] == 1 and numA[1] > 1:
                thisAX = ax[vv]
            elif numA[0] > 1 and numA[1] == 1:
                thisAX = ax[hh]
            else:
                sys.exit(f'shouldnt get here. bad numA {numA}')
            # Setting the ticks to only be at these points
            thisAX.get_xaxis().set_ticks([4, 8, 12, 16, 20, 24, 28])
            # And changing how the ticks look
            thisAX.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='on',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                direction='in',
                labelbottom='off'  # labels along the bottom edge are off)
            )
            # Set the axis labels
            if hh == numA[0] - 1:
                thisAX.set_xlabel('$\\tau/a_\\tau$')
                thisAX.xaxis.set_major_locator(MaxNLocator(integer=True, prune='upper'))
            if vv == 0:
                thisAX.set_ylabel('$\\overline{G}(\\tau)$')

            # Setting some plot details
            thisAX.set_xlim(params['analysis']['GxLim'])
            thisAX.set_yscale('log')
            # Integer only on x-axis markers
            # Get details for this plot
            aName = params['analysis']['anaName'][aa]
            labels = params['analysis']['labelsToAnalyse'][aa]
            print(aa, aName, labels)
            aLab = params['eMass'][aName]['name']
            J = params['eMass'][aName]['J']
            markCount = 0
            for lab in labels:
                if 'markers' in params['analysis'].keys():
                    if markCount > len(params['analysis']['markers']):
                        markCount = 0
                    mark = params['analysis']['markers'][markCount]
                    markCount = markCount + 1
                else:
                    mark = 'd'
                GVD = gv.dataset.avg_data(io.labelNorm(cfDict['data'][lab], params))
                parity = '+'
                if 'F_' in lab:
                    parity = '-'
                # Setting the text label
                text = '${}^{' + J + '^' + parity + '}' + aLab + '$'
                x = np.array(range(0, len(GVD)))
                # Take only up to min of correlator
                minArg = np.argmin(np.fabs(GVD))
                # xPlot = int(np.ceil(len(x)/2)) - 1
                xPlot = int(minArg)  # I don't know why this isn't already an int...
                print(lab, minArg)
                # take only up to the xlimit on plot
                xPlotEnd = params['analysis']['GxLim'][1] + 1
                # essentially selecting up to min of xPlot and xPlotEnd
                x = x[:xPlot][:xPlotEnd]
                # Doing the same for y.
                # Doing as for x doesn't work (cause gvar?)
                if xPlot > xPlotEnd:
                    y = GVD[:xPlotEnd]
                else:
                    y = GVD[:xPlot]
                thisAX = GVP.plot_gvEbar(x, y, thisAX, ma=mark, lab=lab)
                # This sets horizontally adjacent to have same ylimits
                if not sharey:
                    if vv > 0:
                        prevYLim = ax[hh][vv-1].get_ylim()
                        curYLim = thisAX.get_ylim()
                        yLimMin = np.min([prevYLim[0], curYLim[0]])
                        yLimMax = np.max([prevYLim[1], curYLim[1]])
                        thisAX.set_ylim([yLimMin, yLimMax])
                        ax[hh][vv-1].set_ylim([yLimMin, yLimMax])
            print(text)
            # And put the text on the figure
            thisAX.text(0.25, 0.2, text, horizontalalignment='center', verticalalignment='center', transform=thisAX.transAxes)  # noqa: E501
            # Setting plots on RHS to have axis on right
            if vv == numA[1] - 1:
                print('true')
                thisAX.yaxis.tick_right()
                if not setLegend:
                    handles, legendLabels = thisAX.get_legend_handles_labels()
                    # Simplyfing to the Nt only
                    lLeg = []
                    # also the symbol only
                    hLeg = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
                    for ll, leg in enumerate(legendLabels):
                        lLeg.append(leg.split('_')[-1].split('x32')[0])
                    # This offset appropriate for smaller text
                    # thisAX.legend(hLeg, lLeg, bbox_to_anchor=(1.12, 1), loc='upper left', borderaxespad=0, handlelength=0)  # noqa: E501
                    # and this one for bigger
                    # thisAX.legend(hLeg, lLeg, bbox_to_anchor=(1.165, 1), loc='upper left', borderaxespad=0, handlelength=0)  # noqa: E501
                    # and this one for even bigger
                    thisAX.legend(hLeg, lLeg, bbox_to_anchor=(1.19, 1), loc='upper left', borderaxespad=0, handlelength=0)  # noqa: E501
                    setLegend = True

            # Increment analysis counter
            aa = aa + 1
    pdfName = os.path.join(anaDir, 'GPlot_' + params['cfuns']['norm'] + '.pdf')
    print(f'Saving pdf to {pdfName}')
    # plt.show()
    plt.savefig(pdfName)
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    mpl.rcParams['lines.markersize'] = 10.0
    # For Poster/Presentation
    mpl.rcParams['ytick.labelsize'] = 32
    mpl.rcParams['xtick.labelsize'] = 32
    mpl.rcParams['font.size'] = 32
    mpl.rcParams['legend.fontsize'] = 28
    main(sys.argv[1:])
