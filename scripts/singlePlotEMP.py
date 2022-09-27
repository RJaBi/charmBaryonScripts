"""
Plot as function of temperature
each temperature is a separate subplot
Simpler script which loads fits from massChoices.csv
"""

import numpy as np  # type: ignore
import os
import sys
import toml
import gvar as gv  # type: ignore
import pandas as pd  # type: ignore
from typing import Dict, List, Any
import matplotlib.pyplot as plt  # type: ignore
# from matplotlib.container import ErrorbarContainer  # type: ignore
from matplotlib.lines import Line2D  # type: ignore
from matplotlib import rcParams  # type: ignore
# Grabbing some modules from elsewhere
from pathlib import Path
insertPath = os.path.join(Path(__file__).resolve().parent, '..', 'lib')
sys.path.insert(0, insertPath)
import myModules as mo  # type: ignore # noqa: E402
import myGVPlotter as GVP  # type: ignore # noqa: E402


def setYlim(params, ax):
    """
    Sets the ylims on each of the subplots
    such that they have the same scale
    """
    # First iterate over each axis
    # yScale is ymax - ymin
    yScale = 0
    # Middles is the middle point of each subplot
    middles = np.empty(np.shape(ax))
    yLims = np.empty(list(np.shape(ax)) + [2])
    # Determine the yscale
    for vv in range(0, np.shape(ax)[0]):
        for hh in range(0, np.shape(ax)[1]):
            # Need to consider the difference across all temperatures
            if hh == 0:
                ymin, ymax = ax[vv, hh].get_ylim()
            else:
                thisYMin, thisYMax = ax[vv, hh].get_ylim()
                if thisYMax > ymax:
                    ymax = thisYMax
                if thisYMin < ymin:
                    ymin = thisYMin
            hhvvYScale = abs(ymax - ymin)
            if hhvvYScale > yScale:
                yScale = hhvvYScale
            middles[vv, hh] = ymax - hhvvYScale / 2
            yLims[vv, hh, :] = ymin, ymax
            if 'ylim' in params.keys():
                ax[vv, hh].set_ylim(params['ylim'])
    if 'ylim' in params.keys():
        return ax
    # Now set it
    for vv in range(0, np.shape(ax)[0]):
        vertMid = np.median(middles[vv, :])
        yMinVV = vertMid - yScale / 2
        yMaxVV = vertMid + yScale / 2
        # Determine the middle point
        for hh in range(0, np.shape(ax)[1]):
            ymin = vertMid - yScale / 2
            ymax = vertMid + yScale / 2
            # Shift it up/down if necessary
            while ymax < yLims[vv, hh, 1]:
                ymax = ymax + 0.05 * yScale / 2
                ymin = ymin + 0.05 * yScale / 2
            while ymin > yLims[vv, hh, 0]:
                ymax = ymax - 0.05 * yScale / 2
                ymin = ymin - 0.05 * yScale / 2
            if ymin < yMinVV or ymax > yMaxVV:
                yMaxVV = ymax
                yMinVV = ymin
        # and finally set it
        for hh in range(0, np.shape(ax)[1]):
            ax[vv, hh].set_ylim([yMinVV, yMaxVV])
    # and return
    return ax


def main(args: list):
    """
    Reads some parameters in from the command line
    Puts that toml file in
    then reads it. Does a bunch of analysis
    """

    params = mo.GetArgs(args)
    # Printing the input toml back out - easier to read than printing as a dictionary
    toml.dump(params, f=sys.stdout)
    # refining ylims
    params = mo.refineXYLims(params, subDict=None)
    # For conversion to physical units
    hbarc = 0.1973269804  # GeV fm
    a_s = gv.gvar(params['latMass']['as'])
    xi = gv.gvar(params['latMass']['xi'])
    a_t = a_s/xi
    # Load the csv
    print(f"Loading the csv from {params['latMass']['massCSV']}")
    massDF = pd.read_csv(params['latMass']['massCSV'])
    print(massDF.head())
    # Set up where to save
    anaDir = os.path.join(params['latMass']['anaDir'])
    print('Analysis output directory is ', anaDir)
    if not os.path.exists(anaDir):
        os.makedirs(anaDir)

    if 'norm' in params.keys():
        sys.exit('Norm not implemented yet')
        if params['norm'] == 'EP0':
            pdfMod = 'EP0'
        else:
            sys.exit(f'bad normamlisation norm={params["norm"]} selected. Exiting')
        pdfName = os.path.join(anaDir, f'singleEMPPlot_norm_{pdfMod}.pdf')
        pdfNameM = os.path.join(anaDir, f'singleEMPPlot_NegParity_norm_{pdfMod}.pdf')
    else:
        pdfName = os.path.join(anaDir, f'singleEMPPlot.pdf')
        pdfNameM = os.path.join(anaDir, f'singleEMPPlot_NegParity.pdf')
    # now this is hwere savign to
    print(f'Saving pdf to {pdfName} and {pdfNameM}')
    temperatures = params['latMass']['temperatures']
    # A plot for negative and positive parities
    fig, ax = plt.subplots(params['latMass']['vertSplit'], len(temperatures), figsize=(16.6, 11.6), sharey=False, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})  # noqa: E501
    figM, axM = plt.subplots(params['latMass']['vertSplit'], len(temperatures), figsize=(16.6, 11.6), sharey=False, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})  # noqa: E501
    # Iterate over different temperatures
    cols = params['latMass']['colours']
    marks = params['latMass']['markers']
    allHandlesDict: Dict[int, List[Any]] = {}
    allLegendsDict: Dict[int, List[Any]] = {}
    # Setting up the x-axis
    vertDict = {}
    vertCountDict = {}
    uniVert = np.unique(params['latMass']['vertGroup'], return_counts=True)
    xRan = [0, 1]
    for ii in range(0, params['latMass']['vertSplit']):
        vertDict.update({uniVert[0][ii]: uniVert[1][ii]})
        allHandlesDict.update({uniVert[0][ii]: []})
        allLegendsDict.update({uniVert[0][ii]: []})
    for tt, temp in enumerate(temperatures):
        for ii in range(0, params['latMass']['vertSplit']):
            vertCountDict.update({uniVert[0][ii]: 1})
        NT = params['latMass']['NT'][tt]
        # and then over different hadrons
        for aa, aName in enumerate(params['latMass']['anaName']):
            thisVGroup = params['latMass']['vertGroup'][aa]
            thisAX = ax[thisVGroup, tt]
            thisAXM = axM[thisVGroup, tt]
            # Doing the physical mass
            # as lines
            if params['eMass'][aName]['P'] != '':
                physEP = gv.gvar(params['eMass'][aName]['P'])
            else:
                physEP = gv.gvar(None)
            if params['eMass'][aName]['M'] != '':
                physEM = gv.gvar(params['eMass'][aName]['M'])
            else:
                physEM = gv.gvar(None)
            thisAX = GVP.myHSpan(physEP, thisAX, colour=cols[aa], alpha=0.8)
            thisAXM = GVP.myHSpan(physEM, thisAXM, colour=cols[aa], alpha=0.8)
            if NT not in params['latMass']['mALabels'][aa]:
                continue
            # print(aa, aName)
            aDF = massDF.query('Operator == @aName')
            aDF = aDF.rename(columns=lambda x: x.strip())
            # get the hadron name
            had = aDF['Hadron'].values[0][1:-1]  # Strip the $ from it
            thisNT = int(NT)  # noqa: F841
            df = aDF.query('Nt == @thisNT')
            EP = df['EP'].values[0]
            EPSys = df['EPSys'].values[0]
            EM = df['EM'].values[0]
            EMSys = df['EMSys'].values[0]
            # Convert to physical
            EP = gv.gvar(np.asarray(EP)) * hbarc / a_t  # type: ignore
            EPSys = gv.gvar(np.asarray(EPSys)) * hbarc / a_t  # type: ignore
            EM = gv.gvar(np.asarray(EM)) * hbarc / a_t  # type: ignore
            EMSys = gv.gvar(np.asarray(EMSys)) * hbarc / a_t  # type: ignore
            if tt == 0:
                lab = '$' + had + '$'
            else:
                lab = None
            # Put on the plot at some xcoords
            # Width of the box
            width = abs(xRan[0] - xRan[1]) / (2 * vertDict[thisVGroup] + 1)
            xStart = xRan[0] + width * vertCountDict[thisVGroup]
            xStart = xStart - 0.15 * width
            xEnd = xRan[0] + width + width * vertCountDict[thisVGroup]
            xEnd = xEnd + 0.15 * width
            xMid = xStart + abs(xEnd - xStart) / 2
            vertCountDict[thisVGroup] = vertCountDict[thisVGroup] + 2
            # Fill in the legend
            if tt == 0:
                thisPointLabel = Line2D([], [], color=cols[aa], marker=marks[aa], linestyle='', label=lab, markersize=16)  # noqa: E501
                allHandlesDict[thisVGroup].append(thisPointLabel)
                allLegendsDict[thisVGroup].append(lab)
            # and now plot
            thisAX.plot(xMid, gv.mean(EP), marker=marks[aa], linestyle='', color=cols[aa], markeredgecolor='black')  # noqa: E501
            thisAX = GVP.myFill_between([xStart, xEnd], [EP] * 2, thisAX, ls='', colour=cols[aa], alpha=0.5)  # noqa: E501
            thisAX = GVP.myFill_between([xStart, xEnd], [EPSys] * 2, thisAX, ls='', colour=cols[aa], alpha=0.3)  # noqa: E501
            thisAXM.plot(xMid, gv.mean(EM), marker=marks[aa], linestyle='', color=cols[aa], markeredgecolor='black')  # noqa: E501
            thisAXM = GVP.myFill_between([xStart, xEnd], [EM] * 2, thisAXM, ls='', colour=cols[aa], alpha=0.5)  # noqa: E501
            thisAXM = GVP.myFill_between([xStart, xEnd], [EMSys] * 2, thisAXM, ls='', colour=cols[aa], alpha=0.3)  # noqa: E501
            # change y-ticks, etc
            if tt != 0:
                thisAX.get_yaxis().set_visible(False)
                thisAXM.get_yaxis().set_visible(False)
            if tt == len(temperatures) - 1:
                thisAX.yaxis.tick_right()
                thisAXM.yaxis.tick_right()
                thisAX.get_yaxis().set_visible(True)
                thisAXM.get_yaxis().set_visible(True)
        # outside loop
        thisAX.set_xlabel(f'{temp}')
        thisAX.get_xaxis().set_ticks([])
        thisAX.set_xticklabels([])
        fig.supxlabel('Temperature (MeV)')
        fig.supylabel('Mass (GeV)')
        ax[-1, tt].set_xlabel(f'{temp}')
        axM[-1, tt].set_xlabel(f'{temp}')
        for vs in range(0, params['latMass']['vertSplit']):
            charmness = params['latMass']['charmness'][vs]
            ax[vs, 0].text(0.15, 0.8, '$' + f'C={charmness}' + '$', transform=ax[vs, 0].transAxes, usetex=True)  # noqa: E501
            axM[vs, 0].text(0.15, 0.85, '$' + f'C={charmness}' + '$', transform=axM[vs, 0].transAxes, usetex=True)  # noqa: E501
        # Removing xticks
        thisAXM.set_xticklabels([])
        thisAXM.get_xaxis().set_ticks([])
        figM.supxlabel('Temperature (MeV)')
        figM.supylabel('Mass (GeV)')

    # Doing the legend
    # Add a line for experiment
    for k, v in allHandlesDict.items():
        allHandles = v
        allLegends = allLegendsDict[k]
        line_dashed = Line2D([], [], color='black', linestyle='--', label='Exp.')  # noqa: E501
        allHandles.append(line_dashed)
        allLegends.append('Exp.')
        ax[k, len(temperatures) - 1].legend(allHandles, allLegends, bbox_to_anchor=(2.2, 1), borderaxespad=0, handlelength=1.5)  # noqa: E501
        axM[k, len(temperatures) - 1].legend(allHandles, allLegends, bbox_to_anchor=(2.4, 1), borderaxespad=0, handlelength=1.5)  # noqa: E501
    # Now determine  and sety-limits
    ax = setYlim(params, ax)
    axM = setYlim(params, axM)
    ax[0, 0].set_xlim(xRan)
    fig.savefig(pdfName)
    plt.close(fig)
    figM.savefig(pdfNameM)
    plt.close(figM)
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    rcParams['lines.markersize'] = 8
    main(sys.argv[1:])
