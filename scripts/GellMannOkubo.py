"""
Investigate the Gell-Mann-Okubo relations as in Eqnrs 16, 17, 18 of https://arxiv.org/pdf/0910.2419.pdf  # noqa: E501
16) J=1/2+
17) J=3/2+
18) Mix J=1/2+, J=3/2+
Simpler script which loads fits from massChoices.csv
"""

import numpy as np  # type: ignore
import os
import sys
import toml
import gvar as gv  # type: ignore
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.container import ErrorbarContainer  # type: ignore
import matplotlib as mpl  # type: ignore
from typing import List, Tuple
# Grabbing some modules from elsewhere
from pathlib import Path
insertPath = os.path.join(Path(__file__).resolve().parent, '..', 'lib')
sys.path.insert(0, insertPath)
import myModules as mo  # type: ignore # noqa: E402
import myGVPlotter as GVP  # type: ignore # noqa: E402


def eqn16(masses: np.ndarray) -> Tuple[gv.gvar, gv.gvar]:
    # The Gell-Mann-Okubo relation for J=1/2+ baryons
    # as defined in eqn16 of 09102419
    # Assumes masses are ordered left to right
    # 0.5 * (M_\xi + M_N) = (3.0 * M_\lambda + M_\Sigma)/4.0
    return 0.5 * (masses[0] + masses[1]), (3.0 * masses[2] + masses[3]) / 4.0


def eqn18(masses: np.ndarray) -> Tuple[gv.gvar, gv.gvar]:
    # The Gell-Mann-Okubo relation for J=1/2+/J=3/2+ baryons
    # as defined in eqn18 of 09102419
    # Assumes masses are ordered left to right
    # 3.0 * M_\lambda - M_\sigma - 2.0 * M_N = 2.0 * (M_\sigma^* - M_\Delta)
    return 3.0 * masses[0] - masses[1] - 2.0 * masses[2], 2.0 * (masses[3] - masses[4])


def getPhysMasses(params, hadrons: List, spins: List, OPs=False) -> np.ndarray:
    """
    Gets the physical masses of the hadrons with spin specified
    assumes positive parity
    Also gets the operator names optionally and returns that instead
    """
    # print(hadrons, spins)
    masses = np.empty(len(hadrons), dtype=object)
    if OPs:
        opNames = np.empty(len(hadrons), dtype=object)
    for anaName, val in params['eMass'].items():
        if anaName == 'unit':
            continue
        else:
            for hh, name in enumerate(hadrons):
                if val['name'] == name:
                    if val['J'] == spins[hh]:
                        if OPs:
                            opNames[hh] = anaName
                        else:
                            masses[hh] = gv.gvar(val['P'])
    # and return
    if OPs:
        return opNames
    else:
        return masses


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

    # Test the relations for physical (PDG) masses
    # Equation 16 JP=1/2+
    hadrons = ['\Xi_{cc}', 'N', '\Lambda_{c}', '\Sigma_{c}']  # noqa: W605
    spins = ['1/2'] * 4
    masses = getPhysMasses(params, hadrons, spins)
    LHS16, RHS16 = eqn16(masses)
    print('JP=1/2+', LHS16, RHS16, f'{(1 - (LHS16 / RHS16)) * 100} percent difference')

    # Equation 17 JP=3/2+
    # Not doing RHS cause Omega_{ccc} doesn't have physical meas
    hadrons = ['\Sigma_{c}', '\Delta', '\Xi_{cc}', '\Sigma_{c}']  # noqa: W605
    spins = ['3/2'] * len(hadrons)
    masses = getPhysMasses(params, hadrons, spins)
    LHS17 = masses[0] - masses[1]
    CHS17 = masses[2] - masses[3]
    print('JP=3/2+', LHS17, CHS17, f'{(1 - (LHS17 / CHS17)) * 100} percent difference')
    # Equation 18 - mixed JP=1/2+ and JP=3/2+
    hadrons = ['\Lambda_{c}', '\Sigma_{c}', 'N', '\Sigma_{c}', '\Delta']  # noqa: W605
    spins = ['1/2'] * 3 + ['3/2'] * 2
    masses = getPhysMasses(params, hadrons, spins)
    LHS18, RHS18 = eqn18(masses)
    print('JP=1/2+ = JP=3/2+', LHS18, RHS18, f'{(1 - (LHS18 / RHS18)) * 100} percent difference')
    # Now doing the lattice results
    aNames = ['lc16', 'lc17']  # 'lc18
    # Don't have the lattice results to do 18 and only 2/3 of 17
    temperatures = params['latMass']['temperatures']
    pdfName = os.path.join(anaDir, f'GellMannOkuboPlot.pdf')
    # now this is hwere savign to
    print(f'Saving pdf to {pdfName}')
    temperatures = params['latMass']['temperatures']
    # A plot
    fig, ax = plt.subplots(1, len(temperatures), figsize=(16.6, 11.6), sharey=True, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})  # noqa: E501
    # Iterate over different temperatures
    # cols = params['latMass']['colours']
    # marks = params['latMass']['markers']
    for tt, temp in enumerate(temperatures):
        NT = params['latMass']['NT'][tt]
        print('')
        print(NT, temp)
        # print(tt, temp, NT)
        thisAX = ax[tt]
        # and then over different hadrons
        for aa, aName in enumerate(params['latMass']['anaName']):
            if aName not in aNames:
                sys.exit(f'Bad {aName} not in {aNames}')
            if aName == 'lc16':
                # Do the eqn16 with light and charm quarks
                hadrons = ['\Xi_{cc}', 'N', '\Lambda_{c}', '\Sigma_{c}']  # noqa: W605
                spins = ['1/2'] * 4
                # Get the operator names
                operators = getPhysMasses(params, hadrons, spins, OPs=True)
                masses = []  # type: ignore
                # So can get the lattice data
                for op in operators:
                    aDF = massDF.query('Operator == @op')
                    aDF = aDF.rename(columns=lambda x: x.strip())
                    # get the hadron name
                    # had = aDF['Hadron'].values[0][1:-1]  # Strip the $ from it
                    thisNT = int(NT)  # noqa: F841
                    df = aDF.query('Nt == @thisNT')
                    EPSys = df['EPSys'].values[0]
                    masses.append(EPSys)  # type: ignore
                # Convert to physical
                masses = gv.gvar(masses) * hbarc / a_t  # type: ignore
                # Calculate the quantities
                thisNT_LHS16, thisNT_RHS16 = eqn16(masses)
                # legend labels
                LHSLab = '$('+f'{hadrons[0]} + {hadrons[1]}' + ')/2$'
                RHSLab = '$(' + f'3{hadrons[2]} + {hadrons[3]}' + ')/4$'
                # print('JP=1/2+', thisNT_LHS16, thisNT_RHS16, f'{(1 - (thisNT_LHS16 / thisNT_RHS16)) * 100} percent difference')  # noqa: E501
                # Plot LHS
                thisAX.plot(0.3, gv.mean(thisNT_LHS16), marker='d', linestyle='', label=LHSLab)
                col = thisAX.get_lines()[-1].get_color()
                thisAX = GVP.myFill_between([0.2, 0.4], [thisNT_LHS16] * 2, thisAX, ls='', colour=col, alpha=0.4)  # noqa: E501
                # Plot the experimental
                thisAX = GVP.myHSpan(LHS16, thisAX, alpha=0.8, colour=col)
                # and now RHS
                thisAX.plot(0.7, gv.mean(thisNT_RHS16), marker='d', linestyle='', label=RHSLab)
                col = thisAX.get_lines()[-1].get_color()
                thisAX = GVP.myFill_between([0.6, 0.8], [thisNT_RHS16] * 2, thisAX, ma='', ls='', colour=col, alpha=0.4)  # noqa: E501
                # Plot the experimental
                # and now RHS
                thisAX = GVP.myHSpan(RHS16, thisAX, alpha=0.8, colour=col)
                # The difference between RHS and LHS
                thisAX.plot(0.5, gv.mean(thisNT_RHS16 - thisNT_LHS16), marker='o', ls='', label=RHSLab + '\n$-'+LHSLab[1:])  # noqa: E501
                col = thisAX.get_lines()[-1].get_color()
                thisAX = GVP.myFill_between([0.4, 0.6], [thisNT_RHS16 - thisNT_LHS16] * 2, thisAX, ls='', colour=col, alpha=0.4)  # noqa: E501
                # Plot the experimental
                # col = thisAX.get_lines()[-1].get_c()
                thisAX = GVP.myHSpan(RHS16 - LHS16, thisAX, alpha=0.8, colour=col)
            elif aName == 'lc17':
                # Do the eqn17 with light and charm quarks
                hadrons = ['\Xi_{cc}', '\Sigma_{c}', '\Omega_{ccc}', '\Xi_{cc}']  # noqa: W605
                spins = ['3/2'] * len(hadrons)
                # Get the operator names
                operators = getPhysMasses(params, hadrons, spins, OPs=True)
                masses = []  # type: ignore
                # So can get the lattice data
                for op in operators:
                    aDF = massDF.query('Operator == @op')
                    aDF = aDF.rename(columns=lambda x: x.strip())
                    # get the hadron name
                    # had = aDF['Hadron'].values[0][1:-1]  # Strip the $ from it
                    thisNT = int(NT)  # noqa: F841
                    df = aDF.query('Nt == @thisNT')
                    EPSys = df['EPSys'].values[0]
                    masses.append(EPSys)  # type: ignore
                # Convert to physical
                masses = gv.gvar(masses) * hbarc / a_t  # type: ignore
                # Calculate the quantity
                thisNT_CHS17 = masses[0] - masses[1]
                thisNT_RHS17 = masses[2] - masses[3]
                # The label
                CHSLab = '$' + f'{hadrons[0]} - {hadrons[1]}' + '$'
                RHSLab = '$' + f'{hadrons[2]} - {hadrons[3]}' + '$'
                # print('JP=1/2+', thisNT_CHS17, thisNT_RHS17, f'{(1 - (thisNT_CHS17 / thisNT_RHS17)) * 100} percent difference')  # noqa: E501
                thisAX.plot(0.3, gv.mean(thisNT_CHS17), marker='p', ls='', label=CHSLab)
                col = thisAX.get_lines()[-1].get_color()
                thisAX = GVP.myFill_between([0.2, 0.4], [thisNT_CHS17] * 2, thisAX, colour=col, ls='', alpha=0.4)  # noqa: E501
                # Plot the experimental
                thisAX = GVP.myHSpan(CHS17, thisAX, alpha=0.8, colour=col)
                # and now the RHS
                # thisAX = GVP.plot_gvEbar(0.5, thisNT_RHS17, thisAX, ma='X', ls='', lab=RHSLab)
                thisAX.plot(0.7, gv.mean(thisNT_RHS17), marker='X', ls='', label=RHSLab)
                col = thisAX.get_lines()[-1].get_color()
                thisAX = GVP.myFill_between([0.6, 0.8], [thisNT_RHS17] * 2, thisAX, colour=col, ls='', alpha=0.4)  # noqa: E501
                # No experimental RHS
                # thisAX = GVP.myHSpan(RHS17, thisAX, alpha=0.8, colour=col)
            else:
                # Other ones not implemented yet
                continue

        # Get the labels from tt == 0
        # So that have only 1 label for each
        if tt == 0:
            handles, legendLabels = thisAX.get_legend_handles_labels()
            # also the symbol only
            hLeg = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
        if tt == len(temperatures) - 1:
            # Putting a dashed line for experiment in legend
            line_dashed = mpl.lines.Line2D([], [], color='black', linestyle='--', label='Experiment')  # noqa: E501
            hLeg.append(line_dashed)
            legendLabels.append('Experiment')
            thisAX.legend(hLeg, legendLabels, bbox_to_anchor=(1.2, 1), borderaxespad=0, handlelength=1.5)  # noqa: E501
        # Setting the xlim to make plot look nice
        thisAX.set_xlim([0, 1])
        thisAX.set_xlabel(f'{temp}')
        thisAX.set_xticklabels([])
        thisAX.get_xaxis().set_ticks([])
        fig.supxlabel('Temperature (MeV)')
        ax[0].set_ylabel('Mass (GeV)')
    # plt.show()
    if 'ylim' in params.keys():
        ax[0].set_ylim(params['ylim'])
    fig.savefig(pdfName)
    plt.close(fig)
    plt.show()
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    mpl.rcParams['lines.markersize'] = 8
    main(sys.argv[1:])
