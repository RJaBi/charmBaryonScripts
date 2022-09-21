"""
Plot as function of temperature
Simpler script which loads fits from massChoices.csv
"""

import numpy as np  # type: ignore
import os
import sys
import toml
import gvar as gv  # type: ignore
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.backends.backend_pdf import PdfPages  # type: ignore
# Grabbing some modules from elsewhere
from pathlib import Path
insertPath = os.path.join(Path(__file__).resolve().parent, '..', 'lib')
sys.path.insert(0, insertPath)
import myModules as mo  # type: ignore # noqa: E402
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
    # Now iterate over the different operators
    for aa, aName in enumerate(params['latMass']['anaName']):
        print(aa, aName)
        aDF = massDF.query('Operator == @aName')
        aDF = aDF.rename(columns=lambda x: x.strip())
        # get the hadron name
        hadron = aDF['Hadron'].values[0][1:-1]  # Strip the $ from it
        print(hadron)
        # Just some things to hold data
        EP = []
        EPSys = []
        # EP_Sys = []
        EM = []
        EMSys = []
        # EM_Sys = []
        temps = params['latMass']['temperatures'][aa]
        # Loop over the Nt to get the ones we want
        for tt, NT in enumerate(params['latMass']['mALabels'][aa]):
            thisNT = int(NT)  # noqa: F841
            df = aDF.query('Nt == @thisNT')
            EP.append(df['EP'].values[0])
            EPSys.append(df['EPSys'].values[0])
            # EP_Sys.append(df['EP_Sys'].values[0])
            EM.append(df['EM'].values[0])
            EMSys.append(df['EMSys'].values[0])

        # Convert to physical
        EP = gv.gvar(np.asarray(EP)) * hbarc / a_t  # type: ignore
        EPSys = gv.gvar(np.asarray(EPSys)) * hbarc / a_t  # type: ignore
        EM = gv.gvar(np.asarray(EM)) * hbarc / a_t  # type: ignore
        EMSys = gv.gvar(np.asarray(EMSys)) * hbarc / a_t  # type: ignore
        # pdf
        if 'norm' in params.keys():
            if params['norm'] == 'EP0':
                pdfMod = 'EP0'
            else:
                sys.exit(f'bad normamlisation norm={params["norm"]} selected. Exiting')
            pdfName = os.path.join(anaDir, f'simpleEMPPlot_{params["latMass"]["anaName"][aa]}_norm{pdfMod}.pdf')
        else:
            pdfName = os.path.join(anaDir, f'simpleEMPPlot_{params["latMass"]["anaName"][aa]}.pdf')
        print(f'Saving pdf to {pdfName}')
        pdf = PdfPages(pdfName)
        # Getting the name and spin
        J = params['eMass'][aName]['J']
        # Getting the physical mass
        if params['eMass'][aName]['P'] != '':
            physEP = gv.gvar(params['eMass'][aName]['P'])
        else:
            physEP = gv.gvar(None)
        if params['eMass'][aName]['M'] != '':
            physEM = gv.gvar(params['eMass'][aName]['M'])
        else:
            physEM = gv.gvar(None)
        # have now got data
        # print(physEP, physEM, EP, EM, EPSys, EMSys)
        # open a figure
        fig, ax = plt.subplots(figsize=(16.6, 11.6))
        # maybe normmalise
        if 'norm' in params.keys():
            if params['norm'] == 'EP0':
                # Normallise by these values
                normVal = EP[0]
                normValSys = EPSys[0]
                ax.set_ylabel('$M_{{}^{' + J + '}' + hadron + '}$/$M_{{}^{' + J + '}' + hadron + '}(T='+f'{temps[0]}'+')$')  # noqa: W605, E501
            else:
                sys.exit(f'bad normamlisation norm={params["norm"]} selected. Exiting')
            # Do the normallisation
            EP = EP / normVal  # type: ignore
            EM = EM / normVal  # type: ignore
            EPSys = EPSys / normValSys  # type: ignore
            EMSys = EMSys / normValSys  # type: ignore
            physEP = physEP / normVal  # type: ignore
            physEM = physEM / normVal  # type: ignore
        else:
            ax.set_ylabel('$M_{{}^{' + J + '}' + hadron + '}$ (GeV)')  # noqa: W605, E501
        # Plot the data
        # Do stat and sys err on top of each other to show 2 error bars
        ax = GVP.plot_gvEbar(temps, EP, ax, ma='d', ls='', lab='$m_{+}$', col='tab:blue')
        ax = GVP.plot_gvEbar(temps, EPSys, ax, ma='d', ls='', col='tab:blue', alpha=0.5)
        ax = GVP.plot_gvEbar(temps, EM, ax, ma='d', ls='', lab='$m_{-}$', col='tab:red')
        ax = GVP.plot_gvEbar(temps, EMSys, ax, ma='d', ls='', col='tab:red', alpha=0.5)
        # Plot lines across from zero temp
        ax = GVP.myHSpan(EP[0], ax, colour='tab:blue', ls='--')
        ax = GVP.myHSpan(EM[0], ax, colour='tab:red', ls='--')
        # Sometimes don't have physical mass
        # lazy try/except
        try:
            ax = GVP.myFill_between([0, 0.1], [physEP, physEP], ax, alpha=0.8, colour='tab:blue')  # noqa: E501
        except:  # noqa: E722
            print('EP Phys failed')
        try:
            ax = GVP.myFill_between([0, 0.1], [physEM, physEM], ax, alpha=0.8, colour='tab:red')  # noqa: E501
        except:  # noqa: E722
            print('EP Phys failed')
        ax.set_xlabel('$T/T_{c}$')  # noqa: W605
        ax.set_xlim([0, None])
        if 'yLim' in params.keys():
            ax.set_ylim(params['yLim'])
        ax.legend(loc='lower right', ncol=2)
        pdf.savefig(fig)
        plt.close(fig)
        pdf.close()
    # out of loop
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    main(sys.argv[1:])
