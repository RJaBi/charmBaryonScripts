"""
Just plots the inflection points from the csv
"""

import pandas as pd  # type: ignore
import sys
import os
import toml
import gvar as gv  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl
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

    # Load the csv
    inflectDFList = []
    for csv in params['inflectCSV']:
        print(f"Loading the csv from {csv}")
        inflectDF = pd.read_csv(csv)
        # Removing white space in column names
        inflectDF = inflectDF.rename(columns=lambda x: x.strip())
        print(inflectDF.head())
        inflectDFList.append(inflectDF)
    # Set up where to save
    anaDir = os.path.join(params['anaDir'])
    print('Analysis output directory is ', anaDir)
    if not os.path.exists(anaDir):
        os.makedirs(anaDir)
    Tpc = gv.gvar(params['Tpc'])
    fig, ax = plt.subplots(figsize=(16.6, 11.6))
    ax = GVP.myHSpan(Tpc, ax, alpha=0.8, lab=params['TpcLab'])
    for dd, df in enumerate(inflectDFList):
        print(df.columns)
        if params['TpcType'] == 'MeV':
            print(df['MeV'])
            print(df['MeVErr'])
            Tpc = gv.gvar(df['MeV'].values, df['MeVErr'].values)
            yScale = ' (MeV)'
        elif params['TpcType'] == 'TTpc':
            Tpc = gv.gvar(df['TTpc'].values, df['TTpcErr'].values)
            yScale = ''
        else:
            sys.exit(f'bad TpcType = {params["TpcType"]}')
        xAxis = df['ana'].values
        leftScript = params['leftSuperScript'][dd]
        for xx, xV in enumerate(xAxis):
            xAxis[xx] = '${}^{' + leftScript + '}' + xV[1:]
        markers = mo.markers
        for x, y, mark in zip(xAxis, Tpc, markers):
            ax = GVP.plot_gvEbar(x, y, ax, ma=mark, ls='')
    # finalise plot
    ax.set_ylabel('$T_c$' + yScale)
    ax.set_xlabel('Baryon')
    ax.legend(loc='best', ncol=2, fontsize=28)
    plt.savefig(os.path.join(anaDir, 'inflectionPoints.pdf'))
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    mpl.rcParams['lines.markersize'] = 10.0
    # For Poster/Presentation
    mpl.rcParams['ytick.labelsize'] = 32
    mpl.rcParams['xtick.labelsize'] = 32
    mpl.rcParams['font.size'] = 36
    main(sys.argv[1:])
