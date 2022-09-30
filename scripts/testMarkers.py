"""
Just plots the inflection points from the csv
"""

import pandas as pd  # type: ignore
import sys
import os
import gvar as gv  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl
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
    """
    # J=1/2
    csvPath='/home/ryan/Documents/2022/Gen2L/Hawk/charmBaryonData/output/parityRatio/1_200'
    params = {
        'inflectCSV': ['N(uud)_RRatio_TauN4_M_Sdev.txt',
                       'Sigma(uus)_RRatio_TauN4_M_Sdev.txt',
                       'Xi(dss)_RRatio_TauN4_M_Sdev.txt',
                       'Sigma_c(udc)_RRatio_TauN4_M_Sdev.txt',
                       'Xi_cprime(usc)_RRatio_TauN4_M_Sdev.txt',
                       'Omega_c(ssc)_RRatio_TauN4_M_Sdev.txt',
                       'Lambda_c(udc)_RRatio_TauN4_M_Sdev.txt',
                       'Xi_c(usc)_RRatio_TauN4_M_Sdev.txt',
                       'Xi_cc(ccu)_RRatio_TauN4_M_Sdev.txt',
                       'Omega_cc(ccs)_RRatio_TauN4_M_Sdev.txt'
        ],
        'figText': ["$N(uud)$", "$\\Sigma(uus)$", "$\\Xi(dss)$", "$\\Sigma_c(udc)$", "$\\Xi_{c}^{\\prime}(usc)$", "$\\Omega_{c}(ssc)$", "$\\Lambda_{c}(udc)$", "$\\Xi_{c}(usc)$", "$\\Xi_{cc}(ccu)$", "$\\Omega_{cc}(ccs)$" ],  # noqa: E501
        'figTextYPos': [ 0.028, 0.08, 0.132, 0.184, 0.236, 0.288, 0.36, 0.42, 0.64, 0.70]
        }
    """
    # J=3/2
    csvPath = '/home/ryan/Documents/2022/Gen2L/Hawk/charmBaryonData/output/parityRatioJ3_2/1_200'
    params = {
        'inflectCSV': ['Sigma_c(udc)_RRatio_TauN4_M_Sdev.txt',
                       'Xi_c(usc)_RRatio_TauN4_M_Sdev.txt',
                       'Omega_c(ssc)_RRatio_TauN4_M_Sdev.txt',
                       'Xi_cc(ccu)_RRatio_TauN4_M_Sdev.txt',
                       'Omega_cc(ccs)_RRatio_TauN4_M_Sdev.txt',
                       'Omega_ccc(ccc)_RRatio_TauN4_M_Sdev.txt'
        ],
        'figText': ["$\\Sigma_{c}(udc)$", "$\\Xi_{c}(usc)$", "$\\Omega_{c}(ssc)$", "$\\Xi_{cc}(ccu)$", "$\\Omega_{cc}(ccs)$", "$\\Omega_{ccc}(ccc)$"],  # noqa: E501
        'figTextYPos': [0.13, 0.2, 0.29, 0.42, 0.490, 0.68]
        }

    marks = mo.markers
    # Load the csv
    inflectDFList = []

    for csv in params['inflectCSV']:  # type: ignore
        print(f"Loading the csv from {csv}")
        inflectDF = pd.read_csv(os.path.join(csvPath, csv))
        # Removing white space in column names
        inflectDF = inflectDF.rename(columns=lambda x: x.strip())
        inflectDFList.append(inflectDF)
    # Now plot
    fig, ax = plt.subplots(figsize=(16.6, 11.6))
    for dd, df in enumerate(inflectDFList):
        ax = GVP.plot_gvEbar(df['XTemp'], gv.gvar(df['mean'].values, df['err'].values), ax, ma=marks[dd], ls='--')  # noqa: E501
        xt = df['XTemp'].values[-1] + 0.2 * abs(df['XTemp'][0] - df['XTemp'][1])
        yt = params['figTextYPos'][dd]  # type: ignore
        text = ax.text(xt, yt, params['figText'][dd])  # type: ignore
        # Putting a marker in to the right of the text
        text.draw(fig.canvas.get_renderer())
        ex = text.get_window_extent()
        tex = ax.transData.inverted().transform_bbox(ex)
        xM = tex.x0 + 1050.0 * tex.height
        yM = tex.y0 + 0.5 * tex.height
        ax.plot(xM, yM, marker=marks[dd], linestyle='')

    ax.set_ylim([0, 1])
    ax.set_xlabel('Temperature (MeV)')
    ax.set_ylabel('$R$')
    pdf = PdfPages('markerTest.pdf')
    pdf.savefig(fig)
    plt.close(fig)
    pdf.close()


if __name__ == '__main__':
    mo.initBigPlotSettings()
    mpl.rcParams['lines.markersize'] = 16.0
    # For Poster/Presentation
    mpl.rcParams['ytick.labelsize'] = 32
    mpl.rcParams['xtick.labelsize'] = 32
    mpl.rcParams['font.size'] = 36
    main(sys.argv[1:])
