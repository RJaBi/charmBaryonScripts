
import numpy as np  # type: ignore
import os
import sys
import toml
# from inspect import signature
import matplotlib.pyplot as plt  # type: ignore
import matplotlib as mpl  # type: ignore
from matplotlib.backends.backend_pdf import PdfPages  # type: ignore
from typing import Dict, List  # , Any, List, Callable
# from timeit import default_timer as timer
from scipy.interpolate import InterpolatedUnivariateSpline  # type: ignore
import gvar as gv  # type: ignore
import lsqfit as lsq  # type: ignore
# Loading some functions from elsewhere in this repo
from pathlib import Path
modDir = os.path.join(Path(__file__).resolve().parent, '..', 'lib')
sys.path.insert(0, modDir)
import myModules as mo  # type: ignore  # noqa: E402
import myIO as io  # type: ignore  # noqa: E402
import myGVPlotter as GVP  # type: ignore  # noqa: E402


def gvarCorrSum(G: np.ndarray) -> np.ndarray:
    """
    Sums the correlator, i.e.
    GS(t_n) = sum_{t_n}^{Nt} G(t_n)
    but here it's one dim
    """
    GS = np.empty(G.shape[0], dtype=object)
    for it in range(0, G.shape[0]):
        GS[it] = np.sum(G[it:int(G.shape[0]/2)-1])
    return np.asarray(GS)


def RRatio(GJ2: np.ndarray) -> np.ndarray:
    """
    Constructs the sum
    R = (sum_{n} R(t_n)/sigma^2(t_n) )/
        ( sum_{n} sigma^{-2}(t_n))
    """
    sigInv = (1.0 / gv.sdev(GJ2) ** 2.0)
    num = GJ2 * sigInv
    RR = gvarCorrSum(num) / gvarCorrSum(sigInv)
    RR = np.asarray(RR)
    # check for nan as well
    for ii, rii in enumerate(RR):
        if rii != rii:
            RR[ii] = gv.gvar(0, 0)
    return RR


# As suggested in 2007.04188
def arctan(x, p):
    return p['c0'] + p['c1'] * np.arctan(p['c2'] * (x - p['Tpc']))


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def gvIsClose(gv1, gv2, per=0.20):
    """
    returns True if gv1 is within per percent of gv2
    Both ways
    """
    gv1Max = gv.mean(gv1 + 0.2 * gv1)
    gv1Min = gv.mean(gv1 - 0.2 * gv1)
    gv2Max = gv.mean(gv2 + 0.2 * gv2)
    gv2Min = gv.mean(gv2 - 0.2 * gv2)
    close = True
    if gv1 > gv2Max:
        close = False
    elif gv1 < gv2Min:
        close = False
    elif gv2 > gv1Max:
        close = False
    elif gv2 < gv1Min:
        close = False
    return close


def main(args: list):
    """
    Reads some parameters in from the command line
    Puts that toml file in
    then reads it. Does a bunch of analysis
    """
    seed = gv.ranseed(9)  # noqa: F841
    params = mo.GetArgs(args)
    # Printing the input toml back out
    toml.dump(params, f=sys.stdout)
    # Setting x, limits to None if they were 'None' in the toml
    params = mo.refineXYLims(params)

    # Loading correlators
    cfDict, cfLabelsList = io.initCorrelators(params)

    # Now that we've got everything we want.
    # Time to do analysis.

    # the t_n to sum after
    RCurves = params['analysis']['RCurves']
    # The fit window we are considering
    # Here for folder location only
    fMin = params['analysis']['fitMin']
    fMax = params['analysis']['fitMax']

    # Check what scale for Tpc we are in
    # if not present set to 1
    # i.e. are we in T or T/Tpc
    if 'TpcScale' in params['analysis'].keys():
        TpcScale = params['analysis']['TpcScale']
    else:
        TpcScale = 1.0
    if TpcScale == 1.0:
        xScale = 'T/T_c'
        xLeg = ''
    else:
        xLeg = '(MeV)'
        xScale = 'T ' + xLeg
    # where we put the analysis
    anaDir = os.path.join(params['analysis']['anaDir'], f'{int(fMin)}_{int(fMax)}')
    print('Analysis output directory is ', anaDir)
    if not os.path.exists(anaDir):
        os.makedirs(anaDir)

    pdf = PdfPages(os.path.join(anaDir, 'R_gvar_ratios.pdf'))
    # as a function of all t_n
    figTau, axTau = plt.subplots(figsize=(16.6, 11.6))

    RRatios: List[Dict[str, np.ndarray]] = []
    for thisAna in params['analysis']['labelsToAnalyse']:
        thisRR: Dict[str, np.ndarray] = {}
        for ana in thisAna:
            if ana not in cfLabelsList:
                print('label ' + ana + ' is not in cfLabelsList:', cfLabelsList)
                sys.exit('Exiting')
            # Now do the analysis
            G = io.labelNorm(cfDict['data'][ana], params)
            # Taking mean and 2nd order jackknifes
            # GVDO = gv.dataset.avg_data(G)
            # GMean = gv.mean(GVD)
            # Also do the plot for each analysis separately
            thisAnaDir = os.path.join(anaDir, ana)
            if not os.path.exists(thisAnaDir):
                os.makedirs(thisAnaDir)
            # Do a svdcut on the data
            # Removes very small eigenvalues in the cov matrix
            s = gv.dataset.svd_diagnosis(G, nbstrap=1000)
            GVD = gv.svd(s.avgdata, svdcut=s.svdcut)
            GMean = gv.mean(GVD)
            # This does two things
            # Makes sure the svdcut hasn't broken the cov matrix
            # Removes knowledge that it's had an svdcut
            # This makes the RRatio work again
            GVD = gv.gvar(gv.mean(GVD), gv.evalcov(GVD), verify=True)
            GMean = gv.mean(GVD)
            x = np.asarray(range(1, len(GMean)+1))
            writeSteps = False
            if writeSteps:
                saveFile = os.path.join(thisAnaDir, 'mcRRatio_M_Sdev.txt')
                with open(saveFile, 'w') as f:
                    xcount = 0
                    f.write('tau, mean, err \n')
                    for mm, ss in zip(gv.mean(GVD), gv.sdev(GVD)):
                        f.write(f'{x[xcount]}, {mm}, {ss} \n')
                        xcount = xcount + 1
# !#            gvDump = os.path.join(thisAnaDir, 'GVD.gv')
# !#            if os.path.exists(gvDump):
# !#                print(f'LOADING FROM GVDUMP = {gvDump}')
# !#                #print(GVD.shape)
# !#                #GVD = gv.load(gvDump)
# !#                #print(GVD.shape)
# !#                #GMean = gv.mean(GVD)
# !#            else:
# !#                print(f'SAVING TO GVDUMP = {gvDump}')
# !#                gv.dump(GVD, gvDump, add_dependencies=True)
# !#
            # Where to save the pdf
            thisPDF = PdfPages(os.path.join(thisAnaDir, 'R_gvar_ratios.pdf'))
            RR = RRatio(GVD)
            thisRR.update({ana: RR})
            # Plotting the sum to check for excited states
            fig, ax = plt.subplots(figsize=(16.6, 11.6))
            if writeSteps:
                saveFile = os.path.join(thisAnaDir, 'RRatio_M_Sdev.txt')
                with open(saveFile, 'w') as f:
                    xcount = 0
                    f.write('tau_n, mean, err \n')
                    for mm, ss in zip(gv.mean(RR), gv.sdev(RR)):
                        f.write(f'{x[xcount]}, {mm}, {ss} \n')
                        xcount = xcount + 1
            ax = GVP.plot_gvEbar(x, RR, ax, lab=ana)
            ax.set_xlim(params['analysis']['RxLim'])
            ax.set_ylim(params['analysis']['RyLim'])
            ax.set_xlabel('$\\tau_n$')
            ax.set_ylabel('$R(\\tau_n)$')  # noqa: W605
            for tn in RCurves:
                itn = int(tn)
                ax.axvline(x=x[itn], linestyle='--', color='black', alpha=0.25, label=f'$\\tau_n = {itn}$')  # noqa: E501
            ax.legend(loc='best', ncol=2)
            thisPDF.savefig(fig)
            plt.close(fig)
            thisPDF.close()
            # Also plotting them all against each other
            axTau = GVP.plot_gvEbar(x, RR, axTau, lab=ana)
            # and saving that plot
        # update the list
        RRatios.append(thisRR)
    axTau.set_xlim(params['analysis']['RxLim'])
    axTau.set_ylim(params['analysis']['RyLim'])
    axTau.set_xlabel('$\\tau_n$')
    axTau.set_ylabel('$R(\\tau_n)$')  # noqa: W605
    axTau.legend(loc='best', ncol=2)
    pdf.savefig(figTau)
    plt.close(figTau)
    pdf.close()

    # Now plot the curves
    # the curves for single or few t_n
    figCurve, axCurve = plt.subplots(figsize=(16.6, 11.6))
    pdf = PdfPages(os.path.join(anaDir, 'R_gvar_curves.pdf'))
    cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    marks = mo.markers
    markCount = 0
    colCount = 0
    inflectVals = {}
    for ii, thisAna in enumerate(params['analysis']['labelsToAnalyse']):
        if len(RCurves) > len(cols):
            sys.exit(f"Number of curves: {len(RCurves)} is greater than number of colours {len(cols)}")  # noqa: E501
        print(thisAna, f'{params["analysis"]["figText"][ii]}')
        # Data holder for spline
        Rtn = np.empty([len(RCurves), len(RRatios[ii])], dtype=object)
        thisTnCount = 0
        # the x=axis for the points
        # Probably in T/T_c
        XTemps = np.float64(params['analysis']['RXPoints'])
        # Iterating over the tn, and the analysis
        for tn in RCurves:
            itn = int(tn)
            anaCount = 0
            for ana, RR in RRatios[ii].items():
                CV = gv.mean(RR[itn])
                jackErr = gv.sdev(RR[itn])
                # Packaging data for spline
                Rtn[thisTnCount, anaCount] = RR[itn]
                axCurve.errorbar(XTemps[anaCount], y=CV, yerr=jackErr, linestyle='', marker=marks[markCount], color=cols[colCount])  # noqa: E501

                # Increment the counter for labelling
                anaCount = anaCount + 1

            if 'figTextYPos' in params['analysis'].keys():
                # Plotting some text at end of lines
                aC = anaCount - 1
                xt = XTemps[aC] + 0.2 * abs(XTemps[0] - XTemps[1])
                xtP = XTemps[aC] + 0.7 * abs(XTemps[0] - XTemps[1])
                yt = params['analysis']['figTextYPos'][ii]
                axCurve.text(xt, yt, params['analysis']['figText'][ii])
                axCurve.plot(xtP, yt, marker='d', alpha=0)

            if writeSteps:
                saveFile = os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_RRatio_TauN{itn}_M_Sdev.txt')  # noqa: E501
                saveFile = mo.clean_filename(saveFile, whitelist=['/'])
                with open(saveFile, 'w') as f:
                    f.write('XTemp, mean, err \n')
                    for xx in range(0, len(XTemps)):  # type: ignore
                        f.write(f'{XTemps[xx]}, {gv.mean(Rtn[thisTnCount, xx])},  {gv.sdev(Rtn[thisTnCount, xx])} \n')  # noqa: E501
            # Now fit the spline
            tS = np.linspace(XTemps.min(), XTemps.max(), 5000)
            tnSpline = gv.cspline.CSpline(np.float64(XTemps), Rtn[thisTnCount, :], alg='cspline', extrap_order=3)  # noqa: E501

            # The second derivative gives the inflection point
            def derivFunc(x):
                return tnSpline.D2(x) - 0

            # use scipy as a crosscheck
            ncon = 2000
            inflectEval = np.empty(ncon)
            bootIter = gv.bootstrap_iter(Rtn[thisTnCount, :])
            splineCoeffs = np.empty([len(XTemps), ncon])  # type: ignore
            deriv2Coeffs = np.empty([len(XTemps) - 2, ncon])  # type: ignore
            for icon in range(0, ncon):
                RtnIcon = gv.mean(next(bootIter))
                spSpline = InterpolatedUnivariateSpline(np.float64(XTemps), RtnIcon, k=3)
                # Also keep the coeffs
                splineCoeffs[:, icon] = spSpline.get_coeffs()
                deriv2 = spSpline.derivative(2)
                # Also keep the coeffs of the spline
                deriv2Coeffs[:, icon] = deriv2.get_coeffs()
                # Now evaluate it
                deriv2Eval = deriv2.__call__(tS)
                zeroArgMin = np.argmax(abs(deriv2Eval-0))
                zeroClosestVal = tS[zeroArgMin]
                # Putting in a shape 2 thing for the jackknife error code
                inflectEval[icon] = zeroClosestVal
            # Make it into a gvar
            scipyInflectVal = gv.dataset.avg_data(inflectEval)
            splineKnots = spSpline.get_knots()
            deriv2Knots = deriv2.get_knots()
            if writeSteps:
                # Save the inflection point from scipy
                inflectFile = mo.clean_filename(os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_splineInflectVal_TauN{itn}_boot.txt'), whitelist=['/'])  # noqa: E501
                splineKnotsFile = mo.clean_filename(os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_splineKnots_TauN{itn}_boot.txt'), whitelist=['/'])  # noqa: E501
                splineCoeffsFile = mo.clean_filename(os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_splineCoeffs_TauN{itn}_boot.txt'), whitelist=['/'])  # noqa: E501
                deriv2KnotsFile = mo.clean_filename(os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_deriv2Knots_TauN{itn}_boot.txt'), whitelist=['/'])  # noqa: E501
                deriv2CoeffsFile = mo.clean_filename(os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_deriv2Coeffs_TauN{itn}_boot.txt'), whitelist=['/'])  # noqa: E501
                inflectF = open(inflectFile, 'w')
                sKnotsF = open(splineKnotsFile, 'w')
                sCoeffsF = open(splineCoeffsFile, 'w')
                dKnotsF = open(deriv2KnotsFile, 'w')
                dCoeffsF = open(deriv2CoeffsFile, 'w')
                inflectF.write('icon, inflectVal \n')
                sCoeffsF.write('icon, coeffs \n')
                sKnotsF.write('knots \n')
                dCoeffsF.write('icon, coeffs \n')
                dKnotsF.write('knots \n')
                for icon, val in enumerate(inflectEval):
                    inflectF.write(f'{icon}, {val} \n')
                    sCoeffsF.write(f'{icon}, {splineCoeffs[:, icon]} \n')
                    dCoeffsF.write(f'{icon}, {deriv2Coeffs[:, icon]} \n')
                sKnotsF.write(f'{splineKnots} \n')
                dKnotsF.write(f'{deriv2Knots} \n')
                inflectF.close()
                sKnotsF.close()
                sCoeffsF.close()
                dKnotsF.close()
                dCoeffsF.close()
            print(f'scipy inflection point = {scipyInflectVal}')
            try:
                # Use gvar splines
                minSearch = 0.771 * TpcScale
                interval = gv.root.search(derivFunc, minSearch)
                inflectVal = gv.root.refine(derivFunc, interval)
                print(f'minSearch = {minSearch}, inflection Point = {inflectVal}')
                # # for two different interval
                minSearch = 0.925 * TpcScale
                interval = gv.root.search(derivFunc, minSearch)
                inflectVal = gv.root.refine(derivFunc, interval)
                print(f'minSearch = {minSearch}, inflection Point = {inflectVal}')
            except ValueError:
                print('USING SCIPYINFLECTVAL')
                inflectVal = scipyInflectVal

            # Try an arctan fit
            # Priors chosen to give approximately corret shape
            prior = {}
            if TpcScale == 1.0:
                prior['c0'] = gv.gvar(1.0, 0.7)
                prior['c1'] = gv.gvar(-0.5, 1.0)
                prior['c2'] = gv.gvar(10, 10)
                prior['Tpc'] = gv.gvar(1, 0.1)
            else:
                prior['c0'] = gv.gvar(1.0, 0.7)
                prior['c1'] = gv.gvar(-0.5, 1.0)
                prior['c2'] = gv.gvar(0.10, 0.10)
                prior['Tpc'] = gv.gvar(1 * TpcScale, 0.1 * TpcScale)
            try:
                fit = lsq.nonlinear_fit(data=(np.float64(XTemps), Rtn[thisTnCount, :]), fcn=arctan, prior=prior)  # noqa: E501
                print(fit.format(maxline=True))
                fitInflectVal = fit.p['Tpc']
                # Fit only the points of the spline around where expect Tpc to be
                nearX = find_nearest(np.float64(XTemps), TpcScale)
                print(f'nearX = {nearX}')
                if nearX - 2 < 0:
                    xStart = 0
                    xEndMod = abs(nearX-2)
                else:
                    xStart = nearX - 2
                    xEndMod = 0
                if nearX + 3 + xEndMod > len(XTemps):  # type: ignore
                    xEnd = len(XTemps)  # type: ignore
                else:
                    xEnd = nearX + 3 + xEndMod
                fitXRange = np.float64(XTemps)[xStart: xEnd]
                fitYRange = Rtn[thisTnCount, xStart: xEnd]
                fit = lsq.nonlinear_fit(data=(fitXRange, fitYRange), fcn=arctan, prior=prior)  # noqa: E501
                print(fit.format(maxline=True))
                fitInflectVal = fit.p['Tpc']
            except ValueError:
                print('Value Error in fit. Probably residuals not finite')
                print('Using scipyInflectVal for fitInflectVal')
                fitInflectVal = scipyInflectVal

            # Get final number for inflection point
            # Only use the GV value if it's close (34%) to other methods
            close = []
            # for lab, val in zip(['GV', 'SCIPY', 'FIT'], [inflectVal, scipyInflectVal, fitInflectVal]):  # noqa: E501
            for lab, val in zip(['GV'], [inflectVal, scipyInflectVal, fitInflectVal]):  # noqa: E501
                print(lab, val)
                for lab2, val2 in zip(['GV', 'SCIPY', 'FIT'], [inflectVal, scipyInflectVal, fitInflectVal]):  # noqa: E501
                    if (lab == lab2):
                        continue
                    else:
                        print('    ', lab2, val2, gvIsClose(val, val2, per=0.34))
                        thisClose = gvIsClose(val, val2, per=0.34)
                        close.append(thisClose)
            if np.any(np.asarray(close)):
                finalInflectVal = inflectVal
            else:
                finalInflectVal = scipyInflectVal
                for lab, val in zip(['GV', 'SCIPY', 'FIT'], [inflectVal, scipyInflectVal, fitInflectVal]):  # noqa: E501
                    print(lab, val)
            print(f'Using {finalInflectVal}')
            # add inflection point to dictionary for later
            inflectVals.update({params["analysis"]["figText"][ii]: finalInflectVal})
            # Plotting more
            CV = gv.mean(finalInflectVal)
            if CV > 0.8 * TpcScale and CV < 1.2 * TpcScale:
                if xLeg != '':
                    label = f'{params["analysis"]["figText"][ii]}, $\\tau_n = {itn}$' + ', $T_{inf} = '+f'{finalInflectVal}' + '$ $' + xLeg + '$'  # noqa: E501
                else:
                    label = f'{params["analysis"]["figText"][ii]}, $\\tau_n = {itn}$' + ', $T_{inf} = '+f'{finalInflectVal}' + '$'  # noqa: E501
                print(f'LABEL IS {label}')
                axCurve = GVP.myVSpan(finalInflectVal, axCurve, colour=cols[colCount], lab=label)
                # Here don't need a label for spline
                axCurve = GVP.myFill_between(tS, tnSpline(tS), axCurve, colour=cols[colCount], ls='--')  # noqa: E501
                if writeSteps:
                    saveFile = os.path.join(anaDir, f'{params["analysis"]["figText"][ii]}_splineEval_TauN{itn}_M_Sdev.txt')  # noqa: E501
                    saveFile = mo.clean_filename(saveFile, whitelist=['/'])
                    tSave = np.linspace(XTemps.min(), XTemps.max(), 1000)
                    with open(saveFile, 'w') as f:
                        f.write('temp, mean, err \n')
                        for xx in range(0, len(tSave)):
                            f.write(f'{tSave[xx]}, {gv.mean(tnSpline(tSave[xx]))},  {gv.sdev(tnSpline(tSave[xx]))} \n')  # noqa: E501
            else:
                # Here do need one for spline
                splineLabel = f'{params["analysis"]["figText"][ii]}, $\\tau_n = {itn}$'
                axCurve = GVP.myFill_between(tS, tnSpline(tS), axCurve, colour=cols[colCount], lab=splineLabel, alpha=0.25, ls='--')  # noqa: E501
            if markCount == len(marks) - 1:
                markCount = 0
            else:
                markCount = markCount + 1
            if colCount == len(cols) - 1:
                colCount = 0
            else:
                colCount = colCount + 1
            thisTnCount = thisTnCount + 1
    axCurve.set_ylim(params['analysis']['RyLim'])
    axCurve.set_xlabel('$' + xScale + '$')
    axCurve.set_ylabel('$R$')  # noqa: W605

    # axCurve.legend(loc='best', ncol=2)
    pdf.savefig(figCurve)
    plt.close(figCurve)
    pdf.close()

    # Do comparision plot to Tpc from other method
    if 'Tpc' in params['analysis'].keys() and 'TpcLabels' in params['analysis'].keys():
        fig, ax = plt.subplots(figsize=(16.6, 11.6))
        Tpc = gv.gvar(params['analysis']['Tpc'])
        TpcLabs = params['analysis']['TpcLabels']
        if len(Tpc) != len(TpcLabs):
            sys.exit(f'Tpc = {Tpc} and TpcLabs = {TpcLabs} have different lengths')
        # open pdf
        pdf = PdfPages(os.path.join(anaDir, 'R_gvar_Tpc.pdf'))
        # First plot the value from other sources
        for tt, temp in enumerate(Tpc):
            ax = GVP.myHSpan(temp, ax, lab=TpcLabs[tt])
        count = 1
        markCount = 0
        colCount = 0
        anaList = []
        for ana, temp in inflectVals.items():
            if temp < 0.8 * TpcScale or temp > 1.2 * TpcScale:
                continue
            if TpcScale == 1.0:
                tempMeV = temp * Tpc[0]
            else:
                tempMeV = temp
            print(ana, tempMeV)
            anaList.append(ana)
            ax = GVP.plot_gvEbar(ana, tempMeV, ax, ma=marks[markCount], col=cols[colCount], ls='--')
            # and iterate mark/colour counters
            count = count + 1
            if markCount == len(marks) - 1:
                markCount = 0
            else:
                markCount = markCount + 1
            if colCount == len(cols) - 1:
                colCount = 0
            else:
                colCount = colCount + 1
        # change the x-axis
        ax.set_xticklabels(anaList, rotation=45, ha='right')
        # ax.set_xlim([0, count + 1])
        ax.set_ylabel('$T_c$ MeV')
        ax.legend(loc='best', ncol=2)
        pdf.savefig(fig)
        plt.close(fig)
        pdf.close()
        # Write the inflection points to a file
        if writeSteps:
            saveFile = os.path.join(anaDir, 'Inflect_M_Sdev.txt')
            with open(saveFile, 'w') as f:
                f.write(f'ana, MeV, MeVErr, TTpc, TTpcErr \n')
                for ana, temp in inflectVals.items():
                    tempMeV = temp * Tpc[0]
                    f.write(f'{ana}, {gv.mean(tempMeV)}, {gv.sdev(tempMeV)}, {gv.mean(temp)}, {gv.sdev(temp)} \n')  # noqa: E501
    sys.exit('Finished')


if __name__ == '__main__':
    mo.initBigPlotSettings()
    mpl.rcParams['lines.markersize'] = 10.0
    # For Poster/Presentation
    mpl.rcParams['ytick.labelsize'] = 32
    mpl.rcParams['xtick.labelsize'] = 32
    mpl.rcParams['font.size'] = 28
    main(sys.argv[1:])
