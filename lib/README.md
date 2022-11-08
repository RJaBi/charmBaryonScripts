# Python

### Here I keep generic python helper functions
### That apply to multiply smaller programs
---
### myModules
*A couple of very generic functions used by virtually every file*
 * initBigPlotSettings()
   * * Initialise a bunch of plot settings that make the plots look nicer
 * refineXYLims(params, subDict: 'analysis') -> Dict[str, Any]
   * * Toml doesn't support None type. This converts any string 'None' value in the params[subDict] dictionary to be None type. just does params if subDict=None
 * markers
   * * defines a list of marker types for matplotlib
 * GetArgs(args: list) -> MutableMapping[str, Any]
   * * Loads the toml file in the arguments
 * removeZero(data: np.ndarray) -> Tuple[np.ndarray, List[int]]:
   * * Removes any occasions where the data is zero
   * * returns the removed and also indices which were kept
 * replaceZero(data: np.ndarray, rep: float)
   * * Similarly but replaces the zero with the float

### myIO
*More specialised code for loading correlation functions*
 * def strToLambda(string)
   * * Converts a specifically formatted string to a lambda function
   * * I.e.
   * * * outVar:var1,,,...,_N: f(var1,,,...,_N)where f(...) must be an appropriate maths function
 * performVarAnalysis(GList: List[np.ndarray], n: int, t0: int, dt: int, src_t0: int) -> List[np.ndarray]
   * * Performs a variational analysis on the correlators in GList
   * * Assumes order is G11, G12, G21, G22, etc
   * * Uses the correlator at t0 + dt to define the reoccurance relation
 * loadSingleMeson(cfFile):
   * * Loads a single meson correlator in gencf format
   * * * i.e. Nt * double precision complex in big endian. Nt=0 first
 * setParityMask(proj: str) -> np.ndarray:
   * * Sets the appropriate parity mask for the 4x4 dirac matrix in column major order
 * loadGencfMomBaryon(cfFile: str, nt: int, numMom: int, proj: str = '+', mom: int = 0):
   * * Loads a single momenta of a gencf baryon with no lorentz indices
   * * A gencf baryon is complex, double precision, big endian, with a 4x4 dirac matrix at each time step, in column major ordering
 * loadOQCDMeson(cfFile: str, gOff: int = 5):
   * * Loads an openqcd meson correlator which does not have momenta
   * * gOff 5 is gamma 5
   * * gOff 2,3,4 are the vector components
 * def initCorrelators(params):
   * * Loads all correlators in params['cfuns'] in
   * * Does the maths, etc
 * def labelNorm(G: np.ndarray, params):
   * * Specifically for normallising the correlator data [ncon, Nt] at the labelToAnalyse level
   * * Currently only supports by mean at src t0

### myEffE
*Taking an effective mass (energy) of correlators*

 * effE_centre(massdt: np.float64, GJ2: np.ndarray) -> np.ndarray:
   * * The centre finite difference for effective energy
   * * i.e. an arccosh
 * effE_forward(massdt: np.float64, GJ2: np.ndarray) -> np.ndarray:
   * * the foward finite difference effective mass
   * * (1.0/massdt) * log(G(t)/G(t+massdt))
 * effE_solve(massdt: np.float64, GJ2: np.ndarray) -> np.ndarray:
   * * Solves the meson correlator (cosh) for the mass
 * effE_barAna(GJ2: np.ndarray) -> np.ndarray:
   * * Implements the 4 point effective mass for a baryon with periodic BC
   * * Gets the positive parity solution
 * getEffE(params: Dict[str, Any], massdt: np.float64, GJ2: np.ndarray) -> np.ndarray:
   * * Calls the above function based on the effEMethod in params['analysis']
   * * changes any nans to zeros

### myGVPlotter
*Functions for plotting gvar's using matplotlib*
 * plot_gvEbar(x: np.ndarray, y: np.ndarray, ax, ma=None, ls=None, lab=None, col=None, alpha=1):
   * * Makes plotting an errorbar plot of a gvar easier
 * myHSpan(y: np.ndarray, ax, ls='--', colour=None, lab=None, alpha=1.0):
   * * Horizontal spans of gvars
 * myVSpan(y: np.ndarray, ax, ls='--', colour=None, lab=None, alpha=1.0):
   * * Vertical spans of gvars
 * myFill_between(x: np.ndarray, y: np.ndarray, ax, ma=None, ls=None, lab=None, alpha=0.15, colour=None):
   * * Fill_between for gvar
 * colours
   * * A variable containing a list of the 10 tabuleau 10 colours

### myFits
*Functions used in multiple fitting programs*
 * makeFitWins(fitMin: int, fitMax: int, Nt: int, twoSided=True):
   * * Makes all fit windows from fitMin to fitMax
 * plotWinLines(ax, wins: List[str]):
   * * Plots a vertical line where the start of the fit window changes
   * * Assumes wins is sorted
 * getWinTicks(ax, wins: List[str]) -> List[str]:
   * * Gets a list of tick labels where only the first t0 is labeleld
   * * Vertically offsets every second one
   * * Assumes wins is sorted
