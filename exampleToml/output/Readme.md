## EXAMPLE OUTPUTS

These output files are examples only and are not representative of any final analysis. They are often on small subsets of the data.



### mass_sysErr.csv
* anaDir
  * The directory where the fits used are saved
* Nt
  * The number of time slices for the correlator
* Operator
  * The name of the operator in openqcd-hadspec form
* Hadron
  * Which hadron it is
* EP
  * The positive parity energy from the fit from the chosen method (EPMethod)
* EPMethod
  * What the chosen method is
* EP_Sys
  * The systematic error from variation of fit method
* EPSys
  * The positive parity energy from the chosen fit method with the systematic error added in quadrature
* EPFtypes
  * The other fit methods considered for the systematic
* EPFits
  * The values of these other methods. Matches order of EPFtypes
And similarly for the negative parity energy from the fit
* EM
* EMMethod
* EM_Sys
* EMSys
* EMFtypes
* EMFits