AUTHORS

Caroline Barton, Dalhousie University & Nicholas Milson, Dalhousie University

SUPERVISOR

Dr. Philip Bennett, Dalhousie University

INSTALLATION

In your command line, call pip install BinaryStarSolver.
Then, in your python code, write: from binarystarsolve.binarystarsolve import StarSolve

DESCRIPTION

Given a series of radial velocities as a function of time for a star in a binary system, this program solves for various orbital parameters. Namely, it solves for eccentricity (e), argument of periastron (ω), velocity amplitude (K), long term average radial velocity (γ), and orbital period (P).

If the orbital parameters of a primary star are already known, the program can also find the orbital parameters of a companion star, with only a few radial velocity data points.

Note that K = n * a1 * sin(i) / ((1-e)^0.5), where a1 is the primary star's semi-major axis, n is the mean motion (n = 2π / P), and i is the inclination angle of the orbit.

The equation the data is being fitted to is V(t) = γ + K*(cos(v(t) + ω) + e*cos(ω)), where v(t) is the true anomaly (angle from periastron).

For more information on the implementation of this code, refer to the paper "A Numerical Determination of the Orbital Parameters of Spectroscopic Binary Stars" *** (working title) *** . The paper has been submitted for publication. A link to the paper on the arXiv should be posted very soon. Please check back soon.

USER INSTRUCTIONS

All the functionality of this program is accessed via the function StarSolve(). Whether the parameters being solved for are for a primary star or a companion star, StarSolve() is to be used.

The first two arguments of the function StarSolve() are data_file and star. data_file a string, which is either the name of the txt file (if the file is in your working directory), or the file path (if the file is not in the working directory). star is also a string, where the two options are "primary" or "secondary". star is not case sensitive.

PRIMARY STAR APPLICATION INSTRUCTIONS

For solving the parameters of the primary star, this program works for two general types of radial velocity data sets. It works for data sets where all the data is from a singular deduction of the radial velocities. It also works for data sets which are composed of sub data sets, where each sub data set's radial velocities were deduced separately from the other subsets. Having a data set composed of subsets with radial velocities deduced separately may result in discrepancies of the long term average radial velocity (γ) between the subsets. This program can deal with this discrepancy by choosing one of the sub data set's γ as the "true" γ, and shifting the other sub data set's radial velocities to match up with the chosen subset.

For the case where all the data is from a singular deduction of the radial velocities, the user must supply their radial velocity data they wish to fit as a tab separated txt file, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, and the third (optional) value being the weights.

For the case with the RV data set composed of multiple sub data sets, the user must supply their radial velocity data they wish to fit as a tab separated txt file, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, and the third value being the weights, and the fourth value being an integer signifying which sub data set the data comes from. Note, if the user wishes to not assign weights to the data, the third column still must be filled, thus should be a column of ones.

If the period of the orbit is already known, the user may pass in a float with the keyword argument Period (with the period in days). If the period of the orbit is already known, but the user would still like the period to be solved for, the user may pass in a float with the keyword argument Pguess (in days).

If the user would like the estimated covariance matrix returned, the boolean keyword argument covariance may be passed in as True. If the period is previously known (and passed in using Period), the covariance matrix returned will be a 5 by 5 matrix. Otherwise, it will be a 6 by 6. By default, when star = "primary", StarSolve() creates two plots of the fit (one for RV vs RJD, one for RV vs phase). If the user would not like these graphs shown, they may use the optional boolean keyword graphs, and set it as False.

StarSolve will return a list of the solved parameters (along with asini and f(M)), in the order [γ, K, ω, e, T0, P, asin(i), f(M)], followed by a list of uncertainties associated with each parameter (and a*sini and f(M)). As previously stated, the arguement keyword covariance may be passed as True, so that the covariance matrix is returned too

Warning: If no Period or Pguess is provided, StarSolve() cannot find reliable results if the RV data supplied does not span at least 1.5 periods and will fail if the data spans less than one full period. Furthermore, if the data set is composed of multiple sub data sets with different γ velocities, the program may fail to correct the offset in γ velocities if none of the sub data sets span at least 1.5 periods. Using a period determined by other methods (e.g. photometrically) is often preferable to allowing the program to solve for the period. The convergence of the minimization is fairly sensitive to initial period estimates.

Examples:

params, err = StarSolve(data_file = "myRVdata.txt", star = "primary")

params, err, C = StarSolve(data_file = "myRVdata.txt", star = "primary", Period = 3784.3, covariance = True)

params, err = StarSolve(data_file = "myRVdata.txt", star = "primary", Pguess = 7430, graphs = False)

COMPANION STAR APPLICATION INSTRUCTIONS

The RV data must be formatted the same way for star = "secondary" as for star = "primary", i.e. as a tab separated txt file, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, and the third (optional) value being the weights.

To find the parameters of a companion star, a list of the known parameters of the primary star, in the order [γ, K, ω, e, T0, P] (with ω in degrees), must be passed in using the keyword arguement X.

If the user wishes to return the error on the parameters of the companion star, they must also input a list of errors of the known parameters of the primary star (in the order [γ, K, ω, e, T0, P]), using the keyword err.

Due to various observational reasons, there is often a discrepancy between the average radial velocity of the primary's γ1, and the companion star's γ2, though both should be equal. γ1 is assumed to be the "correct" γ, and is the γ returned; but for purposes of fitting the curve to the data, a γ2 is solved for. Using the boolean keyword shift = True, the discrepancy between average velocities (γ2 - γ1) is returned.

By default, when star = "secondary", StarSolve() creates two plots of the fit (one for RV vs RJD, one for RV vs phase). If the user would not like these graphs shown, they may use the optional boolean keyword graphs, and set it as False.

Examples:

params = StarSolve(data_file = "companion.txt", star = "secondary", X = [-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68])

params,err = StarSolve(data_file = "companion.txt", star = "secondary", X = [-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68], err = [0.025298, 0.034264, 0.613298, 0.0020712, 5.80324, 0])

params,shift = StarSolve(data_file = "companion.txt", star = "secondary", X = [-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68], shift = True, graphs = False)

If there are any questions, please contact nick.milson@dal.ca
