## AUTHORS

Nicholas Milson, Caroline Barton, Dr. Philip Bennett

## INSTALLATION

In your command line, call _pip install BinaryStarSolver_.
Then, in your python code, write: _from binarystarsolve.binarystarsolve import StarSolve_

## WEB APPLICATION
binarystarsolver.com

## DESCRIPTION

Given a series of radial velocities as a function of time for a star in a binary system, this program solves for various orbital parameters. Namely, it solves for eccentricity (e), argument of periastron (ω), velocity amplitude (K), long term average radial velocity (γ), time of periastron (T0), and orbital period (P).

If the orbital parameters of a primary star are already known, the program can also find the orbital parameters of a companion star, with only a few radial velocity data points.

In the case of a double-lined binary where both stars have identical orbital coverage, the program can solve for the orbital parameters of both stars at once using a single data file.

The program can also solve for both stars at once when the primary and companion radial velocities were observed at different times or have unequal orbital coverage. In this case, the primary and companion radial velocity data are supplied separately.

Note that K = n * a1 * sin(i) / sqrt(1-e²), where a1 is the primary star's semi-major axis, n is the mean motion (n = 2π / P), and i is the inclination angle of the orbit.

The equation the data is being fitted to is V(t) = γ + K*(cos(v(t) + ω) + e * cos(ω)), where v(t) is the true anomaly (angle from periastron).

For more information on the implementation of this code, refer to the paper "A Python Code to Determine Orbital Parameters of Spectroscopic Binaries": https://arxiv.org/abs/2011.13914

## USER INSTRUCTIONS

All the functionality of this program is accessed via the function StarSolve(). Whether the parameters being solved for are for a primary star or a companion star (or both at once), StarSolve() is to be used.

The first two arguments of the function StarSolve() are data_file and star. data_file may be a string containing either the name of the txt file (if the file is in your working directory) or the file path (if the file is not in the working directory). Alternatively, the radial velocity data may be supplied directly as a two-dimensional array.

star is a string, where the options are "primary", "secondary", "both", or "both_unequal". Note that star is not case sensitive.

The option star = "both" without a second file is used when the primary and companion stars have identical observation times. The unequal-coverage mode may be selected using star = "both_unequal", or by using star = "both" and supplying the companion data with the keyword argument file2.

### PRIMARY STAR APPLICATION INSTRUCTIONS

For solving the parameters of the primary star, this program works for two general types of radial velocity data sets. It works for data sets where all the data is from a singular deduction of the radial velocities. It also works for data sets which are composed of sub data sets, where each sub data set's radial velocities were deduced separately from the other subsets. Having a data set composed of subsets with radial velocities deduced separately may result in discrepancies of the long term average radial velocity (γ) between the subsets. This program can deal with this discrepancy by choosing one of the sub data sets' γ as the "true" γ, and shifting the other sub data set's radial velocities to match up with the chosen subset.

For the case where all the data is from a singular deduction of the radial velocities, the user must supply their radial velocity data as a tab separated txt file or array, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, and the third (optional) value being the weight.

For the case with the RV data set composed of multiple sub data sets, the user must supply their radial velocity data as a tab separated txt file or array, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, the third value being the weight, and the fourth value being an integer signifying which sub data set the data comes from. Note, if the user does not wish to assign weights to the data, the third column still must be filled and thus should be a column of ones.

If the period of the orbit is already known, the user may pass in a float with the keyword argument Period (with the period in days). If an estimate of the period is known, but the user would still like the period to be solved for, the user may pass in a float with the keyword argument Pguess (in days).

If neither Period nor Pguess is supplied, the program makes several initial estimates of the period and selects the estimate that gives the lowest initial sum of squared residuals. If the data contains an exceptionally large gap, the program also attempts to estimate the period using the largest continuous cluster of observations.

If the user would like the estimated covariance matrix returned, the boolean keyword argument covariance may be passed in as True. If the period is previously known (and passed in using Period), the covariance matrix returned will be a 5 by 5 matrix. Otherwise, it will be a 6 by 6 matrix.

By default, when star = "primary", StarSolve() creates two plots of the fit (one for RV vs RJD and one for RV vs phase). If the user would not like these graphs shown, they may use the optional boolean keyword graphs and set it as False.

For circular orbits (e = 0), the argument of periastron (ω) becomes undefined. If the user has reason to believe that the eccentricity is 0, they may pass in the optional boolean keyword parameter zeroEcc as True. This will set e and ω to 0 for the minimization, making the minimization run faster. It is also acceptable to have zeroEcc = False for circular orbits, but just note that the ω returned will not mean anything physically.

StarSolve will return a list of the solved parameters (along with asin(i) and f(M)), in the order [γ, K, ω, e, T0, P, asin(i), f(M)], followed by a list of uncertainties associated with each parameter (and asin(i) and f(M)). As previously stated, the keyword argument covariance may be passed as True so that the covariance matrix is returned too.

Warning: If no Period or Pguess is provided, StarSolve() cannot find reliable results if the RV data supplied does not span at least approximately 1.5 periods and may fail if the data spans less than one full period. Furthermore, if the data set is composed of multiple sub data sets with different γ velocities, the program may fail to correct the offset in γ velocities if none of the sub data sets span at least 1.5 periods. Using a period determined by other methods (e.g. photometrically) is often preferable to allowing the program to solve for the period. The convergence of the minimization is fairly sensitive to initial period estimates.

**Examples**:

```python
params, err = StarSolve(
    data_file="myRVdata.txt",
    star="primary"
)
```

```python
params, err, cov = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    covariance=True
)
```

```python
params, err = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Pguess=7430,
    graphs=False
)
```

### COMPANION STAR APPLICATION INSTRUCTIONS

The RV data must be formatted the same way for star = "secondary" as for star = "primary", i.e. as a tab separated txt file or array, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity in km/s, and the third (optional) value being the weight.

To find the parameters of a companion star, a list of the known parameters of the primary star, in the order [γ, K, ω, e, T0, P] (with ω in degrees), must be passed in using the keyword argument X.

If the user wishes to return the error on the parameters of the companion star, they must also input a list of errors of the known parameters of the primary star (in the order [γ, K, ω, e, T0, P]), using the keyword err.

Due to various observational reasons, there is often a discrepancy between the average radial velocity of the primary star's γ1 and the companion star's γ2, though both should be equal. γ1 is assumed to be the "correct" γ and is the γ returned, but for purposes of fitting the curve to the data, a γ2 is solved for. Using the boolean keyword shift = True, the discrepancy between average velocities (γ2 - γ1) is returned.

By default, when star = "secondary", StarSolve() creates two plots of the fit (one for RV vs RJD and one for RV vs phase). If the user would not like these graphs shown, they may use the optional boolean keyword graphs and set it as False.

**Examples**:

```python
params = StarSolve(
    data_file="companion.txt",
    star="secondary",
    X=[-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68]
)
```

```python
params, err = StarSolve(
    data_file="companion.txt",
    star="secondary",
    X=[-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68],
    err=[0.025298, 0.034264, 0.613298, 0.0020712, 5.80324, 0]
)
```

```python
params, shift = StarSolve(
    data_file="companion.txt",
    star="secondary",
    X=[-6.4354, 13.956, 203.991, 0.20544, 56108.8, 3770.68],
    shift=True,
    graphs=False
)
```

### PRIMARY AND SECONDARY STARS TOGETHER WITH IDENTICAL COVERAGE

Solving the parameters of both stars at once using the keyword star = "both" without supplying file2 only works for data sets where, for every time t, there is both a V1(t) point and a V2(t) point.

Like with the star = "primary" case, the program works for two general types of radial velocity data sets. It works for data sets where all the data is from a singular deduction of the radial velocities. It also works for data sets which are composed of sub data sets, where each sub data set's radial velocities were deduced separately from the other subsets (so the γ's may be different). For the both stars case, this program deals with this discrepancy by finding γ for each sub data set, performing a weighted average to get a value of γ, and then shifting all the data sets to said weighted average value of γ.

For the case where all the data is from a singular deduction of the radial velocities, the user must supply their radial velocity data as a tab separated txt file or array, with the first value being the time in Reduced Julian Date (RJD), the second being the radial velocity of the primary star in km/s, the third being the radial velocity of the secondary star in km/s, and the fourth and fifth (optional) values being the respective weights of the primary and secondary star's radial velocities.

For the case with the RV data set composed of multiple sub data sets, the user must supply their radial velocity data as a tab separated txt file or array, with the first value being the time in Reduced Julian Date (RJD), the second being the primary star's radial velocity in km/s, the third value being the companion star's radial velocity in km/s, the fourth and fifth values being the respective weights of the primary and secondary star's radial velocities, and the sixth value being an integer signifying which sub data set the data comes from. Note, if the user does not wish to assign weights to the data, the fourth and fifth columns still must be filled and thus should be columns of ones.

If the period of the orbit is already known, the user may pass in a float with the keyword argument Period (with the period in days). If an estimate of the period is known, but the user would still like the period to be solved for, the user may pass in a float with the keyword argument Pguess (in days).

If the user would like the estimated covariance matrices returned, the boolean keyword argument covariance may be passed in as True. This returns two covariance matrices (one for each star, with the primary star first). If the period is previously known and passed in using Period, the covariance matrices returned will be 5 by 5 matrices. Otherwise, they will be 6 by 6 matrices.

By default, when star = "both", StarSolve() creates four plots of the fits: one for RV vs RJD of the primary star, one for RV vs phase of the primary, one for RV vs RJD of the companion star, and one for RV vs phase of the companion. If the user would not like these graphs shown, they may use the optional boolean keyword graphs and set it as False.

For circular orbits (e = 0), the argument of periastron (ω) becomes undefined. If the user has reason to believe that the eccentricity is 0, they may pass in the optional boolean keyword parameter zeroEcc as True. This will set e and ω to 0 for the minimization, making the minimization run faster. It is also acceptable to have zeroEcc = False for circular orbits, but just note that the ω returned will not mean anything physically.

When star = "both", StarSolve will return a list of the solved parameters for each star (along with asin(i) and f(M)), in the order [γ, K, ω, e, T0, P, asin(i), f(M)], followed by a list for each star of uncertainties associated with each parameter (and asin(i) and f(M)). As previously stated, the keyword argument covariance may be passed as True so that the covariance matrices are returned too. Note that for all these lists, the primary star's values are returned first.

**Examples**:

```python
params, err, cov = StarSolve(
    data_file="rvData.txt",
    star="both",
    covariance=True,
    zeroEcc=True
)
```

```python
params, err = StarSolve(
    data_file="rvData.txt",
    star="both",
    Pguess=33.8,
    graphs=False
)
```

### PRIMARY AND SECONDARY STARS TOGETHER WITH UNEQUAL COVERAGE

If the primary and companion stars were observed at different times, or contain different numbers of radial velocity measurements, the unequal-coverage mode may be used.

The primary star data are supplied using data_file, while the companion star data are supplied separately using the keyword argument file2. The mode may be selected using star = "both_unequal". Alternatively, using star = "both" while supplying file2 selects the same mode.

The primary and companion data files must each contain the observation time in Reduced Julian Date (RJD), followed by the radial velocity in km/s and an optional weight:

```
time    radial_velocity    weight
```

The two files do not need to contain the same observation times or the same number of measurements.

The fit solves simultaneously for γ1, K1, ω1, e, T0, P, γ2, and K2. The two stars share the same eccentricity, time of periastron, and period. The argument of periastron of the companion is constrained to differ from that of the primary by 180 degrees.

By default, γ1 and γ2 are fitted independently. If the user would like the two average radial velocities to be forced equal, the optional boolean keyword argument gam2 may be set to True.

Period, Pguess, covariance, graphs, and zeroEcc have the same meanings as in the other fitting modes. If covariance=True, a covariance matrix for the simultaneous fit is returned.

StarSolve returns the solved values in the order:

```
[γ1, K1, ω1, e, T0, P, γ2, K2, asin(i), f(M)]
```

This is followed by a list containing the associated uncertainties in the same order.

By default, four plots are created: RV vs RJD and RV vs phase for both the primary and companion stars.

**Examples**:

```python
params, err = StarSolve(
    data_file="primary.txt",
    file2="companion.txt",
    star="both_unequal",
    Period=33.8
)
```

The same mode can also be selected using star = "both":

```python
params, err = StarSolve(
    data_file="primary.txt",
    file2="companion.txt",
    star="both",
    Pguess=33.8,
    graphs=False
)
```

To force γ1 and γ2 to be equal:

```python
params, err, cov = StarSolve(
    data_file="primary.txt",
    file2="companion.txt",
    star="both_unequal",
    Period=33.8,
    gam2=True,
    covariance=True
)
```

### GENERATING MODEL RADIAL VELOCITIES

For any fitting mode, the fitted model radial velocities may optionally be returned using the keyword argument model_output.

If model_output = "input", model velocities are returned at the original observation times:

```python
params, err, model = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    graphs=False,
    model_output="input"
)
```

If model_output = "phase", model velocities are returned on a grid of orbital phases:

```python
params, err, model = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    graphs=False,
    model_output="phase"
)
```

By default, the phase grid runs from 0 to 1 in increments of 0.001. A custom phase grid may be passed using the keyword argument phase_grid as a tuple containing:

```
(minimum phase, maximum phase, phase increment)
```

For example:

```python
params, err, model = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    graphs=False,
    model_output="phase",
    phase_grid=(0, 1, 0.01)
)
```

An explicit array of phases may also be passed:

```python
import numpy as np

phases = np.linspace(0, 1, 101)

params, err, model = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    graphs=False,
    model_output="phase",
    phase_grid=phases
)
```

The returned model table contains three columns:

```
[time, phase, fitted radial velocity]
```

For either of the two-star applications, model is a tuple containing a separate table for the primary and companion stars:

```python
primary_model, companion_model = model
```

The model output is appended to the normal return values. For example, if covariance=True:

```python
params, err, cov, model = StarSolve(
    data_file="myRVdata.txt",
    star="primary",
    Period=3784.3,
    covariance=True,
    graphs=False,
    model_output="input"
)
```

If there are any questions, please contact me at <ins>[nmilson@ualberta.ca](mailto:nmilson@ualberta.ca)</ins>.

