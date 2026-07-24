"""
All of the functionality of this package is through the function StarSolve().
For a more detailed description of how to use this package, please read the ReadMe.
https://github.com/NickMilsonPhysics/BinaryStarSolver/blob/master/README.md
"""

# Import statements
import math
import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np


class ConvergenceError(RuntimeError):
    """Raised when an iterative orbital fit fails to converge."""

def StarSolve(data_file, star="primary", Period=None, Pguess=None, covariance=False, graphs=True, zeroEcc=0, X=None, err=None, shift=False, file2=None, gam2=False, model_output=None, phase_grid=None):
    """
    Please first read the read me at:
   	https://github.com/NickMilsonPhysics/BinaryStarSolver/blob/master/README.md

    Keyword Arguments:

    ** Note some arguments pertain to when parameters of the primary are
    being solved for, some pertain to when parameters of the secondary
    are being solved for, and some pertain to the case when both stars are
    being solved for at once. If an argument doesn't pertain to the situation
    this function is being used for, it may be ignored, as all the case
    specific arguments are None by default. **

    Keyword Arguments mandatory regardless of situation:
    data_file - String holding the name of the tab separated txt file, or an
                array containing the radial velocity data.
    star - Can be "primary", "secondary", "both", or "both_unequal". Is
           "primary" by default. star is a string for if the parameters of
           the primary or secondary are being solved for (or both). Supplying
           file2 while star = "both" also selects the unequal temporal
           coverage application.


    Keyword Arguments pertaining to primary star application:
    Period (optional) - Known value for period (if none provided, by default
                        PrimarySolve makes its own guess).
    Pguess (optional) - If period is not known precisely, but still an
                        estimate is available.
    covariance (optional) - If covariance = True, PrimarySolve will return the
                            entire covariance matrix.
    graphs (optional) - If graphs = true, PrimarySolve creates 2 plots. First, a
                        plot of the radial velocity data as a function of time
                        (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                        where [γ, K, ω, e, T0] are the orbital parameters
                        solved for by PrimarySolve(). Second, a plot of the
                        radial velocity data as a function of phase (between 0
                        and 1) of the star's orbit. Again, fit with
                        V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                        function of phase. On both plots, the long-term average
                        velocity, γ, is also shown.
    zeroEcc (optional) - If zeroEcc = True, eccentricity and arguement of
                         periastron are kept at 0

    Keyword Arguments pertaining to both stars at once application:
    Period (optional) - Known value for period (if none provided, by default
                        PrimarySolve makes its own guess).
    Pguess (optional) - If period is not known precisely, but still an
                        estimate is available.
    covariance (optional) - If covariance = True, PrimarySolve will return the
                            entire covariance matrix.
    graphs (optional) - If graphs = true, BothStars creates 4 plots. First, a
                        plot of the radial velocity data as a function of time
                        (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                        where [γ, K, ω, e, T0] are the orbital parameters
                        solved for by BothStars(). Second, a plot of the
                        radial velocity data as a function of phase (between 0
                        and 1) of the star's orbit. Again, fit with
                        V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                        function of phase. On both plots, the long-term average
                        velocity, γ, is also shown. Then these two plots are
                        made but for the second star. (so K and omega will
                        be different).
    zeroEcc (optional) - If zeroEcc = True, eccentricity and arguement of
                         periastron are kept at 0

    Keyword Arguments pertaining to both stars with unequal temporal coverage:
    file2 - String holding the name of the tab separated txt file containing
            the companion star radial velocity data, or an array containing
            the companion data. data_file contains the primary star data.
            Supplying file2 while star = "both", or using
            star = "both_unequal", selects this application.
    gam2 (optional) - If gam2 = True, the average radial velocities of the
                      primary and companion stars are forced to be equal.
    Period, Pguess, covariance, graphs, and zeroEcc have the same meaning as
    in the both stars at once application. The primary and companion data may
    contain different numbers of observations taken at different times.

    Keyword Arguments pertaining to secondary star application:
    X - List of elements of the primary star, in the order [γ, K, w, e, T0, P].
        (with w in degrees). K does not necessarily need to be known but a
        value must be input.
    err (optional) - List of errors on the elements in X (in the same order)
    shift (optional) - If shift = True, returns discrepancy between γ of primary
    graphs (optional) - If graphs = True, CompanionSolve creates 2 plots. First, a
                        plot of the radial velocity data as a function of time
                        (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                        where [γ, K, ω, e, T0] are the orbital parameters
                        solved for by CompanionSolve(). Second, a plot of the
                        radial velocity data as a function of phase (between 0
                        and 1) of the star's orbit. Again, fit with
                        V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                        function of phase. On both plots, the long-term average
                        velocity, γ, is also shown.

    Keyword Arguments pertaining to generated model data:
    model_output (optional) - Can be "input" or "phase". If model_output =
                              "input", model radial velocities are generated at
                              the original observation times. If model_output =
                              "phase", model radial velocities are generated on
                              a phase grid.
    phase_grid (optional) - Phases used when model_output = "phase". May be an
                            array of phases or a tuple containing
                            (minimum phase, maximum phase, phase increment).
                            By default, phases from 0 to 1 in increments of
                            0.001 are used.


    Returns:

    ** Note which variables are returned depends on the application of this
    function. **

    Returns from the primary star (and both stars) application:
    x - List(s) of solved orbital parameters, the semi major axis and the
        mass function [γ, K, ω, e, T0, P, a, f(M)]
    err - List(s) of standard errors on the determination of the orbital elements,
          in the order [γ, K, ω, e, T0, P, a, f(M)]
    C - Estimated covariance matrix (or matrices) (if argument covariance = True)

    Returns from the both stars with unequal temporal coverage application:
    x - List of solved orbital parameters in the order
        [γ1, K1, ω1, e, T0, P, γ2, K2, a1, f(M1)]
    err - List of standard errors on the determination of the orbital elements,
          in the same order as x
    C - Estimated covariance matrix (if argument covariance = True)

    Returns from the secondary star application:
    x - List of solved elements of the companion star, in the order
        [γ, K, w, e, T0, P, asin(i), f(m)].
    err (optional) - If errors associated with primary star are supplied, list
                     of errors on the elements in x
    shift (optional) - Difference in average radial velocity of primary star, γ1,
                       and of secondary star, γ2 (if argument shift = True). The
                       radial velocity of the primary star is taken to be the
                       'correct' γ in x. Thus shift returned = γ2 - γ1.

    Returns from generated model data:
    model - If model_output is supplied, the generated model data are appended
            to the usual return values. Each model table has columns
            [time, phase, fitted radial velocity]. For a two star application,
            model is a tuple containing the primary and companion tables.
    """

    #%% SET UP

    def load_data(source):
        """Load a numerical table from a path or accept an array directly."""
        if isinstance(source, (str, bytes, os.PathLike)):
            loaded = np.loadtxt(source, ndmin=2)
        else:
            loaded = np.asarray(source, dtype=float)
            if loaded.ndim == 1:
                loaded = np.atleast_2d(loaded)

        if loaded.ndim != 2 or loaded.shape[0] == 0:
            raise ValueError("Radial-velocity data must be a non-empty two-dimensional table.")
        if loaded.shape[1] < 2:
            raise ValueError("Each dataset must contain at least time and radial velocity columns.")
        if not np.isfinite(loaded).all():
            raise ValueError("Radial-velocity data contain NaN or infinite values.")
        return loaded.copy()

    star = star.lower() # making star not case sensitive
    if star in ("both_unequal", "unequal", "secondary_both"):
        star = "secondary_both"
    elif star == "both" and file2 is not None:
        # Backward compatible: the old one-file star="both" mode is unchanged,
        # while supplying file2 naturally selects the new unequal-coverage mode.
        star = "secondary_both"

    # sorting data chronologically
    data = load_data(data_file)
    data_sorted = np.array(sorted(data,key=lambda l:l[0]))
    time = data_sorted[:,0] # RJD times
    V = data_sorted[:,1] # radial velocities

    N = len(time) # number of data points

    if file2 is not None and star == "secondary_both":
        data2 = load_data(file2)
        data_sorted2 = np.array(sorted(data2,key=lambda l:l[0]))
        time2 = data_sorted2[:,0] # RJD times
        N2 = len(time2)
        V2 = data_sorted2[:,1] # radial velocities
        if data2.shape[1] >= 3:
            weight2 = data_sorted2[:,2]
        else:
            weight2 = np.ones(N2)

    if star != "both":
        if X != None:
            Period = X[5]
        # weights
        if data.shape[1] >= 3:
            weight = data_sorted[:,2]
        else:
            weight = np.ones(N)
        if data.shape[1] >= 4:
            setNum = data_sorted[:,3]
    #%% Initial guess

    # estimates of period
    def PGuess(test = False, timeData = None, VData = None):
        """

         Keyword Arguments:
         test (false by default) - If true, PGuess is used as a check to see if
                                   finding a period is possible. If false, PGuess
                                   estimates a period or exits the program.

        When test = False, PGuess returns two estimates of the period, found
        using different methods.

        If less than one period is observed, will call sys.exit signaling that not
        enough data is provided to make an initial estimate of the period. In this
        case the user must provide a known period to StarSolve().

        For both methods, begins by finding all pairs of points where data crosses
        the average radial velocity axis.

        The first method provides an estimate by finding the two points where the
        data crosses the average radial velocity axis, seperated by an integer
        number of periods, that have the smallest difference in radial velocity.

        The period is found from the times (RJD) that these points occur, divided
        by the number of periods apart they are.


        The second method provides an estimate by averaging the separation of
        points crossing the average radial velocity axis by one period, and the
        maxima and minima in between these points.
        """

        if timeData is None:
            timeP = time
            VP = V
        else:
            timeP = np.asarray(timeData)
            VP = np.asarray(VData)
        NP = len(timeP)

        # find all cross points
        # PGuess now uses the actual length of the supplied velocity array, so
        # no special factor-of-two correction is needed for double-lined data.
        Vavg = sum(VP)/NP
        low = None
        crossid = []
        for i in range(1, NP-1):
            if VP[i+1] > Vavg and VP[i-1] < Vavg or VP[i+1] < Vavg and VP[i-1] > Vavg:
                crossid.append(i)

        if NP <=2:
            if test:
                return np.nan
            sys.exit()

        # corrects for double crossings
        l = len(crossid)
        for i in range(l-1):
            a = crossid[i]
            b = crossid[i+1]
            if a+1 != b:
                good = False
                for j in range(a, b):
                    if abs(VP[j] - Vavg) > (max(VP) - Vavg) / 4:
                        good = True
                        break
                if not good:
                    if abs(VP[a] - Vavg) < abs(VP[b] - Vavg):
                        crossid[i+1] = a + 1
                    else:
                        crossid[i] = b - 1

        # if data covers less than one period
        if len(crossid) == 0:
            if test:
                return np.nan
            sys.exit()
        if crossid[0] + 1 != crossid[1]:
            if l < 4:
                if test:
                    return np.nan
                sys.exit()
        else:
            if l < 5:
                if test:
                    return np.nan
                sys.exit()

        crosstime = []
        for i in crossid:
            if (VP[i-1] - Vavg) * (VP[i] - Vavg) <= 0 and VP[i] != VP[i-1]:
                tcross = timeP[i-1] + (Vavg - VP[i-1]) * (timeP[i] - timeP[i-1]) / (VP[i] - VP[i-1])
            elif (VP[i] - Vavg) * (VP[i+1] - Vavg) <= 0 and VP[i+1] != VP[i]:
                tcross = timeP[i] + (Vavg - VP[i]) * (timeP[i+1] - timeP[i]) / (VP[i+1] - VP[i])
            elif VP[i+1] != VP[i-1]:
                tcross = timeP[i-1] + (Vavg - VP[i-1]) * (timeP[i+1] - timeP[i-1]) / (VP[i+1] - VP[i-1])
            else:
                tcross = timeP[i]
            crosstime.append(tcross)

        # excludes last two sets of points from finding period
        if crossid[l-1] - 1 == crossid[l-2]:
            r = l - 4
            check = False
        else:
            r = l - 5
            check = True

        Psum = 0 # sum of periods for method 2
        n = 0 # number of periods measured for method 2

        # for method 1, finds difference in radial velocity using every pair of
        # points seperated by an integer number of periods along Vavg
        # simoultaneously adds time between all pairs of cross points seperated by
        # one period to Tsum for method 2
        first = True
        for i in range(r):
            # method 1: checks difference using every point 4n away
            # method 2: adds time to 4th point ahead to Psum
            Psum += crosstime[i+4] - crosstime[i]
            for j in range(4, l-i, 4):
                if first:
                    low = abs(VP[crossid[4]] - VP[crossid[0]])
                    besta = crosstime[0]
                    bestb = crosstime[4]
                    p = 1 # number of periods between points
                    first = False
                if abs(VP[crossid[i]] - VP[crossid[i+j]]) < low:
                    low = abs(VP[crossid[i]] - VP[crossid[i+j]])
                    besta = crosstime[i]
                    bestb = crosstime[i+j]
                    p = j/4
            # if cross point is first in pair:
            # method 1: checks every point 5 + 4n away
            # method 2: adds time to 5th point ahead to Psum
            if crossid[i] + 1 == crossid[i+1]:
                for j in range(5, l-i, 4):
                    if abs(VP[crossid[i]] - VP[crossid[i+j]]) < low:
                        low = abs(VP[crossid[i]] - VP[crossid[i+j]])
                        besta = crosstime[i]
                        bestb = crosstime[i+j]
                        p = (j-1)/4
                Psum += crosstime[i+5] - crosstime[i]
            # if cross point is second in pair:
            # method 1: checks every point 3 + 4n away
            # method 2: adds time to 3rd point ahead to Psum
            else:
                for j in range(3, l-i, 4):
                    if abs(VP[crossid[i]] - VP[crossid[i+j]]) < low:
                        low = abs(VP[crossid[i]] - VP[crossid[i+j]])
                        besta = crosstime[i]
                        bestb = crosstime[i+j]
                        p = (j+1)/4
                Psum += crosstime[i+3] - crosstime[i]
            n += 2

        # if last point does not belong to pair (will not have been checked with
        # pair of points one period earlier):
        # method 1: checks radial velocity diffence between last point and pair
        # one period before
        # method 2: measures period using last point and both points from one
        # period before, adding to Psum
        if check:
            n += 2
            Psum += crosstime[r+4] - crosstime[r]
            Psum += crosstime[r+4] - crosstime[r+1]
            if low != None:
                if abs(VP[crossid[r]] - VP[crossid[r+4]]) < low:
                    low = abs(VP[crossid[r]] - VP[crossid[r+4]])
                    besta = crosstime[r]
                    bestb = crosstime[r+4]
                    p = 1
                if abs(VP[crossid[r+1]] - VP[crossid[r+4]]) < low:
                    besta = crosstime[r+1]
                    bestb = crosstime[r+4]
                    p = 1

        if low == None:
            P1 = 1000
        else:
            P1 = (bestb - besta)/p # period from first method

        # find maxima and minima for second method
        maxima = []
        minima = []
        # if less than four crossings, neither two maxima nor two minima exist,
        # so no information about period can be gathered, skip finding
        if check and l < 6:
            P2 = Psum/n
        elif not check and l < 7:
            P2 = Psum/n
        else:
            for i in range(l-1):
                if crossid[i] + 1 != crossid[i+1]:
                    if VP[crossid[i]] > Vavg:
                        maximum = VP[crossid[i]]
                        maximumid = crossid[i]
                        for j in range(crossid[i], crossid[i+1]):
                            if VP[j] > maximum:
                                maximum = VP[j]
                                maximumid = j
                        maxima.append(timeP[maximumid])
                    else:
                        minimum = VP[crossid[i]]
                        minimumid = crossid[i]
                        for j in range(crossid[i], crossid[i+1]):
                            if VP[j] < minimum:
                                minimum = VP[j]
                                minimumid = j
                        minima.append(timeP[minimumid])
            # add time between maxima/minima to Psum
            for i in range(len(maxima)-1):
                Psum += maxima[i+1] - maxima[i]
                n += 1
            for i in range(len(minima)-1):
                Psum += minima[i+1] - minima[i]
                n += 1

            P2 = Psum/n # period from second method

        return [P1, P2]

    def PGuessCandidates():
        """Return plausible periods from chronological single-component data."""
        periodCandidates = []

        # In the equal-coverage double-lined mode, ``time`` and ``V`` are later
        # concatenated as [all primary points, all transformed secondary points].
        # That array is useful for the joint residual, but it is not chronological
        # and must not be passed to the crossing-based period guesser. Estimate
        # periods from each original component separately instead.
        if star == "both":
            periodDataSets = [(time1, V1), (time2, V2)]
        else:
            periodDataSets = [(time, V)]

        for periodTime, periodVelocity in periodDataSets:
            periodTime = np.asarray(periodTime)
            periodVelocity = np.asarray(periodVelocity)
            order = np.argsort(periodTime)
            periodTime = periodTime[order]
            periodVelocity = periodVelocity[order]

            fullGuess = PGuess(
                test = True,
                timeData = periodTime,
                VData = periodVelocity,
            )
            if np.isfinite(np.atleast_1d(fullGuess)).any():
                periodCandidates.extend(PGuess(
                    timeData = periodTime,
                    VData = periodVelocity,
                ))

            gaps = np.diff(periodTime) # looking for big breaks in the data collection
            positiveGaps = gaps[gaps > 0]

            if len(positiveGaps) > 0:
                typicalGap = np.median(positiveGaps)
                largeGaps = np.where(gaps > 150 * typicalGap)[0]

                if len(largeGaps) > 0:
                    boundaries = [0] + list(largeGaps + 1) + [len(periodTime)]
                    clusters = [(boundaries[i], boundaries[i+1]) for i in range(len(boundaries)-1)]
                    # finding the largest data set across the clusters
                    start, stop = max(clusters, key=lambda bounds: bounds[1] - bounds[0])
                    clusterTime = periodTime[start:stop]
                    clusterVelocity = periodVelocity[start:stop]
                    clusterGuess = PGuess(
                        test = True,
                        timeData = clusterTime,
                        VData = clusterVelocity,
                    )

                    if np.isfinite(np.atleast_1d(clusterGuess)).any():
                        periodCandidates.extend(PGuess(
                            timeData = clusterTime,
                            VData = clusterVelocity,
                        ))

        periodCandidates = [P for P in periodCandidates if np.isfinite(P) and P > 0]
        periodCandidates = list(dict.fromkeys(periodCandidates))

        if len(periodCandidates) == 0:
            raise ValueError("Unable to estimate period")

        return periodCandidates

    # Mean anomaly at a given time
    def mean(t, T0, P):
        M = 2*np.pi*(t-T0) / P
        while M < 0:
            M = M + 2*np.pi
        while M > 2*np.pi:
            M = M - 2*np.pi
        return M

    # Kepler's Equation rearranged so RHS = 0
    def Kep(t, T0, E, e, P):
        return mean(t, T0, P) - E + e*np.sin(E)

    # Derivative of kep wrt E
    def dKep(E, e):
        return -1 + e*np.cos(E)

    def ecc(t, T0, e, P):
        """
        Uses Newton's Method to find the eccentric anomaly from Kepler equation
        If t is a list, will return array of corresponding E's
        If t is a single value, will return E at time t
        """
        # if time is a scalar:
        if np.isscalar(t):
            E = mean(t, T0, P) + e*np.sin(mean(t, T0, P)) # starting guess
            for _ in range(100):
                residual = Kep(t, T0, E, e, P)
                if abs(residual) <= 0.000005:
                    break
                derivative_value = dKep(E, e)
                if not np.isfinite(derivative_value) or abs(derivative_value) < 1e-12:
                    raise ConvergenceError("Kepler solver encountered a singular derivative.")
                E = E - residual/derivative_value
            else:
                raise ConvergenceError("Kepler solver failed to converge.")
            while E > 2*np.pi:
                E = E - 2*np.pi
            while E < 0:
                E = E + 2*np.pi

            return E
        # if time is a list:
        else:
            EList = []
            for i in range(len(t)):
                 # starting guess (Green p. 145);
                E = mean(t[i], T0, P) + e*np.sin(mean(t[i], T0, P))
                 # if not within error, improves guess until within error:
                for _ in range(100):
                    residual = Kep(t[i], T0, E, e, P)
                    if abs(residual) <= 0.000005:
                        break
                    derivative_value = dKep(E, e)
                    if not np.isfinite(derivative_value) or abs(derivative_value) < 1e-12:
                        raise ConvergenceError("Kepler solver encountered a singular derivative.")
                    E = E - residual/derivative_value
                else:
                    raise ConvergenceError("Kepler solver failed to converge.")
                # simplifies angle to be between 0 and 2pi
                while E > 2*np.pi or E <= 0:
                    if E > 2*np.pi:
                        E = E - 2*np.pi
                    elif E <= 0:
                        E = E + 2 * np.pi
                EList.append(E)

            return np.array(EList)

    def true(t, T0, e, P, ec = False):
        """
        From E, finds true anomaly using tan(v/2) = sqrt((1+e)/(1-e)) * tan(E/2)
        If t is a list, will return array of corresponding v's
        If t is a single value, will return v at time t
        """
        # if t is a scalar:
        if np.isscalar(t):
            E = ecc(t, T0, e, P)
            v = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
            while v > 2*np.pi:
                v = v - 2*np.pi
            while v < 0:
                v = v + 2*np.pi
            if not ec:
                return v
            else:
                return v, E
        # if t is a list:
        else:
            EList = ecc(t, T0, e, P)
            vList = []
            for i in range(len(EList)):
                v = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(EList[i]/2))
                vList.append(v)
            if not ec:
                return np.array(vList)
            else:
                return np.array(vList), np.array(EList)

    def epoch(w, e, P):
        """
        given w, e, P can return an estimate for T0
        """
        vmax = -w
        Emax = 2 * np.arctan((((1+e)/(1-e))**(-0.5)) * np.tan(vmax/2))
        T0 = time[V.argmax()] - (P/(2*np.pi)) * (Emax - e * np.sin(Emax))
        return T0

    def validInitialGuess(guess):
        """Return True when an initial-parameter guess is finite and usable."""
        if guess is None:
            return False
        try:
            guessArray = np.asarray(guess, dtype = float)
        except (TypeError, ValueError):
            return False
        return guessArray.ndim == 1 and len(guessArray) > 0 and np.all(np.isfinite(guessArray))


    def validTrialResult(result, starMode = "primary"):
        """Reject finite-residual fits that are numerically degenerate."""
        try:
            fittedParameters = result[0]
            if starMode == "both":
                fittedParameters = fittedParameters[0]
            fittedParameters = np.asarray(fittedParameters, dtype = float)
        except (TypeError, ValueError, IndexError):
            return False
        if fittedParameters.ndim != 1 or len(fittedParameters) < 6:
            return False
        if not np.all(np.isfinite(fittedParameters[:6])):
            return False
        fittedK = abs(fittedParameters[1])
        fittedE = fittedParameters[3]
        if fittedE < 0 or fittedE >= 1:
            return False
        velocityScale = np.percentile(V, 95) - np.percentile(V, 5)
        if velocityScale > 0 and fittedK < 1e-4 * velocityScale:
            return False
        return True


    def x0(period = None, star = "primary", gam = 0, zeroEcc = False, fixede = None):
        """
        Makes initial approximation for the orbital parameters γ, K, ω, e, T0, & P.

        X0 can find an initial estimate for period, but a user supplied estimate
        can be inputted using the keyword argument 'period'

        fixede can be used to build an initial guess at a chosen eccentricity.

        Note that, in the following two scenerios, x0() may fail to make a good
        guess:
            i) If the RV data does not cover a full period
            ii) If the RV data is very non-uniformly distributed. E.g. there are
                significantly more points at the peak than at the trough.
        """
        if period == None:
            try:
                periodCandidates = PGuessCandidates()
            except:
                periodCandidates = [100,1000]
        else:
            periodCandidates = [period]

        pOnly = False
        if star == "secondary_both":
            pOnly = True

        candidateFits = []
        for P in periodCandidates:

            K = (max(V) - min(V))/2
            if not np.isfinite(K) or K <= 0:
                continue

            if star == "both":
                V0 = gam
            else:
                V0 = sum(V) / len(V)

            lowSS = None # holds lowest sum of squared deviations
            bestX = None # holds best estimate

            if zeroEcc:

                e,w = 0,0
                T0 = epoch(w, e, P)
                bestX = [V0, K, w, e, T0, P]

            elif fixede == None:
                for i in np.arange(0.01, 0.95, 0.01): # iterates through eccentricities
                    e = i
                    temp = (K*e)**-1 * ((max(V) + min(V))* 0.5 - V0)
                    if abs(temp) <= 1:
                        w1 = np.arccos(temp)
                        x1 = [V0, K, w1, e, epoch(w1,e,P), P]
                        w2 = -w1 % (2 * np.pi) # Non principle value
                        x2 = [V0, K, w2, e, epoch(w2,e,P), P]

                        if SumSquared(x1, None, pOnly = pOnly) < SumSquared(x2, None, pOnly = pOnly):
                            w = w1
                        else:
                            w = w2

                        T0 = epoch(w, e, P)
                        x = [V0, K, w, e, T0, P]
                        currentSS = SumSquared(x, None, pOnly = pOnly)
                        if lowSS == None or currentSS < lowSS:
                            lowSS = currentSS
                            bestX = x.copy()
            else:
                e = float(fixede)
                if 0 < e < 1:
                    temp = (K*e)**-1 * ((max(V) + min(V))* 0.5 - V0)
                    if abs(temp) <= 1:
                        w1 = np.arccos(temp)
                        x1 = [V0, K, w1, e, epoch(w1,e,P), P]
                        w2 = -w1 % (2 * np.pi) # Non principle value
                        x2 = [V0, K, w2, e, epoch(w2,e,P), P]

                        if SumSquared(x1, None, pOnly = pOnly) < SumSquared(x2, None, pOnly = pOnly):
                            bestX = x1.copy()
                        else:
                            bestX = x2.copy()

            if validInitialGuess(bestX):
                candidateFits.append(bestX.copy())

        if len(candidateFits) == 0:
            return None

        bestX = candidateFits[0]
        for x in candidateFits[1:]:
            if SumSquared(x, None, pOnly = pOnly) < SumSquared(bestX, None, pOnly = pOnly):
                bestX = x

        return bestX


    def hasPoorPeriodCoverage(period, starMode = "primary"):
        """Return True when the observations span less than 1.5 trial periods."""
        if period is None or not np.isfinite(period) or period <= 0:
            return False

        if starMode == "secondary_both":
            dataSpan = max(max(time), max(time2)) - min(min(time), min(time2))
        else:
            dataSpan = max(time) - min(time)

        return dataSpan < 1.5 * period


    def initialGuessCandidates(period, starMode = "primary", gam = 0):
        """
        Return the normal x0 guess and, for poorly covered supplied periods,
        additional guesses at several fixed eccentricities.
        """
        guesses = []

        def addGuess(guess):
            if not validInitialGuess(guess):
                return
            for existingGuess in guesses:
                if np.allclose(existingGuess, guess, rtol = 1e-10, atol = 1e-10):
                    return
            guesses.append(list(guess))

        addGuess(x0(period, star = starMode, gam = gam, zeroEcc = zeroEcc))

        if zeroEcc or not hasPoorPeriodCoverage(period, starMode = starMode):
            return guesses

        for fixedEccentricity in [0.05, 0.15, 0.25, 0.30, 0.35, 0.50, 0.70, 0.90]:
            addGuess(x0(period, star = starMode, gam = gam,
                        zeroEcc = zeroEcc, fixede = fixedEccentricity))

        return guesses



    #%% FIT OF DATA

    # rounding function. Makes final output neater.
    def sigFig(x, figs):
        # Preserve the original significant-figure formatting while handling
        # exact zero (for example, the uncertainty of a fixed known period).
        if x == 0:
            return 0
        if not np.isfinite(x):
            return x
        return round(x, -int(math.floor(math.log10(abs(x))) - (figs - 1)))

    def lsq_bisect(x,y,wx,wy):
        """
        Parameters
        ----------
        x : x points
        y : y points
        wx : x weights
        wy : y weights

        Returns
        -------
        cov: coefficient vector of the bisector fit line
        """

        cof1= np.polyfit( x, y, 1, w = wy)
        cof2= np.polyfit( y, x, 1, w = wx)

        b1= cof1[0] #slope of normal fit:  y= a1 + b1*x
        b2= 1.0/cof2[0] #slope of inverted fit: x= a2 + b2*y --> y= -a2/b2 + x/b2

        # now compute the slope of the bisector of these two OLS lines

        t1= np.sqrt( 1.0 + b1**2 )
        t2= np.sqrt( 1.0 + b2**2 )

        b = ( b2*t1 + b1*t2 )/( t1 + t2 )     #slope of bisector line

         # the regression line always passes through the centroid (mean(x),mean(y)) and
         # so we use this to find the constant term: a = mean(y) - bs*mean(x)

        xm= np.mean(x)
        ym= np.mean(y)
        a = ym - b*xm
        cof= np.array([b, a])  #define the coefficient vector of the bisector fit line

        return cof

    def derivative(f, params, P, ang, der, h=0.000001, zeroEcc = False, starNum = None, ang2 = None):
        """
        Numerical derivative using central differences

        Keyword Arguments:
        f - The function to differentiate.
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        der - which variable to take derivative w.r.t. The options are "V0", "K",
              "w", "e", and "T0"/
        h - small change variable for taking difference
        starNum - If its for a single star during the star="both" case
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0

        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0

        if star != "secondary_both":

            if der == "V0":
                paramsAb = [V0 + h, K, w, e, T0]
                paramsBe = [V0 - h, K, w, e, T0]
                angAb = ang
                angBe = ang
            elif der == "K":
                paramsAb = [V0, K + h, w, e, T0]
                paramsBe = [V0, K - h, w, e, T0]
                angAb = ang
                angBe = ang
            elif der == "w":
                paramsAb = [V0, K, w + h, e, T0]
                paramsBe = [V0, K, w - h, e, T0]
                angAb = ang
                angBe = ang
            elif der == "e":
                paramsAb = [V0, K, w, e + h, T0]
                paramsBe = [V0, K, w, e - h, T0]
                angAb = true(time, T0, e + h, P, ec = True)
                angBe = true(time, T0, e - h, P, ec = True)
            elif der == "T0":
                paramsAb = [V0, K, w, e, T0 + h]
                paramsBe = [V0, K, w, e, T0 - h]
                angAb = true(time, T0 + h, e, P, ec = True)
                angBe = true(time, T0 - h, e, P, ec = True)
            else:
                params = [V0, K, w, e, T0]
                angAb = true(time, T0, e, P + h, ec = True)
                angBe = true(time, T0, e, P - h, ec = True)
                # Preserve the original package's period-derivative behavior for
                # legacy primary/secondary/same-epoch modes. The unequal-coverage
                # mode below uses the corrected P-h denominator.
                return (f(params, angAb, P+h,starNum) - f(params,angBe, P+h,starNum))/(2*h)

            return (f(paramsAb, angAb, P,starNum) - f(paramsBe, angBe, P, starNum))/(2*h)

        else: # both_secondary
            if der == "V0":
                paramsAb = [V0 + h, K, w, e, T0, V02,K2]
                paramsBe = [V0 - h, K, w, e, T0,V02,K2]
                angAb = ang
                angBe = ang
                angAb2 = ang2
                angBe2 = ang2
            elif der == "K":
                paramsAb = [V0, K+h, w, e, T0, V02,K2]
                paramsBe = [V0, K-h, w, e, T0,V02,K2]
                angAb = ang
                angBe = ang
                angAb2 = ang2
                angBe2 = ang2
            elif der == "w":
                paramsAb = [V0, K, w+h, e, T0, V02,K2]
                paramsBe = [V0, K, w-h, e, T0,V02,K2]
                angAb = ang
                angBe = ang
                angAb2 = ang2
                angBe2 = ang2
            elif der == "e":
                paramsAb = [V0, K, w, e+h, T0, V02,K2]
                paramsBe = [V0, K, w, e-h, T0,V02,K2]
                angAb = true(time, T0, e + h, P, ec = True)
                angBe = true(time, T0, e - h, P, ec = True)
                angAb2 = true(time2, T0, e + h, P, ec = True)
                angBe2 = true(time2, T0, e - h, P, ec = True)

            elif der == "T0":
                paramsAb = [V0, K, w, e, T0+h, V02,K2]
                paramsBe = [V0, K, w, e, T0-h,V02,K2]
                angAb = true(time, T0 + h, e, P, ec = True)
                angBe = true(time, T0 - h, e, P, ec = True)
                angAb2 = true(time2, T0 + h, e, P, ec = True)
                angBe2 = true(time2, T0 - h, e, P, ec = True)

            elif der == "P":
                params = [V0, K, w, e, T0, V02,K2]
                angAb = true(time, T0, e, P + h, ec = True)
                angBe = true(time, T0, e, P - h, ec = True)
                angAb2 = true(time2, T0, e, P+h, ec = True)
                angBe2 = true(time2, T0, e, P-h, ec = True)
                return (f(params, angAb, P+h,starNum, ang2 = angAb2) - f(params,angBe, P-h,starNum, ang2 = angBe2))/(2*h)

            elif der == "V02":
                paramsAb = [V0, K, w, e, T0, V02+h,K2]
                paramsBe = [V0, K, w, e, T0,V02-h,K2]
                angAb = ang
                angBe = ang
                angAb2 = ang2
                angBe2 = ang2
            elif der == "K2":
                paramsAb = [V0, K, w, e, T0, V02,K2+h]
                paramsBe = [V0, K, w, e, T0,V02,K2-h]
                angAb = ang
                angBe = ang
                angAb2 = ang2
                angBe2 = ang2

            return (f(paramsAb, angAb, P,starNum, ang2 = angAb2) - f(paramsBe, angBe, P, starNum, ang2 = angBe2))/(2*h)

    # Sum of squared deviations between the fit and the actual RV data
    def SumSquared(params, P, starNum = None, pOnly = False):

        if star != "secondary_both" or pOnly:
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0

        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0

        if star == "both":
            if starNum == 1:
                t = time1
                epsilon = ((V0 + ((np.cos(true(t, T0, e, P) + w) + e*np.cos(w))*K)) - V1)**2
                SS = sum(epsilon*weight1)

            elif starNum == 2:
                t = time2
                epsilon = ((V0 + ((np.cos(true(t, T0, e, P) + w) + e*np.cos(w))*K)) - V2)**2
                SS = sum(epsilon*weight2)
            else:
                t = timecomb
                epsilon = ((V0 + ((np.cos(true(t, T0, e, P) + w) + e*np.cos(w))*K)) - Vcomb)**2
                SS = sum(epsilon*weightcomb)
        elif star != "secondary_both" or pOnly:
            epsilon = ((V0 + ((np.cos(true(time, T0, e, P) + w) + e*np.cos(w))*K)) - V)**2
            SS = sum(epsilon* weight)
        else:
            w2 = (w + np.pi)%(2*np.pi)
            epsilon1 = ((V0 + ((np.cos(true(time, T0, e, P) + w) + e*np.cos(w))*K)) - V)**2
            epsilon2 = ((V02 + ((np.cos(true(time2, T0, e, P) + w2) + e*np.cos(w2))*K2)) - V2)**2
            SS = sum(epsilon1*weight) + sum(epsilon2*weight2)

        return SS

    def partV0(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt γ.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else: # not nessicary
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2 = ang2[0]
            w2 = (w + np.pi)%(2*np.pi)

        v = ang[0]

        if starNum == 1:
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V1) *weight1
            temp = sum(2*temp)

        elif starNum == 2:
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V2) *weight2
            temp = sum(2*temp)

        elif star == "both":
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - Vcomb) *weightcomb
            temp = sum(2*temp)

        elif star != "secondary_both":
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V) *weight
            temp = sum(2*temp)

        else:
            temp1 = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V) *weight
            temp = sum(2*temp1)

        return temp

    def partK(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt K.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else: # not nessicary
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2 = ang2[0]
            w2 = (w + np.pi)%(2*np.pi)
        v = ang[0]

        if starNum == 1:
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V1) *weight1

        elif starNum == 2:
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V2) *weight2

        elif star == "both":
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - Vcomb) *weightcomb

        else:
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V) *weight

        temp *= (np.cos(v+w) + e*np.cos(w))

        return sum(2*temp)

    def partw(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt w.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2 = ang2[0]
            w2 = (w + np.pi)%(2*np.pi)
        v = ang[0]

        if star != "secondary_both":
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V) *weight
            temp *= (np.sin(v+w)+e*np.sin(w))
            return sum((-2*K*temp))

        else:
            temp1 = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V) *weight
            temp1 *= (np.sin(v+w)+e*np.sin(w))
            temp2 = ((V02 + ((np.cos(v2 + w2) + e*np.cos(w2))*K2)) - V2) *weight2
            temp2 *= (np.sin(v2+w2)+e*np.sin(w2))
            return sum((-2*K*temp1)) + sum((-2*K2*temp2))

    def parte(params, ang, P = None, starNum = None,ang2 = None):
        """
        Derivative of residual wrt e.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies and eccentric anomalies given parameters
              [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2, E2 = ang2
            w2 = (w + np.pi)%(2*np.pi)
        v, E = ang

        if star != "secondary_both":
            dEde = np.sin(E)/(1-e*np.cos(E))
            dvde = np.sqrt((e+1)/(1-e)) * (0.5) * (np.cos(E/2)**-2)
            dvde *= dEde
            dvde += (1-e)**-2 * np.sqrt((1-e)/(1+e)) * np.tan(E/2)
            dvde *= 2 * np.cos(v/2)**2
            temp = K*(-np.sin(v + w)*dvde + np.cos(w))
            temp *= (V0 + (np.cos(v + w) + e*np.cos(w))*K - V) *weight
            return sum(2*temp)
        else:
            dEde = np.sin(E)/(1-e*np.cos(E))
            dvde = np.sqrt((e+1)/(1-e)) * (0.5) * (np.cos(E/2)**-2)
            dvde *= dEde
            dvde += (1-e)**-2 * np.sqrt((1-e)/(1+e)) * np.tan(E/2)
            dvde *= 2 * np.cos(v/2)**2
            temp1 = K*(-np.sin(v + w)*dvde + np.cos(w))
            temp1 *= (V0 + (np.cos(v + w) + e*np.cos(w))*K - V) *weight

            dEde = np.sin(E2)/(1-e*np.cos(E2))
            dvde = np.sqrt((e+1)/(1-e)) * (0.5) * (np.cos(E2/2)**-2)
            dvde *= dEde
            dvde += (1-e)**-2 * np.sqrt((1-e)/(1+e)) * np.tan(E2/2)
            dvde *= 2 * np.cos(v2/2)**2
            temp2 = K2*(-np.sin(v2 + w2)*dvde + np.cos(w2))
            temp2 *= (V02 + (np.cos(v2 + w2) + e*np.cos(w2))*K2 - V2) *weight2

            return sum(2*temp1) +  sum(2*temp2)

    def partT0(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt T0.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies and eccentric anomalies given parameters
              [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2, E2 = ang2
            w2 = (w + np.pi)%(2*np.pi)
        v, E = ang

        if star != "secondary_both":
            temp = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight
            temp *= K * np.sqrt((1+e)/(1-e)) * np.sin(v + w)
            temp *= (1+ np.cos(v)) / (1 + np.cos(E))
            temp *= 2*np.pi / (P*(1-e*np.cos(E)))
            return sum(2*temp)
        else:
            temp1 = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight
            temp1 *= K * np.sqrt((1+e)/(1-e)) * np.sin(v + w)
            temp1 *= (1+ np.cos(v)) / (1 + np.cos(E))
            temp1 *= 2*np.pi / (P*(1-e*np.cos(E)))

            temp2 = ((V02 + ((np.cos(v2 + w2) + e*np.cos(w2))*K2)) - V2)*weight2
            temp2 *= K2 * np.sqrt((1+e)/(1-e)) * np.sin(v2 + w2)
            temp2 *= (1+ np.cos(v2)) / (1 + np.cos(E2))
            temp2 *= 2*np.pi / (P*(1-e*np.cos(E2)))

            return sum(2*temp1) + sum(2*temp2)

    def partP(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt P.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies and eccentric anomalies given parameters
              [γ, K, w, e, T0, P]
        """

        if star != "secondary_both":
            if len(params) == 6:
                V0, K, w, e, T0, P = params
            elif len(params) == 5:
                V0, K, w, e, T0 = params
            elif len(params) == 4:
                V0, K, T0,P = params
                e,w = 0,0
            else:
                V0, K, T0 = params
                e,w = 0,0
        else:
            if len(params) == 8:
                V0, K, w, e, T0, P, V02, K2 = params
            elif len(params) == 7:
                V0, K, w, e, T0, V02, K2 = params
            elif len(params) == 6:
                V0, K, T0,P,V02, K2 = params
                e,w = 0,0
            else:
                V0, K, T0, V02, K2 = params
                e,w = 0,0
            v2, E2 = ang2
            w2 = (w + np.pi)%(2*np.pi)
        v, E = ang
        warnings.filterwarnings("ignore")

        if star != "secondary_both":
            dEdP = ((-2*np.pi)/(P**2)) * (time - T0) * (1-e*np.cos(E))**(-1)
            dvdP = np.sqrt((1+e)/(1-e)) * ((1+ np.cos(v)) / (1 + np.cos(E))) * dEdP
            temp = -K * np.sin(v + w) * dvdP
            temp *= ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight
            return sum(2*temp)

        else:
            dEdP = ((-2*np.pi)/(P**2)) * (time - T0) * (1-e*np.cos(E))**(-1)
            dvdP = np.sqrt((1+e)/(1-e)) * ((1+ np.cos(v)) / (1 + np.cos(E))) * dEdP
            temp1 = -K * np.sin(v + w) * dvdP
            temp1 *= ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight

            dEdP = ((-2*np.pi)/(P**2)) * (time2 - T0) * (1-e*np.cos(E2))**(-1)
            dvdP = np.sqrt((1+e)/(1-e)) * ((1+ np.cos(v2)) / (1 + np.cos(E2))) * dEdP
            temp2 = -K2 * np.sin(v2 + w2) * dvdP
            temp2 *= ((V02 + ((np.cos(v2 + w2) + e*np.cos(w2))*K2)) - V2)*weight2

            return sum(2*temp1) + sum(2*temp2)

    def partV02(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt γ.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """

        if len(params) == 8:
            V0, K, w, e, T0, P, V02, K2 = params
        elif len(params) == 7:
            V0, K, w, e, T0, V02, K2 = params
        elif len(params) == 6:
            V0, K, T0,P,V02, K2 = params
            e,w = 0,0
        else:
            V0, K, T0, V02, K2 = params
            e,w = 0,0
        v2 = ang2[0]
        w2 = (w + np.pi)%(2*np.pi)

        v = ang[0]
        v2 = ang2[0]

        temp1 = ((V02 + ((np.cos(v2 + w2) + e*np.cos(w2)))*K2) - V2) *weight2
        temp = sum(2*temp1)

        return temp

    def partK2(params, ang, P = None, starNum = None, ang2 = None):
        """
        Derivative of residual wrt K.

        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """

        if len(params) == 8:
            V0, K, w, e, T0, P, V02, K2 = params
        elif len(params) == 7:
            V0, K, w, e, T0, V02, K2 = params
        elif len(params) == 6:
            V0, K, T0,P,V02, K2 = params
            e,w = 0,0
        else:
            V0, K, T0, V02, K2 = params
            e,w = 0,0
        v2 = ang2[0]
        w2 = (w + np.pi)%(2*np.pi)
        temp = ((V02 + ((np.cos(v2 + w2) + e*np.cos(w2)))*K2) - V2) *weight2
        temp *= (np.cos(v2+w2) + e*np.cos(w2))

        return sum(2*temp)

    def PrimarySolve(period = None, Pguess = None,covariance = False, graphs = True,correction = False, webgrapher = False, initialX = None, multiStart = True, returnScore = False, trialOnly = False):
        """
        Solves for orbital parameters from radial velocity data

        Keyword Arguments:
        period (optional) - Known value for period (if none provided, by default
                      PrimarySolve makes its own guess)
        covariance (optional) - If covariance = True, PrimarySolve will return the
                                entire covariance matrix.
        graphs (optional) - If graphs = true, PrimarySolve creates 2 plots. First, a
                            plot of the radial velocity data as a function of time
                            (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                            where [γ, K, ω, e, T0] are the orbital parameters
                            solved for by PrimarySolve(). Second, a plot of the
                            radial velocity data as a function of phase (between 0
                            and 1) of the star's orbit. Again, fit with
                            V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                            function of phase. On both plots, the long-term average
                            velocity, γ, is also shown.

        Returns:
        x - List of solved orbital parameters, the semi major axis and the
            mass function [γ, K, ω, e, T0, P, a, f(M)]
        err - List of standard errors on the determination of the orbital elements,
              in the order [γ, K, ω, e, T0, P, a, f(M)]
        C - Estimated covariance matrix (if argument covariance = True)

        """
        requestedPeriod = Pguess if Pguess != None else period
        poorCoverage = hasPoorPeriodCoverage(requestedPeriod)
        if initialX is None and multiStart and poorCoverage:
            trialGuesses = initialGuessCandidates(requestedPeriod)
            bestTrialResult = None
            bestSS = None
            for trialGuess in trialGuesses:
                try:
                    trialResult, trialSS = PrimarySolve(
                        period = period, Pguess = Pguess, covariance = covariance,
                        graphs = graphs, correction = correction,
                        webgrapher = webgrapher,
                        initialX = trialGuess, multiStart = False,
                        returnScore = True)
                except (ConvergenceError, np.linalg.LinAlgError, FloatingPointError,
                        ValueError, TypeError, ZeroDivisionError):
                    continue
                if validTrialResult(trialResult) and np.isfinite(trialSS) and (bestSS is None or trialSS < bestSS):
                    bestTrialResult = trialResult
                    bestSS = trialSS
            if bestTrialResult is not None:
                if returnScore:
                    return bestTrialResult, bestSS
                return bestTrialResult
            raise ConvergenceError("All initial guesses failed to converge.")

        if Pguess != None:
            period = Pguess
        if initialX is None:
            x = x0(period, zeroEcc = zeroEcc) # initial guess
        else:
            x = list(initialX)
        if not validInitialGuess(x):
            raise ConvergenceError("Unable to construct a valid initial guess.")

        if Pguess != None:
            period = None

        if period != None: # if period known and provided by user
            P = period
            x.pop()
            V0, K, w, e, T0 = x
            var = 5
        else:
            V0, K, w, e, T0, P = x
            var = 6

        if zeroEcc:
            if period != None:
                x = [V0, K, T0]
                var = 3
            else:
                x = [V0, K,T0,P]
                var = 4
        l = 3
        ct = 0 # convergence timer
        while True: # perfoming minimization
            if e < 0:
                e = 0  # so some square roots dont cause us trouble
                x[3] = e
            if e >=1:
                e = 0.99
                x[3]= e

            xLast = x
            ang = true(time, T0, e, P, ec = True)

            pVK = derivative(partK, x, P, ang, 'V0')
            pVT0 = derivative(partT0, x, P, ang, 'V0')
            pKT0 = derivative(partT0, x, P, ang, 'K')

            if period == None:

                # Building Hessian matrix

                pVP = derivative(partV0, x, P, ang, 'P')
                pKP = derivative(partK, x, P, ang, 'P')
                pT0P = derivative(partT0, x, P, ang, 'P')

                if not zeroEcc:
                    peP = derivative(parte, x, P, ang, 'P')
                    pVe = derivative(parte, x, P, ang, 'V0')
                    pKe = derivative(parte, x, P, ang, 'K')
                    pwe = derivative(parte, x, P, ang, 'w')
                    peT0 = derivative(partT0, x, P, ang, 'e')
                    pVw = derivative(partw, x, P, ang, 'V0')
                    pKw = derivative(partw, x, P, ang, 'K')
                    pwT0 = derivative(partT0, x, P, ang, 'w')
                    pwP = derivative(partw, x, P, ang, 'P')

                    H1 = [derivative(partV0, x, P, ang, 'V0') * (1+l), pVK, pVw, pVe, pVT0, pVP]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKw, pKe, pKT0, pKP]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w') * (1+l), pwe, pwT0, pwP]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e') * (l+1), peT0, peP]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0') * (1+l), pT0P]
                    H6 = [pVP, pKP, pwP, peP, pT0P, derivative(partP, x, P, ang, 'P') * (1+l)]
                    a = 0.5 * np.array([H1, H2, H3, H4, H5, H6])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P), partw(x, ang, P),
                                         parte(x, ang, P), partT0(x, ang, P), partP(x, ang, P)])

                else:
                    H1 = [derivative(partV0, x, P, ang, 'V0') * (1+l), pVK, pVT0, pVP]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKT0, pKP]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0') * (1+l), pT0P]
                    H4 = [pVP, pKP, pT0P, derivative(partP, x, P, ang, 'P') * (1+l)]
                    a = 0.5 * np.array([H1, H2, H3, H4])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P),
                                         partT0(x, ang, P), partP(x, ang, P)])

            else:

                if not zeroEcc:

                    pVe = derivative(parte, x, P, ang, 'V0')
                    pKe = derivative(parte, x, P, ang, 'K')
                    pwe = derivative(parte, x, P, ang, 'w')
                    peT0 = derivative(partT0, x, P, ang, 'e')
                    pVw = derivative(partw, x, P, ang, 'V0')
                    pKw = derivative(partw, x, P, ang, 'K')
                    pwT0 = derivative(partT0, x, P, ang, 'w')

                    H1 = [derivative(partV0, x, P, ang, "V0") * (1+l), pVK, pVw, pVe, pVT0]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKw, pKe, pKT0]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w') * (1+l), pwe, pwT0]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e') * (l+1), peT0]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0') * (1+l)]
                    a = 0.5 * np.array([H1,H2,H3,H4,H5])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P), partw(x, ang, P),
                                         parte(x, ang, P), partT0(x, ang, P)])

                else:
                    H1 = [derivative(partV0, x, P, ang, "V0") * (1+l), pVK, pVT0]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKT0]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0') * (1+l)]
                    a = 0.5 * np.array([H1,H2,H3])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P),
                                         partT0(x, ang, P)])
            x = xLast + np.linalg.solve(a, b) # Next iterations parameter values

            if period == None:
                if zeroEcc:
                    V0, K, T0, P = x
                else:
                    V0, K, w, e, T0, P = x
            else:
                if zeroEcc:
                    V0, K, T0 = x
                else:
                    V0, K, w, e, T0 = x
            ct +=1
            if ct >= 1100 :
                raise ConvergenceError
            if SumSquared(x, P) < SumSquared(xLast, P): # if params got better
                if abs(SumSquared(xLast, P) - SumSquared(x, P)) < 0.01 or abs(1-(SumSquared(x, P)/SumSquared(xLast, P))) < 1e-3:
                    break # convergence acheived! (probably)
                else:
                    l /= 9 # lowering damping parameter
                    l = max(l, 10e-7)
            else: # if params got worse
                l *= 11 # raising damping parameter
                l = min(l, 10e7)

        if not zeroEcc:
            if x[2] < -np.pi or x[2] >= np.pi:
                while x[2] < -np.pi:
                    x[2] = x[2] + 2*np.pi
                while x[2] >= np.pi:
                    x[2] = x[2] - 2 * np.pi
                w = x[2]

        if period == None:
            if zeroEcc:
                P = x[3]
            else:
                P = x[5]
        x = [V0,K,w,e,T0,P]
        finalSS = SumSquared(x, None)
        if trialOnly:
            return x, finalSS

        # other values to be returned
        nsec = 2*np.pi / (24*3600*P)
        a1 = (x[1]*np.sqrt(1-x[3]**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2
        x.append(a1)
        x.append(fm)

        C = np.linalg.inv(a) # covariance matrix
        err = []

        for q in range(var):
            err.append(np.sqrt(C[q][q]))
        if period != None:
            err.append(0)
        if zeroEcc:
            err.insert(2,0)
            err.insert(3,0)

        # Finding uncertainty in a1
        dade = (-x[1] * x[3] * x[5] * (3600*24)) / (2 * np.pi * np.sqrt(1-x[3]**2))
        dadT = (x[1] * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        dadK = (x[5] * (3600*24) * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        aErr = np.sqrt((dadT * err[5] * (3600*24))**2 + (dadK * err[1])**2 + (dade * err[3])**2)
        aErr *= (10e3 / 10e9)

        # Finding uncertainty in fm
        dfda = 3 * (x[6] * (10e9 / 10e3) )**2 * (3.985e-20) * (x[5])**-2
        dfdT = -2 * (x[6] * (10e9 / 10e3) )**3 * (3.985e-20) * (x[5])**-3
        fErr = np.sqrt((dfda * aErr * (10e9 / 10e3) )**2 + (dfdT * err[5])**2)

        err.append(aErr)
        err.append(fErr)

        # create plots if expected

        # rounding

        webplot = np.empty((1000,4))
        if webgrapher and not correction:

            times = np.linspace(min(time), max(time),1000) # times for fit
            v = true(times, x[4], x[3], P) # array of true anomalies at every time for fit
            fit = x[0] + x[1]*(np.cos(v + x[2])+ x[3]*np.cos(x[2]))

            # RV vs phase data
            phase = ((time - x[4]) / P) % 1 # phase for every RV data point
            # phases = np.linspace(min(phase), max(phase), N)
            phases = ((times - x[4]) / P) % 1 # phase for every point on fit
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            for i in range(len(phases)):
                if i < len(V):
                    phaseV.append([phase[i], V[i]])
                phasesfit.append([phases[i], fit[i]])
            phaseV = np.array(sorted(phaseV,key=lambda l:l[0]))
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))

            webplot[:,0] = times
            webplot[:,1] = fit
            webplot[:,2] = phasesfit[:,0]
            webplot[:,3] = phasesfit[:,1]

        V0 = sigFig(x[0], 5)
        K = sigFig(x[1], 5)
        if not zeroEcc:
            w = sigFig((x[2] * 180/np.pi)%360, 6)
            if e < 0.1:
                e = sigFig(x[3], 4)
            else:
                e = sigFig(x[3], 4)
        T0 = sigFig(x[4], 6)
        a1 = sigFig(x[6], 5)
        fm = sigFig(x[7], 4)

        V0err = sigFig(err[0], 5)
        Kerr = sigFig(err[1], 5)
        werr = err[2]
        eerr = err[3]
        if not zeroEcc:
              werr = sigFig((err[2] * 180/np.pi)%360, 6)
              eerr = sigFig(err[3], 5)
        T0err = sigFig(err[4], 6)
        a1err = sigFig(err[6], 5)
        fmerr = sigFig(err[7], 5)

        if period == None:
            P = sigFig(x[5], 6)
            Perr = sigFig(err[5], 6)
        else:
            P = x[5]
            Perr = err[5]

        x = [V0, K, w, e, T0, P, a1, fm]
        err = [V0err, Kerr, werr, eerr, T0err, Perr, a1err, fmerr]

        # return expected values
        if covariance:
            for i in range(len(C)):
                for j in range(len(C)):
                    C[i][j] = sigFig(C[i][j], 5)
            if webgrapher:
                result = (x,err,C,webplot)
            else:
                result = (x, err, C)
        else:
            result = (x, err)

        if returnScore:
            return result, finalSS
        return result

    def CompanionSolve_both(period = None, Pguess = None,covariance = False, graphs = False, correction = False, webgrapher = False, gam2set = True, initialX = None, multiStart = True, returnScore = False, trialOnly = False):
        """
        Solves for orbital parameters from radial velocity data

        Keyword Arguments:
        period (optional) - Known value for period (if none provided, by default
                      PrimarySolve makes its own guess)
        covariance (optional) - If covariance = True, PrimarySolve will return the
                                entire covariance matrix.
        graphs (optional) - If graphs = true, PrimarySolve creates 2 plots. First, a
                            plot of the radial velocity data as a function of time
                            (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                            where [γ, K, ω, e, T0] are the orbital parameters
                            solved for by PrimarySolve(). Second, a plot of the
                            radial velocity data as a function of phase (between 0
                            and 1) of the star's orbit. Again, fit with
                            V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                            function of phase. On both plots, the long-term average
                            velocity, γ, is also shown.

        Returns:
        x - List of solved orbital parameters, the semi major axis and the
            mass function [γ, K, ω, e, T0, P, a, f(M)]
        err - List of standard errors on the determination of the orbital elements,
              in the order [γ, K, ω, e, T0, P, a, f(M)]
        C - Estimated covariance matrix (if argument covariance = True)

        """

        requestedPeriod = Pguess if Pguess != None else period
        poorCoverage = hasPoorPeriodCoverage(requestedPeriod, starMode = "secondary_both")
        if initialX is None and multiStart and poorCoverage:
            trialGuesses = initialGuessCandidates(requestedPeriod, starMode = "secondary_both")
            bestTrialResult = None
            bestSS = None
            for trialGuess in trialGuesses:
                try:
                    trialResult, trialSS = CompanionSolve_both(
                        period = period, Pguess = Pguess, covariance = covariance,
                        graphs = graphs, correction = correction,
                        webgrapher = webgrapher,
                        gam2set = gam2set,
                        initialX = trialGuess, multiStart = False,
                        returnScore = True)
                except (ConvergenceError, np.linalg.LinAlgError, FloatingPointError,
                        ValueError, TypeError, ZeroDivisionError):
                    continue
                if validTrialResult(trialResult, starMode = "secondary_both") and np.isfinite(trialSS) and (bestSS is None or trialSS < bestSS):
                    bestTrialResult = trialResult
                    bestSS = trialSS
            if bestTrialResult is not None:
                if returnScore:
                    return bestTrialResult, bestSS
                return bestTrialResult
            raise ConvergenceError("All initial guesses failed to converge.")

        if Pguess != None:
            period = Pguess
        if initialX is None:
            x = x0(period, zeroEcc = zeroEcc, star = "secondary_both") # initial guess
        else:
            x = list(initialX)
        if not validInitialGuess(x):
            raise ConvergenceError("Unable to construct a valid initial guess.")

        if Pguess != None:
            period = None

        if period != None: # if period known and provided by user
            P = period
            x.pop()
            V0, K, w, e, T0 = x
            K2 = K
            V02 = V0
            x = [V0,K,w,e,T0,V02,K2]
            var = 7
        else:
            V0, K, w, e, T0, P = x
            K2 = K
            V02 = V0
            x = [V0,K,w,e,T0,P,V02,K2]
            var = 8

        if zeroEcc:
            if period != None:
                x = [V0, K, T0, V02,K2]
                var = 5
            else:
                x = [V0, K,T0,P, V02,K2]
                var = 6
        l = 3
        ct = 0
        while True: # perfoming minimization

            if e < 0:
                e = 0  # so some square roots dont cause us trouble
                x[3] = e
            if e >=1:
                e = 0.99
                x[3]= e

            xLast = x
            ang = true(time, T0, e, P, ec = True)
            ang2 = true(time2, T0, e, P, ec = True)

            pVK = derivative(partK, x, P, ang, 'V0', ang2 = ang2)
            pVT0 = derivative(partT0, x, P, ang, 'V0', ang2 = ang2)
            pVV2 = derivative(partV02, x, P, ang, 'V0', ang2 = ang2)
            pVK2 = derivative(partK2, x, P, ang, 'V0', ang2 = ang2)
            pKT0 = derivative(partT0, x, P, ang, 'K', ang2 = ang2)
            pKV2 = derivative(partV02, x, P, ang, 'K', ang2 = ang2)
            pKK2 = derivative(partK2, x, P, ang, 'K', ang2 = ang2)
            pT0K2 = derivative(partK2, x, P, ang, 'T0', ang2 = ang2)
            pT0V2 = derivative(partV02, x, P, ang, 'T0', ang2 = ang2)
            pV2K2 = derivative(partV02, x, P, ang, 'K2', ang2 = ang2)

            if period == None:

                # Building Hessian matrix

                pVP = derivative(partV0, x, P, ang, 'P', ang2 = ang2)
                pKP = derivative(partK, x, P, ang, 'P', ang2 = ang2)
                pT0P = derivative(partT0, x, P, ang, 'P', ang2 = ang2)
                pPV2 = derivative(partV02, x, P, ang, 'P', ang2 = ang2)
                pPK2 = derivative(partK2, x, P, ang, 'P', ang2 = ang2)

                if not zeroEcc:
                    peP = derivative(parte, x, P, ang, 'P', ang2 = ang2)
                    pVe = derivative(parte, x, P, ang, 'V0', ang2 = ang2)
                    pKe = derivative(parte, x, P, ang, 'K', ang2 = ang2)
                    pwe = derivative(parte, x, P, ang, 'w', ang2 = ang2)
                    peT0 = derivative(partT0, x, P, ang, 'e', ang2 = ang2)
                    peV2 = derivative(partV02, x, P, ang, 'e', ang2 = ang2)
                    peK2 = derivative(partK2, x, P, ang, 'e', ang2 = ang2)

                    pVw = derivative(partw, x, P, ang, 'V0', ang2 = ang2)
                    pKw = derivative(partw, x, P, ang, 'K', ang2 = ang2)
                    pwT0 = derivative(partT0, x, P, ang, 'w', ang2 = ang2)
                    pwP = derivative(partw, x, P, ang, 'P', ang2 = ang2)
                    pwV2 = derivative(partV02, x, P, ang, 'w', ang2 = ang2)
                    pwK2 = derivative(partK2, x, P, ang, 'w', ang2 = ang2)

                    H1 = [derivative(partV0, x, P, ang, 'V0', ang2 = ang2) * (1+l), pVK, pVw, pVe, pVT0, pVP, pVV2, pVK2]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K', ang2 = ang2) * (1+l), pKw, pKe, pKT0, pKP, pKV2, pKK2]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w', ang2 = ang2) * (1+l), pwe, pwT0, pwP, pwV2, pwK2]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e', ang2 = ang2) * (l+1), peT0, peP, peV2, peK2]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0', ang2 = ang2) * (1+l), pT0P, pT0V2, pT0K2]
                    H6 = [pVP, pKP, pwP, peP, pT0P, derivative(partP, x, P, ang, 'P', ang2 = ang2) * (1+l), pPV2, pPK2]
                    H7 = [pVV2, pKV2, pwV2, peV2, pT0V2, pPV2, derivative(partV02, x, P, ang, 'V02', ang2 = ang2) * (1+l), pV2K2]
                    H8 = [pVK2, pKK2, pwK2, peK2, pT0K2, pPK2, pV2K2, derivative(partK2, x, P, ang, 'K2', ang2 = ang2) * (1+l)]

                    a = 0.5 * np.array([H1, H2, H3, H4, H5, H6,H7,H8])
                    b = -0.5 * np.array([partV0(x, ang, P, ang2 = ang2), partK(x, ang, P, ang2 = ang2), partw(x, ang, P, ang2 = ang2),
                                         parte(x, ang, P, ang2 = ang2), partT0(x, ang, P, ang2 = ang2), partP(x, ang, P, ang2 = ang2),
                                         partV02(x, ang, P, ang2 = ang2), partK2(x, ang, P, ang2 = ang2)])

                else:
                    H1 = [derivative(partV0, x, P, ang, 'V0', ang2 = ang2) * (1+l), pVK, pVT0, pVP, pVV2, pVK2]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K', ang2 = ang2) * (1+l), pKT0, pKP, pKV2, pKK2]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0', ang2 = ang2) * (1+l), pT0P, pT0V2, pT0K2]
                    H4 = [pVP, pKP, pT0P, derivative(partP, x, P, ang, 'P', ang2 = ang2) * (1+l), pPV2, pPK2]
                    H5 = [pVV2, pKV2, pT0V2, pPV2, derivative(partV02, x, P, ang, 'V02', ang2 = ang2) * (1+l), pV2K2]
                    H6 = [pVK2, pKK2, pT0K2, pPK2, pV2K2, derivative(partK2, x, P, ang, 'K2', ang2 = ang2) * (1+l)]
                    a = 0.5 * np.array([H1, H2, H3, H4,H5,H6])
                    b = -0.5 * np.array([partV0(x, ang, P, ang2 = ang2), partK(x, ang, P, ang2 = ang2),
                                          partT0(x, ang, P, ang2 = ang2), partP(x, ang, P, ang2 = ang2),
                                         partV02(x, ang, P, ang2 = ang2), partK2(x, ang, P, ang2 = ang2)])

            else:

                if not zeroEcc:

                    pVe = derivative(parte, x, P, ang, 'V0', ang2 = ang2)
                    pKe = derivative(parte, x, P, ang, 'K', ang2 = ang2)
                    pwe = derivative(parte, x, P, ang, 'w', ang2 = ang2)
                    peT0 = derivative(partT0, x, P, ang, 'e', ang2 = ang2)
                    peV2 = derivative(partV02, x, P, ang, 'e', ang2 = ang2)
                    peK2 = derivative(partK2, x, P, ang, 'e', ang2 = ang2)

                    pVw = derivative(partw, x, P, ang, 'V0', ang2 = ang2)
                    pKw = derivative(partw, x, P, ang, 'K', ang2 = ang2)
                    pwT0 = derivative(partT0, x, P, ang, 'w', ang2 = ang2)
                    pwV2 = derivative(partV02, x, P, ang, 'w', ang2 = ang2)
                    pwK2 = derivative(partK2, x, P, ang, 'w', ang2 = ang2)

                    H1 = [derivative(partV0, x, P, ang, 'V0', ang2 = ang2) * (1+l), pVK, pVw, pVe, pVT0, pVV2, pVK2]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K', ang2 = ang2) * (1+l), pKw, pKe, pKT0, pKV2, pKK2]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w', ang2 = ang2) * (1+l), pwe, pwT0, pwV2, pwK2]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e', ang2 = ang2) * (l+1), peT0, peV2, peK2]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0', ang2 = ang2) * (1+l), pT0V2, pT0K2]
                    H6 = [pVV2, pKV2, pwV2, peV2, pT0V2, derivative(partV02, x, P, ang, 'V02', ang2 = ang2) * (1+l), pV2K2]
                    H7 = [pVK2, pKK2, pwK2, peK2, pT0K2, pV2K2, derivative(partK2, x, P, ang, 'K2', ang2 = ang2) * (1+l)]

                    a = 0.5 * np.array([H1,H2,H3,H4,H5,H6,H7])
                    b = -0.5 * np.array([partV0(x, ang, P, ang2 = ang2), partK(x, ang, P, ang2 = ang2), partw(x, ang, P, ang2 = ang2),
                                         parte(x, ang, P, ang2 = ang2), partT0(x, ang, P, ang2 = ang2),
                                         partV02(x, ang, P, ang2 = ang2), partK2(x, ang, P, ang2 = ang2)])

                else:
                    H1 = [derivative(partV0, x, P, ang, 'V0', ang2 = ang2) * (1+l), pVK, pVT0, pVV2, pVK2]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K', ang2 = ang2) * (1+l), pKT0, pKV2, pKK2]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0', ang2 = ang2) * (1+l), pT0V2, pT0K2]
                    H4 = [pVV2, pKV2, pT0V2, derivative(partV02, x, P, ang, 'V02', ang2 = ang2) * (1+l), pV2K2]
                    H5 = [pVK2, pKK2, pT0K2, pV2K2, derivative(partK2, x, P, ang, 'K2', ang2 = ang2) * (1+l)]
                    a = 0.5 * np.array([H1,H2,H3,H4,H5])
                    b = -0.5 * np.array([partV0(x, ang, P, ang2 = ang2), partK(x, ang, P, ang2 = ang2),
                                          partT0(x, ang, P, ang2 = ang2),
                                         partV02(x, ang, P, ang2 = ang2), partK2(x, ang, P, ang2 = ang2)])
            x = xLast + np.linalg.solve(a, b) # Next iterations parameter values

            if period == None:
                if zeroEcc:
                    V0, K, T0, P, V02, K2 = x
                else:
                    V0, K, w, e, T0, P, V02, K2 = x
            else:
                if zeroEcc:
                    V0, K, T0, V02, K2 = x
                else:
                    V0, K, w, e, T0, V02, K2 = x

            ct +=1
            if ct >= 1100 :
                raise ConvergenceError
            if SumSquared(x, P) < SumSquared(xLast, P): # if params got better
                if abs(SumSquared(xLast, P) - SumSquared(x, P)) < 0.01 or abs(1-(SumSquared(x, P)/SumSquared(xLast, P))) < 1e-3:
                    break # convergence acheived! (probably)
                else:
                    l /= 9 # lowering damping parameter
                    l = max(l, 10e-7)
            else: # if params got worse
                l *= 11 # raising damping parameter
                l = min(l, 10e7)

        if not zeroEcc:
            x[2] = x[2]%(2*np.pi)

        if period == None:
            if zeroEcc:
                P = x[3]
            else:
                P = x[5]
        x = [V0,K,w,e,T0,P,V02,K2]
        if not gam2set:
            params = V0,K,w,e,T0,P,V02,K2
        else:
            V02 = V0
            params = V0,K,w,e,T0,P,V0,K2
        ct = 0

        while True:
            ang = true(time, T0, e, P, ec = True)
            ang2 = true(time2, T0, e, P, ec = True)
            paramsLast = params

            if not gam2set:
                varLast = params[-2], params[-1]

                pVK = derivative(partK2, params, P, ang, 'V02',ang2 = ang2 )

                # Building Hessian matrix
                H1 = [derivative(partV02, params, P, ang, 'V02', ang2 = ang2 ), pVK]
                H2 = [pVK, derivative(partK2, params, P, ang, 'K2', ang2 = ang2 )]
                a2 = 0.5 * np.array([H1, H2])
                b2 = -0.5 * np.array([partV02(params, ang, P,ang2 = ang2 ), partK2(params, ang, P,ang2 = ang2 )])

                varnew = varLast + np.linalg.solve(a2, b2)

                params = V0,K, w, e, T0, P,varnew[0], varnew[1]
            else:
                varLast = params[-1]
                a2 = derivative(partK2, params, P, ang, 'K2', ang2 = ang2 )
                b2 = partK2(params, ang, P,ang2 = ang2 )
                varnew = varLast - (b2/a2)
                params = V0,K, w, e, T0, P,V02,varnew

            ct +=1
            if ct >= 1100 :
                raise ConvergenceError
            if abs(SumSquared(paramsLast, None) - SumSquared(params, None)) < 0.01:
                    break
        x = list(params)
        finalSS = SumSquared(x, None)
        if trialOnly:
            return x, finalSS

        # other values to be returned
        nsec = 2*np.pi / (24*3600*P)
        a1 = (x[1]*np.sqrt(1-x[3]**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2
        x.append(a1)
        x.append(fm)

        C = np.linalg.inv(a) # covariance matrix
        err = []

        for q in range(var):
            err.append(np.sqrt((C[q][q])))

        if zeroEcc:
            err.insert(2,0)
            err.insert(3,0)

        if period != None:
                err.insert(5,0  )

        if gam2:
            err[-2] = err[0]
            C[-2][-2] = err[-2]**2

        # Finding uncertainty in a1
        dade = (-x[1] * x[3] * x[5] * (3600*24)) / (2 * np.pi * np.sqrt(1-x[3]**2))
        dadT = (x[1] * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        dadK = (x[5] * (3600*24) * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        aErr = np.sqrt((dadT * err[5] * (3600*24))**2 + (dadK * err[1])**2 + (dade * err[3])**2)
        aErr *= (10e3 / 10e9)

        # Finding uncertainty in fm
        dfda = 3 * (x[8] * (10e9 / 10e3) )**2 * (3.985e-20) * (x[5])**-2
        dfdT = -2 * (x[8] * (10e9 / 10e3) )**3 * (3.985e-20) * (x[5])**-3
        fErr = np.sqrt((dfda * aErr * (10e9 / 10e3) )**2 + (dfdT * err[5])**2)

        err.append(aErr)
        err.append(fErr)

        # create plots if expected
        webplot = np.empty((1000,6))
        if webgrapher:

            times = np.linspace(min(min(time),min(time2)), max(max(time), max(time2)),1000) # times for fit
            v = true(times, x[4], x[3], P) # array of true anomalies at every time for fit
            fit = x[0] + x[1]*(np.cos(v + x[2])+ x[3]*np.cos(x[2]))
            fit2 = x[6] + x[7]*(np.cos(v + x[2] + np.pi)+ x[3]*np.cos(x[2] + np.pi))

            # RV vs phase data
            phase = ((time - x[4]) / P) % 1 # phase for every RV data point
            # phases = np.linspace(min(phase), max(phase), N)
            phases = ((times - x[4]) / P) % 1 # phase for every point on fit
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            phasesfit2 = []
            for i in range(len(phases)):
                phasesfit.append([phases[i], fit[i]])
                phasesfit2.append([phases[i], fit2[i]])
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))
            phasesfit2 = np.array(sorted(phasesfit2,key=lambda l:l[0]))

            webplot[:,0] = times
            webplot[:,1] = fit
            webplot[:,2] = phasesfit[:,0]
            webplot[:,3] = phasesfit[:,1]
            webplot[:,4] = fit2
            webplot[:,5] = phasesfit2[:,1]

        # rounding
        V0 = sigFig(x[0], 5)
        V02 = sigFig(x[6], 5)
        K = sigFig(x[1], 5)
        K2 = sigFig(x[7], 5)
        if not zeroEcc:
            w = sigFig((x[2] * 180/np.pi)%360, 6)
            e = sigFig(x[3], 5)
        else:
            w,e = 0,0
        T0 = sigFig(x[4], 6)
        a1 = sigFig(x[8], 5)
        fm = sigFig(x[9], 4)

        V0err = sigFig(err[0], 5)
        Kerr = sigFig(err[1], 5)
        werr = err[2]
        eerr = err[3]
        if not zeroEcc:
             werr = sigFig((err[2] * 180/np.pi)%360, 6)
             eerr = sigFig(err[3], 5)

        T0err = sigFig(err[4], 6)

        V02err = sigFig(err[6], 5)
        K2err = sigFig(err[7], 5)

        a1err = sigFig(err[8], 5)
        fmerr = sigFig(err[9], 5)

        if period == None:
            P = sigFig(x[5], 6)
            Perr = sigFig(err[5], 6)
        else:
            P = x[5]
            Perr = err[5]

        x = [V0, K, w, e, T0, P, V02, K2, a1,fm]
        err = [V0err, Kerr, werr, eerr, T0err, Perr,V02err, K2err, a1err, fmerr]

        # return expected values
        if covariance:

            for i in range(len(C)):
                for j in range(len(C)):
                    C[i][j] = sigFig(C[i][j], 5)
            if webgrapher:
                result = (x, err, C,webplot)
            else:
                result = (x, err, C)
        else:
            result = (x, err)

        if returnScore:
            return result, finalSS
        return result

    def CompanionSolve(X, err = [None, None, None, None, None, None], shift = True, graphs = True, webgrapher = False):
        """
        Solves orbital parameters for a companion star

        Keyword Arguments:
        X - List of elements of the primary star, in the order [γ, K, w, e, T0, P].
            (with w in degrees). K does not necessarily need to be known but a
            value must be input.
        err (optional) - List of errors on the elements in X (in the same order)
        shift (optional) - If shift = True, returns discrepancy between γ of primary
                           star and of companion star
        graphs (optional) - If graphs = True, companion creates 2 plots. First, a
                            plot of the radial velocity data as a function of time
                            (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                            where [γ, K, ω, e, T0] are the orbital parameters
                            solved for by companion(). Second, a plot of the
                            radial velocity data as a function of phase (between 0
                            and 1) of the star's orbit. Again, fit with
                            V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                            function of phase. On both plots, the long-term average
                            velocity, γ, is also shown.

        Returns:
        x - List of solved elements of the companion star, in the order
            [γ, K, w, e, T0, P, asin(i), f(m)].
        err (optional) - If errors associated with primary star are supplied, list
                         of errors on the elements in x
        shift (optional) - Difference in average radial velocity of primary star, γ1,
                           and of secondary star, γ2 (if argument shift = True). The
                           radial velocity of the primary star is taken to be the
                           'correct' γ in x. Thus shift returned = γ2 - γ1.
        """

        V0, K, wdeg, e, T0, P = X
        wdeg = (wdeg + 180)%360
        w = wdeg*np.pi/180
        if not isinstance(err,list):
            err = [None]

        K = (max(V) - min(V))/2

        params = V0, K, w, e, T0, P
        ang = true(time, T0, e, P, ec = True)

        ct = 0
        while True:
            paramsLast = params
            varLast = params[0], params[1]

            pVK = derivative(partK, params, P, ang, 'V0')

            # Building Hessian matrix
            H1 = [derivative(partV0, params, P, ang, 'V0'), pVK]
            H2 = [pVK, derivative(partK, params, P, ang, 'K') ]
            a = 0.5 * np.array([H1, H2])
            b = -0.5 * np.array([partV0(params, ang, P), partK(params, ang, P)])

            var = varLast + np.linalg.solve(a, b)

            params = var[0], var[1], w, e, T0, P

            ct +=1
            if ct >= 1100 :
                raise ConvergenceError
            if abs(SumSquared(paramsLast, None) - SumSquared(params, None)) < 0.01:
                    break

        K = params[1]

        if shift:
            dV0 = sigFig((params[0] - V0),5)

        nsec = 2*np.pi / (P*24*3600)
        a1 = (K*np.sqrt(1-e**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2

        K = sigFig(K, 5)
        a1 = sigFig(a1, 5)
        fm = sigFig(fm, 4)

        x = [V0, K, w, e, T0, P, a1, fm]

        webplot = np.empty((1000,4))
        if webgrapher:

            times = np.linspace(min(time), max(time),1000) # times for fit
            v = true(times, x[4], x[3], P) # array of true anomalies at every time for fit
            fit = x[0] + x[1]*(np.cos(v + x[2])+ x[3]*np.cos(x[2]))

            # RV vs phase data
            phase = ((time - x[4]) / P) % 1 # phase for every RV data point
            # phases = np.linspace(min(phase), max(phase), N)
            phases = ((times - x[4]) / P) % 1 # phase for every point on fit
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            for i in range(len(phases)):
                if i < len(V):
                    phaseV.append([phase[i], V[i]])
                phasesfit.append([phases[i], fit[i]])
            phaseV = np.array(sorted(phaseV,key=lambda l:l[0]))
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))

            webplot[:,0] = times
            webplot[:,1] = fit
            webplot[:,2] = phasesfit[:,0]
            webplot[:,3] = phasesfit[:,1]

        C = np.linalg.inv(a)

        if err[0] != None: # if user inputted errors
            err[1] += np.sqrt(C[1][1]) # adding error from newtons method onto K error
            err[0] += np.sqrt(C[0][0]) # adding error on shift to gamma error (just for web app)

            dade = (-x[1] * x[3] * x[5] * (3600*24)) / (2 * np.pi * np.sqrt(1-x[3]**2))
            dadT = (x[1] * np.sqrt(1-x[3]**2)) / (2 * np.pi)
            dadK = (x[5] * (3600*24) * np.sqrt(1-x[3]**2)) / (2 * np.pi)
            aErr = np.sqrt((dadT * err[5] * (3600*24))**2 + (dadK * err[1])**2 + (dade * err[3])**2)
            aErr *= (10e3 / 10e9)

            # Finding uncertainty in fm
            dfda = 3 * (x[6] * (10e9 / 10e3) )**2 * (3.985e-20) * (x[5])**-2
            dfdT = -2 * (x[6] * (10e9 / 10e3) )**3 * (3.985e-20) * (x[5])**-3
            fErr = np.sqrt((dfda * aErr * (10e9 / 10e3) )**2 + (dfdT * err[5])**2)

            err.append(aErr)
            err.append(fErr)

            x[2] = sigFig(wdeg, 6)

            V0err = sigFig(err[0], 5)
            Kerr = sigFig(err[1], 5)
            werr = sigFig(err[2],6)
            eerr = sigFig(err[3], 5)
            T0err = sigFig(err[4], 6)
            Perr = sigFig(err[5], 6)
            a1err = sigFig(err[6], 5)
            fmerr = sigFig(err[7], 5)

            if shift:
                output = (x, [V0err,Kerr,werr,eerr,T0err,Perr,a1err,fmerr], dV0)
                if webgrapher:
                    return output + (webplot,)
                return output
            else:
                return x, [V0err,Kerr,werr,eerr,T0err,Perr,a1err,fmerr]

        else:
            x[2] = sigFig(wdeg, 6)
            if shift:
                if webgrapher:
                    return x, dV0, webplot
                return x, dV0
            else:
                return x

    # setting up data for special case of stars = "both"
    if star == "both":
        V1 = np.copy(V) # primary radial velocities
        time1 = np.copy(time)
        V2 = data_sorted[:,2] # companion radial velocities
        time2 = np.copy(time)
        # weights for both stars
        if data.shape[1] >= 4:
            weight1 = data_sorted[:,3]
        else:
            weight1 = np.ones(N)
        if data.shape[1] >= 5:
            weight2 = data_sorted[:,4]
        else:
            weight2 = np.ones(N)
        weight = np.concatenate((weight1,weight2))

        if data.shape[1] >= 6:
            setNum = data_sorted[:,5]
        else:
            linFit = lsq_bisect(V1,V2,weight1,weight2)
            beta, alpha = linFit[0], linFit[1]
            gamma = alpha / (1 - beta)
            time = np.concatenate((time,time))
            V = np.concatenate((V,(V2-gamma)/beta + gamma))
            timecomb = np.copy(time) # holds both stars' times
            Vcomb = np.copy(V) # holds both stars' RV
            weightcomb = np.copy(weight) # holds both stars' weights

    def BothStars(period = None, Pguess = None,covariance = False, graphs = False, correction = False, webgrapher = False, initialX = None, multiStart = True, returnScore = False, trialOnly = False):
        """
        Solves for orbital parameters from radial velocity data, for both stars at once

        Keyword Arguments:
        period (optional) - Known value for period (if none provided, by default
                      PrimarySolve makes its own guess)
        covariance (optional) - If covariance = True, PrimarySolve will return the
                                entire covariance matrix.
        graphs (optional) - If graphs = true, PrimarySolve creates 2 plots. First, a
                            plot of the radial velocity data as a function of time
                            (RJD), fitted with V(t) = γ + K(cos(v(t) + ω) + e*cos(ω)
                            where [γ, K, ω, e, T0] are the orbital parameters
                            solved for by PrimarySolve(). Second, a plot of the
                            radial velocity data as a function of phase (between 0
                            and 1) of the star's orbit. Again, fit with
                            V = γ + K(cos(v + ω) + e*cos(ω) where now V and v are a
                            function of phase. On both plots, the long-term average
                            velocity, γ, is also shown.

        Returns:
        x - List of solved orbital parameters, the semi major axis and the
            mass function (for both stars), in the form:
            [[γ, K1, ω1, e, T0, P, a1, f1(M)], [γ, K2, ω2, e, T0, P, a1, f2(M)]]
        err - List of standard errors on the determination of the orbital elements,
              in the order [γ, K, ω, e, T0, P, a, f(M)]
        C - Estimated covariance matrix (if argument covariance = True)

        # """

        requestedPeriod = Pguess if Pguess != None else period
        poorCoverage = hasPoorPeriodCoverage(requestedPeriod, starMode = "both")
        if initialX is None and multiStart and poorCoverage:
            trialGuesses = initialGuessCandidates(requestedPeriod, starMode = "both", gam = gamma)
            bestTrialResult = None
            bestSS = None
            for trialGuess in trialGuesses:
                try:
                    trialResult, trialSS = BothStars(
                        period = period, Pguess = Pguess, covariance = covariance,
                        graphs = graphs, correction = correction,
                        webgrapher = webgrapher,
                        initialX = trialGuess, multiStart = False,
                        returnScore = True)
                except (ConvergenceError, np.linalg.LinAlgError, FloatingPointError,
                        ValueError, TypeError, ZeroDivisionError):
                    continue
                if validTrialResult(trialResult, starMode = "both") and np.isfinite(trialSS) and (bestSS is None or trialSS < bestSS):
                    bestTrialResult = trialResult
                    bestSS = trialSS
            if bestTrialResult is not None:
                if returnScore:
                    return bestTrialResult, bestSS
                return bestTrialResult
            raise ConvergenceError("All initial guesses failed to converge.")

        if Pguess != None:
            period = Pguess
        if initialX is None:
            x = x0(period, star = "both",gam = gamma, zeroEcc = zeroEcc) # initial guess
        else:
            x = list(initialX)
        if not validInitialGuess(x):
            raise ConvergenceError("Unable to construct a valid initial guess.")

        if Pguess != None:
            period = None

        if period != None: # if period known and provided by user
            P = period
            x.pop()
            V0, K, w, e, T0 = x
            var = 5
        else:
            V0, K, w, e, T0, P = x
            var = 6

        if zeroEcc:
            if period != None:
                x = V0, K, T0
                var = 3
            else:
                x = V0, K,T0,P
                var = 4
        l = 3
        ct = 0
        while True: # perfoming minimization

            if  e < 0 and zeroEcc == False:
                e = 0 # so some sqrts dont cause us troubles
                x[3] = e
            xLast = x
            ang = true(timecomb, T0, e, P, ec = True)

            pVK = derivative(partK, x, P, ang, 'V0')
            pVT0 = derivative(partT0, x, P, ang, 'V0')
            pKT0 = derivative(partT0, x, P, ang, 'K')

            if period == None:

                # Building Hessian matrix

                pVP = derivative(partV0, x, P, ang, 'P')
                pKP = derivative(partK, x, P, ang, 'P')
                pT0P = derivative(partT0, x, P, ang, 'P')

                if not zeroEcc:
                    peP = derivative(parte, x, P, ang, 'P')
                    pVe = derivative(parte, x, P, ang, 'V0')
                    pKe = derivative(parte, x, P, ang, 'K')
                    pwe = derivative(parte, x, P, ang, 'w')
                    peT0 = derivative(partT0, x, P, ang, 'e')
                    pVw = derivative(partw, x, P, ang, 'V0')
                    pKw = derivative(partw, x, P, ang, 'K')
                    pwT0 = derivative(partT0, x, P, ang, 'w')
                    pwP = derivative(partw, x, P, ang, 'P')

                    H1 = [derivative(partV0, x, P, ang, 'V0') * (1+l), pVK, pVw, pVe, pVT0, pVP]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKw, pKe, pKT0, pKP]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w') * (1+l), pwe, pwT0, pwP]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e') * (l+1), peT0, peP]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0') * (1+l), pT0P]
                    H6 = [pVP, pKP, pwP, peP, pT0P, derivative(partP, x, P, ang, 'P') * (1+l)]
                    a = 0.5 * np.array([H1, H2, H3, H4, H5, H6])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P), partw(x, ang, P),
                                         parte(x, ang, P), partT0(x, ang, P), partP(x, ang, P)])

                else:
                    H1 = [derivative(partV0, x, P, ang, 'V0') * (1+l), pVK, pVT0, pVP]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKT0, pKP]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0') * (1+l), pT0P]
                    H4 = [pVP, pKP, pT0P, derivative(partP, x, P, ang, 'P') * (1+l)]
                    a = 0.5 * np.array([H1, H2, H3, H4])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P),
                                         partT0(x, ang, P), partP(x, ang, P)])

            else:

                if not zeroEcc:

                    pVe = derivative(parte, x, P, ang, 'V0')
                    pKe = derivative(parte, x, P, ang, 'K')
                    pwe = derivative(parte, x, P, ang, 'w')
                    peT0 = derivative(partT0, x, P, ang, 'e')
                    pVw = derivative(partw, x, P, ang, 'V0')
                    pKw = derivative(partw, x, P, ang, 'K')
                    pwT0 = derivative(partT0, x, P, ang, 'w')

                    H1 = [derivative(partV0, x, P, ang, "V0") * (1+l), pVK, pVw, pVe, pVT0]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKw, pKe, pKT0]
                    H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w') * (1+l), pwe, pwT0]
                    H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e') * (l+1), peT0]
                    H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0') * (1+l)]
                    a = 0.5 * np.array([H1,H2,H3,H4,H5])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P), partw(x, ang, P),
                                         parte(x, ang, P), partT0(x, ang, P)])

                else:
                    H1 = [derivative(partV0, x, P, ang, "V0") * (1+l), pVK, pVT0]
                    H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKT0]
                    H3 = [pVT0, pKT0, derivative(partT0, x, P, ang, 'T0') * (1+l)]
                    a = 0.5 * np.array([H1,H2,H3])
                    b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P),
                                         partT0(x, ang, P)])

            x = xLast + np.linalg.solve(a, b) # Next iterations parameter values

            if period == None:
                if zeroEcc:
                    V0, K, T0, P = x
                else:
                    V0, K, w, e, T0, P = x
            else:
                if zeroEcc:
                    V0, K, T0 = x
                else:
                    V0, K, w, e, T0 = x

            ct +=1
            if ct >= 1100 :
                raise ConvergenceError
            if SumSquared(x, P) < SumSquared(xLast, P): # if params got better
                if abs(SumSquared(xLast, P) - SumSquared(x, P)) < 0.01 or abs(1-(SumSquared(x, P)/SumSquared(xLast, P))) < 1e-3:
                    break # convergence acheived! (probably)
                else:
                    l /= 9 # lowering damping parameter
                    l = max(l, 10e-7)
            else: # if params got worse
                l *= 11 # raising damping parameter
                l = min(l, 10e7)

        if not zeroEcc:
            if x[2] < -np.pi or x[2] >= np.pi:
                while x[2] < -np.pi:
                    x[2] = x[2] + 2*np.pi
                while x[2] >= np.pi:
                    x[2] = x[2] - 2 * np.pi
                w = x[2]

        if period == None:
            if zeroEcc:
                P = x[3]
            else:
                P = x[5]
        x = [V0,K,w,e,T0,P]
        finalSS = SumSquared(x, None)
        if trialOnly:
            return x, finalSS

        # other values to be returned
        nsec = 2*np.pi / (24*3600*P)
        a1 = (x[1]*np.sqrt(1-x[3]**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2
        x.append(a1)
        x.append(fm)
        K2 = -1*beta * K

        # if period == None:
        #     x = [x[0], x[1], x[2], x[3], x[4], x[5], a1, fm]
        # else:
        #     x = [x[0], x[1], x[2], x[3], x[4], P, a1, fm]

        C = np.linalg.inv(a) # covariance matrix
        err = []

        for q in range(var):
            err.append(np.sqrt(C[q][q]))
        if period != None:
            err.append(0)
        if zeroEcc:
            err.insert(2,0)
            err.insert(3,0)

        # Finding uncertainty in a1
        dade = (-x[1] * x[3] * x[5] * (3600*24)) / (2 * np.pi * np.sqrt(1-x[3]**2))
        dadT = (x[1] * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        dadK = (x[5] * (3600*24) * np.sqrt(1-x[3]**2)) / (2 * np.pi)
        aErr = np.sqrt((dadT * err[5] * (3600*24))**2 + (dadK * err[1])**2 + (dade * err[3])**2)
        aErr *= (10e3 / 10e9)

        # Finding uncertainty in fm
        dfda = 3 * (x[6] * (10e9 / 10e3) )**2 * (3.985e-20) * (x[5])**-2
        dfdT = -2 * (x[6] * (10e9 / 10e3) )**3 * (3.985e-20) * (x[5])**-3
        fErr = np.sqrt((dfda * aErr * (10e9 / 10e3) )**2 + (dfdT * err[5])**2)

        err.append(aErr)
        err.append(fErr)

        # rounding

        V0 = sigFig(x[0], 5)
        K = sigFig(x[1], 5)

        if not zeroEcc:
            w = sigFig((x[2] * 180/np.pi)%360, 6)
            if e < 0.1:
                e = sigFig(x[3], 4)
            else:
                e = sigFig(x[3], 4)
        T0 = sigFig(x[4], 6)
        a1 = sigFig(x[6], 5)
        fm = sigFig(x[7], 4)

        V0err = sigFig(err[0], 5)
        Kerr = sigFig(err[1], 5)
        werr = err[2]
        eerr = err[3]
        if not zeroEcc:
             werr = sigFig((err[2] * 180/np.pi)%360, 6)
             eerr = sigFig(err[3], 5)
        T0err = sigFig(err[4], 6)
        a1err = sigFig(err[6], 5)
        fmerr = sigFig(err[7], 5)

        if period == None:
            P = sigFig(x[5], 6)
            Perr = sigFig(err[5], 6)
        else:
            P = x[5]
            Perr = err[5]

        a2 = (K2*np.sqrt(1-e**2) / nsec) * 10e3 / 10e9
        fm2 = 3.985e-20 * (a2 * 10e9 / 10e3)**3 / P**2

        w2 = (w + 180)%360

        # there are twice as any params to return because stars = "both"!
        x = np.array([[V0, K, w, e, T0, P, a1, fm], [V0, K2, w2, e, T0, P, a2, fm2]])
        err = [V0err, Kerr, werr, eerr, T0err, Perr, a1err, fmerr]
        err = np.array([err,err])
        C = np.array([C,C])

        #now going through one more time to solve Ks and gammas better
        #first star 1
        X = x[0][0:6]
        X[2] = np.radians(X[2])%(2*np.pi)

        ang = true(time1, T0, e, P, ec = True)
        ct = 0
        while True:

            Xlast = X
            varLast = np.array([X[0], X[1]])

            pVK = derivative(partK, X, P, ang, 'V0', starNum = 1)

            # Building Hessian matrix
            H1 = [derivative(partV0, X, P, ang, 'V0',starNum = 1), pVK]
            H2 = [pVK, derivative(partK, X, P, ang, 'K',starNum = 1) ]
            a = 0.5 * np.array([H1, H2])
            b = -0.5 * np.array([partV0(X, ang, P,starNum = 1), partK(X, ang, P, starNum = 1)])

            var = varLast + np.linalg.solve(a, b)

            X = [var[0], var[1], X[2], e, T0, P]

            ct += 1
            if ct >= 1100:
                raise ConvergenceError("Primary-star refinement failed to converge.")
            if abs(SumSquared(Xlast, None, starNum = 1) - SumSquared(X, None, starNum = 1)) < 0.01:
                    break

        K = X[1]
        V0 = X[0]

        nsec = 2*np.pi / (P*24*3600)
        a1 = (K*np.sqrt(1-e**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2

        K = sigFig(K, 5)
        V0 = sigFig(V0,5)
        a1 = sigFig(a1, 5)
        fm = sigFig(fm, 4)

        x[0] = [V0, K,w, e, T0, P, a1, fm]

        C1 = np.linalg.inv(a) # covariance matrix for only star 1 (V0 and K1)
        err[0][0] = sigFig(np.sqrt(C1[0][0]), 5)
        err[0][1] = sigFig(np.sqrt(C1[1][1]), 5)
        C[0][0][0] = (C1[0][0])
        C[0][0][1] = (C1[0][1])
        C[0][1][0] = (C1[1][0])
        C[0][1][1] = (C1[1][1])

        # now doing the other star
        X = x[1][0:6]
        X[2] = np.radians(X[2])%(2*np.pi)

        ang = true(time2, T0, e, P, ec = True)
        ct = 0
        while True:

            Xlast = X
            varLast = np.array([X[0], X[1]])
            pVK = derivative(partK, X, P, ang, 'V0', starNum = 2)

            # Building Hessian matrix
            H1 = [derivative(partV0, X, P, ang, 'V0',starNum = 2), pVK]
            H2 = [pVK, derivative(partK, X, P, ang, 'K',starNum = 2) ]
            a = 0.5 * np.array([H1, H2])
            b = -0.5 * np.array([partV0(X, ang, P,starNum = 2), partK(X, ang, P,starNum = 2)])

            var = varLast + np.linalg.solve(a, b)

            X = [var[0], var[1], X[2], e, T0, P]

            ct += 1
            if ct >= 1100:
                raise ConvergenceError("Companion-star refinement failed to converge.")
            if abs(SumSquared(Xlast, None,starNum = 2) - SumSquared(X, None,starNum = 2)) < 0.01:
                    break

        K2 = X[1]
        V0 = X[0]

        a2 = (K2*np.sqrt(1-e**2) / nsec) * 10e3 / 10e9
        fm2 = 3.985e-20 * (a2 * 10e9 / 10e3)**3 / P**2

        K2 = sigFig(K2, 5)
        V0 = sigFig(V0,5)

        a2 = sigFig(a2, 5)
        fm2 = sigFig(fm2, 5)

        x[1] = [V0, K2, w2, e, T0, P, a2, fm2]

        C2 = np.linalg.inv(a) # covariance matrix for only star 2 (V02 and K2)
        err[1][0] = sigFig(np.sqrt(C2[0][0]), 5)
        err[1][1] = sigFig(np.sqrt(C2[1][1]), 5)
        C[1][0][0] = (C2[0][0])
        C[1][0][1] = (C2[0][1])
        C[1][1][0] = (C2[1][0])
        C[1][1][1] = (C2[1][1])

        webplot = np.empty((1000,6))
        if webgrapher:

            times = np.linspace(min(time1), max(time1),1000) # times for fit
            v = true(times, x[0][4], x[0][3], P) # array of true anomalies at every time for fit
            fit = x[0][0] + x[0][1]*(np.cos(v + np.radians(x[0][2]))+ x[0][3]*np.cos(x[0][2]))
            fit2 = x[1][0] + x[1][1] * (np.cos(v + np.radians(x[1][2]))+ x[1][3]*np.cos(x[1][2]))

            # RV vs phase data
            phase = ((time - x[0][4]) / P) % 1 # phase for every RV data point
            # phases = np.linspace(min(phase), max(phase), N)
            phases = ((times - x[0][4]) / P) % 1 # phase for every point on fit
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            for i in range(len(phases)):
                if i < len(V1):
                    phaseV.append([phase[i], V[i]])
                phasesfit.append([phases[i], fit[i]])
            phaseV = np.array(sorted(phaseV,key=lambda l:l[0]))
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))

            webplot[:,0] = times
            webplot[:,1] = fit
            webplot[:,2] = fit2
            webplot[:,3] = phasesfit[:,0]
            webplot[:,4] = phasesfit[:,1]

            phaseV2 = []
            phasesfit2 = []
            for i in range(len(phases)):
                if i < len(V2):
                    phaseV2.append([phase[i], V2[i]])
                phasesfit2.append([phases[i], fit2[i]])
            phaseV2 = np.array(sorted(phaseV2,key=lambda l:l[0]))
            phasesfit2 = np.array(sorted(phasesfit2,key=lambda l:l[0]))
            webplot[:,5] = phasesfit2[:,1]

        # return requested values
        if covariance: #
            for k in range(2):
                for i in range(np.size(C,1)):
                    for j in range(np.size(C,2)):
                        C[k][i][j] = sigFig(C[k][i][j], 5)
            result = (x, err, C, webplot)
        else:
            result = (x, err)

        if returnScore:
            return result, finalSS
        return result

    def meanMatch(X):
        """
        Finds the γ of an offset dataset

        Keyword Arguments:
        X - List of elements of the primary star, in the order [γ, K, w, e, T0, P].
            (with w in degrees).

        Returns:
        V0 - The corrected value of γ
        """
        V0, K, wdeg, e, T0, P = X
        w = wdeg*np.pi/180
        params = V0, K, w, e, T0, P

        ang = true(time, T0, e, P, ec = True)

        ct = 0
        while True: # single variable Newton's method
            paramsLast = params
            V0 = params[0]
            derivative_value = derivative(partV0, params, P, ang, 'V0')
            if not np.isfinite(derivative_value) or abs(derivative_value) < 1e-12:
                raise ConvergenceError("Dataset-offset correction encountered a singular derivative.")
            V0 = V0 - partV0(params, ang, P)/derivative_value
            params = V0, K, w, e, T0, P
            ct += 1
            if ct >= 1100:
                raise ConvergenceError("Dataset-offset correction failed to converge.")
            if abs(SumSquared(paramsLast, None) - SumSquared(params, None)) < 0.01 :
                    break

        V0 = params[0]

        return V0

    #%% Correcting γ offset bewteen sub data sets

    def split():
        """
        Splits a data set into its constituent data sets.
        """
        sets = np.array(list(set(setNum))) # indices of the different data sets
        origSets = np.empty([len(setNum),data.shape[1],len(sets)],float) # 3D array
        origSets[:,:,:] = np.nan
        # seperating the sets
        for i in range(len(setNum)):
            for j in range(len(sets)):
                if setNum[i] == sets[j]:
                    origSets[i,:,j] = data_sorted[i,:]
        mask = np.isfinite(origSets) # mask for which elements have a value
        return (origSets,mask)

    splitcond = False
    if star == "primary":
        splitcond = data.shape[1] >= 4
    elif star == "both":
        splitcond = data.shape[1] >= 6

    if splitcond: # if γ correcting is required
        splitSets, splitMask = split() # splitting sets
        contenders = [] # holds if a data set could be the best

        if star == "primary":

            for j in range(splitSets.shape[2]):
                V = splitSets[:,1,j][splitMask[:,1,j]]
                time = splitSets[:,0,j][splitMask[:,0,j]]
                N = len(time)
                if np.isfinite(PGuess(test = True)).any(): # if a period can be found
                    contenders.append(j)
            if len(contenders) == 0: # if no sub dataset can estimate period
                big = np.array([0,None]) # simply choose the largest set
                for j in range(splitSets.shape[2]):
                    time = splitSets[:,0,j][splitMask[:,0,j]]
                    N = len(time)
                    if N >= big[0]:
                        big[0], big[1] = N,j
                contenders.append(big[1])

            SS = np.inf # will hold sum of squared deviations
            tempContenders = np.copy(contenders)
            for j in tempContenders: # figuring out which contender is best
                time = splitSets[:,0,int(j)][splitMask[:,0,int(j)]]
                N = len(time)
                V = splitSets[:,1,int(j)][splitMask[:,1,int(j)]]
                weight = np.ones(N)
                if Pguess != None:
                    Period = Pguess
                X = x0(Period)
                if (SumSquared(X,None))/N <= SS:
                    SS = (SumSquared(X,None))/N
                    del contenders[:contenders.index(j)]
                else:
                    contenders.remove(j)

            time = splitSets[:,0,contenders[0]][splitMask[:,0,contenders[0]]]
            V = splitSets[:,1,contenders[0]][splitMask[:,1,contenders[0]]]
            weight = splitSets[:,2,contenders[0]][splitMask[:,2,contenders[0]]]
            N = len(time)
            # finding parameters for "best" data set
            X,err = PrimarySolve(period = Period, covariance = False, graphs = False,correction = True)
            sets = list(set(setNum))
            sets.remove(setNum[contenders[0]]) # removing the "best" set from sets
            for i in sets: # shifting other data sets
                V = splitSets[:,1,int(i-1)][splitMask[:,1,int(i-1)]]
                time = splitSets[:,0,int(i-1)][splitMask[:,0,int(i-1)]]
                weight = splitSets[:,2,int(i-1)][splitMask[:,2,int(i-1)]]
                N = len(time)
                V0 = meanMatch(X[:6]) # finding γ for set i
                shift = V0 - X[0] # the difference in mean RV
                # putting shifted data back together
                for j in range(data_sorted.shape[0]):
                    if data_sorted[j,3] == i:
                        data_sorted[j,1] -= shift # shifting set i

        elif star == "both" :
            # will hold all the different gammas, and their weights
            gammaList = np.empty((splitSets.shape[2],2))
            for j in range(splitSets.shape[2]):
                V = splitSets[:,1,j][splitMask[:,1,j]]
                V2 = splitSets[:,2,j][splitMask[:,2,j]]
                time = splitSets[:,0,j][splitMask[:,0,j]]
                N = len(time)
                linFit = np.polyfit(V,V2,1)
                beta, alpha = linFit[0], linFit[1]
                gamma = alpha / (1 - beta) # this set's gamma
                gammaWeight = np.concatenate((splitSets[:,3,j][splitMask[:,3,j]],splitSets[:,4,j][splitMask[:,4,j]]))
                gammaWeight = np.mean(gammaWeight)
                gammaList[j] = [gamma,gammaWeight]

            #finding weighted mean
            gamma = np.sum(gammaList[:,0] * gammaList[:,1]) / np.sum(gammaList[:,1])
            #shifting data sets to be at weighted mean gamma
            sets = list(set(setNum))
            for i in sets:
                shift = gamma - gammaList[sets.index(i)][0]
                # putting shifted data back together
                for j in range(data_sorted.shape[0]):
                    if data_sorted[j,5] == i:
                        # shifting set i
                        data_sorted[j,1] += shift
                        data_sorted[j,2] += shift

        if Pguess != None:
            Period = None
        time = data_sorted[:,0]
        V = data_sorted[:,1]
        Vplot = np.copy(V)
        tplot = np.copy(time)
        weight = data_sorted[:,2]
        N = len(time)

        if star == "both":
            V1 = np.copy(V)
            V2 = data_sorted[:,2]
            linFit = np.polyfit(V,V2,1)
            beta, alpha = linFit[0], linFit[1]
            gamma = alpha / (1 - beta)
            time = np.concatenate((time,time))
            timecomb = np.copy(time)
            V = np.concatenate((V,(V2-gamma)/beta + gamma))
            Vcomb = np.copy(V)
            weight = np.concatenate((data_sorted[:,3],data_sorted[:,4]))
            weightcomb = np.copy(weight)

    else:
        if Pguess != None:
            Period = None
        if star != "both":
            time = data_sorted[:,0]
            V = data_sorted[:,1]
            Vplot = np.copy(V)
            N = len(time)
            tplot = np.copy(time)

        if star == "secondary_both":
            Vplot2 = np.copy(V2)
            tplot2 = np.copy(time2)

    def fitted_table(params, sample_times):
        """Return [time, phase, fitted_velocity] for one component."""
        sample_times = np.asarray(sample_times, dtype=float)
        V0, K, w_deg, e, T0, P = [float(value) for value in params[:6]]
        w_rad = np.deg2rad(w_deg)
        v = true(sample_times, T0, e, P)
        fitted_velocity = V0 + K*(np.cos(v + w_rad) + e*np.cos(w_rad))
        phase = ((sample_times - T0)/P) % 1
        return np.column_stack((sample_times, phase, fitted_velocity))

    def component_parameters(x_values):
        """Convert each fit mode's return structure into one or two components."""
        if star == "both":
            return np.asarray(x_values[0], dtype=float), np.asarray(x_values[1], dtype=float)
        if star == "secondary_both":
            primary = np.asarray(x_values[:6], dtype=float)
            secondary = np.array([
                x_values[6], x_values[7], (x_values[2] + 180.0) % 360.0,
                x_values[3], x_values[4], x_values[5]
            ], dtype=float)
            return primary, secondary
        return (np.asarray(x_values, dtype=float),)

    def requested_phases():
        if phase_grid is None:
            return np.linspace(0.0, 1.0, 1001)

        values = np.asarray(phase_grid, dtype=float)
        if values.ndim == 1 and values.size == 3 and isinstance(phase_grid, tuple):
            phase_minimum, phase_maximum, phase_increment = values
            if phase_increment <= 0:
                raise ValueError("The phase-grid increment must be positive.")
            if phase_maximum < phase_minimum:
                raise ValueError("The phase-grid maximum must not be below its minimum.")
            count = int(np.floor((phase_maximum - phase_minimum)/phase_increment + 1e-12)) + 1
            phases = phase_minimum + np.arange(count)*phase_increment
            if phases.size == 0 or phases[-1] < phase_maximum - 1e-10:
                phases = np.append(phases, phase_maximum)
            return phases

        phases = np.ravel(values)
        if phases.size == 0 or not np.isfinite(phases).all():
            raise ValueError("phase_grid must contain finite phase values.")
        return phases

    def make_model_output(x_values):
        if model_output is None:
            return None

        output_mode = str(model_output).lower()
        components = component_parameters(x_values)

        if output_mode in ("input", "times", "original"):
            if star == "both":
                sample_times = (data_sorted[:,0], data_sorted[:,0])
            elif star == "secondary_both":
                sample_times = (data_sorted[:,0], data_sorted2[:,0])
            else:
                sample_times = (data_sorted[:,0],)
        elif output_mode in ("phase", "grid", "phase_grid"):
            phases = requested_phases()
            sample_times = tuple(component[4] + phases*component[5] for component in components)
        else:
            raise ValueError('model_output must be None, "input", or "phase".')

        tables = tuple(fitted_table(component, times_) for component, times_ in zip(components, sample_times))
        if output_mode in ("phase", "grid", "phase_grid"):
            # Preserve the requested phase labels exactly, including phase=1.
            for table in tables:
                table[:,1] = phases
        return tables[0] if len(tables) == 1 else tables

    def plot_component(params, observed_times, observed_velocities):
        """Reproduce the package's time and phase plots for one component."""
        observed_times = np.asarray(observed_times, dtype=float)
        observed_velocities = np.asarray(observed_velocities, dtype=float)
        number_fit_points = max(500, len(observed_times))
        fit_times = np.linspace(np.min(observed_times), np.max(observed_times), number_fit_points)
        fit_table = fitted_table(params, fit_times)
        data_phase = ((observed_times - float(params[4]))/float(params[5])) % 1
        phase_order = np.argsort(fit_table[:,1])
        gamma_line = np.zeros(100) + float(params[0])

        fig, ax = plt.subplots()
        ax.scatter(observed_times, observed_velocities, s=18)
        ax.plot(fit_table[:,0], fit_table[:,2], color="orange")
        ax.plot(np.linspace(np.min(observed_times), np.max(observed_times), 100), gamma_line,
                ls="--", color="black")
        ax.set_xlabel("RJD")
        ax.set_ylabel("Velocity (km s$^{-1}$)")
        ax.minorticks_on()
        ax.tick_params(right=True, top=True, which="both")

        fig, ax = plt.subplots()
        ax.scatter(data_phase, observed_velocities, s=18)
        ax.plot(fit_table[phase_order,1], fit_table[phase_order,2], color="orange")
        ax.plot(np.linspace(0, 1, 100), gamma_line, ls="--", color="black")
        ax.set_xlabel("Phase")
        ax.set_ylabel("Velocity (km s$^{-1}$)")
        ax.minorticks_on()
        ax.tick_params(right=True, top=True, which="both")

    def make_plots(x_values):
        components = component_parameters(x_values)
        if star == "both":
            plot_component(components[0], data_sorted[:,0], data_sorted[:,1])
            plot_component(components[1], data_sorted[:,0], data_sorted[:,2])
        elif star == "secondary_both":
            plot_component(components[0], data_sorted[:,0], data_sorted[:,1])
            plot_component(components[1], data_sorted2[:,0], data_sorted2[:,1])
        else:
            plot_component(components[0], data_sorted[:,0], data_sorted[:,1])

    def finish(result):
        x_values = result[0] if isinstance(result, tuple) else result
        if graphs:
            make_plots(x_values)
        generated_model = make_model_output(x_values)
        if generated_model is None:
            return result
        if isinstance(result, tuple):
            return result + (generated_model,)
        return result, generated_model

    if star == "primary":
        if covariance:
            result = PrimarySolve(Period, Pguess, covariance, graphs=False)
        else:
            result = PrimarySolve(Period, Pguess, covariance, graphs=False)
        return finish(result)

    elif star == "secondary":
        if X is None:
            raise ValueError('X must be supplied when star="secondary".')
        companion_errors = err if err is not None else [None, None, None, None, None, None]
        result = CompanionSolve(X, companion_errors, shift, graphs=False)
        return finish(result)

    elif star == "both":
        result = BothStars(Period, Pguess, covariance, graphs=False)
        return finish(result)

    elif star == "secondary_both":
        if file2 is None:
            raise ValueError('file2 must be supplied for an unequal-coverage two-star fit.')
        result = CompanionSolve_both(Period, Pguess, covariance, graphs=False, gam2set=gam2)
        return finish(result)

    raise ValueError(
        'star must be "primary", "secondary", "both", or "both_unequal".'
    )
