"""
All of the functionality of this package is through the function StarSolve(). 
Call help(binarystarsolve.binarystarsolve.StarSolve) for the associated Docstring.

For a more detailed description of how to use this package, please read the ReadMe.
https://github.com/NickMilsonPhysics/BinaryStarSolver/blob/master/README.md
"""

# Import statements
import math
import matplotlib.pyplot as plt
import numpy as np
import sys



def StarSolve(data_file,star = "primary",Period=None,Pguess=None, covariance = False,graphs = True, X = None, err = None,shift = False):
    """
    
        Please first read the read me at:
       	https://github.com/NickMilsonPhysics/binarystarsolve/blob/master/README.md
        
        Keyword Arguments: 
            
        ** Note some arguments pertain to when parameters of the primary are
        being solved for, and some pertain to when parameters of the secondary
        are being solved for. If an argument doesn't pertain to the situation
        this function is being used for, it may be ignored, as all the case
        specific arguments are None by default. **
        
        Keyword Arguments mandatory regardless of situation:
        data_file - String holding the name of the tab separated txt file
        star - Can be either "primary" or "secondary". Is "primary" by default.
               star is a string for if the parameters of the primary or secondary
               are being solved for.
               
               
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
            
        
        Returns:
            
        ** Note which variables are returned depends on the application of this
        function. **
        
        Returns from the primary star application:
        x - List of solved orbital parameters, the semi major axis and the 
            mass function [γ, K, ω, e, T0, P, a, f(M)]
        err - List of standard errors on the determination of the orbital elements,
              in the order [γ, K, ω, e, T0, P, a, f(M)]
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
    """
    
    #%% SET UP
    # sorting data chronologically
    data = np.loadtxt(data_file)
    data_sorted = np.array(sorted(data,key=lambda l:l[0]))
    time = data_sorted[:,0] # RJD times
    V = data_sorted[:,1] # radial velocities
    if X != None:
        Period = X[5]
    
    N = len(time) # number of data points
    
    # weights
    if data.shape[1] >= 3:
        weight = data_sorted[:,2]
    else:
        weight = np.ones(N) 
    if data.shape[1] >= 4:
        setNum = data_sorted[:,3]
        
    star = star.lower() # making star not case sensitive

    
    #%% Initial guess
    
    # estimates of period
    def PGuess(test = False):
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
        
        # find all cross points
        Vavg = sum(V)/N
        low = None
        crossid = []
        for i in range(1, N-1):
            if V[i+1] > Vavg and V[i-1] < Vavg or V[i+1] < Vavg and V[i-1] > Vavg:
                crossid.append(i)
        
        if N <=2:
            if test:
                return np.nan
            print("Not enough data to estimate period.")
            sys.exit()
            
        # corrects for double crossings
        l = len(crossid)
        for i in range(l-1):
            a = crossid[i]
            b = crossid[i+1]
            if a+1 != b:
                good = False
                for j in range(a, b):
                    if abs(V[j] - Vavg) > (max(V) - Vavg) / 4:
                        good = True
                        break
                if not good:
                    if abs(V[a] - Vavg) < abs(V[b] - Vavg):
                        crossid[i+1] = a + 1
                    else:
                        crossid[i] = b - 1
        
        # if data covers less than one period            
        if crossid[0] + 1 != crossid[1]:
            if l < 4:
                if test:
                    return np.nan
                print("Not enough data to estimate period.")
                sys.exit()
        else:
            if l < 5:
                if test:
                    return np.nan
                print("Not enough data to estimate period.")
                sys.exit()
        
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
            Psum += time[crossid[i+4]] - time[crossid[i]]
            for j in range(4, l-i, 4):
                if first:
                    low = abs(V[crossid[4]] - V[crossid[0]])
                    besta = time[crossid[0]]
                    bestb = time[crossid[4]]
                    p = 1 # number of periods between points
                    first = False
                if abs(V[crossid[i]] - V[crossid[i+j]]) < low:
                    low = abs(V[crossid[i]] - V[crossid[i+j]])
                    besta = time[crossid[i]]
                    bestb = time[crossid[i+j]]
                    p = j/4
            # if cross point is first in pair:
            # method 1: checks every point 5 + 4n away
            # method 2: adds time to 5th point ahead to Psum
            if crossid[i] + 1 == crossid[i+1]:
                for j in range(5, l-i, 4):
                    if abs(V[crossid[i]] - V[crossid[i+j]]) < low:
                        low = abs(V[crossid[i]] - V[crossid[i+j]])
                        besta = time[crossid[i]]
                        bestb = time[crossid[i+j]]
                        p = (j-1)/2
                Psum += time[crossid[i+5]] - time[crossid[i]]
            # if cross point is second in pair:
            # method 1: checks every point 3 + 4n away
            # method 2: adds time to 3rd point ahead to Psum
            else:
                for j in range(3, l-i, 4):
                    if abs(V[crossid[i]] - V[crossid[i+j]]) < low:
                        low = abs(V[crossid[i]] - V[crossid[i+j]])
                        besta = time[crossid[i]]
                        bestb = time[crossid[i+j]]
                        p = (j+1)/2
                Psum += time[crossid[i+3]] - time[crossid[i]]
            n += 2
        
        # if last point does not belong to pair (will not have been checked with
        # pair of points one period earlier):
        # method 1: checks radial velocity diffence between last point and pair
        # one period before
        # method 2: measures period using last point and both points from one
        # period before, adding to Psum
        if check:
            n += 2
            Psum += time[crossid[r+4]] - time[crossid[r]]
            Psum += time[crossid[r+4]] - time[crossid[r+1]]
            if low != None:
                if abs(V[crossid[r]] - V[crossid[r+4]]) < low:
                    low = abs(V[crossid[r]] - V[crossid[r+4]])
                    besta = time[crossid[r]]
                    bestb = time[crossid[r+4]]
                    p = 1
                if abs(V[crossid[r+1]] - V[crossid[r+4]]) < low:
                    besta = time[crossid[r+1]]
                    bestb = time[crossid[r+4]]
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
                    if V[crossid[i]] > Vavg:
                        maximum = V[crossid[i]]
                        for j in range(crossid[i], crossid[i+1]):
                            if V[j] > maximum:
                                maximum = V[j]
                        maxima.append(time[j])
                    else:
                        minimum = V[crossid[i]]
                        for j in range(crossid[i], crossid[i+1]):
                            if V[j] < minimum:
                                minimum = V[j]
                        minima.append(time[j])
                else:
                    pass
            # add time between maxima/minima to Psum
            for i in range(len(maxima)-1):
                Psum += maxima[i+1] - maxima[i]
                n += 1
            for i in range(len(minima)-1):
                Psum += minima[i+1] - minima[i]
                n += 1
            
            P2 = Psum/n # period from second method    
        
        return [P1, P2] 
           
    
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
            while abs(Kep(t, T0, E, e, P)) > 0.000005:
                E = E - Kep(t, T0, E, e, P)/dKep(E, e)
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
                while abs(Kep(t[i], T0, E, e, P)) > 0.000005:
                    E = E - Kep(t[i], T0, E, e, P)/dKep(E, e)
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
        
    def x0(period = None):
        """
        Makes initial approximation for the orbital parameters γ, K, ω, e, T0, & P.
        
        X0 can find an initial estimate for period, but a user supplied estimate 
        can be inputted using the keyword argument 'period'
        
        Note that, in the following two scenerios, x0() may fail to make a good 
        guess:
            i) If the RV data does not cover a full period 
            ii) If the RV data is very non-uniformly distributed. E.g. there are 
                significantly more points at the peak than at the trough. 
        """
        if period == None:
            P1, P2 = PGuess()

        
        for j in range(2): # peforms twice as there are two estimates of period
            if period == None:
                if j == 0:
                    P = P1                 
                elif j == 1:
                    P = P2
            else:
                P = period
                j = 2
                
            K = (max(V) - min(V))/2
            V0 = sum(V) / len(V)
            
            lowSS = None # holds lowest sum of squared deviations
            bestX = [None, None, None, None, None, None] # holds best estimate
            
            for i in np.arange(0.01, 0.95, 0.01): # iterates through eccentricities
                e = i
                temp = (K*e)**-1 * ((max(V) + min(V))* 0.5 - V0)
                if abs(temp) <= 1 :
                    w1 = np.arccos(temp)
                    x1 = [V0, K, w1, e, epoch(w1,e,P), P]
                    w2 = -w1 % (2 * np.pi) # Non principle value
                    x2 = [V0, K, w2, e, epoch(w2,e,P), P]
                    
                    if SumSquared(x1, None) < SumSquared(x2, None):
                        w = w1
                    else:
                        w = w2
                    
                    T0 = epoch(w, e, P)
                    x = [V0, K, w, e, T0, P]
                    if lowSS == None:
                        lowSS = SumSquared(x, None)
                        bestX = x.copy()
                        
                    elif SumSquared(x, None) < lowSS:
                        lowSS = SumSquared(x, None)
                        bestX = x.copy()
            if j == 0:
                xa = bestX.copy()
            elif j == 1:
                xb = bestX.copy()
            else:
                return bestX
        
        if SumSquared(xa, None) < SumSquared(xb, None):
            bestX = xa
        elif SumSquared(xb, None) < SumSquared(xa, None): 
            bestX = xb
        
        return bestX
    
    #%% FIT OF DATA
    
    # rounding function. Makes final output neater.
    def sigFig(x, figs):
        return round(x, -int(math.floor(math.log10(abs(x))) - (figs - 1)))
    
    
    def derivative(f, params, P, ang, der, h=0.000001):
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
        """
        
        
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
        
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
            return (f(params, angAb, P+h) - f(params,angBe, P+h))/(2*h)
        
        return (f(paramsAb, angAb, P) - f(paramsBe, angBe, P))/(2*h)
    
    # Sum of squared deviations between the fit and the actual RV data
    def SumSquared(params, P):
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
        v = true(time, T0, e, P)
        epsilon = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)**2
        return sum(epsilon*weight)
    
    def partV0(params, ang, P = None):
        """
        Derivative of residual wrt γ.
        
        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v = ang[0]
        
        temp = ((V0 + ((np.cos(v + w) + e*np.cos(w)))*K) - V) *weight
        return sum(2*temp)
    
    def partK(params, ang, P = None):
        """
        Derivative of residual wrt K.
        
        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v = ang[0]
        
        temp = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*(weight)
        temp *= (np.cos(v+w) + e*np.cos(w))
        return sum(2*temp)
    
    def partw(params, ang, P = None):
        """
        Derivative of residual wrt w.
        
        Keyword Arguments:
        params - value of parameters [γ, K, w, e, T0] or [γ, K, w, e, T0, P] where
                 derivative is to be calculated.
        P - period of orbit (if known). Must be supplied if only [γ, K, w, e, T0]
            are passed in as parameters.
        ang - list of true anomalies given parameters [γ, K, w, e, T0, P]
        """
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v = ang[0]
        
        temp = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V) *weight
        temp *= (np.sin(v+w)+e*np.sin(w))
        return sum((-2*K*temp))
    
    def parte(params, ang, P = None):
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
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v, E = ang
        
        dEde = np.sin(E)/(1-e*np.cos(E))
        dvde = np.sqrt((e+1)/(1-e)) * (0.5) * (np.cos(E/2)**-2)
        dvde *= dEde
        dvde += (1-e)**-2 * np.sqrt((1-e)/(1+e)) * np.tan(E/2)
        dvde *= 2 * np.cos(v/2)**2
        temp = K*(-np.sin(v + w)*dvde + np.cos(w))
        temp *= (V0 + (np.cos(v + w) + e*np.cos(w))*K - V) *weight
        return sum(2*temp)
    
    def partT0(params, ang, P = None):
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
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v, E = ang
        
        temp = ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight
        temp *= K * np.sqrt((1+e)/(1-e)) * np.sin(v + w)
        temp *= (1+ np.cos(v)) / (1 + np.cos(E))
        temp *= 2*np.pi / (P*(1-e*np.cos(E)))
        return sum(2*temp)
    
    def partP(params, ang, P = None):
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
        
        if len(params) == 6:
            V0, K, w, e, T0, P = params
        else:
            V0, K, w, e, T0 = params
            if P == None:
                print('P must be supplied')
                sys.exit()
        v, E = ang
        
        dEdP = ((-2*np.pi)/(P**2)) * (time - T0) * (1-e*np.cos(E))**(-1) 
        dvdP = np.sqrt((1+e)/(1-e)) * ((1+ np.cos(v)) / (1 + np.cos(E))) * dEdP
        temp = -K * np.sin(v + w) * dvdP
        temp *= ((V0 + ((np.cos(v + w) + e*np.cos(w))*K)) - V)*weight
        return sum(2*temp)
    
    def PrimarySolve(period = None, Pguess = None,covariance = False, graphs = True):
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
        
        # """
        if Pguess != None:
            period = Pguess
        x = x0(period) # initial guess
        if x == None:
            print('User must provide an estimate for P.')
            sys.exit()
            
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
    
        l = 3    
        while True: # perfoming minimization
            xLast = x
            ang = true(time, T0, e, P, ec = True)
            
            pVK = derivative(partK, x, P, ang, 'V0')
            pVw = derivative(partw, x, P, ang, 'V0')
            pVe = derivative(parte, x, P, ang, 'V0')
            pVT0 = derivative(partT0, x, P, ang, 'V0')
            pKw = derivative(partw, x, P, ang, 'K')
            pKe = derivative(parte, x, P, ang, 'K')
            pKT0 = derivative(partT0, x, P, ang, 'K')
            pwe = derivative(parte, x, P, ang, 'w')
            pwT0 = derivative(partT0, x, P, ang, 'w')
            peT0 = derivative(partT0, x, P, ang, 'e')
            
            if period == None:
                
                # Building Hessian matrix
                
                pVP = derivative(partV0, x, P, ang, 'P')
                pKP = derivative(partK, x, P, ang, 'P')
                pwP = derivative(partw, x, P, ang, 'P')
                peP = derivative(parte, x, P, ang, 'P')
                pT0P = derivative(partT0, x, P, ang, 'P')
                
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
                H1 = [derivative(partV0, x, P, ang, "V0") * (1+l), pVK, pVw, pVe, pVT0]
                H2 = [pVK, derivative(partK, x, P, ang, 'K') * (1+l), pKw, pKe, pKT0]
                H3 = [pVw, pKw, derivative(partw, x, P, ang, 'w') * (1+l), pwe, pwT0]
                H4 = [pVe, pKe, pwe, derivative(parte, x, P, ang, 'e') * (l+1), peT0]
                H5 = [pVT0, pKT0, pwT0, peT0, derivative(partT0, x, P, ang, 'T0') * (1+l)]
                a = 0.5 * np.array([H1,H2,H3,H4,H5])
                b = -0.5 * np.array([partV0(x, ang, P), partK(x, ang, P), partw(x, ang, P), 
                                     parte(x, ang, P), partT0(x, ang, P)])
            
            x = xLast + np.linalg.solve(a, b) # Next iterations parameter values
            
            if period == None:
                V0, K, w, e, T0, P = x
            else:
                V0, K, w, e, T0 = x
           
            
            if SumSquared(x, P) < SumSquared(xLast, P): # if params got better
                if abs(SumSquared(xLast, P) - SumSquared(x, P)) < 0.01 or abs(1-(SumSquared(x, P)/SumSquared(xLast, P))) < 1e-3:
                    break # convergence acheived! (probably)
                else:
                    l /= 9 # lowering damping parameter
                    l = max(l, 10e-7)
            else: # if params got worse
                l *= 11 # raising damping parameter
                l = min(l, 10e7)
        
        if x[2] < -np.pi or x[2] >= np.pi:
            while x[2] < -np.pi:
                x[2] = x[2] + 2*np.pi
            while x[2] >= np.pi:
                x[2] = x[2] - 2 * np.pi
        
        # other values to be returned
        nsec = 2*np.pi / (24*3600*P)
        a1 = (x[1]*np.sqrt(1-x[3]**2) / nsec) * 10e3 / 10e9
        fm = 3.985e-20 * (a1 * 10e9 / 10e3)**3 / P**2
        
        if period == None:
            x = [x[0], x[1], x[2], x[3], x[4], x[5], a1, fm]
        else:
            x = [x[0], x[1], x[2], x[3], x[4], P, a1, fm]
        
        C = np.linalg.inv(a) # covariance matrix
        err = []
        
        for q in range(var):
            err.append(np.sqrt(C[q][q]))
        if period != None:
            err.append(0)
        
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
        if graphs:
    
            # RV vs time data
            times = np.linspace(min(time), max(time), N) # times for fit
            v = true(times, x[4], x[3], P) # array of true anomalies at every time for fit
            fit = x[0] + x[1]*(np.cos(v + x[2])+ x[3]*np.cos(x[2]))
            
            # RV vs phase data
            phase = ((time - x[4]) / P) % 1 # phase for every RV data point
            phases = ((times - x[4]) / P) % 1 # phase for every point on fit
            
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            for i in range(len(phases)):
                phaseV.append([phase[i], V[i]])
                phasesfit.append([phases[i], fit[i]])
            phaseV = np.array(sorted(phaseV,key=lambda l:l[0]))
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))
            
            # average velocity, γ
            gam = np.zeros(len(phases)) + x[0] # average velocity
            
            # time plot
            plt.figure()
            plt.scatter(time, V, s = 18)
            plt.plot(times, fit, color = "orange")
            plt.xlabel("RJD")
            plt.ylabel("Velocity (km s$^{-1}$)")
            plt.plot(times, gam, ls = "--", color = "black")
            
            # phase plot
            plt.figure()
            plt.scatter(phaseV[:,0], phaseV[:,1], s = 18)
            plt.plot(phasesfit[:,0], phasesfit[:,1], color = "orange")
            plt.xlabel("Phase")
            plt.ylabel("Velocity (km s$^{-1}$)")
            plt.plot(phases, gam, ls = "--", color = "black")
        
        # rounding
        V0 = sigFig(x[0], 5)
        K = sigFig(x[1], 5)
        w = sigFig((x[2] * 180/np.pi)%360, 6)
        e = sigFig(x[3], 5)
        T0 = sigFig(x[4], 6)     
        a1 = sigFig(x[6], 5)
        fm = sigFig(x[7], 5)
            
        V0err = sigFig(err[0], 5)
        Kerr = sigFig(err[1], 5)
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
            return x, err, C
        else:   
            return x, err
        
    def CompanionSolve(X, err = [None, None, None, None, None, None], shift = False, graphs = True):
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
        
        K = (max(V) - min(V))/2
        
        params = V0, K, w, e, T0, P
        ang = true(time, T0, e, P, ec = True)
        
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
        fm = sigFig(fm, 5)
        
        x = [V0, K, w, e, T0, P, a1, fm]
        
        if graphs:
    
            # RV vs time data
            times = np.arange(min(time), max(time)+1, 10)
            v = true(times, x[4], x[3], P) # array of true anomalies at every time for fit
            fit = x[0] + x[1]*(np.cos(v + x[2])+ x[3]*np.cos(x[2]))
            
            # RV vs phase data
            phase = ((time - x[4]) / P) % 1 # phase for every RV data point
            phases = ((times - x[4]) / P) % 1 # phase for every point on fit
            
            
            # rearranging the data to fit as a function of phase
            phaseV = []
            phasesfit = []
            for i in range(N):
                phaseV.append([phase[i], V[i]])
            for i in range(len(phases)):
                phasesfit.append([phases[i], fit[i]])
            phaseV = np.array(sorted(phaseV,key=lambda l:l[0]))
            phasesfit = np.array(sorted(phasesfit,key=lambda l:l[0]))
            
            # average velocity, γ
            gam = np.zeros(len(phases)) + x[0] # average velocity
            
            # time plot
            plt.figure()
            plt.scatter(time,V, s = 18)
            plt.plot(times,fit, color = "orange")
            plt.xlabel("RJD")
            plt.ylabel("Velocity (km s$^{-1}$)")
            plt.plot(times,gam, ls = "--", color = "black")
            
            # phase plot
            plt.figure()
            plt.scatter(phaseV[:,0],phaseV[:,1], s = 18)
            plt.plot(phasesfit[:,0] ,phasesfit[:,1], color = "orange")
            plt.xlabel("Phase")
            plt.ylabel("Velocity (km s$^{-1}$)")
            plt.plot(phases,gam,ls = "--", color = "black")
        
        if err[0] != None: # if user inputted errors   
            err[1] += 0.01 # adding error from newtons method onto K error
            
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
            werr = sigFig((err[2] * 180/np.pi)%360, 6)
            eerr = sigFig(err[3], 5)
            T0err = sigFig(err[4], 6)
            a1err = sigFig(err[6], 5)
            fmerr = sigFig(err[7], 5)
            
            if shift:
                return x, [V0err,Kerr,werr,eerr,T0err,a1err,fmerr], dV0
            else:
                return x, [V0err,Kerr,werr,eerr,T0err,a1err,fmerr]
        
        else:
            x[2] = sigFig(wdeg, 6)
            if shift:
                return x, dV0
            else:
                return x
         
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
    
        while True: # single variable Newton's method
            paramsLast = params
            V0 = params[0]
            V0 = V0 - partV0(params, ang, P)/(derivative(partV0, params, P, ang, 'V0'))
            params = V0, K, w, e, T0, P
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
        origSets = np.empty([len(setNum),4,len(sets)],float) # 3D array
        origSets[:,:,:] = np.nan
        # seperating the sets
        for i in range(len(setNum)):
            for j in range(len(sets)):
                if setNum[i] == sets[j]:
                    origSets[i,:,j] = data_sorted[i,:]
        mask = np.isfinite(origSets) # mask for which elements have a value
        return (origSets,mask)
    
    if data.shape[1] >= 4: # if γ correcting is required
        splitSets, splitMask = split() # splitting sets
        contenders = [] # holds if a data set could be the best
        
    
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
        X,err = PrimarySolve(period = Period, covariance = False, graphs = False)
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
                    
        time = data_sorted[:,0]
        V = data_sorted[:,1] 
        N = len(time) 
        weight = data_sorted[:,2] 
        
        
    if star == "primary":
        if covariance:
            x,err,C = PrimarySolve(Period,Pguess,covariance,graphs)
            return x,err,C
        else:
            x,err = PrimarySolve(Period,Pguess,covariance,graphs)
            return x,err
    
    elif star == "secondary":
        if shift:
            if err != None:
                x,err,dV0 =CompanionSolve(X,err,shift,graphs)
                return x,err,dV0
            else:
                x,dV0 =CompanionSolve(X,err,shift,graphs)
                return x,dV0
        else:
            if err != None:
                x,err =CompanionSolve(X,err,shift,graphs)
                return x,err
            else:
                x =CompanionSolve(X,err,shift,graphs)
                return x
    else:
        print("Please let the parameter, star, be the string \"primary\" or \"secondary\"")
        sys.exit()
