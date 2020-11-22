import numpy as np
import pandas as pd
from scipy import stats


def failurerate(rt, var, conflvlvar, conflvltest):
    ret = -rt
    p = 1-conflvlvar
    T = ret.shape[0]
    ret[rt > var[:, 0]] = 1
    ret[ret != 1] = 0
    N = ret.sum()
    failrate = N/T

    if failrate !=0:
        LR = 2*((N*np.log(failrate)+np.log(1-failrate)*(T-N))-(N*np.log(p))+np.log(1-p)*(T-N))
    else:
        LR=0

    criticalvalue = stats.chi2.ppf(conflvltest, df=1)

    if LR > criticalvalue:
        resultat = 1
    else:
        resultat = 0

    ff = 0
    fs = 0
    sf = 0
    ss = 0

    for i in range(1, T):
        if ret[i] == 1 and ret[i-1] == 1:
            ff = ff + 1
        elif ret[i] == 1 and ret[i-1] == 0:
            fs = fs + 1
        elif ret[i] == 0 and ret[i-1] == 1:
            sf = sf + 1
        elif ret[i] == 0 and ret[i-1] == 0:
            ss = ss + 1

    if failrate == 0:
        LRCCI = 0
    else:
        pi0 = sf/(ss + sf)
        pi1 = ff/(fs + ff)
        pi2 = (sf + ff)/(ss + ff + sf + fs)
        upper = ((1-pi2)**(ss+fs))*(pi2**(sf+ff))
        lower = ((1 - pi0)**ss) * (pi0**sf) * ((1 - pi1)**fs) * (pi1**ff)
        LRCCI = -2*np.log(upper/lower)

    CC = LR + LRCCI
    result = np.array([failrate, LR, resultat, LRCCI, CC])

    return result


def test(rt, var):
    ret = -rt
    T = ret.shape[0]
    ret[rt > var[:, 0]] = 1
    ret[ret != 1] = 0
    N = ret.sum()
    failrate = N/T

    return failrate

