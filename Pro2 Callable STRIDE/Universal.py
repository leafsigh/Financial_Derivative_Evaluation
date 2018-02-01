import numpy as np
import pandas as pd
from scipy import stats
import math


def b_s(cp, s, k, t, v, rf, div):
    """ Price an option using the Black-Scholes model.
    s: initial stock price
    k: strike price
    t: expiration time
    v: volatility
    rf: risk-free rate
    div: dividend
    cp: +1/-1 for call/put
    """

    d1 = (math.log(s / k) + (rf - div + 0.5 * math.pow(v, 2)) * t)/(v * math.sqrt(t))
    d2 = d1 - v * math.sqrt(t)

    optprice = (cp * s * math.exp(- div * t) * stats.norm.cdf(cp * d1))\
               - (cp * k * math.exp(- rf * t) * stats.norm.cdf(cp * d2))
    return optprice

vtrigger = np.zeros((60001, 60001))

principal = 1000
jmax_max = 230
jmax_min = 40
kd = pd.read_csv('.\\read\\keydates.csv', index_col='Date')
yearfrac = kd.yearfrac.tolist()
T = yearfrac[-1]
cpndates = [0, ] + yearfrac[4:]
r = 0.858903322 / 100
sigma = 30.331754 / 100
S0 = 116.62
cpn = 5.833333
barrier = .65 * S0

A = [0, ] * 5000
B = [0, ] * 5000
C = [0, ] * 5000


def calculate(Smax, Smin, jmax, imax):
    dS = 1.0 * (Smax - Smin) / jmax
    ST = np.arange(Smin, Smax + dS, dS)
    dt = T / imax
    time = np.arange(0, T + dt, dt)

    for j in range(1, jmax):
        A[j] = (.5 * sigma ** 2 * j ** 2 + .5 * r * j) * dt
        B[j] = 1 - r * dt - sigma ** 2 * j ** 2 * dt
        C[j] = (.5 * sigma ** 2 * j ** 2 - .5 * r * j) * dt

    for j in range(0, jmax + 1):
        vtrigger[imax][j] = min(principal, ST[j] / S0 * principal) + cpn

# Set upper/lower boundaries.
    cd_sign = 11
    for i in range(imax - 1, -1, -1):
        vtrigger[i][jmax] = vtrigger[i + 1][jmax] * np.exp(-r * dt)
        vtrigger[i][0] = vtrigger[i + 1][0] * np.exp(-r * dt)
        if time[i] < cpndates[cd_sign]:
            vtrigger[i][jmax] += np.exp(-r * (cpndates[cd_sign] - time[i])) * cpn
            vtrigger[i][0] += np.exp(-r * (cpndates[cd_sign] - time[i])) * cpn
            if cd_sign == 9 or cd_sign == 6:
                vtrigger[i][jmax] = min(np.exp(-r * (cpndates[cd_sign] - time[i])) * (principal + cpn),
                                        vtrigger[i][jmax])
                vtrigger[i][0] = min(np.exp(-r * (cpndates[cd_sign] - time[i])) * (principal + cpn),
                                     vtrigger[i][0])
            cd_sign -= 1

    cd_sign = 11
    for i in range(imax - 1, -1, -1):
        for j in range(1, jmax):
            vtrigger[i][j] = A[j] * vtrigger[i + 1][j + 1] + \
                             B[j] * vtrigger[i + 1][j] + \
                             C[j] * vtrigger[i + 1][j - 1]
            if time[i] < cpndates[cd_sign]:
                vtrigger[i][j] += np.exp(-r * (cpndates[cd_sign] - time[i])) * cpn
                if cd_sign == 9 or cd_sign == 6:
                    vtrigger[i][j] = min(np.exp(-r * (cpndates[cd_sign] - time[i])) * (principal + cpn),
                                         vtrigger[i][j])
                cd_sign -= 1

    b = jmax
    p = jmax
    while ST[b] > barrier:
        b -= 1
    while ST[p] > S0:
        p -= 1

    # trigger BS
    vtbs = b_s(-1, S0, S0, T, sigma, r, 0)
    vcb = 0
    for i in range(1, 13):
        vcb += cpn / (1 + r) ** (cpndates[i])
    vcb += 1000 / (1 + r) ** T
    vtbs += vcb
    vt = vtrigger[0][p]


# considering barriers
    for j in range(b, jmax + 1):
        vtrigger[imax][j] = principal + cpn
    for j in range(b + 1, jmax):
        A[j] = (.5 * sigma ** 2 * (j - b) ** 2 + .5 * r * (j - b)) * dt
        B[j] = 1 - r * dt - sigma ** 2 * (j - b) ** 2 * dt
        C[j] = (.5 * sigma ** 2 * (j - b) ** 2 - .5 * r * (j - b)) * dt
    cd_sign = 11
    for i in range(imax - 1, -1, -1):
        for j in range(b + 1, jmax):
            vtrigger[i][j] = A[j] * vtrigger[i + 1][j + 1] + \
                             B[j] * vtrigger[i + 1][j] + \
                             C[j] * vtrigger[i + 1][j - 1]
            if time[i] < cpndates[cd_sign]:
                vtrigger[i][j] += np.exp(-r * (cpndates[cd_sign] - time[i])) * cpn
                if cd_sign == 9 or cd_sign == 6:
                    vtrigger[i][j] = min(np.exp(-r * (cpndates[cd_sign] - time[i])) * (principal + cpn),
                                         vtrigger[i][j])
                cd_sign -= 1
    V = (vtrigger[0][p + 1] - vtrigger[0][p]) / dS * (S0 - ST[p]) + vtrigger[0][p]
    lmbd = (S0 - ST[p]) / dS
    return V, lmbd, vt, vtbs

imax = 20000
for jmax in range(40, 180, 20):
    print jmax, calculate(S0 * 3, 0, jmax, imax)
