#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-body Sun-Earth-Moon evolution with Symplectic integrator
output is a list of Cartesian coordinates for each body at each of Nt
times, starting at J2000
@author: jschnitt
"""
import numpy as np
import re
from CalendarUtilities import mjd, isLeap

# import matplotlib.pyplot as plt


def processData(year):

    # open and read in ephemeris files as strings
    sunData = open("./sun_ephemeris.txt")
    earthData = open("./earth_ephemeris.txt")
    moonData = open("./moon_ephemeris.txt")

    sunTxt = sunData.read()
    earthTxt = earthData.read()
    moonTxt = moonData.read()

    # extract sun/moon/earth positions for the current year
    startPat = str(year) + "-Jan-01"
    endPat = str(year+1) + "-Jan-01"
    cur_yr_regex = fr"(?:SOE[\s\S]*)({startPat}[\s\S]*X=[\s\S]*)(?:{endPat})"

    curSunData = re.search(cur_yr_regex, sunTxt).groups()[0]
    curEarthData = re.search(cur_yr_regex, earthTxt).groups()[0]
    curMoonData = re.search(cur_yr_regex, moonTxt).groups()[0]

    # for the current year's Earth data:
    curSunData = re.sub(r'\n', '', curSunData)  # get rid of newline characters
    curSunDataList = re.split(r' *[A-Z]{1,2} ?= ?', curSunData)[1:]  # split into fields

    # for the current year's Earth data:
    curEarthData = re.sub(r'\n', '', curEarthData)  # get rid of newline characters
    curEarthDataList = re.split(r' *[A-Z]{1,2} ?= ?', curEarthData)[1:]  # split into fields

    # for the current year's Moon data:
    curMoonData = re.sub(r'\n', '', curMoonData)  # get rid of newline characters
    curMoonDataList = re.split(r' *[A-Z]{1,2} ?= ?', curMoonData)[1:]  # split into fields

    return assignVars(curSunDataList, curEarthDataList, curMoonDataList)



def assignVars(curSunDataList, curEarthDataList, curMoonDataList):

    X_1 = (float(curSunDataList[0])) * AU  # -7.139867342351965E-03*AU = -106810895149.63219
    X_2 = (float(curEarthDataList[0])) * AU
    X_3 = (float(curMoonDataList[0])) * AU
    Y_1 = (float(curSunDataList[1])) * AU
    Y_2 = (float(curEarthDataList[1])) * AU
    Y_3 = (float(curMoonDataList[1])) * AU
    Z_1 = (float(curSunDataList[2])) * AU
    Z_2 = (float(curEarthDataList[2])) * AU
    Z_3 = (float(curMoonDataList[2])) * AU
    Vx_1 = (float(curSunDataList[3]) * AU) / DAY_SEC
    Vx_2 = (float(curEarthDataList[3]) * AU) / DAY_SEC
    Vx_3 = (float(curMoonDataList[3]) * AU) / DAY_SEC
    Vy_1 = (float(curSunDataList[4]) * AU) / DAY_SEC
    Vy_2 = (float(curEarthDataList[4]) * AU) / DAY_SEC
    Vy_3 = (float(curMoonDataList[4]) * AU) / DAY_SEC
    Vz_1 = (float(curSunDataList[5]) * AU) / DAY_SEC
    Vz_2 = (float(curEarthDataList[5]) * AU) / DAY_SEC
    Vz_3 = (float(curMoonDataList[5]) * AU) / DAY_SEC

    V_CM = np.zeros(3)
    Vx_CM = (Vx_1 * M_SUN + Vx_2 * M_EARTH + Vx_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    Vy_CM = (Vy_1 * M_SUN + Vy_2 * M_EARTH + Vy_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    Vz_CM = (Vz_1 * M_SUN + Vz_2 * M_EARTH + Vz_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    V_CM[:] = (Vx_CM, Vy_CM, Vz_CM)
    P_CM = np.zeros(3)
    X_CM = (X_1 * M_SUN + X_2 * M_EARTH + X_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    Y_CM = (Y_1 * M_SUN + Y_2 * M_EARTH + Y_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    Z_CM = (Z_1 * M_SUN + Z_2 * M_EARTH + Z_3 * M_MOON) / (M_MOON + M_EARTH + M_SUN)
    P_CM[:] = (X_CM, Y_CM, Z_CM)

    qi = np.zeros((Nbody, 3))
    pi = np.zeros((Nbody, 3))
    qi[0, :] = (X_1, Y_1, Z_1)
    qi[1, :] = (X_2, Y_2, Z_2)
    qi[2, :] = (X_3, Y_3, Z_3)
    qi[0, :] = (qi[0, :] - P_CM)
    qi[1, :] = (qi[1, :] - P_CM)
    qi[2, :] = (qi[2, :] - P_CM)
    pi[0, :] = (Vx_1, Vy_1, Vz_1)
    pi[1, :] = (Vx_2, Vy_2, Vz_2)
    pi[2, :] = (Vx_3, Vy_3, Vz_3)
    pi[0, :] = (pi[0, :] - V_CM)
    pi[1, :] = (pi[1, :] - V_CM)
    pi[2, :] = (pi[2, :] - V_CM)
    pi[0, :] = (pi[0, :]) * M_SUN
    pi[1, :] = (pi[1, :]) * M_EARTH
    pi[2, :] = (pi[2, :]) * M_MOON
    yi = [qi, pi]
    yn = np.array(yi)
    yn[:] = yi[:]

    return yn


def DT(y):
    p_ = y[1]
    dx = np.zeros((Nbody, 3))
    for i in range(Nbody):
        dx[i] = p_[i] / Mi[i]
    return dx


def DV(y):
    x_ = y[0]
    rij = np.zeros((Nbody, Nbody))
    rij_ = np.zeros((Nbody, Nbody, 3))
    Fi_ = np.zeros((Nbody, 3))
    for i in range(Nbody):
        for j in range(Nbody):
            rij_[i, j] = x_[i] - x_[j]
            rij[i, j] = np.linalg.norm(rij_[i, j])
            if (i != j):
                Fi_[i] += -G * Mi[i] * Mi[j] * rij_[i, j] / rij[i, j] ** 3.
    dp = Fi_
    return dp


def advanceSy4(y, dt):
    c1 = 1. / (2. * (2. - 2. ** (1. / 3.)))
    c2 = c1 * (1 - 2. ** (1. / 3.))
    c3 = c2
    c4 = c1
    d1 = 1. / (2. - 2. ** (1. / 3.))
    d2 = -d1 * 2. ** (1. / 3.)
    d3 = d1
    d4 = 0.
    yn = np.array(y)
    dx = DT(yn)
    yn[0] = yn[0] + dt * c1 * dx
    dp = DV(yn)
    yn[1] = yn[1] + dt * d1 * dp
    dx = DT(yn)
    yn[0] = yn[0] + dt * c2 * dx
    dp = DV(yn)
    yn[1] = yn[1] + dt * d2 * dp
    dx = DT(yn)
    yn[0] = yn[0] + dt * c3 * dx
    dp = DV(yn)
    yn[1] = yn[1] + dt * d3 * dp
    dx = DT(yn)
    yn[0] = yn[0] + dt * c4 * dx
    dp = DV(yn)
    yn[1] = yn[1] + dt * d4 * dp
    return (yn)


# T0 = 51544.5  # Epoch for J2000
Nbody = 3
G = 6.674e-8
AU = 14959787070000.0   # cm
DAY_SEC = 86400.0
M_EARTH = 5.97243695e27
M_SUN = 1.988499e33
M_MOON = 7.34922026e25
M_MOON_EARTH = M_MOON + M_EARTH
GM_SUN = 132712440041.93938e15
GM_EARTH = 398600.435436e15
GM_MOON = 4902.800066e15
GMi = np.zeros(Nbody)
Mi = np.zeros(Nbody)
GMi[:] = (GM_SUN, GM_EARTH, GM_MOON)
Mi[:] = (M_SUN, M_EARTH, M_MOON)
Mtotal = np.sum(Mi)


# Returns an array consisting of Sun, Earth, and Moon x, y, z coordinates at Nt times throughout the year
# (each element has shape 4 elements where the last three are themselves each a list of 3)
def get_positions(year):

    T0 = mjd(1, 1, year, (0, 0))  # MJD of the first day of the given year
    tstart = T0
    if isLeap(year): yr_len = 366
    else:            yr_len = 365
    tstop = T0 + yr_len  # days
    Nt = 1001
    tt = np.linspace(tstart, tstop, Nt)
    dt = (tt[1] - tt[0]) * DAY_SEC

    yn = processData(year)
    yi_t = np.zeros((Nt, 2, Nbody, 3))
    yi_t[0] = yn

    positionsArr = []
    for it in range(Nt - 1):

        yn = advanceSy4(yn, dt)
        yi_t[it + 1] = yn
        positionsArr.append([tt[it], yi_t[it, 0, 0], yi_t[it, 0, 1], yi_t[it, 0, 2]])

    return positionsArr


# plt.plot(yi_t[:,0,2,0]-yi_t[:,0,1,0],yi_t[:,0,2,1]-yi_t[:,0,1,1])
# plt.gca().set_aspect('equal')
# plt.show()
# outputFile.close()
get_positions(2023)




