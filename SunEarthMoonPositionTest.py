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


def get_data():
    sunData = open("./sun_ephemeris2023.txt")
    earthData = open("./earth_ephemeris2023.txt")
    moonData = open("./moon_ephemeris2023.txt")

    sunTxt = sunData.read()
    earthTxt = earthData.read()
    moonTxt = moonData.read()

    # strip out newline characters and the file header and footer
    sunTxt = re.sub(r'\n', '', sunTxt)
    sunTxt = re.sub(r'[*]+[\s\S]*SOE', '', sunTxt)
    sunTxt = re.sub(r'EOE[\s\S]*[*]+', '', sunTxt)

    earthTxt = re.sub(r'\n', '', earthTxt)
    earthTxt = re.sub(r'[*]+[\s\S]*SOE', '', earthTxt)
    earthTxt = re.sub(r'EOE[\s\S]*[*]+', '', earthTxt)

    moonTxt = re.sub(r'\n', '', moonTxt)
    moonTxt = re.sub(r'[*]+[\s\S]*SOE', '', moonTxt)
    moonTxt = re.sub(r'EOE[\s\S]*[*]+', '', moonTxt)

    # sunPosList = re.split(r"[0-9]{7}.[0-9]{9} = A.D.", sunTxt)[1:]
    fieldPattern = r'=(\s?\S+)'
    sunXvals = re.findall(r'X ' + fieldPattern, sunTxt)
    sunYvals = re.findall(r'Y ' + fieldPattern, sunTxt)
    sunZvals = re.findall(r'Z ' + fieldPattern, sunTxt)
    sunVXvals = re.findall(r'VX' + fieldPattern, sunTxt)
    sunVYvals = re.findall(r'VY' + fieldPattern, sunTxt)
    sunVZvals = re.findall(r'VZ' + fieldPattern, sunTxt)

    earthXvals = re.findall(r'X ' + fieldPattern, earthTxt)
    earthYvals = re.findall(r'Y ' + fieldPattern, earthTxt)
    earthZvals = re.findall(r'Z ' + fieldPattern, earthTxt)
    earthVXvals = re.findall(r'VX' + fieldPattern, earthTxt)
    earthVYvals = re.findall(r'VY' + fieldPattern, earthTxt)
    earthVZvals = re.findall(r'VZ' + fieldPattern, earthTxt)

    moonXvals = re.findall(r'X ' + fieldPattern, moonTxt)
    moonYvals = re.findall(r'Y ' + fieldPattern, moonTxt)
    moonZvals = re.findall(r'Z ' + fieldPattern, moonTxt)
    moonVXvals = re.findall(r'VX' + fieldPattern, moonTxt)
    moonVYvals = re.findall(r'VY' + fieldPattern, moonTxt)
    moonVZvals = re.findall(r'VZ' + fieldPattern, moonTxt)

    # 3D list (365*3*3) that stores (sun, earth, moon) position data for each day of the year
    dailyPositionsArray = [[[sunXvals[i], sunYvals[i], sunZvals[i], sunVXvals[i], sunVYvals[i], sunVZvals[i]], [earthXvals[i], earthYvals[i], earthZvals[i], earthVXvals[i], earthVYvals[i], earthVZvals[i]], [moonXvals[i], moonYvals[i], moonZvals[i], moonVXvals[i], moonVYvals[i], moonVZvals[i]]] for i in range(len(sunXvals))]

    return dailyPositionsArray


def processData(day_num, dailyPositionsArray):
    curSunData, curEarthData, curMoonData = dailyPositionsArray[day_num]

    X_1 = (float(curSunData[0])) * AU
    X_2 = (float(curEarthData[0])) * AU
    X_3 = (float(curMoonData[0])) * AU
    Y_1 = (float(curSunData[1])) * AU
    Y_2 = (float(curEarthData[1])) * AU
    Y_3 = (float(curMoonData[1])) * AU
    Z_1 = (float(curSunData[2])) * AU
    Z_2 = (float(curEarthData[2])) * AU
    Z_3 = (float(curMoonData[2])) * AU
    Vx_1 = (float(curSunData[3]) * AU) / DAY_SEC
    Vx_2 = (float(curEarthData[3]) * AU) / DAY_SEC
    Vx_3 = (float(curMoonData[3]) * AU) / DAY_SEC
    Vy_1 = (float(curSunData[4]) * AU) / DAY_SEC
    Vy_2 = (float(curEarthData[4]) * AU) / DAY_SEC
    Vy_3 = (float(curMoonData[4]) * AU) / DAY_SEC
    Vz_1 = (float(curSunData[5]) * AU) / DAY_SEC
    Vz_2 = (float(curEarthData[5]) * AU) / DAY_SEC
    Vz_3 = (float(curMoonData[5]) * AU) / DAY_SEC

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
    pi[0, :] = (pi[0, :] - V_CM) * M_SUN
    pi[1, :] = (pi[1, :] - V_CM) * M_EARTH
    pi[2, :] = (pi[2, :] - V_CM) * M_MOON
    yi = [qi, pi]
    yn = np.array(yi)
    yn[:] = yi[:]

    return yn


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


# Returns an array consisting of Sun, Earth, and Moon coordinates at Nt times throughout the year
def get_cur_positions(day_num, dailyPositionsArray):

    yn = processData(day_num, dailyPositionsArray)
    return np.array([yn[0, 0], yn[0, 1], yn[0, 2]])


# print(get_cur_positions(59945.5))




