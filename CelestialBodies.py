import random

import numpy as np
from sympy import *
from sympy import symbols, linsolve, solveset
import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate
import cv2
import glob
import re
from queue import Queue
from DateConversion import year_lists, molad_determination, LEAP_YEARS, hebrew_to_greg, year_of_cycle, length_hebrew_yr, hebrew_year_type
from SunEarthMoonPosition import get_positions
from PIL import Image
from CalendarUtilities import month_len, isLeap, mjd
from SunEarthMoonPositionTest import get_data, get_cur_positions

# user inputs
au_cm = 14959787070000.0
# earth_rad = 637813700   # equatorial radius, in cm
earth_rad = 637100880   # arithmetic mean radius, in cm
moon_rad = 173740000    # volumetric mean radius, in cm
sun_rad = 69570000000   # volumetric mean radius/nominal solar radius, in cm
celestial_rad = au_cm   # 1 AU, in cm
degrees_day = 360.98561     # degrees that Earth rotates in a day
equinox_time = (15, 33)     # GMT
mjd_equinox = 59658.64791666667
mjd_j2000 = 51544
epoch_Rambam = 30, 3, 1178  # Rambam's epoch (Wednesday night, 3 Nisan, 4938 / March 30, 1178)
mjd_Rambam_epoch = mjd(epoch_Rambam[0], epoch_Rambam[1], epoch_Rambam[2], (0, 0))  # epoch mjd
phi_jer = 35.214        # degrees
phi_ny = -74.006        # degrees
phi_texas = -97.7431
# theta_texas = 59.7328
theta_texas = 57.3
theta_jer = 58.232      # degrees
theta_ny = 49.3         # degrees
dayNum_mjdEpoch = 2052005
sunset_angle = 90.8     # degrees
obliquity = 23.43       # degrees
moon_obliquity = 1.5424 # degrees
T_side = 27.321661         # sidereal period of the moon, in days
rotation_phase0 = 217.88   # moon's rotation phase at Epoch 0, in degrees
avg_month_len = 29 + (12 / 24) + ((793 / 1080) / 24)
days_per_century = 36524.2199 # on average
# eccentricity = 0.0167 * 2
eccentricity = 0.0347
# eccentricity = 0.033

## 3/10/24:
psi_cam = 86
phi_cam = 262
## 6/6/24 from 7:30-9:30pm in ny:
# psi_cam = 86
# phi_cam = 300
canvas_size = 256

psi_cam = 90
phi_cam = 270

# convert to radians
radians_day = degrees_day * np.pi/180
phi_jer = phi_jer * np.pi/180
phi_ny = phi_ny * np.pi/180
phi_texas = phi_texas * np.pi/180
theta_jer = theta_jer * np.pi/180
theta_ny = theta_ny * np.pi/180
theta_texas = theta_texas * np.pi/180
sunset_angle = sunset_angle * np.pi/180
obliquity = obliquity * np.pi/180
moon_obliquity = moon_obliquity * np.pi/180
rotation_phase0 = rotation_phase0 * np.pi/180
psi_cam = psi_cam * np.pi/180
phi_cam = phi_cam * np.pi/180

# mjds of rosh chodesh for year 2000
rosh_chodesh_MJD = [(51551, 'Shevat'),  (51580, 'Adar I'), (51581, 'Adar I'), (51610, 'Adar II'), (51611, 'Adar II'),
                    (51640, 'Nisan'), (51669, 'Iyar'), (51670, 'Iyar'), (51699, 'Sivan'), (51728, 'Tammuz'),
                    (51729, 'Tammuz'), (51758, 'Av'), (51787, 'Elul'), (51788, 'Elul'), (51817, 'Tishrei'),
                    (51818, 'Tishrei'), (51846, 'Cheshvan'), (51847, 'Cheshvan'), (51876, 'Kislev'), (51905, 'Teves')]

# # array containing positions of sun, earth, and moon as returned from SunEarthMoonPosition
# positionsArr = SunEarthMoonPosition.get_positions()


def interpolation(xList, yList, order=2):
    system = []
    a, b, c = symbols('a b c')

    # create system of equations using given 4 (x,y) pairs
    for j in range(order+1):
        equation = ((xList[j] ** 2) * a) + ((xList[j]) * b) + c - yList[j]
        system.append(equation)

    # list of coefficients of solution to linear system
    coeff = linsolve(system, [a, b, c])

    X, Y = symbols('X Y')

    # get the one tuple in coefficient FiniteSet (the polynomial)
    polynomial = 0
    for tup in coeff:
        polynomial = (tup[0] * X ** 2) + (tup[1] * X) + tup[2]  # - sunset_angle

    return polynomial


# returns the position of the body in question (sun/earth/moon) at the given mjd
def get_pos(mjd_cur, body, positionsArr):

    index_map = {'s': 1, 'e': 2, 'm': 3}
    i = index_map[body]

    t_start = positionsArr[0][0]  # mjd of Jan 1 of this year
    t_stop = positionsArr[-1][0]
    Nt = len(positionsArr) + 1

    # determine where the given MJD lies relative to the list of mjd/position data points
    daysRecorded = t_stop - t_start
    fractionOfTotal = (mjd_cur - t_start) / daysRecorded
    est_index = int(fractionOfTotal * Nt)  # the index of this MJD in Positions array

    # go to the estimated index and adjust if necessary
    data = positionsArr[est_index]

    while mjd_cur < data[0]:                       # if mjd is lower than current entry, move backward
        est_index -= 1
        data = positionsArr[est_index]

    while mjd_cur > positionsArr[est_index+1][0]:  # if greater than the following entry, move forward
        est_index += 1
        data = positionsArr[est_index]

    # interpolate to get positions at exact time in question
    t_Vals = []
    pos_Vals = []
    for j in range(-1, 2):  # use the prev interval, this one, and the following one
        t_Vals.append(positionsArr[est_index+j][0] - t_start)
        pos_Vals.append(positionsArr[est_index+j][i])

    # isolate the variables
    x_Vals = [coords[0] for coords in pos_Vals]
    y_Vals = [coords[1] for coords in pos_Vals]
    z_Vals = [coords[2] for coords in pos_Vals]

    ## interpolate x, y, and z coordinates separately:
    # create polynomial for each
    p_X = interpolation(t_Vals, x_Vals, 2)
    p_Y = interpolation(t_Vals, y_Vals, 2)
    p_Z = interpolation(t_Vals, z_Vals, 2)

    # now, plug the given t into each polynomial; solve for x, y, and z
    X = symbols('X')
    x = p_X.subs(X, mjd_cur - t_start)
    y = p_Y.subs(X, mjd_cur - t_start)
    z = p_Z.subs(X, mjd_cur - t_start)
    # x = np.interp(mjd_cur-tstart, t_Vals, x_Vals)
    # y = np.interp(mjd_cur-tstart, t_Vals, y_Vals)
    # z = np.interp(mjd_cur-tstart, t_Vals, z_Vals)

    return np.array([x, y, z])


### Utility functions
def normalize(v):
    return v / length(v)


# returns the length of the given vector
def length(v):
    return np.sqrt(float(v.dot(v)))


# returns the euclidian distance between two 2D vectors
def distance(x1, x2, y1, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)


# returns the angular distance between two vectors
def angular_dist(v1, v2):
    return np.arccos(float((np.dot(v1, v2) / (length(v1) * length(v2)))))


# normalized position vector from current Jerusalem position to current position of sun/moon
def n_body_jer(n_jer, p_body):
    p_jer = n_jer * earth_rad
    p_body_jer = p_body - p_jer
    return normalize(p_body_jer)


def ecliptic_to_equatorial(v, mjd_cur, eps=obliquity):
    # x_eq = v[0]
    # y_eq = (np.cos(eps) * v[1]) - (np.sin(eps) * v[2])
    # z_eq = (np.sin(eps) * v[1]) + (np.cos(eps) * v[2])

    r_ec_eq = [[1, 0, 0],
               [0, np.cos(eps), -np.sin(eps)],
               [0, np.sin(eps), np.cos(eps)]]

    ## Add precession of earth's axis:
    cent_since_j2000 = (mjd_cur - mjd_j2000) / days_per_century  # centuries since j2000
    theta_pre = -(5028.796195 * cent_since_j2000) - (1.1054348 * cent_since_j2000 ** 2)

    # precession rotation matrix
    p_ec = [[np.cos(theta_pre), -np.sin(theta_pre), 0],
            [np.sin(theta_pre), np.cos(theta_pre), 0],
            [0, 0, 1]]
    p_eq = np.linalg.inv(p_ec)

    prod_inner = np.dot(p_eq, v)
    eq_p_coords = np.dot(r_ec_eq, prod_inner) * (np.pi / (180*3600))  # radians

    # return np.array([x_eq, y_eq, z_eq])
    return eq_p_coords


def equatorial_to_ecliptic(v, mjd_cur, eps=obliquity):
    x_ec = v[0]
    y_ec = (np.cos(eps) * v[1]) + (np.sin(eps) * v[2])
    z_ec = -(np.sin(eps) * v[1]) + (np.cos(eps) * v[2])

    # rotation matrix
    r_eq_ec = [[1, 0, 0],
               [0, np.cos(eps), np.sin(eps)],
               [0, -np.sin(eps), np.cos(eps)]]

    # matrix-vector multiplication
    ec_coords = np.dot(r_eq_ec, v)

    ## Add precession of earth's axis:
    cent_since_j2000 = (float(mjd_cur) - mjd_j2000) / days_per_century  # centuries since j2000
    theta_pre_arc_sec = -(5028.796195 * cent_since_j2000) - (1.1054348 * cent_since_j2000**2)  # cur precession angle
    theta_pre = theta_pre_arc_sec * (np.pi / (180*3600))    # in radians

    # precession rotation matrix
    p_ec = [[np.cos(theta_pre), -np.sin(theta_pre), 0],
            [np.sin(theta_pre), np.cos(theta_pre), 0],
            [0, 0, 1]]

    ec_p_coords = np.dot(p_ec, ec_coords)

    # return np.array([x_ec, y_ec, z_ec])
    return ec_p_coords



def phiDifAndTheta(time_zone):

    phi_dif = 0  # relative to greenwich
    theta = 0

    if time_zone == 'IST':
        phi_dif = phi_jer
        theta = theta_jer
    elif time_zone == 'EST':
        phi_dif = phi_ny
        theta = theta_ny
    elif time_zone == 'CST':
        phi_dif = phi_texas
        theta = theta_texas

    # elif time_zone == 'GMT':
    #     phi_dif = 0
    #     theta = 0

    return phi_dif, theta


def DST(time_zone, yr):
    dst_start_month = 3  # March

    if time_zone == 'IST':
        dst_start_day_of_week = (0, 'last')  # really starts on Fri before last Sunday

    elif time_zone == 'EST' or time_zone == 'CST':
        dst_start_day_of_week = (0, 2)  # starts on second Sunday in March

    # determine MJD of Friday before last Sunday in March of the given year

    # day of week of March 1
    greg_yrs, hebrew_yrs = year_lists()
    hebrew_yrs = None  # garbage collection

    yr_tuple = greg_yrs[yr + 3760]  # year tuple from year-lists specifying MJD and day of week of Jan 1

    # pointer starting at March 1
    if yr_tuple[3]: cur = 60  # if leap year
    else:           cur = 59  # if not leap year

    # number of Sundays in this March:
    day_of_week_Jan1 = yr_tuple[2]
    day_of_week_March1 = (day_of_week_Jan1 + cur) % 7

    if not day_of_week_March1:
        date_first_sunday = 1
    else:
        date_first_sunday = 8 - day_of_week_March1  # because the first of March is included

    if dst_start_day_of_week[1] == 'last':
        # number of full weeks between first Sunday and last day of month (inclusive)
        weeks = (month_len(dst_start_month, False) - date_first_sunday) // 7

        dst_start_monthDay = date_first_sunday + (7 * weeks)

        # go back to preceding Friday
        if time_zone == 'IST': dst_start_monthDay -= 2

    # if EST, looking for second sunday in March
    elif dst_start_day_of_week[1] == 2:
        dst_start_monthDay = date_first_sunday + 7

    # MJD of daylight-savings start date
    # = MJD of Jan 1 + remaining days in Jan and Feb + days in March until the Friday before its last Sunday
    DST_START = yr_tuple[1] + cur + dst_start_monthDay
    DST_START -= dayNum_mjdEpoch  # convert from anno mundi day count to MJD

    ## DETERMINATION OF DAYLIGHT SAVINGS TIME END DATE:
    # determine MJD of last Sunday in October:

    # days since start of DST (excluding those from October)
    DST_END = DST_START + (month_len(3, False) - dst_start_monthDay + 2)  # days remaining from March
    DST_END += (month_len(4, False) + month_len(5, False) + month_len(6, False) + month_len(7, False)
                + month_len(8, False) + month_len(9, False))  # days from April, May, June, July, August, September

    # date of last Sunday
    day_of_week_Oct1 = (5 + ((DST_END - DST_START + 1) % 7)) % 7  # Friday (day of DST start) + remaining days after complete weeks + 1 to get to Oct 1
    if not day_of_week_Oct1:
        date_first_sunday = 1
    else:
        date_first_sunday = 8 - day_of_week_Oct1  # because the first of Oct is included

    # number of full weeks between first Sunday and last day of month (inclusive)
    weeks = (month_len(10, False) - date_first_sunday) // 7

    date_last_sunday = date_first_sunday + (7 * weeks)
    DST_END += date_last_sunday

    return DST_START, DST_END


# time-difference relative to Greenwich
def time_dif(mjd_cur, time_zone, DST_START, DST_END):
    # DST_START, DST_END = DST(mjd_cur, time_zone, yr)

    if time_zone == 'IST':
        standard_time_dif = 2
        daylight_time_dif = 3

    elif time_zone == 'EST':
        standard_time_dif = -5
        daylight_time_dif = -4

    elif time_zone == 'CST':
        standard_time_dif = -6
        daylight_time_dif = -5

    elif time_zone == 'GMT':
        return 0

    # check if date falls out in DST:
    if DST_START <= mjd_cur < DST_END:
        return daylight_time_dif

    else:
        return standard_time_dif


# # modified julian date
# def mjd(d, m, y, t):
#     epoch = (17, 11, 1858, (0, 0))  # 17 November, 1858, 12 am (midnight)
#
#     days = 0  # days since epoch
#
#     if (d, m, y, t) == epoch:
#         return days  # 0
#
#     if y > epoch[2]:  # this year is past the epoch year
#         start = epoch
#         end = (d, m, y, t)
#
#     elif y < epoch[2]:  # this year precedes the epoch year
#         start = (d, m, y, t)
#         end = epoch
#
#     # complete days remaining from start month (not including start day)
#     days += (month_len(start[1], isLeap(start[2])) - epoch[0])
#
#     # complete months remaining from start year
#     for month in range(start[1] + 1, 13):
#         days += month_len(month, False)  # 1858 (epoch year) was not a leap year
#
#     # complete years between start and end
#     # loop through complete years between epoch year and current year, aggregating their days
#     for year in range(start[2] + 1, end[2]):
#         # is it a leap year?
#         if isLeap(year):
#             days += 366
#         else:
#             days += 365
#
#     # complete months from end year
#     if isLeap(end[2]):
#         leap = True
#     else:
#         leap = False
#
#     for month in range(1, end[1]):
#         days += month_len(month, leap)
#
#     # complete days in end month (not including end day)
#     days += end[0] - 1
#
#     # if date in question is past epoch, must add complete start day (since epoch day is a complete day) and the time that has elapsed in the end day
#     if start == epoch:
#         days += ((end[3][0] + (end[3][1] / 60)) / 24) + 1  # (end[3] = tt)
#
#     # if, however, given date precedes the epoch, only add remainder of start day
#     elif end == epoch:
#         days += ((24 - start[3][0]) + (60 - start[3][1] / 60)) / 24
#
#     return days


# given number of days since equinox, time of equinox in GMT, and phi difference between Greenwich and location in question,
# determine phi value of that location
def phi_val(days_since_eq, phi_dif):

    eq_hours = int(equinox_time[0])
    eq_mins = int(equinox_time[1])

    time_since_noon = ((eq_hours - 12) + (eq_mins / 60)) / 24

    # find phi for greenwich at equinox (by multiplying fraction of day by 360.98 degrees)
    # greenwich_phi_at_equinox = time_since_noon * degrees_day
    greenwich_phi_at_equinox = time_since_noon * radians_day

    # find phi for place in question at equinox (by adding difference in their phi values)
    local_phi_at_equinox = greenwich_phi_at_equinox + phi_dif

    # local_phi_now = (local_phi_at_equinox + (days_since_eq * degrees_day)) % 360
    local_phi_now = (local_phi_at_equinox + (days_since_eq * radians_day)) % (2*np.pi)

    return local_phi_now


# convert normalized vector from spherical to cartesian coordinates
def n_spherical_to_cartesian(theta, phi):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return np.array([x, y, z])


def cartesian_to_spherical(v):
    x, y, z = np.array([float(v[0]), float(v[1]), float(v[2])])
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)

    return np.array([r, phi, theta])

####################


# given day and angles of the location in question, determine position vectors from Earth's center to cur location and from cur location to the sun
def zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr):
    # days since spring equinox
    days_since_equinox = mjd_cur - mjd_equinox

    # location of Jerusalem in terms of theta and phi
    phi = phi_val(days_since_equinox, phi_dif)

    # Jerusalem's normalized position vector
    n_jer_equ = n_spherical_to_cartesian(theta, phi)        # equatorial coords
    n_jer_ecl = equatorial_to_ecliptic(n_jer_equ, mjd_cur)  # ecliptic coords

    p_sun = get_pos(mjd_cur, 's', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_sun = (p_sun-p_earth) * au_cm

    # sun's position relative to Jerusalem
    n_sun_jer = n_body_jer(n_jer_ecl, p_sun)                      # ecliptic coords
    n_sun_jer_equ = ecliptic_to_equatorial(n_sun_jer, mjd_cur)  # equatorial coords

    return n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ


# angle between zenith of location in question and the vector from that location to the sun
def psi_val(mjd_cur, phi_dif, theta, positionsArr):

    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr)

    # dot product of n_jer in ecliptic coords and n_sun_jer
    dot_product = n_jer_ecl.dot(n_sun_jer)

    # resulting angle between normalized vector from center of earth to Jer and normalized vector from Jer to sun
    psi = np.arccos(float(dot_product))

    return psi


# create observer coordinate system
def observer_axes(n_jer_equ):

    # z-axis (zenith)
    n_NP = np.array([0, 0, 1])  # normalized North Pole vector
    e_Z = n_jer_equ

    # East (x) axis
    e_E = np.cross(n_NP, e_Z)
    len_e_E = np.sqrt(e_E.dot(e_E))
    e_E_normal = e_E / len_e_E

    # North (y) axis
    e_N = np.cross(e_Z, e_E_normal)
    len_e_N = np.sqrt(e_N.dot(e_N))
    e_N_normal = e_N / len_e_N

    return np.array([e_N_normal, e_E_normal, e_Z])


# convert position vector from equatorial to observer coords
def n_observer_coords(n, cur_observer_axes):

    e_N, e_E, e_Z = cur_observer_axes
    E_component = n.dot(e_E)
    N_component = n.dot(e_N)
    Z_component = n.dot(e_Z)

    return np.array([E_component, N_component, Z_component])


# convert from observer coords to psi and phi (angles of the sun relative to the observer)
def observer_angles(E_coord, N_coord, Z_coord):
    psi = acos(Z_coord)            # latitude equivalent
    phi = atan2(E_coord, N_coord)  # longitude equivalent

    if phi < 0: phi = phi + 2*np.pi

    # phi = float(phi) * 180/numpy.pi
    # psi = 90 - (psi * 180/numpy.pi)

    return psi, phi


# convert psi and phi angles to position vector in equatorial coords
def n_from_observer_angles(psi, phi, cur_observer_axes):

    cartesian_obs_coords = n_spherical_to_cartesian(psi, phi)       # in observer coordinates
    cartesian_equ_coords = cartesian_obs_coords.dot(cur_observer_axes)  # now in equatorial coordinates

    return cartesian_equ_coords


def psi_vs_time(mjd_cur, time_zone, dst_start, dst_end, positionsArr, time_step=0.01):

    phi_dif, theta = phiDifAndTheta(time_zone)

    # plot psi vs. time
    x = []
    y = []

    for day in range(1):
        if day:  mjd_cur += 1

        time_difference = time_dif(mjd_cur, time_zone, dst_start, dst_end)
        time = 0
        # time = 0-time_difference

        while time < 24:
            # find psi value (angle between epoch and sun position vector) at this point in time
            mjd_frac = mjd_cur + time / 24
            psi = psi_val(mjd_frac, phi_dif, theta, positionsArr)

            acc_time = time + (24 * day) + time_difference
            x.append(acc_time)
            y.append(psi)

            time += time_step

    # plt.xticks(np.arange(0, 25, 1))
    # plt.plot(x, y, color='blue', marker='.', )
    # plt.show()

    # return (psi2 - psi1) / (t2 - t1)
    return x, y


def critical_point(polynomial, hr, interval=1.0):

    X = symbols('X')

    # take derivative of given polynomial; check if x value lies within given hour-range
    derivative = diff(polynomial, X)
    critical_points = solveset(derivative, X)

    for cp in critical_points:
        if hr <= cp <= hr + interval:  # if the time falls within the given hour
            return cp

    return ''


def sunriseSunset(mjd_cur, time_zone, dst_start, dst_end, positionsArr):

    # get list of psi-values for each hour during this day
    x, y = psi_vs_time(mjd_cur, time_zone, dst_start, dst_end, positionsArr, 1)

    # construct polynomial of psi vs. time for each set of 3 hours; check at which times the angle is 91 degrees and at what time it's at a max
    noon = 0
    ans = []
    for i in range(1, 21):

        # get polynomial defined by the angles this hour, the previous hour, and the two following hours
        hours = x[i - 1:i + 2]  # slice of larger hour list
        angles = y[i - 1:i + 2]  # slice of larger angle list

        #### System of Equations
        system = []
        a, b, c = symbols('a b c')

        # create system of equations using given 3 (x,y) pairs
        for j in range(3):
            equation = ((hours[j]**2)*a) + (hours[j]*b) + c - angles[j]
            system.append(equation)

        # list of coefficients of solution to linear system
        coeff = linsolve(system, [a, b, c])

        X, Y = symbols('X Y')

        # get the one tuple in coefficient FiniteSet (the polynomial)
        polynomial = 0
        for tup in coeff:
            polynomial = (tup[0] * X**2) + (tup[1] * X) + tup[2] - sunset_angle
        #####################

        # get time of solar noon
        potential_noon = critical_point(polynomial+sunset_angle, i)
        if potential_noon:  noon = float(potential_noon)

        # for each of the roots, check whether it's real, and it falls within the given hour
        roots = solveset(polynomial, X)
        for r in roots:
            if r.is_real and i <= r <= i + 1:   # if the time falls within the given hour
                ans.append(float(r))

    # return times of sunrise, solar noon, and sunset - in both cleaned-up and raw time
    return ('sunrise: ' + clean_time(ans[0])  + '; solar noon: ' + clean_time(noon) + '; sunset: ' + clean_time(ans[1])), (ans[0], noon, ans[1])


def clean_time(hours):

    m = (hours % 1) * 60
    s = str(int((m % 1) * 60))
    m = str(int(m))

    if len(m) == 1: m = '0' + m
    if len(s) == 1: s = '0' + s

    return str(int(hours)) + ':' + m + ':' + s


def earliest_sunset(yr, mjd_cur, time_zone, dst_start, dst_end, positionsArr):

    # get mjd of first day of given year
    mjd_cur = mjd(1, 1, yr, (0, 0))

    earliest = 24  # hours
    earliestDay = None

    if isLeap(yr): yr_len = 366
    else:          yr_len = 365

    for day in range(50):
    # for day in range(yr_len):

        # find time of sunset
        sunset = sunriseSunset(mjd_cur, time_zone, dst_start, dst_end, positionsArr)[1][2]

        if sunset <= earliest:
            earliest = sunset
            earliestDay = day

        # move to next day
        mjd_cur += 1

    return clean_time(earliest), earliestDay


def timesAndAngle(d, m, y, mjd_cur, time_zone, DST_START, DST_END, positionsArr):

    phi_dif, theta = phiDifAndTheta(time_zone)

    times = sunriseSunset(mjd_cur, time_zone, DST_START, DST_END, positionsArr)

    # time of sunrise and sunset (via interpolations)
    print(str(d) + '/' + str(m) + '/' + str(y) + ': ' + times[0])

    MJD_sunrise = mjd_cur + (times[1][0]/24) - 2/24
    MJD_noon = mjd_cur + (times[1][1]/24) - 2/24
    MJD_sunset = mjd_cur + (times[1][2]/24) - 2/24

    ## phi at SUNRISE in observer coordinates:
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_sunrise, phi_dif, theta, positionsArr)

    # set up observer axes
    axes = observer_axes(n_jer_equ)

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = observer_angles(coords[0], coords[1], coords[2])

    ### phi at SUNSET
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_sunset, phi_dif, theta, positionsArr)
    axes = observer_axes(n_jer_equ)  # observer axes

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = observer_angles(coords[0], coords[1], coords[2])

    ### phi at NOON
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_noon, phi_dif, theta, positionsArr)
    axes = observer_axes(n_jer_equ)  # observer axes

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = observer_angles(coords[0], coords[1], coords[2])
    print('psi, phi at noon:', psi, phi)

    # for hour in range(24):
    #     mjd_cur += 1 / 24
    #     n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(mjd_cur)
    #     axes = observer_axes(n_jer_equ)
    #
    #     # sun-Jerusalem vector in observer coordinates
    #     coords = n_observer_coords(n_sun_jer_equ, axes)
    #     psi, phi = n_observer_angles(coords[0], coords[1], coords[2])
    #     print(hour, psi, phi)
    #     print(psi_vs_time(mjd_cur, time_zone, DST_START, DST_END, 1))


def psi_vs_phi(year, hour, time_zone, time_difference, positionsArr, start_day=1, end_day=365):

    phi_dif, theta = phiDifAndTheta(time_zone)

    # start at given time on the first day of the given year
    mjd_cur = mjd(1, 1, year, (hour - time_difference, 0))

    # lists to hold phi/psi coordinate values
    x = []
    y = []

    equinoxX = []
    equinoxY = []

    # loop through each day in the given year, getting phi and psi values at 12pm on each day
    for day in range(start_day, end_day+1):

        mjd_cur += 1

        # zenith and sun vectors at 12pm on this day (in both ecliptic and equatorial coords)
        n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr)

        # set up observer axes
        axes = observer_axes(n_jer_equ)

        # sun-Jerusalem vector in observer coordinates
        coords = n_observer_coords(n_sun_jer_equ, axes)
        psi, phi = observer_angles(coords[0], coords[1], coords[2])
        psi = psi * 180/np.pi
        phi = phi * 180/np.pi
        print(day, psi, phi)

        if day == 79 or day == 172 or day == 265 or day == 355:
            equinoxX.append(phi)
            equinoxY.append(psi)

        x.append(phi)
        y.append(psi)

    # plot psi vs. phi
    plt.scatter(x, y)

    # plot points for solstices and equinoxes
    plt.scatter(equinoxX, equinoxY, color='red')

    plt.title("Psi vs. Phi at " + clean_time(hour) + time_zone)
    plt.show()


# set up Right and Up axes given the angles of the camera (relative to the observer)
def camera_axes(mjd_cur, p_jer_ecl, positionsArr, obs_axes, center_moon=True):

    e_N_equ, e_E_equ, e_Z_equ = obs_axes
    e_Z_ecl = equatorial_to_ecliptic(e_Z_equ, mjd_cur)
    e_N_ecl = np.array([0, 0, 1])

    p_moon = get_pos(mjd_cur, 'm', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_jer = p_jer_ecl + p_earth

    # camera vector (ecliptic coords)
    # 1. centered at the moon
    if center_moon:
        n_cam = normalize(p_moon-p_jer)
        n_cam = np.array([float(n_cam[0]), float(n_cam[1]), float(n_cam[2])])

    # 2. according to user-defined camera angles
    else:
        n_cam_equ = n_from_observer_angles(psi_cam, phi_cam, obs_axes)
        n_cam = equatorial_to_ecliptic(n_cam_equ, mjd_cur)

    # Right and Up axis vectors (normalized)
    e_R = np.cross(n_cam, e_Z_ecl)      # relative to the zenith
    # e_R = np.cross(n_cam, e_N_ecl)    # relative to the north pole
    e_R_normal = normalize(e_R)

    e_U = np.cross(e_R_normal, n_cam)
    e_U_normal = normalize(e_U)

    return np.array([e_R_normal, e_U_normal, n_cam])


# convert position vector from equatorial coords to camera coords
def n_camera_coords(n, camera_axis_vectors):

    # set up Right and Up axes
    e_R, e_U = camera_axis_vectors[0], camera_axis_vectors[1]

    # determine the given position vector's Right and Up coordinate values
    x = n.dot(e_R)
    y = n.dot(e_U)

    return np.array([x, y])


# plot the position of the sun in camera coordinates at a given time on each day of the given year
def y_vs_x(year, hour, time_zone, time_difference, psi_cam, phi_cam, body, positionsArr, start_day=85, end_day=363):

    phi_dif, theta = phiDifAndTheta(time_zone)

    # lists to hold x/y coordinate values
    xList = []
    yList = []
    equinoxX = []   # coordinates at equinoxes and solstices
    equinoxY = []
    monthsX = []    # coordinates at the first of each month
    monthsY = []
    roshChodeshX = []
    roshChodeshY = []

    sunX = []
    sunY = []

    # month and day-of-month pointers beginning at 0
    month_cur = 1
    day_of_month_cur = 0

    leap = isLeap(year)

    # start at the given hour on the first day of the given year
    mjd_cur = mjd(start_day - 1, 1, year, (hour - time_difference, 0))  # subtract one because add 1 in loop

    # for each day in the given year, get sun position at given time and convert to camera coordinates
    # for day in range(start_day, end_day):
    for day in range(start_day, start_day+1):
        day_of_month_cur += 1
        mjd_cur += 1

        # zenith and sun vectors at given time on this day (in both ecliptic and equatorial coords)
        n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr)

        p_moon = (get_pos(mjd_cur, 'm', positionsArr) - get_pos(mjd_cur, 'e', positionsArr)) * au_cm
        p_moon_equ = ecliptic_to_equatorial(p_moon, mjd_cur)

        # moon's position relative to Jerusalem
        n_moon_jer_equ = n_body_jer(n_jer_equ, p_moon_equ)

        # set up OBSERVER axes according to current position of Jerusalem
        axes = observer_axes(n_jer_equ)

        # now, set up CAMERA axes using observer axes and given position of camera
        cam_axes = camera_axes(mjd_cur, n_jer_ecl*earth_rad, positionsArr, axes)

        # if body == 's':
        # and determine sun position vector in camera coords
        sun_x, sun_y = n_camera_coords(n_sun_jer_equ, cam_axes)

        # elif body == 'm':
        x, y = n_camera_coords(n_moon_jer_equ, cam_axes)

        draw_moon(mjd_cur, cam_axes, n_jer_equ*earth_rad, positionsArr)
        # RA_Dec_lines(cam_axes, n_jer_equ*earth_rad)

        ### store the determined x/y values in the correct list for later plotting
        # days around and including rosh chodesh
        if len(rosh_chodesh_MJD) > 0 and (int(mjd_cur) == rosh_chodesh_MJD[0][0] or int(mjd_cur)-1 == rosh_chodesh_MJD[0][0]
                                          or int(mjd_cur)+1 == rosh_chodesh_MJD[0][0]):

            roshChodeshX.append(x)
            roshChodeshY.append(y)

            sunX.append(sun_x)
            sunY.append(sun_y)
            # plt.text(x, y, rosh_chodesh_MJD[0][1])

            if int(mjd_cur) > rosh_chodesh_MJD[0][0]:
                rosh_chodesh_MJD = rosh_chodesh_MJD[1:]   # slice off the first element

        # equinoxes and solstices (year 2023)
        if day == 79 or day == 172 or day == 265 or day == 355:
            equinoxX.append(x)
            equinoxY.append(y)

            # if day == 79:       plt.text(x, y, " Spring Equinox", color="red", ha='right', fontfamily="times new roman", size=11)
            # elif day == 172:    plt.text(x, y, " Summer Solstice", color="red", va='bottom', fontfamily="times new roman", size=11)
            # elif day == 265:    plt.text(x, y, " Fall Equinox", color="red", fontfamily="times new roman", size=11)
            # else:               plt.text(x, y, " Winter Solstice", color="red", horizontalalignment='center', va='bottom', fontfamily="times new roman", size=11)

        # move to the next month if necessary
        if day_of_month_cur > month_len(month_cur, leap) and month_cur < 12:
            day_of_month_cur = 1
            month_cur += 1

        # if it's the first of the month:
        if day_of_month_cur == 1:
            monthsX.append(x)
            monthsY.append(y)

            # create month labels
            # if jan, feb, mar, apr, jul, or aug, label should lie to the left of the data point
            leftMonths = [1, 2, 3, 4, 7, 8]
            # if month_cur in leftMonths:
            #     plt.text(x, y, DateConversion.MONTHS[month_cur-1], fontfamily="times new roman", va='top', ha='right', color='black')
            #
            # else:   # otherwise, label should lie to the right
            #     plt.text(x, y, DateConversion.MONTHS[month_cur - 1], fontfamily="times new roman", va='top',
            #              color='black')

        # in any case, add this x/y point to the general x and y lists
        xList.append(x)
        yList.append(y)

    # fig = plt.figure(figsize=(5, 5))

    # # plot all y vs. x points
    # plt.scatter(xList, yList, color='blue', linewidths=0)

    # # plot points for solstices and equinoxes
    # plt.scatter(equinoxX, equinoxY, color='red', marker='*')
    #
    # # plot points for the first of each month
    # plt.scatter(monthsX, monthsY, color='yellow', marker='.')

    # plot rosh chodesh points
    plt.scatter(roshChodeshX, roshChodeshY, color='green')

    # plot sun points
    plt.scatter(sunX, sunY, color='red')

    # titles
    plt.title("Up vs. Right at " + clean_time(hour) + ' ' + time_zone, fontfamily="times new roman")
    plt.xlabel("Right", fontfamily="times new roman")
    plt.ylabel("Up", fontfamily="times new roman")

    plt.gca().set_aspect('equal', adjustable='datalim')
    # plt.show()


def plot_sun_moon(mjd_cur, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr):
    n_cam = cam_axes[2]
    p_moon = get_pos(mjd_cur, 'm', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_sun = get_pos(mjd_cur, 's', positionsArr)

    n_sun_earth = normalize(p_sun - p_earth)

    # sun_earth_dist = 1.5e14
    # n_sun_earth1 = normalize(n_jer_ecl*earth_rad + n_sun_jer_ecl*sun_earth_dist)
    # print(n_sun_earth1, n_sun_earth)

    p_moon_earth_ecl = p_moon-p_earth

    # moon's position relative to Jerusalem
    n_moon_jer_ecl = n_body_jer(n_jer_ecl, p_moon_earth_ecl)

    # sun and moon vectors in camera coords
    sun_x, sun_y = n_camera_coords(n_sun_jer_ecl, cam_axes)
    moon_x, moon_y = n_camera_coords(n_moon_jer_ecl, cam_axes)

    # plot
    if n_sun_jer_ecl.dot(n_cam) > 0:  # and np.dot(n_sun_earth, n_jer_ecl) > 0:
        plt.scatter(sun_x, sun_y, color='orange', s=400, marker='*')

    if n_moon_jer_ecl.dot(n_cam) > 0:  # and np.dot(normalize(p_moon_earth_ecl), n_jer_ecl) > 0:
        plt.scatter(moon_x, moon_y, color='blue', s=50, marker='.')

    return [sun_x, sun_y], [moon_x, moon_y]


def drawHorizon(cam_axes, n_jer_equ, mjd_cur):
    n_cam = cam_axes[2]

    e_Z = equatorial_to_ecliptic(n_jer_equ, mjd_cur)
    e_R = np.cross(n_cam, e_Z)
    e_R_normal = normalize(e_R)

    e_F = np.cross(e_Z, e_R_normal)
    e_F_normal = normalize(e_F)

    points = []

    for i in range(0, 360):
        psi = i * np.pi/180
        p = (e_R_normal * np.cos(psi)) + (e_F_normal * np.sin(psi))
        x, y = n_camera_coords(p, cam_axes)  # in camera coords

        if p.dot(n_cam) > 0:  # if visible
            points.append([x, y])

    for i in range(len(points)-1):
        plt.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 'b-', markersize=1)


def RA_Dec_lines(cam_axes, p_jer_equ, mjd_cur):

    # xMin = yMin = -0.999
    # xMax = yMax = 0.999
    # xRange = xMax - xMin
    # yRange = yMax - yMin

    n_cam = cam_axes[2]

    # set plot axis boundaries
    # plt.xlim(-0.25, 0.25)
    # plt.ylim(-0.25, 0.25)
    # plt.xlim(-1, 1)
    # plt.ylim(-1, 1)
    plt.xlim(-0.75, 0.75)
    plt.ylim(-0.75, 0.75)
    plt.gca().set_aspect('equal', adjustable='box')

    # plot horizon
    drawHorizon(cam_axes, normalize(p_jer_equ), mjd_cur)

    # RA LINES
    # theta is held constant; phi varies
    for i in range(0, 180, 30):
        theta = i * np.pi/180
        points = []
        vectors = []
        for j in range(0, 360):
            phi = j * np.pi / 180

            # Project this point onto x/y plane (i.e. convert to camera coordinates)
            n_equ = n_spherical_to_cartesian(theta, phi)  # unit position vector (centered at earth's center)
            p_equ = celestial_rad * n_equ                 # scaled up by radius of celestial sphere
            n_ecl = equatorial_to_ecliptic(n_equ, mjd_cur)# unit position vector in ecliptic coordinates
            x, y = n_camera_coords(n_ecl, cam_axes)       # finally, in camera coords

            if p_equ.dot(p_jer_equ) > 0 and n_equ.dot(n_cam) > 0:  # if visible
                points.append([x, y])
                vectors.append(p_equ)

        # draw line connecting each of the visible points as long as the two points are not on opposite sides of the sphere
        for j in range(len(points) - 1):
            ang_dist = angular_dist(vectors[j], vectors[j+1])   # in radians
            ang_dist = ang_dist * 180/np.pi                     # in degrees
            if ang_dist < 2:
                plt.plot([points[j][0], points[j + 1][0]], [points[j][1], points[j + 1][1]], 'k-')

    # DECLINATION LINES
    # phi is held constant; theta varies
    for i in range(0, 361, 30):
        phi = i * np.pi / 180
        points = []
        vectors = []
        for j in range(0, 181):
            theta = j * np.pi / 180

            # Project this point onto x/y plane (i.e. convert to camera coordinates)
            n_equ = n_spherical_to_cartesian(theta, phi)  # unit position vector (centered at earth's center)
            p_equ = celestial_rad * n_equ                 # scaled up by radius of celestial sphere
            n_ecl = equatorial_to_ecliptic(n_equ, mjd_cur)# unit position vector in ecliptic coordinates
            x, y = n_camera_coords(n_ecl, cam_axes)       # finally, in camera coords

            if p_equ.dot(p_jer_equ) > 0 and n_equ.dot(n_cam) > 0:  # if visible
                points.append([x, y])
                vectors.append(p_equ)

            # draw line connecting each of the visible points as long as the two points are not on opposite sides of the sphere
        for j in range(len(points) - 1):
            ang_dist = angular_dist(vectors[j], vectors[j + 1])  # in radians
            ang_dist = ang_dist * 180 / np.pi  # in degrees
            if ang_dist < 2:
                plt.plot([points[j][0], points[j + 1][0]], [points[j][1], points[j + 1][1]], 'k-')

    # plt.show()


def sky_snapshots(mjd_cur, phi_dif, theta, positionsArr, total_days, steps_per_day):

    N = 100
    # D_fov = 0.05  # field of view, in radians
    D_fov = 0.1  # field of view, in radians
    sun_trail = []
    moon_trail = []

    num_iterations = int(total_days * steps_per_day)
    for i in range(num_iterations):
        print(i)
        n_jer_ecl, n_jer_equ, n_sun_jer_ecl, n_sun_jer_equ = zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr)

        # set up OBSERVER axes according to current position of Jerusalem
        obs_axes = observer_axes(n_jer_equ)

        # now, set up CAMERA axes using observer axes and given position of camera
        cam_axes = camera_axes(mjd_cur, n_jer_ecl*earth_rad, positionsArr, obs_axes, center_moon=False)

        ############
        # now we can plot RA and Dec lines and sun/moon
        f = plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')

        sun_coords, moon_coords = plot_sun_moon(mjd_cur, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr)
        RA_Dec_lines(cam_axes, n_jer_equ * earth_rad, mjd_cur)

        plt.title("MJD:" + str(mjd_cur))
        plt.xlabel("Right", fontfamily="times new roman")
        plt.ylabel("Up", fontfamily="times new roman")

        # plt.show()
        f.savefig("./SkyMaps/"+str(i))
        plt.close()

        ############
        ## Now DRAW sun and moon
        # imageData = draw_sun(mjd_cur, cam_axes, n_jer_ecl*earth_rad, imageData, positionsArr)
        # draw_moon(mjd_cur, cam_axes, n_jer_ecl*earth_rad, moon_coords, imageData, positionsArr)
        imageData, sun_trail = draw_sun2(mjd_cur, cam_axes, n_jer_ecl * earth_rad, positionsArr, N, D_fov, sun_trail, i)
        imageData, moon_trail = draw_moon2(mjd_cur, cam_axes, n_jer_ecl * earth_rad, positionsArr, imageData, N, D_fov, moon_trail, i)

        image = Image.fromarray(imageData)
        file_name = "./SkyImages/" + str(i) + ".png"
        image.save(file_name)

        # now add trails
        for j in range(len(sun_trail)):
            imageData[sun_trail[j]] = [255, 255, 255]
        for j in range(len(moon_trail)):
            imageData[moon_trail[j]] = [255, 255, 255]
        #############

        mjd_cur += (1/steps_per_day)


def plot_ecliptic(cam_axes):
    ecliptic_points = []
    ecliptic_rad = 100 * au_cm

    # plot ecliptic
    for phi in range(0, 360):
        phi_rad = phi * np.pi / 180
        n_star = normalize(np.array([ecliptic_rad * np.cos(phi_rad), ecliptic_rad * np.sin(phi_rad), 0]))
        x, y = n_camera_coords(n_star, cam_axes)

        ecliptic_points.append([x, y])


    for i in range(len(ecliptic_points) - 1):
        plt.plot([ecliptic_points[i][0], ecliptic_points[i + 1][0]], [ecliptic_points[i][1], ecliptic_points[i + 1][1]], 'r-')



def sky_snapshots2(mjd_cur, phi_dif, theta, positionsArr, total_days, steps_per_day):

    sun_trail = []
    moon_trail = []

    num_iterations = int(total_days * steps_per_day)
    for i in range(num_iterations):
        print(i)
        n_jer_ecl, n_jer_equ, n_sun_jer_ecl, n_sun_jer_equ = zenithAndSunVectors(mjd_cur, phi_dif, theta, positionsArr)

        obs_axes = observer_axes(n_jer_equ)
        cam_axes = camera_axes(mjd_cur, n_jer_ecl*earth_rad, positionsArr, obs_axes, center_moon=False)

        # plot RA and Dec lines and sun/moon
        f = plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')

        sun_coords, moon_coords = plot_sun_moon(mjd_cur, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr)
        RA_Dec_lines(cam_axes, n_jer_equ * earth_rad, mjd_cur)

        # plot sun and moon trails
        sun_trail.append(sun_coords)
        moon_trail.append(moon_coords)

        for j in range(len(sun_trail)-1):
            plt.plot([sun_trail[j][0], sun_trail[j + 1][0]], [sun_trail[j][1], sun_trail[j + 1][1]], 'r-')

        for j in range(len(moon_trail)-1):
            plt.plot([moon_trail[j][0], moon_trail[j + 1][0]], [moon_trail[j][1], moon_trail[j + 1][1]], 'k-')

        plot_ecliptic(cam_axes)

        plt.title("MJD:" + str(mjd_cur))
        plt.xlabel("Right", fontfamily="times new roman")
        plt.ylabel("Up", fontfamily="times new roman")

        plt.show()
        # f.savefig("./SkyMaps/"+str(i))
        plt.close()

        mjd_cur += (1/steps_per_day)


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [
        int(text)
        if text.isdigit() else text.lower()
        for text in _nsre.split(s)]


def movie(folder_name, movie_name):

    path = 'C:/Users/tywer/Documents/Celestial Coordinates/' + folder_name + '/*.png'
    img_array = [img for img in glob.glob(path)]
    sorted_images = sorted(img_array, key=natural_sort_key)

    img = cv2.imread(img_array[0])
    height, width, layers = img.shape
    size = (width, height)

    out = cv2.VideoWriter(movie_name+'.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 8, size)

    for image in sorted_images:
        out.write(cv2.imread(image))

    cv2.destroyAllWindows()
    out.release()


# original draw_moon function, looping over points on the moon and projecting each onto the canvas
def draw_moon(mjd_cur, cam_axes, p_jer_ecl, moon_coords, imageData, positionsArr):

    moonColorIm = Image.open('lroc_color_poles_2k.tif')
    im_array = np.asarray(moonColorIm)

    fov = 0.04   # field of view
    # moon_x, moon_y = moon_coords
    moon_x = moon_y = 0
    xMin = moon_x-(fov/2)
    xMax = moon_x+(fov/2)
    xRange = xMax - xMin
    yMin = moon_y-(fov/2)
    yMax = moon_y+(fov/2)
    yRange = yMax - yMin

    # n_cam = np.cross(cam_axes[1], cam_axes[0])

    p_sun = get_pos(mjd_cur, 's', positionsArr)
    p_moon = get_pos(mjd_cur, 'm', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_jer = p_jer_ecl + p_earth

    yellowX = []
    yellowY = []
    blueX = []
    blueY = []

    # loop through each point on the moon
    for i_theta in range(0, 180):
        theta = i_theta * np.pi/180

        for i_phi in range(0, 360):
            phi = i_phi * np.pi / 180

            # this point in ecliptic cartesian coords
            s = n_spherical_to_cartesian(theta, phi) * moon_rad  # vector from center of moon to current point (in moon coords)
            p_s = p_moon + s  # ecliptic coords

            # vector from jer to this point on the moon
            n_moon_jer_ecl = normalize(p_s-p_jer)
            x, y = n_camera_coords(n_moon_jer_ecl, cam_axes)  # camera (x-y) coords

            canvasX = round(((x - xMin) / xRange) * canvas_size)
            canvasY = round(((yMax - y) / yRange) * canvas_size)

            # map this (phi, theta) point onto the 2D moon-color image
            x_val = (2048 * ((i_phi + 180)/360)) % 2048
            y_val = 1024 * (i_theta/180)

            color = im_array[round(y_val)][round(x_val)]

            if s.dot(p_s - p_jer) < 0:
                if s.dot(p_s-p_sun) >= 0:
                    imageData[canvasY - 1:canvasY + 1, canvasX - 1:canvasX + 1] = color/4
                else:
                    imageData[canvasY - 1:canvasY + 1, canvasX - 1:canvasX + 1] = color

            # if visible to earth but not sun, dark
            # if s.dot(p_s-p_jer) < 0 and s.dot(p_s-p_sun) >= 0:
            #     canvasX = round(((x - xMin) / xRange) * canvas_size)
            #     canvasY = round(((yMax - y) / yRange) * canvas_size)
            #     if canvasX == 256: canvasX = 255
            #     if canvasY == 256: canvasY = 255
            #     imageData[canvasY, canvasX] = [47, 79, 79]
            #     imageData[canvasY-1:canvasY+1, canvasX-1:canvasX+1] = [47, 79, 79]
            #
            # # if visible to earth AND sun, illuminated
            # elif s.dot(p_s-p_jer) < 0 and s.dot(p_s-p_sun) < 0:
            #     canvasX = round(((x - xMin) / xRange) * canvas_size)
            #     canvasY = round(((yMax - y) / yRange) * canvas_size)
            #     if canvasX == 256: canvasX = 255
            #     if canvasY == 256: canvasY = 255
                # imageData[canvasY-1:canvasY+1, canvasX-1:canvasX+1] = [255, 255, 255]

    image = Image.fromarray(imageData)
    image.show()

    plt.scatter(yellowX, yellowY, color='yellow')
    plt.scatter(blueX, blueY, color='blue')
    return imageData
    # return image


# given an axis vector of the moon and its rotation phase, rotate that vector
def rot_moon_axis(e, rotation_phase, eps=moon_obliquity):

    # Rotation matrix
    row0 = np.array([np.cos(rotation_phase), -np.cos(eps)*np.sin(rotation_phase), -np.sin(eps)*np.sin(rotation_phase)])
    row1 = np.array([np.cos(eps)*np.sin(rotation_phase),
                     np.cos(rotation_phase) + (np.sin(eps)**2*(1 - np.cos(rotation_phase))),
                     -np.sin(eps)*np.cos(eps)*(1 - np.cos(rotation_phase))])
    row2 = np.array([np.sin(eps)*np.sin(rotation_phase),
                     -np.sin(eps)*np.cos(eps)*(1 - np.cos(rotation_phase)),
                     np.cos(rotation_phase) + (np.cos(eps))**2*(1 - np.cos(rotation_phase))])

    # result of matrix-vector multiplication
    x = np.dot(e, row0)
    y = np.dot(e, row1)
    z = np.dot(e, row2)

    return np.array([x, y, z])


def rotated_moon_axes(mjd_cur, eps=moon_obliquity):

    em_x = np.array([1, 0, 0])
    em_y = np.array([0, np.cos(eps), np.sin(eps)])
    em_z = np.array([0, -np.sin(eps), np.cos(eps)])

    days_since_epoch = mjd_cur - mjd_j2000

    rotation_phase = rotation_phase0 + ((2 * np.pi * days_since_epoch) / T_side)

    # given rotation phase, rotate surface of moon around em_z
    es_x = rot_moon_axis(em_x, rotation_phase)
    es_y = rot_moon_axis(em_y, rotation_phase)
    es_z = rot_moon_axis(em_z, rotation_phase)

    return np.array([es_x, es_y, es_z])


# new draw_moon function, looping through each pixel in the sky and checking if it contains the moon
def draw_moon2(mjd_cur, cam_axes, p_jer_ecl, positionsArr, imageData, N, D_fov, trail, image_num):

    moonColorIm = Image.open('lroc_color_poles_2k.tif')
    im_array = np.asarray(moonColorIm)
    trail_marked = False

    # moon axes rotated around the current phase angle
    es_x, es_y, es_z = rotated_moon_axes(mjd_cur)

    e_R, e_U, n_cam = cam_axes
    p_moon = get_pos(mjd_cur, 'm', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_sun = get_pos(mjd_cur, 's', positionsArr)

    p_jer_solsys = p_earth + p_jer_ecl  # position of jer relative to the center of the solar system
    x_c, y_c, z_c = float(p_moon[0]), float(p_moon[1]), float(p_moon[2])
    x_p, y_p, z_p = float(p_jer_solsys[0]), float(p_jer_solsys[1]), float(p_jer_solsys[2])

    for i in range(-N, N):
        for j in range(-N, N):

            # vector pointing from the observer to the i,j'th pixel
            n_pix = normalize(n_cam + ((i/N) * D_fov/2 * e_R) + ((j/N) * D_fov/2 * e_U))
            n_x, n_y, n_z = n_pix

            a = n_x**2 + n_y**2 + n_z**2
            b = 2*(n_x*(x_p-x_c) + n_y*(y_p-y_c) + n_z*(z_p-z_c))
            c = x_p**2 + x_c**2 - (2*x_p*x_c) + y_p**2 + y_c**2 - (2*y_p*y_c) + z_p**2 + z_c**2 - (2*z_p*z_c) - moon_rad**2

            discrim = b ** 2 - (4 * a * c)

            # if there are real roots (i.e. discriminant is non-negative), then this n_pix must intersect the moon
            if discrim >= 0:
                roots = [(-b + np.sqrt(discrim)) / (2 * a),
                         ((-b - np.sqrt(discrim)) / (2 * a))]

                # points on surface of moon with which n_pix intersects
                p_s1 = p_jer_solsys + (roots[0] * n_pix)
                p_s2 = p_jer_solsys + (roots[1] * n_pix)

                # vectors from center of the moon to those points on its surface
                s1 = p_s1 - p_moon
                s2 = p_s2 - p_moon

                # rotate s1 and s2 according to the moon's rotation
                s1_rotated = np.array([np.dot(s1, es_x), np.dot(s1, es_y), np.dot(s1, es_z)])
                s2_rotated = np.array([np.dot(s2, es_x), np.dot(s2, es_y), np.dot(s2, es_z)])

                # convert them to phi and theta
                r1, phi1, theta1 = cartesian_to_spherical(s1_rotated)
                r2, phi2, theta2 = cartesian_to_spherical(s2_rotated)

                # plot the point that's visible to the observer

                ## the following block never runs
                # if s1.dot(p_s1 - p_jer_solsys) < 0:
                #     # map this (phi, theta) point onto the 2D moon-color image
                #     x1 = (2048 * ((phi1 + 2*np.pi) / (2*np.pi))) % 2048
                #     y1 = 1024 * (theta1 / np.pi)
                #     color = im_array[int(y1)][int(x1)]
                #
                #     # illuminated?
                #     if s1.dot(p_s1 - p_sun) >= 0:
                #         imageData[N-j, i+N] = color / 4
                #     else:
                #         imageData[N-j, i+N] = color
                #
                #     # check if this is the center point - and if so, mark it
                #     if not trail_marked:
                #         c = x_p ** 2 + x_c ** 2 - (2 * x_p * x_c) + y_p ** 2 + y_c ** 2 - \
                #             (2 * y_p * y_c) + z_p ** 2 + z_c ** 2 - (2 * z_p * z_c) - (moon_rad/8)**2
                #         discrim = b ** 2 - (4 * a * c)
                #         if discrim >= 0:
                #             trail.append((N - j, i + N))
                #             trail.append((N - j + 1, i + N - 1))
                #             trail.append((N - j - 1, i + N + 1))
                #             # trail.append((N - j + 1, i + N + 1))
                #             # trail.append((N - j - 1, i + N - 1))
                #             trail_marked = True

                if s2.dot(p_s2 - p_jer_solsys) < 0:# and np.dot(n_pix, normalize(p_jer_ecl)) > 0:
                    x2 = (2048 * ((phi2 + 2*np.pi) / (2 * np.pi))) % 2048
                    y2 = 1024 * (theta2 / np.pi)

                    color = im_array[int(y2)][int(x2)]

                    # illuminated?
                    if s2.dot(p_s2 - p_sun) >= 0:
                        imageData[N-j, i+N] = color / 4
                    else:
                        imageData[N-j, i+N] = color

                    # check if this is the center point - and if so, mark it
                    if not trail_marked:
                        c = x_p ** 2 + x_c ** 2 - (2 * x_p * x_c) + y_p ** 2 + y_c ** 2 - \
                            (2 * y_p * y_c) + z_p ** 2 + z_c ** 2 - (2 * z_p * z_c) - (moon_rad/10)**2
                        discrim = b ** 2 - (4 * a * c)
                        if discrim >= 0:
                            trail.append((N - j, i + N))
                            # trail.append((N - j + 1, i + N - 1))
                            # trail.append((N - j - 1, i + N + 1))
                            trail.append((N - j + 1, i + N + 1))
                            trail.append((N - j - 1, i + N - 1))
                            trail_marked = True

            #
            # else:
            #     # if this point on the sky is visible
            #     if np.dot(n_pix, normalize(p_jer_ecl)) > 0:
            #         imageData[N - j, i + N] = [0, 0, 100]

    # mark trail from previous snapshots
    if trail_marked: items_to_exclude = 3
    else:            items_to_exclude = 0
    for i in range(len(trail) - items_to_exclude):
        imageData[trail[i]] = [255, 255, 255]

    return imageData, trail


# original draw_sun function, looping over points on the sun and projecting each onto the canvas
def draw_sun(mjd_cur, cam_axes, p_jer_ecl, imageData, positionsArr):
    n_cam = np.cross(cam_axes[1], cam_axes[0])

    p_sun = get_pos(mjd_cur, 's', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_jer = p_jer_ecl + p_earth

    yellowX = []
    yellowY = []
    blueX = []
    blueY = []

    # loop through each point on the sun
    for theta in range(0, 180, 5):
        theta = theta * np.pi / 180
        for phi in range(0, 360, 5):
            phi = phi * np.pi / 180

            # this point in ecliptic cartesian coords
            s = n_spherical_to_cartesian(theta,
                                         phi) * sun_rad  # vector from center of sun to current point (in moon coords)
            p_s = p_sun + s  # ecliptic coords
            n_sun_jer_equ = ecliptic_to_equatorial(normalize(p_s - p_jer), mjd_cur)
            x, y = n_camera_coords(n_sun_jer_equ, cam_axes)  # camera (x-y) coords

            # if visible to observer plot as yellow dot...
            if s.dot(p_s - p_jer) < 0:
                yellowX.append(x)
                yellowY.append(y)

            else:
                blueX.append(x)
                blueY.append(y)

    plt.scatter(yellowX, yellowY, color='yellow')
    plt.scatter(blueX, blueY, color='blue')
    return imageData
    # image = Image.fromarray(imageData)
    # return image


# new draw_sun function, looping through each pixel in the sky and checking if it contains the sun
def draw_sun2(mjd_cur, cam_axes, p_jer_ecl, positionsArr, N, D_fov, trail, image_num):

    canvas_size = (2 * N) + 1
    imageData = np.zeros((canvas_size, canvas_size, 3), dtype=np.uint8)
    trail_marked = False

    e_R, e_U, n_cam = cam_axes
    p_earth = get_pos(mjd_cur, 'e', positionsArr)
    p_sun = get_pos(mjd_cur, 's', positionsArr)

    p_jer_solsys = p_earth + p_jer_ecl  # position of jer relative to the center of the solar system
    x_c, y_c, z_c = float(p_sun[0]), float(p_sun[1]), float(p_sun[2])
    x_p, y_p, z_p = float(p_jer_solsys[0]), float(p_jer_solsys[1]), float(p_jer_solsys[2])

    for i in range(-N, N):
        for j in range(-N, N):

            # vector pointing from the observer to the i,j'th pixel
            n_pix = normalize(n_cam + ((i/N) * D_fov/2 * e_R) + ((j/N) * D_fov/2 * e_U))
            n_x, n_y, n_z = n_pix

            a = n_x**2 + n_y**2 + n_z**2
            b = 2*(n_x*(x_p-x_c) + n_y*(y_p-y_c) + n_z*(z_p-z_c))
            c = x_p**2 + x_c**2 - (2*x_p*x_c) + y_p**2 + y_c**2 - (2*y_p*y_c) + z_p**2 + z_c**2 - (2*z_p*z_c) - sun_rad**2

            discrim = b ** 2 - (4 * a * c)

            # if there are real roots (i.e. discriminant is non-negative), then this n_pix must intersect the sun
            if discrim >= 0:
                roots = [(-b + np.sqrt(discrim)) / (2 * a),
                         ((-b - np.sqrt(discrim)) / (2 * a))]

                # points on surface of sun with which n_pix intersects
                p_s1 = p_jer_solsys + (roots[0] * n_pix)
                p_s2 = p_jer_solsys + (roots[1] * n_pix)

                # vectors from center of the sun to those points on its surface
                s1 = p_s1 - p_sun
                s2 = p_s2 - p_sun

                # plot if this either point on the sun is visible to the observer
                if (s1.dot(p_s1 - p_jer_solsys) < 0 or s2.dot(p_s2 - p_jer_solsys)) < 0: # and np.dot(n_pix, normalize(p_jer_ecl)) > 0:
                    imageData[N - j, i + N] = [255, 255, 255]

                    # check if this is the center point - and if so, mark it
                    if not trail_marked:

                        # lambda_x = (x_c-x_p) / n_x
                        # lambda_y = (y_c-y_p) / n_y
                        # lambda_z = (z_c-z_p) / n_z
                        # if lambda_x == lambda_y and lambda_x == lambda_z:
                        #     trail.append((N - j, i + N))
                        #     trail_marked = True

                        c = x_p ** 2 + x_c ** 2 - (2 * x_p * x_c) + y_p ** 2 + y_c ** 2 - \
                            (2 * y_p * y_c) + z_p ** 2 + z_c ** 2 - (2 * z_p * z_c) - (sun_rad/10)**2
                        discrim = b ** 2 - (4 * a * c)
                        if discrim >= 0:
                            trail.append((N - j, i + N))
                            # trail.append((N - j + 2, i + N - 2))
                            # trail.append((N - j - 2, i + N + 2))
                            trail.append((N - j + 1, i + N + 1))
                            trail.append((N - j - 1, i + N - 1))
                            trail_marked = True

            else:
                # if this point on the sky is visible
                if np.dot(n_pix, normalize(p_jer_ecl)) > 0:
                    imageData[N - j, i + N] = [0, 0, 100]

    # mark trail from previous snapshots
    if trail_marked: items_to_exclude = 3
    else:            items_to_exclude = 0
    for i in range(len(trail) - items_to_exclude):
        imageData[trail[i]] = [255, 255, 255]

    return imageData, trail


# returns array consisting of the angles between the sun, earth, and moon
def relative_angles(s, e, m):
    E_M = e - m  # Earth - Moon
    E_S = e - s  # Earth - Sun
    M_S = m - s  # Moon - Sun

    theta1 = np.arccos(float((E_S.dot(M_S)) / (length(E_S) * length(M_S))))  # angle between E-S and M-S
    theta2 = np.arccos(float((E_S.dot(E_M)) / (length(E_S) * length(E_M))))  # angle between E-S and E-M
    theta3 = np.arccos(float((-E_M.dot(M_S)) / (length(E_M) * length(M_S))))  # angle between E-S and E-M
    return np.array([theta1, theta2, theta3])


# test function, checks that interpolated positions using yearly Horizons data agree with the daily position data
def test_positions(year):

    positionsArr = get_positions(year)
    dailyPositionsArray = get_data()

    # starting at Jan 1, get positions at each day of the year
    cur_mjd = mjd(1, 1, year, (0, 0))

    num_days = 100

    expected_positions = np.zeros((num_days, 3, 3))    # as per Horizons daily data
    actual_positions = np.zeros((num_days, 3, 3))      # as per interpolation from yearly Horizons data
    expected_angles = np.zeros((num_days, 3))
    actual_angles = np.zeros((num_days, 3))

    for i in range(num_days):
        expected = get_cur_positions(i, dailyPositionsArray)
        actual = np.array([get_pos(cur_mjd, 's', positionsArr), get_pos(cur_mjd, 'e', positionsArr), get_pos(cur_mjd, 'm', positionsArr)])

        for j in range(3):
            for k in range(3):
                expected_positions[i][j][k] = expected[j][k]
                actual_positions[i][j][k] = actual[j][k]

        ## compare angles between sun, earth, and moon as obtained from the two position lists
        s1, e1, m1 = expected   # from expected data
        s2, e2, m2 = actual     # from actual data

        for j in range(3):
            expected_angles[i][j] = relative_angles(s1, e1, m1)[j]
            actual_angles[i][j] = relative_angles(s2, e2, m2)[j]

        cur_mjd += 1   # move to next day

    ## Plot:
    # # Earth-Sun-Moon Angle
    x = np.array([i for i in range(len(expected_positions))])
    # y_exp = np.array([expected_angles[i][0] for i in range(len(expected_angles))])
    # y_act = np.array([actual_angles[i][0] for i in range(len(actual_angles))])
    # plt.scatter(x, y_act-y_exp)
    # # plt.scatter(x, y_act)
    # plt.title("Earth-Sun-Moon Angle")
    # plt.show()
    #
    # Moon-Earth-Sun Angle
    y_exp = np.array([expected_angles[i][1] for i in range(len(expected_angles))])
    y_act = np.array([actual_angles[i][1] for i in range(len(actual_angles))])
    plt.scatter(x, y_act-y_exp)
    # plt.scatter(x, y_act)
    plt.title("Moon-Earth-Sun Angle")
    plt.show()

    # # Earth-Moon-Sun Angle
    # y_exp = np.array([expected_angles[i][2] for i in range(len(expected_angles))])
    # y_act = np.array([actual_angles[i][2] for i in range(len(actual_angles))])
    # plt.scatter(x, y_act-y_exp)
    # # plt.scatter(x, y_act)
    # plt.title("Earth-Moon-Sun Angle")
    # plt.show()

    # Distance between Moon and Earth
    t = np.array([i for i in range(len(expected_positions))])
    p_exp = np.array([expected_positions[i][2]-expected_positions[i][1] for i in range(len(expected_positions))])
    p_act = np.array([actual_positions[i][2]-actual_positions[i][1] for i in range(len(expected_positions))])

    # really should use cartesian distance / find length of x,y,z difference vector
    xyz_diffs = p_act - p_exp
    x_diffs = np.array([xyz_diffs[i][0] for i in range(len(xyz_diffs))])
    y_diffs = np.array([xyz_diffs[i][1] for i in range(len(xyz_diffs))])
    z_diffs = np.array([xyz_diffs[i][2] for i in range(len(xyz_diffs))])
    diff_magnitude = np.array([length(xyz_diffs[i]) for i in range(len(xyz_diffs))])

    # plt.scatter(t, x_diffs)
    # plt.title("Earth-Moon x-Difference")
    # plt.show()
    # plt.scatter(t, y_diffs)
    # plt.title("Earth-Moon y-Difference")
    # plt.show()
    # plt.scatter(t, z_diffs)
    # plt.title("Earth-Moon z-Difference")
    # plt.show()
    plt.scatter(t, diff_magnitude)
    plt.title("Earth-Moon Distance")
    plt.show()

    moon_period(t, expected_angles, actual_angles)


def moon_period(t, expected_angles, actual_angles):
    # find points where sun-earth-moon angle is at a minimum (i.e. time of the molad)
    # use np.diff() to compute moon-period using actual_angles

    molad_times = []

    # extract only sun-earth-moon angles (actual and expected)
    exp_earth_angles = np.array([expected_angles[i][1] for i in range(len(expected_angles))])
    act_earth_angles = np.array([actual_angles[i][1] for i in range(len(actual_angles))])

    # lists of differences between each angle in expected and actual angle lists
    exp_angles_dif = np.diff(exp_earth_angles)
    act_angles_dif = np.diff(act_earth_angles)

    # identify local minima of list of expected angles --> days on which molad occurs
    for i in range(len(exp_angles_dif)-1):
        if exp_angles_dif[i] <= 0 and exp_angles_dif[i+1] >= 0:
            molad_times.append(t[i+1])

    # number of days between each molad
    month_lens = [molad_times[i]-molad_times[i-1] for i in range(1, len(molad_times))]
    return np.mean(month_lens)


# given a mjd, return angular distance between sun and moon at that time
def cur_sun_moon_dist(mjd_cur, positionsArr):
    p_sun = get_pos(mjd_cur, 's', positionsArr)
    p_moon = get_pos(mjd_cur, 'm', positionsArr)
    p_earth = get_pos(mjd_cur, 'e', positionsArr)

    earth_to_sun = p_sun - p_earth
    earth_to_moon = p_moon - p_earth

    return angular_dist(earth_to_sun, earth_to_moon)

#
# def cur_sun_earth_angle(mjd_cur, positionsArr):
#     p_sun = get_pos(mjd_cur, 's', positionsArr)
#     p_earth = get_pos(mjd_cur, 'e', positionsArr)
#
#     earth_to_sun = p_sun - p_earth
#
#     return angular_dist(p_sun, p_earth)


# returns the mjd of the molad of Tishrei for the given Hebrew year
def molad_date(year):

    # determine molad of first month of given year
    molad = molad_determination(1, year, acc_days=True)
    days, hours, chalakim = molad
    hours += chalakim / 1080
    days += hours / 24

    # calculate the day number (since creation) of the mjd epoch and subtract it from the molad day to get its mjd...
    mjd_epoch_day_num = 2052004.25 + (2 / 24)
    mjd_molad = days - mjd_epoch_day_num

    return mjd_molad


# sky snapshots at each molad over the course of 'total_yrs' years
def molad_snapshots(start_yr, total_yrs, draw=True, plot=False):

    day_of_week, greg_day, greg_month, greg_year = hebrew_to_greg(1, 1, start_yr)

    # build up positionsArr to contain necessary years
    positionsArr = []
    for yr in range(greg_year-1, greg_year + total_yrs):
        positionsArr += get_positions(yr)[:-1]
    positionsArr += get_positions(greg_year + total_yrs)

    # mjd of molad Tishrei of start year
    mjd_molad = molad_date(start_yr)

    molad_mjds = []
    angular_dists = []

    sun_x_list = []
    sun_y_list = []
    moon_x_list = []
    moon_y_list = []

    for month in range(total_yrs * 12):
        print(month)

        # set up
        n_jer_ecl, n_jer_equ, n_sun_jer_ecl, n_sun_jer_equ = zenithAndSunVectors(mjd_molad, phi_ny, theta_ny, positionsArr)
        obs_axes = observer_axes(n_jer_equ)
        cam_axes = camera_axes(mjd_molad, n_jer_ecl * earth_rad, positionsArr, obs_axes)

        if draw:
            # draw the sun and moon at this month's molad
            N = 50
            D_fov = 0.05  # field of view, in radians
            imageData = draw_sun2(mjd_molad, cam_axes, n_jer_ecl * earth_rad, positionsArr, N, D_fov, month)
            draw_moon2(mjd_molad, cam_axes, n_jer_ecl * earth_rad, positionsArr, imageData, N, D_fov, month)

        if plot:
            n_cam = cam_axes[2]
            p_moon = get_pos(mjd_molad, 'm', positionsArr)
            p_earth = get_pos(mjd_molad, 'e', positionsArr)
            p_sun = get_pos(mjd_molad, 's', positionsArr)
            p_moon_earth_ecl = p_moon - p_earth
            n_moon_jer_ecl = n_body_jer(n_jer_ecl, p_moon_earth_ecl)

            # sun and moon vectors in camera coords
            sun_x, sun_y = n_camera_coords(n_sun_jer_ecl, cam_axes)
            moon_x, moon_y = n_camera_coords(n_moon_jer_ecl, cam_axes)
            sun_x_list.append(sun_x)
            sun_y_list.append(sun_y)
            moon_x_list.append(moon_x)
            moon_y_list.append(moon_y)

        else:
            # or, determine the angle between the sun and moon
            dist = cur_sun_moon_dist(mjd_molad, positionsArr)
            print(mjd_molad, dist)

            molad_mjds.append(mjd_molad)
            angular_dists.append(dist)


        ## the following code PLOTS the sun and moon on background ra-dec lines (zoomed-out image)
        # f = plt.figure()
        # plt.gca().set_aspect('equal', adjustable='box')
        #
        # RA_Dec_lines(cam_axes, n_jer_equ * earth_rad, imageData)
        # plot_sun_moon(mjd_molad, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr)
        #
        # plt.title("MJD:" + str(mjd_molad))
        # plt.xlabel("Right", fontfamily="times new roman")
        # plt.ylabel("Up", fontfamily="times new roman")
        #
        # # plt.show()
        # f.savefig("./MoladSnapshotsRef/" + str(month))
        # plt.close(f)

        # move to next molad
        mjd_molad += avg_month_len

    plt.scatter(sun_x_list, sun_y_list, color='yellow')
    plt.scatter(moon_x_list, moon_y_list, color='grey')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


# given a start time, return hourly list of times and corresponding list of sun-moon angular dists
def sun_moon_dist(greg_year, start_mjd, total_days, steps_per_day):

    # build up positionsArr to contain necessary years
    positionsArr = get_positions(greg_year)[:-1] + get_positions(greg_year+1)
    mjd_list = []
    dist_list = []

    mjd_cur = start_mjd
    for hour in range(total_days * steps_per_day):

        # determine the angular distance between the sun and moon
        dist = cur_sun_moon_dist(mjd_cur, positionsArr)

        mjd_list.append(mjd_cur)
        dist_list.append(dist)

        mjd_cur += 1/24  # adjust mjd

    return mjd_list, dist_list


# expects ~3-day long list of times (hourly) as well as corresponding ang_dist list
def act_conjunction(days, ang_dist):
    index = 1
    mjd_cur = days[index]

    while mjd_cur < days[-1]:
        x_list = days[index-1:index+2]
        y_list = ang_dist[index-1:index+2]

        # fit order 2 polynomial to the three time-dist points for the prev, cur, and following hours
        p_X = interpolation(x_list, y_list, 2)
        conjunction_time = critical_point(p_X, mjd_cur, interval=1/24)

        # if we found a minimum sometime within this hour
        if conjunction_time: return float(conjunction_time)

        mjd_cur += 1 / 24
        index += 1

    print("couldn't find conjunction time near", mjd_cur)
    return -1


# for each month, get expected molad, call actual conjunction time function, subtract and plot the dif
def plot_avg_molad_dif(start_yr, total_yrs):
    x = []
    diffs = []
    day_of_week, greg_day, greg_month, greg_year = hebrew_to_greg(1, 1, start_yr)

    mjd_greg_newyear = mjd(1, 1, greg_year - 1, (0, 0))

    # mjd of molad Tishrei of given Hebrew year
    mjd_avg_molad = molad_date(start_yr)

    # mjd of actual molad
    mjd_list, dist_list = sun_moon_dist(greg_year-1, mjd_avg_molad-3, 6, 24)
    mjd_act_molad = act_conjunction(mjd_list, dist_list)

    # for each Hebrew month, add average month len to get average molad
    # and then get actual molad time and plot dif
    for month in range(total_yrs*12):
        print(month)
        if mjd_act_molad != -1:
            x.append(month)
            diffs.append(mjd_avg_molad-mjd_act_molad)

        # move to next month
        mjd_avg_molad += avg_month_len

        # check if we're up to the next greg year
        if isLeap(greg_year): yr_len = 366
        else:                 yr_len = 365
        if mjd_avg_molad - mjd_greg_newyear > yr_len:
            greg_year += 1
            mjd_greg_newyear += yr_len

        mjd_list, dist_list = sun_moon_dist(greg_year-1, mjd_avg_molad - 3, 6, 24)
        mjd_act_molad = act_conjunction(mjd_list, dist_list)

    plt.scatter(x, diffs)
    plt.xlabel('month number')
    plt.ylabel('average molad - actual molad (days)')
    plt.show()


def ecliptic_np_angles(mjd_cur, total_days, time_step, positionsArr):

    times = []  # mjd
    angles = []

    # North Pole unit vector
    n_np = equatorial_to_ecliptic([0, 0, 1], mjd_cur)

    num_iterations = int(total_days/time_step)

    # for each hour, get angle between sun-earth vector and North Pole vector
    for i in range(num_iterations):
        p_sun = get_pos(mjd_cur, 's', positionsArr)
        p_earth = get_pos(mjd_cur, 'e', positionsArr)
        p_earth_sun = p_sun - p_earth

        angle = angular_dist(n_np, p_earth_sun)
        times.append(mjd_cur)
        angles.append(angle)
        # print(mjd_cur, angle-(np.pi/2))
        mjd_cur += time_step

    return times, angles


# interpolate using every 3 consecutive hours and corresponding angles determine when angle == 90 degrees
def get_equinox_time(mjd_cur, total_days, time_step, positionsArr):

    times, angles = ecliptic_np_angles(mjd_cur, total_days, time_step, positionsArr)

    index = 1
    mjd_cur = times[index]

    while mjd_cur < times[-1]:
        x_list = times[index - 1:index + 2]
        y_list = angles[index - 1:index + 2]

        # fit order 2 polynomial to the three time-dist points for the prev, cur, and following hours
        p_X = interpolation(x_list, y_list, 2) - (np.pi/2)

        X = symbols('X')
        roots = solveset(p_X, X)
        for r in roots:
            if r.is_real and mjd_cur <= r <= mjd_cur + 1/24:  # if the time falls within the given hour
                return r

        index += 1
        mjd_cur = times[index]

    print("couldn't find equinox time near", mjd_cur)
    return -1


def get_solstice_time(mjd_cur, total_days, time_step, positionsArr):
    times, angles = ecliptic_np_angles(mjd_cur, total_days, time_step, positionsArr)

    index = 1
    mjd_cur = times[index]

    while mjd_cur < times[-1]:
        x_list = times[index - 1:index + 2]
        y_list = angles[index - 1:index + 2]

        # fit order 2 polynomial to the three time-dist points for the prev, cur, and following hours
        p_X = interpolation(x_list, y_list, 2)
        solstice_time = critical_point(p_X, mjd_cur, interval=1/24)

        # if we found a min or max sometime within this hour
        if solstice_time: return float(solstice_time)

        index += 1
        mjd_cur = times[index]

    print("couldn't find solstice time near", mjd_cur)
    return -1


# get position of sun at moment of equinox/solstice
def sun_pos_at_tekufah(mjd_cur, positionsArr, tekufah_type):

    # positionsArr = get_positions(2000)
    # mjd_cur = mjd(20, 3, 2000, (0, 0))

    if tekufah_type == 'equinox': tekufah_mjd = get_equinox_time(mjd_cur, 1, 1/24, positionsArr)
    else:                         tekufah_mjd = get_solstice_time(mjd_cur, 1, 1/24, positionsArr)

    p_sun = get_pos(tekufah_mjd, 's', positionsArr)
    p_earth = get_pos(tekufah_mjd, 'e', positionsArr)
    p_earth_sun = p_sun - p_earth
    sun_location = np.arctan2(float(p_earth_sun[1]), float(p_earth_sun[0]))

    return sun_location


def season_times(start_tekufah, start_yr, end_yr):

    # build up positionsArr for the requisite number of years
    positionsArr = get_positions(start_yr)
    for i in range(start_yr+1, end_yr):
        positionsArr += get_positions(i)[1:]

    df = pd.DataFrame(columns=['date', 'mjd', 'days since prev equinox/solstice', 'days since this equinox/solstice last year'])

    d, m, y = start_tekufah
    mjd_cur = mjd(d, m, y, (0, 0))
    end_mjd = mjd(0, 0, end_yr, (0, 0))

    prev_season_mjd = None
    prev_seasons_q = Queue(maxsize=4)
    this_season_last_yr_mjd = None

    prev = 'S'
    i = 0
    while mjd_cur <= end_mjd:
        # are we up to an equinox or solstice?
        if prev == 'S':
            tekufah_mjd = get_equinox_time(mjd_cur, 10, 1/24, positionsArr)
            prev = 'E'
        else:
            tekufah_mjd = get_solstice_time(mjd_cur, 10, 1/24, positionsArr)
            prev = 'S'

        if tekufah_mjd == -1: return

        # if queue is full, we can get the mjd of this equinox/solstice last year
        if prev_seasons_q.full():
            this_season_last_yr_mjd = prev_seasons_q.get()

        time_since_prev = tekufah_mjd - prev_season_mjd if prev_season_mjd else None
        time_since_prev_cur_season = tekufah_mjd - this_season_last_yr_mjd if this_season_last_yr_mjd else None

        # function to convert mjd into dd-mm-yyyy and hh:mm:ss
        date = mjd_to_date(tekufah_mjd)

        # add entry to dataframe for this equinox/solstice and corresponding info
        # df = pd.concat([pd.DataFrame([[date, tekufah_mjd, time_since_prev, time_since_prev_cur_season]], columns=df.columns), df], ignore_index=True)
        row = date, tekufah_mjd, time_since_prev, time_since_prev_cur_season
        df.loc[len(df)] = row

        print(date, tekufah_mjd, time_since_prev, time_since_prev_cur_season)
        # note this mjd
        prev_seasons_q.put(tekufah_mjd)
        prev_season_mjd = tekufah_mjd

        i += 1
        mjd_cur = tekufah_mjd//1 + 88

    return df


def mjd_to_date(mjd_cur):

    # count mjd_cur number of days from epoch
    epoch = (17, 11, 1858, (0, 0))  # 17 November, 1858, 12 am (midnight)

    day_cur, month_cur, yr_cur = epoch[:3]
    days_counted = 0

    # move to first of next month for cleaner calculations
    leap = isLeap(yr_cur)
    month_length = month_len(month_cur, leap)
    days_counted += month_length - day_cur

    # if we've reached or passed the number of days to count, return date info
    if days_counted >= mjd_cur:
        day_cur += int(mjd_cur)
        return day_cur, month_cur, yr_cur, clean_time((mjd_cur % 1)*24)

    # otherwise, increment month count and keep going
    month_cur += 1

    while True:
        leap = isLeap(yr_cur)

        # while we're still in this year
        while month_cur <= 12:
            month_length = month_len(month_cur, leap)
            days_counted += month_length

            # if we've reached or passed the number of days to count, return date info
            if days_counted >= mjd_cur:
                day_cur = int(month_length - (days_counted - mjd_cur))
                return day_cur, month_cur, yr_cur, clean_time((mjd_cur % 1)*24)

            month_cur += 1

        # now move to next year
        yr_cur += 1
        month_cur = 1


# calculates solar longitude, lunar longitude, and lunar latitude once a day from start_year til end_year (not inclusive)
# writes datapoints to a file
def write_sun_moon_longitude_and_latitude(start_year, end_year, out_file_path):

    # build up positionsArr
    positionsArr = get_positions(start_year)
    for i in range(start_year+1, end_year):
        positionsArr += get_positions(i)[1:]

    out_file = open(out_file_path, 'w')
    out_file.close()
    out_file = open(out_file_path, 'a')

    mjd_cur = mjd(1, 1, start_year, (0, 0))
    mjd_end = mjd(1, 1, end_year, (0, 0))

    while mjd_cur < mjd_end:
        p_sun = get_pos(mjd_cur, 's', positionsArr)
        p_earth = get_pos(mjd_cur, 'e', positionsArr)
        p_moon = get_pos(mjd_cur, 'm', positionsArr)
        p_earth_sun = p_sun - p_earth
        p_earth_moon = p_moon - p_earth

        # solar longitude
        sun_longitude = np.arctan2(float(p_earth_sun[1]), float(p_earth_sun[0])) * 180/np.pi
        if sun_longitude < 0: sun_longitude += 360

        # lunar longitude and latitude
        moon_longitude = np.arctan2(float(p_earth_moon[1]), float(p_earth_moon[0])) * 180 / np.pi
        moon_latitude = np.arcsin(float(p_earth_moon[2])/length(p_earth_moon)) * 180 / np.pi
        if moon_longitude < 0: moon_longitude += 360

        cur_day_data = f'{mjd_cur} {sun_longitude} {moon_longitude} {moon_latitude}\n'
        out_file.write(cur_day_data)
        mjd_cur += 1

    out_file.close()


def cur_angle2(mjd_cur, epoch_angle, degrees_per_day):

    days_since_epoch = mjd_cur - mjd_Rambam_epoch
    angle_traversed = days_since_epoch * degrees_per_day
    cur_angle = (epoch_angle + angle_traversed) % 360

    return cur_angle


def sun_longitude_Rambam2(mjd_cur, cur_model):

    cur_mean_longitude = cur_angle2(mjd_cur, cur_model.mean_sun_long_at_epoch, cur_model.mean_sun_long_rate) * np.pi/180
    cur_aphelion = cur_angle2(mjd_cur, cur_model.aphelion_at_epoch, cur_model.aphelion_rate) * np.pi/180
    alpha = cur_mean_longitude - cur_aphelion
    if alpha < 0: alpha += 2 * np.pi

    gamma = np.arctan2(np.sin(alpha), np.cos(alpha) + eccentricity)

    return ((gamma + cur_aphelion) * 180 / np.pi) % 360


def get_tosefes(double_elongation):
    tosefes = 0
    if 6 <= double_elongation <= 11:
        tosefes = 1
    elif 12 <= double_elongation <= 18:
        tosefes = 2
    elif 19 <= double_elongation <= 24:
        tosefes = 3
    elif 25 <= double_elongation <= 31:
        tosefes = 4
    elif 32 <= double_elongation <= 38:
        tosefes = 5
    elif 39 <= double_elongation <= 45:
        tosefes = 6
    elif 46 <= double_elongation <= 51:
        tosefes = 7
    elif 52 <= double_elongation <= 59:
        tosefes = 8
    elif 60 <= double_elongation <= 63:
        tosefes = 9

    # e_1 = 10.3167
    # R = 49.6833
    # R_prime = e_1 * np.cos(hamerchak_hakaful) + np.sqrt(R**2 + (e_1**2 * (np.cos(hamerchak_hakaful)**2 - 1)))
    # R_double_prime = np.sqrt(R_prime**2 + e_1**2 - (2 * e_1 * R_prime * np.cos(np.pi-hamerchak_hakaful)))
    # tosefes_2 = np.arccos((R_prime**2 + R_double_prime**2 - e_1**2) / (2 * R_prime * R_double_prime))
    # print(tosefes, tosefes_2)
    return tosefes


def moon_longitude_Rambam2(mjd_cur, cur_model):

    mean_sun_long = cur_angle2(mjd_cur, cur_model.mean_sun_long_at_epoch, cur_model.mean_sun_long_rate)
    mean_moon_long = cur_angle2(mjd_cur, cur_model.mean_moon_long_at_epoch, cur_model.mean_moon_long_rate)
    mean_moon_path = cur_angle2(mjd_cur, cur_model.mean_moon_path_at_epoch, cur_model.mean_moon_path_rate)

    elongation = mean_moon_long - mean_sun_long
    double_elongation = 2 * elongation

    tosefes = get_tosefes(double_elongation)

    correct_moon_course = mean_moon_path + tosefes

    # look up m'nas hamaslul hanachon given maslul hanachon (interpolate)
    mnas_hamaslul_x = [i for i in range(0, 181, 10)]
    mnas_hamaslul_y = np.array([0, 50/60, 1+(38/60), 2+(24/60), 3+(6/60), 3+(44/60), 4+(16/60), 4+(41/60), 5, 5+(5/60),
                        5+(8/60), 4+(59/60), 4+(40/60), 4+(11/60), 3+(33/60), 2+(48/60), 1+(56/60), 59/60, 0]) * cur_model.e1

    if correct_moon_course > 180:
        mnas_hamaslul = np.interp(360 - correct_moon_course, mnas_hamaslul_x, mnas_hamaslul_y)
    else:
        mnas_hamaslul = -np.interp(correct_moon_course, mnas_hamaslul_x, mnas_hamaslul_y)

    true_moon_longitude = mean_moon_long + mnas_hamaslul

    if true_moon_longitude < 0: true_moon_longitude += 360

    return true_moon_longitude


def moon_latitude_Rambam2(mjd_cur, true_moon_longitude, cur_model):
    mean_head_long = cur_angle2(mjd_cur, cur_model.mean_head_long_at_epoch, cur_model.mean_head_long_rate)

    maslul_harochav = true_moon_longitude - mean_head_long
    if maslul_harochav < 0: maslul_harochav += 360

    direction = 'N'
    if 90 <= maslul_harochav <= 180:
        maslul_harochav = 180 - maslul_harochav

    elif 180 < maslul_harochav <= 270:
        maslul_harochav = maslul_harochav - 180
        direction = 'S'

    elif 270 < maslul_harochav <= 360:
        maslul_harochav = 360 - maslul_harochav
        direction = 'S'

    # use maslul harochav to determine declination (given Rambam's table)
    mnas_hamaslul_x = [i for i in range(0, 91, 10)]
    # mnas_hamaslul_y = np.array([0, 52/60, 1+(43/60), 2+(30/60), 3+(13/60), 3+(50/60), 4+(20/60), 4+(42/60), 4+(55/60), 5]) * (5+14/60)/5
    mnas_hamaslul_y = np.array([0, 52/60, 1+(43/60), 2+(30/60), 3+(13/60), 3+(50/60), 4+(20/60), 4+(42/60), 4+(55/60), 5]) * cur_model.e2

    true_moon_latitude = np.interp(maslul_harochav, mnas_hamaslul_x, mnas_hamaslul_y)
    if direction == 'S': true_moon_latitude = -true_moon_latitude

    return true_moon_latitude


def visible(sun_long, moon_long, moon_lat):
    first_longitude = moon_long - sun_long
    first_latitude = moon_lat



class Model(object):
    def __init__(self, mean_sun_long_at_epoch, mean_sun_long_rate, aphelion_at_epoch, aphelion_rate, mean_moon_long_at_epoch, mean_moon_long_rate, mean_moon_path_at_epoch, mean_moon_path_rate, mean_head_long_at_epoch, mean_head_long_rate, e1, e2):
        self.mean_sun_long_at_epoch = mean_sun_long_at_epoch
        self.mean_sun_long_rate = mean_sun_long_rate
        self.aphelion_at_epoch = aphelion_at_epoch
        self.aphelion_rate = aphelion_rate
        self.mean_moon_long_at_epoch = mean_moon_long_at_epoch
        self.mean_moon_long_rate = mean_moon_long_rate
        self.mean_moon_path_at_epoch = mean_moon_path_at_epoch
        self.mean_moon_path_rate = mean_moon_path_rate
        self.mean_head_long_at_epoch = mean_head_long_at_epoch
        self.mean_head_long_rate = mean_head_long_rate
        self.e1 = e1
        self.e2 = e2

    def update_params(self, alpha):
        for attr in self.__dict__:
            attr = attr*alpha


def select_n_rows(file_path, n):
    with open(file_path, 'r') as f:
        df = pd.read_table(file_path, delimiter=' ', names=['mjd', 'sun_longitude', 'moon_longitude', 'moon_latitude'], index_col=0)
        n_random_rows = df.sample(n, random_state=2).sort_index()
    return n_random_rows


# then, read in file and compare angles with Rambam's values
def compare_angles(df, cur_model=None):

    # add columns for Rambam's angles
    df['sun_long_Rambam'] = ''
    df['moon_long_Rambam'] = ''
    df['moon_lat_Rambam'] = ''

    for index, row in df.iterrows():
        mjd_cur = index

        # calculate angles using Rambam for cur mjd
        # sun_long_Rambam = sun_longitude_Rambam(mjd_cur)
        # moon_long_Rambam = moon_longitude_Rambam(mjd_cur)
        # moon_lat_Rambam = moon_latitude_Rambam(mjd_cur, moon_long_Rambam)

        sun_long_Rambam = sun_longitude_Rambam2(mjd_cur, cur_model)
        moon_long_Rambam = moon_longitude_Rambam2(mjd_cur, cur_model)
        moon_lat_Rambam = moon_latitude_Rambam2(mjd_cur, moon_long_Rambam, cur_model)

        # and add to df
        df.loc[mjd_cur, 'sun_long_Rambam'] = sun_long_Rambam
        df.loc[mjd_cur, 'moon_long_Rambam'] = moon_long_Rambam
        df.loc[mjd_cur, 'moon_lat_Rambam'] = moon_lat_Rambam

    # print(tabulate(df.head(10), headers='keys'))

    # error calculations
    residuals = []
    # for i in range(3):
    for i in range(1):

        col_name = df.columns[i]
        expected = np.array(df[df.columns[i]])
        actual = np.array(df[df.columns[i + 3]])
        cur_residuals = np.arcsin(np.sin((expected - actual).astype(float)*np.pi/180)) * 180/np.pi
        residuals += list(cur_residuals)

        plt.scatter(df.index, cur_residuals)
        plt.title(col_name)
        plt.xlabel('mjd')
        plt.ylabel("residuals (data-model, degrees)")
        # plt.xlim(52500, 54000)
        plt.show()

    chi_squared_error = np.sum(np.array(residuals) ** 2)

    return chi_squared_error


# find values of epoch and rate angles that minimize the chi-squared error for sun/moon longitude and moon latitude
def mcmc(file_path, n):

    # data
    df = select_n_rows(file_path, n)

    # prior parameter values
    mean_sun_long_at_epoch = 7 + 3 / 60 + 32 / 3600
    mean_sun_long_rate = (98 + 33 / 60 + 52.3 / 3600) / 100
    # mean_sun_long_rate = (98 + 33 / 60 + 53 / 3600) / 100
    aphelion_at_epoch = 86 + 45 / 60 + 8 / 3600
    aphelion_rate = (18.6/3600) / 100
    # aphelion_rate = (15/3600) / 100
    mean_moon_long_at_epoch = 31 + 14/60 + 43/3600
    mean_moon_long_rate = 13 + 10/60 + 35/(60**2) + 1/(60**3) + 170/(60**4)
    # mean_moon_long_rate = 13 + 10/60 + 35/(60**2) + 1/(60**3) + 48/(60**4)
    mean_moon_path_at_epoch = 84 + 28/60 + 42/3600
    mean_moon_path_rate = 13 + 3/60 + 54/3600
    mean_head_long_at_epoch = -(180 + 57/60 + 28/3600)
    # mean_head_long_at_epoch = -(180 + 57/60 + 10/3600)
    mean_head_long_rate = -((52 + 57/60 + 23/3600) / 1000)
    e1 = 1
    e2 = (5+14/60)/5

    cur_model = Model(mean_sun_long_at_epoch, mean_sun_long_rate, aphelion_at_epoch, aphelion_rate, mean_moon_long_at_epoch, mean_moon_long_rate, mean_moon_path_at_epoch, mean_moon_path_rate, mean_head_long_at_epoch, mean_head_long_rate, e1, e2)
    chi_squared = compare_angles(df, cur_model)     # compare current model against data
    best_chi_squared = chi_squared
    best_model = cur_model

    accepted_list = []
    # param_ranges = [6, 7, 5, 11, 6, 7, 5, 6, 4, 8, 7, 7]
    # param_ranges = [6, 7, 5, 6, 4, 10, 4, 10, 4, 10, 7, 7]
    param_ranges = [6, 6, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]

    for i in range(100):
        print(i)
        # add/subtract random number (normal distribution) from each parameter value
        new_params = []
        cur_model_param_dict = cur_model.__dict__
        j = -1
        for param, val in cur_model_param_dict.items():
            j += 1
            # # scale alpha differently for each paramater
            scale = 10**(-param_ranges[j])
            alpha = np.random.normal(loc=0, scale=scale)
            new_params.append(val + alpha * val)

        new_model = Model(*new_params)
        new_chi_squared = compare_angles(df, new_model)

        print(chi_squared)
        print(new_chi_squared)
        # print(cur_model.__dict__)
        # print(new_model.__dict__)

        # accept?
        temp = 20
        prob_ratio = np.exp(-(new_chi_squared-chi_squared)/(2*temp))
        z = np.random.uniform(0, 1)
        if prob_ratio > z:
            print('accept')
            cur_model = new_model
            chi_squared = new_chi_squared
            accepted_list.append(new_params)

            # save best fit params
            if chi_squared < best_chi_squared:
                best_chi_squared = chi_squared
                best_model = cur_model
        else:
            print('reject')
        print()
    # range of each param
    # col_ranges = np.ptp(accepted_list, axis=0)
    # est_power = [np.ceil(np.abs(np.log10(col_range))) for col_range in col_ranges]
    # print(col_ranges)
    # print(est_power)


def main():

    year = int(input("Year: "))
    month = int(input("Month: "))
    day = int(input("Day: "))
    hour = int(input("Hour (0-23): "))
    minutes = float(input("Minutes: "))
    location = input("Location ('ny', 'jer', 'texas'): ")

    if location == "ny":    time_zone = 'EST'
    elif location == "jer": time_zone = 'IST'
    elif location == "texas": time_zone = 'CST'

    else:
        print("invalid location")
        return

    # daylight savings time beginning and end for this year in this time zone
    DST_START, DST_END = DST(time_zone, year)

    # mjd of date in question (in GMT)
    mjd_cur = mjd(day, month, year, (hour, minutes))  # modified julian date
    mjd_cur -= time_dif(mjd_cur, time_zone, DST_START, DST_END)/24

    # array containing positions of sun, earth, and moon as returned from SunEarthMoonPosition.py
    positionsArr = get_positions(year)

    # print(sunriseSunset(mjd_cur, time_zone, DST_START, DST_END, positionsArr))

    phi_dif, theta = phiDifAndTheta(time_zone)
    # sky_snapshots(mjd_cur, phi_dif, theta, positionsArr, 2/24, 20*24)
    sky_snapshots2(mjd_cur, phi_dif, theta, positionsArr, 3/24, 20*24)
    # movie("SkyImages 3.10.24", "SummerSolsticeMovie")
    # movie("SkyImages 6.6.24", "SpringEquinoxMovie")


# main()


# test_positions(2023)

# molad_snapshots(5764, 40, draw=False, plot=True)

# mjds, angular_dists = sun_moon_dist(2024, 60407, 5, 24)
# print(act_conjunction(mjds, angular_dists))

# plot_avg_molad_dif(5765, 9)


## Chart comparing calculated equinox/solstice times with actual
# start_tekufah = (20, 3, 2000)
# df = season_times(start_tekufah, start_yr=2000, end_yr=2025)
# print(tabulate(df, headers='keys'))

## Compare position of the sun (phi) using Rambam's calculations against those obtained using Symplectic integrator
# mjd_cur = mjd(15, 8, 2024, (0, 0))
# compare_sun_angles(mjd_cur, 2024)

# mjd_cur = mjd(1, 10, 2024, (19, 0))
# mjd_cur = mjd(28, 4, 1178, (0, 0))
# compare_moon_angles(mjd_cur, 2024)


# # employ mcmc to fit params to Rambam's model
out_file_path = './sun_moon_angles.txt'
# # write_sun_moon_longitude_and_latitude(start_year=2000, end_year=2050, out_file_path=out_file_path)
mcmc(out_file_path, n=1000)

