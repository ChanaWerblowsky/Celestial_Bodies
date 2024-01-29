
import numpy as np
from sympy import *
from sympy import symbols, linsolve, solveset
import matplotlib.pyplot as plt
import cv2
import glob
import re
import os
from DateConversion import year_lists
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
phi_jer = 35.214        # degrees
phi_ny = -74.006        # degrees
theta_jer = 58.232      # degrees
theta_ny = 49.3         # degrees
dayNum_mjdEpoch = 2052005
sunset_angle = 90.8     # degrees
obliquity = 23.43       # degrees
moon_obliquity = 1.5424  # degrees
T_side = 27.321661         # sidereal period of the moon, in days
rotation_phase0 = 217.88   # moon's rotation phase at Epoch 0, in degrees
# rotation_phase0 = rotation_phase0 + 45  # test factor, 1/3/24

# psi_cam = 90   # camera position for ~8AM
phi_cam = 259
# psi_cam = 0.1  # camera position for ~12PM
# phi_cam = 180
psi_cam = 76  # camera position for ~6PM
# phi_cam = 205
canvas_size = 256


# convert to radians
radians_day = degrees_day * np.pi/180
phi_jer = phi_jer * np.pi/180
phi_ny = phi_ny * np.pi/180
theta_jer = theta_jer * np.pi/180
theta_ny = theta_ny * np.pi/180
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
    # print(coeff)

    X, Y = symbols('X Y')

    # get the one tuple in coefficient FiniteSet (the polynomial)
    for tup in coeff:
        polynomial = (tup[0] * X ** 2) + (tup[1] * X) + tup[2]  # - sunset_angle

    return polynomial


# returns the position of the body in question (sun/earth/moon) at the given MJD
def get_pos(MJD, body, positionsArr):

    index_map = {'s': 1, 'e': 2, 'm': 3}
    i = index_map[body]

    # tstart, tstop, Nt = SunEarthMoonPosition.tstart, SunEarthMoonPosition.tstop, SunEarthMoonPosition.Nt
    tstart = positionsArr[0][0]  # mjd of Jan 1 of this year
    tstop = tstart + 366
    Nt = len(positionsArr) + 1

    # determine where the given MJD lies relative to the list of mjd/position data points
    daysRecorded = tstop - tstart
    fractionOfTotal = (MJD - tstart) / daysRecorded
    est_index = int(fractionOfTotal * Nt)  # the index of this MJD in Positions array

    # go to the estimated index and adjust if necessary
    data = positionsArr[est_index]
    # print(data)
    while MJD < data[0]:                        # if MJD is lower than current entry, move backward
        est_index -= 1
        data = positionsArr[est_index]

    while MJD > positionsArr[est_index+1][0]:   # if greater than the following entry, move forward
        est_index += 1
        data = positionsArr[est_index]

    # interpolate to get positions at exact time in question
    t_Vals = []
    pos_Vals = []
    for j in range(-1, 2):  # use the prev interval, this one, and the following one
        t_Vals.append(positionsArr[est_index+j][0]-tstart)
        pos_Vals.append(positionsArr[est_index+j][i])

    # isolate the variables
    x_Vals = [coords[0] for coords in pos_Vals]
    y_Vals = [coords[1] for coords in pos_Vals]
    z_Vals = [coords[2] for coords in pos_Vals]

    # x = np.interp(MJD-tstart, t_Vals, x_Vals)
    # y = np.interp(MJD-tstart, t_Vals, y_Vals)
    # z = np.interp(MJD-tstart, t_Vals, z_Vals)


    X, Y = symbols('X Y')

    ## interpolate x, y, and z coordinates separately:
    # create polynomial for each
    p_X = interpolation(t_Vals, x_Vals, 2)
    p_Y = interpolation(t_Vals, y_Vals, 2)
    p_Z = interpolation(t_Vals, z_Vals, 2)

    # now, plug the given t into each polynomial; solve for x, y, and z
    x = p_X.subs(X, MJD-tstart)
    y = p_Y.subs(X, MJD-tstart)
    z = p_Z.subs(X, MJD-tstart)

    return np.array([x, y, z])


### Utility functions
def normalize(v):
    return v / length(v)


def length(v):
    return np.sqrt(float(v.dot(v)))

# normalized position vector from current Jerusalem position to current position of sun/moon

def n_body_jer(n_jer, p_body):
    p_jer = n_jer * earth_rad
    p_body_jer = p_body - p_jer
    return normalize(p_body_jer)


def ecliptic_to_equatorial(v, eps=obliquity):
    x_eq = v[0]
    y_eq = (np.cos(eps) * v[1]) - (np.sin(eps) * v[2])
    z_eq = (np.sin(eps) * v[1]) + (np.cos(eps) * v[2])

    return np.array([x_eq, y_eq, z_eq])


def equatorial_to_ecliptic(v, eps=obliquity):
    x_ec = v[0]
    y_ec = (np.cos(eps) * v[1]) + (np.sin(eps) * v[2])
    z_ec = -(np.sin(eps) * v[1]) + (np.cos(eps) * v[2])

    return np.array([x_ec, y_ec, z_ec])


def phiDifAndTheta(time_zone):

    # GMT
    phi_dif = 0
    theta = 0

    if time_zone == 'IST':
        phi_dif = phi_jer
        theta = theta_jer
    elif time_zone == 'EST':
        phi_dif = phi_ny
        theta = theta_ny
    # elif time_zone == 'GMT':
    #     phi_dif = 0
    #     theta = 0

    return phi_dif, theta


def DST(time_zone, yr):
    dst_start_month = 3  # March

    if time_zone == 'IST':
        dst_start_day_of_week = (0, 'last')  # really starts on Fri before last Sunday

    elif time_zone == 'EST':
        dst_start_day_of_week = (0, 2)  # starts on second Sunday in March

    # determine MJD of Friday before last Sunday in March of the given year

    # day of week of March 1
    greg_yrs, hebrew_yrs = year_lists()
    hebrew_yrs = None  # garbage collection

    yr_tuple = greg_yrs[yr + 3760]  # year tuple from year-lists specifying MJD and day of week of Jan 1

    # pointer starting at March 1
    if yr_tuple[3]:
        cur = 60  # if leap year
    else:
        cur = 59  # if not leap year

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
    day_of_week_Oct1 = (5 + ((
                                         DST_END - DST_START + 1) % 7)) % 7  # Friday (day of DST start) + remaining days after complete weeks + 1 to get to Oct 1
    if not day_of_week_Oct1:
        date_first_sunday = 1
    else:
        date_first_sunday = 8 - day_of_week_Oct1  # because the first of Oct is included

    # number of full weeks between first Sunday and last day of month (inclusive)
    weeks = (month_len(10, False) - date_first_sunday) // 7

    date_last_sunday = date_first_sunday + (7 * weeks)
    DST_END += date_last_sunday

    return DST_START, DST_END


def time_dif(MJD, time_zone, DST_START, DST_END):
    # DST_START, DST_END = DST(MJD, time_zone, yr)

    if time_zone == 'IST':
        standard_time_dif = 2
        daylight_time_dif = 3

    elif time_zone == 'EST':
        standard_time_dif = -5
        daylight_time_dif = -4

    elif time_zone == 'GMT':
        return 0

    # check if date falls out in DST:
    if DST_START <= MJD < DST_END:
        return daylight_time_dif

    else:
        return standard_time_dif


# time in GMT
def mjd(dd, mm, yy, tt):
    epoch = (17, 11, 1858, (0, 0))  # 17 November, 1858, 12 am (midnight)

    days = 0  # days since epoch

    if (dd, mm, yy, tt) == epoch:
        return days  # 0

    if yy > epoch[2]:  # this year is past the epoch year
        start = epoch
        end = (dd, mm, yy, tt)

    elif yy < epoch[2]:  # this year precedes the epoch year
        start = (dd, mm, yy, tt)
        end = epoch

    # complete days remaining from start month (not including start day)
    days += (month_len(start[1], isLeap(start[2])) - epoch[0])

    # complete months remaining from start year
    for month in range(start[1] + 1, 13):
        days += month_len(month, False)  # 1858 (epoch year) was not a leap year

    # complete years between start and end
    # loop through complete years between epoch year and current year, aggregating their days
    for year in range(start[2] + 1, end[2]):
        # is it a leap year?
        if isLeap(year):
            days += 366
        else:
            days += 365

    # complete months from end year
    if isLeap(end[2]):
        leap = True
    else:
        leap = False

    for month in range(1, end[1]):
        days += month_len(month, leap)

    # complete days in end month (not including end day)
    days += end[0] - 1

    # if date in question is past epoch, must add complete start day (since epoch day is a complete day) and the time that has elapsed in the end day
    if start == epoch:
        days += ((end[3][0] + (end[3][1] / 60)) / 24) + 1  # (end[3] = tt)

    # if, however, given date precedes the epoch, only add remainder of start day
    elif end == epoch:
        days += ((24 - start[3][0]) + (60 - start[3][1] / 60)) / 24

    return days


# given number of days since equinox, time of equinox in GMT, and phi difference between Greenwich and location in question,
# determine phi value of that location
def phi_val(days_since_eq, phi_dif):

    eq_hours = int(equinox_time[0])  # am or pm?
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
def zenithAndSunVectors(MJD, phi_dif, theta, positionsArr):
    # days since spring equinox
    days_since_equinox = MJD - mjd_equinox

    # location of Jerusalem in terms of theta and phi
    phi = phi_val(days_since_equinox, phi_dif)

    # Jerusalem's normalized position vector
    n_jer_equ = n_spherical_to_cartesian(theta, phi)            # equatorial coords
    n_jer_ecl = equatorial_to_ecliptic(n_jer_equ)    # equatorial coords

    p_sun = get_pos(MJD, 's', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)
    p_sun = (p_sun-p_earth) * au_cm

    # sun's position relative to Jerusalem
    n_sun_jer = n_body_jer(n_jer_ecl, p_sun)                      # ecliptic coords
    n_sun_jer_equ = ecliptic_to_equatorial(n_sun_jer)  # equatorial coords

    return n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ


# angle between zenith of location in question and the vector from that location to the sun
def psi_val(MJD, phi_dif, theta, positionsArr):

    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD, phi_dif, theta, positionsArr)

    # dot product of n_jer in ecliptic coords and n_sun from jerusalem
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
def n_observer_coords(n, observer_axes):

    e_N, e_E, e_Z = observer_axes
    E_component = n.dot(e_E)
    N_component = n.dot(e_N)
    Z_component = n.dot(e_Z)

    return np.array([E_component, N_component, Z_component])


# convert from observer coords to psi and phi (angles of the sun relative to the observer)
def n_observer_angles(E_coord, N_coord, Z_coord):
    psi = acos(Z_coord)            # latitude equivalent
    phi = atan2(E_coord, N_coord)  # longitude equivalent

    if phi < 0: phi = phi + 2*np.pi

    # phi = float(phi) * 180/numpy.pi
    # psi = 90 - (psi * 180/numpy.pi)

    return psi, phi


# convert psi and phi angles to position vector in equatorial coords
def n_from_observer_angles(psi, phi, observer_axes):

    cartesian_obs_coords = n_spherical_to_cartesian(psi, phi)       # in observer coordinates
    cartesian_equ_coords = cartesian_obs_coords.dot(observer_axes)  # now in equatorial coordinates

    return cartesian_equ_coords


def psi_vs_time(MJD, time_zone, dst_start, dst_end, positionsArr, time_step=0.01):

    phi_dif, theta = phiDifAndTheta(time_zone)

    # plot psi vs. time
    x = []
    y = []

    for day in range(1):

        if day:  MJD += 1

        time_difference = time_dif(MJD, time_zone, dst_start, dst_end)
        time = 0
        # time = 0-time_difference

        while time < 24:
            # find psi value (angle between epoch and sun position vector) at this point in time
            mjd_frac = MJD + time / 24
            psi = psi_val(mjd_frac, phi_dif, theta, positionsArr)

            # x.append(time + time_dif(MJD, time_zone, yr) + (24*day))
            acc_time = time + (24 * day) + time_difference
            x.append(acc_time)
            y.append(psi)

            time += time_step


    # plt.xticks(np.arange(0, 25, 1))
    # plt.plot(x, y, color='blue', marker='.', )
    # plt.show()

    # return (psi2 - psi1) / (t2 - t1)
    return x, y


def solar_noon(polynomial, hr):

    X = symbols('X')

    # take derivative of given polynomial; check if x value lies within given hour-range
    derivative = diff(polynomial, X)

    critical_points = solveset(derivative, X)

    for cp in critical_points:
        if hr <= cp <= hr + 1:  # if the time falls within the given hour
            return cp

    return ''


def sunriseSunset(MJD, time_zone, dst_start, dst_end, positionsArr):

    # get list of psi-values for each hour during this day
    x, y = psi_vs_time(MJD, time_zone, dst_start, dst_end, positionsArr, 1)

    # construct polynomial of psi vs. time for each set of 4 hours; check at which times the angle is 91 degrees and at what time it's at a max
    noon = 0
    ans = []
    for i in range(3, 21):

        # get polynomial defined by the angles this hour, the previous hour, and the two following hours
        hours = x[i - 3:i + 1]  # slice of larger hour list
        angles = y[i - 3:i + 1]  # slice of larger angle list

        #### System of Equations
        system = []
        a, b, c, d = symbols('a b c d')

        # create system of equations using given 4 (x,y) pairs
        for j in range(4):
            equation = ((hours[j]**3)*a) + ((hours[j]**2)*b) + ((hours[j])*c) + d - angles[j]
            system.append(equation)

        # list of coefficients of solution to linear system
        coeff = linsolve(system, [a, b, c, d])

        X, Y = symbols('X Y')

        # get the one tuple in coefficient FiniteSet (the polynomial)
        for tup in coeff:
            polynomial = (tup[0] * X**3) + (tup[1] * X**2) + (tup[2] * X) + tup[3] - sunset_angle
        #####################

        # get time of solar noon
        potential_noon = solar_noon(polynomial+sunset_angle, i)

        if potential_noon:
            noon = float(potential_noon)

        # for each of the roots, check whether it's real and it falls within the given hour
        roots = solveset(polynomial, X)
        for r in roots:
            if r.is_real and i <= r <= i + 1:   # if the time falls within the given hour
                ans.append(float(r))

    # return times of sunrise, solar noon, and sunset - in both cleaned-up and raw time
    return ('sunrise: ' + clean_time(ans[0])  + '; solar noon: ' + clean_time(noon) + '; sunset: ' + clean_time(ans[1])), (ans[0], noon, ans[1])


def clean_time(hours):
    # hour_frac = hours % 1

    # min = str(int(hour_frac * 60))
    min = (hours % 1) * 60
    sec = str(int((min % 1) * 60))

    min = str(int(min))

    if len(min) == 1: min = '0' + min
    if len(sec) == 1: sec = '0' + sec

    return str(int(hours)) + ':' + min + ':' + sec


def earliest_sunset(yr, MJD, time_zone, dst_start, dst_end, positionsArr):

    # get mjd of first day of given year
    mjd_cur = mjd(1, 1, yr, (0, 0))

    earliest = 24  # hours
    earliestDay = None

    if isLeap(yr):
        yr_len = 366
    else:
        yr_len = 365

    for day in range(50):
    # for day in range(yr_len):

        # find time of sunset
        # sunset = sunsetTime(mjd_cur)

        sunset = sunriseSunset(mjd_cur, time_zone, dst_start, dst_end, positionsArr)[1][2]

        if sunset <= earliest:
            earliest = sunset
            earliestDay = day

        # move to next day
        mjd_cur += 1

    return clean_time(earliest), earliestDay


def timesAndAngle(d, m, y, MJD, time_zone, DST_START, DST_END, positionsArr):

    phi_dif, theta = phiDifAndTheta(time_zone)

    times = sunriseSunset(MJD, time_zone, DST_START, DST_END, positionsArr)

    # time of sunrise and sunset (via interpolations)
    print(str(d) + '/' + str(m) + '/' + str(y) + ': ' + times[0])

    #
    MJD_sunrise = MJD + (times[1][0]/24) - 2/24
    MJD_noon = MJD + (times[1][1]/24) - 2/24
    MJD_sunset = MJD + (times[1][2]/24) - 2/24

    ## phi at SUNRISE in observer coordinates:
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_sunrise, phi_dif, theta, positionsArr)

    # set up observer axes
    axes = observer_axes(n_jer_equ)

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = n_observer_angles(coords[0], coords[1], coords[2])

    ### phi at SUNSET
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_sunset, phi_dif, theta, positionsArr)
    # observer axes
    axes = observer_axes(n_jer_equ)

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = n_observer_angles(coords[0], coords[1], coords[2])

    ### phi at NOON
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_noon, phi_dif, theta, positionsArr)
    # observer axes
    axes = observer_axes(n_jer_equ)

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = n_observer_angles(coords[0], coords[1], coords[2])
    print('psi, phi at noon:', psi, phi)

    # for hour in range(24):
    #     MJD += 1 / 24
    #     n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD)
    #     axes = observer_axes(n_jer_equ)
    #
    #     # sun-Jerusalem vector in observer coordinates
    #     coords = n_observer_coords(n_sun_jer_equ, axes)
    #     psi, phi = n_observer_angles(coords[0], coords[1], coords[2])
    #     print(hour, psi, phi)
    #     print(psi_vs_time(MJD, time_zone, DST_START, DST_END, 1))


def psi_vs_phi(year, hour, time_zone, time_difference, positionsArr, start_day=1, end_day=365):

    phi_dif, theta = phiDifAndTheta(time_zone)

    # start at given time on the first day of the given year
    MJD = mjd(1, 1, year, (hour-time_difference, 0))

    # lists to hold phi/psi coordinate values
    x = []
    y = []

    equinoxX = []
    equinoxY = []

    # loop through each day in the given year, getting phi and psi values at 12pm on each day
    for day in range(start_day, end_day+1):

        MJD += 1

        # zenith and sun vectors at 12pm on this day (in both ecliptic and equatorial coords)
        n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD, phi_dif, theta, positionsArr)

        # set up observer axes
        axes = observer_axes(n_jer_equ)

        # sun-Jerusalem vector in observer coordinates
        coords = n_observer_coords(n_sun_jer_equ, axes)
        psi, phi = n_observer_angles(coords[0], coords[1], coords[2])
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


# def equinoxDate(year):

    # spring equinox is either march 19, 20, or 21

    # use polynomial interpolation to determine when sun is exactly 90 degrees above equator
    # psi_vs_phi(year, 11.66, 'GMT', 0)

    # vernal equinox: right between solstices??
    # summer solstice: when psi goes from increasing to decreasing
    # winter solstice: when psi goes from decreasing to increasing

    # method #1: use polynomial interpolation to find solstice dates; equinoxes are directly beteween them
    # method #2: equinox = time at which psi at the equator is 0 and phi is 90


# set up Right and Up axes given the angles of the camera (relative to the observer)
def camera_axes(MJD, p_jer_ecl, positionsArr, obs_axes):

    e_N_equ, e_E_equ, e_Z_equ = obs_axes
    e_Z_ecl = equatorial_to_ecliptic(e_Z_equ)

    p_moon = get_pos(MJD, 'm', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)
    p_jer = p_jer_ecl + p_earth

    # camera vector (ecliptic coords)
    n_cam = normalize(p_moon-p_jer)
    n_cam = np.array([float(n_cam[0]), float(n_cam[1]), float(n_cam[2])])

    # camera's position vector in equatorial coords, given psi and psi angles relative to the observer
    # n_cam = n_from_observer_angles(psi_cam, phi_cam, obs_axes)

    # Right and Up axis vectors (normalized)
    e_R = np.cross(n_cam, e_Z_ecl)
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

    equinoxX = []   # x/y coordinates at equinoxes and solstices
    equinoxY = []

    monthsX = []    # x/y coordinates at the first of each month
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
    MJD = mjd(start_day-1, 1, year, (hour - time_difference, 0))  # subtract one because add 1 in loop

    # for each day in the given year, get sun position at given time and convert to camera coordinates
    # for day in range(start_day, end_day):
    for day in range(start_day, start_day+1):
        day_of_month_cur += 1
        MJD += 1

        # zenith and sun vectors at given time on this day (in both ecliptic and equatorial coords)
        n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD, phi_dif, theta, positionsArr)

        p_moon = (get_pos(MJD, 'm', positionsArr) - get_pos(MJD, 'e', positionsArr)) * au_cm
        p_moon_equ = ecliptic_to_equatorial(p_moon)

        # moon's position relative to Jerusalem
        n_moon_jer_equ = n_body_jer(n_jer_equ, p_moon_equ)

        # set up OBSERVER axes according to current position of Jerusalem
        axes = observer_axes(n_jer_equ)

        # now, set up CAMERA axes using observer axes and given position of camera
        cam_axes = camera_axes(MJD, n_jer_ecl*earth_rad, positionsArr, axes)

        # if body == 's':
        # and determine sun position vector in camera coords
        sun_x, sun_y = n_camera_coords(n_sun_jer_equ, cam_axes)

        # elif body == 'm':
        x, y = n_camera_coords(n_moon_jer_equ, cam_axes)

        draw_moon(MJD, cam_axes, n_jer_equ*earth_rad, positionsArr)
        # RA_Dec_lines(cam_axes, n_jer_equ*earth_rad)

        ### store the determined x/y values in the correct list for later plotting
        # days around and including rosh chodesh
        if len(rosh_chodesh_MJD) > 0 and (int(MJD) == rosh_chodesh_MJD[0][0] or int(MJD)-1 == rosh_chodesh_MJD[0][0]
                                          or int(MJD)+1 == rosh_chodesh_MJD[0][0]):
        # if len(rosh_chodesh_MJD) > 0 and ((int(MJD), 'Cheshvan') in rosh_chodesh_MJD or
        #                                   (int(MJD)-1, 'Cheshvan') in rosh_chodesh_MJD or (int(MJD)+1, 'Cheshvan') in rosh_chodesh_MJD):

            roshChodeshX.append(x)
            roshChodeshY.append(y)

            sunX.append(sun_x)
            sunY.append(sun_y)
            # plt.text(x, y, rosh_chodesh_MJD[0][1])

            if int(MJD) > rosh_chodesh_MJD[0][0]:
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


def plot_sun_moon(MJD, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr):
    # n_cam = np.cross(cam_axes[1], cam_axes[0])
    n_cam = cam_axes[2]
    p_moon = get_pos(MJD, 'm', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)

    p_moon_earth_ecl = p_moon-p_earth

    # moon's position relative to Jerusalem
    n_moon_jer_ecl = n_body_jer(n_jer_ecl, p_moon_earth_ecl)

    # sun and moon vectors in camera coords
    sun_x, sun_y = n_camera_coords(n_sun_jer_ecl, cam_axes)
    moon_x, moon_y = n_camera_coords(n_moon_jer_ecl, cam_axes)

    # plot
    if n_sun_jer_ecl.dot(n_cam) > 0:
        plt.scatter(sun_x, sun_y, color='red', s=400, marker='*')

    if n_moon_jer_ecl.dot(n_cam) > 0:
        plt.scatter(moon_x, moon_y, color='blue', s=100, marker='.')

    return [sun_x, sun_y], [moon_x, moon_y]


def drawHorizon(cam_axes, n_jer_equ):
    # n_cam = np.cross(cam_axes[1], cam_axes[0])
    n_cam = cam_axes[2]
    e_Z = equatorial_to_ecliptic(n_jer_equ)
    e_R = cam_axes[0]
    e_F = np.cross(e_Z, e_R)

    points = []

    for i in range(0, 360, 1):
        psi = i * np.pi/180
        p = (e_R * np.cos(psi)) + (e_F * np.sin(psi))
        x, y = n_camera_coords(p, cam_axes)  # in camera coords

        if p.dot(n_cam) > 0:  # if visible
            points.append([x, y])

    for i in range(len(points)-1):
        plt.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 'b-', markersize=1)


def RA_Dec_lines(cam_axes, p_jer_equ, imageData):

    xMin = -0.999
    xMax = 0.999
    xRange = xMax - xMin
    yMin = -0.999
    yMax = 0.999
    yRange = yMax - yMin
    # xMin = -0.01
    # xMax = 0.01
    # xRange = xMax - xMin
    # yMin = -0.01
    # yMax = 0.01
    # yRange = yMax - yMin

    n_cam = cam_axes[2]

    # set plot axis boundaries
    # plt.xlim(-0.25, 0.25)
    # plt.ylim(-0.25, 0.25)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.gca().set_aspect('equal', adjustable='box')

    # plot horizon
    drawHorizon(cam_axes, normalize(p_jer_equ))

    # RA LINES
    # theta is held constant; phi varies
    for i in range(0, 180, 30):
        theta = i * np.pi/180
        points = []
        for j in range(0, 360):  # starting abt 90 degrees east of phi_jer
            # account for earth's rotation??

            phi = j * np.pi / 180

            # Project this point onto x/y plane (i.e. convert to camera coordinates)
            n_equ = n_spherical_to_cartesian(theta, phi)  # unit position vector (centered at earth's center)
            p_equ = celestial_rad * n_equ                 # scaled up by radius of celestial sphere

            # angleFromCam = acos((n_equ.dot(n_cam))/(np.sqrt(n_equ.dot(n_equ))*np.sqrt(n_cam.dot(n_cam))))
            # print(angleFromCam)
            n_ecl = equatorial_to_ecliptic(n_equ)
            x, y = n_camera_coords(n_ecl, cam_axes)       # in camera coords

            if p_equ.dot(p_jer_equ) > 0 and n_equ.dot(n_cam) > 0:  # if visible
                points.append([x, y])
                # canvasX = round(((x - xMin) / xRange) * canvas_size)
                # canvasY = round(((y - yMin) / yRange) * canvas_size)
                # if canvasX == 256: canvasX = 255
                # if canvasY == 256: canvasY = 255
                # if 0 < canvasX < 256 and 0 < canvasY < 256:
                #     imageData[canvasX, canvasY] = [255, 255, 255]

        for j in range(len(points)-1):
            plt.plot([points[j][0], points[j + 1][0]], [points[j][1], points[j + 1][1]], 'k-')

    # DECLINATION LINES
    # phi is held constant; theta varies
    for i in range(0, 361, 30):
        phi = i * np.pi / 180
        points = []
        for j in range(0, 181):
            theta = j * np.pi / 180

            # Project this point onto x/y plane (i.e. convert to camera coordinates)
            n_equ = n_spherical_to_cartesian(theta, phi)  # unit position vector (centered at earth's center)
            p_equ = celestial_rad * n_equ            # scaled up by radius of celestial sphere

            n_ecl = equatorial_to_ecliptic(n_equ)
            x, y = n_camera_coords(n_ecl, cam_axes)  # in camera coords

            if p_equ.dot(p_jer_equ) > 0 and n_equ.dot(n_cam) > 0:  # if visible
                points.append([x, y])
                # canvasX = round(((x - xMin) / xRange) * canvas_size)
                # canvasY = round(((y - yMin) / yRange) * canvas_size)
                # if canvasX == 256: canvasX = 255
                # if canvasY == 256: canvasY = 255
                # if 0 < canvasX < 256 and 0 < canvasY < 256:
                #     imageData[canvasX, canvasY] = [255, 255, 255]

        for j in range(len(points) - 1):
            plt.plot([points[j][0], points[j + 1][0]], [points[j][1], points[j + 1][1]], 'k-', markersize=1)

    # plt.show()
    return imageData


def moving_RA_Dec(MJD, phi_dif, theta, positionsArr):

    # clear out the folder from previous runs
    dir1 = 'C:/Users/tywer/Documents/Celestial Coordinates/Plots'
    dir2 = 'C:/Users/tywer/Documents/Celestial Coordinates/MoonImages'
    for file in os.scandir(dir1):
        os.remove(file.path)

    for file in os.scandir(dir2):
        os.remove(file.path)

    for interval in range(1):
        # MJD += 1/480  # every 3 min
        MJD += 1
        n_jer_ecl, n_jer_equ, n_sun_jer_ecl, n_sun_jer_equ = zenithAndSunVectors(MJD, phi_dif, theta, positionsArr)

        # set up OBSERVER axes according to current position of Jerusalem
        obs_axes = observer_axes(n_jer_equ)

        # now, set up CAMERA axes using observer axes and given position of camera
        cam_axes = camera_axes(MJD, n_jer_ecl*earth_rad, positionsArr, obs_axes)

        # finally, draw RA and Dec lines and sun/moon
        imageData = np.zeros((canvas_size, canvas_size, 3), dtype=np.uint8)

        f = plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')

        sun_coords, moon_coords = plot_sun_moon(MJD, cam_axes, n_jer_ecl, n_sun_jer_ecl, positionsArr)
        imageData = RA_Dec_lines(cam_axes, n_jer_equ * earth_rad, imageData)
        # imageData = draw_sun(MJD, cam_axes, n_jer_ecl*earth_rad, imageData, positionsArr)
        # draw_moon(MJD, cam_axes, n_jer_ecl*earth_rad, moon_coords, imageData, positionsArr)

        draw_moon2(MJD, cam_axes, n_jer_ecl*earth_rad, positionsArr, interval)

        plt.title("MJD:" + str(MJD))
        plt.xlabel("Right", fontfamily="times new roman")
        plt.ylabel("Up", fontfamily="times new roman")

        plt.show()
        f.savefig("./Plots/"+str(interval))

    # movie("MoonImages", "MoonImage")


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

    out = cv2.VideoWriter(movie_name+'.avi', cv2.VideoWriter_fourcc(*'DIVX'), 8, size)

    for image in sorted_images:
        out.write(cv2.imread(image))

    cv2.destroyAllWindows()
    out.release()


def draw_moon(MJD, cam_axes, p_jer_ecl, moon_coords, imageData, positionsArr):

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

    p_sun = get_pos(MJD, 's', positionsArr)
    p_moon = get_pos(MJD, 'm', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)
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

    em_z = np.array([0, -np.sin(eps), np.cos(eps)])
    em_x = np.array([1, 0, 0])
    em_y = np.array([0, np.cos(eps), np.sin(eps)])

    mjd_epoch = mjd(1, 1, 2000, (0, 0))
    days_since_epoch = mjd_cur - mjd_epoch

    rotation_phase = rotation_phase0 + ((2 * np.pi * days_since_epoch) / T_side)

    # given rotation phase, rotate surface of moon around em_z
    es_x = rot_moon_axis(em_x, rotation_phase)
    es_y = rot_moon_axis(em_y, rotation_phase)
    es_z = rot_moon_axis(em_z, rotation_phase)

    return np.array([es_x, es_y, es_z])


def draw_moon2(MJD, cam_axes, p_jer_ecl, positionsArr, image_num):

    N = 200
    D_fov = 0.015  # field of view, in radians
    canvas_size = (2*N) + 1
    imageData = np.zeros((canvas_size, canvas_size, 3), dtype=np.uint8)

    moonColorIm = Image.open('lroc_color_poles_2k.tif')
    im_array = np.asarray(moonColorIm)

    # moon axes rotated around the current phase angle
    es_x, es_y, es_z = rotated_moon_axes(MJD)

    # e_R, e_U = equatorial_to_ecliptic(cam_axes[0]), equatorial_to_ecliptic(cam_axes[1])
    e_R, e_U = cam_axes[0], cam_axes[1]     # in equatorial coords??
    p_moon = get_pos(MJD, 'm', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)
    p_sun = get_pos(MJD, 's', positionsArr)

    p_jer = p_earth + p_jer_ecl
    x_c, y_c, z_c = float(p_moon[0]), float(p_moon[1]), float(p_moon[2])
    x_p, y_p, z_p = float(p_jer[0]), float(p_jer[1]), float(p_jer[2])
    n_cam = np.cross(e_U, e_R)

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


                p_s1 = p_jer + (roots[0] * n_pix)  # points on surface of moon with which n_pix intersects
                p_s2 = p_jer + (roots[1] * n_pix)

                # vectors from center of the moon to the points on its surface
                s1 = p_s1 - p_moon
                s2 = p_s2 - p_moon

                # rotate s1 and s2 according to the moon's rotation
                s1_rotated = np.array([np.dot(s1, es_x), np.dot(s1, es_y), np.dot(s1, es_z)])
                s2_rotated = np.array([np.dot(s2, es_x), np.dot(s2, es_y), np.dot(s2, es_z)])

                # convert to phi and theta
                r1, phi1, theta1 = cartesian_to_spherical(s1_rotated)
                r2, phi2, theta2 = cartesian_to_spherical(s2_rotated)

                # plot the point that's visible to the observer
                if s1.dot(p_s1 - p_jer) < 0:
                    # map this (phi, theta) point onto the 2D moon-color image
                    x1 = (2048 * ((phi1 + 2*np.pi) / (2*np.pi))) % 2048
                    y1 = 1024 * (theta1 / np.pi)
                    color = im_array[int(y1)][int(x1)]

                    # illuminated?
                    if s1.dot(p_s1 - p_sun) >= 0:
                        imageData[N-j, i+N] = color / 4
                    else:
                        imageData[N-j, i+N] = color

                    # if j == 0:
                    #     print("s1", i, theta1, phi1, x1, y1)

                elif s2.dot(p_s2 - p_jer) < 0:
                    x2 = (2048 * ((phi2 + 2*np.pi) / (2 * np.pi))) % 2048
                    y2 = 1024 * (theta2 / np.pi)

                    # if j == 0:
                    #     print("s2", i, theta2, phi2, x2, y2)

                    color = im_array[int(y2)][int(x2)]

                    # illuminated?
                    if s2.dot(p_s2 - p_sun) >= 0:
                        imageData[N-j, i+N] = color / 4
                    else:
                        imageData[N-j, i+N] = color

    image = Image.fromarray(imageData)
    file_name = "./MoonImages/" + str(image_num) + ".png"
    image.save(file_name)
    image.show()


def draw_sun(MJD, cam_axes, p_jer_ecl, imageData, positionsArr):
    n_cam = np.cross(cam_axes[1], cam_axes[0])

    p_sun = get_pos(MJD, 's', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)
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
            n_sun_jer_equ = ecliptic_to_equatorial(normalize(p_s - p_jer))
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


# returns array consisting of the angles between the sun, earth, and moon
def relative_angles(s, e, m):
    E_M = e - m  # Earth - Moon
    E_S = e - s  # Earth - Sun
    M_S = m - s  # Moon - Sun

    theta1 = np.arccos(float((E_S.dot(M_S)) / (length(E_S) * length(M_S))))  # angle btwn E-S and M-S
    theta2 = np.arccos(float((E_S.dot(E_M)) / (length(E_S) * length(E_M))))  # angle btwn E-S and E-M
    theta3 = np.arccos(float((-E_M.dot(M_S)) / (length(E_M) * length(M_S))))  # angle btwn E-S and E-M
    return np.array([theta1, theta2, theta3])


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

    # numbers of days between each molad
    month_lens = [molad_times[i]-molad_times[i-1] for i in range(1, len(molad_times))]
    return np.mean(month_lens)


def main():
    # day = 13
    # month = 9
    # year = 2023
    # time_zone = 'IST'

    year = int(input("Year: "))
    month = int(input("Month: "))
    day = int(input("Day: "))
    hour = int(input("Hour (0-23): "))
    minutes = float(input("Minutes: "))
    location = input("Location ('ny' or 'jer'): ")

    if location == "ny":    time_zone = 'EST'
    elif location == "jer": time_zone = 'IST'
    else: print("invalid location")

    # daylight savings time beginning and end for this year in this time zone
    DST_START, DST_END = DST(time_zone, year)

    # mjd of date in question (in GMT)
    MJD = mjd(day, month, year, (hour, minutes))  # modified julian date
    MJD -= time_dif(MJD, time_zone, DST_START, DST_END)/24

    # array containing positions of sun, earth, and moon as returned from SunEarthMoonPosition.py
    positionsArr = get_positions(year)

    # print(sunriseSunset(MJD, time_zone, DST_START, DST_END, positionsArr))

    phi_dif, theta = phiDifAndTheta(time_zone)
    moving_RA_Dec(MJD, phi_dif, theta, positionsArr)
    movie("Plots", "Ra_Dec")


main()


# test_positions(2023)

