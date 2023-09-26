
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
phi_jer = 35.214
phi_ny = -74.006
theta_jer = 58.232
theta_ny = 49.3
dayNum_mjdEpoch = 2052005
sunset_angle = 90.8
obliquity = 23.43
# psi_cam = 90   # camera position for ~8AM
# phi_cam = 90
# psi_cam = 0.1  # camera position for ~12PM
# phi_cam = 180
psi_cam = 85  # camera position for ~6PM
phi_cam = 205
canvas_size = 256

# convert to radians
radians_day = degrees_day * np.pi/180
phi_jer = phi_jer * np.pi/180
phi_ny = phi_ny * np.pi/180
theta_jer = theta_jer * np.pi/180
theta_ny = theta_ny * np.pi/180
sunset_angle = sunset_angle * np.pi/180
obliquity = obliquity * np.pi/180
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

    # interpolate to get set of positions at time in question
    t_Vals = []
    pos_Vals = []
    for j in range(-1, 2):  # use the prev interval, this one, and the following one
        t_Vals.append(positionsArr[est_index+j][0]-tstart)
        pos_Vals.append(positionsArr[est_index+j][i])

    # isolate the variables
    x_Vals = [coords[0] for coords in pos_Vals]
    y_Vals = [coords[1] for coords in pos_Vals]
    z_Vals = [coords[2] for coords in pos_Vals]

    X, Y = symbols('X Y')

    # interpolate x, y, and z coordinates separately
    p_X = interpolation(t_Vals, x_Vals, 2)
    p_Y = interpolation(t_Vals, y_Vals, 2)
    p_Z = interpolation(t_Vals, z_Vals, 2)

    # now, plug the given t into each polynomial; solve for x, y, and z
    x = p_X.subs(X, MJD-tstart)
    y = p_Y.subs(X, MJD-tstart)
    z = p_Z.subs(X, MJD-tstart)

    return np.array([x, y, z])


# def moon_pos(MJD):
#
#     p_moon = get_pos(MJD, 'm')
#     p_earth = get_pos(MJD, 'e')
#
#     return p_moon - p_earth


### Utility functions
def normalize(v):
    length = np.sqrt(float(v.dot(v)))
    return v / length


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
    print('psi, phi at sunrise:', psi, phi)

    ### phi at SUNSET
    n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD_sunset, phi_dif, theta, positionsArr)
    # observer axes
    axes = observer_axes(n_jer_equ)

    # sun-Jerusalem vector in observer coordinates
    coords = n_observer_coords(n_sun_jer_equ, axes)
    psi, phi = n_observer_angles(coords[0], coords[1], coords[2])
    print('psi, phi at sunset:', psi, phi)

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
def camera_axes(psi_cam, phi_cam, observer_axes):

    e_N, e_E, e_Z = observer_axes

    # camera's position vector in equatorial coords, given psi and psi angles relative to the observer
    n_cam = n_from_observer_angles(psi_cam, phi_cam, observer_axes)

    # determine Right and Up axis vectors (normalized)
    e_R = np.cross(n_cam, e_Z)
    e_R_normal = normalize(e_R)

    e_U = np.cross(e_R_normal, n_cam)
    e_U_normal = normalize(e_U)

    return np.array([e_R_normal, e_U_normal])


# convert position vector from equatorial coords to camera coords
def n_camera_coords(n, camera_axis_vectors):

    # set up Right and Up axes
    e_R, e_U = camera_axis_vectors

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
        cam_axes = camera_axes(psi_cam, phi_cam, axes)

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


def plot_sun_moon(MJD, cam_axes, n_jer_equ, n_sun_jer_equ, positionsArr):
    n_cam = np.cross(cam_axes[1], cam_axes[0])

    p_moon = get_pos(MJD, 'm', positionsArr)
    p_earth = get_pos(MJD, 'e', positionsArr)

    p_moon_earth = p_moon-p_earth
    p_moon_earth_equ = ecliptic_to_equatorial(p_moon_earth)

    # moon's position relative to Jerusalem
    n_moon_jer_equ = n_body_jer(n_jer_equ, p_moon_earth_equ)

    # sun and moon vectors in camera coords
    sun_x, sun_y = n_camera_coords(n_sun_jer_equ, cam_axes)
    moon_x, moon_y = n_camera_coords(n_moon_jer_equ, cam_axes)

    # plot
    if n_moon_jer_equ.dot(n_cam) > 0:
        plt.scatter(moon_x, moon_y, color='blue', s=100, marker='*')

    if n_sun_jer_equ.dot(n_cam) > 0:
        plt.scatter(sun_x, sun_y, color='red', s=400, marker='*')

    return [sun_x, sun_y], [moon_x, moon_y]


def drawHorizon(cam_axes, n_jer_equ):
    n_cam = np.cross(cam_axes[1], cam_axes[0])
    e_Z = n_jer_equ
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

    n_cam = np.cross(cam_axes[1], cam_axes[0])

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

            x, y = n_camera_coords(n_equ, cam_axes)       # in camera coords

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
            x, y = n_camera_coords(n_equ, cam_axes)  # in camera coords

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

    return imageData


def moving_RA_Dec(MJD, phi_dif, theta, positionsArr):
    dir = 'C:/Users/tywer/Documents/Celestial Coordinates/Plots'
    for file in os.scandir(dir):
        os.remove(file.path)

    # start at 5pm IST  (use conjunction times??)
    # MJD += 14.66/24
    # MJD += 15.66/24
    MJD += 22.66/24

    for interval in range(1):
        MJD += 1/480  # every 3 min

        n_jer_ecl, n_jer_equ, n_sun_jer, n_sun_jer_equ = zenithAndSunVectors(MJD, phi_dif, theta, positionsArr)

        # set up OBSERVER axes according to current position of Jerusalem
        axes = observer_axes(n_jer_equ)

        # now, set up CAMERA axes using observer axes and given position of camera
        cam_axes = camera_axes(psi_cam, phi_cam, axes)

        # finally, draw RA and Dec lines and sun/moon
        imageData = np.zeros((canvas_size, canvas_size, 3), dtype=np.uint8)

        f = plt.figure()
        plt.gca().set_aspect('equal', adjustable='box')

        sun_coords, moon_coords = plot_sun_moon(MJD, cam_axes, n_jer_equ, n_sun_jer_equ, positionsArr)
        imageData = RA_Dec_lines(cam_axes, n_jer_equ * earth_rad, imageData)
        # imageData = draw_sun(MJD, cam_axes, n_jer_ecl*earth_rad, imageData, positionsArr)
        imageData = draw_moon(MJD, cam_axes, n_jer_ecl*earth_rad, moon_coords, imageData, positionsArr)

        image = Image.fromarray(imageData)
        image.show()

        plt.title("MJD:" + str(MJD))
        plt.xlabel("Right", fontfamily="times new roman")
        plt.ylabel("Up", fontfamily="times new roman")

        plt.show()
        f.savefig("./Plots/"+str(interval))


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [
        int(text)
        if text.isdigit() else text.lower()
        for text in _nsre.split(s)]


def movie():

    img_array = [img for img in glob.glob('C:/Users/tywer/Documents/Celestial Coordinates/Plots/*.png')]
    sorted_images = sorted(img_array, key=natural_sort_key)

    img = cv2.imread(img_array[0])
    height, width, layers = img.shape
    size = (width, height)

    out = cv2.VideoWriter('RaDec.avi', cv2.VideoWriter_fourcc(*'DIVX'), 8, size)

    for image in sorted_images:
        out.write(cv2.imread(image))

    cv2.destroyAllWindows()
    out.release()


def draw_moon(MJD, cam_axes, p_jer_ecl, moon_coords, imageData, positionsArr):

    moon_x, moon_y = moon_coords

    xMin = moon_x-0.01
    xMax = moon_x+0.01
    xRange = xMax - xMin
    yMin = moon_y-0.01
    yMax = moon_y+0.01
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
    for theta in range(0, 180, 5):
        theta = theta * np.pi/180
        for phi in range(0, 360, 5):
            phi = phi * np.pi / 180

            # this point in ecliptic cartesian coords
            s = n_spherical_to_cartesian(theta, phi) * moon_rad  # vector from center of moon to current point (in moon coords)
            p_s = p_moon + s  # ecliptic coords

            n_moon_jer_equ = ecliptic_to_equatorial(normalize(p_s-p_jer))
            x, y = n_camera_coords(n_moon_jer_equ, cam_axes)  # camera (x-y) coords

            if s.dot(p_s-p_jer) < 0 and s.dot(p_s-p_sun) >= 0:
                canvasX = round(((x - xMin) / xRange) * canvas_size)
                canvasY = round(((yMax - y) / yRange) * canvas_size)
                if canvasX == 256: canvasX = 255
                if canvasY == 256: canvasY = 255
                imageData[canvasY, canvasX] = [47, 79, 79]
                imageData[canvasY-1:canvasY+1, canvasX-1:canvasX+1] = [47, 79, 79]


            elif s.dot(p_s-p_jer) < 0 and s.dot(p_s-p_sun) < 0:
                canvasX = round(((x - xMin) / xRange) * canvas_size)
                canvasY = round(((yMax - y) / yRange) * canvas_size)
                if canvasX == 256: canvasX = 255
                if canvasY == 256: canvasY = 255
                imageData[canvasY-1:canvasY+1, canvasX-1:canvasX+1] = [255, 255, 255]

    image = Image.fromarray(imageData)
    # image.show()
    plt.scatter(yellowX, yellowY, color='yellow')
    plt.scatter(blueX, blueY, color='blue')
    return imageData
    # return image


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


def main():

    day = 21
    month = 9
    year = 2023
    time_zone = 'EST'

    # mjd of date in question (in GMT)
    MJD = mjd(day, month, year, (0, 0))  # modified julian date

    # MJD = rosh_chodesh_MJD[5][0] + 1.05

    # array containing positions of sun, earth, and moon as returned from SunEarthMoonPosition.py
    positionsArr = get_positions(year)

    # daylight savings time beginning and end for this year in this time zone
    DST_START, DST_END = DST(time_zone, year)

    # time-difference between local time and GMT
    time_difference = time_dif(MJD, time_zone, DST_START, DST_END)

    # print(sunriseSunset(MJD, time_zone, DST_START, DST_END, positionsArr))

    phi_dif, theta = phiDifAndTheta(time_zone)
    moving_RA_Dec(MJD, phi_dif, theta, positionsArr)
    # movie()


main()

