

## Calendar Utility Functions

# Returns length of given gregorian month
def month_len(month, leap):
    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if month == 2 and leap:
        return 29

    return month_lengths[month - 1]


# Returns True if the given year is a Gregorian leap year and False otherwise
def isLeap(year):
    return year % 4 == 0 and \
           (year % 100 != 0 or (year % 100 == 0 and year % 400 == 0))


# Returns the Modified Julian Date for the given day and time
def mjd(dd, mm, yy, tt):
    epoch = (17, 11, 1858, (0, 0))  # 17 November, 1858, 12 am (midnight)

    days = 0  # days since epoch
    start = end = 0

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
