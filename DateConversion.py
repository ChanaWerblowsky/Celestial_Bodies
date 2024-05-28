
# a global set keeps track of which years of lunar cycle are leap years
LEAP_YEARS = {3, 6, 8, 11, 14, 17, 19}

# global array storing days of the week
DAYS = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]

MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
          "October", "November", "December"]


# Returns length of given gregorian month
def len_greg_month(month, leap):

    # index of list is month number-1
    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if month == 2 and leap:
        return 29

    return month_lengths[month-1]


# Returns length of given hebrew month
def len_hebrew_month(month, leap, year_type):

    # if cheshvan and year is meleiah, month has 30 days instead of 29
    if month == 2 and year_type == 'c':
        return 30

    # if kislev and year is chaserah, month has 29 days instead of 30
    if month == 3 and year_type == 'd':
        return 29

    # if leap year, Adar I has 30 days, Adar II has 29
    if leap:
        if month == 6:
            return 30

        if month == 7:
            return 29

    # otherwise, even-number months have 29 days, odd months have 30 (Tishrei = month 1)
    if month % 2 == 0:
        return 29

    return 30


# Returns number of days in hebrew year
def length_hebrew_yr(leap, year_type):
    if not leap:
        if year_type == 'c':    # meleiah (complete)
            return 355

        if year_type == 'r':    # k'sidrah (regular)
            return 354

        if year_type == 'd':    # chaserah (deficient)
            return 353

    # if leap year:
    if year_type == 'c':    # meleiah (complete)
        return 385

    if year_type == 'r':    # k'sidrah (regular)
        return 384

    if year_type == 'd':    # chaserah (deficient)
        return 383


# Multiplies time tuple of form (days, hours, chalakim) by an integer; returns new time tuple (dd always <= 6)
def time_multiplication(time_tuple, factor, acc_days=False):

    # unpack time tuple
    days, hours, chalakim = time_tuple

    # total number of chalakim when time tuple is multiplied by given factor
    chalakim = factor * ((days * 24 * 1080) + (hours * 1080) + chalakim)

    # carry chalakim in excess of 1079
    hours = chalakim // 1080
    chalakim %= 1080

    # carry hours in excess of 23
    days = hours // 24
    hours %= 24

    # disregard days of complete weeks
    if not acc_days: days %= 7

    return days, hours, chalakim


# Adds given time tuple to original molad; returns second molad (days always <= 6)
def time_addition(tuple1, tuple2, acc_days=False):

    # unpack time tuples
    d1, h1, p1 = tuple1
    d2, h2, p2 = tuple2

    chalakim = p1 + p2
    hours = h1 + h2
    days = d1 + d2

    hours += chalakim // 1080
    chalakim %= 1080

    days += hours // 24
    hours %= 24

    if not acc_days: days %= 7

    return days, hours, chalakim


# def tuple_to_hours(time_tuple):
#
#     return (time_tuple[0] * 24) + time_tuple[1] + (time_tuple[2] / 1080)


# Given Hebrew month and year (both ints), returns the month's molad (in form (days, hours, chalakim))
# (Assuming that if leap year, Adar II will count as month #7, etc.)
def molad_determination(month, year, acc_days=False):

    molad_tohu = (1, 5, 204)   # 2d 5h 204p
    month_len = (29, 12, 793)  # 29d 12h 793p

    # First, determine number of months that have passed since Creation

    # number of complete 19-year cycles
    lunar_cycles = (year-1) // 19

    # total number of months in those complete cycles
    months = (lunar_cycles * ((12 * 12) + (7 * 13)))

    # complete years remaining in excess of complete cycles
    remaining_years = (year-1) % 19

    # months from the remaining complete years (# of months per year depends on type of year)
    for i in range(remaining_years):  # for each remaining year
        if i + 1 in LEAP_YEARS:       # if it's a leap year,
            months += 13              # add 13 months to month-count

        else:
            months += 12             # add 12 months otherwise

    # months from current (possibly incomplete) year (up until month in question)
    months += month-1

    # Once know how many months have passed, multiply that # by the length of a month
    # to determine how much the molad has advanced since molad tohu
    molad_advancement = time_multiplication(month_len, months, acc_days)

    # Finally, add the molad's advancement to molad tohu to determine current molad
    molad = time_addition(molad_tohu, molad_advancement, acc_days)

    return molad


# Returns True if time tuple1 is later in the week than (or the same as) time tuple2, and False otherwise
def is_later(tuple1, tuple2):

    # unpack time tuples
    d1, h1, p1 = tuple1
    d2, h2, p2 = tuple2

    molad1 = d1 + ((h1 + (p1/1080))/24)
    molad2 = d2 + ((h2 + (p2/1080))/24)

    return molad1 >= molad2


# Returns tuple indicating day of Rosh Hashanah, whether year is complete/regular/deficient, and first day of Pesach
def hebrew_year_type(molad_tishrei, year_of_lunar_cycle, leap):

    molad_intervals = [(6, 18, 0), (0, 9, 204), (0, 20, 491), (1, 15, 589), (1, 18, 0), (2, 9, 204), (2, 18, 0),
                       (3, 11, 695), (4, 9, 204), (4, 18, 0), (5, 0, 408), (5, 9, 204), (5, 20, 491)]

    # Molad category #1
    if is_later(molad_tishrei, molad_intervals[0]) or not is_later(molad_tishrei, molad_intervals[1]):

        # Both non-leap years and leap years are deficient
        if leap: return 2, 'd', 5
        return 2, 'd', 3

    # Molad category #2
    if is_later(molad_tishrei, molad_intervals[1]) and not is_later(molad_tishrei, molad_intervals[2]):

        # Non-leap years are abundant, leap years are deficient
        if leap: return 2, 'd', 5
        return 2, 'c', 5

    # Molad category #3
    if is_later(molad_tishrei, molad_intervals[2]) and not is_later(molad_tishrei, molad_intervals[3]):

        # Both leap years and non-leap years are abundant
        if leap: return 2, 'c', 7
        return 2, 'c', 5

    # Molad category #4
    if is_later(molad_tishrei, molad_intervals[3]) and not is_later(molad_tishrei, molad_intervals[4]):

        if leap: return 2, 'c', 7
        if year_of_lunar_cycle in [2, 5, 10, 13, 16]: return 2, 'c', 5
        return 3, 'r', 5

    # Molad category #5
    if is_later(molad_tishrei, molad_intervals[4]) and not is_later(molad_tishrei, molad_intervals[5]):

        if leap: return 3, 'r', 7
        return 3, 'r', 5

    # Molad category #6
    if is_later(molad_tishrei, molad_intervals[5]) and not is_later(molad_tishrei, molad_intervals[6]):

        if leap: return 3, 'r', 7
        return 5, 'r', 7

    # Molad category #7
    if is_later(molad_tishrei, molad_intervals[6]) and not is_later(molad_tishrei, molad_intervals[7]):

        if leap: return 5, 'd', 1
        return 5, 'r', 7

    # Molad category #8
    if is_later(molad_tishrei, molad_intervals[7]) and not is_later(molad_tishrei, molad_intervals[8]):

        if leap: return 5, 'c', 3
        return 5, 'r', 7

    # Molad category #9
    if is_later(molad_tishrei, molad_intervals[8]) and not is_later(molad_tishrei, molad_intervals[9]):

        if leap: return 5, 'c', 3
        return 5, 'c', 1

    # Molad category #10
    if is_later(molad_tishrei, molad_intervals[9]) and not is_later(molad_tishrei, molad_intervals[10]):

        if leap: return 7, 'd', 3
        return 7, 'd', 1

    # Molad category #11
    if is_later(molad_tishrei, molad_intervals[10]) and not is_later(molad_tishrei, molad_intervals[11]):

        if leap: return 7, 'd', 3
        if year_of_lunar_cycle in [1, 4, 9, 12, 15]: return 7, 'c', 3
        return 7, 'd', 1

    # Molad category #12
    if is_later(molad_tishrei, molad_intervals[11]) and not is_later(molad_tishrei, molad_intervals[12]):

        if leap: return 7, 'd', 3
        return 7, 'c', 3

    # Molad category #13
    if is_later(molad_tishrei, molad_intervals[12]) and not is_later(molad_tishrei, molad_intervals[0]):

        if leap: return 7, 'c', 5
        return 7, 'c', 3


# Returns dictionary containing each year type (tuple) as keys and their fraction of total time as data
def yr_type_incidence():

    d = {}

    for i in range(19000):

        # which year of 19-year cycle is this?
        year_of_cycle = i % 19
        if year_of_cycle == 0:
            year_of_cycle = 19

        # is it a leap year?
        leap = True if year_of_cycle in {3, 6, 8, 11, 14, 17, 19} else False

        # get year type (tuple)
        year_type = hebrew_year_type(molad_determination(1, i), year_of_cycle, leap)

        # if already encountered this year type, increment its tally
        if year_type in d:
            d[year_type] += 1

        # otherwise, add this year type as a key in the dictionary with data of 1
        else:
            d[year_type] = 1

    # Now should have dictionary with 14 keys (each a different keviyah), where data is tally of years with that keviyah

    # For each year type, convert tally of years into fraction of total years
    for k in d:
        d[k] /= 19000

        # d[k] = (d[k] / 10000) * 100

    return d


def year_lists():

    day_num_hebrew_ny = 1
    day_num_greg_ny = 117

    weekday_greg_ny = 5
    weekday_hebrew_ny = 1

    gregorian_years = [(-3761, -248, 4, True)]
    hebrew_years = []

    # Loop 10,000 times, building up lists of hebrew and gregorian years
    for i in range(1, 6001):

        hebrew_leap = False
        greg_leap = False

        # Which year is this?

        # Gregorian (solar) count starts at 3761 BCE but skips year 0
        if i < 3761:
            gregorian_year = i - 3761

        else:
            gregorian_year = i - 3760

        # Hebrew count starts at 1
        hebrew_year = i

        # Once know year number, determine whether it's a hebrew/gregorian leap year:

        # Given hebrew year, determine which year it is of the current 19-year cycle
        year_of_lunar_cycle = hebrew_year % 19
        if year_of_lunar_cycle == 0:
            year_of_lunar_cycle = 19

        # Is it a hebrew leap year? (depends on year of cycle)
        if year_of_lunar_cycle in LEAP_YEARS:
            hebrew_leap = True

        # Is it a gregorian leap year?
        # (leap year every 4 years - except years that are divisible by 100 and aren't divisible by 400)
        if gregorian_year % 4 == 0 and \
                (gregorian_year % 100 != 0 or (gregorian_year % 100 == 0 and gregorian_year % 400 == 0)):
            greg_leap = True

        # Determine this year's molad tishrei
        molad_tishrei = molad_determination(1, hebrew_year)

        # Then, use molad and year-of-cycle to determine whether year is deficient, regular, or complete
        year_type = hebrew_year_type(molad_tishrei, year_of_lunar_cycle, hebrew_leap)[1]

        # Then, append to list of gregorian/hebrew years:
        # the gregorian/hebrew year number, the day-number and day-of-the-week of its first day, a boolean indicating
        # whether it's a leap year, and the year_type (for hebrew years)

        gregorian_years.append((gregorian_year, day_num_greg_ny, weekday_greg_ny, greg_leap))

        hebrew_years.append((hebrew_year, day_num_hebrew_ny, weekday_hebrew_ny, hebrew_leap, year_of_lunar_cycle, year_type))

        ############

        # Then, figure out how many days THIS year will have to determine the starting date of the NEXT year

        # HEBREW YEAR:

        # Determine how many days this hebrew year has according to year-type and whether it's a leap year
        len_hebrew_yr = length_hebrew_yr(hebrew_leap, year_type)

        # Finally, determine day number and day of the week of next year's Alef Tishrei
        day_num_hebrew_ny += len_hebrew_yr

        weekday_hebrew_ny = (weekday_hebrew_ny + len_hebrew_yr % 7) % 7

        ####

        # GREGORIAN YEAR:

        # Number of days in gregorian year
        len_greg_yr = 366 if greg_leap else 365

        # Day number and day of week of next gregorian new year (January 1)
        day_num_greg_ny += len_greg_yr

        weekday_greg_ny = (weekday_greg_ny + len_greg_yr % 7) % 7

    # print(gregorian_years[1858+3760])
    return gregorian_years, hebrew_years


def hebrew_to_greg(day, month, year):

    gregorian_years, hebrew_years = year_lists()

    # get and unpack tuple of (year number, new-year day number, day of week of new-year, leap (bool), year_type)
    year, ny_day_num, weekday_ny, leap, year_of_cycle, year_type = hebrew_years[year - 1]

    # Error checking:
    # Check for invalid month entries (negative / too high / fractional)
    if month < 1 or ((not leap and month > 12) or (leap and month > 13)) or month % (month // 1):
        return 'Error: invalid month entered'

    # get day from user
    # day = int(input("Enter day number: "))

    # Check for invalid day entries (negative / too high / fractional)
    if day < 1 or day > len_hebrew_month(month, leap, year_type) or day % (day // 1):
        return 'Error: invalid day entered'

    # pointer to begin at new-year (alef Tishrei)
    day_cur = 1
    month_cur = 1
    day_num = ny_day_num
    day_of_week = weekday_ny

    month_len = len_hebrew_month(month_cur, leap, year_type)

    # iterate through days of hebrew calendar (of that year) until reach date in question,
    # keeping track of # days since new year
    while (day_cur, month_cur) != (day, month):
        day_cur += 1
        day_num += 1
        day_of_week = (day_of_week + 1) % 7

        if day_cur > month_len:
            day_cur = 1
            month_cur += 1

            month_len = len_hebrew_month(month_cur, leap, year_type)

    # Now know day-number of hebrew date inputted

    #######

    # look in greg year-list to find year that begins with a day-num <= determined day-num
    # (where next new-year's day-num is greater)
    greg_year = gregorian_years[year - 1]

    if day_num < greg_year[1]: greg_year = gregorian_years[year - 2]
    elif day_num > greg_year[1]: greg_year = gregorian_years[year]

    year, ny_day_num, weekday_ny, leap = greg_year

    # pointer to begin at new-year (1 January)
    day_cur = 1
    month_cur = 1

    month_len = len_greg_month(month_cur, leap)

    for i in range(day_num - ny_day_num):
        # move to next date on gregorian calendar

        day_cur += 1

        if day_cur > month_len:
            day_cur = 1
            month_cur += 1

            month_len = len_greg_month(month_cur, leap)

    # Now know the day of the week and the Gregorian date of the date given
    # Convert numbers into user-friendly output:

    # String for day of the week
    day_of_week = DAYS[day_of_week]

    month_cur = MONTHS[month_cur - 1]

    # Finally, return the day of the week and the Gregorian date
    return day_of_week, day_cur, month_cur, year


def greg_to_hebrew(day, month, year):

    gregorian_years, hebrew_years = year_lists()

    # get and unpack tuple of (year number, new-year day number, day of week of new-year, leap (bool))
    greg_year, ny_day_num, weekday_ny, leap = gregorian_years[year + 3760]

    # Error checking:
    # Check for invalid day entries (negative / too high / fractional)
    if day < 1 or day > len_greg_month(month, leap) or day % (day // 1):
        return 'Error! Invalid day entered.'

    # pointer to begin at new-year (1 January)
    day_cur = 1
    month_cur = 1
    day_num = ny_day_num
    day_of_week = weekday_ny

    month_len = len_greg_month(month_cur, leap)

    # iterate through days of Gregorian calendar (of that year) until reach date in question,
    # keeping track of # days since new year
    while (day_cur, month_cur) != (day, month):
        day_cur += 1
        day_num += 1
        day_of_week = (day_of_week + 1) % 7

        if day_cur > month_len:
            day_cur = 1
            month_cur += 1

            month_len = len_greg_month(month_cur, leap)

    # Now know day-number of Gregorian date inputted

    #######

    # look in Hebrew year-list to find year that begins with a day-num <= determined day-num
    # (where next new-year's day-num is greater)
    hebrew_year = hebrew_years[year + 3760]

    if day_num < hebrew_year[1]: hebrew_year = hebrew_years[year + 3759]

    year, ny_day_num, weekday_ny, leap, yr_of_cycle, year_type = hebrew_year

    # pointer to begin at new-year (Alef Tishrei)
    day_cur = 1
    month_cur = 1

    month_len = len_hebrew_month(month_cur, leap, year_type)

    for i in range(day_num - ny_day_num):

        # move to next date on Hebrew calendar
        day_cur += 1

        if day_cur > month_len:
            day_cur = 1
            month_cur += 1

            # update month-length for new month
            month_len = len_hebrew_month(month_cur, leap, year_type)

    # Now know the day of the week and the Hebrew date of the date given
    # Convert numbers into user-friendly output:

    # Day of the week as a string
    day_of_week = DAYS[day_of_week]

    # Month as a string
    if leap and month_cur == 7:
        month_cur = "Adar II"

    else:
        months = ["Tishrei", "Cheshvan", "Kislev", "Teves", "Shevat", "Adar", "Nisan", "Iyar", "Sivan",
                  "Tammuz", "Av", "Elul"]

        month_cur = months[month_cur - 1]

    # Finally, return the day of the week and the Hebrew date
    return day_of_week, day_cur, month_cur, year


# Given tuple representing either a hebrew or gregorian dd/mm/yy, returns tuple representing date in other calendar
def convert():

    # get calendar type from user
    date_type = (input("Is this a Hebrew or Gregorian date? (Enter 'h' or 'g') ")).lower()

    if date_type != 'g' and date_type != 'h':
        return "Error! Invalid input."

    # if hebrew date, find its year in list of hebrew years --> determine its day-number by adding to new-year's day-num
    if date_type == 'h':

        # get year from user
        year = int(input("Enter year number: "))

        # Check for invalid year entries (negative / too high / fractional)
        if year < 0 or year > 10000 or year % (year // 1):
            return 'Error: invalid year entered'

        # get month and day from user
        month = int(input("Enter month number: "))
        day = int(input("Enter day number: "))

        day_of_week, day_cur, month_cur, year_cur = hebrew_to_greg(day, month, year)

    ####################################

    # if Gregorian date:
    elif date_type == 'g':

        # get year from user
        year = int(input("Enter year number: "))

        # Error checking:
        # Check for invalid year entries (too early / too late / fractional)
        if year < -3761 or year > 6239 or year % (year // 1):
            return "Error! Invalid year entered."

        # get month from user
        month = int(input("Enter month number: "))

        # Check for invalid month entries (negative / too high / fractional)
        if month < 1 or month > 12 or month % (month // 1):
            return 'Error! Invalid month entered.'

        # get day from user
        day = int(input("Enter day number: "))

        # Check that date does not precede world's creation
        if year == -3761 and (month < 9 or (month == 9 and day < 7)):
            return "Error! The date entered precedes the world's creation."

        day_of_week, day_cur, month_cur, year_cur = greg_to_hebrew(day, month, year)

    else:
        return "Invalid date-type entered"

    return day_of_week + ", " + str(day_cur) + " " + month_cur + ", " + str(year_cur)
    ###################################################################################


def main():
    print(convert())


if __name__ == "__main__":
    main()
