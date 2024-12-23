import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

## DATA CUTOFF DATE 2024-12-22 19:26	

def get_station_list(filepath='igra2-station-list.txt'):
    stations = []
    with open(filepath, 'r') as file:
        for line in file:
            stations.append(line[:11])
    return stations

def get_years_and_months(filepath='PFM00059981-data.txt'):
    """
    Extracts available years and their corresponding months from the IGRA file.

    :param filepath: Path to the IGRA sounding data file.
    :return: A dictionary where keys are years and values are lists of available months.
    """
    year_month_data = {}

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    # Extract year and month
                    year = int(line[13:17])  # Year (columns 14-17)
                    month = int(line[18:20])  # Month (columns 19-20)

                    # Add the month to the corresponding year
                    if year not in year_month_data:
                        year_month_data[year] = set()  # Use a set to avoid duplicate months
                    year_month_data[year].add(month)
                except ValueError:
                    continue

    # Convert sets to sorted lists for consistency
    for year in year_month_data:
        year_month_data[year] = sorted(year_month_data[year])

    return year_month_data


def get_available_days_and_times(year, month, filepath='PFM00059981-data.txt'):
    """
    Reads the IGRA file and extracts days and their available times for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A dictionary where keys are days and values are lists of available times.
    """
    available_data = {}

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    line_year = int(line[13:17])  # Extract year (columns 14-17, zero-indexed)
                    line_month = int(line[18:20])  # Extract month (columns 19-20, zero-indexed)
                    line_day = int(line[21:23])  # Extract day (columns 22-23, zero-indexed)
                    line_hour = int(line[24:26])  # Extract hour (columns 25-26, zero-indexed)

                    if line_year == year and line_month == month:
                        if line_day not in available_data:
                            available_data[line_day] = []
                        available_data[line_day].append(line_hour)
                except ValueError:
                    continue

    # Sort the times for each day
    for day in available_data:
        available_data[day] = sorted(available_data[day])

    return available_data


def get_years(filepath='PFM00059981-data.txt'):
    return list(get_years_and_months(filepath).keys())
def get_months(year, filepath='PFM00059981-data.txt'):
    return get_years_and_months(filepath)[year]
def get_days(year, month, filepath='PFM00059981-data.txt'):
    """
    Reads the IGRA file and extracts days where data is available for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A list of days with data available.
    """
    available_days = set()

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    line_year = int(line[13:17])  # Extract year (columns 14-17, zero-indexed)
                    line_month = int(line[18:20])  # Extract month (columns 19-20, zero-indexed)
                    line_day = int(line[21:23])  # Extract day (columns 22-23, zero-indexed)

                    if line_year == year and line_month == month:
                        available_days.add(line_day)
                except ValueError:
                    continue

    return sorted(available_days)
def get_times(year, month, day, filepath='PFM00059981-data.txt'):
    return get_available_days_and_times(year, month, filepath)[day]


def get_daywise_data(year, month, filepath='PFM00059981-data.txt'):
    """
    Reads the IGRA file and extracts geopotential height, pressure, and temperature
    for each available day and time for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A dictionary where each key is a day, and the value is a list of tuples
             (time, geopotential height list, pressure list, temperature list).
    """
    daywise_data = {}

    with open(filepath, 'r') as file:
        current_day = None
        current_time = None

        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    # Extract year, month, day, and time
                    line_year = int(line[13:17])  # Year (columns 14-17)
                    line_month = int(line[18:20])  # Month (columns 19-20)
                    line_day = int(line[21:23])  # Day (columns 22-23)
                    line_hour = int(line[24:26])  # Hour (columns 25-26)
                    if line_year == year and line_month == month:
                        current_day = line_day
                        current_time = line_hour
                        daywise_data[(current_day, current_time)] = {'gph':[], 'pressure':[], 'temp':[]}
                except ValueError:
                    continue

            elif line_year == year and line_month == month:
                gph = int(line[16:21])  # GPH (columns 17-21)
                pressure = int(line[9:15])  # Pressure (columns 10-15)
                temp = int(line[22:27])  # Temperature (columns 23-27)
                if gph != -9999 and pressure != -9999 and temp != -9999:
                    daywise_data[(current_day, current_time)]['gph'].append(gph)
                    daywise_data[(current_day, current_time)]['pressure'].append(pressure / 100) # / 100 convert to hPa
                    daywise_data[(current_day, current_time)]['temp'].append(temp / 10) # / 10 units are in 10ths of a degree celsius
    
    return daywise_data


def interpolate_to_points(x, y, start=5000, end=22000, step=200, kind='linear'):
    """
    Interpolates x to every 200th point between start and end using y values.

    :param x: List or array of x values (independent variable).
    :param y: List or array of y values (dependent variable).
    :param start: Start value for interpolation.
    :param end: End value for interpolation.
    :param step: Step size for interpolation points.
    :return: Interpolated x values and the corresponding points.
    """
    # Ensure x and y are numpy arrays for compatibility
    x = np.array(x)
    y = np.array(y)

    # Define the interpolation function
    interpolation_function = interp1d(
        x, y, kind=kind, bounds_error=False, fill_value="extrapolate"
    )

    # Generate the points for interpolation
    interpolation_points = np.arange(start, end + step, step)

    # Interpolate the x values at the desired points
    interpolated_values = interpolation_function(interpolation_points)

    return interpolation_points, interpolated_values


def plot_graph(x, y, title="X vs Y Graph", xlabel="X-axis", ylabel="Y-axis"):
    """
    Plots a graph of y versus x.

    :param x: List or array of x values.
    :param y: List or array of y values.
    :param title: Title of the graph.
    :param xlabel: Label for the x-axis.
    :param ylabel: Label for the y-axis.
    """
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='b', label="Data", markersize=2.5)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.legend()
    plt.show()


def lapse_rate(gph, temp):
    gph = np.array(gph)
    temp = np.array(temp)
    return gph[:-1], ((-(temp - np.roll(temp, 1)) / (gph - np.roll(gph, 1))) * 1000)[1:], temp[:-1]


def detect_tropopause(gph, lapse):
    for i in range(len(lapse)):
        if lapse[i] < 2:
            if i+10 > len(lapse):
                mean_lapse = np.mean(lapse[i:-1])
            else:
                mean_lapse = np.mean(lapse[i:i+10])
            if mean_lapse > 2:
                continue 
            else:
                return gph[i]





# gph = get_daywise_data(2024, 7)[(31,0)]['gph']
# temp = get_daywise_data(2024, 7)[(31,0)]['temp']
# gph, temp = interpolate_to_points(gph, temp, kind='linear', step=200)
# gph, lapse, temp = lapse_rate(gph, temp)
# tropopause = detect_tropopause(gph, lapse)