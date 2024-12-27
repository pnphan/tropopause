import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geopy.distance import geodesic
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature

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
            

def plot_all_stations(station_data='igra2-station-list.txt'):
    # Read station data and extract latitude and longitude
    with open(station_data, "r") as f:
        stations = f.read().splitlines()[:-101]  # Remove trailing meta-info lines
    
    lat = []
    lon = []
    for station in stations:
        try:
            # Extract latitude and longitude
            latitude = float(station[12:20].strip())
            longitude = float(station[21:30].strip())
            
            # Ignore stations with invalid lat/lon (set to missing values in the file)
            if latitude == -98.8888 or longitude == -998.8888:
                continue
            
            lat.append(latitude)
            lon.append(longitude)
        except ValueError:
            print(f"Skipping invalid line: {station}")
    
    # Plotting
    plt.figure(figsize=(12, 6)) 
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.stock_img()  # Adds a basic map background
    ax.coastlines()  # Adds coastlines for better context
    #ax.add_feature(cfeature.BORDERS, linestyle=':')  # Adds country borders

    ax.add_feature(cfeature.LAND, color='white')
    ax.add_feature(cfeature.OCEAN, color='white')
    ax.add_feature(cfeature.COASTLINE, edgecolor='black')
    #ax.add_feature(cfeature.BORDERS, edgecolor='gray', linestyle=':')
    #ax.add_feature(cfeature.LAKES, color='white', edgecolor='gray')
    #ax.add_feature(cfeature.RIVERS, color='gray')
    
    # Add gridlines for lat/lon
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0, linestyle='--')
    gl.top_labels = False  # Disable labels on the top axis
    gl.right_labels = False  # Disable labels on the right axis
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 60))  # Longitude ticks every 10 degrees
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 30))    # Latitude ticks every 10 degrees
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'}
    
    # Scatter plot of station locations
    plt.scatter(lon, lat, color="red", transform=ccrs.PlateCarree(), marker='.', s=2, label='Stations')
    
    # Title and labels
    plt.title("Station Locations", fontsize=14)
    plt.legend(loc="lower left")
    plt.show()


def plot_stations(station_list=['AEM00041217', 'ACM00078861'], station_data='igra2-station-list.txt'):
    with open(station_data, "r") as f:
        stations = f.read().splitlines()[:-101]  # Remove trailing meta-info lines
    
    lat = []
    lon = []

    target_stations = []
    for i in station_list:
        for j in stations:
            if i == j[:11]:
                target_stations.append(j)

    print(len(target_stations))
    for station in target_stations:
        try:
            # Extract latitude and longitude
            latitude = float(station[12:20].strip())
            longitude = float(station[21:30].strip())
            
            # Ignore stations with invalid lat/lon (set to missing values in the file)
            if latitude == -98.8888 or longitude == -998.8888:
                print(station[:11])
                continue
            
            lat.append(latitude)
            lon.append(longitude)
        except ValueError:
            print(f"Skipping invalid line: {station}")
    
    # Plotting
    plt.figure(figsize=(12, 6)) 
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.stock_img()  # Adds a basic map background
    ax.coastlines()  # Adds coastlines for better context
    #ax.add_feature(cfeature.BORDERS, linestyle=':')  # Adds country borders

    ax.add_feature(cfeature.LAND, color='white')
    ax.add_feature(cfeature.OCEAN, color='white')
    ax.add_feature(cfeature.COASTLINE, edgecolor='black')
    #ax.add_feature(cfeature.BORDERS, edgecolor='gray', linestyle=':')
    #ax.add_feature(cfeature.LAKES, color='white', edgecolor='gray')
    #ax.add_feature(cfeature.RIVERS, color='gray')
    
    # Add gridlines for lat/lon
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0, linestyle='--')
    gl.top_labels = False  # Disable labels on the top axis
    gl.right_labels = False  # Disable labels on the right axis
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 60))  # Longitude ticks every 10 degrees
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 30))    # Latitude ticks every 10 degrees
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'}
    
    # Scatter plot of station locations
    lon.append(-180)
    lon.append(180)
    lat.append(90)
    lat.append(-90)
    plt.scatter(lon, lat, color="red", transform=ccrs.PlateCarree(), marker='o', s=20, label='Stations')
    
    # Title and labels
    plt.title("Station Locations", fontsize=14)
    plt.legend(loc="lower left")
    plt.show()

import numpy as np

def calculate_tropopause_pressure_reichler(temperature, pressure, kappa=0.286, g=9.8, R=287):
    """
    Calculates the tropopause pressure using the Reichler et al. algorithm.
    
    Parameters:
        temperature (np.ndarray): Temperature profile (K) at pressure levels.
        pressure (np.ndarray): Pressure profile (Pa), from surface to upper levels.
        kappa (float): Poisson constant (default is 0.286 for dry air).
        g (float): Acceleration due to gravity (m/s^2), default is 9.8.
        R (float): Gas constant for dry air (J/(kg·K)), default is 287.
    
    Returns:
        p_tropopause (float): Tropopause pressure (Pa).
    """
    # Ensure pressure is descending (surface -> top of atmosphere)
    if not np.all(np.diff(pressure) < 0):
        raise ValueError("Pressure levels must be in descending order.")

    # Calculate half-level pressures using Eq. (3)
    p_half = 0.5 * (pressure[:-1]**kappa + pressure[1:]**kappa)

    # Calculate lapse rate at half levels using Eq. (4)
    lapse_rate_half = (
        (temperature[1:] - temperature[:-1]) *
        (p_half / (pressure[1:]**kappa - pressure[:-1]**kappa)) *
        (kappa * g / R) /
        (temperature[:-1] + temperature[1:])
    )

    # Search for the first level where the lapse rate drops below the critical threshold
    critical_lapse_rate = 2.0  # K/km
    idx = np.where(lapse_rate_half < critical_lapse_rate)[0]

    if len(idx) == 0:
        # No tropopause found
        return np.nan

    # Verify mean lapse rate over the 2 km layer above meets the criteria
    j = idx[0]
    height_diff = np.cumsum(np.diff(pressure) * -R * (temperature[:-1] + temperature[1:]) / (2 * g * pressure[:-1]))
    height_2km_idx = np.where(height_diff >= 2000)[0]
    if j + len(height_2km_idx) >= len(lapse_rate_half) or np.mean(lapse_rate_half[j:j + len(height_2km_idx)]) >= critical_lapse_rate:
        return np.nan

    # Interpolate to find exact tropopause pressure using Eq. (5)
    lapse_rate_j = lapse_rate_half[j]
    lapse_rate_j1 = lapse_rate_half[j + 1]
    p_tropopause = (
        p_half[j] +
        (p_half[j + 1] - p_half[j]) / (lapse_rate_j1 - lapse_rate_j) * (critical_lapse_rate - lapse_rate_j)
    )

    # Restrict results to reasonable range (550 hPa to 75 hPa)
    if p_tropopause < 75000 or p_tropopause > 550000:
        return np.nan

    return p_tropopause


def reichler_tropopause(station, year, month, day, time):
    filepath = station + '-data.txt'
    years = get_years(filepath)
    if year not in years:
        return 'Invalid year'
    else:
        months = get_months(year, filepath)
    if month not in months:
        return 'Invalid month'
    else:
        days = get_days(year, month, filepath)
    if day not in days:
        return 'Invalid day'
    else:
        times = get_times(year, month, day, filepath)
    if time not in times:
        return 'Invalid time'
    
    kappa = 287.053/1005
    pressure = get_daywise_data(year, month, filepath)[(day,time)]['pressure']
    temperature = get_daywise_data(year, month, filepath)[(day,time)]['temp']
    pressure, temperature = interpolate_to_points(pressure, temperature, start=500, end=40, step=-5)
    # print(pressure)
    # print(temperature)

    pressure_half_levels = ((pressure**kappa + np.roll(pressure**kappa,-1))/2)[:-1]
    print(pressure_half_levels)
    lapse_half_levels = ((temperature - np.roll(temperature, 1))[1:] / (pressure**kappa - np.roll(pressure**kappa, 1))[1:]) * ((pressure**kappa + np.roll(pressure**kappa, -1))[:-1] / (temperature + np.roll(temperature, -1))[:-1]) * (kappa * 9.80665/287.053)
    #(pressure - np.roll(pressure, 1))[1:]
    print(lapse_half_levels)
    for i in range(len(lapse_half_levels)):
        if lapse_half_levels[i] < 2:
            if i < len(lapse_half_levels)-9:
                mean_lapse = np.mean(lapse_half_levels[i:i+9])
            else:
                mean_lapse = np.mean(lapse_half_levels[i:])
            if mean_lapse < 2:
                j = i #500 - 5*i
                break

    tropopause = pressure_half_levels[j-1] + (pressure_half_levels[j] - pressure_half_levels[j-1])/(lapse_half_levels[j] - lapse_half_levels[j-1]) * (2 - lapse_half_levels[j-1])
    tropopause = tropopause**(1/kappa)

    return tropopause




#pressure = get_daywise_data(2024, 7)[(31,0)]['pressure']
#print(pressure)
# gph = get_daywise_data(2024, 7)[(31,0)]['gph']
# temp = get_daywise_data(2024, 7)[(31,0)]['temp']
# gph, temp = interpolate_to_points(gph, temp, kind='linear', step=200)
# gph, lapse, temp = lapse_rate(gph, temp)
# tropopause = detect_tropopause(gph, lapse)


seidel2006stationlist = [
    'ARM00087576', 'ASM00094120', 'ASM00094294', 'ASM00094610', 'ASM00094672', 
    'ASM00094998', 'AYM00089009', 'AYM00089050', 'AYM00089532', 'AYM00089542', 
    'AYM00089564', 'AYM00089664', 'BDM00078016', 'BPM00091517', 'BRM00082332', 
    'BRM00083746', 'CAM00071072', 'CAM00071082', 'CAM00071801', 'CAM00071836', 
    'CAM00071926', 'CHM00051709', 'CHM00052681', 'CIM00085442', 'CIM00085469', 
    'CIM00085799', 'COM00080222', 'FIM00002836', 'FJM00091680', 'FMM00091334', 
    'FPM00091938', 'FSM00061996', 'GIM00008495', 'GLM00004360', 'GMM00010868', 
    'HKM00045004', 'ICM00004018', 'IOM00061967', 'ISM00040179', 'IVM00065578', 
    'JAM00047401', 'JAM00047827', 'JAM00047991', 'JNM00001001', 'KEM00063741', 
    'LYM00062010', 'MAM00067083', 'NFM00094996', 'NGM00061052', 'NZM00093844', 
    'NZM00093986', 'POM00008508', 'PSM00091408', 'RMM00091376', 'RQM00078526', 
    'RSM00021504', 'RSM00021965', 'RSM00023415', 'RSM00023472', 'RSM00024266', 
    'RSM00028698', 'RSM00030230', 'RSM00032540', 'RSM00034731', 'RSM00035121', 
    'SAM00041024', 'SFM00068588', 'SFM00068816', 'SFM00068994', 'SGM00061641', 
    'SHM00061902', 'SHM00068906', 'SNM00048698', 'SPM00060020', 'THM00048455', 
    'TXM00038880', 'UKM00003005', 'USM00070026', 'USM00070273', 'USM00070308', 
    'USM00070398', 'USM00072201', 'USM00072208', 'USM00072250', 'USM00072251', 
    'USM00072293', 'USM00072327', 'USM00072403', 'USM00072451', 'USM00072493', 
    'USM00072520', 'USM00072645', 'USM00072694', 'USM00072712', 'USM00072747', 
    'USM00072768', 'USM00072776', 'USM00072797', 'USM00091165', 'USM00091285'
]
#plot_stations(station_list=seidel2006stationlist)