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
            

def temp_gradient_zangl(gph, temp):
    gph = np.array(gph) / 1000
    temp = np.array(temp)
    tempgrad = (temp[1:] - temp[:-1]) / (gph[1:] - gph[:-1])
    return gph, tempgrad
def detect_tropopause_zangl(tempgrad, gph, temp):
    zhalf = (gph[1:] + gph[:-1])/2
    zTP = []
    possible_heights = []
    tempTP = []
    for i in range(len(tempgrad)):
        if tempgrad[i] >= -2:
            if tempgrad[i-1] < -2:
                ztp = gph[i-1] + (-2 - tempgrad[i])/(tempgrad[i] - tempgrad[i-1])*(zhalf[i] - zhalf[i-1])
                zTP.append(ztp)
                possible_heights.append(i)
                if ztp > gph[i]:
                    tempTP.append(temp[i] + tempgrad[i-1] * (ztp - gph[i]))
                else:
                    tempTP.append(temp[i] + tempgrad[i] * (ztp - gph[i]))
    
    for j in range(len(possible_heights)):
        checknext = False
        if possible_heights[j] < len(gph) - 10:
            for k in range(possible_heights[j], possible_heights[j]+10):
                if (temp[k] - tempTP[j]) / (gph[k] - zTP[j]) <= -2:
                    checknext=True
                    break
        else:
            for k in range(possible_heights[j], len(gph)):
                if (temp[k] - tempTP[j]) / (gph[k] - zTP[j]) <= -2:
                    checknext=True
                    break
        if checknext == False:
            return zTP[j]
        else:
            continue
        

def get_all_stations_flask(filepath='igra2-station-list.txt'):
    # Read station data and extract latitude and longitude
    with open(filepath, "r") as f:
        stations = f.read().splitlines()[:-101]  # Remove trailing meta-info lines
    
    stationdata = []
    for station in stations:
        try:
            # Extract latitude and longitude
            latitude = float(station[12:20].strip())
            longitude = float(station[21:30].strip())
            
            # Ignore stations with invalid lat/lon (set to missing values in the file)
            if latitude == -98.8888 or longitude == -998.8888:
                continue

            stationdata.append({'id': station[:11], 'lat': latitude, 'lon': longitude})
        except ValueError:
            print(f"Skipping invalid line: {station}")

    return stationdata
            

def plot_all_stations(station_data='igra2-station-list.txt'):
    # Read station data and extract latitude and longitude
    with open(station_data, "r") as f:
        stations = f.read().splitlines()[:-101]  # Remove trailing meta-info lines
    
    lat = []
    lon = []
    ids = []
    for station in stations:
        try:
            # Extract latitude and longitude
            ids.append(station[:11])
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
    return ids


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


def seek_stations():
    id_list = plot_all_stations()
    good_stations = []
    for id in id_list:
        available_years = get_years(filepath=f'/media/peter/easystore/igra2-extracted/{id}-data.txt')
        if np.max(available_years) >= 2016 and np.min(available_years) < 2000:
            good_stations.append(id)
            print(id)
    return good_stations


def parse_igra_station(station_id):
    """
    Parse an IGRA v2.2 station file of the form <station_id>-data.txt
    and return the available years, months, days, and times in the format:

    {
        "years": [...],
        "data": {
            "YYYY": {
                "MM": {
                    "DD": ["HHZ", ...],
                    ...
                },
                ...
            },
            ...
        }
    }
    """

    # Prepare the output structure
    data_dict = {
        "years": [],
        "data": {}
    }

    # The IGRA v2.2 header record has fixed-column positions:
    # HEADREC columns 1-1   (index 0)
    # ID      columns 2-12  (index 1..11)
    # YEAR    columns 14-17 (index 13..16)
    # MONTH   columns 19-20 (index 18..19)
    # DAY     columns 22-23 (index 21..22)
    # HOUR    columns 25-26 (index 24..25)
    # We will ignore the rest for collecting date/time info.

    file_name = f"data/{station_id}-data.txt"

    try:
        with open(file_name, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                # Only parse lines that start with '#'
                if not line.startswith("#"):
                    continue

                # Extract the columns defined by IGRA
                # Make sure each slice is large enough and handle any short lines gracefully.
                year_str  = line[13:17].strip()
                month_str = line[18:20].strip()
                day_str   = line[21:23].strip()
                hour_str  = line[24:26].strip()

                # If any of these are missing or not digits, skip.
                # (Sometimes lines can be malformed or short.)
                if (not year_str.isdigit() or 
                    not month_str.isdigit() or 
                    not day_str.isdigit() or 
                    not hour_str.isdigit()):
                    continue

                year  = int(year_str)
                month = int(month_str)
                day   = int(day_str)
                hour  = int(hour_str)

                # Hour can be 0..23 or 99 = missing
                if hour == 99:
                    # Skip missing hour
                    continue

                # Convert to strings with zero padding if needed
                year_s  = str(year)
                month_s = f"{month:02d}"
                day_s   = f"{day:02d}"
                hour_s  = f"{hour:02d}Z"

                # Insert into data structure
                if year_s not in data_dict["data"]:
                    data_dict["data"][year_s] = {}
                if month_s not in data_dict["data"][year_s]:
                    data_dict["data"][year_s][month_s] = {}
                if day_s not in data_dict["data"][year_s][month_s]:
                    data_dict["data"][year_s][month_s][day_s] = []

                # Append the hour if not already present
                if hour_s not in data_dict["data"][year_s][month_s][day_s]:
                    data_dict["data"][year_s][month_s][day_s].append(hour_s)

        # Now fill the "years" list. Sort them numerically.
        # The keys in data_dict["data"] are strings, so convert to int for sorting.
        sorted_years = sorted(int(y) for y in data_dict["data"].keys())
        data_dict["years"] = sorted_years

        # Convert them back to strings if you want the top-level "years" 
        # to be integers or strings. The example shows them as integers, 
        # so we can just store them as int.
        # data_dict["years"] = sorted_years

        return data_dict

    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return data_dict
    except Exception as e:
        print(f"Error reading/parsing {file_name}: {e}")
        return data_dict


if __name__ == "__main__":
    # Example usage:
    station_id_example = "RSM00023955"
    output = parse_igra_station(station_id_example)

    import json
    print(json.dumps(output, indent=2))








#pressure = get_daywise_data(2024, 7)[(31,0)]['pressure']
#print(pressure)

# gph, lapse, temp = lapse_rate(gph, temp)
# tropopause = detect_tropopause(gph, lapse)


# seidel2006stationlist = [
#     'ARM00087576', 'ASM00094120', 'ASM00094294', 'ASM00094610', 'ASM00094672', 
#     'ASM00094998', 'AYM00089009', 'AYM00089050', 'AYM00089532', 'AYM00089542', 
#     'AYM00089564', 'AYM00089664', 'BDM00078016', 'BPM00091517', 'BRM00082332', 
#     'BRM00083746', 'CAM00071072', 'CAM00071082', 'CAM00071801', 'CAM00071836', 
#     'CAM00071926', 'CHM00051709', 'CHM00052681', 'CIM00085442', 'CIM00085469', 
#     'CIM00085799', 'COM00080222', 'FIM00002836', 'FJM00091680', 'FMM00091334', 
#     'FPM00091938', 'FSM00061996', 'GIM00008495', 'GLM00004360', 'GMM00010868', 
#     'HKM00045004', 'ICM00004018', 'IOM00061967', 'ISM00040179', 'IVM00065578', 
#     'JAM00047401', 'JAM00047827', 'JAM00047991', 'JNM00001001', 'KEM00063741', 
#     'LYM00062010', 'MAM00067083', 'NFM00094996', 'NGM00061052', 'NZM00093844', 
#     'NZM00093986', 'POM00008508', 'PSM00091408', 'RMM00091376', 'RQM00078526', 
#     'RSM00021504', 'RSM00021965', 'RSM00023415', 'RSM00023472', 'RSM00024266', 
#     'RSM00028698', 'RSM00030230', 'RSM00032540', 'RSM00034731', 'RSM00035121', 
#     'SAM00041024', 'SFM00068588', 'SFM00068816', 'SFM00068994', 'SGM00061641', 
#     'SHM00061902', 'SHM00068906', 'SNM00048698', 'SPM00060020', 'THM00048455', 
#     'TXM00038880', 'UKM00003005', 'USM00070026', 'USM00070273', 'USM00070308', 
#     'USM00070398', 'USM00072201', 'USM00072208', 'USM00072250', 'USM00072251', 
#     'USM00072293', 'USM00072327', 'USM00072403', 'USM00072451', 'USM00072493', 
#     'USM00072520', 'USM00072645', 'USM00072694', 'USM00072712', 'USM00072747', 
#     'USM00072768', 'USM00072776', 'USM00072797', 'USM00091165', 'USM00091285'
# ]
#plot_stations(station_list=seidel2006stationlist)