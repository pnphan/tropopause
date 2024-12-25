import xarray as xr
import numpy as np
import cupy as cp
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
import cdsapi
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def cds_request():
    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": [
            "geopotential",
            "temperature"
        ],
        "pressure_level": [
            "1", "2", "3",
            "5", "7", "10",
            "20", "30", "50",
            "70", "100", "125",
            "150", "175", "200",
            "225", "250", "300",
            "350", "400", "450",
            "500", "550", "600",
            "650", "700", "750",
            "775", "800", "825",
            "850", "875", "900",
            "925", "950", "975",
            "1000"
        ],
        "year": [
            "1940", "1941", "1942",
            "1943", "1944", "1945",
            "1946", "1947", "1948",
            "1949", "1950", "1951",
            "1952", "1953", "1954",
            "1955", "1956", "1957",
            "1958", "1959", "1960",
            "1961", "1962", "1963",
            "1964", "1965", "1966",
            "1967", "1968", "1969",
            "1970", "1971", "1972",
            "1973", "1974", "1975",
            "1976", "1977", "1978",
            "1979", "1980"
        ],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "time": ["00:00"],
        "data_format": "netcdf",
        "download_format": "unarchived"
    }

    client = cdsapi.Client()
    client.retrieve(dataset, request).download()


def get_date_index_1980(month, year):
    return (year - 1980) * 12 + month - 1

def get_date_index_1940(month, year):
    return (year - 1940) * 12 + month - 1

def open_as_np_1940(file_name="monthly-avg-temp-gpheight-1940-1979.nc", variable='t', date_slice=(0,10)):
    """
    Opens the data as a numpy array
    """
    ds = xr.open_dataset(file_name)
    if date_slice != False:
        array = ds[variable].isel(valid_time=slice(date_slice[0],date_slice[1]))
    else:
        array = ds[variable].isel()
    return np.asarray(array)

def open_as_np_1980(file_name="monthly-avg-temp-gpheight-1980-2023.nc", variable='t', date_slice=(0,10)):
    """
    Opens the data as a numpy array
    """
    ds = xr.open_dataset(file_name)
    if date_slice != False:
        array = ds[variable].isel(date=slice(date_slice[0],date_slice[1]))
    else:
        array = ds[variable].isel()
    return np.asarray(array)


def get_geopotential_height_gpu(array):
    return array / 9.80665


def interpolate_to_100m_cubic(temp, geopotential, height_resolution=100):
    """
    Interpolates temperature data to a uniform vertical grid of 100m resolution using cubic interpolation.

    Parameters:
    temp : Temperature array (time, levels, lat, lon).
    geopotential : Geopotential height array (time, levels, lat, lon).
    height_resolution (int): Vertical resolution in meters (default: 100m).

    Returns:
    Interpolated temperature array on a uniform height grid.
    """
    # Ensure arrays are numpy arrays
    temp = np.asarray(temp)
    geopotential = np.asarray(geopotential)
    
    # Define the target vertical grid (100m resolution)
    max_height = np.nanmax(geopotential)  # Maximum height in the dataset
    target_heights = np.arange(0, max_height + height_resolution, height_resolution)
    
    # Prepare an empty array for the interpolated temperature
    interpolated_temp = np.empty((temp.shape[0], len(target_heights), temp.shape[2], temp.shape[3]))
    
    # Loop over each spatial grid point
    for k in range(temp.shape[0]):
        for i in range(temp.shape[2]):  # Latitude
            for j in range(temp.shape[3]):  # Longitude
                # Extract 1D profiles for temperature and geopotential
                temp_profile = temp[k, :, i, j]
                geopotential_profile = geopotential[k, :, i, j]
                
                # Skip if the profile contains NaNs
                if np.any(np.isnan(temp_profile)) or np.any(np.isnan(geopotential_profile)):
                    interpolated_temp[k, :, i, j] = np.nan
                    continue
                
                # Create a cubic spline interpolator
                interpolator = interp1d(
                    geopotential_profile, temp_profile, kind="cubic", bounds_error=False, fill_value="extrapolate"
                )
                
                # Interpolate to the target heights
                interpolated_temp[k, :, i, j] = interpolator(target_heights)
    
    return interpolated_temp #, target_heights


def calculate_lapse_rate(temp):
    """
    Calculate the lapse rate given temperature on 100m levels.
    
    Parameters:
    temp : Temperature array (time, levels, lat, lon) in K.
    geopotential : Geopotential array (time, levels, lat, lon) in m.
    
    Returns:
    Lapse rate array (time, levels-1, lat, lon) in K/km.
    """
    lapse_rate = -np.diff(temp, axis=1) * 10

    return lapse_rate


def convert_to_celsius(temp):
    """
    Parameters:
    temp : (time, levels, lat, lon) in K

    Returns:
    temp : (time, levels, lat, lon) in C
    """
    return temp - 273.15


def calculate_tropopause_height(lapse_rates, height_resolution=100, lapse_rate_threshold=2.0):
    """
    Calculate the tropopause height based on lapse rates.

    Parameters:
    lapse_rates (numpy.ndarray): Lapse rate array of shape (time, height, latitude, longitude).
                                 Units: degrees/km.
    height_resolution (float): Vertical resolution of the height levels in meters (default: 100m).
    lapse_rate_threshold (float): Threshold for the lapse rate to define the tropopause (default: 2 degrees/km).

    Returns:
    numpy.ndarray: Tropopause height array of shape (time, 1, latitude, longitude) in meters.
    """
    lapse_rates = lapse_rates[:,50:231,:,:]

    # Get the shape of the input array
    time, height_levels, lat, lon = lapse_rates.shape

    # Compute heights corresponding to the levels
    heights = np.arange(5000, 23100, height_resolution)

    # Prepare an array to store the tropopause heights
    tropopause_heights = np.full((time, 1, lat, lon), np.nan, dtype=np.float32)

    # Iterate over each grid point and time step
    for t in range(time):
        for i in range(lat):
            for j in range(lon):
                # Extract the lapse rate profile for the grid point
                lapse_rate_profile = lapse_rates[t, :, i, j]
                
                # Initialize the tropopause height
                tropopause_height = np.nan

                # Iterate over height levels
                for k in range(height_levels):
                    # Check if lapse rate is less than or equal to the threshold
                    if lapse_rate_profile[k] <= lapse_rate_threshold:
                        # Calculate the average lapse rate within the next 2 km
                        next_2km_levels = k + int(2000 / height_resolution)
                        if next_2km_levels >= height_levels:
                            next_2km_levels = height_levels - 1
                        
                        average_lapse_rate = np.mean(lapse_rate_profile[k:next_2km_levels + 1])
                        
                        # Check if the average lapse rate is less than or equal to the threshold
                        if average_lapse_rate <= lapse_rate_threshold:
                            tropopause_height = heights[k]
                            break
                
                # Store the tropopause height
                tropopause_heights[t, :, i, j] = tropopause_height

    return tropopause_heights



def detect_second_tropopause(lapse_rates, first_tropopause_heights, height_resolution=100):
    """
    Detect the second tropopause in a 3D lapse rate array.

    Parameters:
    lapse_rates (numpy.ndarray): Lapse rate array of shape (505, 721, 1440) in degrees/km.
    first_tropopause_heights (numpy.ndarray): Array of first tropopause heights (721, 1440) in meters.
    height_resolution (int): Vertical resolution of height levels in meters (default: 100m).

    Returns:
    numpy.ndarray: Second tropopause height array of shape (721, 1440) in meters. 
                   NaN if no second tropopause is found.
    """
    # Get the vertical levels
    num_levels = lapse_rates.shape[0]
    heights = np.arange(0, num_levels * height_resolution, height_resolution)
    
    # Initialize the second tropopause height array with NaN
    second_tropopause_heights = np.full((721, 1440), np.nan, dtype=np.float32)
    
    # Loop over all latitude and longitude points
    for i in range(lapse_rates.shape[1]):  # Latitude
        print(i)
        for j in range(lapse_rates.shape[2]):  # Longitude
            # Extract the lapse rate profile and first tropopause height
            lapse_rate_profile = lapse_rates[:, i, j]
            first_tropopause_height = first_tropopause_heights[i, j]
            
            # Skip grid points where the first tropopause is not defined``
            if np.isnan(first_tropopause_height):
                continue
            
            # Find the index of the first tropopause
            first_tropopause_index = int(first_tropopause_height / height_resolution)
            
            # Check for a region above the first tropopause where the average lapse rate
            # exceeds 3 degrees/km within 1 km
            for k in range(first_tropopause_index + 1, num_levels):
                # Calculate the average lapse rate within 1 km above the current level
                upper_limit = min(k + int(1000 / height_resolution), num_levels)
                avg_lapse_rate = np.mean(lapse_rate_profile[k:upper_limit])
                
                # If the average lapse rate exceeds 3 degrees/km, search for a second tropopause
                if avg_lapse_rate > 3.0:
                    for m in range(k, num_levels):
                        # Apply the same criterion as the first tropopause for the second tropopause
                        if lapse_rate_profile[m] <= 2.0:
                            # Calculate the average lapse rate within 2 km above this level
                            upper_2km_limit = min(m + int(2000 / height_resolution), num_levels)
                            avg_2km_lapse_rate = np.mean(lapse_rate_profile[m:upper_2km_limit])
                            
                            # If the average lapse rate within 2 km is also <= 2.0, record the height
                            if avg_2km_lapse_rate <= 2.0:
                                second_tropopause_heights[i, j] = heights[m]
                                break
                    break  # Stop searching once a second tropopause is found
    
    return second_tropopause_heights



def calculate_second_tropopause_height(lapse_rates, first_tropopause, height_resolution=100, lapse_rate_threshold=2.0):
    """
    Calculate the tropopause height based on lapse rates.

    Parameters:
    lapse_rates (numpy.ndarray): Lapse rate array of shape (time, height, latitude, longitude).
                                 Units: degrees/km.
    height_resolution (float): Vertical resolution of the height levels in meters (default: 100m).
    lapse_rate_threshold (float): Threshold for the lapse rate to define the tropopause (default: 2 degrees/km).

    Returns:
    numpy.ndarray: Tropopause height array of shape (time, 1, latitude, longitude) in meters.
    """
    #lapse_rates = lapse_rates[:,50:231,:,:]

    # Get the shape of the input array
    time, height_levels, lat, lon = lapse_rates.shape

    # Compute heights corresponding to the levels
    heights = np.arange(0, height_levels*100, height_resolution)

    # Prepare an array to store the tropopause heights
    second_tropopause_heights = np.full((time, 1, lat, lon), np.nan, dtype=np.float32)

    for t in range(time):
        for i in range(lat):
            print(i)
            for j in range(lon):
                lapse_rate_profile = lapse_rates[t, :, i, j]
                check=False
                for k in range(int(first_tropopause[t,0,i,j]/100), height_levels):
                    next_1km_levels = k + int(1000 / height_resolution)
                    if next_1km_levels >= height_levels:
                        next_1km_levels = height_levels - 1
                    if np.mean(lapse_rate_profile[k:next_1km_levels + 1]) >= 3:
                        check=k
                        break
                if check != False:
                    for m in range(check, height_levels):
                    # Check if lapse rate is less than or equal to the threshold
                        if lapse_rate_profile[m] <= lapse_rate_threshold:
                            # Calculate the average lapse rate within the next 2 km
                            next_2km_levels = m + int(2000 / height_resolution)
                            if next_2km_levels >= height_levels:
                                next_2km_levels = height_levels - 1
                            
                            average_lapse_rate = np.mean(lapse_rate_profile[m:next_2km_levels + 1])
                            
                            # Check if the average lapse rate is less than or equal to the threshold
                            if average_lapse_rate <= lapse_rate_threshold:
                                second_tropopause_height = heights[m]
                                break
                            else:
                                second_tropopause_height = np.nan
                    
                    # Store the tropopause height
                    second_tropopause_heights[t, :, i, j] = second_tropopause_height

    return second_tropopause_heights



def calculate_rolling_average(data, window=10):
    """
    Calculate the rolling average over the first axis with a specified window size.

    Parameters:
    - data: np.ndarray
        Input array with shape (505, 721, 1440).
    - window: int
        The size of the window for averaging (default is 10).

    Returns:
    - result: np.ndarray
        Output array with shape (505 - window + 1, 721, 1440).
    """
    if data.shape[0] < window:
        raise ValueError("The first dimension of the array must be at least as large as the window size.")
    
    # Use a sliding window to compute the mean
    result = np.empty((data.shape[0] - window + 1, data.shape[1], data.shape[2]))
    
    for i in range(result.shape[0]):
        result[i] = np.mean(data[i:i + window], axis=0)
    
    return result


def plot_tropopause(data, title="Tropopause Height", 
                          xlabel="Longitude", ylabel="Latitude"):
    """
    Plot a contour map from a (1, 721, 1440) array with a world map underlay.

    Parameters:
    - data: np.ndarray
        A (1, 721, 1440) array to be plotted.
    - title: str
        Title of the contour map.
    - xlabel: str
        Label for the x-axis.
    - ylabel: str
        Label for the y-axis.
    """
    # Remove the first dimension
    if data.shape == (1, 721, 1440):
        data_2d = data[0]
    else:
        data_2d = data
    
    # Generate latitude and longitude values
    lat = np.linspace(-90, 90, data_2d.shape[0])  # 721 points
    lon = np.linspace(-180, 180, data_2d.shape[1])  # 1440 points
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6),
                           subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Add the world map
    ax.add_feature(cfeature.LAND, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    
    # Plot the contour map
    contour = ax.contourf(lon, lat, data_2d, levels=50, cmap="viridis", 
                          transform=ccrs.PlateCarree())
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label("Value")
    
    # Add title and labels
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    # Set gridlines
    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
    
    plt.show()


def plot_vertical_temperature_profile(data, lat_index, lon_index, time=0,
                                      title="Vertical Temperature Profile",
                                      xlabel="Temperature (°C)", ylabel="Altitude"):
    """
    Plot the vertical temperature profile for specific latitude and longitude indices.

    Parameters:
    - data: np.ndarray
        A (506, 721, 1440) array representing temperature data.
    - lat_index: int
        The index of the latitude (0 to 720).
    - lon_index: int
        The index of the longitude (0 to 1439).
    - title: str
        Title of the plot.
    - xlabel: str
        Label for the x-axis.
    - ylabel: str
        Label for the y-axis.
    """
    if not (0 <= lat_index < 721):
        raise ValueError("Latitude index must be in the range [0, 720].")
    
    if not (0 <= lon_index < 1440):
        raise ValueError("Longitude index must be in the range [0, 1439].")
    
    # Extract the temperature profile at the specified location
    vertical_profile = data[time, :, lat_index, lon_index]
    
    # Altitude levels (assuming they are represented by index 0-505)
    altitude_levels = np.arange(0, 50600, 100)  # 506 levels, starting from 1
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(vertical_profile, altitude_levels)#, marker='o', linestyle='-')
    #plt.gca().invert_yaxis()  # Invert y-axis (higher altitudes at the top)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()


def plot_vertical_lapse_profile(data, lat_index, lon_index, time=0,
                                      title="Vertical Lapse Rate Profile",
                                      xlabel="Lapse Rate (°C/km)", ylabel="Altitude"):
    """
    Plot the vertical temperature profile for specific latitude and longitude indices.

    Parameters:
    - data: np.ndarray
        A (506, 721, 1440) array representing temperature data.
    - lat_index: int
        The index of the latitude (0 to 720).
    - lon_index: int
        The index of the longitude (0 to 1439).
    - title: str
        Title of the plot.
    - xlabel: str
        Label for the x-axis.
    - ylabel: str
        Label for the y-axis.
    """
    data = calculate_lapse_rate(data)


    if not (0 <= lat_index < 721):
        raise ValueError("Latitude index must be in the range [0, 720].")
    
    if not (0 <= lon_index < 1440):
        raise ValueError("Longitude index must be in the range [0, 1439].")
    
    # Extract the temperature profile at the specified location
    vertical_profile = data[time, :, lat_index, lon_index]
    
    # Altitude levels (assuming they are represented by index 0-505)
    altitude_levels = np.arange(0, 50500, 100)  # 506 levels, starting from 1
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(vertical_profile, altitude_levels)#, marker='o', linestyle='-')
    #plt.gca().invert_yaxis()  # Invert y-axis (higher altitudes at the top)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()


def lat_lon_to_indices(lat, lon, n_lat=721, n_lon=1440):
    """
    Convert latitude and longitude to flat (721, 1440) array indices.

    Parameters:
    - lat: float
        Latitude in degrees (-90 to 90).
    - lon: float
        Longitude in degrees (-180 to 180).
    - n_lat: int
        Number of latitude points (default is 721).
    - n_lon: int
        Number of longitude points (default is 1440).

    Returns:
    - (lat_index, lon_index): tuple of int
        Indices for the latitude and longitude in the array.
    """
    if not (-90 <= lat <= 90):
        raise ValueError("Latitude must be in the range [-90, 90].")
    if not (-180 <= lon <= 180):
        raise ValueError("Longitude must be in the range [-180, 180].")
    
    # Map latitude and longitude to indices
    lat_index = int(round((lat + 90) / 180 * (n_lat - 1)))  # Latitude index
    lon_index = int(round((lon + 180) / 360 * (n_lon - 1)))  # Longitude index
    
    return lat_index, lon_index


def indices_to_lat_lon(lat_index, lon_index, n_lat=721, n_lon=1440):
    """
    Convert flat (721, 1440) array indices to latitude and longitude.

    Parameters:
    - lat_index: int
        Latitude index in the array (0 to n_lat-1).
    - lon_index: int
        Longitude index in the array (0 to n_lon-1).
    - n_lat: int
        Number of latitude points (default is 721).
    - n_lon: int
        Number of longitude points (default is 1440).

    Returns:
    - (lat, lon): tuple of float
        Latitude and longitude in degrees.
    """
    if not (0 <= lat_index < n_lat):
        raise ValueError(f"Latitude index must be in the range [0, {n_lat - 1}].")
    if not (0 <= lon_index < n_lon):
        raise ValueError(f"Longitude index must be in the range [0, {n_lon - 1}].")
    
    # Convert indices to latitude and longitude
    lat = lat_index / (n_lat - 1) * 180 - 90  # Latitude in degrees
    lon = lon_index / (n_lon - 1) * 360 - 180  # Longitude in degrees
    
    return lat, lon


def get_tropopause(month, year):
    """
    Returns global tropopause height map given month (1-12) and year (1940-2023) 
    """
    if year >= 1980:
        index = get_date_index_1980(month, year)
        temp = open_as_np_1980(file_name="monthly-avg-temp-gpheight-1980-2023.nc" , variable='t', date_slice=(index,index+1))
        geopot = open_as_np_1980(file_name="monthly-avg-temp-gpheight-1980-2023.nc", variable='z', date_slice=(index,index+1))
        geopot_height = get_geopotential_height_gpu(geopot)
        temp = interpolate_to_100m_cubic(temp, geopot_height)
        lapse = calculate_lapse_rate(temp)
        tropopause = calculate_tropopause_height(lapse)
    else:
        index = get_date_index_1940(month, year)
        temp = open_as_np_1940(file_name="monthly-avg-temp-gpheight-1940-1979.nc" , variable='t', date_slice=(index,index+1))
        geopot = open_as_np_1940(file_name="monthly-avg-temp-gpheight-1940-1979.nc", variable='z', date_slice=(index,index+1))
        geopot_height = get_geopotential_height_gpu(geopot)
        temp = interpolate_to_100m_cubic(temp, geopot_height)
        lapse = calculate_lapse_rate(temp)
        tropopause = calculate_tropopause_height(lapse)
    return tropopause


def plot_zonal_mean(data):
    """
    Plot the zonal mean (mean across longitudes for each latitude).
    
    Parameters:
    data (numpy.ndarray): Input array of shape (1, 721, 1440).
    latitudes (numpy.ndarray): Array of latitude values of shape (721,).
    """
    # Calculate the zonal mean (average over the longitude dimension)
    zonal_mean = np.mean(data[0,:,:], axis=-1)
    

    lat = np.linspace(-90, 90, data[0].shape[0])  # 721 points
    # Plot the zonal mean
    plt.figure(figsize=(8, 5))
    plt.plot(lat, zonal_mean, label="Zonal Mean")
    plt.xlabel("Latitude")
    plt.ylabel("Mean Value")
    plt.title("Zonal Mean")
    plt.grid()
    plt.legend()
    plt.xlim(-90,90)
    plt.xticks(ticks=np.arange(-90, 91, 30), labels=[f"{lat}°" for lat in np.arange(-90, 91, 30)])
    plt.show()


def plot_meridional_mean(data):
    """
    Plot the meridional mean (mean across latitudes for each longitude).
    
    Parameters:
    data (numpy.ndarray): Input array of shape (1, 721, 1440).
    longitudes (numpy.ndarray): Array of longitude values of shape (1440,).
    """
    # Calculate the meridional mean (average over the latitude dimension)
    meridional_mean = np.mean(data[0, :, :], axis=-2)
    lon = np.linspace(-180, 180, data[0].shape[1])  # 1440 points
    # Plot the meridional mean
    plt.figure(figsize=(8, 5))
    plt.plot(lon, meridional_mean, label="Meridional Mean")
    plt.xlabel("Longitude")
    plt.ylabel("Mean Value")
    plt.title("Meridional Mean")
    plt.grid()
    plt.legend()
    plt.xlim(-180,180)
    plt.xticks(ticks=np.arange(-180, 181, 60), labels=[f"{lon}°" for lon in np.arange(-180, 181, 60)])
    plt.show()






#np.save("ex.npy", get_tropopause(7, 1981))
# tropopause = np.load("ex.npy")
# plot_tropopause(tropopause[0], title="Tropopause Height, July 1981")
# print(np.shape(tropopause))
# plot_meridional_mean(tropopause[0])


# temp = open_as_np_1940(file_name="monthly-avg-temp-gpheight-1940-1979.nc" , variable='t', date_slice=(0,1))
# # print(np.shape(temp))
# # temp = convert_to_celsius(temp)
# # print(temp[0,18,360,720])
# geopot = open_as_np_1940(variable='z', date_slice=(0,1))
# geopot_height = get_geopotential_height_gpu(geopot)
# # print(geopot_height[0,18,360,720])
# temp = interpolate_to_100m_cubic(temp, geopot_height)
temp = np.load("interpolatedtempexample.npy")
lapse = calculate_lapse_rate(temp)
tropopause = calculate_tropopause_height(lapse)
double_tropopause = detect_second_tropopause(lapse[0], tropopause[0,0])
np.save("double.npy", double_tropopause)
plot_tropopause(double_tropopause)

# compute_tropopause_all()

#cds_request()


