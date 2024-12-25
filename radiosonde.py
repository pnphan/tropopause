import numpy as np
import matplotlib.pyplot as plt
from helpers import lat_lon_to_indices, indices_to_lat_lon, lat_lon_to_indices2
import cartopy.crs as ccrs
from geopy.distance import geodesic
import pandas as pd
from scipy.interpolate import interp1d

stations = open("igra2-station-list.txt", "r").read().split("\n")[:-101]
temps = open("temp_00z-mly.txt", "r").read().split("\n")[:-101]
ids = [i[:11] for i in stations]
lat = [float(i[11:20].replace(" ", "")) for i in stations]
lon = [float(i[20:30].replace(" ", "")) for i in stations]
latlon = list(zip(lat,lon))
station_dict = {} # key station ID value (lat lon)
for i in range(len(ids)):
    station_dict[ids[i]] = latlon[i]
latlon_dict = {v: k for k, v in station_dict.items()}


def plot_available_stations(temp_data, station_data, desired_year, desired_month):
    stations = open(station_data, "r").read().split("\n")[:-101]
    temps = open(temp_data, "r").read().split("\n")[:-101]
    ids = [i[:11] for i in stations]
    lat = [float(i[11:20].replace(" ", "")) for i in stations]
    lon = [float(i[20:30].replace(" ", "")) for i in stations]
    latlon = list(zip(lat,lon))
    station_dict = {}
    for i in range(len(ids)):
        station_dict[ids[i]] = latlon[i]
    available_stations = []
    for i in temps:
        if float(i[12:16].replace(" ", "")) == float(desired_year) and float(i[17:19].replace(" ", "")) == float(desired_month):
            available_stations.append(i[:11])
    available_lat_lon = []
    for i in available_stations:
        available_lat_lon.append(station_dict[i])
    lat, lon = zip(*available_lat_lon)
    lat = list(lat)
    lon = list(lon)


    plt.figure(figsize=(10, 5)) 
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()  # Adds a basic map background
    plt.scatter(lon, lat, color="red", transform=ccrs.PlateCarree(), marker='.')
    plt.title("Latitude-Longitude Points on Map")
    plt.show()


def get_station(temp_data, station_data, desired_year, desired_month, desired_lat, desired_lon):
    stations = open(station_data, "r").read().split("\n")[:-101]
    temps = open(temp_data, "r").read().split("\n")[:-101]
    ids = [i[:11] for i in stations]
    lat = [float(i[11:20].replace(" ", "")) for i in stations]
    lon = [float(i[20:30].replace(" ", "")) for i in stations]
    latlon = list(zip(lat,lon))
    station_dict = {}
    for i in range(len(ids)):
        station_dict[ids[i]] = latlon[i]
    available_stations = []
    for i in temps:
        if float(i[12:16].replace(" ", "")) == float(desired_year) and float(i[17:19].replace(" ", "")) == float(desired_month):
            available_stations.append(i[:11])
    available_lat_lon = []
    for i in available_stations:
        available_lat_lon.append(station_dict[i])

    available_stations_dict = {}
    for i in available_stations:
        available_stations_dict[station_dict[i]] = i 

    target_coord = (desired_lat, desired_lon)
    nearest_station = None
    min_distance = float('inf')
    
    for coord, station_id in available_stations_dict.items():
        distance = geodesic(target_coord, coord).kilometers
        if distance < min_distance:
            min_distance = distance
            nearest_station = station_id
    
    return nearest_station

print(get_station("temp_00z-mly.txt", "igra2-station-list.txt", 2008, 1, 27.8, 97.3964))



def plot_all_stations(station_data):
    stations = open(station_data, "r").read().split("\n")[:-101]
    lat = [float(i[11:20].replace(" ", "")) for i in stations]
    lon = [float(i[20:30].replace(" ", "")) for i in stations]

    plt.figure(figsize=(10, 5)) 
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()  # Adds a basic map background
    plt.scatter(lon, lat, color="red", transform=ccrs.PlateCarree(), marker='.')
    plt.title("Latitude-Longitude Points on Map")
    plt.show()



def plot_temperature_vs_height(temp_file, height_file, station_id, year, month):
    # Read temperature data
    temp_columns = ["station_id", "year", "month", "pressure", "temperature", "ignore"]
    temp_df = pd.read_csv(temp_file, delim_whitespace=True, names=temp_columns)
    
    # Read height data
    height_columns = ["station_id", "year", "month", "pressure", "height", "ignore"]
    height_df = pd.read_csv(height_file, delim_whitespace=True, names=height_columns)
    
    # Filter data based on input parameters
    temp_filtered = temp_df[(temp_df["station_id"] == station_id) &
                            (temp_df["year"] == year) &
                            (temp_df["month"] == month)]
    height_filtered = height_df[(height_df["station_id"] == station_id) &
                                (height_df["year"] == year) &
                                (height_df["month"] == month)]
    
    # Check if there is data available for the given inputs
    if temp_filtered.empty or height_filtered.empty:
        print(f"No data available for Station ID {station_id}, Year {year}, Month {month}.")
        return
    
    # Merge temperature and height data on pressure level
    merged_df = pd.merge(temp_filtered, height_filtered, on=["station_id", "year", "month", "pressure"])
    
    # Convert temperature from tenths of degrees to degrees Celsius
    merged_df["temperature"] = merged_df["temperature"] / 10.0
    
    # Plot temperature vs geopotential height
    plt.figure(figsize=(8, 10))
    plt.plot(merged_df["temperature"], merged_df["height"], marker='o')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Geopotential Height (m)")
    plt.title(f"Temperature Profile for Station {station_id} in {year}-{month:02d}")
    #plt.gca().invert_yaxis()  # Invert y-axis to have altitude increasing upwards
    plt.grid(True)
    plt.show()


def plot_temperature_profile_with_lapse_rate(temp_file, height_file, station_id, year, month):
    # Read temperature data
    temp_columns = ["station_id", "year", "month", "pressure", "temperature", "ignore"]
    temp_df = pd.read_csv(temp_file, delim_whitespace=True, names=temp_columns)
    
    # Read height data
    height_columns = ["station_id", "year", "month", "pressure", "height", "ignore"]
    height_df = pd.read_csv(height_file, delim_whitespace=True, names=height_columns)
    
    # Filter data based on input parameters
    temp_filtered = temp_df[(temp_df["station_id"] == station_id) &
                            (temp_df["year"] == year) &
                            (temp_df["month"] == month)]
    height_filtered = height_df[(height_df["station_id"] == station_id) &
                                (height_df["year"] == year) &
                                (height_df["month"] == month)]
    
    # Check if there is data available for the given inputs
    if temp_filtered.empty or height_filtered.empty:
        print(f"No data available for Station ID {station_id}, Year {year}, Month {month}.")
        return
    
    # Merge temperature and height data on pressure level
    merged_df = pd.merge(temp_filtered, height_filtered, on=["station_id", "year", "month", "pressure"])
    
    # Convert temperature from tenths of degrees to degrees Celsius
    merged_df["temperature"] = merged_df["temperature"] / 10.0
    
    # Linear interpolation for temperature profile
    interp_func = interp1d(merged_df["height"], merged_df["temperature"], kind='linear', fill_value='extrapolate')
    height_interp = np.linspace(merged_df["height"].min(), merged_df["height"].max(), 100)
    temperature_interp = interp_func(height_interp)
    
    # Calculate lapse rate (change in temperature per change in height in km)
    lapse_rate = -np.gradient(temperature_interp, height_interp / 1000.0)
    
    # Plot temperature profile and lapse rate side by side
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot temperature profile
    axes[0].plot(merged_df["temperature"], merged_df["height"], 'o', label='Original Data')
    axes[0].plot(temperature_interp, height_interp, '-', label='Interpolated Profile')
    #axes[0].invert_yaxis()  # Invert y-axis to have altitude increasing upwards
    axes[0].set_xlabel("Temperature (°C)")
    axes[0].set_ylabel("Geopotential Height (m)")
    axes[0].set_title(f"Temperature Profile for Station {station_id} in {year}-{month:02d}")
    axes[0].grid(True)
    axes[0].legend()
    
    # Plot lapse rate
    axes[1].plot(lapse_rate, height_interp, '-r')
    #axes[1].invert_yaxis()  # Invert y-axis to have altitude increasing upwards
    axes[1].set_xlabel("Lapse Rate (°C/km)")
    axes[1].set_title(f"Lapse Rate for Station {station_id} in {year}-{month:02d}")
    axes[1].grid(True)
    
    plt.tight_layout()
    plt.show()


def find_nearest_station(latlon_dict, station_dict, latitude, longitude):
    target_coord = (latitude, longitude)
    nearest_station = None
    min_distance = float('inf')
    
    for coord, station_id in latlon_dict.items():
        distance = geodesic(target_coord, coord).kilometers
        if distance < min_distance:
            min_distance = distance
            nearest_station = station_id
    
    return (nearest_station, station_dict[nearest_station])


def get_available_years_months(station_id, temp_file='temp_00z-mly.txt', height_file='ghgt_00z-mly.txt'):
    # Read temperature data
    temp_columns = ["station_id", "year", "month", "pressure", "temperature", "ignore"]
    temp_df = pd.read_csv(temp_file, delim_whitespace=True, names=temp_columns)
    
    # Read height data
    height_columns = ["station_id", "year", "month", "pressure", "height", "ignore"]
    height_df = pd.read_csv(height_file, delim_whitespace=True, names=height_columns)
    
    # Combine both temperature and height data
    combined_df = pd.concat([temp_df, height_df])
    
    # Filter data based on station ID
    filtered_df = combined_df[combined_df["station_id"] == station_id]
    
    # Get unique years and months
    available_data = filtered_df.drop_duplicates(subset=["year", "month"])[["year", "month"]]
    
    available_dict = {}
    for _, row in available_data.iterrows():
        year = row["year"]
        month = row["month"]
        if year not in available_dict:
            available_dict[year] = []
        available_dict[year].append(month)
    
    return available_dict
    #return available_data



print(get_available_years_months("USM00072251"))

# nearest_station = find_nearest_station(latlon_dict, station_dict, 27.8, 97.3964)[0]
# #print(get_available_years_months("temp_00z-mly.txt", "ghgt_00z-mly.txt", nearest_station))
# #plot_temperature_vs_height("temp_00z-mly.txt", "ghgt_00z-mly.txt", "USM00072251", 2008, 2)
# plot_temperature_profile_with_lapse_rate("temp_00z-mly.txt", "ghgt_00z-mly.txt", "USM00072251", 2008, 5)
