from meng2021replication import *


def H_mean_year(stationid, year):
    """
    Compute the mean of tropopause height H for a single station and single year
    """
    with open(f"monthly_means/{stationid}.pkl", "rb") as f:
        data = pickle.load(f)
    years = list(data.keys())
    year = str(year)
    if year in years:
        return np.nanmean(list(data[year].values()))
    else:
        return -9999
    

def H_mean_all_years(stationid):
    """
    Compute the yearly mean of tropopause height H for a single station, for all years
    and saves the file 
    """
    with open(f"monthly_means/{stationid}.pkl", "rb") as f:
        data = pickle.load(f)
    years = list(data.keys())
    years = [str(i) for i in years]

    output = {}
    for year in years:
        data = H_mean_year(stationid, year)
        if data != -9999:
            output[year] = data

    with open(f"yearly_means/{stationid}_yearly_mean.pkl", "wb") as file:
        pickle.dump(output, file)


def H_mean_all_years_all_stations(station_list=station_list):
    """
    Compute the yearly mean of tropopause height H for a single station, for all years, for all stations

    Saves the yearly means 
    """
    for station in station_list:
        H_mean_all_years(station)
        print(station)


def H_mean_yearly_combined(station_list=station_list):
    """
    Compute the global yearly mean for all years

    Returns 
    
    x : list of years
    y : list of means for the corresponding year
    """
    x = []
    y = []

    for year in range(1980, 2025):
        current_y = []
        for station in station_list:
            with open(f"yearly_means/{station}_yearly_mean.pkl", "rb") as f:
                data = pickle.load(f)

            if str(year) in list(data.keys()):
                current_y.append(data[str(year)])

        current_y_avg = np.nanmean(current_y)
        y.append(current_y_avg)
        x.append(year)

    return x, y

        
def H_anomaly_yearly_combined(station_list=station_list):
    """
    Compute the global yearly anomaly over 1980-2024

    Returns 

    x : list of years
    y : list of anomalies for the corresponding year, in meters
    """
    x, y = H_mean_yearly_combined(station_list)
    y = y - np.nanmean(y)
    y = y*1000
    return x, y

        

x, y = H_anomaly_yearly_combined()

plot_values(x, y)
