import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import cartopy.crs as ccrs
from geopy.distance import geodesic
from parameters import *

def H1_yearly_single_month_anomalies(month, station_list=station_list):
    monthly_means = {}
    current_month = []
    for year in range(1980,2024):
        current_year = []
        for station in station_list:
            with open(f"./H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        current_month.append(np.nanmean(current_year))
    monthly_means[month] = np.nanmean(current_month)


    y = []
    x = []
    for year in range(1980,2025):
        current_year = []
        for station in station_list:
            with open(f"./H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        y.append(np.nanmean(current_year) - monthly_means[month])
        x.append(year)

    return x, y

def H2_yearly_single_month_anomalies(month, station_list=station_list):
    monthly_means = {}
    current_month = []
    for year in range(1980,2024):
        current_year = []
        for station in station_list:
            with open(f"./H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        current_month.append(np.nanmean(current_year))
    monthly_means[month] = np.nanmean(current_month)


    y = []
    x = []
    for year in range(1980,2025):
        current_year = []
        for station in station_list:
            with open(f"./H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        y.append(np.nanmean(current_year) - monthly_means[month])
        x.append(year)
    
    return x, y


if __name__ == "__main__":

    for month in months:
        x, y = H1_yearly_single_month_anomalies(month)
        x1, y1 = H2_yearly_single_month_anomalies(month)

        years = range(1980,2025)

        x_ticks = range(1980, 2025)
        x_tick_labels = [str(year) if year%5==0 else '' for year in years]
        

        A = np.vstack([x, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y)[0]

        A1 = np.vstack([x1, np.ones(len(x1))]).T
        m1, c1 = np.linalg.lstsq(A1, y1)[0]


        plt.figure(figsize=(6, 3))
        plt.plot(x, y, label='First Tropopause')
        plt.plot(x1, y1, label='Second Tropopause')
        plt.plot(x, m*np.array(x) + c, label='First Tropopause Trend')
        plt.plot(x1, m1*np.array(x1) + c1, label='Second Tropopause Trend')
        
        plt.xticks(x_ticks, x_tick_labels, rotation=45)

        plt.xlabel('Year')
        plt.ylabel('Anomaly (km)')
        #plt.grid(True)
        #plt.show()
        plt.ylim((-1, 1))
        plt.legend(loc='upper left', fontsize='x-small')
        plt.title(f'Yearly {months_dict[month]} Anomaly of First and Second Tropopause Height')
        plt.savefig(f'./figures/figures_output/both_yearly_{month}_{months_dict[month]}_anomalies.png', dpi=1200, bbox_inches='tight')