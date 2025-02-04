import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import cartopy.crs as ccrs
from geopy.distance import geodesic
from parameters import *

def H1_monthly_anomalies(station_list=station_list):
    monthly_means = {}
    for month in months:
        current_month = []
        for year in range(1980,2025):
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
    counter=0
    for year in range(1980,2025):
        for month in months:
            current_month = []
            for station in station_list:
                with open(f"./H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month) - monthly_means[month])
            counter+=1
            x.append(counter)


    years = range(1980,2025)

    x_ticks = [i * 12 for i in range(len(years))]
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    alpha=0.01
    slope_ci_lower, slope_ci_upper, significant = t_test(years, slope, intercept, r_value, p_value, std_err, alpha)

    plt.figure(figsize=(6, 3))
    plt.plot(x, y)
    plt.plot(x, slope*np.array(x) + intercept, 'r')
    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    plt.xlabel('Year')
    plt.ylabel('Anomaly (km)')
    #plt.grid(True)
    #plt.show()
    plt.title('Monthly Anomaly of First Tropopause Height')
    plt.savefig('./figures/figures_output/H1_monthly_anomalies_v3.png', dpi=1200, bbox_inches='tight')
    return x, y, significant


if __name__ == "__main__":
    print(H1_monthly_anomalies())

