import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import cartopy.crs as ccrs
from geopy.distance import geodesic
from parameters import *


def H1_monthly(station_list=station_list):
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
            y.append(np.nanmean(current_month))
            counter+=1
            x.append(counter)
    
    return x, y

def H2_monthly(station_list=station_list):
    y = []
    x = []
    counter=0
    for year in range(1980,2025):
        for month in months:
            current_month = []
            for station in station_list:
                with open(f"./H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month))
            counter+=1
            x.append(counter)
    
    return x, y


if __name__ == "__main__":
    x, y = H1_monthly()
    x1, y1 = H2_monthly()
    years = range(1980,2025)

    x_ticks = [i * 12 for i in range(len(years))]
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    alpha=0.01
    slope_ci_lower, slope_ci_upper, significant = t_test(years, slope, intercept, r_value, p_value, std_err, alpha)

    plt.figure(figsize=(6, 3))
    fig, ax = plt.subplots()
    plt.plot(x, y, label='First Tropopause')
    plt.plot(x, slope*np.array(x) + intercept, label='First Tropopause Trend')#, 'r')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
    alpha=0.01
    slope_ci_lower, slope_ci_upper, significant = t_test(years, slope, intercept, r_value, p_value, std_err, alpha)
    plt.plot(x1, y1, label='Second Tropopause')
    plt.plot(x1, slope*np.array(x1) + intercept, label='Second Tropopause Trend')#, 'r')

    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    plt.xlabel('Year')
    plt.ylabel('Height (km)')
    #plt.grid(True)
    #plt.show()
    plt.title('Monthly First and Second Tropopause Height 1980-2024')
    plt.legend(fontsize='x-small', bbox_to_anchor=(1, 1))
    fig.subplots_adjust(right=0.75)
    plt.savefig('./figures/figures_output/both_monthly.png', dpi=1200, bbox_inches='tight')










    y2 = np.array(y1) - np.array(y)
    years = range(1980,2025)

    x_ticks = [i * 12 for i in range(len(years))]
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y2)
    alpha=0.01
    slope_ci_lower, slope_ci_upper, significant = t_test(years, slope, intercept, r_value, p_value, std_err, alpha)

    plt.figure(figsize=(6, 3))
    plt.plot(x, y2, label='Monthly Tropopause Thickness')
    plt.plot(x, slope*np.array(x) + intercept, label='Monthly Tropopause Thickness Trend')#, 'r')

    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    plt.xlabel('Year')
    plt.ylabel('Height (km)')
    #plt.grid(True)
    #plt.show()
    plt.title('Monthly Tropopause Thickness 1980-2024')
    plt.savefig('./figures/figures_output/both_monthly_difference.png', dpi=1200, bbox_inches='tight')
