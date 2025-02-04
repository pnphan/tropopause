import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import cartopy.crs as ccrs
from geopy.distance import geodesic
from parameters import *


def H2_yearly_single_month_anomalies_combined(month, station_list=station_list):
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
    years = range(1980,2025)

    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    alpha=0.01
    slope_ci_lower, slope_ci_upper, significant = t_test(years, slope, intercept, r_value, p_value, std_err, alpha)
    # A = np.vstack([x, np.ones(len(x))]).T
    # m, c = np.linalg.lstsq(A, y)[0]
    plt.plot(x, y, label=f'{months_dict[month]}', lw=0.5)
    plt.plot(x, slope*np.array(x) + intercept, label=f'{months_dict[month]} Trend', lw=0.75)
    
    #return x, y


months = ['01', '07']
if __name__ == "__main__":
    plt.figure(figsize=(6, 3))
    fig, ax = plt.subplots()
    for month in months:
        H2_yearly_single_month_anomalies_combined(month)

    years = range(1980,2025)
    x_ticks = range(1980, 2025)
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    plt.xticks(x_ticks, x_tick_labels, rotation=45)

    plt.xlabel('Year')
    plt.ylabel('Anomaly (km)')
    #plt.grid(True)
    #plt.show()
    plt.ylim((-1, 1))
    ax.legend(fontsize='xx-small', bbox_to_anchor=(1, 1))
    fig.subplots_adjust(right=0.75)
    plt.title(f'Yearly Anomaly of Second Tropopause Height by Month')
    plt.savefig(f'./figures/figures_output/H2_yearly_{month[0]}_{month[1]}_anomalies_combined.png', dpi=1200, bbox_inches='tight')