from meng2021replication import *


def daily_H2_H1_mean(station_id):
    """
    Compute the daily means of tropopause thickness (if available) for a single station

    Returns a dict

    updated_data : {year : 
                        {month : 
                            {day : H2} 
                        }
                    }
    """
    with open(f"tropopauses/{station_id}.pkl", "rb") as file:
        data = pickle.load(file)
    data = data['thickness']
    years = list(data.keys())
    updated_data = {}
    for year in years:
        updated_data[year] = {}
        months = list(data[year].keys())
        for month in months:
            updated_data[year][month] = {}
            days_times = list(data[year][month].keys())
            days_times = tuples_to_dict(days_times)
            for day in list(days_times.keys()):
                times = find_pair_in_ranges(days_times[day], (9, 15), (21, 23), (0, 3))
                if times != None:
                    tropopause_time1 = data[year][month][(int(day), times[0])]
                    tropopause_time2 = data[year][month][(int(day), times[1])]
                    updated_data[year][month][day] = (tropopause_time1 + tropopause_time2)/2
                else:
                    updated_data[year][month][day] = data[year][month][(int(day), days_times[day][0])]

    with open(f"H2-H1_daily_means/{station_id}_H2-H1_daily_mean.pkl", "wb") as file:
        pickle.dump(data, file)

    return updated_data


def daily_H2_H1_mean_all_stations(station_list=station_list):
    """
    Compute and save daily means for tropopause thickness for all stations
    """
    for station in station_list:
        daily_H2_H1_mean(station)
        print(station)


def monthly_H2_H1_mean(station_id):
    """
    Compute and save monthly means for tropopause thickness for a single station
    """
    with open(f"H2-H1_daily_means/{station_id}_H2-H1_daily_mean.pkl", "rb") as file:
        data = pickle.load(file)

    updated_data = {}
    years = list(data.keys())
    for year in years:
        updated_data[year] = {}
        months = list(data[year].keys())
        for month in months:
            updated_data[year][month] = np.nanmean(list(data[year][month].values()))

    with open(f"H2-H1_monthly_means/{station_id}_H2-H1_monthly_mean.pkl", "wb") as file:
        pickle.dump(updated_data, file)

    return updated_data


def monthly_H2_H1_mean_all_stations(station_list=station_list):
    """
    Compute the monthly global mean for tropopause thickness

    Returns 
    x : list of years
    y : list of global means for the corresponding year
    """
    for station in station_list:
        monthly_H2_H1_mean(station)
        print(station)


def monthly_H2_H1_mean_global(station_list=station_list):
    """
    Compute monthly global means for tropopause thickness
    
    Returns 
    x : months
    y : global mean for the corresponding month
    """
    x = []
    y = []
    for year in range(1980, 2025):
        for month in [f'0{i}' for i in range(1,10)] + ['10','11','12']:
            current_month = []
            for station in station_list:
                with open(f"H2-H1_monthly_means/{station}_H2-H1_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month))
            x.append((int(year) - 1980) * 12 + (int(month) - 1))

    return x, y


def monthly_H2_H1_mean_global_single_month(month, station_list=station_list):
    """
    Compute yearly trend of global means for tropopause thickness for a single month
    
    Returns 
    x : years
    y : global mean for the year
    """
    x = []
    y = []
    for year in range(1980, 2025):
        current_year = []
        for station in station_list:
            with open(f"H2-H1_monthly_means/{station}_H2-H1_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        y.append(np.nanmean(current_year))
        x.append(year)

    return x, y


def monthly_H2_H1_anomaly_global_single_month(month, station_list=station_list):
    """
    Compute yearly trend of global anomalies for tropopause thickness for a single month

    Returns 
    x : years
    y : global anomaly for the year
    """
    x, y = monthly_H2_H1_mean_global_single_month(month, station_list)
    y = y - np.nanmean(y)
    return x, y


def monthly_H2_H1_anomaly_global(station_list=station_list):
    """
    Compute monthly global anomaly for tropopause thickness

    Returns 
    x : months
    y : global anomaly for the corresponding month
    """
    x, y = monthly_H2_H1_mean_global(station_list)
    y = y - np.nanmean(y)
    return x, y

                
def yearly_H2_H1_mean(station_id):
    """
    Compute and save yearly means for tropopause thickness for a single station
    """
    with open(f"H2-H1_monthly_means/{station_id}_H2-H1_monthly_mean.pkl", "rb") as file:
        data = pickle.load(file)

    updated_data = {}
    years = list(data.keys())
    for year in years:
        updated_data[year] = np.nanmean(list(data[year].values()))

    with open(f"H2-H1_yearly_means/{station_id}_H2-H1_yearly_mean.pkl", "wb") as file:
        pickle.dump(updated_data, file)
    print(updated_data)


def yearly_H2_H1_mean_all_stations(station_list=station_list):
    """
    Compute and save yearly means for tropopause thickness for all stations
    """
    for station in station_list:
        yearly_H2_H1_mean(station)
        print(station)


def yearly_H2_H1_mean_global(station_list=station_list):
    """
    Compute the global mean for tropopause thickness

    Returns 
    x : list of years
    y : list of global means for the corresponding year
    """
    x = []
    y = []
    for year in range(1980, 2025):
        current_y = []
        for station in station_list:
            with open(f"H2-H1_yearly_means/{station}_H2-H1_yearly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()):
                current_y.append(data[str(year)])
        y.append(np.nanmean(current_y))
        x.append(year)

    return x, y


def yearly_H2_H1_anomaly_global(station_list=station_list):
    """
    Compute the yearly global anomaly for tropopause thickness 

    Returns
    x : list of years
    y : list of global anomalies for the corresponding year
    """
    x, y = yearly_H2_H1_mean_global(station_list)
    y = y - np.nanmean(y)
    y = y*1000
    return x, y

x, y = yearly_H2_H1_anomaly_global()
m, b = np.polyfit(x, y, 1)  # 1 means linear fit
x_fit = np.linspace(min(x), max(x), 100)  # Smooth line
y_fit = m * x_fit + b
plt.plot(x_fit, y_fit, color='black', label=f'Fit: y = {m:.2f}x + {b:.2f}', lw=1)
plt.plot(x, y, color='blue', label='Data', lw=1)
plt.title(f'Yearly Northern Hemisphere Tropopause Mean Thickness Anomaly')
plt.xlabel('Year')
plt.ylabel('Thickness Anomaly (m)')
plt.ylim((-700,700))
plt.grid(True)
plt.savefig(f'Yearly Northern Hemisphere Tropopause Mean Thickness Anomaly', dpi=1200)

# import matplotlib
# inputs = [('01', 'January'), ('02', 'February'), ('03', 'March'), ('04', 'April'), ('05', 'May'),
#           ('06', 'June'), ('07', 'July'), ('07', 'August'), ('09', 'September'), ('10', 'October'),
#           ('11', 'November'), ('12', 'December')]
# for (i, j) in inputs:
#     matplotlib.rc('font', family='Times New Roman')
#     x, y = monthly_H2_H1_anomaly_global_single_month(f'{i}')

#     # Fit a linear model y = mx + b
#     m, b = np.polyfit(x, y*1000, 1)  # 1 means linear fit

#     # Generate fitted line
#     x_fit = np.linspace(min(x), max(x), 100)  # Smooth line
#     y_fit = m * x_fit + b

#     # Plot data points
#     plt.plot(x, y*1000, color='blue', label='Data', lw=1)

#     # Plot fitted line
#     plt.plot(x_fit, y_fit, color='black', label=f'Fit: y = {m:.2f}x + {b:.2f}', lw=1)

#     plt.xlabel('Year')
#     plt.ylabel('Thickness Anomaly (m)')
#     plt.title(f'Northern Hemisphere Tropopause Thickness Anomaly for {j}')
#     plt.ylim((-700,700))

#     plt.grid(True)


#     # Show plot
#     plt.savefig(f'Northern Hemisphere Tropopause Thickness Anomaly for {j}.png', dpi=1200)
#     plt.clf()




        
#################################### Compute and save files

# daily_H2_H1_mean_all_stations()
# monthly_H2_H1_mean_all_stations()
# yearly_H2_H1_mean_all_stations()


