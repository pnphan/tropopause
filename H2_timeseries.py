from meng2021replication import *

stationid='USM00091285'
with open(f"tropopauses/{stationid}.pkl", "rb") as f:
    data = pickle.load(f)


def daily_H2_mean(station_id):
    """
    Compute the daily means of second tropopause height (if available) for a single station

    Returns a dict

    updated_data : {year : 
                        {month : 
                            {day : H2} 
                        }
                    }
    """
    with open(f"tropopauses/{station_id}.pkl", "rb") as file:
        data = pickle.load(file)
    data = data['second']
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

    with open(f"H2_daily_means/{station_id}_H2_daily_mean.pkl", "wb") as file:
        pickle.dump(data, file)

    return updated_data


def daily_H2_mean_all_stations(station_list=station_list):
    """
    Compute and save daily means for second tropopause height for all stations
    """
    for station in station_list:
        daily_H2_mean(station)
        print(station)


def monthly_H2_mean(station_id):
    """
    Compute and save monthly means for second tropopause height for a single station
    """
    with open(f"H2_daily_means/{station_id}_H2_daily_mean.pkl", "rb") as file:
        data = pickle.load(file)

    updated_data = {}
    years = list(data.keys())
    for year in years:
        updated_data[year] = {}
        months = list(data[year].keys())
        for month in months:
            updated_data[year][month] = np.nanmean(list(data[year][month].values()))

    with open(f"H2_monthly_means/{station_id}_H2_monthly_mean.pkl", "wb") as file:
        pickle.dump(updated_data, file)

    return updated_data


def monthly_H2_mean_all_stations(station_list=station_list):
    """
    Compute and save monthly means for second tropopause height for all stations
    """
    for station in station_list:
        monthly_H2_mean(station)
        print(station)


def monthly_H2_mean_global(station_list=station_list):
    """
    Compute monthly global means for second tropopause height
    
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
                with open(f"H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month))
            x.append((int(year) - 1980) * 12 + (int(month) - 1))

    return x, y


def yearly_H2_mean(station_id):
    """
    Compute and save yearly means for second tropopause height for a single station
    """
    with open(f"H2_monthly_means/{station_id}_H2_monthly_mean.pkl", "rb") as file:
        data = pickle.load(file)

    updated_data = {}
    years = list(data.keys())
    for year in years:
        updated_data[year] = np.nanmean(list(data[year].values()))

    with open(f"H2_yearly_means/{station_id}_H2_yearly_mean.pkl", "wb") as file:
        pickle.dump(updated_data, file)
    print(updated_data)


def yearly_H2_mean_all_stations(station_list=station_list):
    """
    Compute and save yearly means for second tropopause height for all stations
    """
    for station in station_list:
        yearly_H2_mean(station)
        print(station)


def yearly_H2_mean_global(station_list=station_list):
    """
    Compute the global mean for second tropopause height

    Returns 
    x : list of years
    y : list of global means for the corresponding year
    """
    x = []
    y = []
    for year in range(1980, 2025):
        current_y = []
        for station in station_list:
            with open(f"H2_yearly_means/{station}_H2_yearly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()):
                current_y.append(data[str(year)])
        y.append(np.nanmean(current_y))
        x.append(year)

    return x, y


def yearly_H2_anomaly_global(station_list=station_list):
    """
    Compute the yearly global anomaly for second tropopause height 

    Returns
    x : list of years
    y : list of global anomalies for the corresponding year
    """
    x, y = yearly_H2_mean_global()
    y = y - np.nanmean(y)
    y = y*1000
    return x, y


def monthly_H2_mean_global_single_month(month, station_list=station_list):
    """
    Compute yearly trend of global means for first tropopause height for a single month
    
    Returns 
    x : years
    y : global mean for the year
    """
    x = []
    y = []
    for year in range(1980, 2025):
        current_year = []
        for station in station_list:
            with open(f"H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                data = pickle.load(file)
            if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                current_year.append(data[str(year)][month])
        y.append(np.nanmean(current_year))
        x.append(year)

    return x, y




# station_id = 'USM00091285'
# with open(f"H2_monthly_means/{station_id}_H2_monthly_mean.pkl", "rb") as file:
#     data = pickle.load(file)


# print(data.keys())
def H1_monthly_anomalies(station_list=station_list):
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    monthly_means = {}
    for month in months:
        current_month = []
        for year in range(1980,2024):
            current_year = []
            for station in station_list:
                with open(f"H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_year.append(data[str(year)][month])
            current_month.append(np.nanmean(current_year))
        monthly_means[month] = np.nanmean(current_month)


    y = []
    x = []
    counter = 0
    for year in range(1980,2024):
        for month in months:
            current_month = []
            for station in station_list:
                with open(f"H1_monthly_means/{station}_H1_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month) - monthly_means[month])
            counter+=1 
            x.append(counter)
    years = range(1980,2025)

    x_ticks = [i * 12 for i in range(len(years))]
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    

    # A = np.vstack([x, np.ones(len(x))]).T
    # m, c = np.linalg.lstsq(A, y)[0]
    # plt.plot(x, y)
    # plt.plot(x, m*np.array(x) + c, 'r')
    # plt.xticks(x_ticks, x_tick_labels, rotation=45)

    # plt.xlabel('Year')
    # plt.ylabel('Anomaly (km)')
    # #plt.grid(True)
    # plt.show()
    return x, y

def H2_monthly_anomalies(station_list=station_list):
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    monthly_means = {}
    for month in months:
        current_month = []
        for year in range(1980,2024):
            current_year = []
            for station in station_list:
                with open(f"H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_year.append(data[str(year)][month])
            current_month.append(np.nanmean(current_year))
        monthly_means[month] = np.nanmean(current_month)


    y = []
    x = []
    counter = 0
    for year in range(1980,2024):
        for month in months:
            current_month = []
            for station in station_list:
                with open(f"H2_monthly_means/{station}_H2_monthly_mean.pkl", "rb") as file:
                    data = pickle.load(file)
                if str(year) in list(data.keys()) and month in list(data[str(year)].keys()):
                    current_month.append(data[str(year)][month])
            y.append(np.nanmean(current_month) - monthly_means[month])
            counter+=1 
            x.append(counter)
    years = range(1980,2025)

    x_ticks = [i * 12 for i in range(len(years))]
    x_tick_labels = [str(year) if year%5==0 else '' for year in years]
    

    # A = np.vstack([x, np.ones(len(x))]).T
    # m, c = np.linalg.lstsq(A, y)[0]
    # plt.plot(x, y)
    # plt.plot(x, m*np.array(x) + c, 'r')
    # plt.xticks(x_ticks, x_tick_labels, rotation=45)

    # plt.xlabel('Year')
    # plt.ylabel('Anomaly (km)')
    # #plt.grid(True)
    # plt.show()
    return x, y



years = range(1980,2025)
x, y = H1_monthly_anomalies()
A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A, y)[0]
plt.plot(x, y, label='H1')
plt.plot(x, m*np.array(x) + c, 'r')

x1, y1 = H2_monthly_anomalies()
A1 = np.vstack([x1, np.ones(len(x1))]).T
m1, c1 = np.linalg.lstsq(A1, y1)[0]
plt.plot(x1, y1, label = 'H2')
plt.plot(x1, m1*np.array(x1) + c1, 'r')


x_ticks = [i * 12 for i in range(len(years))]
x_tick_labels = [str(year) if year%5==0 else '' for year in years]
plt.xticks(x_ticks, x_tick_labels, rotation=45)

plt.xlabel('Year')
plt.ylabel('Anomaly (km)')
#plt.grid(True)
plt.legend()
plt.show()

    

        

        

        


                


        


