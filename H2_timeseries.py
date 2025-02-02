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




x, y = yearly_H2_anomaly_global()
plot_values(x, y)



        


