from functions import get_station_list
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

station_list = ['GQM00091212', 'RQM00078526', 'USM00091285', 'USM00091165', 'HKM00045004', 'JAM00047991', 'USM00072201', 'JAM00047945',
 'USM00072250', 'JAM00047971', 'USM00072210', 'USM00072261', 'USM00072240', 'CHM00058457', 'JAM00047827', 'USM00072364',
 'CHM00058238', 'USM00072265', 'ISM00040179', 'USM00072274', 'USM00072235', 'USM00072293', 'USM00072208', 'CHM00057127', 
 'JAM00047678', 'JAM00047778', 'JAM00047807', 'CHM00058150', 'CHM00058027', 'JAM00047681', 'USM00072340', 'KSM00047158', 
 'USM00072363', 'CHM00053915', 'JAM00047646', 'CHM00054857', 'USM00072317', 'USM00072327', 'CHM00053845', 'CHM00051828', 
 'JAM00047600', 'CHM00053772', 'USM00072493', 'USM00072451', 'USM00072402', 'CHM00053614', 'CHM00052681', 'USM00072403', 
 'CHM00051777', 'USM00072456', 'CHM00051709', 'JAM00047582', 'CHM00054511', 'CHM00052418', 'SPM00008221', 'USM00072520', 
 'ITM00016320', 'JAM00047580', 'USM00072572', 'USM00072562', 'USM00072558', 'USM00074494', 'CHM00054218', 'USM00072597', 
 'USM00072518', 'CHM00054292', 'USM00072528', 'JAM00047412', 'USM00072681', 'CHM00054135', 'CHM00051463', 'USM00074389', 
 'CAM00071600', 'CHM00054102', 'USM00072662', 'USM00072645', 'FRM00007510', 'USM00072694', 'JAM00047401', 'ITM00016080', 
 'CHM00050953', 'CAM00071722', 'USM00072764', 'SZM00006610', 'USM00072712', 'HUM00012843', 'USM00072776', 'USM00072786', 
 'CHM00051076', 'USM00072797', 'USM00072768', 'GMM00010868', 'AUM00011035', 'USM00072747', 'CAM00071815', 'GMM00010739', 
 'LOM00011952', 'CHM00050557', 'EZM00011520', 'CAM00071811', 'GMM00010548', 'CAM00071109', 'CAM00071836', 'GMM00010410', 
 'EIM00003953', 'GMM00010393', 'PLM00012374', 'RSM00032540', 'CAM00071816', 'CAM00071119', 'CAM00071867', 'GMM00010184', 
 'GMM00010035', 'USM00070398', 'USM00070316', 'USM00070308', 'USM00070350', 'CAM00071906', 'CAM00071907', 'USM00070326', 
 'CAM00071913', 'USM00070361', 'CAM00071934', 'UKM00003005', 'USM00070219', 'FIM00002963', 'USM00070273', 'GLM00004270', 
 'USM00070231', 'NOM00001241', 'CAM00071909', 'ICM00004018', 'USM00070200', 'RSM00022543', 'USM00070261', 'CAM00071043', 
 'CAM00071957', 'CAM00071081', 'CAM00071925', 'JNM00001001', 'USM00070026', 'SVM00001028', 'CAM00071924', 'CAM00071917', 
 'SPM00008001', 'UKM00003808', 'NOM00001415', 'FIM00002836', 'GLM00004220']

def get_daywise_data(year, month, station_id='PFM00059981'):
    """
    Reads the IGRA file and extracts geopotential height, pressure, and temperature
    for each available day and time for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A dictionary where each key is a day, and the value is a list of tuples
             (time, geopotential height list, pressure list, temperature list).
    """
    daywise_data = {}

    filepath = f"data/{station_id}-data.txt"

    with open(filepath, 'r') as file:
        current_day = None
        current_time = None

        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    # Extract year, month, day, and time
                    line_year = int(line[13:17])  # Year (columns 14-17)
                    line_month = int(line[18:20])  # Month (columns 19-20)
                    line_day = int(line[21:23])  # Day (columns 22-23)
                    line_hour = int(line[24:26])  # Hour (columns 25-26)
                    if line_year == year and line_month == month:
                        current_day = line_day
                        current_time = line_hour
                        daywise_data[(current_day, current_time)] = {'gph':[], 'pressure':[], 'temp':[]}
                except ValueError:
                    continue

            elif line_year == year and line_month == month:
                gph = int(line[16:21])  # GPH (columns 17-21)
                #pressure = int(line[9:15])  # Pressure (columns 10-15)
                temp = int(line[22:27])  # Temperature (columns 23-27)
                if gph != -9999 and gph != -8888 and temp != -9999: #and pressure != -9999 and
                    daywise_data[(current_day, current_time)]['gph'].append(gph)
                    #daywise_data[(current_day, current_time)]['pressure'].append(pressure / 100) # / 100 convert to hPa
                    daywise_data[(current_day, current_time)]['temp'].append(temp / 10) # / 10 units are in 10ths of a degree celsius
    
    return daywise_data


def get_availability(station_id):
# def parse_igra_station(station_id):
    """
    Parse an IGRA v2.2 station file of the form <station_id>-data.txt
    and return the available years, months, days, and times in the format:

    {
        "years": [...],
        "data": {
            "YYYY": {
                "MM": {
                    "DD": ["HHZ", ...],
                    ...
                },
                ...
            },
            ...
        }
    }
    """

    # Prepare the output structure
    data_dict = {
        "years": [],
        "data": {}
    }

    # The IGRA v2.2 header record has fixed-column positions:
    # HEADREC columns 1-1   (index 0)
    # ID      columns 2-12  (index 1..11)
    # YEAR    columns 14-17 (index 13..16)
    # MONTH   columns 19-20 (index 18..19)
    # DAY     columns 22-23 (index 21..22)
    # HOUR    columns 25-26 (index 24..25)
    # We will ignore the rest for collecting date/time info.

    file_name = f"data/{station_id}-data.txt"

    try:
        with open(file_name, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                # Only parse lines that start with '#'
                if not line.startswith("#"):
                    continue

                # Extract the columns defined by IGRA
                # Make sure each slice is large enough and handle any short lines gracefully.
                year_str  = line[13:17].strip()
                month_str = line[18:20].strip()
                day_str   = line[21:23].strip()
                hour_str  = line[24:26].strip()

                # If any of these are missing or not digits, skip.
                # (Sometimes lines can be malformed or short.)
                if (not year_str.isdigit() or 
                    not month_str.isdigit() or 
                    not day_str.isdigit() or 
                    not hour_str.isdigit()):
                    continue

                year  = int(year_str)
                month = int(month_str)
                day   = int(day_str)
                hour  = int(hour_str)

                # Hour can be 0..23 or 99 = missing
                if hour == 99:
                    # Skip missing hour
                    continue

                # Convert to strings with zero padding if needed
                year_s  = str(year)
                month_s = f"{month:02d}"
                day_s   = f"{day:02d}"
                hour_s  = f"{hour:02d}" #Z

                # Insert into data structure
                if year_s not in data_dict["data"]:
                    data_dict["data"][year_s] = {}
                if month_s not in data_dict["data"][year_s]:
                    data_dict["data"][year_s][month_s] = {}
                if day_s not in data_dict["data"][year_s][month_s]:
                    data_dict["data"][year_s][month_s][day_s] = []

                # Append the hour if not already present
                if hour_s not in data_dict["data"][year_s][month_s][day_s]:
                    data_dict["data"][year_s][month_s][day_s].append(hour_s)

        # Now fill the "years" list. Sort them numerically.
        # The keys in data_dict["data"] are strings, so convert to int for sorting.
        sorted_years = sorted(int(y) for y in data_dict["data"].keys())
        data_dict["years"] = sorted_years

        # Convert them back to strings if you want the top-level "years" 
        # to be integers or strings. The example shows them as integers, 
        # so we can just store them as int.
        # data_dict["years"] = sorted_years

        return data_dict

    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return data_dict
    except Exception as e:
        print(f"Error reading/parsing {file_name}: {e}")
        return data_dict
    

def temp_gradient_zangl(gph, temp):
    gph = np.array(gph) / 1000
    temp = np.array(temp)
    tempgrad = (temp[1:] - temp[:-1]) / (gph[1:] - gph[:-1])
    return gph, tempgrad

def detect_tropopause_zangl(tempgrad, gph, temp):
    zhalf = (gph[1:] + gph[:-1])/2
    zTP = []
    possible_heights = []
    tempTP = []
    for i in range(len(tempgrad)):
        if tempgrad[i] >= -2:
            if tempgrad[i-1] < -2:
                ztp = gph[i-1] + (-2 - tempgrad[i])/(tempgrad[i] - tempgrad[i-1])*(zhalf[i] - zhalf[i-1])
                zTP.append(ztp)
                possible_heights.append(i)
                if ztp > gph[i]:
                    tempTP.append(temp[i] + tempgrad[i-1] * (ztp - gph[i]))
                else:
                    tempTP.append(temp[i] + tempgrad[i] * (ztp - gph[i]))
    
    for j in range(len(possible_heights)):
        checknext = False
        if possible_heights[j] < len(gph) - 10:
            for k in range(possible_heights[j], possible_heights[j]+10):
                if (temp[k] - tempTP[j]) / (gph[k] - zTP[j]) <= -2:
                    checknext=True
                    break
        else:
            for k in range(possible_heights[j], len(gph)):
                if (temp[k] - tempTP[j]) / (gph[k] - zTP[j]) <= -2:
                    checknext=True
                    break
        if checknext == False:
            return zTP[j]
        else:
            continue

def plot_values(x_values, y_values):
    """
    Plots a graph using the given x and y values.

    Parameters:
        x_values (list): List of x values.
        y_values (list): List of y values.
    """
    if len(x_values) != len(y_values):
        raise ValueError("x_values and y_values must have the same length.")
    
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, marker='o', linestyle='-', color='b')
    plt.xlabel('X Values')
    plt.ylabel('Y Values')
    plt.title('Plot of X vs. Y')
    plt.grid(True)
    plt.show()

def interpolate_gph_temp(gph, temp, kind='linear'):
    interpolation_func = interp1d(gph, temp, kind=kind, fill_value="extrapolate")
    gph = np.arange(5000, 25200, 200)
    temp = interpolation_func(gph)
    return gph, temp

def detect_tropopause(gph, temp):
    gph = gph/1000
    tropopause_found = False
    start_level = 0

    while tropopause_found == False:
        gph_half = (gph[:-1] + gph[1:])/2
        temp_gradient = (temp[1:] - temp[:-1]) / (gph[1:] - gph[:-1])
        i = -9999
        for j in range(start_level, len(temp_gradient)):
            if temp_gradient[j] >= -2:
                i = j
                break
        if i == -9999:
            return -9999
        zTP = gph_half[i-1] + ((0.2)/(temp_gradient[i] - temp_gradient[i-1])) * (-2 - temp_gradient[i])
        if zTP > gph[i]:
            tTP = temp[i] + temp_gradient[i-1] * (zTP - gph[i])
        else:
            tTP = temp[i] + temp_gradient[i] * (zTP - gph[i])
        
        restart = False
        for j in range(len(gph)):
            if 0 < gph[j] - zTP <= 2 and (temp[j] - tTP)/(gph[j] - zTP) <= -2:
                restart = True
                break
        
        if restart == True:
            if start_level == len(temp_gradient):
                return -9999
            else:
                start_level = i+1
        else:
            return zTP
        

def compute_all_tropopause(stationid):
    availability = get_availability(stationid)['data']
    years = list(availability.keys())
    for j in range(len(years)):
        if int(years[j]) >= 1980:
            i=j
            break
    years = years[i:]
    tropopause_data = {}
    for year in years:
        print(year)
        tropopause_data[year] = {}
        months = list(availability[year].keys())
        for month in months:
            tropopause_data[year][month] = {}
            data_month = get_daywise_data(int(year), int(month), 'JAM00047991')
            days_times = list(data_month.keys())
            for day_time in days_times:
                gph = data_month[day_time]['gph']
                temp = data_month[day_time]['temp']
                if len(gph) == 0 or len(temp) == 0:
                    continue
                gph, temp = interpolate_gph_temp(gph, temp, 'linear')
                tropopause = detect_tropopause(gph, temp)
                if tropopause != -9999 and abs(tropopause) <= 20:
                    tropopause_data[year][month][day_time] = tropopause

    return tropopause_data
    

print(compute_all_tropopause('JAM00047991'))


# data = get_daywise_data(2016, 6, 'JAM00047991')[(2,0)]
# gph = data['gph']
# temp = data['temp']
# gph, temp = interpolate_gph_temp(gph, temp, 'cubic')
# print(gph)
# print(temp)
# plot_values(temp, gph)
# print(detect_tropopause(gph, temp))


############ STATION LIST ###################

# ['GQM00091212', 'RQM00078526', 'USM00091285', 'USM00091165', 'HKM00045004', 'JAM00047991', 'USM00072201', 'JAM00047945',
#  'USM00072250', 'JAM00047971', 'USM00072210', 'USM00072261', 'USM00072240', 'CHM00058457', 'JAM00047827', 'USM00072364',
#  'CHM00058238', 'USM00072265', 'ISM00040179', 'USM00072274', 'USM00072235', 'USM00072293', 'USM00072208', 'CHM00057127', 
#  'JAM00047678', 'JAM00047778', 'JAM00047807', 'CHM00058150', 'CHM00058027', 'JAM00047681', 'USM00072340', 'KSM00047158', 
#  'USM00072363', 'CHM00053915', 'JAM00047646', 'CHM00054857', 'USM00072317', 'USM00072327', 'CHM00053845', 'CHM00051828', 
#  'JAM00047600', 'CHM00053772', 'USM00072493', 'USM00072451', 'USM00072402', 'CHM00053614', 'CHM00052681', 'USM00072403', 
#  'CHM00051777', 'USM00072456', 'CHM00051709', 'JAM00047582', 'CHM00054511', 'CHM00052418', 'SPM00008221', 'USM00072520', 
#  'ITM00016320', 'JAM00047580', 'USM00072572', 'USM00072562', 'USM00072558', 'USM00074494', 'CHM00054218', 'USM00072597', 
#  'USM00072518', 'CHM00054292', 'USM00072528', 'JAM00047412', 'USM00072681', 'CHM00054135', 'CHM00051463', 'USM00074389', 
#  'CAM00071600', 'CHM00054102', 'USM00072662', 'USM00072645', 'FRM00007510', 'USM00072694', 'JAM00047401', 'ITM00016080', 
#  'CHM00050953', 'CAM00071722', 'USM00072764', 'SZM00006610', 'USM00072712', 'HUM00012843', 'USM00072776', 'USM00072786', 
#  'CHM00051076', 'USM00072797', 'USM00072768', 'GMM00010868', 'AUM00011035', 'USM00072747', 'CAM00071815', 'GMM00010739', 
#  'LOM00011952', 'CHM00050557', 'EZM00011520', 'CAM00071811', 'GMM00010548', 'CAM00071109', 'CAM00071836', 'GMM00010410', 
#  'EIM00003953', 'GMM00010393', 'PLM00012374', 'RSM00032540', 'CAM00071816', 'CAM00071119', 'CAM00071867', 'GMM00010184', 
#  'GMM00010035', 'USM00070398', 'USM00070316', 'USM00070308', 'USM00070350', 'CAM00071906', 'CAM00071907', 'USM00070326', 
#  'CAM00071913', 'USM00070361', 'CAM00071934', 'UKM00003005', 'USM00070219', 'FIM00002963', 'USM00070273', 'GLM00004270', 
#  'USM00070231', 'NOM00001241', 'CAM00071909', 'ICM00004018', 'USM00070200', 'RSM00022543', 'USM00070261', 'CAM00071043', 
#  'CAM00071957', 'CAM00071081', 'CAM00071925', 'JNM00001001', 'USM00070026', 'SVM00001028', 'CAM00071924', 'CAM00071917', 
#  'SPM00008001', 'UKM00003808', 'NOM00001415', 'FIM00002836', 'GLM00004220']

# meng_stations = [
#     '91212', '78526', '91285', '91165', '45004', '47991', '72201', '47945', '72250', '47971', 
#     '72210', '72261', '72240', '58457', '47827', '72364', '58238', '72265', '40179', '72274', 
#     '72235', '72293', '72208', '57127', '47678', '47778', '47807', '58150', '58027', '47681', 
#     '72340', '47158', '72363', '53915', '47646', '54857', '72317', '72327', '53845', '51828', 
#     '47600', '53772', '72493', '72451', '72402', '53614', '52681', '72403', '51777', '72456', 
#     '51709', '47582', '54511', '52418', '8221', '72520', '16320', '47580', '72572', '72562', 
#     '72558', '74494', '54218', '72597', '72518', '54292', '72528', '47412', '8001', '72681', 
#     '54135', '51463', '74389', '71600', '54102', '72662', '72645', '7510', '72694', '47401', '16080', 
#     '50953', '71722', '72764', '6610', '72712', '12843', '72776', '72786', '51076', '72797', '72768', 
#     '10868', '11035', '72747', '71815', '10739', '11952', '50557', '11520', '3808', '71811', '10548', '71109', 
#     '71836', '10410', '3953', '10393', '12374', '32540', '71816', '71119', '71867', '10184', '10035', '70398', 
#     '70316', '70308', '70350', '71906', '71907', '70326', '71913', '1415', '70361', '71934', '3005', '70219', 
#     '2963', '70273', '4270', '70231', '1241', '71909', '4018', '70200', '22543', '70261', '71043', '2836', '71957', 
#     '4220', '71081', '71925', '1001', '70026', '1028', '71924', '71917'
# ]

# duplicates = [] # SPM00008001, UKM00003808, NOM00001415, FIM00002836, GLM00004220 