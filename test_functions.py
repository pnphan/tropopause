from functions import *


# month = 5
# day = 15
# stepsize=50
# differences=[]
# for day in range(1,30):
#     gph = get_daywise_data(2024, month)[(day,0)]['gph']
#     temp = get_daywise_data(2024, month)[(day,0)]['temp']
#     gph, temp = interpolate_to_points(gph, temp, kind='linear', step=stepsize)
#     gph, tempgrad = temp_gradient_zangl(gph, temp)
#     tropopause1 = detect_tropopause_zangl(tempgrad, gph, temp)

#     gph = get_daywise_data(2024, month)[(day,0)]['gph']
#     temp = get_daywise_data(2024, month)[(day,0)]['temp']
#     gph, temp = interpolate_to_points(gph, temp, kind='linear', step=stepsize)
#     gph, lapse, temp = lapse_rate(gph, temp)
#     tropopause2 = detect_tropopause(gph, lapse)

#     if tropopause1 != None and tropopause2 != None:
#         differences.append(tropopause1*1000 - tropopause2)

# print(np.mean(np.sqrt(np.array(differences) ** 2)))


# print(get_years_and_months())


# def get_availability_flask(station_id='ZZM00099013'):
#     """
#     {
#     "years": [1980, 1981, ..., 2024],
#     "data": {
#         "1980": {
#         "01": {
#             "01": ["00Z", "12Z"],
#             "02": ["06Z"],
#             "15": ["00Z", "06Z", "12Z"],
#             ...
#         },
#         "02": {
#             "03": ["00Z", "12Z"],
#             ...
#         }
#         ...
#         },
#         ...
#         "2024": {
#         ...
#         }
#     }
#     }
#     """
#     output_dict = {}
#     years_months = get_years_and_months(filepath=station_id+'-data.txt')
#     years = list(years_months.keys())
#     output_dict['years'] = years
#     output_dict['data'] = {}
#     for year in years:
#         print(year)
#         output_dict['data'][year] = {}
#         for month in years_months[year]:
#             output_dict['data'][year][month] = {}
#             for day in get_days(year, month, filepath=station_id+'-data.txt'):
#                 output_dict['data'][year][month][day] = get_times(year, month, day, filepath=station_id+'-data.txt')

#     print(output_dict)

# get_availability_flask()


import matplotlib.pyplot as plt

def plot_igra_sounding(station_id, year, month, day, hour):
    """
    Reads an IGRA v2.2 station file named <station_id>-data.txt and plots:
      1) Temperature vs Geopotential Height
      2) Temperature vs Pressure
    for the specified sounding (year, month, day, hour).

    Parameters
    ----------
    station_id : str
        The station identifier (e.g., 'PFM00059981').
    year : int
        The 4-digit year of the sounding (e.g. 1992).
    month : int
        The 2-digit month of the sounding (1..12).
    day : int
        The 2-digit day of the month (1..31).
    hour : int
        The nominal observation hour in UTC (0..23). (99 = missing, which is skipped)

    Returns
    -------
    None
        Displays two matplotlib figures side by side (Temperature vs GPH, Temperature vs Pressure).
    """

    file_name = f"data/{station_id}-data.txt"

    # Data structure to hold the (pressure, height, temperature) for the chosen sounding
    data_rows = []
    
    # Internal state tracking
    found_header = False
    needed_numlev = 0  # We'll read this from the header

    try:
        with open(file_name, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")

                # Header record lines always start with '#'
                if line.startswith('#'):
                    # Parse header record according to IGRA fixed columns
                    # YEAR is columns 14-17 (line[13:17])
                    # MONTH is columns 19-20 (line[18:20])
                    # DAY is columns 22-23 (line[21:23])
                    # HOUR is columns 25-26 (line[24:26])
                    try:
                        hyear  = int(line[13:17])
                        hmonth = int(line[18:20])
                        hday   = int(line[21:23])
                        hhour  = int(line[24:26])
                    except ValueError:
                        # If parsing fails (malformed line), skip
                        continue

                    # Compare to our requested date/time
                    if (hyear == year and 
                        hmonth == month and 
                        hday == day and 
                        hhour == hour):
                        
                        # Mark that we've found the target header
                        found_header = True

                        # NUMLEV is columns 33-36 (line[32:36])
                        try:
                            needed_numlev = int(line[32:36])
                        except ValueError:
                            needed_numlev = 0
                        
                        # Reset data rows for new sounding
                        data_rows = []
                    else:
                        found_header = False

                else:
                    # If we are in the block of the desired header, parse data records
                    if found_header and needed_numlev > 0:
                        # IGRA data record layout (fixed columns):
                        # PRESS (columns 10-15)    => row[9:15]
                        # GPH   (columns 17-21)    => row[16:21]
                        # TEMP  (columns 23-27)    => row[22:27]
                        row = line
                        if len(row) < 27:
                            continue  # skip malformed lines

                        try:
                            press_str = row[9:15].strip()
                            gph_str   = row[16:21].strip()
                            temp_str  = row[22:27].strip()

                            press = int(press_str)
                            gph   = int(gph_str)
                            temp  = int(temp_str)
                        except ValueError:
                            # If any field is non-integer or missing
                            continue

                        data_rows.append((press, gph, temp))
                        
                        # If we've read all levels in the sounding, we can break early
                        if len(data_rows) >= needed_numlev:
                            found_header = False  # done reading this sounding
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return
    except Exception as e:
        print(f"Error reading/parsing {file_name}: {e}")
        return

    # Now we have the raw data rows for the specified sounding in data_rows.
    # Next step: convert raw IGRA values to physical units and filter out missing data.

    valid_press = []
    valid_gph   = []
    valid_temp  = []

    for (p, z, t) in data_rows:
        # Pressure: -9999 / -8888 => missing or removed
        # If p is e.g. 100000 => 1000.0 hPa
        if p not in (-9999, -8888) and p >= 0:
            p_hpa = p / 100.0  # Convert from Pa to hPa
        else:
            p_hpa = None

        # Geopotential Height: -9999 / -8888 => missing or removed
        # Typically in meters
        if z not in (-9999, -8888) and z >= 0:
            z_m = z
        else:
            z_m = None

        # Temperature: -9999 / -8888 => missing or removed
        # Values are tenths of a degree C => dividing by 10
        if t not in (-9999, -8888):
            temp_c = t / 10.0
        else:
            temp_c = None

        # Store
        valid_press.append(p_hpa)
        valid_gph.append(z_m)
        valid_temp.append(temp_c)

    # We'll build separate lists for each plot (filter out any None).
    # 1) Temperature vs Geopotential Height
    # 2) Temperature vs Pressure
    # Note: We might have some levels that have a valid pressure but missing GPH, etc.

    # For Temperature vs GPH:
    temp_gph = []
    gph_vals = []
    for t, z in zip(valid_temp, valid_gph):
        if t is not None and z is not None:
            temp_gph.append(t)
            gph_vals.append(z)

    # For Temperature vs Pressure:
    temp_pres = []
    pres_vals = []
    for t, p_ in zip(valid_temp, valid_press):
        if t is not None and p_ is not None:
            temp_pres.append(t)
            pres_vals.append(p_)

    if not temp_gph and not temp_pres:
        print("No valid data found for the specified sounding.")
        return

    # Create the subplots
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(5, 5))

    # Plot 1: Temperature vs Geopotential Height
    if temp_gph and gph_vals:
        ax1.plot(temp_gph, gph_vals, marker='o', linestyle='-')
        ax1.set_xlabel("Temperature (°C)")
        ax1.set_ylabel("Geopotential Height (m)")
        ax1.set_title(f"{station_id}\n{year}-{month:02d}-{day:02d} {hour:02d}\nTemp vs Height")
    else:
        ax1.text(0.5, 0.5, "No T-GPH data", ha='center', va='center')
        ax1.set_title("No T vs GPH data")
        ax1.set_xlabel("Temperature (°C)")
        ax1.set_ylabel("Geopotential Height (m)")

    # Plot 2: Temperature vs Pressure
    if temp_pres and pres_vals:
        ax2.plot(temp_pres, pres_vals, marker='o', linestyle='-')
        ax2.set_xlabel("Temperature (°C)")
        ax2.set_ylabel("Pressure (hPa)")
        ax2.set_title(f"{station_id}\n{year}-{month:02d}-{day:02d} {hour:02d}\nTemp vs Pressure")
        # Invert the Pressure axis so that higher altitude is higher on the plot
        ax2.invert_yaxis()
    else:
        ax2.text(0.5, 0.5, "No T-P data", ha='center', va='center')
        ax2.set_title("No T vs Pressure data")
        ax2.set_xlabel("Temperature (°C)")
        ax2.set_ylabel("Pressure (hPa)")

    plt.tight_layout()
    plt.show()


# Example usage (uncomment to try):
plot_igra_sounding("PFM00059981", 1992, 5, 23, 12)
