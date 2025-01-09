from flask import Flask, jsonify, send_file
import matplotlib.pyplot as plt
import io

def get_station_data(filepath='igra2-station-list.txt'):
    # Read station data and extract latitude and longitude
    with open(filepath, "r") as f:
        stations = f.read().splitlines()[:-101]  # Remove trailing meta-info lines
    
    stationdata = []
    for station in stations:
        try:
            # Extract latitude and longitude
            latitude = float(station[12:20].strip())
            longitude = float(station[21:30].strip())
            
            # Ignore stations with invalid lat/lon (set to missing values in the file)
            if latitude == -98.8888 or longitude == -998.8888:
                continue

            stationdata.append({'id': station[:11], 'lat': latitude, 'lon': longitude})
        except ValueError:
            print(f"Skipping invalid line: {station}")

    return stationdata


def get_years_and_months(filepath='PFM00059981-data.txt'):
    """
    Extracts available years and their corresponding months from the IGRA file.

    :param filepath: Path to the IGRA sounding data file.
    :return: A dictionary where keys are years and values are lists of available months.
    """
    year_month_data = {}

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    # Extract year and month
                    year = int(line[13:17])  # Year (columns 14-17)
                    month = int(line[18:20])  # Month (columns 19-20)

                    # Add the month to the corresponding year
                    if year not in year_month_data:
                        year_month_data[year] = set()  # Use a set to avoid duplicate months
                    year_month_data[year].add(month)
                except ValueError:
                    continue

    # Convert sets to sorted lists for consistency
    for year in year_month_data:
        year_month_data[year] = sorted(year_month_data[year])

    return year_month_data


def get_available_days_and_times(year, month, filepath='PFM00059981-data.txt'):
    """
    Reads the IGRA file and extracts days and their available times for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A dictionary where keys are days and values are lists of available times.
    """
    available_data = {}

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    line_year = int(line[13:17])  # Extract year (columns 14-17, zero-indexed)
                    line_month = int(line[18:20])  # Extract month (columns 19-20, zero-indexed)
                    line_day = int(line[21:23])  # Extract day (columns 22-23, zero-indexed)
                    line_hour = int(line[24:26])  # Extract hour (columns 25-26, zero-indexed)

                    if line_year == year and line_month == month:
                        if line_day not in available_data:
                            available_data[line_day] = []
                        available_data[line_day].append(line_hour)
                except ValueError:
                    continue

    # Sort the times for each day
    for day in available_data:
        available_data[day] = sorted(available_data[day])

    return available_data


def get_years(filepath='PFM00059981-data.txt'):
    return list(get_years_and_months(filepath).keys())
def get_months(year, filepath='PFM00059981-data.txt'):
    return get_years_and_months(filepath)[year]
def get_days(year, month, filepath='PFM00059981-data.txt'):
    """
    Reads the IGRA file and extracts days where data is available for a specific year and month.

    :param filepath: Path to the IGRA sounding data file.
    :param year: The year to filter the data.
    :param month: The month to filter the data.
    :return: A list of days with data available.
    """
    available_days = set()

    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Header line
                try:
                    line_year = int(line[13:17])  # Extract year (columns 14-17, zero-indexed)
                    line_month = int(line[18:20])  # Extract month (columns 19-20, zero-indexed)
                    line_day = int(line[21:23])  # Extract day (columns 22-23, zero-indexed)

                    if line_year == year and line_month == month:
                        available_days.add(line_day)
                except ValueError:
                    continue

    return sorted(available_days)
def get_times(year, month, day, filepath='PFM00059981-data.txt'):
    return get_available_days_and_times(year, month, filepath)[day]


def get_availability_flask(station_id):
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


app = Flask(__name__, static_folder='static')

@app.route("/stations")
def get_stations():
    data = get_station_data()
    return jsonify(data)

@app.route('/availability/<station_id>')
def get_availability(station_id):
    data = get_availability_flask(station_id)
    return jsonify(data)

@app.route('/plot-graph/<input>')
def plot_igra_sounding(input):
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
    parts = input.split('split')

    # Extract each part into separate variables
    station_id = parts[0]
    year = int(parts[1])
    month = int(parts[2])
    day = int(parts[3])
    hour = int(parts[4])

    print(station_id)
    print(year)
    print(month)
    print(day)
    print(hour)


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
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))

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
    # plt.show()

    buf = io.BytesIO()
    fig.savefig(buf, format='png')  # can also use 'jpg' etc.
    buf.seek(0)  # Important: go to the start of the buffer

    # 3. Optionally close the figure to free memory
    plt.close(fig)

    # 4. Return the buffer as a file-like object with 'image/png' MIME type
    return send_file(buf, mimetype='image/png')



@app.route("/")
def index():
    return app.send_static_file("index_calendar.html") #copy

if __name__ == "__main__":
    app.run(debug=True)
