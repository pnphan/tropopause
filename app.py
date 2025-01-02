from flask import Flask, jsonify

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

app = Flask(__name__, static_folder='static')

@app.route("/stations")
def get_stations():
    data = get_station_data()
    return jsonify(data)

@app.route("/")
def index():
    return app.send_static_file("index copy.html")

if __name__ == "__main__":
    app.run(debug=True)
