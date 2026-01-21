import json
import csv
import os
from datetime import datetime
from skyfield.api import load, Loader
from skyfield import almanac

def load_config(path="config.json"):
    with open(path, "r") as f:
        return json.load(f)

def ensure_de440(filename):
    loader = Loader('.') 
    if not os.path.exists(filename):
        print(f"Downloading de440.bsp...")
        planets = loader(filename)
    else:
        planets = loader(filename)
    return planets

def calculate_seasons():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    ts = load.timescale()
    
    # Parse UTC times from config
    t_start = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["start"].replace('Z', '+00:00'))
    )
    t_end = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["end"].replace('Z', '+00:00'))
    )
    
    print(f"Calculating Seasons...")

    # skyfield.almanac.seasons identifies when the Sun crosses 0, 90, 180, and 270 degrees
    # of geocentric apparent ecliptic longitude.
    f = almanac.seasons(eph)
    
    # find_discrete is the precise method for state-change events
    t_times, y_events = almanac.find_discrete(t_start, t_end, f)

    season_names = {
        0: "March Equinox",
        1: "June Solstice",
        2: "September Equinox",
        3: "December Solstice"
    }

    # Output to CSV
    output_filename = "seasons.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        
        for t, event_code in zip(t_times, y_events):
            writer.writerow([
                t.utc_iso(),
                season_names[event_code],
                ""  # Third column is empty
            ])
            
    print(f"Success. Wrote {len(t_times)} events to {output_filename}.")

if __name__ == "__main__":
    calculate_seasons()
