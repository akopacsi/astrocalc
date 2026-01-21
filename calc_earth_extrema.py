import json
import csv
import os
from datetime import datetime
from skyfield.api import load, Loader
from skyfield import searchlib
import numpy as np

def load_config(path="config.json"):
    with open(path, "r") as f:
        return json.load(f)

def ensure_de440(filename):
    loader = Loader('.') 
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        planets = loader(filename)
    else:
        planets = loader(filename)
    return planets

def calculate_earth_extrema():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    
    ts = load.timescale()
    t_start = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["start"].replace('Z', '+00:00'))
    )
    t_end = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["end"].replace('Z', '+00:00'))
    )
    
    print(f"Calculating Earth Perihelion/Aphelion (Velocity Flip Method)...")

    # --- THE CORRECT ELEGANT LOGIC ---
    # We define a function that returns an INTEGER (0 or 1).
    # 1 = Moving Away (Distance Increasing)
    # 0 = Moving Closer (Distance Decreasing)
    
    def is_moving_away(t):
        p = (earth - sun).at(t)
        
        # Dot product of position and velocity
        # pos · vel > 0 means distance is increasing
        pos = p.position.au
        vel = p.velocity.au_per_d
        
        # We use numpy sum for efficient dot product across arrays
        dot_product = np.sum(pos * vel, axis=0)
        
        # Return 1 if increasing, 0 if decreasing
        return dot_product > 0

    # Set the step size (period hint).
    # searchlib will scan this function.
    is_moving_away.rough_period = 365.0

    # find_discrete finds exactly when the return value changes
    t_times, values = searchlib.find_discrete(t_start, t_end, is_moving_away)
    
    events = []

    for t, state_after_flip in zip(t_times, values):
        # Calculate distance for the report
        p = (earth - sun).at(t)
        dist_km = p.distance().km
        
        # Logic:
        # If we flipped TO 1 (True), we started moving away.
        # That means we were just at the closest point -> Perihelion.
        if state_after_flip == 1:
            name = "Earth Perihelion"
        
        # If we flipped TO 0 (False), we started moving closer.
        # That means we were just at the furthest point -> Aphelion.
        else:
            name = "Earth Aphelion"
            
        events.append({
            "datetime": t.utc_iso(),
            "event": name,
            "details": f"{dist_km:.2f} km"
        })

    # Output to CSV
    output_filename = "earth_extrema.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        
        for event in events:
            writer.writerow([
                event["datetime"],
                event["event"],
                event["details"]
            ])
            
    print(f"Success. Wrote {len(events)} precise events to {output_filename}.")

if __name__ == "__main__":
    calculate_earth_extrema()
