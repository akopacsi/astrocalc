import json
import csv
import os
import itertools
from datetime import datetime
from skyfield.api import load, Loader, wgs84
from skyfield import searchlib

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

def calculate_conjunctions():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    moon = eph['moon']
    
    # Define Planets (using Barycenters for safety with de440)
    # de440 usually has:
    # 1=Mercury, 2=Venus, 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune (Barycenters)
    # 199/299/399/499 are specific planet centers but often missing for outer ones.
    # Barycenters are accurate enough for visual conjunctions.
    planets = {
        'Mercury': eph['MERCURY BARYCENTER'],
        'Venus': eph['VENUS BARYCENTER'],
        'Mars': eph['MARS BARYCENTER'],
        'Jupiter': eph['JUPITER BARYCENTER'],
        'Saturn': eph['SATURN BARYCENTER'],
        'Uranus': eph['URANUS BARYCENTER'],
        'Neptune': eph['NEPTUNE BARYCENTER']
    }
    
    # We define the observer on Earth
    loc = config["location"]
    observer = earth + wgs84.latlon(
        loc["latitude"], loc["longitude"], elevation_m=loc["elevation_m"]
    )

    ts = load.timescale()
    t_start = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["start"].replace('Z', '+00:00'))
    )
    t_end = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["end"].replace('Z', '+00:00'))
    )

    events = []
    
    print("Calculating Close Approaches...")

    # ==========================================
    # 1. MOON vs PLANETS (Threshold < 6.0 deg)
    # ==========================================
    print("  > Scanning Moon-Planet conjunctions...")
    
    for p_name, p_obj in planets.items():
        # Define Separation Function
        def moon_planet_sep(t):
            obs = observer.at(t)
            # Apparent position (accounting for light time, aberration, deflection)
            m = obs.observe(moon).apparent()
            p = obs.observe(p_obj).apparent()
            return m.separation_from(p).degrees
        
        # Period hint: Moon circles Earth every ~27.3 days
        moon_planet_sep.rough_period = 27.3
        
        # Find Minima
        t_mins, val_mins = searchlib.find_minima(t_start, t_end, moon_planet_sep)
        
        for t, sep in zip(t_mins, val_mins):
            if sep < 6.0:
                events.append({
                    "datetime": t.utc_iso(),
                    "event": f"Conjunction: Moon & {p_name}",
                    "details": f"Separation: {sep:.2f}°"
                })

    # ==========================================
    # 2. PLANET vs PLANETS (Threshold < 3.0 deg)
    # ==========================================
    print("  > Scanning Planet-Planet conjunctions...")
    
    # Get all unique pairs of planets
    planet_pairs = list(itertools.combinations(planets.keys(), 2))
    
    for p1_name, p2_name in planet_pairs:
        p1_obj = planets[p1_name]
        p2_obj = planets[p2_name]
        
        def planet_planet_sep(t):
            obs = observer.at(t)
            p1 = obs.observe(p1_obj).apparent()
            p2 = obs.observe(p2_obj).apparent()
            return p1.separation_from(p2).degrees
        
        # Period hint: Inner planets move fast. 
        # 40 days is small enough to catch Mercury/Venus events 
        # without excessive computation.
        planet_planet_sep.rough_period = 40.0
        
        t_mins, val_mins = searchlib.find_minima(t_start, t_end, planet_planet_sep)
        
        for t, sep in zip(t_mins, val_mins):
            if sep < 3.0:
                # Solar Proximity Check:
                # If both planets are too close to the Sun, the event is invisible.
                # Let's check Sun separation for one of them.
                # (Optional, but useful. We will list all, but maybe note if invisible?)
                # For now, just listing all geometry as requested.
                
                events.append({
                    "datetime": t.utc_iso(),
                    "event": f"Conjunction: {p1_name} & {p2_name}",
                    "details": f"Separation: {sep:.2f}°"
                })

    # Sort and Output
    events.sort(key=lambda x: x["datetime"])

    output_filename = "conjunctions.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        
        for event in events:
            writer.writerow([
                event["datetime"],
                event["event"],
                event["details"]
            ])
            
    print(f"Success. Wrote {len(events)} conjunctions to {output_filename}.")

if __name__ == "__main__":
    calculate_conjunctions()