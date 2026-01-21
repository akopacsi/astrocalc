import json
import csv
import os
from datetime import datetime
import numpy as np
from skyfield.api import load, Loader, wgs84
from skyfield import almanac
from skyfield import searchlib
from skyfield.magnitudelib import planetary_magnitude

# Planet radii (km) for diameter calculations
PLANET_RADII = {
    'mars': 3389.5,
    'jupiter': 69911.0,
    'saturn': 58232.0,
    'uranus': 25362.0,
    'neptune': 24622.0
}

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

def get_visual_details(observer, t, planet_obj, planet_name_key):
    """
    Returns a formatted string with Magnitude and Apparent Diameter.
    Uses 'arcsec' to prevent CSV errors.
    """
    obs = observer.at(t).observe(planet_obj)
    app = obs.apparent()
    
    # 1. Magnitude
    try:
        mag = planetary_magnitude(app)
        mag_str = f"{mag:.1f}m"
    except ValueError:
        mag_str = "?m"

    # 2. Apparent Diameter (arcseconds)
    dist_km = obs.distance().km
    radius_km = PLANET_RADII[planet_name_key.lower()]
    # angular diameter = 2 * arctan(r/d)
    diam_deg = 2 * np.degrees(np.arctan(radius_km / dist_km))
    diam_arcsec = diam_deg * 3600.0
    
    return f"Mag: {mag_str}, Dia: {diam_arcsec:.1f} arcsec"

def calculate_outer_planet_events():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    loc = config["location"]
    # Define observer on Earth surface
    observer = eph['earth'] + wgs84.latlon(
        loc["latitude"], loc["longitude"], elevation_m=loc["elevation_m"]
    )
    # Earth center for pure distance calculations
    earth_center = eph['earth'] 
    
    ts = load.timescale()
    t_start = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["start"].replace('Z', '+00:00'))
    )
    t_end = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["end"].replace('Z', '+00:00'))
    )

    events = []
    
    # Map display names to keys in de440.bsp
    outer_planets_map = {
        'Mars': 'MARS BARYCENTER',
        'Jupiter': 'JUPITER BARYCENTER',
        'Saturn': 'SATURN BARYCENTER',
        'Uranus': 'URANUS BARYCENTER',
        'Neptune': 'NEPTUNE BARYCENTER'
    }

    for display_name, eph_key in outer_planets_map.items():
        print(f"Calculating events for {display_name}...")
        planet = eph[eph_key]
        
        periods = { 'Mars': 780.0, 'Jupiter': 399.0, 'Saturn': 378.0, 'Uranus': 370.0, 'Neptune': 367.0 }
        rough_p = periods[display_name]

        # ==========================================
        # 1. OPPOSITION & CONJUNCTION
        # ==========================================
        f = almanac.oppositions_conjunctions(eph, planet)
        t_oc, y_oc = almanac.find_discrete(t_start, t_end, f)

        for t, val in zip(t_oc, y_oc):
            event_type = "Opposition" if val == 1 else "Conjunction with Sun"
            vis = get_visual_details(observer, t, planet, display_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{display_name} {event_type}",
                "details": vis
            })

        # ==========================================
        # 2. PERIGEE & APOGEE
        # ==========================================
        def distance_func(t):
            return earth_center.at(t).observe(planet).distance().au
        distance_func.rough_period = rough_p

        # Perigee (Closest)
        t_peri, dist_peri = searchlib.find_minima(t_start, t_end, distance_func)
        for t, d in zip(t_peri, dist_peri):
            vis = get_visual_details(observer, t, planet, display_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{display_name} Perigee",
                "details": f"Dist: {d:.4f} AU, {vis}"
            })

        # Apogee (Farthest)
        t_apo, dist_apo = searchlib.find_maxima(t_start, t_end, distance_func)
        for t, d in zip(t_apo, dist_apo):
            vis = get_visual_details(observer, t, planet, display_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{display_name} Apogee",
                "details": f"Dist: {d:.4f} AU, {vis}"
            })

    events.sort(key=lambda x: x["datetime"])

    output_filename = "outer_planets.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        for e in events:
            writer.writerow([e["datetime"], e["event"], e["details"]])
            
    print(f"Success. Wrote {len(events)} events to {output_filename}.")

if __name__ == "__main__":
    calculate_outer_planet_events()