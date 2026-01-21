import json
import csv
import os
from datetime import datetime
import numpy as np
from skyfield.api import load, Loader, wgs84
from skyfield import almanac
from skyfield import searchlib
from skyfield.magnitudelib import planetary_magnitude

# Physical constants
SUN_RADIUS_KM = 696340.0
# Planet radii (km) for diameter/transit calculations
PLANET_RADII = {
    'mercury': 2439.7,
    'venus': 6051.8
}

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

def get_angular_radius_deg(distance_km, radius_km):
    return np.degrees(np.arcsin(radius_km / distance_km))

def get_visual_details(observer, t, planet_obj, planet_name_key):
    """
    Returns a formatted string with Magnitude and Apparent Diameter.
    Uses 'arcsec' instead of double quotes to prevent CSV formatting issues.
    """
    obs = observer.at(t).observe(planet_obj)
    app = obs.apparent()
    
    # 1. Magnitude
    try:
        mag = planetary_magnitude(app)
        mag_str = f"{mag:.1f}m"
    except ValueError:
        mag_str = "?m"

    # 2. Apparent Diameter (in arcseconds)
    dist_km = obs.distance().km
    radius_km = PLANET_RADII[planet_name_key.lower()]
    # angular diameter = 2 * arctan(r/d)
    diam_deg = 2 * np.degrees(np.arctan(radius_km / dist_km))
    diam_arcsec = diam_deg * 3600.0
    
    return f"Mag: {mag_str}, Dia: {diam_arcsec:.1f} arcsec"

def calculate_inner_planet_events():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    
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
    target_planets = ['mercury', 'venus']

    for planet_name in target_planets:
        print(f"Calculating events for {planet_name.capitalize()}...")
        planet = eph[planet_name]
        
        if planet_name == 'mercury':
            synodic_days = 115.0
        else:
            synodic_days = 584.0

        # ==========================================
        # 1. PERIGEE & APOGEE
        # ==========================================
        def earth_planet_distance(t):
            return earth.at(t).observe(planet).distance().km
        
        earth_planet_distance.rough_period = synodic_days

        # Perigee
        t_min, val_min = searchlib.find_minima(t_start, t_end, earth_planet_distance)
        for t, dist in zip(t_min, val_min):
            vis = get_visual_details(observer, t, planet, planet_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{planet_name.capitalize()} Perigee (Closest)",
                "details": f"Dist: {dist:.0f} km, {vis}"
            })

        # Apogee
        t_max, val_max = searchlib.find_maxima(t_start, t_end, earth_planet_distance)
        for t, dist in zip(t_max, val_max):
            vis = get_visual_details(observer, t, planet, planet_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{planet_name.capitalize()} Apogee (Farthest)",
                "details": f"Dist: {dist:.0f} km, {vis}"
            })

        # ==========================================
        # 2. GREATEST ELONGATIONS
        # ==========================================
        def sun_planet_separation(t):
            e = earth.at(t)
            s = e.observe(sun)
            p = e.observe(planet)
            return s.separation_from(p).degrees

        sun_planet_separation.rough_period = synodic_days

        t_elong, val_elong = searchlib.find_maxima(t_start, t_end, sun_planet_separation)
        
        for t, sep in zip(t_elong, val_elong):
            e = earth.at(t)
            _, lon_sun, _ = e.observe(sun).apparent().ecliptic_latlon()
            _, lon_planet, _ = e.observe(planet).apparent().ecliptic_latlon()
            
            diff = (lon_planet.degrees - lon_sun.degrees) % 360
            if 0 < diff < 180:
                direction = "Eastern (Evening)"
            else:
                direction = "Western (Morning)"

            vis = get_visual_details(observer, t, planet, planet_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{planet_name.capitalize()} Greatest {direction} Elongation",
                "details": f"Sep: {sep:.1f}°, {vis}"
            })

        # ==========================================
        # 3. CONJUNCTIONS & TRANSITS
        # ==========================================
        def is_planet_ahead(t):
            e = earth.at(t)
            _, lon_sun, _ = e.observe(sun).apparent().ecliptic_latlon()
            _, lon_planet, _ = e.observe(planet).apparent().ecliptic_latlon()
            diff = (lon_planet.degrees - lon_sun.degrees) % 360
            return diff < 180 

        is_planet_ahead.rough_period = synodic_days
        
        t_conjs, vals = searchlib.find_discrete(t_start, t_end, is_planet_ahead)
        
        for t in t_conjs:
            d_sun = earth.at(t).observe(sun).distance().au
            d_planet = earth.at(t).observe(planet).distance().au
            
            # Inferior Conjunction
            if d_planet < d_sun:
                type_name = "Inferior Conjunction"
                
                # Check for Transit
                def transit_sep(tx):
                    e = earth.at(tx)
                    return e.observe(sun).separation_from(e.observe(planet)).degrees
                
                transit_sep.rough_period = 0.1
                curr_sep = transit_sep(t)
                
                if curr_sep < 0.27:
                    # Detailed Transit Check
                    def local_transit_overlap(tx):
                        obs = observer.at(tx)
                        s = obs.observe(sun).apparent()
                        p = obs.observe(planet).apparent()
                        sep = s.separation_from(p).degrees
                        r_sun = get_angular_radius_deg(s.distance().km, SUN_RADIUS_KM)
                        r_p = get_angular_radius_deg(p.distance().km, PLANET_RADII[planet_name])
                        return sep < (r_sun + r_p)

                    t_scan_start = ts.tt_jd(t.tt - 0.5)
                    t_scan_end = ts.tt_jd(t.tt + 0.5)
                    local_transit_overlap.rough_period = 0.1
                    
                    t_cont, vals_cont = searchlib.find_discrete(t_scan_start, t_scan_end, local_transit_overlap)
                    
                    # If we found local overlap events (entry/exit)
                    if 1 in vals_cont:
                        t_beg = None
                        t_end_local = None
                        for tc, val in zip(t_cont, vals_cont):
                            if val == 1: t_beg = tc
                            if val == 0: t_end_local = tc
                        
                        if t_beg is not None and t_end_local is not None:
                            t_max_list, min_seps = searchlib.find_minima(t_beg, t_end_local, transit_sep)
                            if t_max_list:
                                tm = t_max_list[0]
                                alt_sun = observer.at(tm).observe(sun).apparent().altaz()[0].degrees
                                
                                # VISIBILITY CHECK
                                if alt_sun > 0:
                                    vis = get_visual_details(observer, tm, planet, planet_name)
                                    details = (
                                        f"Start: {t_beg.utc_strftime('%H:%M')}, "
                                        f"End: {t_end_local.utc_strftime('%H:%M')}, "
                                        f"{vis}"
                                    )
                                    events.append({
                                        "datetime": tm.utc_iso(),
                                        "event": f"{planet_name.capitalize()} Transit of Sun",
                                        "details": details
                                    })
                                else:
                                    # Not Visible -> No Details
                                    events.append({
                                        "datetime": tm.utc_iso(),
                                        "event": f"{planet_name.capitalize()} Transit of Sun (Not visible)",
                                        "details": "" 
                                    })

            else:
                type_name = "Superior Conjunction"

            vis = get_visual_details(observer, t, planet, planet_name)
            events.append({
                "datetime": t.utc_iso(),
                "event": f"{planet_name.capitalize()} {type_name}",
                "details": vis
            })

    events.sort(key=lambda x: x["datetime"])

    output_filename = "inner_planets.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        for e in events:
            writer.writerow([e["datetime"], e["event"], e["details"]])
            
    print(f"Success. Wrote {len(events)} inner planet events to {output_filename}.")

if __name__ == "__main__":
    calculate_inner_planet_events()