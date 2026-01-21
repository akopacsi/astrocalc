import json
import csv
import os
from datetime import datetime
from skyfield.api import load, Loader, Topos
from skyfield import almanac
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

def find_sun_at_minus_6_deg(ts, t_start, t_end, observer, sun):
    """
    Finds exact times when Sun crosses -6 degrees altitude (Civil Twilight).
    Returns list of (time, type) tuples where type is "Dawn" or "Dusk".
    """
    def is_sun_above_minus_6(t):
        astrometric = observer.at(t).observe(sun)
        app = astrometric.apparent()
        alt, az, distance = app.altaz()
        return alt.degrees > -6.0

    is_sun_above_minus_6.rough_period = 0.5 
    
    t_times, values = searchlib.find_discrete(t_start, t_end, is_sun_above_minus_6)
    
    results = []
    for t, state_after_change in zip(t_times, values):
        if state_after_change == 1: # False -> True (Rising/Dawn)
            event_type = "Dawn"
        else: # True -> False (Setting/Dusk)
            event_type = "Dusk"
        results.append((t, event_type))
        
    return results

def calculate_moon_events():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    moon = eph['moon']
    
    loc = config["location"]
    observer = earth + Topos(
        latitude_degrees=loc["latitude"],
        longitude_degrees=loc["longitude"],
        elevation_m=loc["elevation_m"]
    )

    ts = load.timescale()
    t_start = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["start"].replace('Z', '+00:00'))
    )
    t_end = ts.from_datetime(
        datetime.fromisoformat(config["time_range"]["end"].replace('Z', '+00:00'))
    )
    
    print("Calculating Moon Phases, Supermoons, and Crescent Visibility...")
    
    # 1. Calculate Standard Phases
    t_phases, phase_codes = almanac.find_discrete(
        t_start, t_end, almanac.moon_phases(eph)
    )
    
    phase_names = {
        0: "New Moon",
        1: "First Quarter",
        2: "Full Moon",
        3: "Last Quarter"
    }

    events = []

    for t, code in zip(t_phases, phase_codes):
        name = phase_names[code]
        details = ""
        
        # --- SUPERMOON LOGIC (Full Moon Only) ---
        if code == 2: # Full Moon
            dist_km = (earth - moon).at(t).distance().km
            
            # CRITERIA UPDATE:
            # The Jan 3, 2026 Full Moon is 362,312 km.
            # To include this (and align with Espenak/Media), we set the limit to 363,000 km.
            if dist_km < 363_000:
                details = f"Supermoon (Dist: {dist_km:.0f} km)"
            else:
                details = "" 
                
            events.append({
                "datetime": t.utc_iso(),
                "event": name,
                "details": details
            })

        # --- CRESCENT VISIBILITY LOGIC (New Moon Only) ---
        elif code == 0: # New Moon
            events.append({
                "datetime": t.utc_iso(),
                "event": name,
                "details": ""
            })
            
            # A. Find Last OLD Moon Crescent (Morning)
            for day_offset in range(1, 5):
                t_check_start = ts.tt_jd(t.tt - day_offset - 0.5)
                t_check_end = ts.tt_jd(t.tt - day_offset + 0.5)
                
                sun_events = find_sun_at_minus_6_deg(ts, t_check_start, t_check_end, observer, sun)
                dawns = [x[0] for x in sun_events if x[1] == "Dawn"]
                
                if not dawns: continue
                t_dawn = dawns[0]
                
                moon_at = observer.at(t_dawn).observe(moon).apparent()
                alt, az, _ = moon_at.altaz()
                
                if alt.degrees > 0:
                    age_days = t_dawn - t # Negative
                    events.append({
                        "datetime": t_dawn.utc_iso(),
                        "event": "Old Moon Crescent Visible",
                        "details": f"Alt: {alt.degrees:.1f}°, Az: {az.degrees:.0f}°, Age: {age_days:.2f}d"
                    })
                    break 
            
            # B. Find First YOUNG Moon Crescent (Evening)
            for day_offset in range(1, 5):
                t_check_start = ts.tt_jd(t.tt + day_offset - 0.5)
                t_check_end = ts.tt_jd(t.tt + day_offset + 0.5)
                
                sun_events = find_sun_at_minus_6_deg(ts, t_check_start, t_check_end, observer, sun)
                dusks = [x[0] for x in sun_events if x[1] == "Dusk"]
                
                if not dusks: continue
                t_dusk = dusks[0]
                
                moon_at = observer.at(t_dusk).observe(moon).apparent()
                alt, az, _ = moon_at.altaz()
                
                if alt.degrees > 0:
                    age_days = t_dusk - t # Positive
                    events.append({
                        "datetime": t_dusk.utc_iso(),
                        "event": "Young Moon Crescent Visible",
                        "details": f"Alt: {alt.degrees:.1f}°, Az: {az.degrees:.0f}°, Age: +{age_days:.2f}d"
                    })
                    break 

        else:
            events.append({
                "datetime": t.utc_iso(),
                "event": name,
                "details": ""
            })

    events.sort(key=lambda x: x["datetime"])

    output_filename = "moon_phases_detailed.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        
        for event in events:
            writer.writerow([
                event["datetime"],
                event["event"],
                event["details"]
            ])
            
    print(f"Success. Wrote {len(events)} events to {output_filename}.")

if __name__ == "__main__":
    calculate_moon_events()
