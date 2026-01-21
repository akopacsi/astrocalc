import json
import csv
import os
from datetime import datetime, timedelta
import numpy as np
from skyfield.api import load, Loader, wgs84
from skyfield import almanac

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

def format_time_iso(t):
    """Returns HH:MM string from Skyfield Time object."""
    return t.utc_strftime('%H:%M')

def format_gst_vector(gast_hours):
    """
    Vectorized conversion of GAST (hours) to HH:MM:SS strings.
    """
    h = gast_hours.astype(int)
    remainder_h = (gast_hours - h) * 60.0
    m = remainder_h.astype(int)
    s = ((remainder_h - m) * 60.0).astype(int)
    # Handle rounding edge cases (e.g. 59.999s)
    return [f"{hh:02d}:{mm:02d}:{ss:02d}" for hh, mm, ss in zip(h, m, s)]

def calculate_calendar_final():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    moon = eph['moon']
    
    loc = config["location"]
    
    # 1. Location & Observer
    # Location (Topos) is needed for almanac search functions
    location = wgs84.latlon(
        loc["latitude"], loc["longitude"], elevation_m=loc["elevation_m"]
    )
    # Observer (Vector) is needed for altitude calculations
    observer = earth + location

    ts = load.timescale()
    
    # 2. Time Range Setup
    start_str = config["time_range"]["start"].replace('Z', '+00:00')
    end_str = config["time_range"]["end"].replace('Z', '+00:00')
    start_dt = datetime.fromisoformat(start_str).replace(hour=0, minute=0, second=0)
    end_dt = datetime.fromisoformat(end_str).replace(hour=0, minute=0, second=0)

    # We add a buffer of +/- 1 day to the search range to catch events
    # that occur near midnight boundaries.
    t_search_start = ts.from_datetime(start_dt - timedelta(days=1))
    t_search_end = ts.from_datetime(end_dt + timedelta(days=1))

    print(f"Generating Calendar from {start_dt.date()} to {end_dt.date()}...")

    # =================================================================
    # PHASE 1: BATCH ALMANAC CALCULATIONS (Events)
    # =================================================================
    
    # A. Sun Rise/Set
    print("  > Calculating Sun Rise/Set...")
    f_sun = almanac.risings_and_settings(eph, sun, location)
    t_sun, y_sun = almanac.find_discrete(t_search_start, t_search_end, f_sun)

    # B. Moon Rise/Set
    print("  > Calculating Moon Rise/Set...")
    f_moon = almanac.risings_and_settings(eph, moon, location)
    t_moon, y_moon = almanac.find_discrete(t_search_start, t_search_end, f_moon)

    # C. Sun Transits (Upper Only)
    print("  > Calculating Sun Transits...")
    f_trans_sun = almanac.meridian_transits(eph, sun, location)
    t_tr_sun_raw, y_tr_sun_raw = almanac.find_discrete(t_search_start, t_search_end, f_trans_sun)
    # Filter: y=1 is Upper Transit (Noon), y=0 is Lower Transit (Midnight)
    mask_sun = (y_tr_sun_raw == 1)
    t_tr_sun = t_tr_sun_raw[mask_sun]

    # D. Moon Transits (Upper Only)
    print("  > Calculating Moon Transits...")
    f_trans_moon = almanac.meridian_transits(eph, moon, location)
    t_tr_moon_raw, y_tr_moon_raw = almanac.find_discrete(t_search_start, t_search_end, f_trans_moon)
    mask_moon = (y_tr_moon_raw == 1)
    t_tr_moon = t_tr_moon_raw[mask_moon]

    # E. Moon Phases
    print("  > Calculating Moon Phases...")
    f_phases = almanac.moon_phases(eph)
    t_phases, y_phases = almanac.find_discrete(t_search_start, t_search_end, f_phases)

    # =================================================================
    # PHASE 2: DAILY ITERATION & MAPPING
    # =================================================================
    
    total_days = (end_dt.date() - start_dt.date()).days + 1
    date_list = [start_dt.date() + timedelta(days=i) for i in range(total_days)]
    
    # Pre-calculate Vectorized Times for Noon and Midnight
    noon_dts = [datetime(d.year, d.month, d.day, 12, 0, 0, tzinfo=start_dt.tzinfo) for d in date_list]
    midnight_dts = [datetime(d.year, d.month, d.day, 0, 0, 0, tzinfo=start_dt.tzinfo) for d in date_list]
    
    t_noons = ts.from_datetimes(noon_dts)
    t_midnights = ts.from_datetimes(midnight_dts)
    
    # --- Vectorized Calculations ---
    # 1. Julian Date at 12:00 UTC (Terrestrial Time)
    jd_array = t_noons.tt 
    
    # 2. Greenwich Sidereal Time at 00:00 UTC
    gst_array = t_midnights.gast
    gst_strings = format_gst_vector(gst_array)
    
    # 3. Equation of Time at 12:00 UTC
    # Eq = GHA_Sun - 12h
    gast_noons = t_noons.gast
    ra_noons = earth.at(t_noons).observe(sun).apparent().radec(epoch='date')[0].hours
    gha_hours = (gast_noons - ra_noons) % 24.0
    # Normalize to [-12, +12]
    gha_hours = np.where(gha_hours > 12.0, gha_hours - 24.0, gha_hours)
    # Skyfield GHA is 'how far past the meridian'. 
    # If GHA is 0.1h, Sun passed 0.1h ago (Sun is fast/early). EoT is positive.
    eot_minutes = gha_hours * 60.0

    # 4. Sun Altitude at Transit
    # We observe the sun at the EXACT moment of transit found in Phase 1
    # But we need to map these altitudes to the specific day later.
    if len(t_tr_sun) > 0:
        sun_trans_altitudes = observer.at(t_tr_sun).observe(sun).apparent().altaz()[0].degrees
    else:
        sun_trans_altitudes = []

    # --- Data Bucketing ---
    # We create a map: Date_String -> Event_Data
    daily_events = {d.isoformat(): {} for d in date_list}

    def bucket_event(t_objs, values, key_map):
        iso_strs = t_objs.utc_iso()
        for i, t_str in enumerate(iso_strs):
            d_key = t_str[:10] # YYYY-MM-DD
            if d_key in daily_events:
                time_val = t_str[11:16]
                if values is not None:
                    # Map code (1=Rise) to string key
                    code = values[i]
                    if code in key_map:
                        daily_events[d_key][key_map[code]] = time_val
                else:
                    # Direct assignment (Transit)
                    daily_events[d_key][key_map] = time_val
                    # Special case for Sun Transit Altitude
                    if key_map == "sun_transit":
                        daily_events[d_key]["sun_alt"] = f"{sun_trans_altitudes[i]:.1f}"

    # Map Events to Dates
    bucket_event(t_sun, y_sun, {1: "sunrise", 0: "sunset"})
    bucket_event(t_moon, y_moon, {1: "moonrise", 0: "moonset"})
    bucket_event(t_tr_sun, None, "sun_transit")
    bucket_event(t_tr_moon, None, "moon_transit")

    # Map Phases
    phase_names = {0: "New Moon", 1: "First Q", 2: "Full Moon", 3: "Last Q"}
    iso_ph = t_phases.utc_iso()
    for i, t_str in enumerate(iso_ph):
        d_key = t_str[:10]
        if d_key in daily_events:
            p_name = phase_names[y_phases[i]]
            p_time = t_str[11:16]
            daily_events[d_key]["moon_phase"] = f"{p_name} {p_time}"

    # =================================================================
    # PHASE 3: WRITE CSV
    # =================================================================
    rows = []
    days_abbr = ["Mo", "Tu", "We", "Th", "Fr", "Sa", "Su"]

    for i, d in enumerate(date_list):
        d_str = d.isoformat()
        evt = daily_events[d_str]
        
        # Day of Week & Day of Year
        dow = days_abbr[d.weekday()]
        doy = d.timetuple().tm_yday
        
        # Pre-computed vectors
        jd = f"{jd_array[i]:.2f}"
        gst = gst_strings[i]
        eot = f"{eot_minutes[i]:+.1f}m"

        rows.append([
            d_str,
            dow,
            doy,
            evt.get("sunrise", "--:--"),
            evt.get("sun_transit", "--:--"),
            evt.get("sunset", "--:--"),
            evt.get("sun_alt", ""),
            eot,
            evt.get("moonrise", "--:--"),
            evt.get("moon_transit", "--:--"),
            evt.get("moonset", "--:--"),
            evt.get("moon_phase", ""),
            jd,
            gst
        ])

    filename = "daily_calendar.csv"
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "date", "dow", "doy", 
            "sunrise", "sun_transit", "sunset", "sun_alt_trans", "eq_of_time",
            "moonrise", "moon_transit", "moonset", "moon_phase",
            "julian_date_12h", "greenwich_st_00h"
        ])
        writer.writerows(rows)

    print(f"Success. Wrote {len(rows)} days to {filename}.")

if __name__ == "__main__":
    calculate_calendar_final()