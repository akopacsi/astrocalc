import json
import csv
import os
from datetime import datetime
import numpy as np
from skyfield.api import load, Loader, wgs84
from skyfield import almanac
from skyfield import searchlib
from skyfield.functions import angle_between

# Physical constants for angular size calculation
SUN_RADIUS_KM = 696340.0
MOON_RADIUS_KM = 1737.1

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

def get_angular_radius(distance_km, radius_km):
    """
    Calculates the angular radius of a body in degrees.
    """
    return np.degrees(np.arcsin(radius_km / distance_km))

def calculate_eclipses():
    config = load_config()
    eph_filename = config.get("ephemeris_file", "de440.bsp")
    eph = ensure_de440(eph_filename)
    
    sun = eph['sun']
    earth = eph['earth']
    moon = eph['moon']
    
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
    
    print("Scanning for Solar and Lunar Eclipses...")

    # 1. Global Scan using Moon Phases
    t_phases, y_phases = almanac.find_discrete(
        t_start, t_end, almanac.moon_phases(eph)
    )

    events = []

    for t, phase in zip(t_phases, y_phases):
        
        # ==========================================
        # SOLAR ECLIPSE CHECK (New Moon)
        # ==========================================
        if phase == 0:
            # 1. Global Pre-check
            g_sun = earth.at(t).observe(sun)
            g_moon = earth.at(t).observe(moon)
            sep = g_sun.separation_from(g_moon).degrees
            
            if sep > 1.6: 
                continue

            # 2. Local Visibility Check
            def local_overlap(time_obj):
                obs = observer.at(time_obj)
                s = obs.observe(sun).apparent()
                m = obs.observe(moon).apparent()
                
                curr_sep = s.separation_from(m).degrees
                r_sun = get_angular_radius(s.distance().km, SUN_RADIUS_KM)
                r_moon = get_angular_radius(m.distance().km, MOON_RADIUS_KM)
                
                return curr_sep < (r_sun + r_moon)

            t_check_start = ts.tt_jd(t.tt - 0.17)
            t_check_end = ts.tt_jd(t.tt + 0.17)
            
            local_overlap.rough_period = 0.1
            t_contacts, values = searchlib.find_discrete(t_check_start, t_check_end, local_overlap)
            
            # --- VISIBLE SOLAR ECLIPSE ---
            if 1 in values:
                t_beg = None
                t_end_local = None
                for tc, val in zip(t_contacts, values):
                    if val == 1: t_beg = tc
                    if val == 0: t_end_local = tc
                
                if t_beg is not None and t_end_local is not None:
                    def separation_func(tx):
                        obs = observer.at(tx)
                        return obs.observe(sun).separation_from(obs.observe(moon)).degrees
                    
                    separation_func.rough_period = 0.1
                    t_max_list, min_seps = searchlib.find_minima(t_beg, t_end_local, separation_func)
                    
                    if t_max_list:
                        tm = t_max_list[0]
                        sep_min = min_seps[0]
                        
                        obs = observer.at(tm)
                        s = obs.observe(sun).apparent()
                        m = obs.observe(moon).apparent()
                        r_sun = get_angular_radius(s.distance().km, SUN_RADIUS_KM)
                        r_moon = get_angular_radius(m.distance().km, MOON_RADIUS_KM)
                        
                        mag = (r_moon + r_sun - sep_min) / (2 * r_sun)
                        
                        if mag >= 1.0:
                            event_name = "Solar Eclipse (Total)"
                        else:
                            if (r_moon < r_sun) and (sep_min < (r_sun - r_moon)):
                                event_name = "Solar Eclipse (Annular)"
                            else:
                                event_name = "Solar Eclipse (Partial)"

                        details_str = (
                            f"Begin: {t_beg.utc_strftime('%H:%M')}, "
                            f"Max: {tm.utc_strftime('%H:%M')}, "
                            f"End: {t_end_local.utc_strftime('%H:%M')}, "
                            f"Mag: {mag:.3f}"
                        )
                        
                        events.append({
                            "datetime": tm.utc_iso(),
                            "event": event_name,
                            "details": details_str
                        })
            
            # --- NOT VISIBLE SOLAR ECLIPSE ---
            else:
                events.append({
                    "datetime": t.utc_iso(),
                    "event": "Solar Eclipse (Not visible)",
                    "details": "" 
                })

        # ==========================================
        # LUNAR ECLIPSE CHECK (Full Moon)
        # ==========================================
        elif phase == 2:
            p_sun = earth.at(t).observe(sun).apparent()
            p_moon = earth.at(t).observe(moon).apparent()
            pos_sun = p_sun.position.au
            pos_moon = p_moon.position.au
            
            sep_shadow = angle_between(pos_moon, -pos_sun) * 57.2958
            
            if sep_shadow > 1.0:
                continue

            # Umbra Contact Check
            def moon_in_umbra(tx):
                sun_p = earth.at(tx).observe(sun).apparent().position.au
                moon_p = earth.at(tx).observe(moon).apparent().position.au
                sep = angle_between(moon_p, -sun_p) * 57.2958
                return sep < (0.65 + 0.26)

            t_check_start = ts.tt_jd(t.tt - 0.17)
            t_check_end = ts.tt_jd(t.tt + 0.17)
            moon_in_umbra.rough_period = 0.1
            
            t_contacts, values = searchlib.find_discrete(t_check_start, t_check_end, moon_in_umbra)
            
            t_u1 = None
            t_u4 = None
            for tc, val in zip(t_contacts, values):
                if val == 1: t_u1 = tc
                if val == 0: t_u4 = tc
            
            if t_u1 is not None and t_u4 is not None:
                def shadow_sep_func(tx):
                    sun_p = earth.at(tx).observe(sun).apparent().position.au
                    moon_p = earth.at(tx).observe(moon).apparent().position.au
                    return angle_between(moon_p, -sun_p)
                
                shadow_sep_func.rough_period = 0.1
                t_max_list, min_seps = searchlib.find_minima(t_u1, t_u4, shadow_sep_func)
                
                if t_max_list:
                    tm = t_max_list[0]
                    sep_deg = min_seps[0] * 57.2958
                    
                    r_moon = 0.26
                    r_umbra = 0.65 
                    mag = (r_umbra + r_moon - sep_deg) / (2 * r_moon)
                    
                    # --- VISIBILITY CHECK ---
                    # Check altitude at the moment of maximum eclipse
                    obs_at_max = observer.at(tm).observe(moon).apparent()
                    alt, _, _ = obs_at_max.altaz()
                    
                    if alt.degrees < 0:
                        # NOT VISIBLE -> HIDE DATA
                        events.append({
                            "datetime": tm.utc_iso(),
                            "event": "Lunar Eclipse (Not visible)",
                            "details": ""
                        })
                    else:
                        # VISIBLE -> SHOW DATA
                        if mag >= 1.0:
                            event_name = "Lunar Eclipse (Total)"
                        else:
                            event_name = "Lunar Eclipse (Partial)"
                        
                        details_str = (
                            f"Begin: {t_u1.utc_strftime('%H:%M')}, "
                            f"Max: {tm.utc_strftime('%H:%M')}, "
                            f"End: {t_u4.utc_strftime('%H:%M')}, "
                            f"Mag: {mag:.2f}"
                        )
                        
                        events.append({
                            "datetime": tm.utc_iso(),
                            "event": event_name,
                            "details": details_str
                        })

    events.sort(key=lambda x: x["datetime"])

    output_filename = "eclipses.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["datetime", "event", "details"])
        
        for event in events:
            writer.writerow([
                event["datetime"],
                event["event"],
                event["details"]
            ])
            
    print(f"Success. Wrote {len(events)} eclipse events to {output_filename}.")

if __name__ == "__main__":
    calculate_eclipses()
