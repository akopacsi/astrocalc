import csv
import os
import sys

# Import your calculation modules
import calc_seasons
import calc_earth_extrema
import calc_moon_phases
import calc_eclipses
import calc_inner_planets
import calc_outer_planets
import calc_conjunctions
import calc_daily_calendar

def main():
    print("==========================================")
    print("      ASTROCALC: MASTER CALCULATION       ")
    print("==========================================")

    # ---------------------------------------------------------
    # 1. EXECUTE EVENT SCRIPTS
    # ---------------------------------------------------------
    print("\n[1/7] Calculating Seasons...")
    calc_seasons.calculate_seasons()
    
    print("\n[2/7] Calculating Earth Perihelion/Aphelion...")
    calc_earth_extrema.calculate_earth_extrema()
    
    print("\n[3/7] Calculating Moon Phases & Visibility...")
    calc_moon_phases.calculate_moon_events()
    
    print("\n[4/7] Calculating Eclipses...")
    calc_eclipses.calculate_eclipses()
    
    print("\n[5/7] Calculating Inner Planets...")
    calc_inner_planets.calculate_inner_planet_events()
    
    print("\n[6/7] Calculating Outer Planets...")
    calc_outer_planets.calculate_outer_planet_events()
    
    print("\n[7/7] Calculating Conjunctions...")
    calc_conjunctions.calculate_conjunctions()

    # ---------------------------------------------------------
    # 2. ASK USER: GENERATE DAILY CALENDAR?
    # ---------------------------------------------------------
    print("\n------------------------------------------")
    response = input("Generate detailed Daily Calendar (daily_calendar.csv)? [y/N]: ").strip().lower()

    if response in ['y', 'yes']:
        print("\n[8/8] Generating Daily Data Grid...")
        try:
            # Tries to call the latest function name from your final script
            calc_daily_calendar.calculate_calendar_final()
        except AttributeError:
            try:
                calc_daily_calendar.calculate_calendar_precise()
            except AttributeError:
                calc_daily_calendar.calculate_calendar()
    else:
        print("\n[Skipping Daily Calendar generation]")

    # ---------------------------------------------------------
    # 3. MERGE EVENTS INTO MASTER CSV
    # ---------------------------------------------------------
    print("\n------------------------------------------")
    print("Merging discrete events into 'astrocalc.csv'...")

    # List of files to merge
    event_files = [
        "seasons.csv",
        "earth_extrema.csv",
        "moon_phases_detailed.csv",
        "eclipses.csv",
        "inner_planets.csv",
        "outer_planets.csv",
        "conjunctions.csv"
    ]

    master_events = []

    for filename in event_files:
        if os.path.exists(filename):
            with open(filename, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                if 'datetime' in reader.fieldnames:
                    for row in reader:
                        master_events.append(row)
                else:
                    print(f"  Warning: Skipping {filename} (missing 'datetime' column)")
        else:
            print(f"  Warning: {filename} not found (script failed?)")

    # Sort by datetime
    master_events.sort(key=lambda x: x['datetime'])

    # Write Master File
    output_file = "astrocalc.csv"
    headers = ["datetime", "event", "details"]

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        
        count = 0
        for row in master_events:
            cleaned_row = {
                "datetime": row.get("datetime", ""),
                "event": row.get("event", ""),
                "details": row.get("details", "")
            }
            writer.writerow(cleaned_row)
            count += 1

    print("------------------------------------------")
    print(f"SUCCESS! {count} events merged into '{output_file}'.")
    if response in ['y', 'yes']:
        print(f"Daily Ephemeris Grid saved to 'daily_calendar.csv'.")
    print("==========================================")

if __name__ == "__main__":
    main()