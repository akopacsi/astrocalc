# Astrocalc 🌌

A professional-grade astronomical almanac generator written in Python. This tool calculates precise astronomical events, planetary phenomena, and generates detailed daily ephemeris calendars.

It was designed to match the precision of printed astronomical almanacs, taking into account observer location, timezone, and visual magnitude.

## 🚀 Features & Calculation Criteria

### 1. Master Event List (`astrocalc.csv`)
Generates a chronological list of major astronomical phenomena for a time range.

#### 🌍 Earth & Seasons
* **Apsides:** Calculates **Perihelion** (closest to Sun) and **Aphelion** (farthest) using precise velocity-vector flip detection.
* **Seasons:** Calculates precise times for March/September Equinoxes and June/December Solstices using geocentric apparent ecliptic longitude.

#### 🌙 Moon Phenomena
* **Phases:** Standard New Moon, First Quarter, Full Moon, and Last Quarter.
* **Supermoons:** Detects Full Moons occurring when the Moon is closer than **363,000 km** to Earth.
* **Crescent Visibility:**
    * Calculates the first visibility of the "Young Moon" (Evening) and last visibility of the "Old Moon" (Morning).
    * **Criteria:** The Moon must be above the horizon (Altitude > 0°) at the moment the Sun reaches **Civil Twilight (Altitude -6°)**.
    * Reports the crescent's Age (days), Altitude, and Azimuth.

#### 🌑 Eclipses
* **Solar Eclipses:**
    * **Global Criteria:** Moon center passes within **1.6°** of Sun center.
    * **Local Visibility:** Checks if the Moon's disk overlaps the Sun's disk at the specific observer location.
    * **Types:** Total (Mag ≥ 1.0), Annular, and Partial.
* **Lunar Eclipses:**
    * **Criteria:** Moon passes within **~1.0°** of the anti-solar point.
    * **Umbra Check:** Detects entry into Earth's umbral shadow radius (~0.9°).
    * **Visibility:** Filters out events where the Moon is below the horizon at maximum eclipse.
    * **Types:** Total and Partial.

#### 🪐 Planetary Events
* **Inner Planets (Mercury & Venus):**
    * **Extrema:** Perigee (Closest) and Apogee (Farthest).
    * **Elongations:** Greatest Eastern (Evening visibility) and Western (Morning visibility) Elongations.
    * **Conjunctions:** Inferior (between Earth/Sun) and Superior (behind Sun).
    * **Transits of the Sun:**
        * **Criteria:** Separation < **0.27°** during Inferior Conjunction.
        * **Visibility:** Checks if Sun is above horizon during transit.
* **Outer Planets (Mars, Jupiter, Saturn, Uranus, Neptune):**
    * **Oppositions:** When the planet is opposite the Sun (best visibility).
    * **Solar Conjunctions:** When the planet is behind the Sun (unobservable).
    * **Extrema:** Perigee and Apogee distances.
* **Visual Data:** All planetary events include estimated **Magnitude** (brightness) and **Apparent Diameter** (in arcseconds).

#### 💫 Conjunctions
* **Moon & Planets:** Detects close approaches where separation is **< 6.0°**.
* **Planet & Planet:** Detects close approaches where separation is **< 3.0°**.

### 2. Daily Calendar Grid (`daily_calendar.csv`)
Generates a dense daily data grid ideal for printing or planning observations.

* **Sun:** Sunrise, Sunset, **Upper Meridian Transit (Noon)**, Altitude at Transit, and Equation of Time.
* **Moon:** Moonrise, Transit, Moonset, and daily Phase status.
* **Time Standards:**
    * **Julian Date:** Terrestrial Time (TT) at 12:00 UTC.
    * **Sidereal Time:** Greenwich Apparent Sidereal Time (GAST) at 00:00 UTC.

## 🛠️ Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/akopacsi/astrocalc.git](https://github.com/akopacsi/astrocalc.git)
    cd astrocalc
    ```

2.  **Install dependencies:**
    This project relies on the high-precision [Skyfield](https://rhodesmill.org/skyfield/) library.
    ```bash
    pip install skyfield numpy
    ```

3.  **First Run:**
    The first time you run the script, it will automatically download the required NASA JPL ephemeris file (`de440.bsp`, ~114 MB).

## ⚙️ Configuration

Open `config.json` to set your observer location and time range.

Example for Hungary (CET):
```json
{
    "location": {
        "latitude": 47.5,
        "longitude": 19.0,
        "elevation_m": 0
    },
    "time_range": {
        "start": "2024-01-01",
        "end": "2024-12-31"
    },
    "timezone_offset": 1.0,
    "ephemeris_file": "de440.bsp"
}
```

## Copyright (c) 2026 Antal Kopácsi

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
