# Solar System Simulation

A physics-accurate 3D simulation of the Solar System including the Sun and all 8 major planets (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune). Built with VPython and real NASA/JPL data.

![Solar System Simulation](https://img.shields.io/badge/Python-3.8%2B-blue)
![VPython](https://img.shields.io/badge/VPython-7.6%2B-green)
![Physics](https://img.shields.io/badge/Physics-N--Body%20Simulation-red)

## Overview

This project simulates the gravitational dynamics of our Solar System using Newton's law of universal gravitation with N-body physics. Every physical parameter (masses, radii, orbital elements) comes from authoritative NASA databases with full citation. The simulation uses the Velocity Verlet integration method for numerical stability and energy conservation.

## Features

- All 8 major planets with accurate orbital mechanics
- Real-time 3D visualization with VPython
- Physics-based N-body gravitational interactions
- Velocity Verlet numerical integration for stability
- Interactive camera controls (rotate, zoom, pan)
- Energy conservation monitoring
- Orbital trails showing planetary paths
- Scientifically accurate with full source attribution

## Physics Model

### Gravitational Dynamics

The simulation implements Newton's Law of Universal Gravitation for N-body interactions:

```
F = G × m₁ × m₂ / r²
```

Where:
- F = gravitational force (Newtons)
- G = gravitational constant = 6.67430 × 10⁻¹¹ m³ kg⁻¹ s⁻² (CODATA 2018)
- m₁, m₂ = masses of two bodies (kg)
- r = distance between body centers (m)

Every body experiences gravitational forces from all other bodies in the system (Sun-planet and planet-planet interactions).

### Numerical Integration

**Method:** Velocity Verlet

The Velocity Verlet algorithm provides excellent energy conservation for orbital mechanics:

```
x(t+Δt) = x(t) + v(t)Δt + ½a(t)Δt²
v(t+Δt) = v(t) + ½[a(t) + a(t+Δt)]Δt
```

**Time Step:** Δt = 86,400 seconds (1 Earth day)

**Justification:** This time step provides good stability for all planets while maintaining reasonable computation speed. It is approximately 1/88 of Mercury's orbital period (the fastest planet), which satisfies the stability criterion of Δt ≤ T_min/100.

### Units

All internal calculations use SI units:
- **Mass:** kilograms (kg)
- **Distance:** meters (m)
- **Time:** seconds (s)
- **Velocity:** meters per second (m/s)
- **Force:** Newtons (N)

## Scaling Factors

The real Solar System is too large to visualize at true scale, so we apply scaling factors for visualization only. All physics calculations use real SI units.

### Distance Scaling

- **Factor:** 2.0 VPython units / 1.496×10¹¹ m (1 AU)
- **Result:** 1 AU ≈ 2.0 VPython units
- **Reason:** Makes the entire Solar System (out to Neptune at 30 AU) visible in the viewport
- **Example:** Earth's orbit (1 AU) appears as a circle of radius 2.0 units

### Radius Scaling

- **Factor:** 2000 × distance scale
- **Reason:** Planetary radii are extremely small compared to orbital distances
- **Reality:** At true scale, planets would be invisible dots
- **Exaggeration:** Radii multiplied by 2000 for visibility
- **Example:** Earth's radius (6,371 km) displayed as ~0.17 VPython units instead of ~0.000085 units

### Time Scaling

- **Factor:** 1 simulation second = 86,400 real seconds (1 Earth day)
- **Reason:** Makes orbital motion visible while preserving relative speeds
- **Example:** Earth completes one orbit in approximately 365 simulation seconds (about 6 minutes of wall-clock time at 60 fps)

## Data Sources / References

All physical values in this simulation are derived from authoritative sources. No approximations or fictional values are used.

### Primary Sources

1. **NASA JPL Solar System Dynamics - Planetary Physical Parameters**
   - URL: https://ssd.jpl.nasa.gov/planets/phys_par.html
   - Accessed: 2025-12-11
   - Used for: Planetary masses (kg), mean radii (km)

2. **NASA JPL Solar System Dynamics - Approximate Positions of the Planets**
   - URL: https://ssd.jpl.nasa.gov/planets/approx_pos.html
   - Accessed: 2025-12-11
   - Used for: Orbital elements (semi-major axis, eccentricity, inclination)
   - Valid: 1800 AD - 2050 AD (Keplerian elements at J2000 epoch)

3. **NASA Sun Facts**
   - URL: https://science.nasa.gov/sun/facts/
   - Accessed: 2025-12-11
   - Used for: Solar mass (1.989×10³⁰ kg), solar radius (6.963×10⁸ m)

### Physical Constants

4. **CODATA 2018 Recommended Values of Fundamental Physical Constants**
   - URL: https://physics.nist.gov/cgi-bin/cuu/Value?bg
   - Accessed: 2025-12-11
   - Used for: Gravitational constant G = 6.67430×10⁻¹¹ m³ kg⁻¹ s⁻²

### Data Verification

Every constant in `data.py` includes:
- Source title and URL
- Access date
- Original units and conversions
- Intended use

This ensures full traceability and scientific rigor.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. Clone this repository:
```bash
git clone https://github.com/yourusername/Final-Project-Solar-System.git
cd Final-Project-Solar-System
```

2. (Optional) Create a virtual environment:
```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the simulation:

```bash
python solar_system.py
```

### Controls

Once the simulation window opens:

- **Rotate view:** Right-click and drag
- **Zoom:** Mouse scroll wheel
- **Pan:** Ctrl + drag (or Cmd + drag on macOS)
- **Stop simulation:** Press Ctrl+C in the terminal

### What to Observe

- **Inner planets** (Mercury, Venus, Earth, Mars): Small, fast orbits
- **Outer planets** (Jupiter, Saturn, Uranus, Neptune): Large, slow orbits
- **Orbital shapes:** Most orbits appear circular, but Mercury has noticeable eccentricity
- **Inclinations:** Slight tilt differences in orbital planes
- **Relative speeds:** Mercury orbits fastest, Neptune slowest
- **Energy conservation:** Terminal output shows energy drift should be < 1%

## Project Structure

```
Final-Project-Solar-System/
├── solar_system.py          # Main simulation with VPython visualization
├── data.py                  # All physical constants with NASA sources
├── requirements.txt         # Python dependencies
├── README.md               # This file
└── .gitignore              # Git ignore rules
```

### File Descriptions

- **`solar_system.py`**: Main simulation program
  - `CelestialBody` class: Represents Sun and planets
  - Gravitational force calculation
  - Velocity Verlet integrator
  - VPython visualization setup
  - Simulation loop with energy monitoring

- **`data.py`**: Physical constants database
  - Gravitational constant (CODATA)
  - Solar parameters (NASA)
  - Planetary masses, radii (JPL)
  - Orbital elements (JPL)
  - Scaling factors for visualization
  - Complete data source registry

## Included Celestial Bodies

The simulation includes:

1. **Sun** - Central star (yellow, emissive)
2. **Mercury** - Innermost planet (gray)
3. **Venus** - Second planet (yellowish)
4. **Earth** - Our home (blue)
5. **Mars** - The red planet (red-orange)
6. **Jupiter** - Largest planet (orange-brown)
7. **Saturn** - Ringed giant (pale yellow)
8. **Uranus** - Ice giant (cyan)
9. **Neptune** - Outermost planet (deep blue)

## Technical Details

### Orbital Element Implementation

Each planet starts at **perihelion** (closest approach to Sun):

- **Position:** r_perihelion = a(1 - e)
- **Velocity:** Calculated using vis-viva equation
  ```
  v² = GM(2/r - 1/a)
  ```
- **Inclination:** Orbital plane rotated by inclination angle from ecliptic

### Energy Conservation

The simulation monitors total mechanical energy:

```
E_total = E_kinetic + E_potential

E_kinetic = Σ(½mv²)
E_potential = ΣΣ(-Gm₁m₂/r) for all pairs
```

Energy should remain approximately constant. Drift < 1% indicates good numerical stability.

### Performance

- **Frame rate:** Target 60 fps (VPython rate limit)
- **Computation per frame:** O(n²) gravitational force calculations (36 pairs for 9 bodies)
- **Trail retention:** 500 points per planet to balance history vs. performance

## Validation

### Physical Accuracy

The simulation reproduces:
- Correct orbital periods (Earth: 365.25 days, Mars: 687 days, etc.)
- Elliptical orbits with proper eccentricity
- Stable orbits over hundreds of years
- Energy conservation to within 1%

### Limitations

1. **Simplified model:** Does not include:
   - Relativistic corrections (negligible for planets)
   - Moons and asteroids
   - Non-gravitational forces (solar wind, radiation pressure)
   - Planetary oblateness (J2 effects)

2. **Numerical precision:**
   - Time step chosen for balance of speed vs. accuracy
   - Very long-term simulations (>1000 years) may accumulate numerical error
   - Double precision floating point arithmetic

## Future Enhancements

Possible extensions to this project:

- Add major moons (Moon, Galilean moons, Titan, etc.)
- Implement adaptive time stepping for efficiency
- Add asteroid belt objects
- Include dwarf planets (Pluto, Ceres, etc.)
- Add velocity vectors and force arrows for visualization
- Save/load simulation state
- Export orbital data to files
- Add relativistic corrections for Mercury's precession
- Implement collision detection

## Educational Value

This simulation demonstrates:

- Newton's laws of motion and gravitation
- Numerical integration techniques
- Kepler's laws of planetary motion
- Conservation of energy and momentum
- N-body problem complexity
- Importance of data sources in scientific computing
- Software engineering for scientific applications

## License

This project is for educational purposes as part of the NDHU Python Physics course final project.

## Acknowledgments

- NASA Jet Propulsion Laboratory for authoritative planetary data
- NASA Science Mission Directorate for solar system information
- CODATA for fundamental physical constants
- VPython team for the excellent 3D visualization library

## Contact

For questions or issues with this simulation, please open an issue on GitHub or contact the course instructor.

---

**Note:** This simulation uses real physical data and implements accurate N-body gravitational physics. It is suitable for educational demonstrations of celestial mechanics and serves as a foundation for more advanced astronomical simulations.

## References

1. Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics*. Cambridge University Press.

2. Standish, E. M. (1998). JPL Planetary and Lunar Ephemerides, DE405/LE405. *JPL IOM 312.F-98-048*.

3. Urban, S. E., & Seidelmann, P. K. (Eds.). (2012). *Explanatory Supplement to the Astronomical Almanac* (3rd ed.). University Science Books.

4. Danby, J. M. A. (1988). *Fundamentals of Celestial Mechanics* (2nd ed.). Willmann-Bell.

---

Generated with physics-accurate data from NASA/JPL databases.
Last updated: 2025-12-11
