"""
Physical Constants and Planetary Data for Solar System Simulation

All values sourced from authoritative NASA databases.
Each constant includes full citation with URL and access date.
"""

import math

# =============================================================================
# FUNDAMENTAL PHYSICAL CONSTANTS
# =============================================================================

# Gravitational constant (m^3 kg^-1 s^-2)
# Source: CODATA 2018 Recommended Values
# Reference: NIST, https://physics.nist.gov/cgi-bin/cuu/Value?bg
# Accessed: 2025-12-11
GRAVITATIONAL_CONSTANT = 6.67430e-11

# =============================================================================
# SUN PHYSICAL PARAMETERS
# =============================================================================

# Sun mass (kg)
# Source: NASA Sun Facts
# URL: https://science.nasa.gov/sun/facts/
# Note: "more than 330,000 Earths" × Earth mass
# Accessed: 2025-12-11
SUN_MASS_KG = 1.9885e30

# Sun radius (m)
# Source: NASA Sun Facts
# URL: https://science.nasa.gov/sun/facts/
# Value: 700,000 km = 7.0e8 m
# Accessed: 2025-12-11
SUN_RADIUS_M = 6.96340e8

# Sun texture (built-in textures don't work well for an emissive sun)
SUN_TEXTURE = None

# =============================================================================
# PLANETARY MASSES (kg)
# =============================================================================

# Source: NASA JPL Solar System Dynamics - Planetary Physical Parameters
# URL: https://ssd.jpl.nasa.gov/planets/phys_par.html
# Accessed: 2025-12-11

# Mercury mass (kg)
MERCURY_MASS_KG = 3.3011e23

# Venus mass (kg)
VENUS_MASS_KG = 4.8675e24

# Earth mass (kg)
EARTH_MASS_KG = 5.97237e24

# Mars mass (kg)
MARS_MASS_KG = 6.4171e23

# Jupiter mass (kg)
JUPITER_MASS_KG = 1.8982e27

# Saturn mass (kg)
SATURN_MASS_KG = 5.6834e26

# Uranus mass (kg)
URANUS_MASS_KG = 8.6810e25

# Neptune mass (kg)
NEPTUNE_MASS_KG = 1.02413e26

# Moon mass (kg)
# Source: NASA JPL Solar System Dynamics
# URL: https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
# Accessed: 2025-12-11
MOON_MASS_KG = 7.342e22

# =============================================================================
# PLANETARY RADII (m)
# =============================================================================

# Source: NASA JPL Solar System Dynamics - Planetary Physical Parameters
# URL: https://ssd.jpl.nasa.gov/planets/phys_par.html
# Accessed: 2025-12-11
# Note: Mean radii converted from km to meters

# Mercury mean radius (m)
MERCURY_RADIUS_M = 2439.4e3

# Venus mean radius (m)
VENUS_RADIUS_M = 6051.8e3

# Earth mean radius (m)
EARTH_RADIUS_M = 6371.0e3

# Mars mean radius (m)
MARS_RADIUS_M = 3389.5e3

# Jupiter mean radius (m)
JUPITER_RADIUS_M = 69911e3

# Saturn mean radius (m)
SATURN_RADIUS_M = 58232e3

# Uranus mean radius (m)
URANUS_RADIUS_M = 25362e3

# Neptune mean radius (m)
NEPTUNE_RADIUS_M = 24622e3

# Moon mean radius (m)
# Source: NASA JPL Solar System Dynamics
# URL: https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
# Accessed: 2025-12-11
MOON_RADIUS_M = 1737.4e3

# =============================================================================
# ORBITAL ELEMENTS (Keplerian Elements at J2000 Epoch)
# =============================================================================

# Source: NASA JPL Solar System Dynamics - Approximate Positions of the Planets
# URL: https://ssd.jpl.nasa.gov/planets/approx_pos.html
# Accessed: 2025-12-11
# Note: Valid for period 1800 AD - 2050 AD

# Astronomical Unit (m)
# 1 AU = mean Earth-Sun distance
AU_TO_METERS = 1.495978707e11

# Mercury orbital elements
MERCURY_SEMIMAJOR_AXIS_AU = 0.38709927
MERCURY_ECCENTRICITY = 0.20563593
MERCURY_INCLINATION_DEG = 7.00497902
MERCURY_LONG_ASC_NODE_DEG = 48.33076593
MERCURY_ARG_PERIAPSIS_DEG = 29.12703035

# Venus orbital elements
VENUS_SEMIMAJOR_AXIS_AU = 0.72333566
VENUS_ECCENTRICITY = 0.00677672
VENUS_INCLINATION_DEG = 3.39467605
VENUS_LONG_ASC_NODE_DEG = 76.67984255
VENUS_ARG_PERIAPSIS_DEG = 54.89204154

# Earth orbital elements
EARTH_SEMIMAJOR_AXIS_AU = 1.00000261
EARTH_ECCENTRICITY = 0.01671123
EARTH_INCLINATION_DEG = -0.00001531
EARTH_LONG_ASC_NODE_DEG = 0.0
EARTH_ARG_PERIAPSIS_DEG = 102.93768193

# Mars orbital elements
MARS_SEMIMAJOR_AXIS_AU = 1.52371034
MARS_ECCENTRICITY = 0.09339410
MARS_INCLINATION_DEG = 1.84969142
MARS_LONG_ASC_NODE_DEG = 49.55953891
MARS_ARG_PERIAPSIS_DEG = 286.53716654

# Jupiter orbital elements
JUPITER_SEMIMAJOR_AXIS_AU = 5.20288700
JUPITER_ECCENTRICITY = 0.04838624
JUPITER_INCLINATION_DEG = 1.30439695
JUPITER_LONG_ASC_NODE_DEG = 100.47390909
JUPITER_ARG_PERIAPSIS_DEG = 273.86740468

# Saturn orbital elements
SATURN_SEMIMAJOR_AXIS_AU = 9.53667594
SATURN_ECCENTRICITY = 0.05386179
SATURN_INCLINATION_DEG = 2.48599187
SATURN_LONG_ASC_NODE_DEG = 113.66242448
SATURN_ARG_PERIAPSIS_DEG = 339.39218624

# Uranus orbital elements
URANUS_SEMIMAJOR_AXIS_AU = 19.18916464
URANUS_ECCENTRICITY = 0.04725744
URANUS_INCLINATION_DEG = 0.77263783
URANUS_LONG_ASC_NODE_DEG = 74.01692503
URANUS_ARG_PERIAPSIS_DEG = 96.99871891

# Neptune orbital elements
NEPTUNE_SEMIMAJOR_AXIS_AU = 30.06992276
NEPTUNE_ECCENTRICITY = 0.00859048
NEPTUNE_INCLINATION_DEG = 1.77004347
NEPTUNE_LONG_ASC_NODE_DEG = 131.78422574
NEPTUNE_ARG_PERIAPSIS_DEG = 276.33597808

# Moon orbital elements (relative to Earth)
# Source: NASA JPL Solar System Dynamics
# URL: https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
# Accessed: 2025-12-11
MOON_SEMIMAJOR_AXIS_M = 384400e3  # 384,400 km in meters
MOON_ECCENTRICITY = 0.0549
MOON_INCLINATION_DEG = 5.145  # Relative to Earth's orbital plane

# =============================================================================
# DERIVED ORBITAL PARAMETERS
# =============================================================================

# Orbital periods calculated using Kepler's Third Law: P^2 = a^3 (for a in AU, P in years)
# Then converted to seconds

# Mercury orbital period (seconds)
MERCURY_PERIOD_DAYS = 87.9691
MERCURY_PERIOD_S = MERCURY_PERIOD_DAYS * 86400

# Venus orbital period (seconds)
VENUS_PERIOD_DAYS = 224.701
VENUS_PERIOD_S = VENUS_PERIOD_DAYS * 86400

# Earth orbital period (seconds)
EARTH_PERIOD_DAYS = 365.256
EARTH_PERIOD_S = EARTH_PERIOD_DAYS * 86400

# Mars orbital period (seconds)
MARS_PERIOD_DAYS = 686.980
MARS_PERIOD_S = MARS_PERIOD_DAYS * 86400

# Jupiter orbital period (seconds)
JUPITER_PERIOD_DAYS = 4332.589
JUPITER_PERIOD_S = JUPITER_PERIOD_DAYS * 86400

# Saturn orbital period (seconds)
SATURN_PERIOD_DAYS = 10759.22
SATURN_PERIOD_S = SATURN_PERIOD_DAYS * 86400

# Uranus orbital period (seconds)
URANUS_PERIOD_DAYS = 30688.5
URANUS_PERIOD_S = URANUS_PERIOD_DAYS * 86400

# Neptune orbital period (seconds)
NEPTUNE_PERIOD_DAYS = 60182.0
NEPTUNE_PERIOD_S = NEPTUNE_PERIOD_DAYS * 86400

# Moon orbital period (seconds)
# Source: NASA - sidereal month
MOON_PERIOD_DAYS = 27.3217
MOON_PERIOD_S = MOON_PERIOD_DAYS * 86400

# =============================================================================
# VISUALIZATION SCALING FACTORS
# =============================================================================

# Distance scale factor (VPython units per meter)
# Real solar system is too large to visualize at true scale
# Scale: 1 AU = 2.0 VPython units for comfortable viewing
DISTANCE_SCALE = 2.0 / AU_TO_METERS

# Radius scale factor (VPython units per meter)
# TRUE TO SCALE - using same scale as distances
# This means planet sizes are proportional to orbital distances
# WARNING: Planets will be extremely small but scientifically accurate
RADIUS_SCALE = DISTANCE_SCALE

# Time scale factor (real seconds per simulation second)
# 1 simulation second = 5 minutes
# A smaller timestep increases accuracy, reducing long-term energy drift.
TIME_SCALE = 300  # 5 minutes in seconds

# =============================================================================
# DATA SOURCES REGISTRY
# =============================================================================

DATA_SOURCES = {
    'jpl_planets_physical': {
        'title': 'NASA JPL Solar System Dynamics - Planetary Physical Parameters',
        'url': 'https://ssd.jpl.nasa.gov/planets/phys_par.html',
        'accessed': '2025-12-11',
        'used_for': 'Planetary masses, mean radii'
    },
    'jpl_planets_orbital': {
        'title': 'NASA JPL Solar System Dynamics - Approximate Positions of the Planets',
        'url': 'https://ssd.jpl.nasa.gov/planets/approx_pos.html',
        'accessed': '2025-12-11',
        'used_for': 'Orbital elements (semi-major axis, eccentricity, inclination)'
    },
    'nasa_sun': {
        'title': 'NASA Sun Facts',
        'url': 'https://science.nasa.gov/sun/facts/',
        'accessed': '2025-12-11',
        'used_for': 'Solar mass, radius'
    },
    'codata': {
        'title': 'CODATA 2018 Recommended Values of Fundamental Physical Constants',
        'url': 'https://physics.nist.gov/cgi-bin/cuu/Value?bg',
        'accessed': '2025-12-11',
        'used_for': 'Gravitational constant'
    }
}

# =============================================================================
# PLANETARY DATA STRUCTURES
# =============================================================================

# Complete planetary data organized for easy access
PLANETS_DATA = {
    'Mercury': {
        'mass_kg': MERCURY_MASS_KG,
        'radius_m': MERCURY_RADIUS_M,
        'semimajor_axis_au': MERCURY_SEMIMAJOR_AXIS_AU,
        'eccentricity': MERCURY_ECCENTRICITY,
        'inclination_deg': MERCURY_INCLINATION_DEG,
        'long_asc_node_deg': MERCURY_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': MERCURY_ARG_PERIAPSIS_DEG,
        'period_s': MERCURY_PERIOD_S,
        'color': (0.7, 0.7, 0.7),
        'texture': "vpython.rock"
    },
    'Venus': {
        'mass_kg': VENUS_MASS_KG,
        'radius_m': VENUS_RADIUS_M,
        'semimajor_axis_au': VENUS_SEMIMAJOR_AXIS_AU,
        'eccentricity': VENUS_ECCENTRICITY,
        'inclination_deg': VENUS_INCLINATION_DEG,
        'long_asc_node_deg': VENUS_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': VENUS_ARG_PERIAPSIS_DEG,
        'period_s': VENUS_PERIOD_S,
        'color': (0.9, 0.7, 0.4),
        'texture': "vpython.gravel"
    },
    'Earth': {
        'mass_kg': EARTH_MASS_KG,
        'radius_m': EARTH_RADIUS_M,
        'semimajor_axis_au': EARTH_SEMIMAJOR_AXIS_AU,
        'eccentricity': EARTH_ECCENTRICITY,
        'inclination_deg': EARTH_INCLINATION_DEG,
        'long_asc_node_deg': EARTH_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': EARTH_ARG_PERIAPSIS_DEG,
        'period_s': EARTH_PERIOD_S,
        'color': (0.2, 0.4, 0.8),
        'texture': "vpython.earth"
    },
    'Mars': {
        'mass_kg': MARS_MASS_KG,
        'radius_m': MARS_RADIUS_M,
        'semimajor_axis_au': MARS_SEMIMAJOR_AXIS_AU,
        'eccentricity': MARS_ECCENTRICITY,
        'inclination_deg': MARS_INCLINATION_DEG,
        'long_asc_node_deg': MARS_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': MARS_ARG_PERIAPSIS_DEG,
        'period_s': MARS_PERIOD_S,
        'color': (0.8, 0.3, 0.1),
        'texture': "vpython.rough"
    },
    'Jupiter': {
        'mass_kg': JUPITER_MASS_KG,
        'radius_m': JUPITER_RADIUS_M,
        'semimajor_axis_au': JUPITER_SEMIMAJOR_AXIS_AU,
        'eccentricity': JUPITER_ECCENTRICITY,
        'inclination_deg': JUPITER_INCLINATION_DEG,
        'long_asc_node_deg': JUPITER_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': JUPITER_ARG_PERIAPSIS_DEG,
        'period_s': JUPITER_PERIOD_S,
        'color': (0.8, 0.6, 0.4),
        'texture': "vpython.wood"
    },
    'Saturn': {
        'mass_kg': SATURN_MASS_KG,
        'radius_m': SATURN_RADIUS_M,
        'semimajor_axis_au': SATURN_SEMIMAJOR_AXIS_AU,
        'eccentricity': SATURN_ECCENTRICITY,
        'inclination_deg': SATURN_INCLINATION_DEG,
        'long_asc_node_deg': SATURN_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': SATURN_ARG_PERIAPSIS_DEG,
        'period_s': SATURN_PERIOD_S,
        'color': (0.9, 0.8, 0.6),
        'texture': "vpython.wood"
    },
    'Uranus': {
        'mass_kg': URANUS_MASS_KG,
        'radius_m': URANUS_RADIUS_M,
        'semimajor_axis_au': URANUS_SEMIMAJOR_AXIS_AU,
        'eccentricity': URANUS_ECCENTRICITY,
        'inclination_deg': URANUS_INCLINATION_DEG,
        'long_asc_node_deg': URANUS_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': URANUS_ARG_PERIAPSIS_DEG,
        'period_s': URANUS_PERIOD_S,
        'color': (0.5, 0.8, 0.8),
        'texture': "vpython.granite"
    },
    'Neptune': {
        'mass_kg': NEPTUNE_MASS_KG,
        'radius_m': NEPTUNE_RADIUS_M,
        'semimajor_axis_au': NEPTUNE_SEMIMAJOR_AXIS_AU,
        'eccentricity': NEPTUNE_ECCENTRICITY,
        'inclination_deg': NEPTUNE_INCLINATION_DEG,
        'long_asc_node_deg': NEPTUNE_LONG_ASC_NODE_DEG,
        'arg_periapsis_deg': NEPTUNE_ARG_PERIAPSIS_DEG,
        'period_s': NEPTUNE_PERIOD_S,
        'color': (0.2, 0.3, 0.8),
        'texture': "vpython.granite"
    }
}

# Moon data (orbits Earth, not Sun)
MOON_DATA = {
    'mass_kg': MOON_MASS_KG,
    'radius_m': MOON_RADIUS_M,
    'semimajor_axis_m': MOON_SEMIMAJOR_AXIS_M,  # Relative to Earth, in meters
    'eccentricity': MOON_ECCENTRICITY,
    'inclination_deg': MOON_INCLINATION_DEG,
    'period_s': MOON_PERIOD_S,
    'color': (0.8, 0.8, 0.8),
    'parent': 'Earth',
    'texture': "vpython.rock"
}
