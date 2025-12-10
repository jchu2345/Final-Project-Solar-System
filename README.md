# Final-Project-Solar-System
NDHU python physics course Final project

## Setting up a Virtual Environment

It is recommended to use a virtual environment to manage project dependencies. Follow these steps to set up and activate your virtual environment:

1.  **Create the virtual environment:**
    ```bash
    python -m venv .venv
    ```

2.  **Activate the virtual environment:**

    *   **On Windows:**
        ```bash
        .venv\Scripts\activate
        ```
    *   **On macOS/Linux:**
        ```bash
        source .venv/bin/activate
        ```

After activation, your terminal prompt should show `(.venv)` indicating that the virtual environment is active.

## Solar System Simulation

This project contains a 3D simulation of the solar system, built using VPython. The simulation uses real-world physical data for celestial bodies to model their orbits.

### Features

*   **Real-world data:** The simulation uses real physical data for planets, moons, and the sun, including mass, radius, and orbital parameters.
*   **All 8 planets + Pluto:** The simulation includes all 8 planets in our solar system, plus the dwarf planet Pluto.
*   **Earth's Moon:** The simulation includes the Earth's Moon, orbiting the Earth.
*   **Saturn's Rings:** Saturn is rendered with its iconic rings.
*   **3D Visualization:** The simulation is visualized in 3D using VPython, with planets leaving trails to show their orbital paths.
*   **Realistic Orbits:** The simulation accounts for orbital eccentricity and inclination, resulting in more realistic elliptical and tilted orbits.

### Running the Simulation

1.  **Activate your virtual environment** (see "Setting up a Virtual Environment" above).

2.  **Install VPython:**
    ```bash
    pip install vpython
    ```

3.  **Run the simulation:**
    ```bash
    python main.py
    ```

The VPython window will open and the simulation will start.
