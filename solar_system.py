"""
Solar System Simulation with VPython

A physics-accurate N-body simulation of the Sun and all 8 major planets.
Uses real NASA data for all physical parameters and orbital elements.
Implements Velocity Verlet integration for numerical stability.

Author: Solar System Simulation Project
Date: 2025-12-11
"""

from vpython import (
    vector, sphere, canvas, rate, color, mag, norm, keysdown, button, wtext, curve,
    graph, gcurve, textures
)
import math
from data import (
    GRAVITATIONAL_CONSTANT, SUN_MASS_KG, SUN_RADIUS_M, SUN_TEXTURE,
    PLANETS_DATA, MOON_DATA, AU_TO_METERS, DISTANCE_SCALE, RADIUS_SCALE, TIME_SCALE
)


class CelestialBody:
    """
    Represents a celestial body (Sun or planet) with mass, position, velocity.

    Attributes:
        name: Body name (e.g., "Earth")
        mass: Mass in kilograms
        position: Position vector in meters (VPython vector)
        velocity: Velocity vector in m/s (VPython vector)
        sphere: VPython sphere object for visualization
        trail: List of previous positions for orbit visualization
        parent_body: The body this object orbits (e.g., Sun for Earth)
        initial_specific_angular_momentum: h = |r x v| at start
        initial_semimajor_axis: Orbital semi-major axis 'a' at start
        error_label: VPython wtext object for displaying error metrics
    """

    def __init__(self, name, mass_kg, radius_m, position_m, velocity_m_s,
                 body_color, texture=None, make_trail=True):
        """
        Initialize a celestial body.

        Args:
            name: Name of the body
            mass_kg: Mass in kilograms
            radius_m: Physical radius in meters
            position_m: Initial position vector in meters
            velocity_m_s: Initial velocity vector in m/s
            body_color: RGB tuple for body color
            texture: URL or path to a texture file
            make_trail: Whether to show orbital trail
        """
        self.name = name
        self.mass = mass_kg
        self.position = position_m  # meters
        self.velocity = velocity_m_s  # m/s
        self.acceleration = vector(0, 0, 0)  # m/s^2

        # Physics validation attributes
        self.parent_body = None
        self.initial_specific_angular_momentum = None
        self.initial_semimajor_axis = None
        self.error_label = None

        # A mapping from string names to VPython texture objects
        vpython_textures = {
            "vpython.earth": textures.earth,
            "vpython.rock": textures.rock,
            "vpython.gravel": textures.gravel,
            "vpython.rough": textures.rough,
            "vpython.wood": textures.wood,
            "vpython.granite": textures.granite
        }

        # Determine texture
        body_texture = vpython_textures.get(texture)

        # Create VPython sphere for visualization
        self.sphere = sphere(
            pos=position_m * DISTANCE_SCALE,
            radius=radius_m * RADIUS_SCALE,
            color=vector(*body_color),
            make_trail=False,  # We'll use custom curve trails instead
            texture=body_texture
        )

        # Create orbital path curve if needed
        if make_trail:
            self.trail = curve(color=vector(*body_color))  # No radius = thin line
            self.trail_positions = []
            self.trail_interval = 10  # Add point every N frames
            self.trail_counter = 0
        else:
            self.trail = None

        # Special rendering for the Sun
        if name == "Sun":
            self.sphere.emissive = True
            # For emissive objects, color should be white to show texture properly
            if body_texture:
                 self.sphere.color = color.white

    def update_position(self, dt):
        """
        Update position using Velocity Verlet algorithm (first step).

        x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2

        Args:
            dt: Time step in seconds
        """
        self.position = self.position + self.velocity * dt + 0.5 * self.acceleration * dt**2
        self.sphere.pos = self.position * DISTANCE_SCALE

        # Update trail
        if self.trail is not None:
            self.trail_counter += 1
            if self.trail_counter >= self.trail_interval:
                self.trail.append(pos=self.sphere.pos)
                self.trail_counter = 0

    def update_velocity(self, new_acceleration, dt):
        """
        Update velocity using Velocity Verlet algorithm (second step).

        v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt

        Args:
            new_acceleration: Acceleration at new time step (m/s^2)
            dt: Time step in seconds
        """
        self.velocity = self.velocity + 0.5 * (self.acceleration + new_acceleration) * dt
        self.acceleration = new_acceleration


def calculate_gravitational_force(body1, body2):
    """
    Calculate gravitational force on body1 due to body2 using Newton's law.

    F = G * m1 * m2 / r^2 (direction: from body1 toward body2)

    Args:
        body1: CelestialBody object (force is calculated on this body)
        body2: CelestialBody object (force is due to this body)

    Returns:
        vector: Gravitational force on body1 in Newtons

    Physics:
        Newton's Law of Universal Gravitation
        F = G * m1 * m2 / |r|^2 * r_hat
        where r_hat is the unit vector from body1 to body2
    """
    # Vector from body1 to body2
    r_vec = body2.position - body1.position
    r_magnitude = mag(r_vec)

    # Avoid division by zero if bodies are at same position
    if r_magnitude < 1e3:  # Less than 1 km
        return vector(0, 0, 0)

    # Calculate force magnitude: F = G * m1 * m2 / r^2
    force_magnitude = (GRAVITATIONAL_CONSTANT * body1.mass * body2.mass) / (r_magnitude**2)

    # Force direction: unit vector from body1 to body2
    force_direction = norm(r_vec)

    # Force vector
    force = force_magnitude * force_direction

    return force


def calculate_total_acceleration(body, all_bodies):
    """
    Calculate total gravitational acceleration on a body from all other bodies.

    a = sum(F_i / m) for all forces F_i

    Args:
        body: CelestialBody to calculate acceleration for
        all_bodies: List of all CelestialBody objects in system

    Returns:
        vector: Total acceleration in m/s^2
    """
    total_force = vector(0, 0, 0)

    for other_body in all_bodies:
        if other_body is not body:
            force = calculate_gravitational_force(body, other_body)
            total_force = total_force + force

    # a = F / m
    acceleration = total_force / body.mass

    return acceleration


def initialize_planet_position_velocity(semimajor_axis_au, eccentricity,
                                        inclination_deg, long_asc_node_deg,
                                        arg_periapsis_deg, sun_mass_kg):
    """
    Initialize planet position and velocity at perihelion using proper 3D orbital mechanics.

    Uses Keplerian orbital elements to correctly position planet in 3D space.

    Args:
        semimajor_axis_au: Semi-major axis in AU
        eccentricity: Orbital eccentricity (0 = circle, <1 = ellipse)
        inclination_deg: Orbital inclination in degrees (tilt of orbital plane)
        long_asc_node_deg: Longitude of ascending node (Ω) in degrees
        arg_periapsis_deg: Argument of periapsis (ω) in degrees
        sun_mass_kg: Mass of Sun in kg

    Returns:
        tuple: (position_vector_m, velocity_vector_m_s)

    Physics:
        - Perihelion distance: r_p = a(1 - e)
        - Vis-viva equation: v = sqrt(G*M*(2/r - 1/a))
        - 3D rotation: Apply Ω, i, ω rotations to transform from orbital plane to ecliptic
    """
    # Convert to SI units and radians
    semimajor_axis_m = semimajor_axis_au * AU_TO_METERS
    i = math.radians(inclination_deg)      # Inclination
    omega = math.radians(long_asc_node_deg)  # Longitude of ascending node (Ω)
    w = math.radians(arg_periapsis_deg)    # Argument of periapsis (ω)

    # Perihelion distance (closest approach)
    perihelion_distance = semimajor_axis_m * (1 - eccentricity)

    # Position at perihelion in orbital plane (periapsis direction = x-axis)
    x_orb = perihelion_distance
    y_orb = 0.0
    z_orb = 0.0

    # Velocity at perihelion using vis-viva equation
    v_magnitude = math.sqrt(
        GRAVITATIONAL_CONSTANT * sun_mass_kg *
        (2.0 / perihelion_distance - 1.0 / semimajor_axis_m)
    )

    # Velocity perpendicular to position in orbital plane (y-direction)
    vx_orb = 0.0
    vy_orb = v_magnitude
    vz_orb = 0.0

    # Transform from orbital plane to ecliptic plane using rotation matrices
    # Combined rotation: R = R_z(Ω) * R_x(i) * R_z(ω)

    # Precompute trig functions
    cos_omega, sin_omega = math.cos(omega), math.sin(omega)
    cos_i, sin_i = math.cos(i), math.sin(i)
    cos_w, sin_w = math.cos(w), math.sin(w)

    # Apply rotation matrix to position
    # First rotate by ω (argument of periapsis)
    x1 = x_orb * cos_w - y_orb * sin_w
    y1 = x_orb * sin_w + y_orb * cos_w
    z1 = z_orb

    # Then rotate by i (inclination)
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i

    # Finally rotate by Ω (longitude of ascending node)
    x_ecl = x2 * cos_omega - y2 * sin_omega
    y_ecl = x2 * sin_omega + y2 * cos_omega
    z_ecl = z2

    position = vector(x_ecl, y_ecl, z_ecl)

    # Apply same rotation to velocity
    vx1 = vx_orb * cos_w - vy_orb * sin_w
    vy1 = vx_orb * sin_w + vy_orb * cos_w
    vz1 = vz_orb

    vx2 = vx1
    vy2 = vy1 * cos_i - vz1 * sin_i
    vz2 = vy1 * sin_i + vz1 * cos_i

    vx_ecl = vx2 * cos_omega - vy2 * sin_omega
    vy_ecl = vx2 * sin_omega + vy2 * cos_omega
    vz_ecl = vz2

    velocity = vector(vx_ecl, vy_ecl, vz_ecl)

    return position, velocity


def calculate_total_energy(bodies):
    """
    Calculate total mechanical energy (kinetic + potential) of the system.

    Used for validation - should remain approximately constant for stable orbit.

    Args:
        bodies: List of all CelestialBody objects

    Returns:
        float: Total energy in Joules

    Physics:
        E_total = E_kinetic + E_potential
        E_kinetic = sum(0.5 * m * v^2)
        E_potential = sum(sum(-G * m1 * m2 / r)) for all pairs
    """
    kinetic_energy = 0
    potential_energy = 0

    # Kinetic energy
    for body in bodies:
        v_squared = mag(body.velocity)**2
        kinetic_energy += 0.5 * body.mass * v_squared

    # Potential energy (sum over all pairs)
    for i, body1 in enumerate(bodies):
        for body2 in bodies[i+1:]:
            r = mag(body2.position - body1.position)
            if r > 1e3:  # Avoid division by zero
                potential_energy -= (GRAVITATIONAL_CONSTANT * body1.mass *
                                   body2.mass / r)

    return kinetic_energy + potential_energy


def initialize_moon_position_velocity(earth_position, earth_velocity,
                                      semimajor_axis_m, eccentricity,
                                      inclination_deg, earth_mass_kg):
    """
    Initialize Moon position and velocity relative to Earth.

    Args:
        earth_position: Earth's position vector in meters
        earth_velocity: Earth's velocity vector in m/s
        semimajor_axis_m: Moon's semi-major axis around Earth in meters
        eccentricity: Moon's orbital eccentricity
        inclination_deg: Moon's orbital inclination in degrees
        earth_mass_kg: Earth's mass in kg

    Returns:
        tuple: (position_vector_m, velocity_vector_m_s)
    """
    inclination_rad = math.radians(inclination_deg)

    # Periapsis distance (closest approach to Earth)
    periapsis_distance = semimajor_axis_m * (1 - eccentricity)

    # Moon's position relative to Earth (in Earth's orbital plane + inclination)
    moon_relative_pos = vector(
        periapsis_distance * math.cos(inclination_rad),
        0,
        periapsis_distance * math.sin(inclination_rad)
    )

    # Moon's absolute position (Earth's position + relative position)
    moon_position = earth_position + moon_relative_pos

    # Moon's orbital velocity around Earth using vis-viva equation
    # v^2 = G*M*(2/r - 1/a)
    v_magnitude = math.sqrt(
        GRAVITATIONAL_CONSTANT * earth_mass_kg *
        (2.0 / periapsis_distance - 1.0 / semimajor_axis_m)
    )

    # Moon's velocity relative to Earth (perpendicular to position)
    moon_relative_vel = vector(
        0,
        v_magnitude * math.cos(inclination_rad),
        -v_magnitude * math.sin(inclination_rad)
    )

    # Moon's absolute velocity (Earth's velocity + relative velocity)
    moon_velocity = earth_velocity + moon_relative_vel

    return moon_position, moon_velocity


def create_solar_system():
    """
    Create and initialize all celestial bodies in the solar system.

    Returns:
        list: List of CelestialBody objects [Sun, Mercury, Venus, ..., Moon]
    """
    bodies = []

    # Create the Sun at origin
    sun = CelestialBody(
        name="Sun",
        mass_kg=SUN_MASS_KG,
        radius_m=SUN_RADIUS_M,
        position_m=vector(0, 0, 0),
        velocity_m_s=vector(0, 0, 0),
        body_color=(1, 1, 0),  # Base color if texture fails
        texture=SUN_TEXTURE,
        make_trail=False
    )
    bodies.append(sun)

    # Create all 8 planets
    earth_body = None
    for planet_name, data in PLANETS_DATA.items():
        position, velocity = initialize_planet_position_velocity(
            semimajor_axis_au=data['semimajor_axis_au'],
            eccentricity=data['eccentricity'],
            inclination_deg=data['inclination_deg'],
            long_asc_node_deg=data['long_asc_node_deg'],
            arg_periapsis_deg=data['arg_periapsis_deg'],
            sun_mass_kg=SUN_MASS_KG
        )

        planet = CelestialBody(
            name=planet_name,
            mass_kg=data['mass_kg'],
            radius_m=data['radius_m'],
            position_m=position,
            velocity_m_s=velocity,
            body_color=data['color'],
            texture=data.get('texture'),
            make_trail=True
        )
        bodies.append(planet)

        # Save Earth reference for Moon initialization
        if planet_name == 'Earth':
            earth_body = planet

    # Create the Moon orbiting Earth
    if earth_body is not None:
        moon_position, moon_velocity = initialize_moon_position_velocity(
            earth_position=earth_body.position,
            earth_velocity=earth_body.velocity,
            semimajor_axis_m=MOON_DATA['semimajor_axis_m'],
            eccentricity=MOON_DATA['eccentricity'],
            inclination_deg=MOON_DATA['inclination_deg'],
            earth_mass_kg=earth_body.mass
        )

        moon = CelestialBody(
            name="Moon",
            mass_kg=MOON_DATA['mass_kg'],
            radius_m=MOON_DATA['radius_m'],
            position_m=moon_position,
            velocity_m_s=moon_velocity,
            body_color=MOON_DATA['color'],
            texture=MOON_DATA.get('texture'),
            make_trail=True
        )
        bodies.append(moon)
    
    # Initialize data needed for Kepler's laws validation
    initialize_keplerian_validation(bodies)

    return bodies


def initialize_keplerian_validation(bodies):
    """
    Calculate and store initial orbital parameters for validation against Kepler's Laws.
    - 2nd Law: Conservation of angular momentum (equal areas in equal time).
    - 3rd Law: Conservation of orbital energy / semi-major axis.
    """
    # Find Sun and Earth to use as parent bodies
    sun = next((b for b in bodies if b.name == 'Sun'), None)
    earth = next((b for b in bodies if b.name == 'Earth'), None)

    for body in bodies:
        if body.name == 'Sun':
            continue  # Sun doesn't orbit anything in this model

        parent = None
        if body.name == 'Moon':
            parent = earth
        else:
            parent = sun
        
        if parent is None:
            continue

        body.parent_body = parent
        
        # Calculate initial parameters relative to the parent body
        r_vec = body.position - parent.position
        v_vec = body.velocity - parent.velocity
        r_mag = mag(r_vec)
        v_mag_sq = mag(v_vec)**2

        # Kepler's 2nd Law: Specific Angular Momentum h = |r x v|
        # This value should be conserved for a perfect two-body orbit.
        h_vector = r_vec.cross(v_vec)
        body.initial_specific_angular_momentum = mag(h_vector)

        # Kepler's 3rd Law is tied to orbital energy. We check conservation of the semi-major axis 'a'.
        # Specific Orbital Energy: E = v^2/2 - GM/r
        # Semi-major axis: a = -GM / (2E)
        mu = GRAVITATIONAL_CONSTANT * (parent.mass + body.mass)
        specific_energy = v_mag_sq / 2 - mu / r_mag
        
        if abs(specific_energy) > 1e-6:
             body.initial_semimajor_axis = -mu / (2 * specific_energy)
        else:
            body.initial_semimajor_axis = float('inf') # Parabolic/hyperbolic case


def setup_scene():
    """
    Set up the VPython canvas with appropriate settings and controls.

    Returns:
        canvas: VPython canvas object
    """
    scene = canvas(
        title='Solar System Simulation - TRUE SCALE (Real NASA Data)',
        width=1400,
        height=900,
        center=vector(0, 0, 0),
        background=color.black,
        range=0.02  # Start zoomed in on Sun (true scale is tiny!)
    )

    # Add physics validation title
    scene.append_to_caption("\n<b>PHYSICS VALIDATION (KEPLER'S LAWS)</b>\n")

    # Add instructions
    scene.append_to_caption("""
    Solar System - 100% TRUE SCALE + 3D ORBITAL MECHANICS

    ALL VALUES 100% REAL NASA DATA:
    ✓ Gravitational constant: 6.67430×10⁻¹¹ m³/(kg·s²)
    ✓ All masses: Real NASA JPL values (kg)
    ✓ All radii: Real NASA JPL values (m) - SAME SCALE AS DISTANCES
    ✓ All orbital distances: Real semi-major axes (AU)
    ✓ All velocities: Calculated from vis-viva equation
    ✓ All eccentricities: Real values (Mercury=0.206, etc.)
    ✓ ALL INCLINATIONS: Real 3D orbital planes implemented!
    ✓ Longitude of ascending node (Ω): Real values
    ✓ Argument of periapsis (ω): Real values
    ✓ N-body gravity: F = G×m₁×m₂/r² for all pairs
    ✓ Integration: Velocity Verlet
    ✓ Moon: 384,400 km from Earth, proper 3D orbit

    3D ORBITAL MECHANICS:
    - Mercury tilted 7° (most inclined planet)
    - Venus tilted 3.4°
    - Mars tilted 1.8°
    - All planets in their real orbital planes!
    - Proper rotation matrices for 3D space

    CAMERA CONTROLS:
    - Press 0-9: Lock to body (3=Earth, 9=Moon)
    - Press C: Free camera
    - Scroll OUT to see 3D orbital structure
    - Right-drag: Rotate | Ctrl+drag: Pan

    Time: 1 sim sec = 10 minutes | Thin lines = orbital paths
    """)

    return scene


def setup_plots():
    """
    Create a graph window for plotting real-time simulation data.
    """
    g = graph(
        title="<b>Real-time Physics Validation</b>",
        xtitle="Time (years)",
        ytitle="Error (%)",
        width=800,
        height=400,
        align='right'
    )
    
    total_energy_curve = gcurve(
        graph=g,
        color=color.red,
        label="Total System Energy Drift (%)"
    )
    
    earth_h_error_curve = gcurve(
        graph=g,
        color=color.blue,
        label="Earth Ang. Momentum Error (2nd Law) (%)"
    )

    earth_a_error_curve = gcurve(
        graph=g,
        color=color.green,
        label="Earth Semi-Major Axis Error (3rd Law) (%)"
    )
    
    return {
        "total_energy": total_energy_curve,
        "earth_h_error": earth_h_error_curve,
        "earth_a_error": earth_a_error_curve,
    }


def simulation_loop(bodies, dt_seconds, scene, plots):
    """
    Main simulation loop using Velocity Verlet integration.

    Args:
        bodies: List of CelestialBody objects
        dt_seconds: Time step in real seconds
        scene: VPython canvas object for camera control
        plots: Dictionary of gcurve objects for plotting
    """
    # Calculate initial energy for monitoring
    initial_energy = calculate_total_energy(bodies)

    # Setup for Kepler's Laws validation display
    # A single wtext object is used and updated dynamically.
    validation_wtext = wtext(text="", pos=scene.caption_anchor)
    scene.append_to_caption("\n") # Add spacing

    simulation_time = 0  # seconds
    frame_count = 0

    # Camera lock variables
    camera_lock_index = None  # None = free camera, otherwise index in bodies list
    camera_lock_body = None

    print("Starting simulation...")
    print(f"Time step: {dt_seconds / 60:.1f} minutes ({dt_seconds / 3600:.3f} hours, {dt_seconds / 86400:.5f} days)")
    print(f"Initial total energy: {initial_energy:.3e} J")
    print("\n=== CAMERA CONTROLS ===")
    print("Press 0-9 to lock camera on bodies:")
    for i, body in enumerate(bodies):
        if i <= 9:
            print(f"  {i}: {body.name}")
    print("  C: Free camera (unlock)")
    print("\nPress Ctrl+C to stop.\n")

    try:
        while True:
            rate(60)  # 60 frames per second max

            # Check for keyboard input for camera lock
            keys = keysdown()
            if 'c' in keys or 'C' in keys:
                camera_lock_index = None
                camera_lock_body = None
                print("Camera unlocked - free camera mode")
            else:
                # Check for number keys 0-9
                for i in range(min(10, len(bodies))):
                    if str(i) in keys:
                        camera_lock_index = i
                        camera_lock_body = bodies[i]
                        print(f"Camera locked to: {bodies[i].name}")
                        break

            # Update camera position if locked to a body
            if camera_lock_body is not None:
                scene.center = camera_lock_body.sphere.pos

            # Velocity Verlet integration
            # Step 1: Calculate current accelerations
            for body in bodies:
                body.acceleration = calculate_total_acceleration(body, bodies)

            # Step 2: Update positions using current velocity and acceleration
            for body in bodies:
                body.update_position(dt_seconds)

            # Step 3: Calculate new accelerations at new positions
            new_accelerations = []
            for body in bodies:
                new_acc = calculate_total_acceleration(body, bodies)
                new_accelerations.append(new_acc)

            # Step 4: Update velocities using average of old and new accelerations
            for body, new_acc in zip(bodies, new_accelerations):
                body.update_velocity(new_acc, dt_seconds)

            # Update simulation time
            simulation_time += dt_seconds
            frame_count += 1

            # Update status and physics validation every 20 frames for performance
            if frame_count % 20 == 0:
                # --- Total System Energy (Integrator Stability) ---
                current_energy = calculate_total_energy(bodies)
                energy_change = abs((current_energy - initial_energy) / initial_energy) * 100
                days_elapsed = simulation_time / 86400
                years_elapsed = days_elapsed / 365.25
                
                # Update console output less frequently
                if frame_count % 100 == 0:
                    print(f"Time: {years_elapsed:.4f} years ({days_elapsed:.2f} days) | "
                          f"Total Energy Drift: {energy_change:.6f}%")
                
                plots["total_energy"].plot(years_elapsed, energy_change)
                
                # --- Per-Body Keplerian Validation ---
                max_error_body = None
                max_h_error = -1

                # First, calculate error for all bodies
                all_body_errors = {}
                for body in bodies:
                    if body.parent_body is not None:
                        r_vec = body.position - body.parent_body.position
                        v_vec = body.velocity - body.parent_body.velocity
                        r_mag = mag(r_vec)
                        v_mag_sq = mag(v_vec)**2

                        h_current = mag(r_vec.cross(v_vec))
                        h_error = 100 * abs(h_current - body.initial_specific_angular_momentum) / body.initial_specific_angular_momentum

                        mu = GRAVITATIONAL_CONSTANT * (body.parent_body.mass + body.mass)
                        specific_energy = v_mag_sq / 2 - mu / r_mag
                        a_current = -mu / (2 * specific_energy) if abs(specific_energy) > 1e-6 else float('inf')
                        a_error = 100 * abs(a_current - body.initial_semimajor_axis) / body.initial_semimajor_axis if body.initial_semimajor_axis > 0 else 0
                        
                        all_body_errors[body.name] = (h_error, a_error)

                        if h_error > max_h_error:
                            max_h_error = h_error
                            max_error_body = body
                        
                        # Update plots for specific bodies
                        if body.name == "Earth":
                            plots["earth_h_error"].plot(years_elapsed, h_error)
                            plots["earth_a_error"].plot(years_elapsed, a_error)

                
                # --- Update On-Screen Text ---
                body_to_display = None
                if camera_lock_body is not None and camera_lock_body.parent_body is not None:
                    # If camera is locked, show info for that body
                    body_to_display = camera_lock_body
                elif max_error_body is not None:
                    # Otherwise, show the body with the max error
                    body_to_display = max_error_body

                if body_to_display and body_to_display.name in all_body_errors:
                    h_err, a_err = all_body_errors[body_to_display.name]
                    validation_wtext.text = (f"<b>Displaying: {body_to_display.name}</b> | "
                                             f"Ang. Momentum Error (2nd Law): {h_err:.5f}% | "
                                             f"Semi-Major Axis Error (3rd Law): {a_err:.5f}%")
                else:
                    validation_wtext.text = "Free camera. No body selected."


    except KeyboardInterrupt:
        print("\nSimulation stopped by user.")
        print(f"Total simulation time: {simulation_time / (86400 * 365.25):.2f} years")


def main():
    """
    Main entry point for the solar system simulation.
    """
    print("=" * 70)
    print("SOLAR SYSTEM SIMULATION")
    print("Real NASA/JPL Data - Physics-Accurate N-Body Simulation")
    print("=" * 70)
    print()

    # Set up visualization
    scene = setup_scene()
    plots = setup_plots()

    # Create all celestial bodies
    bodies = create_solar_system()

    print(f"Created {len(bodies)} celestial bodies:")
    for body in bodies:
        print(f"  - {body.name}")
    print()

    # Time step selection
    # 10 minute timestep provides smooth animation and excellent accuracy
    # Much smaller than orbital periods (Mercury: 88 days = 12,672 ten-minute intervals)
    dt_seconds = TIME_SCALE  # 10 minutes per simulation step

    print(f"Physics parameters:")
    print(f"  Gravitational constant: {GRAVITATIONAL_CONSTANT:.5e} m^3 kg^-1 s^-2")
    print(f"  Time step: {dt_seconds / 60:.1f} minutes ({dt_seconds / 3600:.3f} hours)")
    print(f"  Integration method: Velocity Verlet")
    print()

    # Start simulation
    simulation_loop(bodies, dt_seconds, scene, plots)


if __name__ == "__main__":
    main()