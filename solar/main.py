from vpython import *
import math

# Constants
G = 6.67430e-11  # Gravitational constant
SCALE_FACTOR_PLANET_RADIUS = 1000
SCALE_FACTOR_SUN_RADIUS = 10
SCALE_FACTOR_MOON_RADIUS = 1000

# Data for celestial bodies
# Data format: [mass (kg), radius (m), semi-major axis (m), eccentricity, orbital inclination (deg), axial tilt (deg), color]
celestial_data = {
    "Sun": [1.989e30, 6.9634e8, 0, 0, 0, 0, color.yellow],
    "Mercury": [3.301e23, 2.440e6, 5.79e10, 0.205, 7.0, 0.01, color.gray(0.5)],
    "Venus": [4.867e24, 6.052e6, 1.082e11, 0.007, 3.4, 177.4, color.orange],
    "Earth": [5.972e24, 6.371e6, 1.496e11, 0.017, 0.0, 23.5, color.blue],
    "Mars": [6.417e23, 3.389e6, 2.279e11, 0.094, 1.9, 25.2, color.red],
    "Jupiter": [1.898e27, 6.991e7, 7.786e11, 0.049, 1.3, 3.1, color.brown],
    "Saturn": [5.683e26, 5.823e7, 1.4335e12, 0.057, 2.5, 26.7, color.yellow],
    "Uranus": [8.681e25, 2.536e7, 2.8725e12, 0.046, 0.8, 97.8, color.cyan],
    "Neptune": [1.024e26, 2.462e7, 4.4951e12, 0.011, 1.8, 28.3, color.blue],
    "Pluto": [1.309e22, 1.188e6, 5.9064e12, 0.244, 17.2, 122.5, color.gray(0.7)]
}

moon_data = {
    "Moon": [7.342e22, 1.737e6, 3.844e8, 0.0549, 5.1, color.gray(0.8)]
}

class CelestialBody(sphere):
    def __init__(self, name, mass, radius, color, pos, velocity, axis):
        super().__init__(pos=pos, radius=radius, color=color, make_trail=True)
        self.mass = mass
        self.velocity = velocity
        self.name = name
        self.axis = axis

def main():
    # Scene setup
    scene.caption = "Solar System Simulation"
    scene.background = color.black

    # Create celestial bodies
    bodies = {}
    sun_mass = celestial_data["Sun"][0]
    for name, data in celestial_data.items():
        mass, radius, a, e, incl, tilt, col = data
        
        axis = vector(0,1,0).rotate(angle=radians(tilt), axis=vector(0,0,1))

        if name == "Sun":
            radius *= SCALE_FACTOR_SUN_RADIUS
            pos = vector(0,0,0)
            vel = vector(0,0,0)
        else:
            radius *= SCALE_FACTOR_PLANET_RADIUS
            # Position at perihelion
            x = a * (1 - e)
            # Velocity at perihelion
            v = math.sqrt(G * sun_mass * (1 + e) / (a * (1 - e)))
            
            pos = vector(x, 0, 0)
            vel = vector(0, v, 0)
            
            # Apply orbital inclination
            incl_rad = radians(incl)
            pos = pos.rotate(angle=incl_rad, axis=vector(1,0,0))
            vel = vel.rotate(angle=incl_rad, axis=vector(1,0,0))

        bodies[name] = CelestialBody(name=name,
                                     mass=mass,
                                     radius=radius,
                                     color=col,
                                     pos=pos,
                                     velocity=vel,
                                     axis=axis)

    # Create moons
    for name, data in moon_data.items():
        mass, radius, a, e, incl, col = data
        radius *= SCALE_FACTOR_MOON_RADIUS
        
        # Position and velocity relative to Earth
        earth = bodies["Earth"]
        x = a * (1 - e)
        v = math.sqrt(G * earth.mass * (1 + e) / (a * (1 - e)))
        
        pos_rel = vector(x, 0, 0)
        vel_rel = vector(0, v, 0)
        
        # Apply orbital inclination
        incl_rad = radians(incl)
        pos_rel = pos_rel.rotate(angle=incl_rad, axis=vector(1,0,0))
        vel_rel = vel_rel.rotate(angle=incl_rad, axis=vector(1,0,0))
        
        bodies[name] = CelestialBody(name=name,
                                     mass=mass,
                                     radius=radius,
                                     color=col,
                                     pos=earth.pos + pos_rel,
                                     velocity=earth.velocity + vel_rel,
                                     axis=vector(0,1,0))


    sun = bodies.pop("Sun")
    planets = [p for p in bodies.values() if p.name not in moon_data]
    moons = [m for m in bodies.values() if m.name in moon_data]
    
    # Add rings to Saturn
    saturn = bodies["Saturn"]
    saturn_rings = ring(pos=saturn.pos,
                        axis=saturn.axis,
                        radius=saturn.radius * 2.5,
                        thickness=10,
                        color=saturn.color)


    # Simulation loop
    dt = 3600  # 1 hour in seconds
    while True:
        rate(200)  # Limit the loop rate

        # Update forces
        for body1 in bodies.values():
            body1.force = vector(0, 0, 0)
            # Force of Sun
            r = body1.pos - sun.pos
            body1.force += -G * sun.mass * body1.mass * r.hat / mag2(r)

            # Force of other bodies
            for body2 in bodies.values():
                if body1 is not body2:
                    r = body1.pos - body2.pos
                    body1.force += -G * body2.mass * body1.mass * r.hat / mag2(r)

        # Update positions
        for body in bodies.values():
            # Euler-Cromer integration
            body.velocity = body.velocity + body.force / body.mass * dt
            body.pos = body.pos + body.velocity * dt
            
        # Update Saturn's rings position
        saturn_rings.pos = saturn.pos


if __name__ == "__main__":
    main()