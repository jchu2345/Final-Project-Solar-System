from vpython import *
import math
import json

# Constants
G = 6.67430e-11  # Gravitational constant
SCALE_FACTOR_PLANET_RADIUS = 1000
SCALE_FACTOR_SUN_RADIUS = 10
SCALE_FACTOR_MOON_RADIUS = 5000 # Increased for visibility

# Load data from JSON
with open('solar/sun.json') as f:
    data = json.load(f)
    celestial_data = data['celestial_data']
    moon_data = data['moon_data']


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
    sun_mass = celestial_data["Sun"]['mass']
    for name, data in celestial_data.items():
        mass = data['mass']
        radius = data['radius']
        a = data['semi_major_axis']
        e = data['eccentricity']
        incl = data['orbital_inclination']
        tilt = data['axial_tilt']
        col = vector(data['color'][0], data['color'][1], data['color'][2])
        
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
        mass = data['mass']
        radius = data['radius']
        a = data['semi_major_axis']
        e = data['eccentricity']
        incl = data['orbital_inclination']
        col = vector(data['color'][0], data['color'][1], data['color'][2])

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
    
    # Add rings to Saturn
    saturn = bodies["Saturn"]
    saturn_rings = ring(pos=saturn.pos,
                        axis=saturn.axis,
                        radius=saturn.radius * 2.5,
                        thickness=10,
                        color=saturn.color)

    # Camera control
    def follow_planet(menu):
        planet_name = menu.selected
        if planet_name == "Free":
            scene.camera.follow(None)
        elif planet_name in bodies:
            scene.camera.follow(bodies[planet_name])

    planet_names = ["Free"] + list(bodies.keys())
    menu(choices=planet_names, bind=follow_planet, pos=scene.title_anchor)

    # Simulation loop
    dt = 600  # 10 minutes in seconds
    while True:
        rate(1000)  # Increased rate for smoother animation

        # Update forces
        all_bodies = list(bodies.values()) + [sun]
        for i in range(len(all_bodies)):
            body1 = all_bodies[i]
            if body1.name == 'Sun': continue
            body1.force = vector(0, 0, 0)

            for j in range(len(all_bodies)):
                if i == j: continue
                body2 = all_bodies[j]
                r = body1.pos - body2.pos
                if mag(r) > 0: # Avoid division by zero if bodies collide
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