from orbit_sim import OrbitalSim, Body, circular_orbit_velocity, G
sim = OrbitalSim()
sun = sim.add(Body("Sun", 1e30, 0, 0))
r = 1e10
v = circular_orbit_velocity(sun.mass, r)
earth = sim.add(Body("Earth", 1e24, r, 0, 0, v))
e0 = sim.total_energy()
sim.run(100, 1000)
e1 = sim.total_energy()
assert abs(e1 - e0) / abs(e0) < 0.01, f"Energy conservation violated: {e0} vs {e1}"
print("Orbit sim tests passed")