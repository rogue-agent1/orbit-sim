#!/usr/bin/env python3
"""Orbital mechanics simulator."""
import math

G = 6.674e-11  # gravitational constant

class Body:
    def __init__(self, mass, x, y, vx=0, vy=0):
        self.mass = mass
        self.x, self.y = x, y
        self.vx, self.vy = vx, vy

def gravity_force(a, b):
    dx, dy = b.x - a.x, b.y - a.y
    r = math.sqrt(dx**2 + dy**2)
    if r < 1: return 0, 0
    f = G * a.mass * b.mass / r**2
    return f * dx/r, f * dy/r

def step(bodies, dt):
    forces = [(0, 0)] * len(bodies)
    for i in range(len(bodies)):
        fx, fy = 0, 0
        for j in range(len(bodies)):
            if i != j:
                dfx, dfy = gravity_force(bodies[i], bodies[j])
                fx += dfx; fy += dfy
        forces[i] = (fx, fy)
    for i, b in enumerate(bodies):
        ax, ay = forces[i][0]/b.mass, forces[i][1]/b.mass
        b.vx += ax * dt; b.vy += ay * dt
        b.x += b.vx * dt; b.y += b.vy * dt

def orbital_period(mass_central, radius):
    return 2 * math.pi * math.sqrt(radius**3 / (G * mass_central))

def orbital_velocity(mass_central, radius):
    return math.sqrt(G * mass_central / radius)

def energy(bodies):
    ke = sum(0.5 * b.mass * (b.vx**2 + b.vy**2) for b in bodies)
    pe = 0
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            dx = bodies[j].x - bodies[i].x
            dy = bodies[j].y - bodies[i].y
            r = math.sqrt(dx**2 + dy**2)
            if r > 0: pe -= G * bodies[i].mass * bodies[j].mass / r
    return ke + pe

if __name__ == "__main__":
    sun = Body(1.989e30, 0, 0)
    earth = Body(5.972e24, 1.496e11, 0, 0, 29780)
    for _ in range(365):
        for _ in range(24*3600//3600):
            step([sun, earth], 3600)
    print(f"Earth position after ~1yr: ({earth.x:.2e}, {earth.y:.2e})")

def test():
    # Orbital velocity
    v = orbital_velocity(1.989e30, 1.496e11)
    assert abs(v - 29780) < 100  # ~29.78 km/s
    # Period (Earth ~365.25 days)
    T = orbital_period(1.989e30, 1.496e11)
    assert abs(T / 86400 - 365.25) < 1
    # Energy conservation (short sim)
    sun = Body(1e30, 0, 0)
    planet = Body(1e24, 1e11, 0, 0, orbital_velocity(1e30, 1e11))
    e0 = energy([sun, planet])
    for _ in range(100):
        step([sun, planet], 1000)
    e1 = energy([sun, planet])
    assert abs((e1 - e0) / e0) < 0.01  # <1% energy drift
    print("  orbit_sim: ALL TESTS PASSED")
