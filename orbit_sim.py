#!/usr/bin/env python3
"""Orbital mechanics simulator. Zero dependencies."""
import math, sys

G = 6.674e-11

class Body:
    def __init__(self, name, mass, x, y, vx=0, vy=0):
        self.name = name
        self.mass = mass
        self.x, self.y = float(x), float(y)
        self.vx, self.vy = float(vx), float(vy)
        self.trail = []

    def dist_to(self, other):
        return math.sqrt((self.x-other.x)**2 + (self.y-other.y)**2)

    def kinetic_energy(self):
        return 0.5 * self.mass * (self.vx**2 + self.vy**2)

class OrbitalSim:
    def __init__(self):
        self.bodies = []

    def add(self, body):
        self.bodies.append(body); return body

    def step(self, dt):
        forces = [(0, 0)] * len(self.bodies)
        for i in range(len(self.bodies)):
            fx, fy = 0, 0
            for j in range(len(self.bodies)):
                if i == j: continue
                dx = self.bodies[j].x - self.bodies[i].x
                dy = self.bodies[j].y - self.bodies[i].y
                r = math.sqrt(dx**2 + dy**2)
                if r < 1e-10: continue
                f = G * self.bodies[i].mass * self.bodies[j].mass / (r*r)
                fx += f * dx/r; fy += f * dy/r
            forces[i] = (fx, fy)
        for i, b in enumerate(self.bodies):
            ax = forces[i][0] / b.mass
            ay = forces[i][1] / b.mass
            b.vx += ax * dt; b.vy += ay * dt
            b.x += b.vx * dt; b.y += b.vy * dt
            b.trail.append((b.x, b.y))

    def run(self, steps, dt):
        for _ in range(steps):
            self.step(dt)

    def total_energy(self):
        ke = sum(b.kinetic_energy() for b in self.bodies)
        pe = 0
        for i in range(len(self.bodies)):
            for j in range(i+1, len(self.bodies)):
                r = self.bodies[i].dist_to(self.bodies[j])
                if r > 0:
                    pe -= G * self.bodies[i].mass * self.bodies[j].mass / r
        return ke + pe

def circular_orbit_velocity(M, r):
    return math.sqrt(G * M / r)

if __name__ == "__main__":
    sim = OrbitalSim()
    sun = sim.add(Body("Sun", 1.989e30, 0, 0))
    v = circular_orbit_velocity(sun.mass, 1.496e11)
    earth = sim.add(Body("Earth", 5.972e24, 1.496e11, 0, 0, v))
    sim.run(365, 86400)
    print(f"Earth after 1 year: ({earth.x:.3e}, {earth.y:.3e})")
    print(f"Distance from Sun: {earth.dist_to(sun):.3e} m")
