#!/usr/bin/env python3
"""orbit_sim - Orbital mechanics simulator."""
import argparse, math, sys, time

G = 6.674e-11

class Body:
    def __init__(self, name, mass, x, y, vx=0, vy=0):
        self.name=name; self.mass=mass; self.x=x; self.y=y; self.vx=vx; self.vy=vy
        self.trail = []

def simulate(bodies, dt=3600, steps=1000):
    for _ in range(steps):
        for b in bodies: b.trail.append((b.x, b.y))
        # Calculate forces
        for i, a in enumerate(bodies):
            ax, ay = 0, 0
            for j, b in enumerate(bodies):
                if i == j: continue
                dx, dy = b.x-a.x, b.y-a.y
                r = math.sqrt(dx*dx+dy*dy)
                if r < 1e3: continue
                f = G * b.mass / (r*r)
                ax += f * dx/r; ay += f * dy/r
            a.vx += ax*dt; a.vy += ay*dt
        for b in bodies:
            b.x += b.vx*dt; b.y += b.vy*dt

def render_ascii(bodies, w=70, h=35, scale=None):
    all_x = [b.x for b in bodies]; all_y = [b.y for b in bodies]
    if not scale:
        rx = max(abs(x) for x in all_x) or 1; ry = max(abs(y) for y in all_y) or 1
        scale = max(rx, ry) * 1.2
    grid = [["."]*w for _ in range(h)]
    for b in bodies:
        # Trail
        for tx, ty in b.trail[-100:]:
            px = int((tx/scale+1)*w/2); py = int((ty/scale+1)*h/2)
            if 0<=px<w and 0<=py<h: grid[py][px] = "·"
        # Body
        px = int((b.x/scale+1)*w/2); py = int((b.y/scale+1)*h/2)
        if 0<=px<w and 0<=py<h: grid[py][px] = b.name[0]
    return "\n".join("".join(row) for row in grid)

def main():
    p = argparse.ArgumentParser(description="Orbital mechanics")
    p.add_argument("--demo", choices=["earth-moon","solar","binary"], default="earth-moon")
    p.add_argument("-s","--steps", type=int, default=200)
    a = p.parse_args()
    if a.demo == "earth-moon":
        bodies = [
            Body("Earth", 5.972e24, 0, 0),
            Body("Moon", 7.342e22, 3.844e8, 0, 0, 1022),
        ]
        dt = 3600
    elif a.demo == "binary":
        bodies = [
            Body("Star_A", 2e30, -5e10, 0, 0, 2e4),
            Body("Star_B", 2e30, 5e10, 0, 0, -2e4),
        ]
        dt = 36000
    else:
        bodies = [
            Body("Sun", 1.989e30, 0, 0),
            Body("Mercury", 3.301e23, 5.79e10, 0, 0, 47400),
            Body("Venus", 4.867e24, 1.082e11, 0, 0, 35020),
            Body("Earth", 5.972e24, 1.496e11, 0, 0, 29780),
        ]
        dt = 86400
    print(f"Simulating {a.demo} ({len(bodies)} bodies, {a.steps} steps)")
    for step in range(a.steps):
        simulate(bodies, dt, 10)
        if step % 10 == 0:
            sys.stdout.write(f"\033[2J\033[HStep {step}\n")
            print(render_ascii(bodies))
            time.sleep(0.05)
    print(f"\nFinal positions:")
    for b in bodies: print(f"  {b.name}: ({b.x:.2e}, {b.y:.2e})")

if __name__ == "__main__": main()
