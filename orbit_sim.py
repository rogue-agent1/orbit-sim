#!/usr/bin/env python3
"""Keplerian orbits and orbital mechanics."""
import sys, math

def orbit(r0, v0, G=1, M=1000, dt=0.001, steps=10000):
    x, y = r0; vx, vy = v0; traj = [(x,y)]
    for _ in range(steps):
        r = math.sqrt(x**2+y**2)+1e-10
        ax = -G*M*x/r**3; ay = -G*M*y/r**3
        vx += ax*dt; vy += ay*dt; x += vx*dt; y += vy*dt
        traj.append((x,y))
    return traj

def orbital_elements(r, v, mu=1000):
    x,y = r; vx,vy = v; r_mag = math.sqrt(x**2+y**2)
    v_mag = math.sqrt(vx**2+vy**2)
    energy = 0.5*v_mag**2 - mu/r_mag
    a = -mu/(2*energy) if energy != 0 else float('inf')
    h = x*vy - y*vx
    e = math.sqrt(1 + 2*energy*h**2/mu**2) if mu > 0 else 0
    return {"semi_major": a, "eccentricity": e, "energy": energy, "angular_momentum": h}

def main():
    # Circular orbit
    r = 10; v_circ = math.sqrt(1000/r)
    traj = orbit((r,0),(0,v_circ), steps=20000)
    elem = orbital_elements((r,0),(0,v_circ))
    print(f"Circular: a={elem['semi_major']:.2f}, e={elem['eccentricity']:.4f}")
    # Elliptical
    traj2 = orbit((r,0),(0,v_circ*0.7), steps=20000)
    elem2 = orbital_elements((r,0),(0,v_circ*0.7))
    print(f"Elliptical: a={elem2['semi_major']:.2f}, e={elem2['eccentricity']:.4f}")

if __name__ == "__main__": main()
