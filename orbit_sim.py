#!/usr/bin/env python3
"""Orbital Mechanics - N-body gravitational simulation with Verlet integration."""
import sys, math

G = 6.674e-11

class Body:
    def __init__(self, name, mass, x, y, vx=0, vy=0):
        self.name=name;self.mass=mass;self.x=x;self.y=y;self.vx=vx;self.vy=vy;self.ax=0;self.ay=0

def compute_forces(bodies):
    for b in bodies: b.ax=0;b.ay=0
    for i,a in enumerate(bodies):
        for b in bodies[i+1:]:
            dx=b.x-a.x;dy=b.y-a.y;r2=dx*dx+dy*dy;r=math.sqrt(r2)
            if r<1e6: continue
            f=G*a.mass*b.mass/r2;fx=f*dx/r;fy=f*dy/r
            a.ax+=fx/a.mass;a.ay+=fy/a.mass
            b.ax-=fx/b.mass;b.ay-=fy/b.mass

def step(bodies, dt):
    for b in bodies:
        b.x+=b.vx*dt+0.5*b.ax*dt*dt
        b.y+=b.vy*dt+0.5*b.ay*dt*dt
    old_a = [(b.ax,b.ay) for b in bodies]
    compute_forces(bodies)
    for i,b in enumerate(bodies):
        b.vx+=0.5*(old_a[i][0]+b.ax)*dt
        b.vy+=0.5*(old_a[i][1]+b.ay)*dt

def energy(bodies):
    ke=sum(0.5*b.mass*(b.vx**2+b.vy**2) for b in bodies)
    pe=0
    for i,a in enumerate(bodies):
        for b in bodies[i+1:]:
            r=math.sqrt((a.x-b.x)**2+(a.y-b.y)**2)
            if r>0: pe-=G*a.mass*b.mass/r
    return ke+pe

def main():
    AU=1.496e11
    sun=Body("Sun",1.989e30,0,0)
    earth=Body("Earth",5.972e24,AU,0,0,29780)
    mars=Body("Mars",6.39e23,1.524*AU,0,0,24077)
    bodies=[sun,earth,mars]
    compute_forces(bodies)
    dt=86400;days=365
    print(f"=== Orbital Mechanics ({days} days) ===\n")
    print(f"Initial energy: {energy(bodies):.3e} J")
    for day in range(days):
        step(bodies, dt)
        if day % 90 == 0:
            print(f"  Day {day:3d}:", end="")
            for b in bodies[1:]:
                r=math.sqrt(b.x**2+b.y**2)/AU
                print(f" {b.name}={r:.3f}AU", end="")
            print()
    print(f"\nFinal energy: {energy(bodies):.3e} J")
    drift=abs(energy(bodies)-energy([Body("Sun",1.989e30,0,0),Body("Earth",5.972e24,AU,0,0,29780),Body("Mars",6.39e23,1.524*AU,0,0,24077)]))

if __name__ == "__main__":
    main()
