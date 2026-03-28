#!/usr/bin/env python3
"""Orbital mechanics simulator (Kepler orbits)."""
import math
def orbit(a=1.0,e=0.3,mu=1.0,steps=360):
    points=[]
    for i in range(steps+1):
        theta=2*math.pi*i/steps;r=a*(1-e*e)/(1+e*math.cos(theta))
        x=r*math.cos(theta);y=r*math.sin(theta)
        v=math.sqrt(mu*(2/r-1/a))
        points.append((x,y,r,v,theta))
    return points
def hohmann(r1,r2,mu=1.0):
    a_t=(r1+r2)/2;dv1=math.sqrt(mu/r1)*(math.sqrt(2*r2/(r1+r2))-1)
    dv2=math.sqrt(mu/r2)*(1-math.sqrt(2*r1/(r1+r2)));t=math.pi*math.sqrt(a_t**3/mu)
    return dv1,dv2,t
if __name__=="__main__":
    pts=orbit(a=1.0,e=0.5)
    rmin=min(p[2] for p in pts);rmax=max(p[2] for p in pts)
    print(f"Orbit: periapsis={rmin:.3f}, apoapsis={rmax:.3f}")
    dv1,dv2,t=hohmann(1.0,1.524)
    print(f"Hohmann Earth->Mars: dv1={dv1:.4f}, dv2={dv2:.4f}, time={t:.2f}")
    print("Orbit sim OK")
