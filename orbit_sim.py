#!/usr/bin/env python3
"""orbit_sim - Kepler orbit simulation."""
import sys,argparse,json,math
def simulate_orbit(mass_central=1e6,mass_orbiter=1,r0=100,v0=None,G=1,dt=0.01,steps=10000):
    if v0 is None:v0=math.sqrt(G*mass_central/r0)
    x,y=r0,0.0;vx,vy=0.0,v0;history=[]
    for i in range(steps):
        r=math.sqrt(x**2+y**2)
        if r<1:break
        a=-G*mass_central/r**3
        ax=a*x;ay=a*y
        vx+=ax*dt;vy+=ay*dt;x+=vx*dt;y+=vy*dt
        if i%(steps//100)==0:
            ke=0.5*mass_orbiter*(vx**2+vy**2)
            pe=-G*mass_central*mass_orbiter/r
            history.append({"t":round(i*dt,2),"x":round(x,2),"y":round(y,2),"r":round(r,2),"energy":round(ke+pe,2)})
    return history
def main():
    p=argparse.ArgumentParser(description="Orbit simulator")
    p.add_argument("--radius",type=float,default=100);p.add_argument("--velocity",type=float)
    p.add_argument("--steps",type=int,default=50000);p.add_argument("--mass",type=float,default=1e6)
    args=p.parse_args()
    history=simulate_orbit(args.mass,r0=args.radius,v0=args.velocity,steps=args.steps)
    min_r=min(h["r"] for h in history) if history else 0
    max_r=max(h["r"] for h in history) if history else 0
    print(json.dumps({"initial_radius":args.radius,"periapsis":round(min_r,2),"apoapsis":round(max_r,2),"eccentricity":round((max_r-min_r)/(max_r+min_r),4) if max_r+min_r>0 else 0,"trajectory_points":len(history),"sample":history[:20]},indent=2))
if __name__=="__main__":main()
