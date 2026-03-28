#!/usr/bin/env python3
"""Orbital mechanics — two-body simulation, Kepler elements."""
import sys, math
def simulate(m1=1e10, m2=1, r0=100, v0=None, G=6.674e-11, dt=0.1, steps=1000):
    if v0 is None: v0=math.sqrt(G*m1/r0)
    x,y=r0,0.0; vx,vy=0.0,v0; points=[]
    for _ in range(steps):
        r=math.sqrt(x*x+y*y); a=-G*m1/(r*r*r)
        vx+=a*x*dt; vy+=a*y*dt; x+=vx*dt; y+=vy*dt
        points.append((x,y))
    return points
def ascii_orbit(points, w=60, h=30):
    xs=[p[0] for p in points]; ys=[p[1] for p in points]
    mnx,mxx=min(xs),max(xs); mny,mxy=min(ys),max(ys)
    rx,ry=mxx-mnx or 1, mxy-mny or 1
    grid=[[" "]*w for _ in range(h)]
    grid[h//2][w//2]="*"
    for x,y in points:
        c=int((x-mnx)/rx*(w-1)); r=int((y-mny)/ry*(h-1))
        if 0<=r<h and 0<=c<w: grid[r][c]="."
    for row in grid: print("".join(row))
def cli():
    steps=int(sys.argv[1]) if len(sys.argv)>1 else 2000
    points=simulate(steps=steps)
    ascii_orbit(points)
    xs=[p[0] for p in points]; ys=[p[1] for p in points]
    print(f"  Steps: {steps}, X range: [{min(xs):.1f}, {max(xs):.1f}]")
if __name__=="__main__": cli()
