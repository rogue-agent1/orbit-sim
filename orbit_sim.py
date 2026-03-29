import argparse, math

def simulate(altitude_km=400, steps=5000, dt=10):
    G = 6.674e-11
    M = 5.972e24  # Earth mass
    R = 6.371e6   # Earth radius
    r0 = R + altitude_km * 1000
    v_orbital = math.sqrt(G * M / r0)
    x, y = r0, 0
    vx, vy = 0, v_orbital
    history = []
    for step in range(steps):
        r = math.sqrt(x*x + y*y)
        a = -G * M / (r * r)
        ax = a * x / r
        ay = a * y / r
        vx += ax * dt
        vy += ay * dt
        x += vx * dt
        y += vy * dt
        alt = (math.sqrt(x*x + y*y) - R) / 1000
        if step % (steps//20) == 0:
            history.append((step*dt, x/1e6, y/1e6, alt, math.sqrt(vx*vx+vy*vy)))
    return history

def plot_orbit(history, w=50, h=25):
    xs = [p[1] for p in history]
    ys = [p[2] for p in history]
    xmin, xmax = min(xs)-0.5, max(xs)+0.5
    ymin, ymax = min(ys)-0.5, max(ys)+0.5
    grid = [[" "]*w for _ in range(h)]
    # Draw Earth
    cx = int((0 - xmin) / (xmax - xmin) * (w-1))
    cy = int((0 - ymin) / (ymax - ymin) * (h-1))
    if 0 <= cy < h and 0 <= cx < w: grid[cy][cx] = "⊕"
    for x, y in zip(xs, ys):
        col = int((x - xmin) / (xmax - xmin) * (w-1))
        row = int((y - ymin) / (ymax - ymin) * (h-1))
        if 0 <= row < h and 0 <= col < w: grid[row][col] = "·"
    for r in grid: print("".join(r))

def main():
    p = argparse.ArgumentParser(description="Orbit simulator")
    p.add_argument("-a", "--altitude", type=float, default=400, help="Altitude in km")
    p.add_argument("-n", "--steps", type=int, default=5000)
    args = p.parse_args()
    history = simulate(args.altitude, args.steps)
    print(f"Orbit at {args.altitude} km altitude")
    plot_orbit(history)
    print(f"\n{'Time(s)':>8} {'X(Mm)':>8} {'Y(Mm)':>8} {'Alt(km)':>10} {'V(m/s)':>10}")
    for t, x, y, alt, v in history:
        print(f"{t:8.0f} {x:8.2f} {y:8.2f} {alt:10.1f} {v:10.1f}")

if __name__ == "__main__":
    main()
