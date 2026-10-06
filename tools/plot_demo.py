#!/usr/bin/env python3
"""绘制 demo_replan 的输出（经纬度平面示意 + 曲率剖面）。用法: python3 plot_demo.py <out_dir>"""
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 6371000.0  # 与 demo 中一致（地球平均半径）


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def cap_boundary(lat, lon, r_m, n=180):
    la, lo, rr = math.radians(lat), math.radians(lon), r_m / R
    c = (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))
    north = (-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la))
    east = (-math.sin(lo), math.cos(lo), 0.0)
    pts = []
    for i in range(n + 1):
        a = 2 * math.pi * i / n
        p = [math.cos(rr) * c[k] + math.sin(rr) * (math.cos(a) * north[k] + math.sin(a) * east[k]) for k in range(3)]
        pts.append((math.degrees(math.atan2(p[1], p[0])), math.degrees(math.asin(max(-1, min(1, p[2]))))))
    return pts


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6.5), gridspec_kw={"width_ratios": [1.25, 1]})
    for o in read(os.path.join(d, "obstacles.csv")):
        b = cap_boundary(float(o["lat"]), float(o["lon"]), float(o["radius_m"]))
        dyn = o["dynamic"] == "1"
        ax.fill([p[0] for p in b], [p[1] for p in b], color="tab:red" if dyn else "0.6", alpha=0.45)
    # 黑色实飞轨迹先画，绿色“突现后规划”最后画：两者事后几乎重合（跟踪误差 km 级），
    # 若绿色在下方会被黑色完全遮住
    styles = [("plan_initial.csv", "tab:blue", "-", 1.6, "initial plan"),
              ("plan_before_dynamic.csv", "tab:orange", "--", 1.6, "plan before new obstacle"),
              ("vehicle_track.csv", "k", "-", 1.0, "driven track"),
              ("plan_after_dynamic.csv", "tab:green", "-.", 1.3, "plan after new obstacle")]
    for fn, col, ls, lw, lab in styles:
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            continue
        rows = read(p)
        ax.plot([float(r["lon"]) for r in rows], [float(r["lat"]) for r in rows], color=col, lw=lw, ls=ls, label=lab)
    ax.set_aspect("equal")
    ax.set_ylim(-10, 60)
    ax.set_xlabel("lon [deg]")
    ax.set_ylabel("lat [deg]")
    ax.legend(loc="lower right")
    ax.set_title("Earth R=6371 km, no-fly zones (geodesic disks)")

    rows = read(os.path.join(d, "plan_initial.csv"))
    ax2.plot([float(r["s"]) / 1000.0 for r in rows], [float(r["kappa"]) * 1000.0 for r in rows],
             label="initial plan kappa_g")
    tr = read(os.path.join(d, "vehicle_track.csv"))
    # 累计飞行距离（t 轴 -> 航程轴），用于画随时间变化的转弯半径模型
    dist, xs, kc, kl = 0.0, [], [], []
    for i, r in enumerate(tr):
        if i > 0:
            dist += float(r["v"]) * (float(r["t"]) - float(tr[i - 1]["t"]))
        xs.append(dist / 1000.0)
        kc.append(float(r["kappa_cmd"]) * 1000.0)
        kl.append(float(r["kappa_lim"]) * 1000.0)
    klim0 = kl[0] if kl else 0.0
    ax2.axhline(klim0, color="r", ls=":", lw=1)
    ax2.axhline(-klim0, color="r", ls=":", lw=1, label="kappa_lim at entry")
    ax2.plot(xs, kl, color="r", ls="--", lw=1.0, label="kappa_lim(t) (linear turn radius)")
    ax2.plot(xs, [-x for x in kl], color="r", ls="--", lw=1.0)
    ax2.plot(xs, kc, color="k", lw=0.7, alpha=0.8, label="flown kappa_cmd")
    ax2b = ax2.twinx()
    ax2b.plot(xs, [float(r["solve_ms"]) for r in tr], color="0.7", lw=0.6,
              label="replan time [ms]")
    ax2.set_xlabel("arc length (plan s / flown distance) [km]")
    ax2.set_ylabel("geodesic curvature [1/km]")
    ax2b.set_ylabel("replan time [ms]")
    ax2.legend(loc="upper left")
    ax2b.legend(loc="upper right")
    fig.tight_layout()
    out = os.path.join(d, "demo.png")
    fig.savefig(out, dpi=110)
    print("saved", out)


if __name__ == "__main__":
    main()
