// 闭环演示：从零规划 → 每个控制周期重规划 + 球面纯追踪跟踪 → 中途出现新障碍。
// 输出统计信息，并写出 CSV（可用 tools/plot_demo.py 画图）。
#include <algorithm>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "sphtraj/sphere_planner.hpp"

using namespace sphtraj;

// 以物理曲率 kappa 在球面上前进 ds（沿小圆精确积分）
static void moveVehicle(Vec3& p, Vec3& h, double kappa, double ds, double R) {
  const double ku = kappa * R;
  const double beta = std::atan2(1.0, ku);  // cot(beta) = ku，beta∈(0,π)
  const Vec3 axis = p * std::cos(beta) + p.cross(h) * std::sin(beta);
  const double th = (ds / R) / std::sin(beta);
  p = rotate(p, axis, th).normalized();
  h = projectTangent(rotate(h, axis, th), p).normalized();
}

static void writeCurve(const std::string& file, const SphereSpline& sp) {
  FILE* f = std::fopen(file.c_str(), "w");
  std::fprintf(f, "s,lat,lon,kappa\n");
  const int n = 400;
  for (int i = 0; i <= n; ++i) {
    const double s = sp.length() * i / n;
    const auto g = sp.geom(sp.uOfS(s));
    double la, lo;
    latLonOf(g.pos, la, lo);
    std::fprintf(f, "%.4f,%.8f,%.8f,%.8f\n", s, la, lo, g.kappaU / sp.radius());
  }
  std::fclose(f);
}

int main(int argc, char** argv) {
  const std::string outDir = argc > 1 ? argv[1] : ".";
  // ---------------- 场景 ----------------
  PlannerParams prm;
  prm.R = 1000.0;               // 球半径 1000 m
  prm.kappaMax = 1.0 / 40.0;    // 最小转弯半径 40 m
  prm.spanLength = 50.0;        // 每个 span 约 50 m
  prm.optMargin = 10.0;         // 优化膨胀 10 m
  prm.checkMargin = 4.0;        // 校验膨胀 4 m（给跟踪误差留余量）
  prm.pinStartCurvature = true;   // 起点曲率 = 车辆当前曲率（转向指令跨周期连续）
  prm.pinGoalCurvature = true;    // 以零曲率（直线）到达终点
  prm.goalCurvature = 0.0;
  SphereTrajPlanner planner(prm);

  const double R = prm.R;
  auto cap = [&](double lat, double lon, double rDeg) { return Obstacle{fromLatLon(lat, lon), deg2rad(rDeg) * R}; };
  std::vector<Obstacle> obs = {cap(8, 10, 5),   cap(17, 22, 4.5), cap(24, 31, 5), cap(13, 29, 4),
                               cap(27, 16, 4),  cap(32, 40, 3),   cap(5, 24, 3.5), cap(20, 40, 3)};
  planner.setObstacles(obs);
  const Vec3 ps = fromLatLon(0, 0), pg = fromLatLon(35, 45);
  const SurfaceState start{ps, headingFromAzimuth(ps, 30), 0.0};
  planner.setGoal({pg, headingFromAzimuth(pg, 90), 0.0});

  std::printf("great-circle distance start->goal: %.1f m\n", R * angleBetween(ps, pg));
  std::printf("gradient check (rel err): %.2e\n", planner.gradientCheck(start));

  // ---------------- 从零规划 ----------------
  PlanStats st = planner.planFromScratch(start);
  std::printf("[scratch] valid=%d  len=%.1f m  spans=%d  iters=%d  evals=%d  minClr=%.1f m  maxK=%.4f (lim %.4f)  %.2f ms\n",
              st.valid, st.length, st.spans, st.iters, st.evals, st.minClearance, st.maxKappa, prm.kappaMax,
              st.solveMs);
  if (!st.valid) return 1;
  writeCurve(outDir + "/plan_initial.csv", planner.curve());
  {
    FILE* f = std::fopen((outDir + "/obstacles.csv").c_str(), "w");
    std::fprintf(f, "lat,lon,radius_m,dynamic\n");
    for (const auto& o : obs) {
      double la, lo;
      latLonOf(o.center, la, lo);
      std::fprintf(f, "%.8f,%.8f,%.4f,0\n", la, lo, o.radius);
    }
    std::fclose(f);
  }

  // ---------------- 闭环仿真 ----------------
  const double v = 5.0, dt = 0.05, lookahead = 20.0, drift = 0.25;  // 速度 m/s, 周期 s, 前视 m, 侧向漂移 m/s
  std::mt19937 rng(42);
  std::normal_distribution<double> headNoise(0.0, deg2rad(0.3));
  Vec3 p = start.p, h = start.h;
  double kCmd = 0.0, prevKCmd = 0.0;
  bool dynAdded = false;
  std::map<std::string, int> statusCount;
  std::vector<double> times;
  long iterSum = 0, iterCnt = 0;
  double maxCross = 0, maxDK = 0, minClrDriven = 1e9, devAtEvent = 0;
  std::vector<double> devs;  // 常规周期的相邻规划偏差
  int eventCycle = -1;
  int classSwitches = 0;
  FILE* ft = std::fopen((outDir + "/vehicle_track.csv").c_str(), "w");
  std::fprintf(ft, "t,lat,lon,cross_track,kappa_cmd,solve_ms,status\n");

  int cycle = 0;
  for (double t = 0.0; t < 600.0; t += dt, ++cycle) {
    const SurfaceState meas{p, h, kCmd};

    // 中途出现新障碍：放在当前路径前方 150 m 处
    if (!dynAdded && t >= 60.0) {
      const SphereSpline& c = planner.curve();
      const double sNow = c.sOfU(c.closestU(p, planner.lastU(), 0.5, 3.0));
      const auto g = c.geom(c.uOfS(sNow + 150.0));
      const Vec3 center = geodesicStep(g.pos, g.left, 8.0 / R);
      planner.addObstacle({center, 30.0});
      FILE* f = std::fopen((outDir + "/obstacles.csv").c_str(), "a");
      double la, lo;
      latLonOf(center, la, lo);
      std::fprintf(f, "%.8f,%.8f,%.4f,1\n", la, lo, 30.0);
      std::fclose(f);
      writeCurve(outDir + "/plan_before_dynamic.csv", c);
      dynAdded = true;
      eventCycle = cycle;
      std::printf("[t=%.1f] dynamic obstacle added 150 m ahead\n", t);
    }

    // 记录重规划前旧轨迹的前 200 m，用于衡量相邻周期一致性
    const SphereSpline oldCurve = planner.curve();
    const double sOld = oldCurve.sOfU(oldCurve.closestU(p, planner.lastU(), 0.5, 3.0));

    // 每周期重规划
    st = planner.replan(meas, 0.0);
    statusCount[toString(st.status)]++;
    times.push_back(st.solveMs);
    if (st.status != PlanStatus::Frozen) { iterSum += st.iters; ++iterCnt; }
    if (st.status == PlanStatus::ClassChangedForced || st.status == PlanStatus::FrontEndReplan)
      writeCurve(outDir + "/plan_after_dynamic.csv", planner.curve());

    // 低频（1 Hz）检查是否值得主动切换同伦类
    if (cycle % 20 == 0 && planner.considerClassSwitch(meas)) ++classSwitches;

    // 相邻周期偏差（仅统计“同类且已接受”的周期）
    if (st.status == PlanStatus::Ok) {
      const SphereSpline& nc = planner.curve();
      double dev = 0;
      for (double d = 0; d <= 200.0; d += 20.0) {
        if (sOld + d > oldCurve.length() || d > nc.length()) break;
        const Vec3 a = oldCurve.pos(oldCurve.uOfS(sOld + d)), b = nc.pos(nc.uOfS(d));
        dev = std::max(dev, R * angleBetween(a, b));
      }
      if (cycle == eventCycle) devAtEvent = dev;
      else devs.push_back(dev);
    }

    // 跟踪：球面纯追踪
    const TrackingRef ref = planner.trackingReference(meas, lookahead);
    kCmd = clampd(ref.kappaPurePursuit, -prm.kappaMax, prm.kappaMax);
    if (cycle > 0) maxDK = std::max(maxDK, std::fabs(kCmd - prevKCmd));
    prevKCmd = kCmd;
    maxCross = std::max(maxCross, std::fabs(ref.crossTrack));
    for (const auto& o : planner.obstacles())
      minClrDriven = std::min(minClrDriven, R * angleBetween(p, o.center) - o.radius);

    double la, lo;
    latLonOf(p, la, lo);
    std::fprintf(ft, "%.3f,%.8f,%.8f,%.4f,%.6f,%.4f,%s\n", t, la, lo, ref.crossTrack, kCmd, st.solveMs,
                 toString(st.status));

    if (ref.remaining < 0.5) {
      std::printf("[t=%.1f] goal reached: pos err %.2f m, heading err %.2f deg\n", t, R * angleBetween(p, pg),
                  rad2deg(angleBetween(h, headingFromAzimuth(pg, 90))));
      break;
    }
    // 车辆运动 + 扰动
    moveVehicle(p, h, kCmd, v * dt, R);
    h = rotate(h, p, headNoise(rng));
    p = geodesicStep(p, p.cross(h), drift * dt / R);
    h = projectTangent(h, p).normalized();
  }
  std::fclose(ft);
  writeCurve(outDir + "/plan_final.csv", planner.curve());

  std::sort(times.begin(), times.end());
  double mean = 0;
  for (double x : times) mean += x;
  mean /= times.size();
  std::printf("\nreplan cycles: %zu   solve time: mean %.3f ms, p50 %.3f, p99 %.3f, max %.3f\n", times.size(), mean,
              times[times.size() / 2], times[times.size() * 99 / 100], times.back());
  std::printf("mean L-BFGS iterations per replan: %.1f\n", double(iterSum) / std::max(1L, iterCnt));
  for (const auto& kv : statusCount) std::printf("  status %-28s %d\n", kv.first.c_str(), kv.second);
  std::printf("voluntary class switches: %d\n", classSwitches);
  std::printf("max |cross-track| = %.2f m, max |dKappa_cmd| per cycle = %.5f 1/m\n", maxCross, maxDK);
  std::sort(devs.begin(), devs.end());
  std::printf("consecutive-plan deviation over first 200 m: p50 %.3f m, p99 %.3f m, max %.3f m (event cycle: %.2f m)\n",
              devs[devs.size() / 2], devs[devs.size() * 99 / 100], devs.back(), devAtEvent);
  std::printf("min clearance of driven track to obstacle boundary = %.2f m\n", minClrDriven);
  return 0;
}
