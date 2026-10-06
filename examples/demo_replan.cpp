// 闭环演示：高超声速滑翔再入场景（地球尺度）。
// 从再入界面从零规划 → 每个制导周期（1 Hz）重规划 + 球面纯追踪跟踪 → 中途出现新禁飞区。
// 速度模型：从再入界面 7.0 km/s 线性衰减到末段 2.0 km/s（历时 2500 s），之后保持。
// 转弯半径模型（简化）：等侧向过载下 r = v²/(n·g) 随速度平方减小，这里简化为与速度同相位的线性减小，
// 即 kappaMax(t) = 1/rMin(t) 随时间单调放宽，机动能力随减速增强。
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
  // ---------------- 场景：地球尺度滑翔再入 ----------------
  // 速度与转弯半径模型（线性）：v 从 7000 → 2000 m/s（2500 s）；转弯半径与速度同相位线性减小
  // 1200 km → 650 km（等效侧向过载 v²/(r·g)：再入界面 ~4 g，末段 ~2.7 g），即 κmax 从 8.3e-7 → 6.7e-6 1/m。
  const double vEntry = 7000.0, vTerm = 2000.0, tDecay = 2500.0;
  const double rMinEntry = 1200000.0, rMinTerm = 650000.0;
  auto speedAt = [&](double t) { return vEntry + (vTerm - vEntry) * std::min(t / tDecay, 1.0); };
  auto rMinAt = [&](double t) { return rMinEntry + (rMinTerm - rMinEntry) * std::min(t / tDecay, 1.0); };
  PlannerParams prm;
  prm.R = 6371000.0;             // 地球平均半径
  prm.kappaMax = 1.0 / rMinEntry; // 初始转弯半径 1200 km，之后每周期按线性减小模型更新（setKappaMax）
  prm.spanLength = 300000.0;     // 每个 span 约 300 km（≈ 末段转弯半径的 2 倍，span 角约 2.7°）
  prm.optMargin = 60000.0;       // 优化膨胀 60 km
  prm.checkMargin = 25000.0;     // 校验膨胀 25 km（给跟踪误差留余量）
  prm.pinStartCurvature = true;   // 起点曲率 = 飞行器当前曲率（转向指令跨周期连续）
  prm.pinGoalCurvature = true;    // 以零曲率（直线）进入末段交接
  prm.goalLeadLength = 0.6 * prm.spanLength;  // 末段切向臂约 spanLength/2，末段“近似直线”长度随之增长
  prm.goalCurvature = 0.0;
  SphereTrajPlanner planner(prm);

  const double R = prm.R;
  auto cap = [&](double lat, double lon, double rDeg) { return Obstacle{fromLatLon(lat, lon), deg2rad(rDeg) * R}; };
  // 静态禁飞区：角半径 3°~5°，即物理半径约 334~556 km（再入路径规划的典型量级）
  std::vector<Obstacle> obs = {cap(27, 37, 3.1),   cap(30, 58, 3.5), cap(21, 76, 3.2), cap(18, 52, 3.0)};
  planner.setObstacles(obs);
  // 再入界面 (15,8) → 末段目标 (27,100)：大圆航程约 10000 km（洲际滑翔的典型量级）
  const Vec3 ps = fromLatLon(15, 8), pg = fromLatLon(27, 100);
  const SurfaceState start{ps, headingFromAzimuth(ps, 60), 0.0};
  planner.setGoal({pg, headingFromAzimuth(pg, 70), 0.0});

  std::printf("great-circle distance start->goal: %.1f km\n", R * angleBetween(ps, pg) / 1000.0);
  std::printf("gradient check (rel err): %.2e\n", planner.gradientCheck(start));

  // ---------------- 从零规划 ----------------
  PlanStats st = planner.planFromScratch(start);
  std::printf("[scratch] valid=%d  len=%.1f km  spans=%d  iters=%d  evals=%d  minClr=%.1f km  maxK=%.3e (lim %.3e)  %.2f ms\n",
              st.valid, st.length / 1000.0, st.spans, st.iters, st.evals, st.minClearance / 1000.0, st.maxKappa,
              prm.kappaMax, st.solveMs);
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
  // 5 s 制导/重规划周期，前视 125 km，侧向漂移 30 m/s（急流量级），航向噪声 1°/步。
  // 0.2 s 计算/通信延迟用 advanceDist = v*0.2 补偿。转弯半径每周期按线性减小模型收紧 planner 限幅。
  const double dt = 5.0, lookahead = 125000.0, drift = 30.0;
  std::mt19937 rng(42);
  std::normal_distribution<double> headNoise(0.0, deg2rad(1.0));
  Vec3 p = start.p, h = start.h;
  double kCmd = 0.0, prevKCmd = 0.0;
  bool dynAdded = false;
  std::map<std::string, int> statusCount;
  std::vector<double> times;
  long iterSum = 0, iterCnt = 0;
  double maxCross = 0, maxDK = 0, minClrDriven = 1e9, devAtEvent = 0, driven = 0, kapUtil = 0;
  std::vector<double> devs;  // 常规周期的相邻规划偏差
  int eventCycle = -1;
  int classSwitches = 0;
  FILE* ft = std::fopen((outDir + "/vehicle_track.csv").c_str(), "w");
  std::fprintf(ft, "t,lat,lon,s,cross_track,v,kappa_cmd,kappa_lim,r_min_m,solve_ms,status\n");

  int cycle = 0;
  for (double t = 0.0; t < 2500.0; t += dt, ++cycle) {
    const double v = speedAt(t);
    // 转弯半径线性减小：更新 planner/跟踪的曲率限幅（限幅单调放宽，旧轨迹始终可行）
    const double rMin = rMinAt(t), kapNow = 1.0 / rMin;
    planner.setKappaMax(kapNow);
    const SurfaceState meas{p, h, kCmd};

    // 中途出现新禁飞区：放在当前路径前方 1500 km 处
    if (!dynAdded && t >= 600.0) {
      const SphereSpline& c = planner.curve();
      const double sNow = c.sOfU(c.closestU(p, planner.lastU(), 0.5, 3.0));
      const auto g = c.geom(c.uOfS(sNow + 1.5e6));
      const Vec3 center = geodesicStep(g.pos, g.left, 150000.0 / R);
      planner.addObstacle({center, 280000.0});
      FILE* f = std::fopen((outDir + "/obstacles.csv").c_str(), "a");
      double la, lo;
      latLonOf(center, la, lo);
      std::fprintf(f, "%.8f,%.8f,%.4f,1\n", la, lo, 280000.0);
      std::fclose(f);
      writeCurve(outDir + "/plan_before_dynamic.csv", c);
      dynAdded = true;
      eventCycle = cycle;
      std::printf("[t=%.1f] dynamic no-fly zone added 1500 km ahead\n", t);
    }

    // 记录重规划前旧轨迹的前 1000 km，用于衡量相邻周期一致性
    const SphereSpline oldCurve = planner.curve();
    const double sOld = oldCurve.sOfU(oldCurve.closestU(p, planner.lastU(), 0.5, 3.0));

    // 每周期重规划（advanceDist 补偿计算/通信延迟）
    st = planner.replan(meas, v * 0.2);
    statusCount[toString(st.status)]++;
    times.push_back(st.solveMs);
    if (st.status != PlanStatus::Frozen) { iterSum += st.iters; ++iterCnt; }
    // 突现禁飞区当周期无论状态如何（Ok/ClassChangedForced/FrontEndReplan）都写出新规划
    if (cycle == eventCycle || st.status == PlanStatus::ClassChangedForced || st.status == PlanStatus::FrontEndReplan)
      writeCurve(outDir + "/plan_after_dynamic.csv", planner.curve());

    // 低频（0.2 Hz）检查是否值得主动切换同伦类
    if (cycle % 5 == 0 && planner.considerClassSwitch(meas)) ++classSwitches;

    // 相邻周期偏差（仅统计“同类且已接受”的周期）
    if (st.status == PlanStatus::Ok) {
      const SphereSpline& nc = planner.curve();
      double dev = 0;
      for (double d = 0; d <= 1.0e6; d += 1.0e5) {
        if (sOld + d > oldCurve.length() || d > nc.length()) break;
        const Vec3 a = oldCurve.pos(oldCurve.uOfS(sOld + d)), b = nc.pos(nc.uOfS(d));
        dev = std::max(dev, R * angleBetween(a, b));
      }
      if (cycle == eventCycle) devAtEvent = dev;
      else devs.push_back(dev);
    }

    // 跟踪：球面纯追踪（限幅取当前转弯半径模型）
    const TrackingRef ref = planner.trackingReference(meas, lookahead);
    kCmd = clampd(ref.kappaPurePursuit, -kapNow, kapNow);
    if (cycle > 0) maxDK = std::max(maxDK, std::fabs(kCmd - prevKCmd));
    prevKCmd = kCmd;
    kapUtil = std::max(kapUtil, std::fabs(kCmd) / kapNow);
    maxCross = std::max(maxCross, std::fabs(ref.crossTrack));
    for (const auto& o : planner.obstacles())
      minClrDriven = std::min(minClrDriven, R * angleBetween(p, o.center) - o.radius);

    double la, lo;
    latLonOf(p, la, lo);
    std::fprintf(ft, "%.3f,%.8f,%.8f,%.2f,%.4f,%.2f,%.6e,%.6e,%.0f,%.4f,%s\n", t, la, lo, ref.s, ref.crossTrack, v,
                 kCmd, kapNow, rMin, st.solveMs, toString(st.status));

    if (ref.remaining < 3000.0) {
      std::printf("[t=%.1f] goal reached: pos err %.3f km, heading err %.2f deg, final speed %.1f m/s, r_min %.0f km, flown %.1f km\n",
                  t, R * angleBetween(p, pg) / 1000.0,
                  rad2deg(angleBetween(h, headingFromAzimuth(pg, 90))), v, rMin / 1000.0, driven / 1000.0);
      break;
    }
    // 飞行器运动 + 扰动
    moveVehicle(p, h, kCmd, v * dt, R);
    driven += v * dt;
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
  std::printf("turn-radius model: %.0f -> %.0f km linear (kappaLim %.3e -> %.3e 1/m), max kappa_cmd utilization = %.0f%%\n",
              rMinEntry / 1000.0, rMinTerm / 1000.0, 1.0 / rMinEntry, 1.0 / rMinTerm, 100.0 * kapUtil);
  std::printf("max |cross-track| = %.2f km, max |dKappa_cmd| per cycle = %.3e 1/m\n", maxCross / 1000.0, maxDK);
  std::sort(devs.begin(), devs.end());
  std::printf("consecutive-plan deviation over first 1000 km: p50 %.3f km, p99 %.3f km, max %.3f km (event cycle: %.2f km)\n",
              devs[devs.size() / 2] / 1000.0, devs[devs.size() * 99 / 100] / 1000.0, devs.back() / 1000.0,
              devAtEvent / 1000.0);
  std::printf("min clearance of flown track to no-fly-zone boundary = %.2f km\n", minClrDriven / 1000.0);
  return 0;
}
