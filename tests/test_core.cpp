// 核心数学的自检：基函数导数、几何闭式公式、边界条件精确性、梯度、弧长互查、纯追踪公式。
#include <cstdio>
#include <random>

#include "sphtraj/sphere_planner.hpp"

using namespace sphtraj;

static int g_fail = 0;
#define CHECK(cond, ...)                         \
  do {                                           \
    if (!(cond)) {                               \
      ++g_fail;                                  \
      std::printf("[FAIL] " __VA_ARGS__);        \
      std::printf("   (%s:%d)\n", __FILE__, __LINE__); \
    } else {                                     \
      std::printf("[ OK ] " __VA_ARGS__);        \
      std::printf("\n");                         \
    }                                            \
  } while (0)

static void testQuadAndBasis() {
  std::vector<double> x, w;
  gaussLegendre01(5, x, w);
  double sw = 0, s4 = 0;
  for (int i = 0; i < 5; ++i) { sw += w[i]; s4 += w[i] * std::pow(x[i], 8); }
  CHECK(std::fabs(sw - 1) < 1e-14 && std::fabs(s4 - 1.0 / 9) < 1e-14, "Gauss-Legendre 权重和/9阶精确");

  ClampedCubicBasis B(7);
  double maxPU = 0, maxD1 = 0, maxD2 = 0;
  for (double u = 0.0; u <= 7.0; u += 0.0137) {
    double uu = u;
    const int k = B.spanOf(uu);
    double D[3][4];
    B.eval(uu, k, D);
    maxPU = std::max(maxPU, std::fabs(D[0][0] + D[0][1] + D[0][2] + D[0][3] - 1));
    const double h = 1e-6;
    if (u > h && u < 7 - h && std::fabs(u - std::round(u)) > 2 * h) {
      double Dp[3][4], Dm[3][4];
      B.eval(u + h, k, Dp);
      B.eval(u - h, k, Dm);
      for (int j = 0; j < 4; ++j) {
        maxD1 = std::max(maxD1, std::fabs((Dp[0][j] - Dm[0][j]) / (2 * h) - D[1][j]));
        maxD2 = std::max(maxD2, std::fabs((Dp[1][j] - Dm[1][j]) / (2 * h) - D[2][j]));
      }
    }
  }
  CHECK(maxPU < 1e-14 && maxD1 < 1e-7 && maxD2 < 1e-6, "B 样条单位分解与导数 (err %.1e %.1e %.1e)", maxPU, maxD1, maxD2);
}

static SphereSpline randomSpline(int M, unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> nd(0, 0.05);
  std::vector<Vec3> P;
  for (int j = 0; j < M + 3; ++j) {
    Vec3 p = fromLatLon(10.0 * j / (M + 2) + 0, 30.0 * j / (M + 2)) + Vec3(nd(rng), nd(rng), nd(rng));
    P.push_back(p * (1.0 + nd(rng)));
  }
  SphereSpline sp;
  sp.set(M, P, 2.0);
  return sp;
}

static void testGeometry() {
  const SphereSpline sp = randomSpline(6, 3);
  double eS = 0, eK = 0;
  for (double u = 0.31; u < 5.9; u += 0.37) {
    const auto g = sp.geom(u);
    const double h = 1e-4;
    const Vec3 a = sp.pos(u - h), b = sp.pos(u), c = sp.pos(u + h);
    const Vec3 d1 = (c - a) / (2 * h), d2 = (c - b * 2.0 + a) / (h * h);
    const double kNum = b.cross(d1).dot(d2) / std::pow(d1.norm(), 3);
    eS = std::max(eS, std::fabs(d1.norm() - g.speedU));
    eK = std::max(eK, std::fabs(kNum - g.kappaU) / std::max(1.0, std::fabs(g.kappaU)));
  }
  CHECK(eS < 1e-7 && eK < 1e-5, "闭式 |γ'| 与 κ_g 对比差分 (err %.1e %.1e)", eS, eK);

  double eInv = 0;
  for (double u = 0.0; u <= 6.0; u += 0.173) eInv = std::max(eInv, std::fabs(sp.uOfS(sp.sOfU(u)) - u));
  CHECK(eInv < 1e-9, "弧长<->参数互查 (err %.1e)", eInv);

  // 简单 Simpson 交叉验证总长
  double Ls = 0;
  const int n = 6000;
  for (int i = 0; i < n; ++i) {
    const double a = 6.0 * i / n, b = 6.0 * (i + 1) / n;
    Ls += (b - a) / 6 * (sp.speedU(a) + 4 * sp.speedU(0.5 * (a + b)) + sp.speedU(b));
  }
  CHECK(std::fabs(Ls * 2.0 - sp.length()) < 1e-9 * sp.length(), "总长 GL vs Simpson (rel %.1e)", std::fabs(Ls * 2.0 - sp.length()) / sp.length());

  const Vec3 target = (sp.pos(2.4) + Vec3(0, 0, 0.01)).normalized();
  const double uc = sp.closestU(target, 2.0, 1.0, 1.0);
  CHECK(std::fabs(sp.geom(uc).d1.dot(target)) < 1e-10, "最近点一阶条件 (u=%.4f)", uc);
}

static void testPurePursuit() {
  // 沿已知曲率的小圆运动后，公式应还原该曲率
  const double R = 50.0, kap = -0.07;
  const Vec3 p = fromLatLon(10, 20), h = headingFromAzimuth(p, 40);
  const double ku = kap * R, beta = std::atan2(1.0, ku);
  const Vec3 axis = p * std::cos(beta) + p.cross(h) * std::sin(beta);
  const Vec3 q = rotate(p, axis, (8.0 / R) / std::sin(beta));
  const double k2 = purePursuitCurvature(p, h, q, R);
  CHECK(std::fabs(k2 - kap) < 1e-12, "球面纯追踪曲率公式 (%.6f vs %.6f)", k2, kap);
}

static void testTrackingRefModes() {
  const SphereSpline sp = randomSpline(6, 3);
  const double R = sp.radius(), L = sp.length();
  const auto gEnd = sp.geom(sp.uMax());
  // 车辆取曲线上的点且航向与切向一致，最近点/弧长可精确对照
  const double u0 = 2.0;
  const Vec3 p = sp.pos(u0), h = sp.geom(u0).tangent;

  // 固定前视 (默认参数)：超出终点时沿终点切向大圆外延
  const TrackingRef rFix = computeTrackingRef(sp, p, h, u0, L);
  CHECK(std::fabs(rFix.lookaheadS - (rFix.s + L)) < 1e-12 && rFix.lookaheadS > L,
        "固定前视弧长外推 (s=%.3f, sLd=%.3f)", rFix.s, rFix.lookaheadS);
  CHECK(std::fabs(R * angleBetween(rFix.lookaheadPos, gEnd.pos) - (rFix.lookaheadS - L)) < 1e-9,
        "外延点沿终点切向大圆 (err %.1e)",
        std::fabs(R * angleBetween(rFix.lookaheadPos, gEnd.pos) - (rFix.lookaheadS - L)));

  // 比例前视远场：ratio*remaining > cap 时取 cap（与固定前视一致）
  const double ratio = 0.3, cap = 0.1 * L;
  const TrackingRef rFar = computeTrackingRef(sp, p, h, u0, cap, ratio);
  CHECK(ratio * rFar.remaining > cap && std::fabs(rFar.lookaheadS - (rFar.s + cap)) < 1e-9,
        "比例前视远场取上限 (rem=%.3f, cap=%.3f, sLd=%.3f)", rFar.remaining, cap, rFar.lookaheadS);

  // 比例前视近场：Ld = ratio*remaining < cap，前视点在曲线上且平滑收缩
  const double u1 = 5.5;
  const Vec3 p1 = sp.pos(u1), h1 = sp.geom(u1).tangent;
  const TrackingRef rNear = computeTrackingRef(sp, p1, h1, u1, cap, ratio);
  const double Ld1 = std::min(ratio * rNear.remaining, cap);
  CHECK(Ld1 < cap && std::fabs(rNear.lookaheadS - (rNear.s + Ld1)) < 1e-9 && rNear.lookaheadS <= L,
        "比例前视近场按剩余弧长收缩 (rem=%.3f, Ld=%.3f)", rNear.remaining, Ld1);
  CHECK(angleBetween(rNear.lookaheadPos, sp.pos(sp.uOfS(rNear.lookaheadS))) < 1e-12,
        "比例前视点落在曲线上");
}

static void testPlannerBoundaryAndGradient() {
  PlannerParams prm;
  prm.R = 1000.0;
  prm.kappaMax = 1.0 / 40.0;
  prm.spanLength = 50.0;
  prm.optMargin = 10.0;
  prm.checkMargin = 4.0;
  prm.pinStartCurvature = true;
  prm.wKappaVar = 1e3;
  SphereTrajPlanner pl(prm);
  std::vector<Obstacle> obs = {{fromLatLon(10, 12), deg2rad(5) * prm.R}, {fromLatLon(22, 26), deg2rad(4) * prm.R}};
  pl.setObstacles(obs);
  const Vec3 ps = fromLatLon(0, 0), pg = fromLatLon(30, 40);
  pl.setGoal({pg, headingFromAzimuth(pg, 90), 0.0});
  const SurfaceState s0{ps, headingFromAzimuth(ps, 30), 0.012};

  const double gerr = pl.gradientCheck(s0, 7);
  CHECK(gerr >= 0 && gerr < 1e-5, "解析梯度 vs 中心差分 (rel %.1e)", gerr);

  const PlanStats st = pl.planFromScratch(s0);
  CHECK(st.valid, "从零规划有效 (len %.1f, iters %d, %.2f ms)", st.length, st.iters, st.solveMs);
  const SphereSpline& sp = pl.curve();
  const auto g0 = sp.geom(0.0), g1 = sp.geom(sp.uMax());
  const double hs = angleBetween(g0.tangent, s0.h);
  CHECK(hs < 1e-10, "起点航向精确 (%.1e rad)", hs);
  const double k0 = g0.kappaU / prm.R;
  CHECK(std::fabs(k0 - 0.012) < 1e-10, "起点曲率精确 (%.6f)", k0);
  CHECK(angleBetween(g0.pos, ps) < 1e-14 && angleBetween(g1.pos, pg) < 1e-14, "起终点位置精确");
  // 终端航向不约束：终点切向不应被钉死到 setGoal 传入的航向
  CHECK(angleBetween(g1.tangent, headingFromAzimuth(pg, 90)) > 0.05,
        "终点切向不钉死到 setGoal 航向 (%.2f rad)", angleBetween(g1.tangent, headingFromAzimuth(pg, 90)));
}

static void testBlockedReplan() {
  PlannerParams prm;
  prm.R = 1000.0;
  prm.kappaMax = 1.0 / 40.0;
  prm.spanLength = 50.0;
  prm.optMargin = 10.0;
  prm.checkMargin = 4.0;
  SphereTrajPlanner pl(prm);
  pl.setObstacles({{fromLatLon(10, 12), deg2rad(5) * prm.R}});
  const Vec3 ps = fromLatLon(0, 0), pg = fromLatLon(25, 30);
  pl.setGoal({pg, headingFromAzimuth(pg, 60), 0.0});
  const SurfaceState s0{ps, headingFromAzimuth(ps, 45), 0.0};
  PlanStats st = pl.planFromScratch(s0);
  CHECK(st.valid, "阻塞测试：初始规划有效");
  // 在路径中段放置一个大障碍，完全挡住当前路径
  const SphereSpline& c = pl.curve();
  pl.addObstacle({c.pos(c.uOfS(0.55 * c.length())), 120.0});
  st = pl.replan(s0, 0.0);
  const auto v = pl.validate(pl.curve(), 0.0);
  CHECK(v.ok && (st.status == PlanStatus::Ok || st.status == PlanStatus::FrontEndReplan),
        "阻塞后重规划 (无同伦守卫, 可行即接受)：status=%s, 新路径有效 (clr %.1f m)",
        toString(st.status), v.minClearance);
}

static void testRetainedStatsAndEmptyCurve() {
  PlannerParams prm;
  prm.R = 1000.0;
  prm.kappaMax = 1.0 / 40.0;
  prm.spanLength = 50.0;
  prm.optMargin = 10.0;
  prm.checkMargin = 4.0;
  const Vec3 ps = fromLatLon(0, 0), pg = fromLatLon(25, 30);
  const SurfaceState s0{ps, headingFromAzimuth(ps, 45), 0.0};
  SphereTrajPlanner pl(prm);
  pl.setGoal({pg, headingFromAzimuth(pg, 60), 0.0});
  PlanStats st = pl.planFromScratch(s0);
  CHECK(st.valid && pl.hasPlan(), "保留测试：初始规划有效");
  const double len0 = pl.curve().length();
  const int spans0 = pl.curve().numSpans();
  const Vec3 end0 = pl.curve().pos(pl.curve().uMax());
  // 新终点放在障碍中心：任何到达它的曲线都穿障。旧曲线远离该点，应被留下。
  const Vec3 bad = fromLatLon(40, 50);
  pl.addObstacle({bad, 80.0});
  pl.setGoal({bad, headingFromAzimuth(bad, 60), 0.0});
  st = pl.replan(s0, 0.0);
  CHECK(st.status == PlanStatus::KeptPreviousInfeasible, "新解穿障时保留旧解 status=%s", toString(st.status));
  CHECK(st.valid, "保留的旧解仍然有效");
  CHECK(std::fabs(st.length - len0) < 1e-6 && st.spans == spans0, "统计量属于保留曲线");
  CHECK(st.minClearance >= 0.0, "保留曲线间隔 %.3f m", st.minClearance);
  CHECK(angleBetween(pl.curve().pos(pl.curve().uMax()), end0) < 1e-9, "曲线终点未被坏目标替换");

  PlannerParams badP = prm;
  badP.checkMargin = 1.0e9;
  SphereTrajPlanner plBad(badP);
  plBad.setObstacles({{fromLatLon(10, 10), 50.0}});
  plBad.setGoal({pg, headingFromAzimuth(pg, 60), 0.0});
  st = plBad.planFromScratch(s0);
  CHECK(!st.valid && !plBad.hasPlan() && st.spans == 0 && st.length == 0.0, "失败的从零规划不留下曲线");
  const TrackingRef tr = plBad.trackingReference(s0, 10.0);
  CHECK(tr.u == 0.0 && tr.kappaPurePursuit == 0.0, "无曲线时跟踪返回空参考");
  SphereSpline empty;
  const auto g = empty.geom(0.3);
  CHECK(g.cnorm == 0.0 && empty.uOfS(1.0) == 0.0 && empty.closestU(ps, 0.2) == 0.0, "空样条几何查询返回空值");
  CHECK(empty.samplePositions(0.0, 0.1).empty(), "空样条采样为空");
  const auto badVal = plBad.validate(empty, 0.0);
  CHECK(!badVal.ok, "空曲线校验为无效");
}

// 曲率变化率软罚：解析梯度、有效性、无硬上限、平滑效果。
// 曲率总变化 TV = Σ|Δκ| 沿等参采样序列累计。
static double kappaTotalVariation(const SphereSpline& c, int n, double R) {
  double tv = 0, kapPrev = 0;
  for (int i = 0; i <= n; ++i) {
    const double u = c.uMax() * i / n;
    const double kap = c.geom(u).kappaU / R;
    if (i > 0) tv += std::fabs(kap - kapPrev);
    kapPrev = kap;
  }
  return tv;
}

static void testCurvatureVariation() {
  const double R = 6371000.0;
  const Vec3 ps = fromLatLon(15, 8), pg = fromLatLon(27, 100);
  const SurfaceState s0{ps, headingFromAzimuth(ps, 60.0), 0.0};
  auto cap = [&](double lat, double lon, double rd) { return Obstacle{fromLatLon(lat, lon), deg2rad(rd) * R}; };
  const std::vector<Obstacle> obs = {cap(27, 37, 3.1), cap(30, 58, 3.5), cap(21, 76, 3.2)};

  auto make = [&](double wKv) {
    PlannerParams p;
    p.R = R;
    p.kappaMax = 1.0 / 1.0e6;
    p.spanLength = 600e3;
    p.optMargin = 80e3;
    p.checkMargin = 25e3;
    p.wKappaVar = wKv;
    return p;
  };

  // 梯度（含变化罚项）
  PlannerParams gp = make(1e3);
  SphereTrajPlanner plg(gp);
  plg.setObstacles(obs);
  plg.setGoal({pg, headingFromAzimuth(pg, 70.0), 0.0});
  const double gerr = plg.gradientCheck(s0, 11);
  CHECK(gerr >= 0 && gerr < 1e-5, "曲率变化罚：解析梯度 vs 中心差分 (rel %.1e)", gerr);

  auto run = [&](double wKv) {
    SphereTrajPlanner pl(make(wKv));
    pl.setObstacles(obs);
    pl.setGoal({pg, headingFromAzimuth(pg, 70.0), 0.0});
    const PlanStats st = pl.planFromScratch(s0);
    const double tv = st.valid ? kappaTotalVariation(pl.curve(), 200, R) : -1.0;
    return std::make_pair(st, tv);  // (stats, 曲率总变化)
  };
  auto [s0r, tv0] = run(0.0);
  auto [sKr, tvK] = run(1e3);
  CHECK(s0r.valid && sKr.valid, "曲率变化罚：开/关均有效 (len %.0f / %.0f km)", s0r.length / 1e3, sKr.length / 1e3);
  CHECK(sKr.minClearance > 0.0 && s0r.minClearance > 0.0, "曲率变化罚：障碍间隔为正 (%.1f / %.1f km)",
        sKr.minClearance / 1e3, s0r.minClearance / 1e3);
  CHECK(tvK < tv0, "曲率变化罚：开启后曲率总变化下降 (TV %.3e -> %.3e, len +%.1f%%)",
        tv0, tvK, 100.0 * (sKr.length / s0r.length - 1.0));

  // 无硬上限：把 kappaMax 收紧到远小于绕障所需曲率, 校验仍通过且 maxKappa 超过参考值
  PlannerParams tp = make(1e3);
  tp.kappaMax = 1.0 / 8.0e6;
  tp.leadLength = 600e3;  // 避免 leadLength 自动值 = 1/kappaMax 过大
  SphereTrajPlanner plt(tp);
  plt.setObstacles(obs);
  plt.setGoal({pg, headingFromAzimuth(pg, 70.0), 0.0});
  const PlanStats stT = plt.planFromScratch(s0);
  CHECK(stT.valid && stT.maxKappa > tp.kappaMax,
        "无曲率硬上限：收紧参考曲率后仍有效 (maxKappa/kappaMax = %.2f)", stT.maxKappa / tp.kappaMax);
}

int main() {
  testQuadAndBasis();
  testGeometry();
  testPurePursuit();
  testTrackingRefModes();
  testPlannerBoundaryAndGradient();
  testBlockedReplan();
  testRetainedStatsAndEmptyCurve();
  testCurvatureVariation();
  std::printf("\n%s (%d failures)\n", g_fail ? "SOME TESTS FAILED" : "ALL TESTS PASSED", g_fail);
  return g_fail ? 1 : 0;
}
