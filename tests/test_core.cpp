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

static void testPlannerBoundaryAndGradient() {
  PlannerParams prm;
  prm.R = 1000.0;
  prm.kappaMax = 1.0 / 40.0;
  prm.spanLength = 50.0;
  prm.optMargin = 10.0;
  prm.checkMargin = 4.0;
  prm.pinStartCurvature = true;
  prm.pinGoalCurvature = true;
  prm.goalCurvature = -0.004;
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
  const double hs = angleBetween(g0.tangent, s0.h), hg = angleBetween(g1.tangent, headingFromAzimuth(pg, 90));
  CHECK(hs < 1e-10 && hg < 1e-10, "起终点航向精确 (%.1e, %.1e rad)", hs, hg);
  const double k0 = g0.kappaU / prm.R, k1 = g1.kappaU / prm.R;
  CHECK(std::fabs(k0 - 0.012) < 1e-10 && std::fabs(k1 + 0.004) < 1e-10, "起终点曲率精确 (%.6f, %.6f)", k0, k1);
  CHECK(angleBetween(g0.pos, ps) < 1e-14 && angleBetween(g1.pos, pg) < 1e-14, "起终点位置精确");
}

static std::vector<Vec3> densify(const std::vector<Vec3>& pts, int n) {
  std::vector<Vec3> out;
  for (size_t i = 1; i < pts.size(); ++i)
    for (int k = 0; k < n; ++k) out.push_back(slerp(pts[i - 1], pts[i], double(k) / n));
  out.push_back(pts.back());
  return out;
}

static void testHomotopySignature() {
  PlannerParams prm;
  prm.R = 1.0;
  prm.kappaMax = 10.0;
  prm.spanLength = 0.05;
  SphereTrajPlanner pl(prm);
  pl.setObstacles({{fromLatLon(0, 10), deg2rad(2)}, {fromLatLon(0, 30), deg2rad(2)}});
  pl.setHomotopyPole(fromLatLon(-10, -160));
  const Vec3 A = fromLatLon(0, 0), B = fromLatLon(0, 20);
  const auto north = pl.signatureOfPoints(densify({A, fromLatLon(5, 10), B}, 50));
  const auto north2 = pl.signatureOfPoints(densify({A, fromLatLon(8, 7), fromLatLon(3, 14), B}, 50));
  const auto south = pl.signatureOfPoints(densify({A, fromLatLon(-5, 10), B}, 50));
  const double dSame = std::fabs(north[0] - north2[0]), dDiff = std::fabs(north[0] - south[0]);
  CHECK(dSame < 0.5 && std::fabs(dDiff - 2 * kPi) < 1e-6, "同伦签名：同侧差 %.3f，异侧差 %.6f (应为 2π)", dSame, dDiff);
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
  CHECK(v.ok && (st.status == PlanStatus::ClassChangedForced || st.status == PlanStatus::FrontEndReplan ||
                 st.status == PlanStatus::Ok),
        "阻塞后重规划：status=%s, 新路径有效 (clr %.1f m)", toString(st.status), v.minClearance);
}

int main() {
  testQuadAndBasis();
  testGeometry();
  testPurePursuit();
  testPlannerBoundaryAndGradient();
  testHomotopySignature();
  testBlockedReplan();
  std::printf("\n%s (%d failures)\n", g_fail ? "SOME TESTS FAILED" : "ALL TESTS PASSED", g_fail);
  return g_fail ? 1 : 0;
}
