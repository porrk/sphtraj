#pragma once
// 方案 A：球面中心投影 B 样条轨迹优化 + 切线图前端 + 每周期重规划（同伦一致性）。
//
// 曲线：c(u) = Σ N_j(u) P_j（R^3 夹持三次 B 样条），γ = c/|c|。
// 边界条件（构造式精确满足）：
//   P0 = ps, P1 = ps + a_s hs                  → 起点航向精确
//   P_{M+2} = pg, P_{M+1} = pg - a_g hg         → 终点航向精确
//   ℓs·P2 = 3 a_s² κ0  (ℓs = ps × hs)           → 起点曲率 = 当前曲率（可选，保证转向指令连续）
//   ℓg·P_M = 3 a_g² κg (ℓg = pg × hg)           → 终点曲率（可选）
// 代价：J = wE·E + wB·B + wK·K + wO·O + wR·Rad + wP·Prox
//   E   = M ∫|γ'|² du / Lref²         （≥ (L/Lref)²，等号当且仅当匀速 → 最短且近似弧长参数化）
//   B   = ∫ κ² ds / (κlim² Lref)      （弯曲能量正则）
//   K   = (1/M) ∫ [max(0, |κ|/κlim - 1)³ + max(0, 1 - |γ'|/smin)³] du   （曲率上限 + 速度下限防尖点）
//   O   = (1/M) ∫ Σ_i max(0, (γ·n_i - cos r_i)/(sin r_i σ))³ du
//   Rad = Σ (|P_j|² - 1)²             （消除径向冗余自由度）
//   Prox= Σ |P_j - P_j^warm|² / Δ²    （近端项：相邻周期一致性）
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include "bspline_basis.hpp"
#include "lbfgs.hpp"
#include "sphere_math.hpp"
#include "sphere_spline.hpp"
#include "tangent_graph.hpp"

namespace sphtraj {

// 车辆/边界状态：p 单位位置向量，h 单位切向航向，kappa 物理测地曲率（1/长度，左转为正）
struct SurfaceState {
  Vec3 p;
  Vec3 h;
  double kappa = 0.0;
};

// 障碍：中心单位向量 + 物理测地半径（应已包含车体尺寸）
struct Obstacle {
  Vec3 center;
  double radius = 0.0;
};

struct PlannerParams {
  // ---- 几何 ----
  double R = 1.0;             // 球半径（物理单位）
  double kappaMax = 1.0;      // 最大测地曲率（1/物理长度）
  double spanLength = 0.1;    // 每个 knot span 期望弧长（物理）
  int minSpans = 4, maxSpans = 80;
  int quadNodes = 5;          // 每 span 的 Gauss-Legendre 节点数
  // ---- 边界条件 ----
  bool pinStartCurvature = true;
  bool pinGoalCurvature = false;
  double goalCurvature = 0.0;  // 1/物理长度
  // ---- 安全裕度（物理长度）----
  double optMargin = 0.0;            // 优化时额外膨胀
  double checkMargin = 0.0;          // 校验时额外膨胀（应小于 optMargin）
  double frontEndExtraMargin = -1;   // 前端额外膨胀，<0 自动 = 0.5*spanLength
  // ---- 代价权重 ----
  double wEnergy = 1.0, wBend = 0.02, wCurv = 1e3, wObs = 1e4, wRadial = 1.0, wProx = 0.05;
  double curvSafety = 0.9;    // 优化中使用 curvSafety*kappaMax 作为曲率上限
  double minSpeedFrac = 0.3;  // 参数速度下限 = minSpeedFrac * L/M（防止尖点/打结）
  double maxStepFrac = 0.5;   // L-BFGS 每步控制点最大位移 = maxStepFrac * spanLength
  double obsScale = -1;       // 穿透量归一化尺度（物理），<0 自动 = 0.2*spanLength
  // ---- 求解器 ----
  LbfgsParams lbfgsFull;      // 从零规划
  LbfgsParams lbfgsReplan;    // 每周期重规划（热启动）
  // ---- 重规划 / 拼接 ----
  double stitchPosTol = -1;               // <0 自动 = 0.5*spanLength
  double stitchHeadingTol = deg2rad(15);  // rad
  double spanRatioLo = 0.7, spanRatioHi = 1.4;  // span 数迟滞
  double freezeDistance = -1;             // 剩余长度小于此值停止重规划，<0 自动 = 1/kappaMax
  // ---- 校验 ----
  double checkStep = -1;      // 校验采样步长（物理），<0 自动 = 0.1*spanLength
  double curvTolerance = 1.05;
  // ---- 同伦切换迟滞 ----
  double switchGain = 0.15;   // 新同伦类须短于当前剩余长度的 (1-gain)
  int switchHoldCount = 3;    // 且连续满足的次数
  // ---- 前端引导段（让初值满足起终点航向）----
  double leadLength = -1;     // <0 自动 = 1/kappaMax

  void finalize() {
    if (frontEndExtraMargin < 0) frontEndExtraMargin = 0.5 * spanLength;
    if (obsScale < 0) obsScale = 0.2 * spanLength;
    if (stitchPosTol < 0) stitchPosTol = 0.5 * spanLength;
    if (freezeDistance < 0) freezeDistance = 1.0 / kappaMax;
    if (checkStep < 0) checkStep = 0.1 * spanLength;
    if (leadLength < 0) leadLength = 1.0 / kappaMax;
    minSpans = std::max(minSpans, 3);
    maxSpans = std::max(maxSpans, minSpans);
  }
  PlannerParams() {
    lbfgsFull.maxIter = 300;
    lbfgsReplan.maxIter = 40;
    lbfgsReplan.fRelTol = 1e-6;
  }
};

enum class PlanStatus {
  Ok,                          // 接受新解，同伦类不变
  ClassChangedForced,          // 旧路径已不可行，接受了不同同伦类的新解
  KeptPreviousHomotopyGuard,   // 新解改变了同伦类而旧解仍可行 → 保留旧解
  KeptPreviousInfeasible,      // 新解校验失败 → 保留旧解
  FrontEndReplan,              // 新旧都不可行 → 前端重新规划
  Frozen,                      // 接近终点，冻结
  Failed
};

inline const char* toString(PlanStatus s) {
  switch (s) {
    case PlanStatus::Ok: return "Ok";
    case PlanStatus::ClassChangedForced: return "ClassChangedForced";
    case PlanStatus::KeptPreviousHomotopyGuard: return "KeptPreviousHomotopyGuard";
    case PlanStatus::KeptPreviousInfeasible: return "KeptPreviousInfeasible";
    case PlanStatus::FrontEndReplan: return "FrontEndReplan";
    case PlanStatus::Frozen: return "Frozen";
    default: return "Failed";
  }
}

struct PlanStats {
  PlanStatus status = PlanStatus::Failed;
  bool valid = false;
  double solveMs = 0.0;
  int iters = 0, evals = 0;
  double cost = 0.0;
  double length = 0.0;          // 物理
  double minClearance = 0.0;    // 物理（相对 radius+checkMargin）
  double maxKappa = 0.0;        // 物理
  double homotopyDiff = 0.0;    // 新旧签名最大差（rad）
  int spans = 0;
  bool stitched = false;
};

class SphereTrajPlanner {
 public:
  struct Validation {
    bool ok = true;
    double minClearance = std::numeric_limits<double>::infinity();
    double maxKappa = 0.0;
    double minCnorm = std::numeric_limits<double>::infinity();
  };

  explicit SphereTrajPlanner(const PlannerParams& prm) : prm_(prm) { prm_.finalize(); }

  void setObstacles(const std::vector<Obstacle>& obs) { obs_ = obs; }
  void addObstacle(const Obstacle& o) { obs_.push_back(o); }
  const std::vector<Obstacle>& obstacles() const { return obs_; }
  void setGoal(const SurfaceState& g) {
    goal_ = normalizeState(g);
    hasGoal_ = true;
  }
  void setHomotopyPole(const Vec3& z) {
    pole_ = z.normalized();
    hasPole_ = true;
  }
  bool hasPlan() const { return curve_.valid(); }
  const SphereSpline& curve() const { return curve_; }
  const PlannerParams& params() const { return prm_; }
  double lastU() const { return lastU_; }

  // ---------------- 从零规划：前端切线图 → 拟合 → 优化 → 校验 ----------------
  PlanStats planFromScratch(const SurfaceState& start0) {
    const auto t0 = std::chrono::steady_clock::now();
    PlanStats st;
    if (!hasGoal_) return st;
    const SurfaceState start = normalizeState(start0);
    const double R = prm_.R;
    std::vector<Vec3> poly;
    if (!frontEnd(start, goal_, poly)) {
      st.solveMs = msSince(t0);
      return st;
    }
    const double LU = polylineLength(poly);
    const int M = std::clamp(static_cast<int>(std::ceil(LU * R / prm_.spanLength)), prm_.minSpans, prm_.maxSpans);
    const std::vector<Vec3> pts = resamplePolyline(poly, 4 * M + 1);
    Problem pb = makeProblem(M, start, goal_, LU);
    const std::vector<Vec3> P0 = fitCtrl(M, pts, start.p, goal_.p);
    std::vector<double> x = toVars(pb.bd, P0, 0.3 * LU / (3.0 * M));
    SolveOut so = solve(pb, x, prm_.lbfgsFull);
    Validation v = validate(so.sp, 0.0);
    for (int retry = 0; retry < 2 && !v.ok; ++retry) {  // 加大障碍/曲率权重再解
      pb.wO *= 10.0;
      pb.wK *= 10.0;
      so = solve(pb, so.x, prm_.lbfgsFull);
      v = validate(so.sp, 0.0);
    }
    if (!hasPole_) setDefaultPole(start.p, goal_.p);
    if (v.ok || !curve_.valid()) {
      curve_ = so.sp;
      lastU_ = 0.0;
    }
    fillStats(st, so, v, M);
    st.status = v.ok ? PlanStatus::Ok : PlanStatus::Failed;
    st.solveMs = msSince(t0);
    return st;
  }

  // ---------------- 每个控制周期调用：拼接 + 热启动 + 近端优化 + 校验 + 同伦守卫 ----------------
  // advanceDist：计算/通信延迟期间车辆预计前进的距离（物理），拼接点取在其前方。
  PlanStats replan(const SurfaceState& meas0, double advanceDist = 0.0) {
    const auto t0 = std::chrono::steady_clock::now();
    if (!curve_.valid()) return planFromScratch(meas0);
    PlanStats st;
    const SurfaceState meas = normalizeState(meas0);
    const double R = prm_.R;
    const SphereSpline& old = curve_;
    const double uV = old.closestU(meas.p, lastU_, 0.5, 3.0);
    lastU_ = uV;
    const double sV = old.sOfU(uV), L = old.length();
    if (L - sV < prm_.freezeDistance) {
      st.status = PlanStatus::Frozen;
      st.valid = true;
      st.length = L;
      st.spans = old.numSpans();
      st.solveMs = msSince(t0);
      return st;
    }
    // 1) 拼接：误差小时从旧轨迹上取起点状态（含曲率），否则用实测状态
    const SphereSpline::Geom gv = old.geom(uV);
    const double posErr = R * angleBetween(meas.p, gv.pos);
    const Vec3 Tp = projectTangent(gv.tangent, meas.p).normalized();
    const double hErr = std::atan2(Tp.cross(meas.h).dot(meas.p), Tp.dot(meas.h));
    SurfaceState start;
    double u0;
    st.stitched = posErr <= prm_.stitchPosTol && std::fabs(hErr) <= prm_.stitchHeadingTol;
    if (st.stitched) {
      u0 = old.uOfS(std::min(sV + advanceDist, L - 0.5 * prm_.freezeDistance));
      const SphereSpline::Geom g0 = old.geom(u0);
      start = {g0.pos, g0.tangent, g0.kappaU / R};
    } else {
      u0 = uV;
      start = meas;
    }
    const double s0 = old.sOfU(u0);
    const double remU = (L - s0) / R;
    // 2) 热启动：旧轨迹剩余部分按参数线性重采样 → 最小二乘拟合新控制点
    //    （旧解由能量项保证近似匀速，按参数采样可保留其参数化，拟合误差极小）
    const int M = chooseSpans(remU * R, old.numSpans());
    std::vector<Vec3> pts = old.sampleUniformU(u0, old.uMax(), 4 * M + 1);
    pts.front() = start.p;
    Problem pb = makeProblem(M, start, goal_, remU);
    const std::vector<Vec3> Pw = fitCtrl(M, pts, start.p, goal_.p);
    std::vector<double> x = toVars(pb.bd, Pw, 0.3 * remU / (3.0 * M));
    buildCtrl(pb.bd, x, pb.anchors);  // 近端锚点
    // 3) 优化
    const SolveOut so = solve(pb, x, prm_.lbfgsReplan);
    // 4) 校验 + 同伦签名比较
    const Validation vNew = validate(so.sp, 0.0);
    st.homotopyDiff = maxAbsDiff(signature(so.sp, 0.0), signature(old, u0));
    const bool sameClass = st.homotopyDiff < kPi;
    fillStats(st, so, vNew, M);
    // 旧轨迹剩余段的校验只在需要做决策时才计算（常规周期省掉这部分开销）
    const bool needOld = !(vNew.ok && sameClass);
    const Validation vOld = needOld ? validate(old, u0) : Validation{};
    if (vNew.ok && (sameClass || !vOld.ok)) {
      curve_ = so.sp;
      lastU_ = 0.0;
      st.status = sameClass ? PlanStatus::Ok : PlanStatus::ClassChangedForced;
    } else if (vNew.ok) {
      st.status = PlanStatus::KeptPreviousHomotopyGuard;
    } else if (vOld.ok) {
      st.status = PlanStatus::KeptPreviousInfeasible;
    } else {
      const PlanStats fs = planFromScratch(start);
      st = fs;
      st.status = fs.valid ? PlanStatus::FrontEndReplan : PlanStatus::Failed;
    }
    st.solveMs = msSince(t0);
    return st;
  }

  // ---------------- 低频调用（如 1~5 Hz）：带迟滞的同伦类主动切换 ----------------
  bool considerClassSwitch(const SurfaceState& meas0) {
    if (!curve_.valid()) return false;
    const SurfaceState meas = normalizeState(meas0);
    const double uV = curve_.closestU(meas.p, lastU_, 0.5, 3.0);
    const double remaining = curve_.length() - curve_.sOfU(uV);
    if (remaining < 2.0 * prm_.freezeDistance) return false;
    std::vector<Vec3> poly;
    if (!frontEnd(meas, goal_, poly)) {
      switchHold_ = 0;
      return false;
    }
    const double Lfe = polylineLength(poly) * prm_.R;
    const double diff = maxAbsDiff(signatureOfPoints(poly), signature(curve_, uV));
    if (diff > kPi && Lfe < (1.0 - prm_.switchGain) * remaining) ++switchHold_;
    else switchHold_ = 0;
    if (switchHold_ < prm_.switchHoldCount) return false;
    switchHold_ = 0;
    const PlanStats st = planFromScratch(meas);  // 仅在有效时替换 curve_
    return st.valid;
  }

  // ---------------- 跟踪参考量（最近点 / 横向误差 / 前视点 / 纯追踪曲率） ----------------
  TrackingRef trackingReference(const SurfaceState& meas, double lookahead) {
    TrackingRef r = computeTrackingRef(curve_, meas.p.normalized(), meas.h, lastU_, lookahead);
    lastU_ = r.u;
    return r;
  }

  // ---------------- 校验：碰撞（障碍半径 + checkMargin）、曲率、退化 ----------------
  Validation validate(const SphereSpline& sp, double u0) const {
    Validation v;
    const double R = prm_.R, stepU = prm_.checkStep / R;
    const int M = sp.numSpans();
    for (int k = std::clamp(static_cast<int>(std::floor(u0)), 0, M - 1); k < M; ++k) {
      const int n = std::max(2, static_cast<int>(std::ceil(sp.spanLengthU(k) / stepU)));
      for (int i = 0; i <= n; ++i) {
        const double u = k + double(i) / n;
        if (u < u0) continue;
        const SphereSpline::Geom g = sp.geom(u);
        v.maxKappa = std::max(v.maxKappa, std::fabs(g.kappaU) / R);
        v.minCnorm = std::min(v.minCnorm, g.cnorm);
        for (const Obstacle& o : obs_) {
          const double clr = R * angleBetween(g.pos, o.center) - o.radius - prm_.checkMargin;
          v.minClearance = std::min(v.minClearance, clr);
        }
      }
    }
    v.ok = v.minClearance >= 0.0 && v.maxKappa <= prm_.kappaMax * prm_.curvTolerance && v.minCnorm > 0.2;
    return v;
  }

  // ---------------- 同伦签名（球面版 H-signature） ----------------
  // 以 pole 为极点做球极投影，累计路径绕每个障碍像点的转角。
  std::vector<double> signatureOfPoints(const std::vector<Vec3>& pts) const {
    const StereoChart chart(pole_);
    std::vector<std::pair<double, double>> oc;
    for (const Obstacle& o : obs_) {
      double a = 0, b = 0;
      if (angleBetween(o.center, pole_) * prm_.R < o.radius || !chart.project(o.center, a, b)) {
        oc.push_back({std::numeric_limits<double>::quiet_NaN(), 0.0});  // 含极点的障碍：作为外边界，不计
      } else {
        oc.push_back({a, b});
      }
    }
    std::vector<double> w(obs_.size(), 0.0);
    bool havePrev = false;
    double pu = 0, pv = 0;
    for (const Vec3& x : pts) {
      double u, v;
      if (!chart.project(x, u, v)) continue;
      if (havePrev) {
        for (size_t i = 0; i < oc.size(); ++i) {
          if (std::isnan(oc[i].first)) continue;
          const double ax = pu - oc[i].first, ay = pv - oc[i].second;
          const double bx = u - oc[i].first, by = v - oc[i].second;
          w[i] += std::atan2(ax * by - ay * bx, ax * bx + ay * by);
        }
      }
      pu = u;
      pv = v;
      havePrev = true;
    }
    return w;
  }
  std::vector<double> signature(const SphereSpline& sp, double u0) const {
    const double stepU = std::min(prm_.checkStep / prm_.R, deg2rad(0.5));
    return signatureOfPoints(sp.samplePositions(u0, stepU));
  }

  // ---------------- 诊断：解析梯度 vs 中心差分（返回最大误差 / 最大梯度分量） ----------------
  double gradientCheck(const SurfaceState& start0, unsigned seed = 1) {
    const SurfaceState start = normalizeState(start0);
    std::vector<Vec3> poly;
    if (!frontEnd(start, goal_, poly)) return -1.0;
    const double LU = polylineLength(poly);
    const int M = std::clamp(static_cast<int>(std::ceil(LU * prm_.R / prm_.spanLength)), prm_.minSpans, prm_.maxSpans);
    Problem pb = makeProblem(M, start, goal_, LU);
    for (CapOpt& c : pb.caps) {  // 放大障碍、收紧曲率，让所有罚项都处于激活状态
      c.r *= 1.6;
      c.cosr = std::cos(c.r);
      c.sinr = std::sin(c.r);
    }
    pb.kapLimU *= 0.3;
    std::vector<Vec3> P0 = fitCtrl(M, resamplePolyline(poly, 4 * M + 1), start.p, goal_.p);
    std::vector<double> x = toVars(pb.bd, P0, 0.3 * LU / (3.0 * M));
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, 0.01);
    buildCtrl(pb.bd, x, pb.anchors);
    for (double& xi : x) xi += nd(rng);
    std::vector<double> g, gtmp;
    costGrad(pb, x, g);
    double maxErr = 0.0, maxG = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
      const double h = 1e-6 * std::max(1.0, std::fabs(x[i]));
      std::vector<double> xp = x, xm = x;
      xp[i] += h;
      xm[i] -= h;
      const double fd = (costGrad(pb, xp, gtmp) - costGrad(pb, xm, gtmp)) / (2.0 * h);
      maxErr = std::max(maxErr, std::fabs(fd - g[i]));
      maxG = std::max(maxG, std::fabs(g[i]));
    }
    return maxErr / std::max(maxG, 1e-12);
  }

 private:
  // ======================= 内部数据结构 =======================
  struct Boundary {
    int M = 0;
    Vec3 ps, hs, ls, pg, hg, lg;
    double kapS = 0.0, kapG = 0.0;  // 单位球曲率
    bool pinS = true, pinG = false;
    int numVars() const { return 2 + (pinS ? 2 : 3) + 3 * (M - 3) + (pinG ? 2 : 3); }
  };
  struct CapOpt {
    Vec3 n;
    double r = 0, cosr = 1, sinr = 0;
  };
  // 采样表：每个 span 含 nq 个 Gauss-Legendre 节点（用于积分项与罚项）
  // + 节点 u=k（以及最后一个 span 的 u=M）处的附加点（仅用于曲率/障碍罚项，防止端点处曲率尖峰漏检）
  struct SampleTable {
    int M = 0, nq = 0;
    std::vector<int> offset;                // span k 的样本为 [offset[k], offset[k+1])
    std::vector<double> wInt, wPen;         // 积分项权重 / 罚项权重（u 测度）
    std::vector<std::array<double, 12>> N;  // N[d*4+j]
  };
  struct Problem {
    Boundary bd;
    const SampleTable* tab = nullptr;
    std::vector<CapOpt> caps;
    double LrefU = 1, kapLimU = 1, obsScaleU = 1, proxScaleU = 1, sMinU = 0;
    double wE = 1, wB = 0, wK = 0, wO = 0, wR = 0, wP = 0;
    std::vector<Vec3> anchors;  // 为空则不加近端项
  };
  struct SolveOut {
    std::vector<double> x;
    LbfgsResult res;
    SphereSpline sp;
  };
  struct FitCache {
    int M = 0, K = 0, n = 0;
    std::vector<int> first;
    std::vector<std::array<double, 4>> N;
    std::vector<double> L;  // Cholesky 下三角 n×n
  };

  // ======================= 变量 <-> 控制点 =======================
  static void buildCtrl(const Boundary& b, const std::vector<double>& x, std::vector<Vec3>& P) {
    const int M = b.M;
    P.resize(M + 3);
    const double as = std::exp(x[0]), ag = std::exp(x[1]);
    int k = 2;
    P[0] = b.ps;
    P[1] = b.ps + b.hs * as;
    if (b.pinS) {
      P[2] = b.ps * x[k] + b.hs * x[k + 1] + b.ls * (3.0 * as * as * b.kapS);
      k += 2;
    } else {
      P[2] = Vec3(x[k], x[k + 1], x[k + 2]);
      k += 3;
    }
    for (int j = 3; j <= M - 1; ++j, k += 3) P[j] = Vec3(x[k], x[k + 1], x[k + 2]);
    if (b.pinG) {
      P[M] = b.pg * x[k] + b.hg * x[k + 1] + b.lg * (3.0 * ag * ag * b.kapG);
    } else {
      P[M] = Vec3(x[k], x[k + 1], x[k + 2]);
    }
    P[M + 1] = b.pg - b.hg * ag;
    P[M + 2] = b.pg;
  }

  static void chainGrad(const Boundary& b, const std::vector<double>& x, const std::vector<Vec3>& gP,
                        std::vector<double>& g) {
    const int M = b.M;
    g.assign(x.size(), 0.0);
    const double as = std::exp(x[0]), ag = std::exp(x[1]);
    g[0] = gP[1].dot(b.hs) * as;
    g[1] = -gP[M + 1].dot(b.hg) * ag;
    int k = 2;
    auto put3 = [&](const Vec3& v) {
      g[k] = v.x;
      g[k + 1] = v.y;
      g[k + 2] = v.z;
      k += 3;
    };
    if (b.pinS) {
      g[0] += gP[2].dot(b.ls) * 6.0 * as * as * b.kapS;
      g[k] = gP[2].dot(b.ps);
      g[k + 1] = gP[2].dot(b.hs);
      k += 2;
    } else {
      put3(gP[2]);
    }
    for (int j = 3; j <= M - 1; ++j) put3(gP[j]);
    if (b.pinG) {
      g[1] += gP[M].dot(b.lg) * 6.0 * ag * ag * b.kapG;
      g[k] = gP[M].dot(b.pg);
      g[k + 1] = gP[M].dot(b.hg);
    } else {
      put3(gP[M]);
    }
  }

  static std::vector<double> toVars(const Boundary& b, const std::vector<Vec3>& P, double aMin) {
    const int M = b.M;
    std::vector<double> x(b.numVars());
    x[0] = std::log(std::max((P[1] - b.ps).dot(b.hs), aMin));
    x[1] = std::log(std::max((b.pg - P[M + 1]).dot(b.hg), aMin));
    int k = 2;
    auto put3 = [&](const Vec3& v) {
      x[k] = v.x;
      x[k + 1] = v.y;
      x[k + 2] = v.z;
      k += 3;
    };
    if (b.pinS) {
      x[k] = P[2].dot(b.ps);
      x[k + 1] = P[2].dot(b.hs);
      k += 2;
    } else {
      put3(P[2]);
    }
    for (int j = 3; j <= M - 1; ++j) put3(P[j]);
    if (b.pinG) {
      x[k] = P[M].dot(b.pg);
      x[k + 1] = P[M].dot(b.hg);
    } else {
      put3(P[M]);
    }
    return x;
  }

  // ======================= 代价与梯度 =======================
  // 单个采样点的形状项（能量、弯曲、曲率上限）及其对 (c, c', c'') 的梯度（手推解析式）
  static double shapeTerm(const Vec3& c, const Vec3& c1, const Vec3& c2, double aE, double aB, double aK,
                          double kapLim, double sMin, Vec3& gc, Vec3& gc1, Vec3& gc2) {
    const double r2 = c.squaredNorm();
    const Vec3 wv = c.cross(c1);
    const double W2 = wv.squaredNorm();
    double J = aE * W2 / (r2 * r2);       // s² = W2 / r2²
    double dW2 = aE / (r2 * r2);
    double dr2 = -2.0 * aE * W2 / (r2 * r2 * r2);
    double dd = 0.0;
    if (W2 > 1e-24) {
      const double r = std::sqrt(r2), W = std::sqrt(W2);
      const double d = wv.dot(c2);
      const double s = W / r2;
      const double kap = r2 * r * d / (W2 * W);
      const double dk_dr2 = 1.5 * kap / r2, dk_dW2 = -1.5 * kap / W2, dk_dd = r2 * r / (W2 * W);
      const double ds_dW2 = 0.5 / (W * r2), ds_dr2 = -W / (r2 * r2);
      J += aB * kap * kap * s;
      double gk = aB * 2.0 * kap * s, gs = aB * kap * kap;
      const double v = (std::fabs(kap) - kapLim) / kapLim;
      if (v > 0) {
        J += aK * v * v * v;
        gk += aK * 3.0 * v * v * (kap > 0 ? 1.0 : -1.0) / kapLim;
      }
      const double vs = (sMin - s) / sMin;  // 速度下限
      if (vs > 0) {
        J += aK * vs * vs * vs;
        gs -= aK * 3.0 * vs * vs / sMin;
      }
      dW2 += gk * dk_dW2 + gs * ds_dW2;
      dr2 += gk * dk_dr2 + gs * ds_dr2;
      dd += gk * dk_dd;
    }
    // ∂r2/∂c = 2c；∂W2/∂c = 2 c'×w；∂W2/∂c' = 2 w×c；∂d/∂c = c'×c''；∂d/∂c' = c''×c；∂d/∂c'' = w
    gc = c * (2.0 * dr2) + c1.cross(wv) * (2.0 * dW2) + c1.cross(c2) * dd;
    gc1 = wv.cross(c) * (2.0 * dW2) + c2.cross(c) * dd;
    gc2 = wv * dd;
    return J;
  }

  double costGrad(const Problem& pb, const std::vector<double>& x, std::vector<double>& g) {
    const int M = pb.bd.M;
    std::vector<Vec3>& P = bufP_;
    std::vector<Vec3>& gP = bufGP_;
    buildCtrl(pb.bd, x, P);
    gP.assign(M + 3, Vec3());
    double J = 0.0;
    // 径向正则
    for (int j = 2; j <= M; ++j) {
      const double e = P[j].squaredNorm() - 1.0;
      J += pb.wR * e * e;
      gP[j] += P[j] * (4.0 * pb.wR * e);
    }
    // 近端项
    if (!pb.anchors.empty() && pb.wP > 0) {
      const double s = pb.wP / (pb.proxScaleU * pb.proxScaleU);
      for (int j = 1; j <= M + 1; ++j) {
        const Vec3 d = P[j] - pb.anchors[j];
        J += s * d.squaredNorm();
        gP[j] += d * (2.0 * s);
      }
    }
    // 采样点项
    const SampleTable& tab = *pb.tab;
    const double aE0 = pb.wE * M / (pb.LrefU * pb.LrefU);
    const double aB0 = pb.wB / (pb.kapLimU * pb.kapLimU * pb.LrefU);
    const double aK0 = pb.wK / M, aO0 = pb.wO / M;
    std::vector<int>& nearIdx = bufNear_;
    for (int k = 0; k < M; ++k) {
      // 凸包性质：该 span 落在其 4 个控制点的包围球冠内 → 精确剔除远处障碍
      nearIdx.clear();
      if (!pb.caps.empty()) {
        Vec3 ph[4], m;
        for (int j = 0; j < 4; ++j) {
          ph[j] = P[k + j].normalized();
          m += ph[j];
        }
        m = m.normalized();
        double cphi = 1.0;
        for (int j = 0; j < 4; ++j) cphi = std::min(cphi, m.dot(ph[j]));
        const double phi = std::acos(clampd(cphi, -1.0, 1.0));
        for (int i = 0; i < static_cast<int>(pb.caps.size()); ++i)
          if (angleBetween(m, pb.caps[i].n) < phi + pb.caps[i].r + 1e-6) nearIdx.push_back(i);
      }
      for (int q = tab.offset[k]; q < tab.offset[k + 1]; ++q) {
        const std::array<double, 12>& N = tab.N[q];
        const double w = tab.wInt[q], wp = tab.wPen[q];
        Vec3 c, c1, c2;
        for (int j = 0; j < 4; ++j) {
          c += P[k + j] * N[j];
          c1 += P[k + j] * N[4 + j];
          c2 += P[k + j] * N[8 + j];
        }
        Vec3 gc, gc1, gc2;
        J += shapeTerm(c, c1, c2, aE0 * w, aB0 * w, aK0 * wp, pb.kapLimU, pb.sMinU, gc, gc1, gc2);
        if (!nearIdx.empty()) {
          const double r = c.norm();
          const Vec3 gam = c / r;
          for (int i : nearIdx) {
            const CapOpt& cp = pb.caps[i];
            const double scale = 1.0 / (std::max(cp.sinr, 1e-3) * pb.obsScaleU);
            const double nd = gam.dot(cp.n);
            const double v = (nd - cp.cosr) * scale;
            if (v <= 0) continue;
            J += aO0 * wp * v * v * v;
            gc += (cp.n - gam * nd) * (3.0 * aO0 * wp * v * v * scale / r);  // ∂(γ·n)/∂c = (n-(γ·n)γ)/|c|
          }
        }
        for (int j = 0; j < 4; ++j) gP[k + j] += gc * N[j] + gc1 * N[4 + j] + gc2 * N[8 + j];
      }
    }
    chainGrad(pb.bd, x, gP, g);
    return J;
  }

  SolveOut solve(Problem& pb, std::vector<double> x, const LbfgsParams& lp) {
    auto fg = [&](const std::vector<double>& xx, std::vector<double>& gg) { return costGrad(pb, xx, gg); };
    SolveOut out;
    LbfgsParams lp2 = lp;
    lp2.maxStep = std::min(lp.maxStep, prm_.maxStepFrac * prm_.spanLength / prm_.R);
    out.res = lbfgsMinimize(fg, x, lp2);
    std::vector<Vec3> P;
    buildCtrl(pb.bd, x, P);
    out.sp.set(pb.bd.M, P, prm_.R);
    out.x = std::move(x);
    return out;
  }

  // ======================= 问题构造 =======================
  const SampleTable& table(int M, int nq) {
    auto it = tables_.find({M, nq});
    if (it != tables_.end()) return it->second;
    SampleTable t;
    t.M = M;
    t.nq = nq;
    std::vector<double> gx, gw;
    gaussLegendre01(nq, gx, gw);
    ClampedCubicBasis B(M);
    auto push = [&](int k, double u, double wInt, double wPen) {
      double D[3][4];
      B.eval(u, k, D);
      std::array<double, 12> a;
      for (int d = 0; d < 3; ++d)
        for (int j = 0; j < 4; ++j) a[d * 4 + j] = D[d][j];
      t.N.push_back(a);
      t.wInt.push_back(wInt);
      t.wPen.push_back(wPen);
    };
    for (int k = 0; k < M; ++k) {
      t.offset.push_back(static_cast<int>(t.N.size()));
      push(k, double(k), 0.0, 1.0 / nq);
      for (int q = 0; q < nq; ++q) push(k, k + gx[q], gw[q], gw[q]);
      if (k == M - 1) push(k, double(M), 0.0, 1.0 / nq);
    }
    t.offset.push_back(static_cast<int>(t.N.size()));
    return tables_.emplace(std::make_pair(M, nq), std::move(t)).first->second;
  }

  Problem makeProblem(int M, const SurfaceState& s, const SurfaceState& g, double LrefU) {
    const double R = prm_.R, kmaxU = prm_.kappaMax * R;
    Problem pb;
    Boundary& b = pb.bd;
    b.M = M;
    b.ps = s.p;
    b.hs = s.h;
    b.ls = s.p.cross(s.h);
    b.kapS = clampd(s.kappa * R, -kmaxU, kmaxU);
    b.pg = g.p;
    b.hg = g.h;
    b.lg = g.p.cross(g.h);
    b.kapG = clampd(prm_.goalCurvature * R, -kmaxU, kmaxU);
    b.pinS = prm_.pinStartCurvature;
    b.pinG = prm_.pinGoalCurvature;
    pb.tab = &table(M, prm_.quadNodes);
    for (const Obstacle& o : obs_) {
      CapOpt c;
      c.n = o.center.normalized();
      c.r = (o.radius + prm_.optMargin) / R;
      c.cosr = std::cos(c.r);
      c.sinr = std::sin(c.r);
      pb.caps.push_back(c);
    }
    pb.LrefU = std::max(LrefU, 1e-9);
    pb.kapLimU = prm_.curvSafety * kmaxU;
    pb.obsScaleU = prm_.obsScale / R;
    pb.proxScaleU = prm_.spanLength / R;
    pb.sMinU = prm_.minSpeedFrac * pb.LrefU / M;
    pb.wE = prm_.wEnergy;
    pb.wB = prm_.wBend;
    pb.wK = prm_.wCurv;
    pb.wO = prm_.wObs;
    pb.wR = prm_.wRadial;
    pb.wP = prm_.wProx;
    return pb;
  }

  // 等参数最小二乘拟合：样本 u_k = M k/(K-1)，P0 与 P_{M+2} 固定；Cholesky 因子按 (M,K) 缓存
  std::vector<Vec3> fitCtrl(int M, const std::vector<Vec3>& pts, const Vec3& p0, const Vec3& pEnd) {
    const int K = static_cast<int>(pts.size());
    FitCache& fc = fitCache(M, K);
    const int n = fc.n;
    std::vector<Vec3> rhs(n);
    for (int k = 0; k < K; ++k) {
      const int f = fc.first[k];
      const auto& N = fc.N[k];
      Vec3 y = pts[k];
      for (int j = 0; j < 4; ++j) {
        const int idx = f + j;
        if (idx == 0) y -= p0 * N[j];
        if (idx == M + 2) y -= pEnd * N[j];
      }
      for (int j = 0; j < 4; ++j) {
        const int idx = f + j;
        if (idx >= 1 && idx <= M + 1) rhs[idx - 1] += y * N[j];
      }
    }
    // 前代 / 回代
    const std::vector<double>& L = fc.L;
    for (int i = 0; i < n; ++i) {
      Vec3 s = rhs[i];
      for (int j = 0; j < i; ++j) s -= rhs[j] * L[i * n + j];
      rhs[i] = s / L[i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
      Vec3 s = rhs[i];
      for (int j = i + 1; j < n; ++j) s -= rhs[j] * L[j * n + i];
      rhs[i] = s / L[i * n + i];
    }
    std::vector<Vec3> P(M + 3);
    P[0] = p0;
    P[M + 2] = pEnd;
    for (int i = 0; i < n; ++i) P[i + 1] = rhs[i];
    return P;
  }

  FitCache& fitCache(int M, int K) {
    auto it = fitCaches_.find({M, K});
    if (it != fitCaches_.end()) return it->second;
    FitCache fc;
    fc.M = M;
    fc.K = K;
    fc.n = M + 1;
    const int n = fc.n;
    ClampedCubicBasis B(M);
    std::vector<double> A(n * n, 0.0);
    for (int k = 0; k < K; ++k) {
      double u = M * double(k) / (K - 1);
      const int f = B.spanOf(u);
      double D[3][4];
      B.eval(u, f, D);
      fc.first.push_back(f);
      fc.N.push_back({D[0][0], D[0][1], D[0][2], D[0][3]});
      for (int a = 0; a < 4; ++a) {
        const int ia = f + a - 1;
        if (ia < 0 || ia >= n) continue;
        for (int b = 0; b < 4; ++b) {
          const int ib = f + b - 1;
          if (ib < 0 || ib >= n) continue;
          A[ia * n + ib] += D[0][a] * D[0][b];
        }
      }
    }
    for (int i = 0; i < n; ++i) A[i * n + i] += 1e-10;
    fc.L.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = A[i * n + j];
        for (int k = 0; k < j; ++k) s -= fc.L[i * n + k] * fc.L[j * n + k];
        fc.L[i * n + j] = (i == j) ? std::sqrt(std::max(s, 1e-300)) : s / fc.L[j * n + j];
      }
    }
    return fitCaches_.emplace(std::make_pair(M, K), std::move(fc)).first->second;
  }

  // ======================= 前端 =======================
  bool frontEnd(const SurfaceState& s, const SurfaceState& g, std::vector<Vec3>& poly) const {
    const double R = prm_.R;
    std::vector<Cap> caps;
    for (const Obstacle& o : obs_) {
      Cap c{o.center.normalized(), (o.radius + prm_.optMargin + prm_.frontEndExtraMargin) / R};
      // 起/终点落在（膨胀后的）障碍内时，仅对前端收缩该障碍
      for (const Vec3& e : {s.p, g.p}) c.radius = std::min(c.radius, angleBetween(e, c.center) - 1e-4);
      if (c.radius > 0) caps.push_back(c);
    }
    auto leadPoint = [&](const Vec3& p, const Vec3& dir) {
      for (double d = prm_.leadLength / R; d > 1e-9; d *= 0.5) {
        const Vec3 q = geodesicStep(p, dir, d);
        bool ok = true;
        for (const Cap& c : caps)
          if (detail::pointInCap(q, c) || detail::arcHitsCap(p, q, c)) { ok = false; break; }
        if (ok) return q;
      }
      return p;
    };
    const Vec3 ps2 = leadPoint(s.p, s.h), pg2 = leadPoint(g.p, -g.h);
    TangentGraph tg;
    const FrontEndPath fe = tg.plan(ps2, pg2, caps);
    if (!fe.ok) return false;
    poly.clear();
    if (angleBetween(ps2, s.p) > 1e-12) poly.push_back(s.p);
    poly.insert(poly.end(), fe.points.begin(), fe.points.end());
    if (angleBetween(pg2, g.p) > 1e-12) poly.push_back(g.p);
    return true;
  }

  static double polylineLength(const std::vector<Vec3>& poly) {
    double L = 0;
    for (size_t i = 1; i < poly.size(); ++i) L += angleBetween(poly[i - 1], poly[i]);
    return L;
  }

  static std::vector<Vec3> resamplePolyline(const std::vector<Vec3>& poly, int K) {
    std::vector<double> cum(poly.size(), 0.0);
    for (size_t i = 1; i < poly.size(); ++i) cum[i] = cum[i - 1] + angleBetween(poly[i - 1], poly[i]);
    std::vector<Vec3> out(K);
    size_t seg = 1;
    for (int k = 0; k < K; ++k) {
      const double s = cum.back() * k / double(K - 1);
      while (seg + 1 < poly.size() && cum[seg] < s) ++seg;
      const double len = cum[seg] - cum[seg - 1];
      const double t = len > 0 ? clampd((s - cum[seg - 1]) / len, 0.0, 1.0) : 0.0;
      out[k] = slerp(poly[seg - 1], poly[seg], t);
    }
    return out;
  }

  // ======================= 杂项 =======================
  int chooseSpans(double remLen, int Mcur) const {
    const double ratio = remLen / (Mcur * prm_.spanLength);
    int M = Mcur;
    if (ratio < prm_.spanRatioLo || ratio > prm_.spanRatioHi) M = static_cast<int>(std::lround(remLen / prm_.spanLength));
    return std::clamp(M, prm_.minSpans, prm_.maxSpans);
  }

  void setDefaultPole(const Vec3& a, const Vec3& b) {
    Vec3 m = a + b;
    if (m.norm() < 1e-6) {
      Vec3 e1, e2;
      tangentFrame(a, e1, e2);
      m = e1;
    }
    pole_ = -m.normalized();  // 作业区域中点的对跖点
    hasPole_ = true;
  }

  static SurfaceState normalizeState(const SurfaceState& s) {
    SurfaceState o = s;
    o.p = s.p.normalized();
    o.h = projectTangent(s.h, o.p).normalized();
    return o;
  }
  static double maxAbsDiff(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (size_t i = 0; i < std::min(a.size(), b.size()); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
  }
  static double msSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  }
  void fillStats(PlanStats& st, const SolveOut& so, const Validation& v, int M) const {
    st.valid = v.ok;
    st.iters = so.res.iters;
    st.evals = so.res.evals;
    st.cost = so.res.f;
    st.length = so.sp.length();
    st.minClearance = v.minClearance;
    st.maxKappa = v.maxKappa;
    st.spans = M;
  }

  PlannerParams prm_;
  std::vector<Obstacle> obs_;
  SurfaceState goal_;
  bool hasGoal_ = false;
  Vec3 pole_;
  bool hasPole_ = false;
  SphereSpline curve_;
  double lastU_ = 0.0;
  int switchHold_ = 0;
  std::map<std::pair<int, int>, SampleTable> tables_;
  std::map<std::pair<int, int>, FitCache> fitCaches_;
  std::vector<Vec3> bufP_, bufGP_;
  std::vector<int> bufNear_;
};

}  // namespace sphtraj
