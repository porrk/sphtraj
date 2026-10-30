#pragma once
// 简化版：球面中心投影 B 样条轨迹优化 + 切线图前端 + 每周期重规划。
// 目标只有两个 —— 规避禁飞区、到达终点区域。无终端航向约束、无同伦守卫、
// 无近端冻结、无曲率硬上限；曲率以软罚正则（幅值 + 变化率）。
//
// 曲线：c(u) = Σ N_j(u) P_j（R^3 夹持三次 B 样条），γ = c/|c|。
// 边界条件（构造式精确满足）：
//   P0 = ps, P1 = ps + a_s hs                  → 起点航向精确（转向连续）
//   ℓs·P2 = 3 a_s² κ0  (ℓs = ps × hs)           → 起点曲率 = 当前曲率（可选）
//   P_{M+2} = pg；P_M、P_{M+1} 自由              → 终端航向由优化决定
// 代价：J = wE·E + wB·B + wKv·Kv + wSp·Sp + wO·O + wR·Rad
//   E   = M ∫|γ'|² du / Lref²         （≥ (L/Lref)²，等号当且仅当匀速 → 最短且近似弧长参数化）
//   B   = ∫ κ² ds / (κref² Lref)      （曲率幅值软罚：κref = kappaMax·R 仅作归一化尺度）
//   Kv  = (1/M) Σ_q (Δκ_q/κref)²·(spanU/Δs_q)   （曲率变化率软罚：相邻采样点曲率差）
//   Sp  = (1/M) ∫ max(0, 1 - |γ'|/smin)³ du     （参数速度下限，防尖点/打结）
//   O   = (1/M) ∫ Σ_i max(0, (γ·n_i - cos r_i)/(sin r_i σ))³ du
//   Rad = Σ (|P_j|² - 1)²             （消除径向冗余自由度）
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
  double kappaMax = 1.0;      // 参考测地曲率（1/物理长度）。不是约束：仅作曲率罚项的
                              // 归一化尺度与前端引导段长度 (leadLength 自动值) 的基准。
  double spanLength = 0.1;    // 每个 knot span 期望弧长（物理）
  int minSpans = 4, maxSpans = 80;
  int quadNodes = 5;          // 每 span 的 Gauss-Legendre 节点数
  // ---- 边界条件 ----
  bool pinStartCurvature = true;  // 起点曲率 = 实测曲率（replan 转向连续）
  // ---- 安全裕度（物理长度）----
  double optMargin = 0.0;            // 优化时额外膨胀
  double checkMargin = 0.0;          // 校验时额外膨胀（应小于 optMargin）
  double frontEndExtraMargin = -1;   // 前端额外膨胀，<0 自动 = 0.5*spanLength
  // ---- 代价权重 ----
  double wEnergy = 1.0;        // 能量（最短 + 匀速参数化）
  double wBend = 0.02;         // 曲率幅值软罚 ∫κ² ds
  double wKappaVar = 1e3;      // 曲率变化率软罚 Σ(Δκ)²/Δs（相邻采样点曲率差）
  double wSpeed = 1e3;         // 参数速度下限罚权重（防尖点/打结）
  double wObs = 1e4;           // 障碍穿透罚
  double wRadial = 1.0;        // 径向正则
  double minSpeedFrac = 0.3;   // 参数速度下限 = minSpeedFrac * L/M（防止尖点/打结）
  double maxStepFrac = 0.5;    // L-BFGS 每步控制点最大位移 = maxStepFrac * spanLength
  double obsScale = -1;        // 穿透量归一化尺度（物理），<0 自动 = 0.2*spanLength
  // ---- 求解器 ----
  LbfgsParams lbfgsFull;      // 从零规划
  LbfgsParams lbfgsReplan;    // 每周期重规划（热启动）
  // ---- 重规划 / 拼接 ----
  double stitchPosTol = -1;               // <0 自动 = 0.5*spanLength
  double stitchHeadingTol = deg2rad(15);  // rad
  double spanRatioLo = 0.7, spanRatioHi = 1.4;  // span 数迟滞
  // ---- 校验 ----
  double checkStep = -1;      // 校验采样步长（物理），<0 自动 = 0.1*spanLength
  // ---- 前端引导段（让初值满足起点航向）----
  double leadLength = -1;     // <0 自动 = 1/kappaMax
  // ---- 跟踪参考（前视点选取）----
  double lookaheadRatio = 0.0;  // >0 时按剩余弧长比例选前视点：Ld = min(ratio*remaining, lookahead)

  void finalize() {
    if (frontEndExtraMargin < 0) frontEndExtraMargin = 0.5 * spanLength;
    if (obsScale < 0) obsScale = 0.2 * spanLength;
    if (stitchPosTol < 0) stitchPosTol = 0.5 * spanLength;
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
  Ok,                      // 接受新解
  KeptPreviousInfeasible,  // 新解校验失败 → 保留旧解
  FrontEndReplan,          // 新旧都不可行 → 前端重新规划
  Failed
};

inline const char* toString(PlanStatus s) {
  switch (s) {
    case PlanStatus::Ok: return "Ok";
    case PlanStatus::KeptPreviousInfeasible: return "KeptPreviousInfeasible";
    case PlanStatus::FrontEndReplan: return "FrontEndReplan";
    default: return "Failed";
  }
}

struct PlanStats {
  PlanStatus status = PlanStatus::Failed;
  // valid / length / minClearance / maxKappa / spans 描述规划器里保留的那条曲线。
  // 新候选被拒绝时，这些字段仍属于继续跟踪的旧曲线。
  // iters / evals / cost 描述本周期尝试的那次求解。
  bool valid = false;
  double solveMs = 0.0;
  int iters = 0, evals = 0;
  double cost = 0.0;
  double length = 0.0;          // 保留曲线的物理总长；没有保留曲线时为 0
  double minClearance = 0.0;    // 保留曲线受检段相对 radius+checkMargin 的间隔（物理）
  double maxKappa = 0.0;        // 保留曲线受检段的最大物理曲率（统计上报，不设上限）
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
  // 有一条曾通过校验并被采纳的曲线。从零规划失败不会留下落选的初值。
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
    for (int retry = 0; retry < 2 && !v.ok; ++retry) {  // 加大障碍权重再解
      pb.wO *= 10.0;
      so = solve(pb, so.x, prm_.lbfgsFull);
      v = validate(so.sp, 0.0);
    }
    if (v.ok) {
      curve_ = so.sp;
      lastU_ = 0.0;
      fillStats(st, so, v, M);
    } else {
      // 落选的曲线不写入 curve_。已有旧曲线时，几何量改记旧曲线，valid 仍表示本次没有采纳新解。
      st.iters = so.res.iters;
      st.evals = so.res.evals;
      st.cost = so.res.f;
      st.valid = false;
      if (curve_.valid()) {
        const Validation vKeep = validate(curve_, 0.0);
        st.length = curve_.length();
        st.minClearance = vKeep.minClearance;
        st.maxKappa = vKeep.maxKappa;
        st.spans = curve_.numSpans();
      }
    }
    st.status = v.ok ? PlanStatus::Ok : PlanStatus::Failed;
    st.solveMs = msSince(t0);
    return st;
  }

  // ---------------- 每个控制周期调用：拼接 + 热启动 + 优化 + 校验 ----------------
  // advanceDist：计算/通信延迟期间车辆预计前进的距离（物理），拼接点取在其前方。
  // 到达终点前持续重规划（无近端冻结）；新解可行即接受（无同伦守卫）。
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
    // 1) 拼接：误差小时从旧轨迹上取起点状态（含曲率），否则用实测状态
    const SphereSpline::Geom gv = old.geom(uV);
    const double posErr = R * angleBetween(meas.p, gv.pos);
    const Vec3 Tp = projectTangent(gv.tangent, meas.p).normalized();
    const double hErr = std::atan2(Tp.cross(meas.h).dot(meas.p), Tp.dot(meas.h));
    SurfaceState start;
    double u0;
    st.stitched = posErr <= prm_.stitchPosTol && std::fabs(hErr) <= prm_.stitchHeadingTol;
    if (st.stitched) {
      u0 = old.uOfS(std::min(sV + advanceDist, L));
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
    // 3) 优化
    const SolveOut so = solve(pb, x, prm_.lbfgsReplan);
    // 4) 校验：可行即接受
    const Validation vNew = validate(so.sp, 0.0);
    fillStats(st, so, vNew, M);
    if (vNew.ok) {
      curve_ = so.sp;
      lastU_ = 0.0;
      st.status = PlanStatus::Ok;
    } else {
      const Validation vOld = validate(old, u0);
      if (vOld.ok) {
        st.status = PlanStatus::KeptPreviousInfeasible;
        describeRetained(st, old, vOld);
      } else {
        const PlanStats fs = planFromScratch(start);
        const bool accepted = fs.valid;
        st = fs;
        st.status = accepted ? PlanStatus::FrontEndReplan : PlanStatus::Failed;
        if (!accepted) describeRetained(st, curve_, curve_.valid() ? validate(curve_, u0) : Validation{});
      }
    }
    st.solveMs = msSince(t0);
    return st;
  }

  // ---------------- 跟踪参考量（最近点 / 横向误差 / 前视点 / 纯追踪曲率） ----------------
  TrackingRef trackingReference(const SurfaceState& meas, double lookahead) {
    if (!curve_.valid()) return {};
    TrackingRef r = computeTrackingRef(curve_, meas.p.normalized(), meas.h, lastU_, lookahead,
                                        prm_.lookaheadRatio);
    lastU_ = r.u;
    return r;
  }

  // ---------------- 校验：碰撞（障碍半径 + checkMargin）、退化 ----------------
  Validation validate(const SphereSpline& sp, double u0) const {
    Validation v;
    if (!sp.valid()) {
      v.ok = false;
      v.minClearance = -std::numeric_limits<double>::infinity();
      v.minCnorm = 0.0;
      return v;
    }
    const double R = prm_.R, stepU = prm_.checkStep / R;
    const int M = sp.numSpans();
    int k0 = static_cast<int>(std::floor(u0));
    if (k0 < 0) k0 = 0;
    if (k0 >= M) k0 = M - 1;
    for (int k = k0; k < M; ++k) {
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
    v.ok = v.minClearance >= 0.0 && v.minCnorm > 0.2;  // maxKappa 仅统计上报（无曲率硬上限）
    return v;
  }

  // ---------------- 诊断：解析梯度 vs 中心差分（返回最大误差 / 最大梯度分量） ----------------
  double gradientCheck(const SurfaceState& start0, unsigned seed = 1) {
    const SurfaceState start = normalizeState(start0);
    std::vector<Vec3> poly;
    if (!frontEnd(start, goal_, poly)) return -1.0;
    const double LU = polylineLength(poly);
    const int M = std::clamp(static_cast<int>(std::ceil(LU * prm_.R / prm_.spanLength)), prm_.minSpans, prm_.maxSpans);
    Problem pb = makeProblem(M, start, goal_, LU);
    for (CapOpt& c : pb.caps) {  // 放大障碍，让穿透罚项处于激活状态
      c.r *= 1.6;
      c.cosr = std::cos(c.r);
      c.sinr = std::sin(c.r);
    }
    std::vector<Vec3> P0 = fitCtrl(M, resamplePolyline(poly, 4 * M + 1), start.p, goal_.p);
    std::vector<double> x = toVars(pb.bd, P0, 0.3 * LU / (3.0 * M));
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, 0.01);
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
    Vec3 ps, hs, ls, pg;          // 起点 (位置/航向/侧向) 与终点位置；终点航向不约束
    double kapS = 0.0;            // 单位球曲率
    bool pinS = true;
    int numVars() const { return 1 + (pinS ? 2 : 3) + 3 * (M - 3) + 6; }  // P_M 与 P_{M+1} 恒自由
  };
  struct CapOpt {
    Vec3 n;
    double r = 0, cosr = 1, sinr = 0;
  };
  // 采样表：每个 span 含 nq 个 Gauss-Legendre 节点（用于积分项与罚项）
  // + 节点 u=k（以及最后一个 span 的 u=M）处的附加点（仅用于罚项，防止端点处尖峰漏检）
  struct SampleTable {
    int M = 0, nq = 0;
    std::vector<int> offset;                // span k 的样本为 [offset[k], offset[k+1])
    std::vector<double> wInt, wPen;         // 积分项权重 / 罚项权重（u 测度）
    std::vector<double> u;                  // 每个样本的参数 u（曲率变化罚项用，全局按 u 递增）
    std::vector<std::array<double, 12>> N;  // N[d*4+j]
  };
  struct Problem {
    Boundary bd;
    const SampleTable* tab = nullptr;
    std::vector<CapOpt> caps;
    double LrefU = 1, kapRefU = 1, obsScaleU = 1, sMinU = 0;
    double wE = 1, wB = 0, wKv = 0, wSp = 0, wO = 0, wR = 0;
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
    const double as = std::exp(x[0]);
    int k = 1;
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
    P[M] = Vec3(x[k], x[k + 1], x[k + 2]);
    k += 3;
    P[M + 1] = Vec3(x[k], x[k + 1], x[k + 2]);
    P[M + 2] = b.pg;
  }

  static void chainGrad(const Boundary& b, const std::vector<double>& x, const std::vector<Vec3>& gP,
                        std::vector<double>& g) {
    const int M = b.M;
    g.assign(x.size(), 0.0);
    const double as = std::exp(x[0]);
    g[0] = gP[1].dot(b.hs) * as;
    int k = 1;
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
    put3(gP[M]);
    put3(gP[M + 1]);
  }

  static std::vector<double> toVars(const Boundary& b, const std::vector<Vec3>& P, double aMin) {
    const int M = b.M;
    std::vector<double> x(b.numVars());
    x[0] = std::log(std::max((P[1] - b.ps).dot(b.hs), aMin));
    int k = 1;
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
    put3(P[M]);
    put3(P[M + 1]);
    return x;
  }

  // ======================= 代价与梯度 =======================
  // 每采样点几何量与偏导（κ = 单位球测地曲率, s = 单位球弧长速率 ds/du）
  struct GeomPart {
    double kap = 0.0, s = 0.0;
    double dk_dr2 = 0, dk_dW2 = 0, dk_dd = 0;  // ∂κ/∂(r2, W2, d)
    double ds_dW2 = 0, ds_dr2 = 0;             // ∂s/∂(r2, W2)
    bool ok = false;                           // W2 > 0（非退化）
  };
  static GeomPart geomPart(const Vec3& c, const Vec3& c1, const Vec3& c2) {
    GeomPart gp;
    const double r2 = c.squaredNorm();
    const Vec3 wv = c.cross(c1);
    const double W2 = wv.squaredNorm();
    if (W2 <= 1e-24) return gp;
    const double r = std::sqrt(r2), W = std::sqrt(W2);
    const double d = wv.dot(c2);
    gp.kap = r2 * r * d / (W2 * W);
    gp.s = W / r2;
    gp.dk_dr2 = 1.5 * gp.kap / r2;
    gp.dk_dW2 = -1.5 * gp.kap / W2;
    gp.dk_dd = r2 * r / (W2 * W);
    gp.ds_dW2 = 0.5 / (W * r2);
    gp.ds_dr2 = -W / (r2 * r2);
    gp.ok = true;
    return gp;
  }

  // 点状代价（能量 + 曲率幅值软罚 + 参数速度下限）及其对 (c, c', c'') 的梯度（手推解析式）
  static double pointwiseTerm(const Vec3& c, const Vec3& c1, const Vec3& c2, const GeomPart& gp,
                              double aE, double aB, double aSp, double sMin,
                              Vec3& gc, Vec3& gc1, Vec3& gc2) {
    const double r2 = c.squaredNorm();
    const Vec3 wv = c.cross(c1);
    const double W2 = wv.squaredNorm();
    double J = aE * W2 / (r2 * r2);       // s² = W2 / r2²
    double dW2 = aE / (r2 * r2);
    double dr2 = -2.0 * aE * W2 / (r2 * r2 * r2);
    double dd = 0.0;
    if (gp.ok) {
      J += aB * gp.kap * gp.kap * gp.s;   // ∫κ² ds 的离散项
      double gk = aB * 2.0 * gp.kap * gp.s, gs = aB * gp.kap * gp.kap;
      const double vs = (sMin - gp.s) / sMin;  // 速度下限
      if (vs > 0) {
        J += aSp * vs * vs * vs;
        gs -= aSp * 3.0 * vs * vs / sMin;
      }
      dW2 += gk * gp.dk_dW2 + gs * gp.ds_dW2;
      dr2 += gk * gp.dk_dr2 + gs * gp.ds_dr2;
      dd += gk * gp.dk_dd;
    }
    // ∂r2/∂c = 2c；∂W2/∂c = 2 c'×w；∂W2/∂c' = 2 w×c；∂d/∂c = c'×c''；∂d/∂c' = c''×c；∂d/∂c'' = w
    gc = c * (2.0 * dr2) + c1.cross(wv) * (2.0 * dW2) + c1.cross(c2) * dd;
    gc1 = wv.cross(c) * (2.0 * dW2) + c2.cross(c) * dd;
    gc2 = wv * dd;
    return J;
  }

  // 曲率变化罚项的 (κ, s) 梯度回传到 (c, c', c'')
  static void kapSGrad(const Vec3& c, const Vec3& c1, const Vec3& c2, const Vec3& wv, const GeomPart& gp,
                       double gkap, double gs, Vec3& gc, Vec3& gc1, Vec3& gc2) {
    if (!gp.ok || (gkap == 0.0 && gs == 0.0)) {
      gc = gc1 = gc2 = Vec3();
      return;
    }
    const double gr2 = gkap * gp.dk_dr2 + gs * gp.ds_dr2;
    const double gW2 = gkap * gp.dk_dW2 + gs * gp.ds_dW2;
    const double gd = gkap * gp.dk_dd;
    gc = c * (2.0 * gr2) + c1.cross(wv) * (2.0 * gW2) + c1.cross(c2) * gd;
    gc1 = wv.cross(c) * (2.0 * gW2) + c2.cross(c) * gd;
    gc2 = wv * gd;
  }

  double costGrad(const Problem& pb, const std::vector<double>& x, std::vector<double>& g) {
    const int M = pb.bd.M;
    std::vector<Vec3>& P = bufP_;
    std::vector<Vec3>& gP = bufGP_;
    buildCtrl(pb.bd, x, P);
    gP.assign(M + 3, Vec3());
    double J = 0.0;
    // 径向正则（P_{M+1} 现为自由控制点, 一并拉回单位球; P_1 由边界式保持良态）
    for (int j = 2; j <= M + 1; ++j) {
      const double e = P[j].squaredNorm() - 1.0;
      J += pb.wR * e * e;
      gP[j] += P[j] * (4.0 * pb.wR * e);
    }
    // 采样点项（点状部分），并记录每采样点的 (κ, s) 供曲率变化罚项使用
    const SampleTable& tab = *pb.tab;
    const int nS = static_cast<int>(tab.N.size());
    std::vector<double>& kap = bufKap_;
    std::vector<double>& sp = bufS_;
    kap.assign(nS, 0.0);
    sp.assign(nS, 0.0);
    const double aE0 = pb.wE * M / (pb.LrefU * pb.LrefU);
    const double aB0 = pb.wB / (pb.kapRefU * pb.kapRefU * pb.LrefU);
    const double aSp0 = pb.wSp / M, aO0 = pb.wO / M;
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
        const GeomPart gp = geomPart(c, c1, c2);
        kap[q] = gp.kap;
        sp[q] = gp.s;
        Vec3 gc, gc1, gc2;
        J += pointwiseTerm(c, c1, c2, gp, aE0 * w, aB0 * w, aSp0 * wp, pb.sMinU, gc, gc1, gc2);
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
    // 曲率变化率软罚：J += (wKv/M)·Σ_q (Δκ_q/κref)²·(spanU/Δs_q)，
    // Δs_q = 0.5(s_q + s_{q+1})·Δu（单位球弧长增量）。κ 在一个 span 内变化一个 κref
    // 量级时整体贡献 ≈ wKv，与其它 O(1) 代价项同尺度。
    if (pb.wKv > 0 && nS > 1) {
      std::vector<double>& gk = bufGK_;
      std::vector<double>& gs = bufGS_;
      gk.assign(nS, 0.0);
      gs.assign(nS, 0.0);
      const double spanU = pb.LrefU / M;
      const double wVar = pb.wKv / M / (pb.kapRefU * pb.kapRefU);
      for (int q = 0; q + 1 < nS; ++q) {
        const double du = tab.u[q + 1] - tab.u[q];
        if (du <= 0.0) continue;
        const double dkap = kap[q + 1] - kap[q];
        const double dsU = 0.5 * (sp[q] + sp[q + 1]) * du;
        if (dsU <= 1e-12) continue;
        const double c0 = wVar * spanU / dsU;
        J += c0 * dkap * dkap;
        const double gd = 2.0 * c0 * dkap;               // ∂J/∂(Δκ)
        gk[q] -= gd;
        gk[q + 1] += gd;
        const double ges = -c0 * dkap * dkap / dsU;      // ∂J/∂(Δs)
        gs[q] += ges * 0.5 * du;
        gs[q + 1] += ges * 0.5 * du;
      }
      // (κ, s) 梯度回传到控制点
      for (int k = 0; k < M; ++k) {
        for (int q = tab.offset[k]; q < tab.offset[k + 1]; ++q) {
          if (gk[q] == 0.0 && gs[q] == 0.0) continue;
          const std::array<double, 12>& N = tab.N[q];
          Vec3 c, c1, c2;
          for (int j = 0; j < 4; ++j) {
            c += P[k + j] * N[j];
            c1 += P[k + j] * N[4 + j];
            c2 += P[k + j] * N[8 + j];
          }
          const GeomPart gp = geomPart(c, c1, c2);
          Vec3 gc, gc1, gc2;
          kapSGrad(c, c1, c2, c.cross(c1), gp, gk[q], gs[q], gc, gc1, gc2);
          for (int j = 0; j < 4; ++j) gP[k + j] += gc * N[j] + gc1 * N[4 + j] + gc2 * N[8 + j];
        }
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
      t.u.push_back(u);
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
    const double R = prm_.R, kapRefU = std::max(prm_.kappaMax * R, 1e-9);
    Problem pb;
    Boundary& b = pb.bd;
    b.M = M;
    b.ps = s.p;
    b.hs = s.h;
    b.ls = s.p.cross(s.h);
    b.kapS = clampd(s.kappa * R, -kapRefU, kapRefU);
    b.pg = g.p;
    b.pinS = prm_.pinStartCurvature;
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
    pb.kapRefU = kapRefU;
    pb.obsScaleU = prm_.obsScale / R;
    pb.sMinU = prm_.minSpeedFrac * pb.LrefU / M;
    pb.wE = prm_.wEnergy;
    pb.wB = prm_.wBend;
    pb.wKv = prm_.wKappaVar;
    pb.wSp = prm_.wSpeed;
    pb.wO = prm_.wObs;
    pb.wR = prm_.wRadial;
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
    const Vec3 ps2 = leadPoint(s.p, s.h);
    // 终端航向不约束：目标引导方向取起点指向目标的自然切向（仅前端初值）
    Vec3 gBack = s.p - g.p;
    gBack -= g.p * gBack.dot(g.p);  // 投影到目标处切平面
    const Vec3 gLead = gBack.squaredNorm() < 1e-18 ? -g.h : gBack.normalized();
    const Vec3 pg2 = leadPoint(g.p, gLead);
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

  static SurfaceState normalizeState(const SurfaceState& s) {
    SurfaceState o = s;
    o.p = s.p.normalized();
    o.h = projectTangent(s.h, o.p).normalized();
    return o;
  }
  static double msSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  }
  void fillStats(PlanStats& st, const SolveOut& so, const Validation& v, int M) const {
    st.iters = so.res.iters;
    st.evals = so.res.evals;
    st.cost = so.res.f;
    describeRetained(st, so.sp, v);
    st.spans = M;
  }
  // 几何字段改记实际保留的曲线。求解次数和代价留给本周期的候选解。
  static void describeRetained(PlanStats& st, const SphereSpline& sp, const Validation& v) {
    st.valid = sp.valid() && v.ok;
    st.length = sp.valid() ? sp.length() : 0.0;
    st.minClearance = sp.valid() ? v.minClearance : 0.0;
    st.maxKappa = sp.valid() ? v.maxKappa : 0.0;
    st.spans = sp.valid() ? sp.numSpans() : 0;
  }

  PlannerParams prm_;
  std::vector<Obstacle> obs_;
  SurfaceState goal_;
  bool hasGoal_ = false;
  SphereSpline curve_;
  double lastU_ = 0.0;
  std::map<std::pair<int, int>, SampleTable> tables_;
  std::map<std::pair<int, int>, FitCache> fitCaches_;
  std::vector<Vec3> bufP_, bufGP_;
  std::vector<int> bufNear_;
  std::vector<double> bufKap_, bufS_, bufGK_, bufGS_;  // 曲率变化罚项的每采样点缓冲
};

}  // namespace sphtraj
