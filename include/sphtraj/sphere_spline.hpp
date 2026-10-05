#pragma once
// 方案 A 的曲线表示：R^3 中的夹持三次 B 样条 c(u)，经中心投影 γ(u) = c/|c| 落到球面。
// 提供：点/切向/测地曲率、弧长表、弧长<->参数互查、最近点投影、跟踪参考量、球面纯追踪。
#include <algorithm>
#include <cmath>
#include <vector>

#include "bspline_basis.hpp"
#include "sphere_math.hpp"

namespace sphtraj {

class SphereSpline {
 public:
  struct Ambient {
    Vec3 c, c1, c2;  // c, dc/du, d2c/du2
  };
  struct Geom {
    Vec3 pos;           // γ（单位向量）
    Vec3 tangent;       // 单位切向 T
    Vec3 left;          // 左法向 γ × T
    Vec3 d1, d2;        // dγ/du, d²γ/du²（单位球）
    double speedU = 0;  // |dγ/du|（单位球）
    double kappaU = 0;  // 测地曲率（单位球，左转为正）；物理曲率 = kappaU / R
    double cnorm = 0;   // |c|，用于退化检测
  };

  void set(int M, const std::vector<Vec3>& ctrl, double R) {
    M_ = M;
    P_ = ctrl;
    R_ = R;
    basis_.reset(M);
    gaussLegendre01(8, glx_, glw_);
    cum_.assign(M + 1, 0.0);
    for (int k = 0; k < M; ++k) cum_[k + 1] = cum_[k] + speedIntegral(k, k, k + 1.0);
  }

  bool valid() const { return M_ > 0; }
  int numSpans() const { return M_; }
  double uMax() const { return double(M_); }
  double radius() const { return R_; }
  const std::vector<Vec3>& ctrl() const { return P_; }
  double lengthU() const { return cum_.empty() ? 0.0 : cum_.back(); }
  double length() const { return R_ * lengthU(); }  // 物理总长
  double spanLengthU(int k) const { return cum_[k + 1] - cum_[k]; }

  Ambient ambient(double u) const {
    const int k = basis_.spanOf(u);
    double N[3][4];
    basis_.eval(u, k, N);
    Ambient a;
    for (int j = 0; j < 4; ++j) {
      a.c += P_[k + j] * N[0][j];
      a.c1 += P_[k + j] * N[1][j];
      a.c2 += P_[k + j] * N[2][j];
    }
    return a;
  }

  // 由 (c, c', c'') 闭式计算球面几何量
  static Geom geomFromAmbient(const Ambient& a) {
    Geom g;
    const double r = a.c.norm();
    g.cnorm = r;
    g.pos = a.c / r;
    const double rp = g.pos.dot(a.c1);                                       // |c|'
    g.d1 = (a.c1 - g.pos * rp) / r;                                           // γ'
    const double rpp = (a.c1.squaredNorm() + a.c.dot(a.c2) - rp * rp) / r;  // |c|''
    g.d2 = (a.c2 - g.d1 * (2.0 * rp) - g.pos * rpp) / r;                    // γ''
    const Vec3 w = a.c.cross(a.c1);
    const double W = w.norm();
    g.speedU = W / (r * r);                                    // |γ'| = |c×c'|/|c|²
    g.tangent = g.speedU > 0 ? g.d1 / g.d1.norm() : Vec3();
    g.left = g.pos.cross(g.tangent);
    g.kappaU = W > 1e-15 ? r * r * r * w.dot(a.c2) / (W * W * W) : 0.0;  // |c|³ (c×c')·c'' / |c×c'|³
    return g;
  }
  Geom geom(double u) const { return geomFromAmbient(ambient(u)); }
  Vec3 pos(double u) const { return ambient(u).c.normalized(); }
  double speedU(double u) const {
    const Ambient a = ambient(u);
    return a.c.cross(a.c1).norm() / a.c.squaredNorm();
  }

  // 参数 -> 物理弧长
  double sOfU(double u) const {
    const int k = basis_.spanOf(u);
    return R_ * (cum_[k] + speedIntegral(k, double(k), u));
  }

  // 物理弧长 -> 参数（二分定位 span + Newton）
  double uOfS(double s) const {
    const double sU = clampd(s / R_, 0.0, lengthU());
    int k = static_cast<int>(std::upper_bound(cum_.begin(), cum_.end(), sU) - cum_.begin()) - 1;
    k = std::clamp(k, 0, M_ - 1);
    const double t = sU - cum_[k], L = spanLengthU(k);
    double u = k + (L > 0 ? t / L : 0.0);
    for (int it = 0; it < 8; ++it) {
      const double f = speedIntegral(k, double(k), u) - t;
      const double sp = speedU(u);
      if (sp < 1e-14) break;
      const double un = clampd(u - f / sp, double(k), double(k + 1));
      const bool done = std::fabs(un - u) < 1e-13;
      u = un;
      if (done) break;
    }
    return u;
  }

  // 最近点投影：在窗口 [uGuess-back, uGuess+fwd] 内最大化 γ(u)·p（粗搜 + 带保护的 Newton）
  double closestU(const Vec3& p, double uGuess, double back = 0.5, double fwd = 2.0) const {
    const double lo = std::max(0.0, uGuess - back), hi = std::min(uMax(), uGuess + fwd);
    const int n = std::max(2, static_cast<int>(std::ceil((hi - lo) / 0.1)));
    double u = lo, best = -2.0;
    for (int i = 0; i <= n; ++i) {
      const double ui = lo + (hi - lo) * i / n;
      const double v = pos(ui).dot(p);
      if (v > best) { best = v; u = ui; }
    }
    for (int it = 0; it < 12; ++it) {
      const Geom g = geom(u);
      const double f1 = g.d1.dot(p), f2 = g.d2.dot(p), cur = g.pos.dot(p);
      double du = (f2 < -1e-14) ? -f1 / f2 : (f1 > 0 ? 0.05 : -0.05);
      du = clampd(du, -0.25, 0.25);
      double un = clampd(u + du, lo, hi);
      for (int h = 0; h < 12 && pos(un).dot(p) < cur - 1e-15; ++h) {
        du *= 0.5;
        un = clampd(u + du, lo, hi);
      }
      const bool done = std::fabs(un - u) < 1e-11;
      u = un;
      if (done) break;
    }
    return u;
  }

  // 在物理弧长区间 [s0, s1] 上等弧长采样 K 个点
  std::vector<Vec3> sampleUniformS(double s0, double s1, int K) const {
    std::vector<Vec3> out(K);
    for (int i = 0; i < K; ++i) out[i] = pos(uOfS(s0 + (s1 - s0) * i / double(K - 1)));
    return out;
  }

  // 在参数区间 [u0, u1] 上等参数采样 K 个点（热启动用：保留旧解的参数化，拟合几乎无损）
  std::vector<Vec3> sampleUniformU(double u0, double u1, int K) const {
    std::vector<Vec3> out(K);
    for (int i = 0; i < K; ++i) out[i] = pos(u0 + (u1 - u0) * i / double(K - 1));
    return out;
  }

  // 从 u0 到终点，按角步长不超过 maxStepU 采样位置
  std::vector<Vec3> samplePositions(double u0, double maxStepU) const {
    std::vector<Vec3> out;
    u0 = clampd(u0, 0.0, uMax());
    out.push_back(pos(u0));
    for (int k = std::clamp(static_cast<int>(std::floor(u0)), 0, M_ - 1); k < M_; ++k) {
      const int n = std::max(1, static_cast<int>(std::ceil(spanLengthU(k) / maxStepU)));
      for (int i = 1; i <= n; ++i) {
        const double u = k + double(i) / n;
        if (u > u0) out.push_back(pos(u));
      }
    }
    return out;
  }

 private:
  // span k 内 ∫_a^b |γ'(u)| du（单位球）。8 点 Gauss-Legendre + 自适应二分：
  // 对接近匀速的优化结果通常一次通过；速度近零（近尖点）的曲线会自动细分。
  double speedIntegral(int /*k*/, double a, double b) const {
    if (b <= a) return 0.0;
    return adaptiveGL(a, b, gl8(a, b), 0);
  }
  double gl8(double a, double b) const {
    double sum = 0.0;
    for (size_t q = 0; q < glx_.size(); ++q) sum += glw_[q] * speedU(a + (b - a) * glx_[q]);
    return sum * (b - a);
  }
  double adaptiveGL(double a, double b, double whole, int depth) const {
    const double m = 0.5 * (a + b);
    const double l = gl8(a, m), r = gl8(m, b);
    if (depth >= 12 || std::fabs(l + r - whole) <= 1e-13 * std::max(1e-3, std::fabs(whole))) return l + r;
    return adaptiveGL(a, m, l, depth + 1) + adaptiveGL(m, b, r, depth + 1);
  }

  int M_ = 0;
  double R_ = 1.0;
  std::vector<Vec3> P_;
  ClampedCubicBasis basis_;
  std::vector<double> cum_, glx_, glw_;
};

// 球面纯追踪：过 p、在 p 处与 h 相切、且经过前视点 q 的小圆的测地曲率（物理单位，左转为正）。
// 平面极限下退化为 2 sin(alpha) / Ld。
inline double purePursuitCurvature(const Vec3& p, const Vec3& h, const Vec3& q, double R) {
  const double den = 1.0 - p.dot(q);
  if (den < 1e-14) return 0.0;
  return q.dot(p.cross(h)) / (R * den);
}

struct TrackingRef {
  double u = 0, s = 0;          // 最近点参数与物理弧长
  Vec3 refPos, refTangent;      // 最近点及其切向
  double crossTrack = 0;        // 横向误差（物理），车辆在路径左侧为正
  double headingError = 0;      // 车辆航向相对路径切向（rad），左偏为正
  double kappaRef = 0;          // 路径曲率前馈（1/物理长度）
  Vec3 lookaheadPos;            // 前视点
  double lookaheadS = 0;
  double kappaPurePursuit = 0;  // 纯追踪曲率指令
  double remaining = 0;         // 剩余弧长（物理）
};

inline TrackingRef computeTrackingRef(const SphereSpline& sp, const Vec3& p, const Vec3& h, double uGuess,
                                      double lookahead) {
  TrackingRef r;
  const double R = sp.radius();
  r.u = sp.closestU(p, uGuess, 0.5, 2.0);
  const SphereSpline::Geom g = sp.geom(r.u);
  r.s = sp.sOfU(r.u);
  r.refPos = g.pos;
  r.refTangent = g.tangent;
  r.crossTrack = R * std::asin(clampd(p.dot(g.left), -1.0, 1.0));
  const Vec3 Tp = projectTangent(g.tangent, p).normalized();
  r.headingError = std::atan2(Tp.cross(h).dot(p), Tp.dot(h));
  r.kappaRef = g.kappaU / R;
  // 前视点超出终点时沿终点切向大圆外延，使车辆以目标航向到达终点
  r.lookaheadS = r.s + lookahead;
  if (r.lookaheadS <= sp.length()) {
    r.lookaheadPos = sp.pos(sp.uOfS(r.lookaheadS));
  } else {
    const SphereSpline::Geom ge = sp.geom(sp.uMax());
    r.lookaheadPos = geodesicStep(ge.pos, ge.tangent, (r.lookaheadS - sp.length()) / R);
  }
  r.kappaPurePursuit = purePursuitCurvature(p, h, r.lookaheadPos, R);
  r.remaining = sp.length() - r.s;
  return r;
}

}  // namespace sphtraj
