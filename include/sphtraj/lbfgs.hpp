#pragma once
// 精简 L-BFGS：弱 Wolfe 线搜索（Lewis-Overton 二分/扩张，与 LBFGS-Lite 同类），无外部依赖。
// 代价函数签名：double fg(const std::vector<double>& x, std::vector<double>& grad)
// 若工程中已有 Eigen，可直接换成 LBFGS-Lite (ZJU-FAST-Lab) 或 LBFGSpp，代价函数无需改动逻辑。
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace sphtraj {

struct LbfgsParams {
  int mem = 8;              // 历史对数
  int maxIter = 200;        // 最大迭代
  double gTol = 1e-6;       // |g| <= gTol * max(1,|x|)
  double fRelTol = 1e-7;    // past 次迭代内相对下降量阈值
  int past = 3;
  int maxLineSearch = 40;
  double c1 = 1e-4;         // Armijo
  double c2 = 0.9;          // 弱曲率条件
  double maxStep = std::numeric_limits<double>::infinity();  // 每次迭代 |Δx|∞ 上限（防止大步跳入尖点/打结形状）
};

struct LbfgsResult {
  double f = 0.0;
  int iters = 0;
  int evals = 0;
  int status = 2;  // 0: 梯度收敛, 1: 函数值收敛, 2: 达到迭代上限, -1: 线搜索失败
};

template <class F>
LbfgsResult lbfgsMinimize(F&& fg, std::vector<double>& x, const LbfgsParams& prm) {
  const int n = static_cast<int>(x.size());
  const int m = std::max(1, prm.mem);
  const int past = std::max(1, prm.past);
  LbfgsResult res;
  std::vector<double> g(n), xt(n), gt(n), d(n), xlo(n), glo(n);
  std::vector<std::vector<double>> S(m, std::vector<double>(n)), Y(m, std::vector<double>(n));
  std::vector<double> rho(m, 0.0), alpha(m, 0.0), hist(past, 0.0);
  auto dot = [n](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
  };

  double f = fg(x, g);
  res.evals = 1;
  res.f = f;
  hist[0] = f;
  double gn = std::sqrt(dot(g, g));
  if (gn <= prm.gTol * std::max(1.0, std::sqrt(dot(x, x)))) {
    res.status = 0;
    res.iters = 0;
    return res;
  }
  int stored = 0, head = 0;
  for (int i = 0; i < n; ++i) d[i] = -g[i];
  double step = 1.0 / gn;

  for (int it = 1; it <= prm.maxIter; ++it) {
    res.iters = it;
    double gd = dot(g, d);
    if (!(gd < 0.0)) {  // 非下降方向：重置为最速下降
      for (int i = 0; i < n; ++i) d[i] = -g[i];
      gd = dot(g, d);
      stored = 0;
      step = 1.0 / std::sqrt(-gd);
    }
    // ---- 弱 Wolfe 线搜索（步长受 maxStep 限制）----
    double dInf = 0.0;
    for (int i = 0; i < n; ++i) dInf = std::max(dInf, std::fabs(d[i]));
    const double tMax = prm.maxStep / std::max(dInf, 1e-300);
    double lo = 0.0, hi = std::numeric_limits<double>::infinity(), t = std::min(step, tMax), ft = f, flo = f;
    bool ok = false, haveLo = false;
    for (int ls = 0; ls < prm.maxLineSearch; ++ls) {
      for (int i = 0; i < n; ++i) xt[i] = x[i] + t * d[i];
      ft = fg(xt, gt);
      ++res.evals;
      if (!std::isfinite(ft) || ft > f + prm.c1 * t * gd) {
        hi = t;
      } else if (dot(gt, d) < prm.c2 * gd) {
        lo = t;
        haveLo = true;
        xlo = xt;
        glo = gt;
        flo = ft;
        if (t >= tMax) break;  // 已达步长上限且满足 Armijo：接受
      } else {
        ok = true;
        break;
      }
      t = std::isfinite(hi) ? 0.5 * (lo + hi) : std::min(2.0 * t, tMax);
    }
    if (!ok) {
      if (!haveLo) {
        res.status = -1;
        break;
      }
      xt = xlo;
      gt = glo;
      ft = flo;
    }
    // ---- 更新历史 ----
    std::vector<double>& s = S[head];
    std::vector<double>& y = Y[head];
    for (int i = 0; i < n; ++i) {
      s[i] = xt[i] - x[i];
      y[i] = gt[i] - g[i];
    }
    const double sy = dot(s, y), yy = dot(y, y);
    x.swap(xt);
    g.swap(gt);
    f = ft;
    res.f = f;
    if (sy > 1e-12 * yy && yy > 0.0) {
      rho[head] = 1.0 / sy;
      head = (head + 1) % m;
      stored = std::min(stored + 1, m);
    }
    // ---- 收敛判据 ----
    gn = std::sqrt(dot(g, g));
    if (gn <= prm.gTol * std::max(1.0, std::sqrt(dot(x, x)))) {
      res.status = 0;
      break;
    }
    if (it >= past) {
      const double fp = hist[it % past];
      if ((fp - f) / std::max(1.0, std::fabs(f)) < prm.fRelTol) {
        res.status = 1;
        break;
      }
    }
    hist[it % past] = f;
    // ---- 双循环递推：d = -H g ----
    for (int i = 0; i < n; ++i) d[i] = g[i];
    int idx = head;
    for (int k = 0; k < stored; ++k) {
      idx = (idx - 1 + m) % m;
      alpha[idx] = rho[idx] * dot(S[idx], d);
      for (int i = 0; i < n; ++i) d[i] -= alpha[idx] * Y[idx][i];
    }
    if (stored > 0) {
      const int last = (head - 1 + m) % m;
      const double gamma = 1.0 / (rho[last] * dot(Y[last], Y[last]));
      for (int i = 0; i < n; ++i) d[i] *= gamma;
    }
    for (int k = 0; k < stored; ++k) {
      const double beta = rho[idx] * dot(Y[idx], d);
      for (int i = 0; i < n; ++i) d[i] += S[idx][i] * (alpha[idx] - beta);
      idx = (idx + 1) % m;
    }
    for (int i = 0; i < n; ++i) d[i] = -d[i];
    step = stored > 0 ? 1.0 : 1.0 / std::max(gn, 1e-300);
  }
  return res;
}

}  // namespace sphtraj
