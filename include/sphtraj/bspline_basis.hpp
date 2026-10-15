#pragma once
// 夹持(clamped)均匀三次 B 样条基函数（含 1、2 阶导数）与 Gauss-Legendre 求积节点。
// 节点向量：{0,0,0,0,1,2,...,M-1,M,M,M,M}，共 M 个 span、M+3 个控制点，参数 u ∈ [0, M]。
#include <algorithm>
#include <cmath>
#include <vector>

namespace sphtraj {

// [0,1] 上的 n 点 Gauss-Legendre 节点与权重（权重和为 1）
inline void gaussLegendre01(int n, std::vector<double>& x, std::vector<double>& w) {
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5));
    double dp = 1.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (std::fabs(dz) < 1e-15) break;
    }
    x[i] = 0.5 * (1.0 - z);
    w[i] = 1.0 / ((1.0 - z * z) * dp * dp);
  }
}

class ClampedCubicBasis {
 public:
  ClampedCubicBasis() = default;
  explicit ClampedCubicBasis(int M) { reset(M); }

  void reset(int M) {
    M_ = M;
    knots_.assign(M + 7, 0.0);
    for (int i = 0; i < M + 7; ++i) knots_[i] = std::clamp(double(i - 3), 0.0, double(M));
  }
  int numSpans() const { return M_; }
  int numCtrl() const { return M_ + 3; }

  // 返回 u 所在 span k（作用于控制点 k..k+3），并把 u 夹持到 [0, M]。
  // M==0 时没有合法 span，返回 0 且不调用上下界颠倒的 clamp。
  int spanOf(double& u) const {
    if (M_ <= 0) {
      u = 0.0;
      return 0;
    }
    u = std::clamp(u, 0.0, double(M_));
    int k = static_cast<int>(std::floor(u));
    if (k < 0) k = 0;
    if (k >= M_) k = M_ - 1;
    return k;
  }

  // ders[d][j]：第 d 阶导数(d=0,1,2)的基函数值，对应控制点 k+j
  void eval(double u, int k, double ders[3][4]) const { dersBasis(k + 3, u, ders); }

 private:
  // The NURBS Book, Algorithm A2.3（p=3, n=2）
  void dersBasis(int i, double u, double ders[3][4]) const {
    constexpr int p = 3, nd = 2;
    const std::vector<double>& U = knots_;
    double ndu[p + 1][p + 1], left[p + 1], right[p + 1], a[2][p + 1];
    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
      left[j] = u - U[i + 1 - j];
      right[j] = U[i + j] - u;
      double saved = 0.0;
      for (int r = 0; r < j; ++r) {
        ndu[j][r] = right[r + 1] + left[j - r];
        const double temp = ndu[r][j - 1] / ndu[j][r];
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];
    for (int r = 0; r <= p; ++r) {
      int s1 = 0, s2 = 1;
      a[0][0] = 1.0;
      for (int k = 1; k <= nd; ++k) {
        double d = 0.0;
        const int rk = r - k, pk = p - k;
        if (r >= k) {
          a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
          d = a[s2][0] * ndu[rk][pk];
        }
        const int j1 = (rk >= -1) ? 1 : -rk;
        const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
        for (int j = j1; j <= j2; ++j) {
          a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
          d += a[s2][j] * ndu[rk + j][pk];
        }
        if (r <= pk) {
          a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
          d += a[s2][k] * ndu[r][pk];
        }
        ders[k][r] = d;
        std::swap(s1, s2);
      }
    }
    int r = p;
    for (int k = 1; k <= nd; ++k) {
      for (int j = 0; j <= p; ++j) ders[k][j] *= r;
      r *= (p - k);
    }
  }

  int M_ = 0;
  std::vector<double> knots_;
};

}  // namespace sphtraj
