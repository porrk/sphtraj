#pragma once
// 球面几何基础工具：轻量 Vec3、旋转、slerp、球冠、球极投影。
// 约定：内部一律在单位球上计算，位置用单位向量表示；物理长度 = 角度 * R。
#include <algorithm>
#include <cmath>

namespace sphtraj {

constexpr double kPi = 3.14159265358979323846;
inline double deg2rad(double d) { return d * kPi / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / kPi; }
inline double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 cross(const Vec3& o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
  double squaredNorm() const { return dot(*this); }
  double norm() const { return std::sqrt(squaredNorm()); }
  Vec3 normalized() const {
    const double n = norm();
    return n > 0 ? *this / n : Vec3();
  }
};
inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

// 两向量夹角（atan2 形式，在 0 与 pi 附近都数值稳健）
inline double angleBetween(const Vec3& a, const Vec3& b) { return std::atan2(a.cross(b).norm(), a.dot(b)); }

// v 在 p 处切平面上的投影
inline Vec3 projectTangent(const Vec3& v, const Vec3& p) { return v - p * v.dot(p); }

// 与单位向量 n 正交的右手标架 (e1, e2, n)
inline void tangentFrame(const Vec3& n, Vec3& e1, Vec3& e2) {
  const Vec3 t = std::fabs(n.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
  e1 = projectTangent(t, n).normalized();
  e2 = n.cross(e1);
}

// Rodrigues 旋转：绕单位轴 k 旋转 ang
inline Vec3 rotate(const Vec3& v, const Vec3& k, double ang) {
  const double c = std::cos(ang), s = std::sin(ang);
  return v * c + k.cross(v) * s + k * (k.dot(v) * (1.0 - c));
}

// 从 p 出发沿初始方向 h（单位切向量）走大圆角距离 ang
inline Vec3 geodesicStep(const Vec3& p, const Vec3& h, double ang) {
  return (p * std::cos(ang) + h * std::sin(ang)).normalized();
}

inline Vec3 slerp(const Vec3& a, const Vec3& b, double t) {
  const double th = angleBetween(a, b);
  if (th < 1e-12) return a;
  const double s = std::sin(th);
  return (a * (std::sin((1.0 - t) * th) / s) + b * (std::sin(t * th) / s)).normalized();
}

// 单位球上的测地圆盘（球冠）：中心单位向量 + 角半径 (rad)
struct Cap {
  Vec3 center;
  double radius = 0.0;
};

// 以 pole 为投影极点的球极投影（用于同伦签名）
struct StereoChart {
  Vec3 pole, e1, e2;
  explicit StereoChart(const Vec3& z) : pole(z.normalized()) { tangentFrame(pole, e1, e2); }
  bool project(const Vec3& x, double& u, double& v) const {
    const double c = x.dot(pole);
    const double denom = 1.0 - c;
    if (denom < 1e-12) return false;
    const Vec3 xp = x - pole * c;
    u = xp.dot(e1) / denom;
    v = xp.dot(e2) / denom;
    return true;
  }
};

// ---- 经纬度工具（仅演示/测试用）----
inline Vec3 fromLatLon(double latDeg, double lonDeg) {
  const double la = deg2rad(latDeg), lo = deg2rad(lonDeg);
  return {std::cos(la) * std::cos(lo), std::cos(la) * std::sin(lo), std::sin(la)};
}
inline void latLonOf(const Vec3& p, double& latDeg, double& lonDeg) {
  latDeg = rad2deg(std::asin(clampd(p.z, -1.0, 1.0)));
  lonDeg = rad2deg(std::atan2(p.y, p.x));
}
// 方位角（自北向东，度）→ 单位切向量；极点处无定义
inline Vec3 headingFromAzimuth(const Vec3& p, double azDeg) {
  const Vec3 north = projectTangent(Vec3(0, 0, 1), p).normalized();
  const Vec3 east = north.cross(p);
  const double a = deg2rad(azDeg);
  return north * std::cos(a) + east * std::sin(a);
}

}  // namespace sphtraj
