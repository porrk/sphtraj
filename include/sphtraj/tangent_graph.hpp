#pragma once
// 前端：球面切线图（tangent graph）。在若干球冠障碍之间求起终点的最短无碰路径
// （由大圆切线段与障碍边界小圆弧组成，G1 连续，不考虑曲率约束）。
// 用途：确定同伦类 + 为后端 B 样条优化提供初值。
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <vector>

#include "sphere_math.hpp"

namespace sphtraj {

struct FrontEndPath {
  bool ok = false;
  std::vector<Vec3> points;  // 折线（相邻点间按大圆插值）
  double lengthU = 0.0;      // 单位球长度
};

namespace detail {

inline bool pointInCap(const Vec3& p, const Cap& c, double eps = 1e-7) {
  return angleBetween(p, c.center) < c.radius - eps;
}

// 小弧 a→b 是否进入球冠内部（半径略缩，允许相切）
inline bool arcHitsCap(const Vec3& a, const Vec3& b, const Cap& cap, double eps = 1e-7) {
  const double rEff = cap.radius - eps;
  if (rEff <= 0) return false;
  const Vec3& n = cap.center;
  Vec3 w = a.cross(b);
  const double wn = w.norm();
  double dist;
  if (wn < 1e-12) {
    dist = angleBetween(n, a);
  } else {
    w = w / wn;
    Vec3 q = n - w * n.dot(w);
    const double qn = q.norm();
    if (qn < 1e-12) {
      dist = kPi / 2;
    } else {
      q = q / qn;
      const bool onArc = a.cross(q).dot(w) >= 0 && q.cross(b).dot(w) >= 0;
      dist = onArc ? angleBetween(n, q) : std::min(angleBetween(n, a), angleBetween(n, b));
    }
  }
  return dist < rEff;
}

// 求单位向量 a，满足 a·n1 = A, a·n2 = B（最多两个解）。a 为某大圆的极点。
inline int solvePole(const Vec3& n1, const Vec3& n2, double A, double B, Vec3 out[2]) {
  const double c12 = n1.dot(n2);
  const Vec3 nx = n1.cross(n2);
  const double det = nx.squaredNorm();
  if (det < 1e-14) return 0;
  const double al = (A - B * c12) / det, be = (B - A * c12) / det;
  const Vec3 base = n1 * al + n2 * be;
  const double t2 = (1.0 - base.squaredNorm()) / det;
  if (t2 < 0) return 0;
  const double t = std::sqrt(t2);
  out[0] = (base + nx * t).normalized();
  out[1] = (base - nx * t).normalized();
  return 2;
}

// 极点为 a 的大圆与球冠的切点
inline Vec3 tangentPoint(const Cap& c, const Vec3& a) { return (c.center - a * a.dot(c.center)).normalized(); }

}  // namespace detail

class TangentGraph {
 public:
  FrontEndPath plan(const Vec3& s, const Vec3& g, const std::vector<Cap>& caps, double arcStep = deg2rad(1.0)) {
    using namespace detail;
    struct Node {
      Vec3 p;
      int cap;
      double phi;
    };
    struct Edge {
      int to;
      double w;
      int cap;  // <0：大圆段；>=0：沿该球冠边界的圆弧
      double dphi;
    };
    const int K = static_cast<int>(caps.size());
    std::vector<Vec3> e1(K), e2(K);
    for (int i = 0; i < K; ++i) tangentFrame(caps[i].center, e1[i], e2[i]);

    std::vector<Node> nodes;
    std::vector<std::vector<Edge>> adj;
    auto addNode = [&](const Vec3& p, int cap) -> int {
      for (int i = 0; i < K; ++i)
        if (i != cap && pointInCap(p, caps[i])) return -1;
      double phi = 0.0;
      if (cap >= 0) phi = std::atan2(p.dot(e2[cap]), p.dot(e1[cap]));
      nodes.push_back({p, cap, phi});
      adj.emplace_back();
      return static_cast<int>(nodes.size()) - 1;
    };
    auto segFree = [&](const Vec3& a, const Vec3& b) {
      for (int i = 0; i < K; ++i)
        if (arcHitsCap(a, b, caps[i])) return false;
      return true;
    };
    auto addSeg = [&](int i, int j) {
      if (i < 0 || j < 0) return;
      if (!segFree(nodes[i].p, nodes[j].p)) return;
      const double w = angleBetween(nodes[i].p, nodes[j].p);
      adj[i].push_back({j, w, -1, 0.0});
      adj[j].push_back({i, w, -1, 0.0});
    };

    FrontEndPath out;
    const int S = addNode(s, -1), G = addNode(g, -1);
    if (S < 0 || G < 0) return out;  // 起点或终点在障碍内（调用方应先处理）
    addSeg(S, G);

    // 起/终点到各球冠的切线
    for (int e : {S, G}) {
      for (int i = 0; i < K; ++i) {
        Vec3 a[2];
        const int ns = solvePole(nodes[e].p, caps[i].center, 0.0, std::sin(caps[i].radius), a);
        for (int t = 0; t < ns; ++t) addSeg(e, addNode(tangentPoint(caps[i], a[t]), i));
      }
    }
    // 球冠之间的公切线（外公切 + 内公切）
    for (int i = 0; i < K; ++i) {
      for (int j = i + 1; j < K; ++j) {
        for (double sj : {1.0, -1.0}) {
          Vec3 a[2];
          const int ns = solvePole(caps[i].center, caps[j].center, std::sin(caps[i].radius),
                                   sj * std::sin(caps[j].radius), a);
          for (int t = 0; t < ns; ++t) {
            const int ni = addNode(tangentPoint(caps[i], a[t]), i);
            const int nj = addNode(tangentPoint(caps[j], a[t]), j);
            addSeg(ni, nj);
          }
        }
      }
    }
    // 每个球冠上相邻切点之间的边界圆弧
    for (int c = 0; c < K; ++c) {
      std::vector<int> ids;
      for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        if (nodes[i].cap == c) ids.push_back(i);
      if (ids.size() < 2) continue;
      std::sort(ids.begin(), ids.end(), [&](int a, int b) { return nodes[a].phi < nodes[b].phi; });
      const double sr = std::sin(caps[c].radius), cr = std::cos(caps[c].radius);
      for (size_t t = 0; t < ids.size(); ++t) {
        const int a = ids[t], b = ids[(t + 1) % ids.size()];
        double dphi = nodes[b].phi - nodes[a].phi;
        if (dphi <= 0) dphi += 2.0 * kPi;
        const int ns = std::max(2, static_cast<int>(std::ceil(dphi * sr / arcStep)));
        bool free = true;
        for (int q = 1; q < ns && free; ++q) {
          const double ph = nodes[a].phi + dphi * q / ns;
          const Vec3 pt = caps[c].center * cr + (e1[c] * std::cos(ph) + e2[c] * std::sin(ph)) * sr;
          for (int o = 0; o < K; ++o)
            if (o != c && pointInCap(pt, caps[o])) { free = false; break; }
        }
        if (!free) continue;
        adj[a].push_back({b, sr * dphi, c, dphi});
        adj[b].push_back({a, sr * dphi, c, -dphi});
      }
    }

    // Dijkstra
    const int Nn = static_cast<int>(nodes.size());
    std::vector<double> dist(Nn, std::numeric_limits<double>::infinity());
    std::vector<int> par(Nn, -1), parEdge(Nn, -1);
    using QE = std::pair<double, int>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
    dist[S] = 0.0;
    pq.push({0.0, S});
    while (!pq.empty()) {
      const auto [d, u] = pq.top();
      pq.pop();
      if (d > dist[u]) continue;
      if (u == G) break;
      for (int ei = 0; ei < static_cast<int>(adj[u].size()); ++ei) {
        const Edge& e = adj[u][ei];
        const double nd = d + e.w;
        if (nd < dist[e.to]) {
          dist[e.to] = nd;
          par[e.to] = u;
          parEdge[e.to] = ei;
          pq.push({nd, e.to});
        }
      }
    }
    if (!std::isfinite(dist[G])) return out;

    // 回溯并展开圆弧
    std::vector<int> seq;
    for (int v = G; v >= 0; v = par[v]) seq.push_back(v);
    std::reverse(seq.begin(), seq.end());
    out.points.push_back(nodes[S].p);
    for (size_t t = 1; t < seq.size(); ++t) {
      const int u = seq[t - 1], v = seq[t];
      const Edge& e = adj[u][parEdge[v]];
      if (e.cap >= 0) {
        const int c = e.cap;
        const double sr = std::sin(caps[c].radius), cr = std::cos(caps[c].radius);
        const int ns = std::max(2, static_cast<int>(std::ceil(std::fabs(e.dphi) * sr / arcStep)));
        for (int q = 1; q < ns; ++q) {
          const double ph = nodes[u].phi + e.dphi * q / ns;
          out.points.push_back(caps[c].center * cr + (e1[c] * std::cos(ph) + e2[c] * std::sin(ph)) * sr);
        }
      }
      out.points.push_back(nodes[v].p);
    }
    out.lengthU = dist[G];
    out.ok = true;
    return out;
  }
};

}  // namespace sphtraj
