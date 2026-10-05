# sphtraj：球面中心投影 B 样条轨迹优化（方案 A）

纯头文件 C++17 实现，**无第三方依赖**（自带精简 L-BFGS）。包含：切线图前端、B 样条后端优化、每周期重规划（拼接 + 热启动 + 近端项 + 同伦守卫）、跟踪接口（弧长、最近点、前视点、球面纯追踪）。

## 目录

```
include/sphtraj/
  sphere_math.hpp      Vec3、旋转、slerp、球冠、球极投影
  bspline_basis.hpp    夹持均匀三次 B 样条基函数(含 1/2 阶导) + Gauss-Legendre
  lbfgs.hpp            精简 L-BFGS（弱 Wolfe 线搜索，带最大步长）
  sphere_spline.hpp    曲线 γ=c/|c|：几何闭式、弧长表、s<->u、最近点、纯追踪
  tangent_graph.hpp    前端：球面切线图 + Dijkstra
  sphere_planner.hpp   后端：代价/解析梯度、拟合、重规划、校验、同伦签名
examples/demo_replan.cpp  闭环仿真演示
tests/test_core.cpp       自检（梯度、边界条件精确性、几何公式等）
tools/plot_demo.py        画图
```

## 编译运行

```bash
g++ -std=c++17 -O2 -Iinclude tests/test_core.cpp -o test_core && ./test_core
g++ -std=c++17 -O2 -Iinclude examples/demo_replan.cpp -o demo_replan && ./demo_replan out
python3 tools/plot_demo.py out        # 生成 out/demo.png
```

也提供了 CMakeLists.txt（`cmake -B build && cmake --build build`）。

## 用法

```cpp
#include "sphtraj/sphere_planner.hpp"
using namespace sphtraj;

PlannerParams prm;
prm.R = 1000.0;            // 球半径
prm.kappaMax = 1.0 / 40;   // 最大测地曲率 = 1/最小转弯半径
prm.spanLength = 50.0;     // 每个 knot span 的期望弧长（≈ 1~1.5 倍转弯半径）
prm.optMargin = 10.0;      // 优化膨胀
prm.checkMargin = 4.0;     // 校验膨胀（< optMargin，留给跟踪误差）
SphereTrajPlanner planner(prm);

planner.setObstacles(obstacles);                  // {中心单位向量, 物理测地半径}
planner.setGoal({pg, hg, 0.0});                   // 终点位置/航向/曲率
planner.planFromScratch({p0, h0, kappa0});        // 初始化（前端 + 优化）

// 每个控制周期：
PlanStats st = planner.replan(state, v * latency);       // state = {p, h, 当前曲率}
TrackingRef ref = planner.trackingReference(state, Ld);  // 横向误差/航向误差/前视点/曲率指令
// 低频（1~5 Hz，可放后台线程）：
planner.considerClassSwitch(state);

// 轨迹查询：
const SphereSpline& c = planner.curve();
double L = c.length();                     // 总长
auto g = c.geom(c.uOfS(s));                // 弧长 s 处的位置 g.pos / 切向 g.tangent / 曲率 g.kappaU / R
```

位置一律用单位向量表示，长度和曲率用物理单位（内部除以 / 乘以 R）。

## 与需求的对应关系

1. **起终点航向约束**：P0=ps, P1=ps+a_s·hs；P_{M+2}=pg, P_{M+1}=pg−a_g·hg。航向对任意 a_s, a_g>0 精确成立，测试中误差约 1e-15 rad。可选地固定 ℓs·P2 = 3a_s²κ0，使起点曲率等于车辆当前曲率，转向指令跨周期连续；终点曲率同理。
2. **最短**：能量项 E = M∫|γ'|²du/Lref² ≥ (L/Lref)²，取等号当且仅当匀速。因此它同时最小化长度，并使参数近似弧长参数化。
3. **平滑**：曲线 C² 连续（曲率连续）。曲率上限以罚项实现（优化用 0.9·κmax，校验用 1.05·κmax）。另有弯曲能量正则，以及速度下限罚项（防尖点）。罚项样本除 Gauss-Legendre 点外，还包含每个节点和两个端点，避免端点处曲率尖峰漏检。
4. **同伦一致**：
   - 拼接：误差小时，起点状态取自旧轨迹（含曲率），而不是噪声测量。
   - 热启动：旧轨迹剩余段按参数重采样，再用缓存的 Cholesky 做最小二乘拟合。
   - 近端项：约束相邻周期之间的控制点变化。
   - 每步位移上限：L-BFGS 每次迭代的控制点位移有上限。
   - 同伦守卫：用球面 H-signature 比较新解与旧轨迹剩余段，类别改变且旧解仍可行时保留旧解。
   - 主动切换：前端低频评估，带收益阈值和连续次数迟滞。
5. **实时**：障碍剔除利用凸包性质（span 落在其 4 个控制点的包围球冠内，剔除是精确的）。梯度是手推的解析式，只依赖每个样本的 (c, c', c'')，再散射到 4 个控制点。
6. **跟踪**：
   - 弧长：8 点 Gauss-Legendre 加自适应细分。
   - s→u：二分加 Newton。
   - 最近点：窗口内粗搜，再用带保护的 Newton。
   - 横向误差和航向误差带符号。
   - 前视点越过终点时沿终点切向外延。
   - 球面纯追踪曲率闭式：κ = q·(p×h) / (R(1−p·q))。

## 演示结果（本仓库沙箱，单核）

场景：R = 1000 m，κmax = 1/40 m⁻¹，8 个静态障碍，中途在路径前方 150 m 处突然加入 1 个障碍，侧向漂移 0.25 m/s，航向噪声 0.3°/步，20 Hz 每周期重规划。

| 指标 | 结果 |
|---|---|
| 从零规划（22 span） | 约 3.5 ms |
| 每周期重规划 | 均值 0.31 ms，p99 0.95 ms，平均 8 次 L-BFGS 迭代 |
| 相邻周期规划偏差（前 200 m） | p50 1 cm，p99 0.21 m；新障碍出现那一周期为 31.9 m（预期内的绕行） |
| 横向误差最大值 | 1.9 m |
| 到达终点 | 位置误差 1.4 m，航向误差 3.9° |

实际耗时取决于 CPU、span 数和附近障碍数，请在目标平台上测。

## 调参建议

- `spanLength`：约 1~1.5 倍转弯半径。过小会导致变量多、收敛慢；过大则曲率分辨率不足。每个 span 对应的角度应 ≲ 10°，保证 |c| 远离 0。
- `wObs` / `wCurv`：默认 1e4 / 1e3（已归一化）。如果校验经常失败，先调大它们；`planFromScratch` 内部也会自动 ×10 重试两次。
- `wProx`：调大更"稳"（相邻周期更一致），但对新障碍的反应更慢。
- `lbfgsReplan.maxIter`：40 次以内通常足够；它也是最坏情况耗时的上限。
- `checkStep`：校验采样步长，应小于最小障碍直径的一半，防止"穿隧"。

## 替换为 LBFGS-Lite / LBFGSpp

`costGrad(pb, x, g)` 的签名就是"给定 x 返回 f 并写入梯度"，可以直接包一层 Eigen::VectorXd 适配器换成 LBFGS-Lite（ZJU-FAST-Lab）或 LBFGSpp。建议保留"最大步长"限制：LBFGS-Lite 对应 `max_linesearch` / 步长上界相关参数，或者在回调里自行限制。

## 已知局限与扩展点

- **前端不考虑曲率**：前端是无曲率约束的切线图，只用起终点"引导段"（长度 = 转弯半径）近似处理航向。若终点航向需要大角度掉头，可将起终点左右转弯圆作为有向虚拟球冠加入切线图（Dubins 风格）。
- **碰撞约束是软约束**：障碍为罚函数，最终由 `validate()` 硬校验把关。如果需要可证明的安全性，可改为"每个 span 的半球约束 a·P_j ≤ 0"加 QP（OSQP）的 SCP 版本。
- **同伦签名的定义域**：签名在"球面减去极点"上定义。默认极点取初始起终点中点的对跖点，作业区域不应接近该点；也可以用 `setHomotopyPole()` 把极点放进某个障碍内部。
- **重规划起点**：`replan` 的 `advanceDist` 用于补偿计算和通信延迟，通常取 v·延迟。
