> [← 返回总目录](MCM_A_备赛学习指南.md)　·　本册：附录 A—F

## 附录 A 符号总表

> 全书统一使用下表符号。遇到与章节内特殊记号冲突时，以该章节的"符号说明"为准。

| 符号 | 含义 | 常用单位 | 首次出现 |
|---|---|---|---|
| $t$ | 时间（自变量） | s, h, d, yr | 2.1 |
| $x,y,z$ | 空间坐标 | m, km | 3.5 |
| $\mathbf{x}$ | 状态向量 | 依系统 | 2.8 |
| $x_i$ | 第 $i$ 个状态分量 | 依系统 | 2.8 |
| $u$ | 控制/决策变量 | 依系统 | 3.9, 4.1 |
| $\boldsymbol\theta$ | 参数向量 | 依模型 | 3.7 |
| $y$ | 观测量 | 依实验 | 3.7 |
| $N$ | 种群数量（个体数） | 个 | 3.2 |
| $M$ | 生物量/体重 | kg | 3.1, 9.1 |
| $S,I,R$ | 易感/感染/恢复 | 人 | 5.2 |
| $P$ | 功率 | W | 3.1, 9.4 |
| $v$ | 速度 | m/s | 3.1 |
| $s$ | 荷电状态 SOC | 无量纲 | 2.13, 9.8 |
| $T$ | 温度 | °C, K | 2.13 |
| $W$ | 含水量 | 无量纲, mm | 2.14, 9.5 |
| $L$ | 枯落物质量/长度 | kg, m | 9.3 |
| $h$ | 深度/厚度/步长 | m, — | 2.14, 3.6 |
| $r$ | 内禀增长率 | 1/time | 3.2 |
| $K$ | 环境容量 | 个体 | 3.2 |
| $q$ | 捕捞系数 | 1/(effort·time) | 3.2 |
| $E$ | 捕捞努力量/能量 | h, J | 3.2, 9.1 |
| $\alpha,\beta,\gamma$ | 模型系数（含义随章节说明） | — | 各处 |
| $\lambda$ | 特征值/拉格朗日乘子/失效率 | — | 3.3, 3.9, 4.1, 5.2 |
| $\sigma$ | 标准差/噪声强度 | 依量 | 3.7 |
| $\tau$ | 时间常数/延迟 | time | 3.1, 4.9 |
| $\rho$ | 密度/相关系数/利用率 | kg/m³, — | 2.9, 3.7 |
| $\eta$ | 效率/尺度参数 | — | 3.1, 5.2 |
| $\delta$ | 磨损深度/扰动 | m | 2.14 |
| $k$ | 速率常数 | 1/time | 3.2 |
| $D$ | 扩散系数 | m²/s | 3.5 |
| $J$ | 目标函数 | 依问题 | 3.9, 4.1 |
| $H$ | 哈密顿函数/硬度 | —, Pa | 4.1, 2.14 |
| $\pi_i$ | 无量纲群 | — | 2.9 |
| $\mathrm{Re},\mathrm{Pe},\mathrm{Da},\mathrm{Bi}$ | 雷诺/佩克莱/达姆科勒/毕渥数 | — | 3.13 |

---

## 附录 B 核心公式速查

**守恒与输运**

$$
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0,\qquad
\frac{d}{dt}\int_V \rho\,\phi\,dV=-\oint_{\partial V}\mathbf J\cdot\mathbf n\,dS+\int_V s\,dV.
$$

**一维守恒型 PDE（对流—扩散—反应）**

$$
\frac{\partial c}{\partial t}+v\frac{\partial c}{\partial x}=D\frac{\partial^2 c}{\partial x^2}+R(c).
$$

**Logistic + 捕捞**

$$
\dot N=rN\left(1-\frac NK\right)-qEN,\quad N^*=K\left(1-\frac{qE}{r}\right),\quad E_c=\frac rq .
$$

**Ricker / Beverton–Holt**

$$
N_{t+1}=N_t e^{r(1-N_t/K)},\qquad N_{t+1}=\frac{RN_t}{1+(R-1)N_t/K}.
$$

**Leslie 矩阵**

$$
\mathbf n_{t+1}=L\mathbf n_t,\qquad \lambda_1=\text{谱半径}(L).
$$

**热方程（一维）与显式格式**

$$
\frac{\partial T}{\partial t}=\alpha\frac{\partial^2 T}{\partial x^2},\qquad
T_i^{n+1}=T_i^n+r\left(T_{i+1}^n-2T_i^n+T_{i-1}^n\right),\quad r=\frac{\alpha\Delta t}{\Delta x^2}\le\frac12 .
$$

**PMP**

$$
H=L+\lambda^{\mathsf T}f,\quad \dot x=H_\lambda,\quad \dot\lambda=-H_x,\quad H_u=0,\quad \lambda(t_f)=S_x.
$$

**Kalman（标量）**

$$
K=\frac{P^-}{P^-+R},\quad \hat x=\hat x^-+K(z-\hat x^-),\quad P=(1-K)P^-.
$$

**普通 Kriging**

$$
\begin{bmatrix}\Gamma&\mathbf 1\\ \mathbf 1^{\mathsf T}&0\end{bmatrix}
\begin{bmatrix}\mathbf w\\ \mu\end{bmatrix}=
\begin{bmatrix}\boldsymbol\gamma_0\\ 1\end{bmatrix},\qquad
\hat Z(\mathbf x_0)=\mathbf w^{\mathsf T}\mathbf Z,\quad
\sigma^2=\mathbf w^{\mathsf T}\boldsymbol\gamma_0+\mu .
$$

**Sobol（Saltelli 估计）**

$$
\hat S_i=\frac{\frac1N\sum y_A y_{A_B^{(i)}}-f_0^2}{\hat V},\qquad
\hat S_{Ti}=1-\frac{\frac1N\sum y_A y_{A_B^{(-i)}}-f_0^2}{\hat V}.
$$

**贝叶斯反演**

$$
p(\boldsymbol\theta\mid \mathbf y)\propto p(\mathbf y\mid\boldsymbol\theta)\,p(\boldsymbol\theta),\qquad
\hat{\boldsymbol\theta}_{\text{MAP}}=\arg\max \ln p(\mathbf y\mid\boldsymbol\theta)+\ln p(\boldsymbol\theta).
$$

**可辨识性（Fisher 信息）**

$$
\mathcal I(\boldsymbol\theta)=\sum_i\frac{1}{\sigma_i^2}\left(\frac{\partial G(\mathbf x_i;\boldsymbol\theta)}{\partial\boldsymbol\theta}\right)^{\mathsf T}
\left(\frac{\partial G(\mathbf x_i;\boldsymbol\theta)}{\partial\boldsymbol\theta}\right).
$$

**量纲分析（Buckingham π）**

$$
n_{\pi}=n_{\text{var}}-\operatorname{rank}(\text{量纲矩阵}).
$$

**Archard 磨损**

$$
\delta=\frac{K F_n s}{H A_c}.
$$

---

## 附录 C MATLAB 函数索引（建模手视角）

| 函数 | 用途 | 数学对应 | 章节 |
|---|---|---|---|
| `ode45` | 非刚性 ODE 初值问题 | RK45 自适应 | 3.2, 8.2 |
| `ode15s` | 刚性 ODE/DAE | BDF | 3.6 |
| `odeset` | 事件、容差、雅可比 | 事件检测 | 3.2, 8.2 |
| `bvp4c` | 两点边值问题 | 打靶/配置法 | 4.1 |
| `pdepe` | 一维抛物/椭圆 PDE | 有限体积/线法 | 3.5 |
| `pdepe` 的 `bcfun` | 边界条件 | 三类边界 | 3.5 |
| `fmincon` | 约束非线性优化 | SQP/内点 | 3.9 |
| `linprog` | 线性规划 | 单纯形/内点 | 3.9 |
| `intlinprog` | 混合整数线性规划 | 分支定界 | 5.3 |
| `ga` / `particleswarm` | 全局启发式 | 元启发式 | 5.5 |
| `lsqcurvefit` | 非线性最小二乘拟合 | Gauss–Newton/LM | 3.7 |
| `lsqnonlin` | 非线性最小二乘残差 | 同上 | 4.3 |
| `nlparci` | 参数置信区间 | 线性化协方差 | 3.7 |
| `fsolve` | 非线性方程组 | Newton/信赖域 | 3.6 |
| `fzero` | 单变量求根 | 混合法 | 3.6 |
| `integral` | 数值积分 | 自适应 Simpson | 3.6 |
| `interp1` | 一维插值 | 线性/样条 | 4.7 |
| `scatteredInterpolant` | 散点插值 | Delaunay/自然邻域 | 4.7 |
| `rand` / `randn` | 均匀/正态随机数 | MC 采样 | 3.11 |
| `lhsdesign` | 拉丁超立方 | 分层采样 | 3.11 |
| `sobolset` / `qmc` | 低差异序列 | 准 MC | 3.11 |
| `gamultiobj` | 多目标遗传算法 | Pareto 前沿 | 4.8 |
| `paretosearch` | 多目标模式搜索 | 前沿 | 4.8 |
| `fitrgp` | 高斯过程回归 | GP/代理 | 5.5 |
| `bayesopt` | 贝叶斯优化 | EI/UCB | 5.5 |
| `graph` / `shortestpath` | 图与最短路 | Dijkstra | 5.3 |
| `maxflow` | 最大流 | 最小割 | 5.3 |
| `fitdist` | 分布拟合 | MLE | 5.2 |
| `arima` | 时间序列模型 | ARIMA | 4.10 |
| `eig` | 特征值 | 稳定性/谱半径 | 3.3 |
| `svd` | 奇异值分解 | POD/PCA/条件数 | 5.4 |
| `cond` | 条件数 | 可辨识性诊断 | 3.7 |
| `jacobian`（符号） | 雅可比 | 灵敏度 | 3.10 |

---

## 附录 D 教材与可靠资料

**只列真实存在的书与官方入口；不提供可能失效的第三方链接。**

**数学建模总论**

- 姜启源、谢金星、叶俊.《数学模型》（第五版）. 高等教育出版社.
- 司守奎、孙玺菁.《数学建模算法与应用》. 国防工业出版社.

**微分方程与动力系统**

- Boyce, W. E., & DiPrima, R. C. *Elementary Differential Equations and Boundary Value Problems*. Wiley.
- Strogatz, S. H. *Nonlinear Dynamics and Chaos*. Westview Press.
- Murray, J. D. *Mathematical Biology I: An Introduction*. Springer.

**PDE 与数值方法**

- Logan, J. D. *Applied Partial Differential Equations*. Springer.
- LeVeque, R. J. *Finite Difference Methods for Ordinary and Partial Differential Equations*. SIAM.
- Hairer, E., & Wanner, G. *Solving Ordinary Differential Equations I/II*. Springer.
- Nocedal, J., & Wright, S. J. *Numerical Optimization*. Springer.

**统计、参数估计与 UQ**

- Montgomery, D. C., & Runger, G. C. *Applied Statistics and Probability for Engineers*. Wiley.
- Saltelli, A., et al. *Global Sensitivity Analysis: The Primer*. Wiley.
- Efron, B., & Tibshirani, R. J. *An Introduction to the Bootstrap*. Chapman & Hall.
- Rasmussen, C. E., & Williams, C. K. I. *Gaussian Processes for Machine Learning*. MIT Press（作者主页有合法公开版）.

**优化与控制**

- Boyd, S., & Vandenberghe, L. *Convex Optimization*. Cambridge University Press（作者主页有合法公开版）.
- Bertsekas, D. P. *Dynamic Programming and Optimal Control*. Athena Scientific.
- Liberzon, D. *Calculus of Variations and Optimal Control Theory*. Princeton University Press.
- Ben-Tal, A., El Ghaoui, L., & Nemirovski, A. *Robust Optimization*. Princeton University Press.

**随机与时间序列**

- Øksendal, B. *Stochastic Differential Equations*. Springer.
- Shumway, R. H., & Stoffer, D. S. *Time Series Analysis and Its Applications*. Springer.
- Hyndman, R. J., & Athanasopoulos, G. *Forecasting: Principles and Practice*. OTexts（官方免费在线版）.

**复杂系统与仿真**

- Railsback, S. F., & Grimm, V. *Agent-Based and Individual-Based Modeling*. Princeton University Press.
- Forrester, J. W. *Industrial Dynamics*. MIT Press.

**官方入口（可自行检索，不依赖具体链接）**

- COMAP（MCM/ICM 主办方）：`comap.com`，历年题目与优秀论文官方入口。
- MathWorks 文档：`mathworks.com/help`，函数用法与示例。
- 本仓库：`MCM_A_2020-2026_历年分析报告.md`，2020—2026 年题目的逐年分析。

---

## 附录 E 建模交底模板（可直接复制）

```markdown
## 建模交底：<题目简称>

### 1. 问题与输出
- 题目要的输出（数值/曲线/场/策略/建议）：
- 我建的模型回答的是哪一问：

### 2. 变量表
| 符号 | 类型 | 含义 | 单位 | 量级/来源 |
|---|---|---|---|---|
| | 状态/控制/参数 | | | |

### 3. 方程
- 主方程（每条附一句"从哪条定律来"）：
- 初始条件：
- 边界条件/终端条件：
- 辅助代数关系：

### 4. 假设与理由
| 假设 | 理由 | 若被违反的后果 |
|---|---|---|

### 5. 参数来源
| 参数 | 取值 | 来源（文献/拟合/假设） | 不确定度 |
|---|---|---|---|

### 6. 求解方案
- 解析/数值：
- 离散格式与步长：
- 需要编程实现的部分（输入/输出接口）：

### 7. 验证清单
- [ ] 量纲一致
- [ ] 守恒（质量/能量/概率）
- [ ] 极限行为
- [ ] 网格/步长收敛
- [ ] 与解析解或文献对比
- [ ] 参数可辨识性

### 8. 可信度与边界
- 关键参数及其区间：
- 结论在什么范围内成立：
- 最可能出错的地方：
```

**给建模手的一句话：** 这张表填不满，就不要让编程手开始写代码——**空着的格子最后都会变成论文里的漏洞**。

---

## 附录 F 原文档 36 个模型去向索引

| 原编号 | 原模型 | 新位置 | 新级别 | 处理说明 |
|---|---|---|---|---|
| 1 | 机理建模与守恒 | 3.1 | P0 | 保留并加 Reynolds 输运推导 |
| 2 | ODE | 3.2, 3.3 | P0 | 拆出稳定性与分岔 |
| 3 | 概率统计与回归 | 3.8, 3.7 | P0 | 回归归 3.8，估计归 3.7 |
| 4 | 数值方法与蒙特卡洛 | 3.6, 4.4 | P0/P1 | 数值归 P0，MC 归 P1 |
| 5 | 优化（LP/NLP/DP） | 3.9, 4.2 | P0/P1 | DP 独立为 P1 |
| 6 | 量纲分析与相似 | 2.9, 3.13 | P0 | 基础进方法论，进阶留 3.13 |
| 7 | 敏感性与 UQ | 3.10, 3.11 | P0 | 拆成敏感性与 UQ 两节 |
| 8 | PDE | 3.5 | **P0（升级）** | 原文列 P1，与 A 题权重不符 |
| 9 | 时间序列与随机过程 | 4.10, 4.4 | P1 | 拆分 |
| 10 | 空间插值与 Kriging | 4.7 | P1 | 保留 |
| 11 | 降维 PCA/POD | 5.4 | P2 | 降级：A 题少单独考 |
| 12 | 评价 AHP/TOPSIS/模糊 | 5.1 | P2 | 降级：非机理主线 |
| 13 | CA/ABM/DES | 4.5, 4.6 | P1 | CA/ABM 保留，DES 并入 5.3 |
| 14 | Leslie | 3.4, 5.2 | P0 例/P2 | 作为差分例 + P2 条目 |
| 15 | Beverton–Holt/Ricker | 3.4, 5.2 | P0 例/P2 | 同上 |
| 16 | 最优控制 PMP | 4.1 | **P1（升级）** | 原文标"选学"，实为 2022 核心 |
| 17 | 逆问题与贝叶斯反演 | 4.3 | **P1（升级）** | 2025 核心能力 |
| 18 | 生存分析 Weibull | 5.2 | P2 | 降级：A 题偶用 |
| 19 | 排队论/DES | 5.3 | P2 | 降级：与 A 题弱相关 |
| 20 | 图论/最短路/网络流 | 5.3 | P2 | 降级：工具性，用时查 |
| 21 | 阈值/分段/混合动态 | 5.3 | P2 | 保留思想，事件在 3.2 |
| 22 | 多目标 Pareto | 4.8 | P1 | 保留 |
| 23 | SIR/SEIR | 5.2 | **P2（降级）** | 用户指定 P2，A 题少直接考 |
| 24 | 指数平滑/Holt–Winters | 5.1 | P2 | 降级：预测非主线 |
| 25 | 聚类 K-means/GMM | 5.4 | P2 | 降级，修正轮廓系数 |
| 26 | GM(1,1) | 5.1 | P2 | 降级，修正数值 |
| 27 | MDP | 4.2 | P1 | 保留（与 DP 合并） |
| 28 | 博弈 Nash/演化 | 5.3 | P2 | 降级 |
| 29 | 分层贝叶斯 | 5.4 | P2 | 降级，修正收缩公式 |
| 30 | 机器学习 | 5.1 | P2 | 降级为"只作基线/代理" |
| 31 | 系统动力学 | 4.9 | P1 | 保留 |
| 32 | 元启发式 GA/SA/PSO/DE | 5.5 | P2 | 降级：备用 |
| 33 | Kalman/数据同化 | 5.4 | P2 | 保留 P2 |
| 34 | SDE | 5.5 | P2 | 保留 P2 |
| 35 | 鲁棒/随机规划 | 5.3 | P2 | 保留 P2 |
| 36 | 代理模型/贝叶斯优化 | 5.5 | P2 | 保留 P2 |

**新增（原文缺，A 题需要）**

| 新增内容 | 位置 | 理由 |
|---|---|---|
| 稳定性与分岔 | 3.3 | 判断平衡点是否可信，是建模手的基本功 |
| 差分方程独立成节 | 3.4 | 原文散落在 Leslie/Ricker 中，缺统一框架 |
| 参数估计系统化 | 3.7 | 原文只有零散 nlinfit 片段 |
| 模型验证 | 3.12 | 原文验证只有 5 行清单 |
| 不确定性分析 | 3.11 | 原文与敏感性混在一起 |
| 选择/比较/组合/改进 | 第六章 | 原文只有 3 页速查 |
| 工作流与可信度声明 | 第七章 | 把"验证"变成可执行流程 |
| 真题逐年拆解 | 第九章 | 原文只有零散"应用"段 |
| 四周计划 | 第十章 | 原文是 12—14 周，与一个月冲突 |
| P2 深做路线与五个示范 | 5.7—5.8 | 补上"P2 太浅、需另找教材"的缺口：21 项各有推导/实现/读物/验收，SDE、分层贝叶斯、鲁棒、生存、Kalman 给完整推导 |

---

*文档版本：2026-11 重构版。原版（6695 行）见本仓库 git 历史中的上一提交，可用 `git show HEAD:MCM_A_备赛学习指南.md` 查看。*

